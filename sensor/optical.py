"""A ground-based optical sensor for the geostationary belt.

Hardware (docs/01-physical.md): a 200 mm f/2.8 lens on an astronomy camera
(3008 x 3008 px, 3.76 um pixels: 3.88 arcsec/px, 3.2 deg field) on a mount
tracking the stars. In a 2 s exposure stars are round and a GEO satellite,
which stays fixed against the ground, draws a short streak (about 30 arcsec).

Pipeline for one frame:
  1. detect sources (background, threshold, connected components, moments)
  2. split round sources (stars) from elongated ones along the expected
     satellite direction (streaks)
  3. astrometry: match stars to the catalogue from the mount's approximate
     pointing (translation vote, then nearest neighbour) and fit a gnomonic
     plate solution by least squares
  4. measure each streak's centre (mid-exposure position) in RA/Dec
  5. correlate with catalogue predictions; residuals; flag uncorrelated tracks
     (UCTs: possibly unknown objects) and large residuals (manoeuvre or stale
     elements)

render() makes synthetic frames for rehearsal: catalogue stars plus satellites
at positions the caller supplies (in the rehearsal: the operators' own orbits),
with a Gaussian PSF, sky background, photon and read noise.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import ndimage

ARCSEC = math.pi / (180 * 3600)


@dataclass
class Camera:
    width: int = 3008
    height: int = 3008
    scale_arcsec: float = 3.88          # arcsec per pixel
    psf_sigma_px: float = 1.1
    exposure_s: float = 2.0
    zero_point: float = 21.5            # magnitude giving 1 e-/s
    sky_e_per_px: float = 120.0         # suburban sky, 2 s
    read_noise_e: float = 3.5
    full_well_e: float = 50000.0
    distortion_k1: float = 2.2e-10        # radial: r' = r (1 + k1 r^2); about 0.05% (0.75 px) at the corner
    seeing_jitter_px: float = 0.25        # image motion from the atmosphere in a 2 s exposure (about 1 arcsec)
    catalogue_error_arcsec: float = 0.3   # star catalogue position error (Tycho-class; Gaia is far better)


# ------------------------------------------------------------------ tangent plane

def project(ra, dec, ra0, dec0):
    """Gnomonic projection: (ra, dec) radians -> standard coordinates (xi, eta) radians."""
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return xi, eta


def deproject(xi, eta, ra0, dec0):
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / np.where(rho > 0, rho, 1), 0))
        ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return ra % (2 * np.pi), dec


# ------------------------------------------------------------------ synthetic stars

def star_field(ra0, dec0, radius_deg, seed, density_per_deg2=160.0, mag_limit=12.5):
    """A deterministic synthetic star catalogue around a pointing (the rehearsal's sky).
    Counts follow log N ~ 0.45 m, roughly the real sky away from the galactic plane."""
    rng = np.random.default_rng(seed)
    area = math.pi * radius_deg ** 2
    n = rng.poisson(density_per_deg2 * area)
    u = rng.random(n)
    mag = mag_limit + np.log10(np.maximum(u, 1e-6)) / 0.45
    r = np.radians(radius_deg) * np.sqrt(rng.random(n))
    th = rng.random(n) * 2 * np.pi
    ra, dec = deproject(r * np.cos(th), r * np.sin(th), ra0, dec0)
    keep = mag > 3
    return ra[keep], dec[keep], mag[keep]


# ------------------------------------------------------------------ rendering

def render(cam: Camera, ra0, dec0, roll, stars, sats, seed=0):
    """stars: (ra, dec, mag); sats: list of dict(ra, dec, mag, rate_ra, rate_dec) with rates in rad/s
    of apparent motion relative to the stars. Returns image (float32) and the true pixel positions."""
    rng = np.random.default_rng(seed)
    img = np.full((cam.height, cam.width), cam.sky_e_per_px, dtype=np.float32)
    cr, sr = math.cos(roll), math.sin(roll)

    def to_px(ra, dec, jitter=True):
        xi, eta = project(ra, dec, ra0, dec0)
        u = (cr * xi - sr * eta) / (cam.scale_arcsec * ARCSEC)
        v = (sr * xi + cr * eta) / (cam.scale_arcsec * ARCSEC)
        k = 1 + cam.distortion_k1 * (u * u + v * v)              # lens distortion
        x, y = u * k + cam.width / 2, v * k + cam.height / 2
        if jitter:
            x = x + rng.normal(0, cam.seeing_jitter_px, np.shape(x))
            y = y + rng.normal(0, cam.seeing_jitter_px, np.shape(y))
        return x, y

    def stamp(x, y, flux, dx=0.0, dy=0.0, steps=1):
        s = cam.psf_sigma_px
        r = int(4 * s + abs(dx) + abs(dy)) + 2
        xi, yi = int(round(x)), int(round(y))
        x0, x1, y0, y1 = max(0, xi - r), min(cam.width, xi + r + 1), max(0, yi - r), min(cam.height, yi + r + 1)
        if x0 >= x1 or y0 >= y1:
            return
        gx, gy = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1))
        acc = np.zeros_like(gx, dtype=np.float64)
        for k in range(steps):
            f = (k + 0.5) / steps - 0.5
            acc += np.exp(-((gx - (x + f * dx)) ** 2 + (gy - (y + f * dy)) ** 2) / (2 * s * s))
        acc *= flux / (2 * math.pi * s * s * steps)
        img[y0:y1, x0:x1] += acc.astype(np.float32)

    sx, sy = to_px(stars[0], stars[1])
    flux = 10 ** (-0.4 * (stars[2] - cam.zero_point)) * cam.exposure_s
    for x, y, f in zip(sx, sy, flux):
        if -20 < x < cam.width + 20 and -20 < y < cam.height + 20:
            stamp(x, y, f)
    truth = []
    for s in sats:
        x, y = to_px(np.array([s["ra"]]), np.array([s["dec"]]))
        x2, y2 = to_px(np.array([s["ra"] + s["rate_ra"] * cam.exposure_s]), np.array([s["dec"] + s["rate_dec"] * cam.exposure_s]), jitter=False)
        x2, y2 = x2 + (x - to_px(np.array([s["ra"]]), np.array([s["dec"]]), jitter=False)[0]), y2 + (y - to_px(np.array([s["ra"]]), np.array([s["dec"]]), jitter=False)[1])
        dx, dy = float(x2[0] - x[0]), float(y2[0] - y[0])
        cx, cy = float(x[0] + dx / 2), float(y[0] + dy / 2)      # mid-exposure
        f = 10 ** (-0.4 * (s["mag"] - cam.zero_point)) * cam.exposure_s
        if 0 <= cx < cam.width and 0 <= cy < cam.height:
            stamp(cx, cy, f, dx, dy, steps=max(1, int(math.hypot(dx, dy) * 2)))
            truth.append({**s, "x": cx, "y": cy, "dx": dx, "dy": dy})
    img = rng.poisson(np.clip(img, 0, None)).astype(np.float32) + rng.normal(0, cam.read_noise_e, img.shape).astype(np.float32)
    return np.minimum(img, cam.full_well_e), truth


# ------------------------------------------------------------------ detection

def detect(img, nsigma=5.0, min_px=4):
    med = float(np.median(img))
    mad = float(np.median(np.abs(img - med))) * 1.4826
    sm = ndimage.gaussian_filter(img - med, 1.0)
    lab, n = ndimage.label(sm > nsigma * mad / 2.0)
    out = []
    objs = ndimage.find_objects(lab)
    for k, sl in enumerate(objs, start=1):
        m = lab[sl] == k
        if m.sum() < min_px:
            continue
        sub = (img[sl] - med) * m
        sub = np.clip(sub, 0, None)
        tot = sub.sum()
        if tot <= 0:
            continue
        yy, xx = np.mgrid[sl[0], sl[1]]
        cx, cy = (sub * xx).sum() / tot, (sub * yy).sum() / tot
        mxx = (sub * (xx - cx) ** 2).sum() / tot
        myy = (sub * (yy - cy) ** 2).sum() / tot
        mxy = (sub * (xx - cx) * (yy - cy)).sum() / tot
        tr, det = mxx + myy, mxx * myy - mxy ** 2
        l1 = tr / 2 + math.sqrt(max(tr * tr / 4 - det, 0))
        l2 = tr / 2 - math.sqrt(max(tr * tr / 4 - det, 0))
        out.append({"x": float(cx), "y": float(cy), "flux": float(tot), "snr": float(tot / (mad * math.sqrt(m.sum()) + 1e-9)),
                    "elong": float(math.sqrt(l1 / max(l2, 1e-6))), "angle": float(0.5 * math.atan2(2 * mxy, mxx - myy)),
                    "length_px": float(math.sqrt(12 * max(l1 - l2, 0)))})
    return out, mad


# ------------------------------------------------------------------ astrometry

def _terms(u, v, order):
    """Polynomial terms in normalised pixel offsets (1, u, v, u2, uv, v2, ...)."""
    u, v = np.asarray(u, float) / 1500.0, np.asarray(v, float) / 1500.0
    cols = [np.ones_like(u)]
    for n in range(1, order + 1):
        for i in range(n + 1):
            cols.append(u ** (n - i) * v ** i)
    return np.column_stack(cols)


@dataclass
class Plate:
    ra0: float
    dec0: float
    a_xi: np.ndarray        # polynomial coefficients: pixel offsets -> standard coordinates (radians)
    a_eta: np.ndarray
    cx: float
    cy: float
    order: int
    rms_arcsec: float
    n_stars: int

    def to_radec(self, x, y):
        T = _terms(np.atleast_1d(x) - self.cx, np.atleast_1d(y) - self.cy, self.order)
        ra, dec = deproject(T @ self.a_xi, T @ self.a_eta, self.ra0, self.dec0)
        return (ra, dec) if np.ndim(x) else (float(ra[0]), float(dec[0]))

    def to_px(self, ra, dec):
        """Inverse by Newton iteration on the polynomial (for star subtraction)."""
        xi, eta = project(np.asarray(ra), np.asarray(dec), self.ra0, self.dec0)
        lin = np.array([[self.a_xi[1], self.a_xi[2]], [self.a_eta[1], self.a_eta[2]]]) / 1500.0
        inv = np.linalg.inv(lin)
        uv = inv @ np.vstack([xi - self.a_xi[0], eta - self.a_eta[0]])
        for _ in range(4):
            T = _terms(uv[0], uv[1], self.order)
            r = np.vstack([T @ self.a_xi - xi, T @ self.a_eta - eta])
            uv = uv - inv @ r
        return uv[0] + self.cx, uv[1] + self.cy


def solve(dets, cat_ra, cat_dec, cat_mag, ra_guess, dec_guess, cam: Camera, roll_guess=0.0, tol_px=3.0, order=3):
    """Plate solution from round sources and the catalogue, starting from the mount's pointing."""
    stars = sorted([d for d in dets if d["elong"] < 1.6], key=lambda d: -d["flux"])[:250]
    if len(stars) < 8:
        return None
    s = cam.scale_arcsec * ARCSEC
    cr, sr = math.cos(roll_guess), math.sin(roll_guess)
    xi, eta = project(cat_ra, cat_dec, ra_guess, dec_guess)
    px = (cr * xi - sr * eta) / s + cam.width / 2
    py = (sr * xi + cr * eta) / s + cam.height / 2
    bright = np.argsort(cat_mag)[:250]
    # 1. translation vote: offsets between every bright detection and every bright predicted star
    dx = np.array([d["x"] for d in stars])[:, None] - px[bright][None, :]
    dy = np.array([d["y"] for d in stars])[:, None] - py[bright][None, :]
    lim = 600
    m = (np.abs(dx) < lim) & (np.abs(dy) < lim)
    h, xe, ye = np.histogram2d(dx[m], dy[m], bins=240, range=[[-lim, lim], [-lim, lim]])
    i, j = np.unravel_index(np.argmax(h), h.shape)
    ox, oy = (xe[i] + xe[i + 1]) / 2, (ye[j] + ye[j + 1]) / 2
    # 2. nearest-neighbour matches, then least squares, twice
    det_xy = np.array([[d["x"], d["y"]] for d in dets if d["elong"] < 1.6])
    cat_xy = np.stack([px + ox, py + oy], 1)
    from scipy.spatial import cKDTree
    plate = None
    for it in range(3):
        tree = cKDTree(cat_xy)
        dist, idx = tree.query(det_xy, distance_upper_bound=tol_px * (4 if it == 0 else 1.5))
        ok = np.isfinite(dist)
        if ok.sum() < 8:
            return None
        dpx, cidx = det_xy[ok], idx[ok]
        ordr = 1 if it == 0 else order
        T = _terms(dpx[:, 0] - cam.width / 2, dpx[:, 1] - cam.height / 2, ordr)
        a_xi, *_ = np.linalg.lstsq(T, xi[cidx], rcond=None)
        a_eta, *_ = np.linalg.lstsq(T, eta[cidx], rcond=None)
        # move the tangent point to the frame centre, then re-express the fit there
        ra_c, dec_c = deproject(np.array([a_xi[0]]), np.array([a_eta[0]]), ra_guess, dec_guess)
        ra_guess, dec_guess = float(ra_c[0]), float(dec_c[0])
        xi, eta = project(cat_ra, cat_dec, ra_guess, dec_guess)
        a_xi, *_ = np.linalg.lstsq(T, xi[cidx], rcond=None)
        a_eta, *_ = np.linalg.lstsq(T, eta[cidx], rcond=None)
        res = np.column_stack([T @ a_xi - xi[cidx], T @ a_eta - eta[cidx]])
        rms = float(np.sqrt(np.mean((res ** 2).sum(1)))) / ARCSEC
        plate = Plate(ra_guess, dec_guess, a_xi, a_eta, cam.width / 2, cam.height / 2, ordr, rms, int(ok.sum()))
        cat_xy = np.column_stack(plate.to_px(cat_ra, cat_dec))
    return plate


def subtract_stars(img, plate: Plate, cat_ra, cat_dec, cat_mag, dets, cam: Camera):
    """Remove every catalogue star from the image using the plate solution and a photometric
    zero point fitted on matched stars, so streaks that overlap stars can be measured."""
    from scipy.spatial import cKDTree
    sx, sy = plate.to_px(cat_ra, cat_dec)
    inside = (sx > -10) & (sx < cam.width + 10) & (sy > -10) & (sy < cam.height + 10)
    sx, sy, mag = sx[inside], sy[inside], cat_mag[inside]
    stars = [d for d in dets if d["elong"] < 1.6]
    tree = cKDTree(np.column_stack([sx, sy]))
    zps = []
    for d in stars:
        dist, k = tree.query([d["x"], d["y"]])
        if dist < 1.5 and d["flux"] > 0:
            zps.append(2.5 * math.log10(d["flux"] / cam.exposure_s) + mag[k])
    zp = float(np.median(zps)) if zps else cam.zero_point
    model = np.zeros_like(img)
    s = cam.psf_sigma_px
    r = int(4 * s) + 2
    for x, y, m in zip(sx, sy, mag):
        f = 10 ** (-0.4 * (m - zp)) * cam.exposure_s
        xi, yi = int(round(x)), int(round(y))
        x0, x1, y0, y1 = max(0, xi - r), min(cam.width, xi + r + 1), max(0, yi - r), min(cam.height, yi + r + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        gx, gy = np.meshgrid(np.arange(x0, x1), np.arange(y0, y1))
        model[y0:y1, x0:x1] += f / (2 * math.pi * s * s) * np.exp(-((gx - x) ** 2 + (gy - y) ** 2) / (2 * s * s))
    return img - model, zp


# ------------------------------------------------------------------ satellites

def streaks(dets, expected_angle=None, min_elong=1.8, expected_len_px=None, star_xy=None, star_mag=None,
            bright_mag=8.0, bright_mask_px=12.0, star_mask_px=2.5, star_flux=None, residual_ratio=0.35):
    """Satellite candidates: elongated along the expected direction, about the expected length
    (known from the tracking mode), and not explained by a star (saturated bright stars leave
    residuals after subtraction)."""
    out = []
    for d in dets:
        if d["elong"] < min_elong:
            continue
        if expected_angle is not None:
            da = abs((d["angle"] - expected_angle + math.pi / 2) % math.pi - math.pi / 2)
            if da > math.radians(15):
                continue
        if expected_len_px is not None and not (0.5 * expected_len_px <= d["length_px"] <= 1.8 * expected_len_px):
            continue
        if star_xy is not None and len(star_xy):
            dist = np.hypot(star_xy[:, 0] - d["x"], star_xy[:, 1] - d["y"])
            near_bright = (dist < bright_mask_px) & (star_mag < bright_mag)
            if near_bright.any() or (dist < star_mask_px).any():
                continue
            # residual of imperfectly subtracted stars: a small fraction of the starlight around it
            if star_flux is not None:
                around = star_flux[dist < 6.0].sum()
                if around > 0 and d["flux"] < residual_ratio * around:
                    continue
        out.append(d)
    return out


def correlate(measured, predicted, gate_arcsec=300.0):
    """measured: [(ra, dec)]; predicted: [{norad, ra, dec}]. Nearest prediction inside the gate."""
    res = []
    used = set()
    for k, (ra, dec) in enumerate(measured):
        best, bd = None, 1e9
        for p in predicted:
            d = math.degrees(math.acos(min(1.0, math.sin(dec) * math.sin(p["dec"]) + math.cos(dec) * math.cos(p["dec"]) * math.cos(ra - p["ra"])))) * 3600
            if d < bd:
                best, bd = p, d
        if best is not None and bd <= gate_arcsec and best["norad"] not in used:
            used.add(best["norad"])
            dra = ((ra - best["ra"] + math.pi) % (2 * math.pi) - math.pi) * math.cos(dec) / ARCSEC
            res.append({"k": k, "norad": best["norad"], "sep_arcsec": bd, "dra_arcsec": dra, "ddec_arcsec": (dec - best["dec"]) / ARCSEC})
        else:
            res.append({"k": k, "norad": None, "sep_arcsec": bd})
    return res


def solve_checked(dets, cat_ra, cat_dec, cat_mag, ra_guess, dec_guess, cam: Camera, max_rms=3.0, min_fraction=0.5):
    """Plate solution with a quality gate. A frame whose solution matches too few of its stars or fits
    them badly would put every satellite measurement in it off by the same error, so it is retried
    with a coarser start and otherwise rejected (returns None, reason)."""
    n_round = sum(1 for d in dets if d["elong"] < 1.6)
    tries = [dict(), dict(tol_px=6.0), dict(tol_px=2.0)]
    best = None
    for kw in tries:
        p = solve(dets, cat_ra, cat_dec, cat_mag, ra_guess, dec_guess, cam, **kw)
        if p is None:
            continue
        if best is None or p.rms_arcsec < best.rms_arcsec:
            best = p
        if p.rms_arcsec <= max_rms and p.n_stars >= min_fraction * n_round:
            return p, "ok"
    if best is None:
        return None, "no solution"
    return None, f"rejected: {best.n_stars} of {n_round} stars, {best.rms_arcsec:.1f} arcsec"
