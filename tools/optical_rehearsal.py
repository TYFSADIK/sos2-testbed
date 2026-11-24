#!/usr/bin/env python3
"""A night of GEO-belt observing from Toronto, rehearsed with synthetic frames.

For each field (centred on a Canadian GEO slot):
  * the satellites in the frame are drawn where they really are: at the
    operators' own orbits (CelesTrak SupGP: Telesat, Intelsat, SES, EUMETSAT)
    when available, otherwise at the catalogue position plus an error drawn
    from the measured uncertainty model (labelled "modelled truth");
  * two tests are injected: an object that is not in the catalogue, and a
    catalogued satellite displaced 0.06 deg along the belt (about 40 km) as if
    it had manoeuvred;
  * the pipeline (sensor/optical.py) detects, solves, measures and correlates
    against the catalogue's predictions, knowing nothing of the truth.

    python tools/optical_rehearsal.py data/raw results/optical
"""
import json
import math
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sensor.optical import ARCSEC, Camera, correlate, detect, render, solve_checked, star_field, streaks, subtract_stars  # noqa: E402
from sos.catalog import Catalog, load_snapshots, read_gp_csv, read_satcat  # noqa: E402
from sos.uncertainty import Model  # noqa: E402

from skyfield.api import load, wgs84  # noqa: E402
from skyfield.sgp4lib import EarthSatellite  # noqa: E402

raw, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
ts = load.timescale()
SITE = wgs84.latlon(43.66, -79.40, elevation_m=150)          # Toronto
cam = Camera()
model = Model(Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).resolve().parent.parent / "results" / "uncertainty.json")
rng = np.random.default_rng(42)

cat = Catalog(load_snapshots(raw))
sup = {}
for f in ("intelsat", "ses", "telesat", "eumetsat"):
    p = raw / f"sup-{f}.csv"
    if p.exists():
        for e in read_gp_csv(p.read_text()):
            if e.norad in cat.index and cat.elements[cat.index[e.norad]].intdes == e.intdes:
                sup[e.norad] = e
satcat = read_satcat((raw / "satcat.csv").read_text()) if (raw / "satcat.csv").exists() else {}
geo = [e for e in cat.elements if e.regime == "GEO"]
print("GEO objects", len(geo), "with operator truth", sum(1 for e in geo if e.norad in sup))


def radec(sat, t):
    a = (sat - SITE).at(t)
    ra, dec, _ = a.radec()
    alt, _, dist = a.altaz()
    return ra.radians, dec.radians, alt.degrees, dist.km


def magnitude(kind, rng_km, phase_deg):
    area = {"PAY": 20.0, "R/B": 8.0, "DEB": 0.8}.get(kind, 4.0)            # m^2
    ph = math.radians(phase_deg)
    f = (2 / (3 * math.pi ** 2)) * ((math.pi - ph) * math.cos(ph) + math.sin(ph))
    return -26.74 - 2.5 * math.log10(0.175 * area * f / (rng_km * 1000) ** 2)


# observing time: 23:00 local (04:00 UTC), late November
t0 = datetime(2026, 11, 25, 4, 0, tzinfo=timezone.utc)
canadian = ["ANIK F1R", "ANIK F2", "NIMIQ 6", "NIMIQ 4", "NIMIQ 5", "TELSTAR 19V", "ANIK F3", "NIMIQ 2"]
# plus a survey of the visible belt: the catalogued GEO object nearest every 4.5 degrees from 30 W to 133.5 W
from sos.catalog import ecef_to_geodetic, jd_array, teme_to_ecef  # noqa: E402
_jd, _fr = jd_array([t0])
_r, _, _ = cat.propagate(_jd, _fr, subset=[cat.index[e.norad] for e in geo])
_lat, _lon, _h = ecef_to_geodetic(teme_to_ecef(_r[:, 0, :], _jd, _fr))
survey = []
for lon in np.arange(-30.0, -134.0, -4.5):
    k = int(np.argmin(np.abs(_lon - lon)))
    if abs(_lon[k] - lon) < 2.0 and geo[k].name not in canadian and geo[k].name not in survey:
        survey.append(geo[k].name)
targets = canadian + survey
eph = load("de421.bsp") if Path("de421.bsp").exists() else None
fields, measures = [], []
injected = {"uct": [], "manoeuvre": []}
for fi, name in enumerate(targets):
    tgt = next((e for e in geo if e.name.startswith(name)), None)
    if tgt is None:
        continue
    t = ts.from_datetime(t0 + timedelta(minutes=4 * fi))
    tsat = EarthSatellite.from_satrec(tgt.sat, ts)
    ra0, dec0, alt, _ = radec(tsat, t)
    if alt < 15:
        continue
    point_ra = ra0 + math.radians(rng.normal(0, 0.12))       # mount pointing error
    point_dec = dec0 + math.radians(rng.normal(0, 0.12))
    roll = math.radians(rng.normal(0, 0.3))
    stars = star_field(point_ra, point_dec, 2.6, seed=1000 + fi)
    # the pipeline's catalogue has its own small errors relative to the sky that is rendered
    cat_err = np.radians(cam.catalogue_error_arcsec / 3600)
    stars_cat = (stars[0] + rng.normal(0, cat_err, len(stars[0])) / np.cos(stars[1]), stars[1] + rng.normal(0, cat_err, len(stars[1])), stars[2])
    sats, predicted = [], []
    for e in geo:
        s_gp = EarthSatellite.from_satrec(e.sat, ts)
        ra, dec, a, d = radec(s_gp, t)
        if math.degrees(math.acos(min(1, math.sin(dec) * math.sin(point_dec) + math.cos(dec) * math.cos(point_dec) * math.cos(ra - point_ra)))) > 2.4:
            continue
        predicted.append({"norad": e.norad, "name": e.name, "ra": ra, "dec": dec, "age_days": (t0 - e.epoch).total_seconds() / 86400})
        # truth: operator orbit if available, else catalogue + modelled error
        if e.norad in sup:
            tr = EarthSatellite.from_satrec(sup[e.norad].sat, ts)
            rt, dt_, _, dist = radec(tr, t)
            rt2, dt2, _, _ = radec(tr, ts.from_datetime(t.utc_datetime() + timedelta(seconds=1)))
            source = "operator"
        else:
            sig = model.sigma("GEO", (t0 - e.epoch).total_seconds() / 86400)
            off = rng.normal(0, 1, 3) * sig                          # km in R, I, C ~ angles at range
            rt = ra + off[1] / d * 1 / max(math.cos(dec), 0.2)
            dt_ = dec + off[2] / d
            ra2, dec2, _, _ = radec(s_gp, ts.from_datetime(t.utc_datetime() + timedelta(seconds=1)))
            rt2, dt2 = rt + (ra2 - ra), dt_ + (dec2 - dec)
            dist = d
            source = "modelled"
        kind = satcat.get(e.norad, {}).get("OBJECT_TYPE", "PAY")
        mag = magnitude(kind, dist, 25.0) + rng.normal(0, 0.3)
        # observation time t is mid-exposure: the streak starts half an exposure earlier
        h = cam.exposure_s / 2
        sats.append({"norad": e.norad, "name": e.name, "ra": rt - (rt2 - rt) * h, "dec": dt_ - (dt2 - dt_) * h, "mag": mag,
                     "rate_ra": rt2 - rt, "rate_dec": dt2 - dt_, "truth": source, "kind": kind, "ra_mid": rt, "dec_mid": dt_})
    # injected tests
    if sats:
        base = sats[0]
        # an uncatalogued object at a spot at least 0.1 deg from every satellite
        for off in ((0.7, -0.4), (-0.8, 0.5), (0.5, 0.9), (-0.6, -0.9), (1.1, 0.2)):
            u_ra, u_dec = point_ra + math.radians(off[0]) / math.cos(point_dec), point_dec + math.radians(off[1])
            if all(math.hypot((o["ra_mid"] - u_ra) * math.cos(u_dec), o["dec_mid"] - u_dec) > math.radians(0.1) for o in sats):
                break
        uct = {"norad": None, "name": "UNCATALOGUED", "ra": u_ra, "dec": u_dec,
               "mag": 13.0, "rate_ra": base["rate_ra"], "rate_dec": base["rate_dec"], "truth": "injected", "kind": "UCT"}
        uct["ra_mid"] = uct["ra"] + uct["rate_ra"] * cam.exposure_s / 2
        uct["dec_mid"] = uct["dec"] + uct["rate_dec"] * cam.exposure_s / 2
        sats.append(uct)
        injected["uct"].append(fi)
        def clear(s0, others, shift=0.0):
            ra_s = s0["ra_mid"] + shift
            return all(math.hypot((o["ra_mid"] - ra_s) * math.cos(o["dec_mid"]), o["dec_mid"] - s0["dec_mid"]) > math.radians(0.1)
                       for o in others if o is not s0)
        shift_try = math.radians(0.06)
        def central(x):
            return math.hypot((x["ra_mid"] - point_ra) * math.cos(point_dec), x["dec_mid"] - point_dec) < math.radians(1.2)
        cands = [x for x in sats[:-1] if central(x) and clear(x, sats, shift_try / max(math.cos(x["dec_mid"]), 0.2))]
        mv = min(cands or sats[:-1], key=lambda s: s["mag"])   # brightest satellite with clear sky after the shift
        shift = math.radians(0.06) / max(math.cos(mv["dec"]), 0.2)
        mv["ra"] += shift
        mv["ra_mid"] += shift
        mv["truth"] += "+manoeuvre"
        injected["manoeuvre"].append(mv["norad"])
    img, truth = render(cam, point_ra, point_dec, roll, stars, sats, seed=fi)
    dets, noise = detect(img)
    plate, quality = solve_checked(dets, *stars_cat, point_ra, point_dec, cam)
    if plate is None:
        fields.append({"target": name, "solved": False, "quality": quality, "satellites_in_frame": len(truth)})
        print("REJECTED", name, quality)
        continue
    exp_angle = math.atan2(truth[0]["dy"], truth[0]["dx"]) if truth else None
    clean, zp = subtract_stars(img, plate, *stars_cat, dets, cam)
    dets2, _ = detect(clean)
    # what the pipeline knows: the expected streak from the catalogue predictions in the field, the star positions
    exp_len = math.hypot(truth[0]["dx"], truth[0]["dy"]) if truth else None
    sx, sy = plate.to_px(stars_cat[0], stars_cat[1])
    found = streaks(dets2, exp_angle, expected_len_px=exp_len, star_xy=np.column_stack([sx, sy]), star_mag=stars_cat[2],
                    star_flux=10 ** (-0.4 * (stars_cat[2] - zp)) * cam.exposure_s)
    false_streaks = [d for d in found if all(math.hypot(d["x"] - t["x"], d["y"] - t["y"]) >= 4 for t in truth)]
    for d in false_streaks:
        k = int(np.argmin(np.hypot(sx - d["x"], sy - d["y"])))
        print("FALSE", name, {q: round(v, 2) for q, v in d.items()}, "nearest star px", round(float(np.hypot(sx[k] - d["x"], sy[k] - d["y"])), 1),
              "mag", round(float(stars_cat[2][k]), 2), "stars within 12 px", [(round(float(np.hypot(sx[q] - d["x"], sy[q] - d["y"])), 1), round(float(stars_cat[2][q]), 2)) for q in np.nonzero(np.hypot(sx - d["x"], sy - d["y"]) < 12)[0]], "peak", float(img[int(d["y"]) - 3:int(d["y"]) + 4, int(d["x"]) - 3:int(d["x"]) + 4].max()), "nearest sat px", round(min([math.hypot(d["x"] - t["x"], d["y"] - t["y"]) for t in truth] or [0]), 1))
    meas = [plate.to_radec(d["x"], d["y"]) for d in found]
    meas = [(float(a), float(b)) for a, b in meas]
    corr = correlate(meas, predicted, gate_arcsec=400)
    # score each truth satellite
    for tr in truth:
        best, bd = None, 1e9
        for k, d in enumerate(found):
            dd = math.hypot(d["x"] - tr["x"], d["y"] - tr["y"])
            if dd < bd:
                best, bd = k, dd
        detected = best is not None and bd < 4
        if not detected:
            near = sorted(dets2, key=lambda d: math.hypot(d["x"] - tr["x"], d["y"] - tr["y"]))[:2]
            print("MISS", name, tr["name"], round(tr["x"]), round(tr["y"]), "streak dx,dy", round(tr["dx"], 1), round(tr["dy"], 1),
                  [(round(math.hypot(d["x"] - tr["x"], d["y"] - tr["y"]), 1), round(d["elong"], 2), round(math.degrees(d["angle"]), 1), round(d["snr"])) for d in near],
                  "expected angle", round(math.degrees(exp_angle), 1) if exp_angle is not None else None)
        rec = {"field": name, "norad": tr["norad"], "name": tr["name"], "mag": round(tr["mag"], 2), "truth": tr["truth"],
               "kind": tr["kind"], "detected": bool(detected)}
        if detected:
            ra_m, dec_m = meas[best]
            err = math.degrees(math.acos(min(1, math.sin(dec_m) * math.sin(tr["dec_mid"]) + math.cos(dec_m) * math.cos(tr["dec_mid"]) * math.cos(ra_m - tr["ra_mid"])))) * 3600
            c = corr[best]
            rec.update({"astrometric_error_arcsec": round(err, 2), "correlated_to": c["norad"], "residual_arcsec": round(c["sep_arcsec"], 1)})
            if c["norad"] is not None:
                pred = next(p for p in predicted if p["norad"] == c["norad"])
                sig_km = float(np.linalg.norm(model.sigma("GEO", pred["age_days"])[1:]))
                exp_arcsec = sig_km / 37000 / ARCSEC
                rec["expected_arcsec_1sigma"] = round(exp_arcsec, 1)
                rec["flag"] = "residual over 4 sigma: manoeuvre or stale elements" if c["sep_arcsec"] > 4 * exp_arcsec else ""
            else:
                rec["flag"] = "uncorrelated track"
        measures.append(rec)
    fields.append({"target": name, "solved": True, "stars_matched": plate.n_stars, "plate_rms_arcsec": round(plate.rms_arcsec, 2),
                   "detections": len(dets), "streaks": len(found), "satellites_in_frame": len(truth),
                   "false_streaks": len(false_streaks),
                   "false_uct_alarms": sum(1 for d in false_streaks if corr[found.index(d)]["norad"] is None),
                   "pointing_error_arcmin": round(math.degrees(math.acos(min(1, math.sin(plate.dec0) * math.sin(point_dec) + math.cos(plate.dec0) * math.cos(point_dec) * math.cos(plate.ra0 - point_ra)))) * 60, 2)})
    if fi == 0:
        from PIL import Image                     # a stretched crop of the first frame, for the docs
        crop = img[1000:2000, 1000:2000].astype(float)
        lo, hi = np.percentile(crop, [5, 99.7])
        Image.fromarray((np.clip((crop - lo) / (hi - lo), 0, 1) ** 0.5 * 255).astype("uint8")).save(out / "frame0.jpg", quality=88)

cat_meas = [m for m in measures if m["kind"] != "UCT"]
det = [m for m in cat_meas if m["detected"]]
summary = {
    "fields": len(fields), "solved": sum(f.get("solved", False) for f in fields),
    "rejected_by_quality_gate": [f"{f['target']}: {f['quality']}" for f in fields if not f.get("solved")],
    "plate_rms_arcsec_median": float(np.median([f["plate_rms_arcsec"] for f in fields if f.get("solved")])),
    "satellites_in_frames": len(cat_meas), "detected": len(det),
    "completeness_by_mag": {f"{lo}-{lo + 1}": [sum(1 for m in cat_meas if lo <= m["mag"] < lo + 1 and m["detected"]), sum(1 for m in cat_meas if lo <= m["mag"] < lo + 1)] for lo in range(9, 18)},
    "astrometric_error_arcsec_median": float(np.median([m["astrometric_error_arcsec"] for m in det])) if det else None,
    "astrometric_error_arcsec_p95": float(np.percentile([m["astrometric_error_arcsec"] for m in det], 95)) if det else None,
    "correct_correlation": sum(1 for m in det if m["correlated_to"] == m["norad"]), "wrong_correlation": sum(1 for m in det if m["correlated_to"] not in (None, m["norad"])),
    "uncorrelated_catalogued": sum(1 for m in det if m["correlated_to"] is None),
    "uct_injected": len(injected["uct"]), "uct_detected_and_flagged": sum(1 for m in measures if m["kind"] == "UCT" and m.get("flag") == "uncorrelated track"),
    "manoeuvre_injected": sum(1 for m in cat_meas if "manoeuvre" in m["truth"]),
    "manoeuvre_flagged": sum(1 for m in det if "manoeuvre" in m["truth"] and m.get("flag", "").startswith("residual")),
    "false_manoeuvre_flags": sum(1 for m in det if "manoeuvre" not in m["truth"] and m.get("flag", "").startswith("residual")),
    "false_streaks": sum(f.get("false_streaks", 0) for f in fields),
    "false_uct_alarms": sum(f.get("false_uct_alarms", 0) for f in fields),
    "operator_truth_objects": sum(1 for m in cat_meas if m["truth"].startswith("operator")),
}
(out / "optical.json").write_text(json.dumps({"summary": summary, "fields": fields, "measurements": measures}, indent=1))
print(json.dumps(summary, indent=1))
