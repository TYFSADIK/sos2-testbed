"""The orbital catalogue: public GP element sets, SATCAT metadata, propagation.

GP ("general perturbations", the data behind TLEs) element sets come from
CelesTrak in OMM CSV form. Each is turned into an SGP4 satellite record;
SGP4 is the only propagator these elements are valid with. Positions are in
the TEME frame (what SGP4 produces) unless converted.
"""
from __future__ import annotations

import csv
import io
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone

import numpy as np
from sgp4 import omm
from sgp4.api import SatrecArray, Satrec, jday

MU = 398600.4418            # km3/s2
RE = 6378.137               # km
NOT_UPDATED = "has not updated since your last successful"


@dataclass
class Element:
    norad: int
    name: str
    intdes: str
    epoch: datetime
    fields: dict
    sat: Satrec = field(repr=False, default=None)
    source: str = "gp"             # gp | supgp (operator-derived)

    @property
    def mean_motion(self):          # rev/day
        return float(self.fields["MEAN_MOTION"])

    @property
    def ecc(self):
        return float(self.fields["ECCENTRICITY"])

    @property
    def inc(self):
        return float(self.fields["INCLINATION"])

    @property
    def sma(self):                  # km
        n = self.mean_motion * 2 * math.pi / 86400
        return (MU / n ** 2) ** (1 / 3)

    @property
    def perigee(self):
        return self.sma * (1 - self.ecc) - RE

    @property
    def apogee(self):
        return self.sma * (1 + self.ecc) - RE

    @property
    def period_min(self):
        return 1440.0 / self.mean_motion

    @property
    def regime(self):
        return regime(self.perigee, self.apogee, self.ecc, self.period_min, self.inc)


def regime(perigee, apogee, ecc, period_min, inc):
    """LEO / MEO / GEO / HEO / GTO, by the usual altitude and period bands."""
    if 1300 < period_min < 1580 and ecc < 0.03 and inc < 25:
        return "GEO"
    if ecc > 0.25:
        return "GTO" if perigee < 2000 and 30000 < apogee < 45000 else "HEO"
    if apogee < 2000:
        return "LEO"
    if perigee > 2000 and apogee < 35000:
        return "MEO"
    return "HEO"


def parse_epoch(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "")).replace(tzinfo=timezone.utc)


def read_gp_csv(text: str) -> list[Element]:
    if NOT_UPDATED in text[:300]:
        raise ValueError("CelesTrak: group not updated since the last download")
    out = []
    for f in omm.parse_csv(io.StringIO(text)):
        s = Satrec()
        omm.initialize(s, f)
        out.append(Element(int(f["NORAD_CAT_ID"]), f["OBJECT_NAME"].strip(), f["OBJECT_ID"].strip(),
                           parse_epoch(f["EPOCH"]), dict(f), s))
    return out


def read_satcat(text: str) -> dict[int, dict]:
    out = {}
    for r in csv.DictReader(io.StringIO(text)):
        try:
            out[int(r["NORAD_CAT_ID"])] = r
        except (KeyError, ValueError):
            continue
    return out


def jd_array(times: list[datetime] | np.ndarray):
    jd, fr = zip(*(jday(t.year, t.month, t.day, t.hour, t.minute, t.second + t.microsecond / 1e6) for t in times))
    return np.array(jd), np.array(fr)


class Catalog:
    """A set of element sets propagated together (vectorised SGP4)."""

    def __init__(self, elements: list[Element]):
        # one element set per object: the newest
        best: dict[int, Element] = {}
        for e in elements:
            if e.norad not in best or e.epoch > best[e.norad].epoch:
                best[e.norad] = e
        self.elements = sorted(best.values(), key=lambda e: e.norad)
        self.norad = np.array([e.norad for e in self.elements])
        self.index = {n: i for i, n in enumerate(self.norad)}
        self._array = SatrecArray([e.sat for e in self.elements]) if self.elements else None

    def __len__(self):
        return len(self.elements)

    def propagate(self, jd, fr, subset=None):
        """TEME position (km) and velocity (km/s), shape (n_objects, n_times, 3); error codes."""
        arr = self._array if subset is None else SatrecArray([self.elements[i].sat for i in subset])
        e, r, v = arr.sgp4(np.atleast_1d(jd), np.atleast_1d(fr))
        return r, v, e

    def epochs(self) -> np.ndarray:
        return np.array([e.epoch.timestamp() for e in self.elements])


def teme_to_ecef(r_teme: np.ndarray, jd_ut1: np.ndarray, fr: np.ndarray) -> np.ndarray:
    """Rotate TEME to an Earth-fixed frame by Greenwich mean sidereal time (polar motion
    ignored: under 15 m, far below GP accuracy)."""
    t = ((jd_ut1 - 2451545.0) + fr) / 36525.0
    gmst = (67310.54841 + (876600.0 * 3600 + 8640184.812866) * t + 0.093104 * t ** 2 - 6.2e-6 * t ** 3)
    gmst = np.radians((gmst % 86400) / 240.0)
    c, s = np.cos(gmst), np.sin(gmst)
    x, y, z = r_teme[..., 0], r_teme[..., 1], r_teme[..., 2]
    return np.stack([c * x + s * y, -s * x + c * y, z], axis=-1)


def ecef_to_geodetic(r: np.ndarray):
    """WGS84 latitude, longitude (degrees) and height (km), vectorised."""
    a, f = 6378.137, 1 / 298.257223563
    e2 = f * (2 - f)
    x, y, z = r[..., 0], r[..., 1], r[..., 2]
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)
    lat = np.arctan2(z, p * (1 - e2))
    for _ in range(5):
        n = a / np.sqrt(1 - e2 * np.sin(lat) ** 2)
        h = p / np.cos(lat) - n
        lat = np.arctan2(z, p * (1 - e2 * n / (n + h)))
    return np.degrees(lat), np.degrees(lon), h


def ric_frame(r: np.ndarray, v: np.ndarray):
    """Unit vectors radial, in-track, cross-track for each state (..., 3)."""
    rh = r / np.linalg.norm(r, axis=-1, keepdims=True)
    c = np.cross(r, v)
    ch = c / np.linalg.norm(c, axis=-1, keepdims=True)
    ih = np.cross(ch, rh)
    return rh, ih, ch


def load_snapshots(raw) -> list[Element]:
    """Every element set from every snapshot under raw/snapshots/<time>/gp-*.csv (skipping
    CelesTrak's "not updated" notices). A Catalog built from these keeps the newest per object."""
    import glob
    from pathlib import Path
    out = []
    for f in sorted(glob.glob(str(Path(raw) / "snapshots" / "*" / "gp-*.csv"))):
        try:
            out += read_gp_csv(open(f).read())
        except ValueError:
            pass
    return out
