"""How wrong is a public element set, and how fast does it get worse?

Measured, not assumed: tools/calibrate.py compares catalogue (GP) predictions
with operators' own orbit data (CelesTrak SupGP: Starlink, GPS, Intelsat, SES,
Telesat, EUMETSAT, laser-ranging targets) and NASA's ISS ephemeris, and fits
position error against element-set age per orbit regime, in the radial (R),
in-track (I) and cross-track (C) directions. This module turns that fit into a
covariance for a given object, regime and age, and into collision probability.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from .catalog import ric_frame

DEFAULT = Path(__file__).resolve().parent.parent / "results" / "uncertainty.json"
# used only if no calibration exists yet (rough literature values, km at age 0 and growth per day)
FALLBACK = {
    "LEO": {"R": [0.1, 0.05], "I": [0.4, 1.0], "C": [0.1, 0.05]},
    "MEO": {"R": [0.5, 0.1], "I": [1.0, 0.5], "C": [0.5, 0.1]},
    "GEO": {"R": [1.0, 0.2], "I": [3.0, 1.0], "C": [1.0, 0.2]},
    "HEO": {"R": [2.0, 1.0], "I": [5.0, 5.0], "C": [2.0, 1.0]},
}


class Model:
    def __init__(self, path: Path | str = DEFAULT):
        self.fit = FALLBACK
        self.source = "fallback"
        p = Path(path)
        if p.exists():
            d = json.loads(p.read_text())
            self.fit = {k: v for k, v in d["fit"].items()}
            self.source = d.get("source", str(p))

    def sigma(self, regime: str, age_days: float) -> np.ndarray:
        """1-sigma position uncertainty (km) in R, I, C."""
        f = self.fit.get(regime) or self.fit.get({"GTO": "HEO"}.get(regime, "HEO")) or FALLBACK["HEO"]
        a = abs(age_days)
        return np.array([f[k][0] + f[k][1] * a + (f[k][2] * a * a if len(f[k]) > 2 else 0.0) for k in ("R", "I", "C")])

    @staticmethod
    def sigma_operator(rms_km, age_days) -> np.ndarray:
        """Operator-derived orbits (SupGP): CelesTrak publishes each fit's RMS against the operator's
        ephemeris. Assumption, not measured here: that RMS at epoch, growing 0.3 km/day in-track."""
        r = max(float(rms_km or 0.2), 0.05)
        a = abs(age_days)
        return np.array([r / 2, r + 0.3 * a, r / 2])

    def element_sigma(self, element, age_days) -> np.ndarray:
        """Operator orbits get their fit RMS only within a day of their epoch; after that their real
        accuracy is not known from the data here, so the measured public-data model applies (taking the
        larger of the two). Without this, operator orbits propagated days ahead produced tiny
        uncertainties and implausibly high collision probabilities."""
        gp = self.sigma(element.regime, age_days)
        if getattr(element, "source", "gp") == "supgp" and abs(age_days) <= 1.0:
            return self.sigma_operator(element.fields.get("RMS"), age_days)
        if getattr(element, "source", "gp") == "supgp":
            return np.maximum(gp, self.sigma_operator(element.fields.get("RMS"), age_days))
        return gp

    def covariance_eci(self, r, v, regime, age_days, element=None) -> np.ndarray:
        rh, ih, ch = ric_frame(np.asarray(r), np.asarray(v))
        m = np.stack([rh, ih, ch], axis=1)          # columns R, I, C
        s = self.element_sigma(element, age_days) if element is not None else self.sigma(regime, age_days)
        return m @ np.diag(s ** 2) @ m.T


def encounter_plane(ra, va, rb, vb):
    """Miss vector and 2x3 projection onto the plane normal to the relative velocity."""
    dr = np.asarray(rb) - np.asarray(ra)
    dv = np.asarray(vb) - np.asarray(va)
    x = dv / np.linalg.norm(dv)
    y = np.cross(dr, x)
    if np.linalg.norm(y) < 1e-12:
        y = np.cross(x, [0, 0, 1.0])
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    P = np.stack([y, z])                          # 2 x 3
    return P @ dr, P


def pc_2d(miss2: np.ndarray, cov2: np.ndarray, hbr_km: float, n: int = 120) -> float:
    """Probability that the miss vector falls inside the hard-body circle (Foster's method),
    integrated numerically on a polar grid."""
    try:
        inv = np.linalg.inv(cov2)
        det = np.linalg.det(cov2)
    except np.linalg.LinAlgError:
        return 0.0
    if det <= 0:
        return 0.0
    rr = (np.arange(n) + 0.5) / n * hbr_km
    th = (np.arange(n) + 0.5) / n * 2 * np.pi
    R, T = np.meshgrid(rr, th)
    pts = np.stack([miss2[0] + R * np.cos(T), miss2[1] + R * np.sin(T)], axis=-1)
    q = np.einsum("...i,ij,...j->...", pts, inv, pts)
    dens = np.exp(-0.5 * q) / (2 * np.pi * math.sqrt(det))
    return float((dens * R).sum() * (hbr_km / n) * (2 * np.pi / n))


def pc_max(miss_km: float, hbr_km: float) -> float:
    """Worst-case probability over the size of a circular covariance (Alfano): R^2 / (e d^2)."""
    if miss_km <= hbr_km:
        return 1.0
    return min(1.0, hbr_km ** 2 / (math.e * miss_km ** 2))


HBR = {"PAYLOAD": 0.005, "ROCKET BODY": 0.008, "DEBRIS": 0.002, "UNKNOWN": 0.003}   # km, combined with the other object


def hard_body_radius(kind_a: str, kind_b: str, rcs_a: str = "", rcs_b: str = "") -> float:
    size = {"SMALL": 0.0005, "MEDIUM": 0.0015, "LARGE": 0.004}
    ra = size.get(rcs_a.upper(), HBR.get(kind_a, 0.003)) if rcs_a else HBR.get(kind_a, 0.003)
    rb = size.get(rcs_b.upper(), HBR.get(kind_b, 0.003)) if rcs_b else HBR.get(kind_b, 0.003)
    return ra + rb
