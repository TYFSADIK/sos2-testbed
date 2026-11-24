"""Deep space: characterise the geostationary belt and watch Canadian assets.

For every GEO object: where it sits (longitude), how fast it drifts, how
inclined it is, and so what it is doing:
  station-kept           holds its slot: drift under 0.05 deg/day, inclination under 0.2 deg
  inclined station-kept  holds longitude but lets inclination grow (saves fuel late in life)
  drifting               moving along the belt: dead, being relocated, or in a graveyard orbit

For each Canadian GEO satellite (SATCAT owner CA), list the objects that come
closest over the next day, by actual propagated distance.
"""
from __future__ import annotations

import numpy as np

from .catalog import Catalog, ecef_to_geodetic, teme_to_ecef

N_GEO = 1.00273791          # sidereal rev/day


def drift_deg_per_day(mean_motion: float) -> float:
    return (mean_motion - N_GEO) * 360.0


def classify(drift: float, inc: float) -> str:
    if abs(drift) < 0.05:
        return "station-kept" if inc < 0.2 else "inclined station-kept"
    return "drifting"


def belt(cat: Catalog, jd: float, fr: float, satcat: dict) -> list[dict]:
    idx = [i for i, e in enumerate(cat.elements) if e.regime == "GEO"]
    if not idx:
        return []
    r, _, err = cat.propagate(np.array([jd]), np.array([fr]), subset=idx)
    lat, lon, h = ecef_to_geodetic(teme_to_ecef(r[:, 0, :], np.array([jd]), np.array([fr])))
    out = []
    for k, i in enumerate(idx):
        e = cat.elements[i]
        sc = satcat.get(e.norad, {})
        d = drift_deg_per_day(e.mean_motion)
        out.append({"norad": e.norad, "name": e.name, "owner": sc.get("OWNER", ""), "type": sc.get("OBJECT_TYPE", ""),
                    "status": sc.get("OPS_STATUS_CODE", ""), "lon": round(float(lon[k]), 3), "lat": round(float(lat[k]), 3),
                    "alt_km": round(float(h[k]), 1), "drift_deg_day": round(d, 4), "inc": round(e.inc, 3),
                    "state": classify(d, e.inc)})
    return out


def canadian_neighbours(cat: Catalog, jd: float, satcat: dict, hours: float = 24.0, step_min: float = 10.0, top: int = 5):
    geo = [i for i, e in enumerate(cat.elements) if e.regime == "GEO"]
    ca = [i for i in geo if satcat.get(cat.elements[i].norad, {}).get("OWNER") == "CA"]
    if not ca:
        return []
    fr = np.arange(0, hours * 60 + 0.1, step_min) / 1440.0
    r, _, _ = cat.propagate(np.full(len(fr), jd), fr, subset=geo)
    pos = {g: r[k] for k, g in enumerate(geo)}
    out = []
    for i in ca:
        d = {g: float(np.min(np.linalg.norm(pos[g] - pos[i], axis=1))) for g in geo if g != i}
        near = sorted(d.items(), key=lambda kv: kv[1])[:top]
        e = cat.elements[i]
        out.append({"norad": e.norad, "name": e.name, "neighbours": [
            {"norad": cat.elements[g].norad, "name": cat.elements[g].name,
             "owner": satcat.get(cat.elements[g].norad, {}).get("OWNER", ""),
             "type": satcat.get(cat.elements[g].norad, {}).get("OBJECT_TYPE", ""), "min_km_24h": round(dist, 1)}
            for g, dist in near]})
    return out
