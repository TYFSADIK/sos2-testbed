#!/usr/bin/env python3
"""Measure how wrong public element sets are, and fit an uncertainty model.

Truth sources (all public):
  * CelesTrak SupGP: orbits fitted to the operators' own ephemerides (Starlink,
    GPS, Intelsat, SES, Telesat, EUMETSAT) and to laser-ranging predictions (CPF).
    Used only within 12 h after their own epoch, where they track the operator
    data to their published fit RMS.
  * NASA's ISS ephemeris (CCSDS OEM, EME2000), 15 days ahead.

For every object present in both, the catalogue (GP) element set is propagated
to the truth times and the difference is resolved into radial, in-track and
cross-track. Error is then fitted against element-set age per orbit regime.

    python tools/calibrate.py data/raw results/
"""
import glob
import json
import math
import sys
from collections import defaultdict
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sos.catalog import Catalog, jd_array, load_snapshots, read_gp_csv, ric_frame  # noqa: E402

raw, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)

cat = Catalog(load_snapshots(raw))
print("catalogue objects", len(cat))

samples = defaultdict(list)          # regime -> [(age_days, dR, dI, dC, source)]
per_source = defaultdict(lambda: {"objects": 0, "samples": 0})
mismatched = defaultdict(int)
for f in sorted(glob.glob(str(raw / "sup-*.csv"))):
    src = Path(f).stem.replace("sup-", "")
    if src == "iss":
        continue                     # the ISS is checked against NASA's own ephemeris below
    try:
        sup = read_gp_csv(open(f).read())
    except ValueError:
        continue
    for s in sup:
        i = cat.index.get(s.norad)
        if i is None:
            continue
        g = cat.elements[i]
        if g.fields.get("EPOCH") == s.fields.get("EPOCH"):
            continue
        if g.intdes != s.intdes:
            mismatched[src] += 1         # temporary catalogue numbers are reused across sources
            continue
        times = [s.epoch + timedelta(hours=h) for h in np.arange(0, 12.01, 1.0)]
        jd, fr = jd_array(times)
        e1, rt, vt = s.sat.sgp4_array(jd, fr)
        e2, rg, _ = g.sat.sgp4_array(jd, fr)
        ok = (e1 == 0) & (e2 == 0)
        if not ok.any():
            continue
        rh, ih, ch = ric_frame(rt[ok], vt[ok])
        d = rg[ok] - rt[ok]
        ages = np.array([(t - g.epoch).total_seconds() / 86400 for t, k in zip(times, ok) if k])
        for a, dr, di, dc in zip(ages, (d * rh).sum(1), (d * ih).sum(1), (d * ch).sum(1)):
            samples[g.regime].append((float(a), float(dr), float(di), float(dc), src))
        per_source[src]["objects"] += 1
        per_source[src]["samples"] += int(ok.sum())

# ---- ISS against NASA's ephemeris
iss_curve = []
oem = raw / "iss-oem.txt"
if oem.exists() and 25544 in cat.index:
    from skyfield.api import load as sf_load
    from skyfield.sgp4lib import EarthSatellite
    ts = sf_load.timescale()
    g = cat.elements[cat.index[25544]]
    sat = EarthSatellite.from_satrec(g.sat, ts)
    rows = []
    for line in oem.read_text().splitlines():
        p = line.split()
        if len(p) == 7 and p[0][:4].isdigit() and "T" in p[0]:
            rows.append((datetime.fromisoformat(p[0]).replace(tzinfo=timezone.utc), *map(float, p[1:])))
    rows = rows[::8]                 # every 32 minutes
    t = ts.from_datetimes([r[0] for r in rows])
    pos = sat.at(t)
    rg = pos.position.km.T
    rt = np.array([r[1:4] for r in rows])
    vt = np.array([r[4:7] for r in rows])
    rh, ih, ch = ric_frame(rt, vt)
    d = rg - rt
    for (tt, *_), dr, di, dc in zip(rows, (d * rh).sum(1), (d * ih).sum(1), (d * ch).sum(1)):
        age = (tt - g.epoch).total_seconds() / 86400
        iss_curve.append({"age_days": round(age, 3), "R": round(float(dr), 3), "I": round(float(di), 3), "C": round(float(dc), 3),
                          "total": round(float(math.sqrt(dr * dr + di * di + dc * dc)), 3)})
        samples["ISS"].append((age, float(dr), float(di), float(dc), "nasa-oem"))


# ---- fit sigma(age) = a + b * age per regime and axis, from RMS in age bins
def rms(x):
    """Robust 1-sigma: 1.4826 x median absolute error (manoeuvring objects are a separate population,
    caught by manoeuvre detection; they would otherwise dominate a root-mean-square)."""
    return float(1.4826 * np.median(np.abs(x))) if len(x) else None


fit, table = {}, {}
for reg, s in samples.items():
    a = np.array([x[0] for x in s])
    d = np.array([x[1:4] for x in s])
    bins = [0, 0.25, 0.5, 1, 2, 3, 5, 8, 15]
    rows = []
    for lo, hi in zip(bins, bins[1:]):
        m = (a >= lo) & (a < hi)
        if m.sum() >= 5:
            rows.append({"age_from": lo, "age_to": hi, "n": int(m.sum()), "R": rms(d[m, 0]), "I": rms(d[m, 1]), "C": rms(d[m, 2]),
                         "total_median": float(np.median(np.linalg.norm(d[m], axis=1))),
                         "total_p95": float(np.percentile(np.linalg.norm(d[m], axis=1), 95)),
                         "share_over_50km": float(np.mean(np.linalg.norm(d[m], axis=1) > 50))})
    f = {}
    for k, axis in (("R", 0), ("I", 1), ("C", 2)):
        if len(rows) >= 2:
            x = np.array([(r["age_from"] + r["age_to"]) / 2 for r in rows])
            y = np.array([r[k] for r in rows])
            w = np.sqrt([r["n"] for r in rows])
            if k == "I" and len(rows) >= 4:
                # along-track error from a small mismatch in predicted decay grows faster than linearly
                c2, b, a0 = np.polyfit(x, y, 2, w=w)
                f[k] = [max(float(a0), 0.01), max(float(b), 0.0), max(float(c2), 0.0)]
            else:
                b, a0 = np.polyfit(x, y, 1, w=w)
                f[k] = [max(float(a0), 0.01), max(float(b), 0.0)]
        elif rows:
            f[k] = [rows[0][k], 0.0]
    if f:
        fit[reg] = f
    table[reg] = {"samples": len(s), "objects_by_source": dict(per_source), "bins": rows,
                  "median_total_km": float(np.median(np.linalg.norm(d, axis=1))),
                  "p95_total_km": float(np.percentile(np.linalg.norm(d, axis=1), 95)),
                  "median_age_days": float(np.median(a))}

model = {"source": "tools/calibrate.py: GP vs operator SupGP and NASA ISS ephemeris", "created": datetime.now(timezone.utc).isoformat(),
         "fit": {k: v for k, v in fit.items() if k != "ISS"}, "fit_iss": fit.get("ISS")}
(out / "uncertainty.json").write_text(json.dumps(model, indent=1))
(out / "calibration.json").write_text(json.dumps({"regimes": table, "sources": dict(per_source), "iss_vs_nasa": iss_curve}, indent=1))
for reg, t in table.items():
    print(f"{reg:4} samples {t['samples']:6}  median {t['median_total_km']:.3f} km  p95 {t['p95_total_km']:.3f} km  median age {t['median_age_days']:.2f} d  fit {fit.get(reg)}")
if iss_curve:
    for c in iss_curve[:: max(1, len(iss_curve) // 12)]:
        print("ISS age", c["age_days"], "d: total", c["total"], "km (in-track", c["I"], ")")
print(json.dumps(dict(per_source)), "skipped (designator mismatch):", dict(mismatched))
