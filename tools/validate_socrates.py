#!/usr/bin/env python3
"""Check the conjunction screener against CelesTrak SOCRATES, on identical inputs.

SOCRATES (STK/CAT + SGP4) lists every conjunction under 5 km between active
payloads and the whole public catalogue over 7 days. Sorted by time of
closest approach, its first entries are a complete time slice. For each pair
CelesTrak publishes the exact element sets it used.

  A. Pairwise: for every SOCRATES conjunction, refine TCA and miss distance from
     the same element sets; compare.
  B. Screening: run the full screener over the same time slice on every object
     that appears in it; does it find every SOCRATES conjunction, and anything
     SOCRATES did not list?

Docked or co-located objects (relative speed under 10 m/s, such as ISS modules
and visiting spacecraft) are separated out on both sides.

    python tools/validate_socrates.py data/raw results/socrates
"""
import glob
import json
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from sos.catalog import Catalog, jd_array, read_gp_csv, read_satcat  # noqa: E402
from sos.conjunction import refine, screen  # noqa: E402
from sos.socrates import parse  # noqa: E402

raw, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
soc = parse((raw / "socrates-tca.html").read_text())
items = soc["conjunctions"]
print("SOCRATES:", soc["meta"], len(items), "conjunctions listed")

els = []
for f in sorted(glob.glob(str(raw / "socrates-gp" / "*.csv"))):
    try:
        els += read_gp_csv(open(f).read())
    except ValueError:
        pass
cat = Catalog(els)
have = set(cat.index)
items = [c for c in items if c["norad1"] in have and c["norad2"] in have]
print("conjunctions with both element sets available:", len(items), "objects:", len(cat))
EPOCH0 = datetime(2000, 1, 1, 12, tzinfo=timezone.utc)


def jd_of(dt):
    jd, fr = jd_array([dt])
    return jd[0], fr[0]


def to_dt(jd):
    return EPOCH0 + timedelta(days=jd - 2451545.0)


# ---- A: pairwise
pair = []
for c in items:
    tca = datetime.fromisoformat(c["tca"])
    jd, fr = jd_of(tca)
    i, j = cat.index[c["norad1"]], cat.index[c["norad2"]]
    jd_tca, miss, vrel, *_ = refine(cat, i, j, jd, fr, 60 / 86400)
    pair.append({"a": c["norad1"], "b": c["norad2"], "dt_s": (to_dt(jd_tca) - tca).total_seconds(), "socrates_km": c["range_km"],
                 "mine_km": miss, "d_miss_m": (miss - c["range_km"]) * 1000, "socrates_rel": c["rel_speed_kms"], "mine_rel": vrel,
                 "colocated": c["rel_speed_kms"] < 0.01})
real = [p for p in pair if not p["colocated"]]


def stats(xs):
    a = np.abs(np.array(xs))
    return {"median": round(float(np.median(a)), 4), "p95": round(float(np.percentile(a, 95)), 4), "max": round(float(a.max()), 4)} if len(a) else None


A = {"compared": len(real), "colocated_skipped": len(pair) - len(real),
     "abs_dtca_s": stats([p["dt_s"] for p in real]), "abs_dmiss_m": stats([p["d_miss_m"] for p in real]),
     "abs_drel_speed_ms": stats([(p["mine_rel"] - p["socrates_rel"]) * 1000 for p in real])}
print("A pairwise:", json.dumps(A))

# ---- B: screening over the same slice
tcas = sorted(datetime.fromisoformat(c["tca"]) for c in items)
start = tcas[0] - timedelta(seconds=30)
stop = tcas[-1] - timedelta(seconds=60)           # the list is cut at its last entry: stay inside it
satcat = read_satcat((raw / "satcat.csv").read_text()) if (raw / "satcat.csv").exists() else {}
# SOCRATES primaries are active payloads: use its own status flags where given, else SATCAT
status = {}
for c in soc["conjunctions"]:
    status[c["norad1"]] = c["status1"]
    status[c["norad2"]] = c["status2"]
ACTIVE = {"+", "P", "B", "S", "X"}
prim = np.array([status.get(int(n), satcat.get(int(n), {}).get("OPS_STATUS_CODE", "")) in ACTIVE
                 and satcat.get(int(n), {}).get("OBJECT_TYPE", "PAY") == "PAY" for n in cat.norad])
jd, fr = jd_of(start)
days = (stop - start).total_seconds() / 86400
t0 = time.time()
mine, colo, st = screen(cat, jd + fr, days, 5.0, 20.0, primaries=prim)
runtime = time.time() - t0
# score only the interior: one minute in from SOCRATES's own window start and from the end of the slice
win_start = datetime.strptime(soc["meta"]["start"], "%Y %b %d %H:%M:%S").replace(tzinfo=timezone.utc)
inner = (max(start, win_start) + timedelta(seconds=60), stop - timedelta(seconds=60))
soc_real = [c for c in items if c["rel_speed_kms"] >= 0.01 and inner[0] <= datetime.fromisoformat(c["tca"]) <= inner[1]]
mine = [m for m in mine if inner[0] <= to_dt(m.tca_jd) <= inner[1]]
used = set()
matched = []
for c in soc_real:
    key = (min(c["norad1"], c["norad2"]), max(c["norad1"], c["norad2"]))
    t = datetime.fromisoformat(c["tca"])
    best = None
    for k, m in enumerate(mine):
        if m.key() == key and k not in used and abs((to_dt(m.tca_jd) - t).total_seconds()) < 120:
            best = k
            break
    if best is not None:
        used.add(best)
        matched.append(c)
missed = [c for c in soc_real if c not in matched]
extra = [m for k, m in enumerate(mine) if k not in used]
B = {"slice_utc": [start.isoformat(), stop.isoformat()], "scored_utc": [inner[0].isoformat(), inner[1].isoformat()], "slice_minutes": round(days * 1440, 1), "objects": len(cat),
     "primaries": int(prim.sum()), "socrates_conjunctions": len(soc_real), "found": len(matched), "missed": len(missed),
     "extra_not_in_socrates": len(extra), "runtime_s": round(runtime, 1), "screen_stats": st,
     "missed_detail": [{k: c[k] for k in ("norad1", "name1", "norad2", "name2", "tca", "range_km", "rel_speed_kms")} for c in missed[:20]],
     "extra_detail": [{"a": m.a, "b": m.b, "tca": to_dt(m.tca_jd).isoformat(), "miss_km": round(m.miss_km, 3),
                       "rel_kms": round(m.rel_speed_kms, 3), "a_primary": bool(prim[cat.index[m.a]]), "b_primary": bool(prim[cat.index[m.b]])}
                      for m in extra[:20]]}
print("B screening:", json.dumps({k: v for k, v in B.items() if not k.endswith("detail")}))
(out / "socrates-validation.json").write_text(json.dumps({"socrates_meta": soc["meta"], "pairwise": A, "screening": B,
                                                          "pairs": pair}, indent=1, default=str))
