"""Space situational awareness API and 3D viewer.

/api/positions        every catalogued object now (or at ?t=), compact arrays
/api/object/<norad>   metadata, freshness, uncertainty now, element history, one orbit
/api/conjunctions     the latest screening run, by collision probability
/api/conjunction/<a>/<b>  both tracks around the time of closest approach
/api/geo              the GEO belt now, and the neighbours of Canadian satellites
/api/search?q=        by name or number
/api/stats            catalogue, freshness and ingest health
"""
from __future__ import annotations


import os
import sys
import threading
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sos.catalog import ecef_to_geodetic, jd_array, teme_to_ecef          # noqa: E402
from sos.geo import canadian_neighbours                                    # noqa: E402
from sos.jobs import load_catalog, satcat_map                              # noqa: E402
from sos.store import connect                                              # noqa: E402
from sos.uncertainty import Model                                          # noqa: E402

DB = os.environ.get("DATABASE_URL", "postgresql://sos:sos@postgis/sos")
STATIC = Path(__file__).resolve().parent / "static"
conn = connect(DB)
model = Model(os.environ.get("UNCERTAINTY", "/app/results/uncertainty.json"))
state = {"cat": None, "satcat": {}, "loaded": 0.0, "lock": threading.Lock(), "geo_nb": None}
REG = {"LEO": 0, "MEO": 1, "GEO": 2, "HEO": 3, "GTO": 4}
TYP = {"PAY": 0, "R/B": 1, "DEB": 2, "UNK": 3}


def catalog():
    with state["lock"]:
        stale = time.time() - state["loaded"] > (300 if state["cat"] is not None and len(state["cat"]) else 10)
        if state["cat"] is None or stale:
            state["cat"] = load_catalog(conn)
            state["satcat"] = satcat_map(conn)
            state["loaded"] = time.time()
            state["regime"] = np.array([REG.get(e.regime, 3) for e in state["cat"].elements], dtype=np.int8)
            state["type"] = np.array([TYP.get(state["satcat"].get(e.norad, {}).get("OBJECT_TYPE", "UNK"), 3) for e in state["cat"].elements], dtype=np.int8)
            state["epochs"] = state["cat"].epochs()
        return state["cat"]


def when(t: str | None) -> datetime:
    return datetime.fromisoformat(t.replace("Z", "+00:00")) if t else datetime.now(timezone.utc)


app = FastAPI(title="Surveillance of Space testbed")


@app.get("/api/health")
def health():
    return {"ok": True, "objects": len(catalog())}


@app.get("/api/positions")
def positions(t: str | None = None):
    cat = catalog()
    now = when(t)
    jd, fr = jd_array([now])
    r, _, err = cat.propagate(jd, fr)
    lat, lon, h = ecef_to_geodetic(teme_to_ecef(r[:, 0, :], jd, fr))
    ok = (err[:, 0] == 0) & np.isfinite(h)
    age = (now.timestamp() - state["epochs"]) / 86400
    return {"t": now.isoformat(), "n": int(ok.sum()),
            "norad": cat.norad[ok].tolist(), "lat": np.round(lat[ok], 3).tolist(), "lon": np.round(lon[ok], 3).tolist(),
            "alt": np.round(h[ok], 0).tolist(), "regime": state["regime"][ok].tolist(), "type": state["type"][ok].tolist(),
            "age": np.round(age[ok], 2).tolist()}


@app.get("/api/object/{norad}")
def obj(norad: int, t: str | None = None):
    cat = catalog()
    i = cat.index.get(norad)
    if i is None:
        raise HTTPException(404, "not in catalogue")
    e = cat.elements[i]
    now = when(t)
    age = (now - e.epoch).total_seconds() / 86400
    sc = conn.execute("SELECT name, intdes, object_type, ops_status, owner, launch_date, rcs FROM satcat WHERE norad=%s", (norad,)).fetchone()
    hist = conn.execute("SELECT epoch, fetched_at, sma_km, inc, source FROM gp WHERE norad=%s ORDER BY epoch", (norad,)).fetchall()
    man = conn.execute("SELECT epoch_before, epoch_after, d_sma_km, d_inc_deg FROM manoeuvre WHERE norad=%s ORDER BY epoch_after", (norad,)).fetchall()
    # one orbit, Earth-fixed, for drawing
    n = 180
    span = min(e.period_min, 1440) / 1440
    times = [now + timedelta(days=span * k / n) for k in range(n + 1)]
    jd, fr = jd_array(times)
    r, _, err = cat.propagate(jd, fr, subset=[i])
    lat, lon, h = ecef_to_geodetic(teme_to_ecef(r[0], jd, fr))
    sig = model.element_sigma(e, age)
    return {"norad": norad, "name": e.name, "intdes": e.intdes, "regime": e.regime, "orbit_source": e.source,
            "satcat": dict(zip(["name", "intdes", "type", "status", "owner", "launch", "rcs_m2"], [str(x) if x is not None else None for x in sc])) if sc else None,
            "orbit": {"perigee_km": round(e.perigee, 1), "apogee_km": round(e.apogee, 1), "inc_deg": e.inc, "period_min": round(e.period_min, 2), "ecc": e.ecc},
            "epoch": e.epoch.isoformat(), "age_days": round(age, 2),
            "uncertainty_1sigma_km": {"radial": round(float(sig[0]), 2), "in_track": round(float(sig[1]), 2), "cross_track": round(float(sig[2]), 2),
                                      "model": model.source},
            "history": [{"epoch": a.isoformat(), "fetched": b.isoformat(), "sma_km": round(c, 3), "inc": d, "source": s} for a, b, c, d, s in hist],
            "manoeuvres": [{"before": a.isoformat(), "after": b.isoformat(), "d_sma_km": round(c, 3), "d_inc_deg": round(d, 4)} for a, b, c, d in man],
            "track": [[round(float(x), 3), round(float(y), 3), round(float(z), 1)] for x, y, z in zip(lon, lat, h)]}


@app.get("/api/conjunctions")
def conjunctions(limit: int = 50, order: str = "pc"):
    run = conn.execute("SELECT id, created, window_start, days, threshold_km, objects, primaries, runtime_s, stats FROM conjunction_run ORDER BY id DESC LIMIT 1").fetchone()
    if not run:
        return {"run": None, "items": []}
    col = {"pc": "pc DESC NULLS LAST", "miss": "miss_km", "tca": "tca"}.get(order, "pc DESC NULLS LAST")
    rows = conn.execute(f"""SELECT c.a, sa.name, sa.object_type, c.b, sb.name, sb.object_type, c.tca, c.miss_km, c.rel_speed_kms, c.pc, c.pc_max,
                                   c.age_a_days, c.age_b_days, c.sigma_a_km, c.sigma_b_km, ST_Y(c.geom::geometry), ST_X(c.geom::geometry), ST_Z(c.geom::geometry)
                            FROM conjunction c LEFT JOIN satcat sa ON sa.norad=c.a LEFT JOIN satcat sb ON sb.norad=c.b
                            WHERE c.run_id=%s AND NOT c.colocated ORDER BY {col} LIMIT %s""", (run[0], limit)).fetchall()
    keys = ["a", "name_a", "type_a", "b", "name_b", "type_b", "tca", "miss_km", "rel_speed_kms", "pc", "pc_max", "age_a_days", "age_b_days",
            "sigma_a_km", "sigma_b_km", "lat", "lon", "alt_m"]
    total = conn.execute("SELECT count(*) FILTER (WHERE NOT colocated), count(*) FILTER (WHERE colocated) FROM conjunction WHERE run_id=%s", (run[0],)).fetchone()
    return {"run": dict(zip(["id", "created", "window_start", "days", "threshold_km", "objects", "primaries", "runtime_s", "stats"], run)),
            "conjunctions": total[0], "colocated": total[1],
            "items": [dict(zip(keys, r)) for r in rows]}


@app.get("/api/conjunction/{a}/{b}")
def conjunction(a: int, b: int):
    row = conn.execute("""SELECT tca, miss_km, rel_speed_kms, pc FROM conjunction WHERE run_id=(SELECT max(id) FROM conjunction_run)
                          AND ((a=%s AND b=%s) OR (a=%s AND b=%s)) ORDER BY pc DESC NULLS LAST LIMIT 1""", (a, b, b, a)).fetchone()
    if not row:
        raise HTTPException(404)
    cat = catalog()
    tca = row[0]
    times = [tca + timedelta(seconds=s) for s in range(-600, 601, 10)]
    jd, fr = jd_array(times)
    out = {}
    for n in (a, b):
        i = cat.index[n]
        r, _, _ = cat.propagate(jd, fr, subset=[i])
        lat, lon, h = ecef_to_geodetic(teme_to_ecef(r[0], jd, fr))
        out[str(n)] = [[round(float(x), 4), round(float(y), 4), round(float(z), 3)] for x, y, z in zip(lon, lat, h)]
    return {"tca": tca.isoformat(), "miss_km": row[1], "rel_speed_kms": row[2], "pc": row[3], "tracks": out}


@app.get("/api/geo")
def geo():
    snap = conn.execute("SELECT max(snapshot) FROM geo_state").fetchone()[0]
    rows = conn.execute("""SELECT norad, name, owner, object_type, lon, lat, drift_deg_day, inc, state FROM geo_state WHERE snapshot=%s
                           ORDER BY lon""", (snap,)).fetchall() if snap else []
    if state["geo_nb"] is None or state["geo_nb"][0] != snap:
        cat = catalog()
        jd, fr = jd_array([datetime.now(timezone.utc)])
        state["geo_nb"] = (snap, canadian_neighbours(cat, jd[0] + fr[0], state["satcat"]))
    keys = ["norad", "name", "owner", "type", "lon", "lat", "drift_deg_day", "inc", "state"]
    return {"snapshot": snap, "belt": [dict(zip(keys, r)) for r in rows], "canadian": state["geo_nb"][1]}


@app.get("/api/search")
def search(q: str, limit: int = 20):
    if q.isdigit():
        rows = conn.execute("SELECT norad, name FROM latest_gp WHERE norad=%s", (int(q),)).fetchall()
    else:
        rows = conn.execute("SELECT norad, name FROM latest_gp WHERE name ILIKE %s ORDER BY name LIMIT %s", (f"%{q}%", limit)).fetchall()
    return [{"norad": a, "name": b} for a, b in rows]


@app.get("/api/stats")
def stats():
    cat = catalog()
    now = datetime.now(timezone.utc).timestamp()
    age = (now - state["epochs"]) / 86400
    by = {}
    for k, v in REG.items():
        m = state["regime"] == v
        if m.any():
            by[k] = {"objects": int(m.sum()), "median_age_days": round(float(np.median(age[m])), 2),
                     "older_than_7d": int((age[m] > 7).sum())}
    log = conn.execute("SELECT source, finished, status, rows, new_rows, message FROM fetch_log ORDER BY id DESC LIMIT 20").fetchall()
    return {"objects": len(cat), "element_sets": conn.execute("SELECT count(*) FROM gp").fetchone()[0],
            "by_regime": by, "manoeuvres": conn.execute("SELECT count(*) FROM manoeuvre").fetchone()[0],
            "fetch_log": [dict(zip(["source", "finished", "status", "rows", "new", "message"], r)) for r in log],
            "uncertainty_model": model.source}


app.mount("/static", StaticFiles(directory=STATIC), name="static")


@app.get("/")
def index():
    return FileResponse(STATIC / "index.html")


def _warm():
    try:
        catalog()
    except Exception as e:  # database still empty
        print("warm-up:", e, flush=True)


threading.Thread(target=_warm, daemon=True).start()
