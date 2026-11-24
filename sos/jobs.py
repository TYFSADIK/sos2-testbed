"""Derived products computed after each ingest: GEO belt snapshot, manoeuvre
detection from element-set history, and a conjunction screening run."""
from __future__ import annotations

import json
import math
import time
from datetime import datetime, timedelta, timezone

import numpy as np
from sgp4 import omm
from sgp4.api import Satrec

from .catalog import Catalog, Element, jd_array
from .conjunction import screen
from .geo import belt
from .uncertainty import Model, encounter_plane, hard_body_radius, pc_2d, pc_max

# manoeuvre thresholds: a change in semi-major axis (km) or inclination (deg) between consecutive
# element sets that drag (from the element set's own mean-motion derivative) does not explain
SMA_JUMP = {"LEO": 0.6, "MEO": 1.5, "GEO": 2.0, "HEO": 5.0, "GTO": 5.0}
INC_JUMP = 0.02


def _element(norad, name, intdes, fields, source="gp"):
    s = Satrec()
    omm.initialize(s, fields)
    return Element(norad, name, intdes, datetime.fromisoformat(fields["EPOCH"]).replace(tzinfo=timezone.utc), fields, s, source)


def load_catalog(conn, prefer_operator: bool = True, max_operator_age_days: float = 3.0) -> Catalog:
    """Best available orbit per object: the operator-derived SupGP when it is fresh and is the same
    object (international designator), otherwise the newest public GP element set. This is what
    SOCRATES does; inside manoeuvring constellations public GP alone produces many false conjunctions."""
    best = {}
    for norad, name, intdes, fields in conn.execute("SELECT norad, name, intdes, fields FROM latest_gp"):
        best[norad] = _element(norad, name, intdes, fields)
    if prefer_operator:
        rows = conn.execute("""SELECT DISTINCT ON (norad) norad, name, intdes, fields FROM supgp
                               WHERE epoch > now() - make_interval(days => %s) ORDER BY norad, epoch DESC""",
                            (int(max_operator_age_days),)).fetchall()
        for norad, name, intdes, fields in rows:
            g = best.get(norad)
            if g is not None and g.intdes == intdes:
                best[norad] = _element(norad, name, intdes, fields, "supgp")
    return Catalog(list(best.values()))


def satcat_map(conn) -> dict:
    return {n: {"OBJECT_TYPE": t or "", "OWNER": o or "", "OPS_STATUS_CODE": s or "", "RCS": r}
            for n, t, o, s, r in conn.execute("SELECT norad, object_type, owner, ops_status, rcs FROM satcat")}


def run_geo(conn, cat: Catalog, satcat: dict, when: datetime) -> int:
    jd, fr = jd_array([when])
    rows = belt(cat, jd[0], fr[0], satcat)
    with conn.cursor() as cur:
        cur.executemany("""INSERT INTO geo_state VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s, ST_SetSRID(ST_MakePoint(%s,%s),4326)::geography)
                           ON CONFLICT DO NOTHING""",
                        [(when, r["norad"], r["name"], r["owner"], r["type"], r["lon"], r["lat"], r["drift_deg_day"], r["inc"], r["state"],
                          r["lon"], r["lat"]) for r in rows])
    return len(rows)


def is_manoeuvre(regime, dt_days, a0, a1, inc0, inc1, ndot, n) -> tuple[bool, float]:
    """Between two element sets: did the orbit change more than drag explains?
    MEAN_MOTION_DOT is ndot/2 (rev/day^2); da = -2/3 a dn/n. Returns (flag, expected da)."""
    if dt_days <= 0 or dt_days > 10:
        return False, 0.0
    exp_da = -2.0 / 3.0 * a0 * (2 * ndot * dt_days) / n if n else 0.0
    return (abs((a1 - a0) - exp_da) > SMA_JUMP.get(regime, 2.0) or abs(inc1 - inc0) > INC_JUMP), exp_da


def run_manoeuvres(conn, since: datetime | None = None) -> int:
    """Compare each object's consecutive element sets."""
    q = """SELECT norad, epoch, sma_km, inc, regime, (fields->>'MEAN_MOTION_DOT')::float, mean_motion,
                  lag(epoch) OVER w, lag(sma_km) OVER w, lag(inc) OVER w
           FROM gp WINDOW w AS (PARTITION BY norad ORDER BY epoch)"""
    found = 0
    for norad, ep, a, inc, reg, ndot, n, ep0, a0, inc0 in conn.execute(q):
        if ep0 is None or (since and ep < since):
            continue
        flag, exp_da = is_manoeuvre(reg, (ep - ep0).total_seconds() / 86400, a0, a, inc0, inc, ndot, n)
        if flag:
            conn.execute("""INSERT INTO manoeuvre VALUES (%s,%s,%s,%s,%s,%s,%s,now()) ON CONFLICT DO NOTHING""",
                         (norad, ep0, ep, a - a0, inc - inc0, exp_da, reg))
            found += 1
    return found


def run_conjunctions(conn, cat: Catalog, satcat: dict, model: Model, start: datetime, days: float, threshold_km: float = 5.0,
                     step_s: float = 20.0, progress=None) -> dict:
    """Active payloads against everything, like SOCRATES; probability from the calibrated covariance."""
    t0 = time.time()
    prim = np.array([satcat.get(int(n), {}).get("OBJECT_TYPE") == "PAY" and satcat.get(int(n), {}).get("OPS_STATUS_CODE") in ("+", "P", "B", "S", "X")
                     for n in cat.norad])
    jd, fr = jd_array([start])
    conj, colo, stats = screen(cat, jd[0] + fr[0], days, threshold_km, step_s, primaries=prim, progress=progress)
    runtime = time.time() - t0
    run_id = conn.execute("""INSERT INTO conjunction_run (created, window_start, days, threshold_km, objects, primaries, runtime_s, stats)
                             VALUES (now(), %s, %s, %s, %s, %s, %s, %s) RETURNING id""",
                          (start, days, threshold_km, len(cat), int(prim.sum()), runtime, json.dumps(stats))).fetchone()[0]
    rows = []
    for c, colocated in [(x, False) for x in conj] + [(x, True) for x in colo]:
        ea, eb = cat.elements[cat.index[c.a]], cat.elements[cat.index[c.b]]
        tca = datetime(2000, 1, 1, 12, tzinfo=timezone.utc) + timedelta(days=c.tca_jd - 2451545.0)
        age_a = (tca - ea.epoch).total_seconds() / 86400
        age_b = (tca - eb.epoch).total_seconds() / 86400
        sa, sb = satcat.get(c.a, {}), satcat.get(c.b, {})
        hbr = hard_body_radius(sa.get("OBJECT_TYPE", ""), sb.get("OBJECT_TYPE", ""))
        pc = pmax = None
        if not colocated:
            ca = model.covariance_eci(c.r_a, c.v_a, ea.regime, age_a, ea)
            cb = model.covariance_eci(c.r_b, c.v_b, eb.regime, age_b, eb)
            miss2, P = encounter_plane(c.r_a, c.v_a, c.r_b, c.v_b)
            pc = pc_2d(miss2, P @ (ca + cb) @ P.T, hbr)
            pmax = pc_max(c.miss_km, hbr)
        # sub-satellite point of the encounter for the map
        from .catalog import ecef_to_geodetic, teme_to_ecef
        jdw = np.array([math.floor(c.tca_jd - 0.5) + 0.5])
        frw = np.array([c.tca_jd - jdw[0]])
        lat, lon, h = ecef_to_geodetic(teme_to_ecef(np.array([c.r_a]), jdw, frw))
        rows.append((run_id, c.a, c.b, tca, c.miss_km, c.rel_speed_kms, pc, pmax, age_a, age_b,
                     float(np.linalg.norm(model.element_sigma(ea, age_a))), float(np.linalg.norm(model.element_sigma(eb, age_b))),
                     hbr * 1000, colocated, float(lon[0]), float(lat[0]), float(h[0]) * 1000))
    with conn.cursor() as cur:
        cur.executemany("""INSERT INTO conjunction VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,
                           ST_SetSRID(ST_MakePoint(%s,%s,%s),4326)::geography)""", rows)
    sources = {}
    for c in conj:
        k = "-".join(sorted((cat.elements[cat.index[c.a]].source, cat.elements[cat.index[c.b]].source)))
        sources[k] = sources.get(k, 0) + 1
    conn.execute("UPDATE conjunction_run SET stats = stats || %s::jsonb WHERE id = %s",
                 (json.dumps({"by_orbit_source": sources, "operator_orbits": int(sum(e.source == "supgp" for e in cat.elements))}), run_id))
    return {"run_id": run_id, "runtime_s": round(runtime, 1), **stats, "by_orbit_source": sources}
