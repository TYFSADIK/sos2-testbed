"""Conjunction screening: who comes within a few kilometres of whom, and when.

1. Propagate every object on a fixed time grid (vectorised SGP4).
2. At each step, a k-d tree finds pairs closer than the threshold plus the
   distance two objects can close in half a step at the largest relative speed.
3. Linear relative motion inside the step estimates the closest approach; only
   pairs that could get within the threshold are kept.
4. Each surviving candidate is refined with SGP4 itself (bounded minimisation
   of the true separation) to a time of closest approach (TCA), miss distance
   and relative speed.

Pairs that never separate (docked modules, objects flying in formation) are
reported separately as co-located, not as conjunctions.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.spatial import cKDTree

from .catalog import Catalog

V_REL_MAX = 16.0          # km/s: head-on LEO crossing


@dataclass
class Conjunction:
    a: int                # NORAD ids
    b: int
    tca_jd: float         # Julian date (whole + fraction)
    miss_km: float
    rel_speed_kms: float
    r_a: np.ndarray
    v_a: np.ndarray
    r_b: np.ndarray
    v_b: np.ndarray

    def key(self):
        return (min(self.a, self.b), max(self.a, self.b))


def _sep(sat_a, sat_b, jd0, t):
    e1, r1, _ = sat_a.sgp4(jd0, t)
    e2, r2, _ = sat_b.sgp4(jd0, t)
    if e1 or e2:
        return 1e9
    return math.dist(r1, r2)


def refine(cat: Catalog, i: int, j: int, jd0: float, t_guess: float, half: float):
    """Minimise the true SGP4 separation within +-half (days) of t_guess (fraction of jd0)."""
    sa, sb = cat.elements[i].sat, cat.elements[j].sat
    res = minimize_scalar(lambda t: _sep(sa, sb, jd0, t), bounds=(t_guess - half, t_guess + half),
                          method="bounded", options={"xatol": 1e-9})
    t = float(res.x)
    e1, r1, v1 = sa.sgp4(jd0, t)
    e2, r2, v2 = sb.sgp4(jd0, t)
    r1, r2, v1, v2 = map(np.array, (r1, r2, v1, v2))
    return jd0 + t, float(np.linalg.norm(r1 - r2)), float(np.linalg.norm(v1 - v2)), r1, v1, r2, v2


def screen(cat: Catalog, jd0: float, days: float, threshold_km: float = 5.0, step_s: float = 20.0,
           primaries: np.ndarray | None = None, chunk_steps: int = 180, colocated_kms: float = 0.01,
           progress=None):
    """Screen [jd0, jd0 + days). primaries: boolean mask; a pair needs at least one primary."""
    n = len(cat)
    prim = np.ones(n, bool) if primaries is None else primaries
    steps = int(round(days * 86400 / step_s))
    dt = step_s / 86400.0
    pad = threshold_km + V_REL_MAX * step_s / 2
    candidates: dict = {}          # pair -> list of (t_frac, est_miss)
    stats = {"steps": steps, "objects": n, "primaries": int(prim.sum()), "kd_pairs": 0, "linear_kept": 0}
    for c0 in range(0, steps, chunk_steps):
        ks = np.arange(c0, min(steps, c0 + chunk_steps))
        fr = ks * dt
        r, v, err = cat.propagate(np.full(len(ks), jd0), fr)
        bad = err != 0
        for s_i, k in enumerate(ks):
            pos = r[:, s_i, :]
            ok = ~bad[:, s_i] & np.isfinite(pos).all(axis=1)
            idx = np.nonzero(ok)[0]
            tree = cKDTree(pos[idx])
            pairs = tree.query_pairs(pad, output_type="ndarray")
            if len(pairs) == 0:
                continue
            i, j = idx[pairs[:, 0]], idx[pairs[:, 1]]
            keep = prim[i] | prim[j]
            i, j = i[keep], j[keep]
            stats["kd_pairs"] += len(i)
            dr = pos[j] - pos[i]
            dv = v[j, s_i, :] - v[i, s_i, :]
            vv = (dv * dv).sum(1)
            tstar = np.where(vv > 1e-12, -(dr * dv).sum(1) / np.maximum(vv, 1e-12), 0.0)
            tstar = np.clip(tstar, -step_s / 2, step_s / 2)
            dmin = np.linalg.norm(dr + dv * tstar[:, None], axis=1)
            m = dmin < threshold_km * 1.5 + 0.5
            stats["linear_kept"] += int(m.sum())
            for a, b, ts, d in zip(i[m], j[m], tstar[m], dmin[m]):
                key = (int(min(a, b)), int(max(a, b)))
                candidates.setdefault(key, []).append((float(fr[s_i] + ts / 86400.0), float(d)))
        if progress:
            progress(min(steps, c0 + chunk_steps), steps)

    conj, colocated = [], []
    for (i, j), hits in candidates.items():
        hits.sort()
        # group hits closer than 10 minutes into one encounter
        groups, cur = [], [hits[0]]
        for h in hits[1:]:
            if h[0] - cur[-1][0] < 600 / 86400:
                cur.append(h)
            else:
                groups.append(cur)
                cur = [h]
        groups.append(cur)
        for g in groups:
            t_guess = min(g, key=lambda h: h[1])[0]
            jd_tca, miss, vrel, ra, va, rb, vb = refine(cat, i, j, jd0, t_guess, dt)
            if miss > threshold_km or not (0 <= jd_tca - jd0 < days):
                continue
            c = Conjunction(int(cat.norad[i]), int(cat.norad[j]), jd_tca, miss, vrel, ra, va, rb, vb)
            (colocated if vrel < colocated_kms else conj).append(c)
    stats["candidate_pairs"] = len(candidates)
    stats["conjunctions"] = len(conj)
    stats["colocated"] = len(colocated)
    return conj, colocated, stats
