# Architecture reference

```
CelesTrak GP + SATCAT ──(every 2 h, polite)──► ingest ──► PostgreSQL / PostGIS
                                                  │         gp            every element set ever fetched (history)
                                                  │         satcat        names, owners, types, sizes
                                                  │         fetch_log     every download, including "not updated" refusals
                                                  │         geo_state     GEO belt snapshots
                                                  │         manoeuvre     jumps between consecutive element sets
                                                  │         conjunction_run / conjunction (with encounter geography)
                                                  └──► derived products after every ingest
API (FastAPI) ──► CesiumJS viewer: positions, orbits, metadata, freshness, conjunctions, Canadian GEO
Grafana ──(SQL)──► catalogue health
tools/calibrate.py        GP vs operator SupGP and NASA ISS OEM  ──► results/uncertainty.json
tools/validate_socrates.py   screener vs SOCRATES on identical inputs
tools/optical_rehearsal.py   a night of GEO-belt frames, synthetic, scored against operator orbits
```

## Frames and time

SGP4 works in TEME. Earth-fixed positions use Greenwich mean sidereal time
(polar motion and the equation of the equinoxes are ignored: tens of metres,
far below GP accuracy). The optical pipeline uses Skyfield for topocentric
astrometric RA/Dec, the same frame as the star catalogue.

## Conjunction screening

Grid step 20 s; k-d tree radius 5 km + 16 km/s x 10 s; linear closest approach
inside the step; candidates with an estimated miss under 8 km refined with
SGP4 by bounded minimisation over +-1 step; encounters of the same pair within
10 minutes merged. Primaries are active payloads (SATCAT operational status
+, P, B, S, X); a pair needs one primary. Relative speed under 10 m/s is
co-located (docked, formation), not a conjunction.

Probability: each object's covariance from the uncertainty model (regime and
element-set age at TCA, radial / in-track / cross-track), combined and
projected on the plane normal to the relative velocity, integrated over a
hard-body circle whose radius comes from object type and radar size. The worst
case over covariance size (Alfano) is stored alongside.

## Uncertainty model

Per regime and axis, sigma(age) = a + b age (+ c age^2 in-track), fitted to
robust per-bin sigmas (1.4826 x median absolute error) of catalogue minus
truth. Manoeuvring objects are a separate population (their share over 50 km
is reported) and are caught by manoeuvre detection.

## Manoeuvre detection

For consecutive element sets of one object: the change in semi-major axis
minus the change drag predicts from the element set's own mean-motion
derivative, and the change in inclination. Thresholds: 0.6 km (LEO), 1.5 km
(MEO), 2 km (GEO), 5 km (HEO/GTO), 0.02 deg.

## Optical pipeline

Detect (5 sigma, connected components, moments) -> match round sources to the
star catalogue from the mount's pointing (translation vote, nearest
neighbour) -> cubic plate solution -> subtract every catalogue star with a
fitted zero point -> streaks of the expected length and direction, not
explained by a bright star or by star residuals -> RA/Dec at mid-exposure ->
correlate with catalogue predictions (400 arcsec gate) -> flag uncorrelated
tracks and residuals over 4 sigma of the uncertainty model.
