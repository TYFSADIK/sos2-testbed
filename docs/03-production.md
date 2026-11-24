# 3. Before anything real

| Area | Here | For real |
| --- | --- | --- |
| Catalogue | public CelesTrak groups (about 17,000 objects) | the full public catalogue (Space-Track, about 30,000+) or, for operations, the Space Surveillance Network's own catalogue with covariance |
| Element sets | GP/TLE, SGP4 | special-perturbations orbits with covariance from the observations themselves |
| Uncertainty | fitted from operator data per orbit regime | per-object covariance; manoeuvring objects get their own treatment |
| Collision probability | 2D Foster method with the fitted covariance and size from object type | operator-provided covariance, hard-body radius from real dimensions, CARA-style screening volumes |
| Screening | 20 s steps, linear pre-filter, SGP4 refinement | the same idea; validate after any change (tools/validate_socrates.py) |
| Optical sensor | rehearsed on synthetic frames | real frames with Gaia astrometry, photometric calibration, observation formats the network accepts (e.g. CCSDS TDM) |
| Security | Grafana anonymous viewer, API without login, database password in compose | behind an identity provider (see my ICAM testbed), secrets, TLS |
| Availability | one machine | the database replicated (see my ITI C2 testbed) |
