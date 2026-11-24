# 2. Deploy

## On one machine

```bash
./local/fetch.sh            # public data into data/raw (slow: CelesTrak serves about 15 KB/s and updates every 2 h)
./local/rehearse.sh up      # viewer http://127.0.0.1:8000, Grafana http://127.0.0.1:3000
./local/rehearse.sh test    # the full rehearsal (about 25 minutes)
```

Docker with 4 CPUs and 8 GB. `BASE=<image>` builds the Python image on top of
an existing one (saves downloading NumPy and SciPy on a slow link).

## Live mode

`MODE=live docker compose up -d` makes the ingest fetch CelesTrak itself every
`INTERVAL` seconds (default 7,200, CelesTrak's update cycle) and rerun the
derived products. Every element set is kept, so history (and manoeuvre
detection) grows with every cycle. `CONJ_DAYS` sets the screening window.

Please keep `INTERVAL` at 2 hours or more: CelesTrak is a free service, run by
one person, and refuses repeat downloads inside the update cycle anyway.

## Data sources and terms

| Source | What | Notes |
| --- | --- | --- |
| CelesTrak GP | element sets for active satellites, the GEO protected zone, debris clouds, recent launches, analyst objects | the full public catalogue needs a Space-Track account; with one, add it as another source |
| CelesTrak SATCAT | names, owners, launch dates, object types, sizes | |
| CelesTrak SupGP | orbits fitted to operators' ephemerides and laser-ranging predictions | used only as truth for calibration |
| CelesTrak SOCRATES | conjunction reports | used only for validation |
| NASA ISS OEM | the ISS ephemeris NASA publishes | used only for validation |

Credit CelesTrak and NASA when you publish results.
