#!/usr/bin/env bash
# Rehearse the space surveillance testbed on one machine with Docker.
#
#   ./local/fetch.sh            download public data (CelesTrak, NASA) into data/raw (run first; slow and rate-limited)
#   ./local/rehearse.sh up      build and start: PostgreSQL/PostGIS, ingest, API + 3D viewer, Grafana
#   ./local/rehearse.sh test    the full rehearsal; results in results/$RUN
#   ./local/rehearse.sh down    stop and delete everything
set -euo pipefail
cd "$(dirname "$0")/../stacks"
RUN=${RUN:-rehearsal}
OUT=../results/$RUN
PW=${PLAYWRIGHT:-}                 # path to playwright-core, for the screenshots (optional)
API=http://127.0.0.1:8000
mkdir -p "$OUT"
tools() { docker compose --profile tools run --rm -T tools "$@"; }

case "${1:-}" in
  up)
    docker compose up -d --build
    echo "viewer $API   grafana http://127.0.0.1:3000" ;;
  down)
    docker compose --profile tools down -v --remove-orphans ;;
  test)
    docker compose --profile tools down -v --remove-orphans >/dev/null 2>&1 || true
    docker compose build -q   # set BASE=<image> to build on an existing Python image

    echo "== 1. uncertainty: catalogue vs operator orbits and NASA's ISS ephemeris"
    tools python tools/calibrate.py data/raw "results/$RUN" | tee "$OUT/calibrate.txt"
    cp "$OUT/uncertainty.json" ../results/uncertainty.json

    echo "== 2. ingest every snapshot (history), GEO belt, manoeuvres, 7-day conjunction screening"
    CONJ_DAYS=${CONJ_DAYS:-7} docker compose up -d >/dev/null
    for _ in $(seq 1 360); do
      docker compose logs ingest 2>&1 | grep -q "'conjunctions'" && break
      docker compose logs ingest 2>&1 | grep -q Traceback && { docker compose logs ingest | tail -30; exit 1; }
      sleep 5
    done
    docker compose logs ingest 2>&1 | grep -E "snapshot|conjunctions|satcat" | sed 's/^ingest-1  | //' | tee "$OUT/ingest.txt"
    for _ in $(seq 1 60); do curl -fsS "$API/api/stats" 2>/dev/null | grep -q '"objects":[1-9]' && break; sleep 2; done

    echo "== 3. conjunction screening against CelesTrak SOCRATES (identical inputs)"
    tools python tools/validate_socrates.py data/raw "results/$RUN/socrates" | tail -2 | tee "$OUT/socrates.txt"

    echo "== 4. optical sensor: a night of GEO-belt observing from Toronto"
    tools python tools/optical_rehearsal.py data/raw "results/$RUN/optical" "results/$RUN/uncertainty.json" | tail -40 > "$OUT/optical.txt"

    echo "== 5. platform"
    curl -fsS "$API/api/stats" > "$OUT/stats.json"
    curl -fsS "$API/api/geo" > "$OUT/geo.json"
    curl -fsS "$API/api/conjunctions?limit=100" > "$OUT/conjunctions.json"
    curl -fsS "$API/api/object/25544" > "$OUT/object-iss.json"
    curl -fsS "$API/api/object/28868" > "$OUT/object-anik-f1r.json"
    python3 - "$API" > "$OUT/api-latency.json" <<'PY'
import json, statistics, sys, time, urllib.request
api = sys.argv[1]
def t(path, n):
    xs = []
    for _ in range(n):
        s = time.perf_counter(); urllib.request.urlopen(api + path).read(); xs.append((time.perf_counter() - s) * 1000)
    xs.sort()
    return {"n": n, "median_ms": round(statistics.median(xs), 1), "p95_ms": round(xs[int(0.95 * (n - 1))], 1)}
print(json.dumps({"positions_all_objects": t("/api/positions", 30), "object_detail": t("/api/object/25544", 30),
                  "conjunctions_top50": t("/api/conjunctions?limit=50", 30), "search": t("/api/search?q=ANIK", 30)}))
PY
    docker compose exec -T postgis psql -U sos -tA -F, -c \
      "select (select count(*) from gp), (select count(*) from latest_gp), (select count(*) from manoeuvre), (select count(distinct snapshot) from geo_state), (select count(*) from fetch_log where status='not-updated')" > "$OUT/db.csv"
    # conjunctions per hour by object family, mine (7-day run) against the SOCRATES hour
    docker compose exec -T postgis psql -U sos -tA -F, -c "
      select case when sa.name like 'STARLINK%' and sb.name like 'STARLINK%' then 'starlink-starlink'
                  when sa.name like 'STARLINK%' or sb.name like 'STARLINK%' then 'starlink-other' else 'other' end,
             round(count(*) / (select days * 24 from conjunction_run order by id desc limit 1)::numeric, 1)
      from conjunction c left join satcat sa on sa.norad=c.a left join satcat sb on sb.norad=c.b
      where not colocated and run_id=(select max(id) from conjunction_run) group by 1 order by 1" > "$OUT/rate-per-hour.csv"
    tools python -c "
import sys, collections; sys.path.insert(0, '.')
from sos.socrates import parse
s = [x for x in parse(open('data/raw/socrates-tca.html').read())['conjunctions'] if x['rel_speed_kms'] >= 0.01]
f = collections.Counter('starlink-starlink' if 'STARLINK' in x['name1'] and 'STARLINK' in x['name2'] else 'starlink-other' if 'STARLINK' in x['name1'] + x['name2'] else 'other' for x in s)
print('\n'.join(f'{k},{round(v / 62.2 * 60, 1)}' for k, v in sorted(f.items())))" > "$OUT/rate-per-hour-socrates.csv"
    if [ -n "$PW" ]; then
      node ../tools/shots.cjs "$PW" "$OUT" | tee "$OUT/shots.json"
    fi
    echo "results: $OUT" ;;
  *)
    echo "usage: $0 up|test|down" >&2; exit 2 ;;
esac
