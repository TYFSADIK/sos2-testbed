"""PostgreSQL/PostGIS store: every element set ever fetched, the satellite
catalogue, fetch log, conjunction runs, GEO belt snapshots and manoeuvres."""
from __future__ import annotations

import json
import time
from datetime import datetime, timezone

import psycopg

SCHEMA = """
CREATE EXTENSION IF NOT EXISTS postgis;
CREATE TABLE IF NOT EXISTS gp (
  norad int, epoch timestamptz, fetched_at timestamptz, source text, name text, intdes text,
  mean_motion double precision, ecc double precision, inc double precision, raan double precision,
  argp double precision, mean_anomaly double precision, bstar double precision, sma_km double precision,
  perigee_km double precision, apogee_km double precision, regime text, fields jsonb,
  PRIMARY KEY (norad, epoch));
CREATE INDEX IF NOT EXISTS gp_norad_epoch ON gp (norad, epoch DESC);
CREATE TABLE IF NOT EXISTS supgp (norad int, epoch timestamptz, fetched_at timestamptz, source text, name text, intdes text,
  rms_km double precision, data_source text, fields jsonb, PRIMARY KEY (norad, epoch));
CREATE TABLE IF NOT EXISTS satcat (norad int PRIMARY KEY, name text, intdes text, object_type text, ops_status text,
  owner text, launch_date date, decay_date date, period_min double precision, inclination double precision,
  apogee_km double precision, perigee_km double precision, rcs double precision, row jsonb);
CREATE TABLE IF NOT EXISTS fetch_log (id serial PRIMARY KEY, source text, started timestamptz, finished timestamptz,
  status text, rows int, new_rows int, message text);
CREATE TABLE IF NOT EXISTS conjunction_run (id serial PRIMARY KEY, created timestamptz, window_start timestamptz,
  days double precision, threshold_km double precision, objects int, primaries int, runtime_s double precision, stats jsonb);
CREATE TABLE IF NOT EXISTS conjunction (run_id int REFERENCES conjunction_run(id) ON DELETE CASCADE, a int, b int,
  tca timestamptz, miss_km double precision, rel_speed_kms double precision, pc double precision, pc_max double precision,
  age_a_days double precision, age_b_days double precision, sigma_a_km double precision, sigma_b_km double precision,
  hbr_m double precision, colocated bool DEFAULT false, geom geography(PointZ, 4326));
CREATE INDEX IF NOT EXISTS conjunction_run_pc ON conjunction (run_id, pc DESC);
CREATE TABLE IF NOT EXISTS geo_state (snapshot timestamptz, norad int, name text, owner text, object_type text,
  lon double precision, lat double precision, drift_deg_day double precision, inc double precision, state text,
  geom geography(Point, 4326), PRIMARY KEY (snapshot, norad));
CREATE TABLE IF NOT EXISTS manoeuvre (norad int, epoch_before timestamptz, epoch_after timestamptz, d_sma_km double precision,
  d_inc_deg double precision, expected_d_sma_km double precision, regime text, detected_at timestamptz,
  PRIMARY KEY (norad, epoch_after));
CREATE OR REPLACE VIEW latest_gp AS
  SELECT DISTINCT ON (norad) * FROM gp ORDER BY norad, epoch DESC;
"""


def connect(url: str, tries: int = 60):
    for _ in range(tries):
        try:
            conn = psycopg.connect(url, autocommit=True)
            conn.execute("SELECT pg_advisory_lock(7342)")      # one service at a time creates the schema
            try:
                conn.execute(SCHEMA)
            finally:
                conn.execute("SELECT pg_advisory_unlock(7342)")
            return conn
        except psycopg.OperationalError:
            time.sleep(1)
    raise SystemExit("database not reachable")


def insert_elements(conn, elements, source: str, fetched_at: datetime | None = None) -> int:
    fetched_at = fetched_at or datetime.now(timezone.utc)
    rows = [(e.norad, e.epoch, fetched_at, source, e.name, e.intdes, e.mean_motion, e.ecc, e.inc,
             float(e.fields["RA_OF_ASC_NODE"]), float(e.fields["ARG_OF_PERICENTER"]), float(e.fields["MEAN_ANOMALY"]),
             float(e.fields["BSTAR"]), e.sma, e.perigee, e.apogee, e.regime, json.dumps(e.fields)) for e in elements]
    before = conn.execute("SELECT count(*) FROM gp").fetchone()[0]
    with conn.cursor() as cur:
        cur.executemany("""INSERT INTO gp VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
                           ON CONFLICT (norad, epoch) DO NOTHING""", rows)
    return conn.execute("SELECT count(*) FROM gp").fetchone()[0] - before


def insert_supgp(conn, elements, source: str, fetched_at: datetime) -> int:
    """Operator-derived orbits (CelesTrak SupGP), kept apart from the public GP history."""
    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None
    rows = [(e.norad, e.epoch, fetched_at, source, e.name, e.intdes, num(e.fields.get("RMS")), e.fields.get("DATA_SOURCE"),
             json.dumps(e.fields)) for e in elements]
    before = conn.execute("SELECT count(*) FROM supgp").fetchone()[0]
    with conn.cursor() as cur:
        cur.executemany("INSERT INTO supgp VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s) ON CONFLICT DO NOTHING", rows)
    return conn.execute("SELECT count(*) FROM supgp").fetchone()[0] - before


def insert_satcat(conn, satcat: dict) -> int:
    def num(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None

    def date(x):
        return x or None
    rows = [(n, r.get("OBJECT_NAME"), r.get("OBJECT_ID"), r.get("OBJECT_TYPE"), r.get("OPS_STATUS_CODE"), r.get("OWNER"),
             date(r.get("LAUNCH_DATE")), date(r.get("DECAY_DATE")), num(r.get("PERIOD")), num(r.get("INCLINATION")),
             num(r.get("APOGEE")), num(r.get("PERIGEE")), num(r.get("RCS")), json.dumps(r)) for n, r in satcat.items()]
    with conn.cursor() as cur:
        cur.executemany("""INSERT INTO satcat VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)
                           ON CONFLICT (norad) DO UPDATE SET ops_status = excluded.ops_status, decay_date = excluded.decay_date,
                           row = excluded.row""", rows)
    return len(rows)


def log_fetch(conn, source, started, status, rows=0, new_rows=0, message=""):
    conn.execute("INSERT INTO fetch_log (source, started, finished, status, rows, new_rows, message) VALUES (%s,%s,now(),%s,%s,%s,%s)",
                 (source, started, status, rows, new_rows, message[:500]))
