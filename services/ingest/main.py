"""Ingest: public orbital data into PostgreSQL, then the derived products.

MODE=live   fetch from CelesTrak every INTERVAL seconds (default 2 h, which is
            how often CelesTrak updates its groups; a group that has not
            updated is answered with a short notice, which is logged as
            "not-updated" and never stored as data)
MODE=files  load snapshots from /data/raw/snapshots/<time>/ in time order
            (rehearsal: reproducible, and no load on CelesTrak)

After each ingest: GEO belt snapshot, manoeuvre detection, and (every
CONJ_EVERY ingests) a conjunction screening run.
"""
import glob
import os
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from sos.catalog import read_gp_csv, read_satcat              # noqa: E402
from sos.jobs import load_catalog, run_conjunctions, run_geo, run_manoeuvres, satcat_map  # noqa: E402
from sos.store import connect, insert_elements, insert_satcat, insert_supgp, log_fetch  # noqa: E402
from sos.uncertainty import Model                              # noqa: E402

DB = os.environ.get("DATABASE_URL", "postgresql://sos:sos@postgis/sos")
MODE = os.environ.get("MODE", "files")
RAW = Path(os.environ.get("RAW_DIR", "/data/raw"))
INTERVAL = int(os.environ.get("INTERVAL", "7200"))
CONJ_DAYS = float(os.environ.get("CONJ_DAYS", "1"))
CONJ_EVERY = int(os.environ.get("CONJ_EVERY", "1"))
MODEL = Path(os.environ.get("UNCERTAINTY", "/app/results/uncertainty.json"))
UA = "sos2-testbed/0.1 (research; github.com/TYFSADIK/sos2-testbed)"
GROUPS = [("GROUP", "active"), ("SPECIAL", "gpz-plus"), ("GROUP", "cosmos-2251-debris"), ("GROUP", "iridium-33-debris"),
          ("GROUP", "fengyun-1c-debris"), ("GROUP", "last-30-days"), ("GROUP", "analyst")]
SUPGP = ["starlink", "oneweb", "planet", "iridium", "gps", "glonass", "intelsat", "ses", "telesat", "eumetsat", "iss", "cpf"]

conn = connect(DB)
model = Model(MODEL)
print(f"ingest mode={MODE} uncertainty={model.source}", flush=True)


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=1800) as r:
        return r.read().decode("utf-8", "replace")


def ingest_text(source, text, fetched_at):
    started = datetime.now(timezone.utc)
    try:
        els = read_gp_csv(text)
    except ValueError as e:
        log_fetch(conn, source, started, "not-updated", message=str(e))
        return 0
    new = insert_elements(conn, els, source, fetched_at)
    log_fetch(conn, source, started, "ok", len(els), new)
    return new


def derived(when, cycle):
    cat = load_catalog(conn)
    sc = satcat_map(conn)
    g = run_geo(conn, cat, sc, when)
    m = run_manoeuvres(conn)
    out = {"catalogue": len(cat), "geo": g, "manoeuvre_rows": m}
    if CONJ_EVERY and cycle % CONJ_EVERY == 0:
        out["conjunctions"] = run_conjunctions(conn, cat, sc, model, when, CONJ_DAYS)
    print(when.isoformat(), out, flush=True)


if MODE == "files":
    sat = RAW / "satcat.csv"
    if sat.exists():
        print("satcat rows", insert_satcat(conn, read_satcat(sat.read_text())), flush=True)
    for f in sorted(glob.glob(str(RAW / "sup-*.csv"))):
        started = datetime.now(timezone.utc)
        try:
            els = read_gp_csv(Path(f).read_text())
            n = insert_supgp(conn, els, Path(f).stem, datetime.fromtimestamp(Path(f).stat().st_mtime, timezone.utc))
            log_fetch(conn, Path(f).stem, started, "ok", len(els), n)
        except ValueError as e:
            log_fetch(conn, Path(f).stem, started, "not-updated", message=str(e))
    snaps = sorted(glob.glob(str(RAW / "snapshots" / "*")))
    for k, d in enumerate(snaps):
        fetched = datetime.strptime(Path(d).name, "%Y%m%dT%H%M%SZ").replace(tzinfo=timezone.utc)
        new = sum(ingest_text(Path(f).stem, Path(f).read_text(), fetched) for f in sorted(glob.glob(f"{d}/gp-*.csv")))
        print(f"snapshot {Path(d).name}: {new} new element sets", flush=True)
        if k == len(snaps) - 1 or os.environ.get("DERIVE_EACH") == "1":
            derived(fetched, k)
    conn.execute("INSERT INTO fetch_log (source, started, finished, status, message) VALUES ('files', now(), now(), 'done', %s)",
                 (f"{len(snaps)} snapshots",))
    while True:
        time.sleep(3600)
else:
    cycle = 0
    while True:
        now = datetime.now(timezone.utc)
        for key, g in GROUPS:
            try:
                ingest_text(g, fetch(f"https://celestrak.org/NORAD/elements/gp.php?{key}={g}&FORMAT=csv"), now)
            except Exception as e:                       # network trouble: log and carry on
                log_fetch(conn, g, now, "error", message=repr(e))
            time.sleep(2)
        for f in SUPGP:
            try:
                els = read_gp_csv(fetch(f"https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE={f}&FORMAT=csv"))
                log_fetch(conn, f"sup-{f}", now, "ok", len(els), insert_supgp(conn, els, f"sup-{f}", now))
            except ValueError as e:
                log_fetch(conn, f"sup-{f}", now, "not-updated", message=str(e))
            except Exception as e:
                log_fetch(conn, f"sup-{f}", now, "error", message=repr(e))
            time.sleep(2)
        if cycle % 12 == 0:
            try:
                insert_satcat(conn, read_satcat(fetch("https://celestrak.org/pub/satcat.csv")))
            except Exception as e:
                log_fetch(conn, "satcat", now, "error", message=repr(e))
        derived(now, cycle)
        cycle += 1
        time.sleep(INTERVAL)
