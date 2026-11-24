#!/usr/bin/env bash
# Download the public data the testbed uses into data/raw, politely (one request at a time).
# CelesTrak updates each group every 2 hours and answers a repeat download inside that time with a
# short notice instead of data; the files keep that notice, and the ingest logs it as "not-updated".
# Each run makes a new snapshot of the catalogue under data/raw/snapshots/<UTC time>/ (history).
set -euo pipefail
cd "$(dirname "$0")/../data/raw" 2>/dev/null || { mkdir -p "$(dirname "$0")/../data/raw"; cd "$(dirname "$0")/../data/raw"; }
UA="sos2-testbed/0.1 (research)"
get() {
  local out=$1 url=$2
  for _ in 1 2 3; do
    if curl -sS -m 1800 -A "$UA" -o "$out.part" "$url"; then mv "$out.part" "$out"; echo "ok $out $(wc -c < "$out")"; return; fi
    sleep 5
  done
  echo "FAIL $out"
}
SNAP=snapshots/$(date -u +%Y%m%dT%H%M%SZ)
mkdir -p "$SNAP"
for g in active gpz-plus cosmos-2251-debris iridium-33-debris fengyun-1c-debris last-30-days analyst; do
  key=GROUP; [ "$g" = gpz-plus ] && key=SPECIAL
  get "$SNAP/gp-$g.csv" "https://celestrak.org/NORAD/elements/gp.php?$key=$g&FORMAT=csv"
done
if [ "${SNAPSHOT_ONLY:-0}" = 1 ]; then exit 0; fi
for f in iss intelsat ses telesat eumetsat gps cpf starlink; do
  get "sup-$f.csv" "https://celestrak.org/NORAD/elements/supplemental/sup-gp.php?FILE=$f&FORMAT=csv"
done
get satcat.csv "https://celestrak.org/pub/satcat.csv"
get iss-oem.txt "https://nasa-public-data.s3.amazonaws.com/iss-coords/current/ISS_OEM/ISS.OEM_J2K_EPH.txt"
# SOCRATES: the first 1,000 conjunctions by time of closest approach, and the element sets it used for each pair
get socrates-tca.html "https://celestrak.org/SOCRATES/table-socrates.php?NAME=,&ORDER=TCA&MAX=1000"
mkdir -p socrates-gp
grep -oE "data\.php\?CATNR=[0-9]+,[0-9]+" socrates-tca.html | sort -u | while read -r q; do
  f="socrates-gp/${q#data.php?CATNR=}"; f="${f/,/-}.csv"
  [ -s "$f" ] && continue
  curl -sS -m 120 -A "$UA" -o "$f" "https://celestrak.org/SOCRATES/$q" || rm -f "$f"
  sleep 0.5
done
echo "done: $SNAP"
