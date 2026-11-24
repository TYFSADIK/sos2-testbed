"""Parse CelesTrak SOCRATES Plus result pages (the independent reference for
conjunction screening)."""
from __future__ import annotations

import html
import re
from datetime import datetime, timezone

ROW = re.compile(r"<tr[^>]*align=center[^>]*>(.*?)</tr>", re.S)   # alternate rows carry shading attributes
CELL = re.compile(r"<td[^>]*>(.*?)</td>", re.S)


def _txt(s):
    return html.unescape(re.sub(r"<[^>]+>", "", s)).strip()


def parse(page: str) -> dict:
    meta = {}
    m = re.search(r"Data current as of ([^<]+) UTC", page)
    meta["data_as_of"] = m.group(1).strip() if m else None
    m = re.search(r"Start = ([^,]+) UTC, Stop = ([^<]+) UTC", page)
    if m:
        meta["start"], meta["stop"] = m.group(1).strip(), m.group(2).strip()
    m = re.search(r"Considering: ([\d,]+) Primaries, ([\d,]+) Secondaries \(([\d,]+) Conjunctions\)", page)
    if m:
        meta["primaries"], meta["secondaries"], meta["conjunctions_total"] = (int(x.replace(",", "")) for x in m.groups())
    rows = [c for c in (CELL.findall(r) for r in ROW.findall(page)) if len(c) >= 6]
    out = []
    for a, b in zip(rows[0::2], rows[1::2]):
        if len(a) < 7 or len(b) < 6:
            continue
        a, b = [_txt(x) for x in a], [_txt(x) for x in b]
        name1, status1 = _split_status(a[2])
        name2, status2 = _split_status(b[2])
        out.append({
            "norad1": int(a[1]), "name1": name1, "status1": status1, "dse1": float(a[3]),
            "tca": datetime.strptime(a[4], "%Y-%m-%d %H:%M:%S.%f").replace(tzinfo=timezone.utc).isoformat(),
            "range_km": float(a[5]), "rel_speed_kms": float(a[6]),
            "norad2": int(b[1]), "name2": name2, "status2": status2, "dse2": float(b[3]),
            "max_prob": float(b[4]), "dilution_km": float(b[5]) if b[5] not in ("", "N/A") else None,
        })
    return {"meta": meta, "conjunctions": out}


def _split_status(s):
    m = re.match(r"(.*?)\s*\[(.)\]\s*$", s)
    return (m.group(1), m.group(2)) if m else (s, "")
