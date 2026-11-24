import math
from datetime import datetime, timezone

import numpy as np
import pytest

from sos.catalog import Catalog, ecef_to_geodetic, jd_array, read_gp_csv, regime, teme_to_ecef
from sos.conjunction import screen
from sos.geo import classify, drift_deg_per_day
from sos.socrates import parse
from sos.uncertainty import Model, encounter_plane, pc_2d, pc_max

HDR = "OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT\n"
ISS = "ISS (ZARYA),1998-067A,2026-11-23T04:53:05.576064,15.49560177,.0007419,51.6320,233.1644,56.4417,303.7279,0,U,25544,999,58665,.13296E-3,.7028E-4,0\n"
GEO = "NIMIQ 5,2009-050A,2026-11-23T10:00:00.000000,1.00271000,.0002000,0.0500,90.0,10.0,200.0,0,U,35873,999,6000,0,-.1E-6,0\n"


def test_parse_and_regimes():
    els = read_gp_csv(HDR + ISS + GEO)
    assert [e.regime for e in els] == ["LEO", "GEO"]
    assert 410 < els[0].perigee < 430 and abs(els[1].sma - 42164) < 5
    assert regime(200, 35786, 0.73, 630, 27) == "GTO"


def test_not_updated_notice_is_refused():
    with pytest.raises(ValueError):
        read_gp_csv("GP data has not updated since your last successful\ndownload of GROUP=active")


def test_geodetic_round_trip_on_equator():
    lat, lon, h = ecef_to_geodetic(np.array([[42164.0, 0.0, 0.0]]))
    assert abs(lat[0]) < 1e-9 and abs(lon[0]) < 1e-9 and abs(h[0] - (42164 - 6378.137)) < 1e-6


def test_iss_position_is_sensible():
    cat = Catalog(read_gp_csv(HDR + ISS))
    jd, fr = jd_array([datetime(2026, 11, 23, 12, tzinfo=timezone.utc)])
    r, v, e = cat.propagate(jd, fr)
    lat, lon, h = ecef_to_geodetic(teme_to_ecef(r[:, 0, :], jd, fr))
    assert e[0, 0] == 0 and 400 < h[0] < 440 and abs(lat[0]) <= 51.7
    assert abs(np.linalg.norm(v[0, 0]) - 7.66) < 0.05


def test_geo_classification():
    assert abs(drift_deg_per_day(1.00273791)) < 1e-9
    assert classify(0.01, 0.05) == "station-kept" and classify(0.01, 2.0) == "inclined station-kept"
    assert classify(1.5, 0.05) == "drifting"


def test_screen_finds_a_constructed_crossing():
    # two circular orbits at the same altitude in different planes cross twice a revolution
    a = "SAT A,2020-001A,2026-11-23T00:00:00.000000,15.20000000,.0001,53.0,10.0,0.0,0.0,0,U,90001,999,1,0,0,0\n"
    b = "SAT B,2020-002A,2026-11-23T00:00:00.000000,15.20000000,.0001,97.0,10.0,0.0,0.0,0,U,90002,999,1,0,0,0\n"
    cat = Catalog(read_gp_csv(HDR + a + b))
    jd, fr = jd_array([datetime(2026, 11, 23, tzinfo=timezone.utc)])
    conj, colo, st = screen(cat, jd[0] + fr[0], 0.2, threshold_km=50.0, step_s=20.0)
    assert conj, st
    assert min(c.miss_km for c in conj) < 50 and all(c.rel_speed_kms > 5 for c in conj)


def test_probability_behaves():
    m = np.array([0.0, 0.0])
    assert pc_2d(m, np.eye(2) * 1e-6, 0.01) > 0.99          # tiny covariance, dead centre: certain
    far = pc_2d(np.array([5.0, 0.0]), np.eye(2) * 0.01, 0.01)
    assert far < 1e-12
    assert abs(pc_max(1.0, 0.01) - 0.01 ** 2 / (math.e * 1.0)) < 1e-12
    miss, P = encounter_plane([7000, 0, 0], [0, 7.5, 0], [7000.5, 0.1, 0], [0, 0, 7.5])
    assert P.shape == (2, 3) and abs(np.linalg.norm(miss) - math.hypot(0.5, 0.1)) < 0.01


def test_uncertainty_grows_with_age():
    m = Model("/nonexistent")
    assert (m.sigma("LEO", 5) > m.sigma("LEO", 0)).all()


def test_socrates_parse():
    page = """Data current as of 2026 Nov 23 18:15:04 UTC <p>Computation Interval: Start = 2026 Nov 23 19:00:00 UTC, Stop = 2026 Nov 30 19:00:00 UTC<br>
    Considering: 16,902 Primaries, 32,971 Secondaries (196,692 Conjunctions)<br>
    <tr align=center><td>x</td><td>25544</td><td>ISS (ZARYA) [+]</td><td>1.334</td><td rowspan=2>2026-11-24 12:54:15.791</td><td>0.432</td><td>11.200</td></tr>
    <tr align=center><td>y</td><td>12345</td><td>DEB [-]</td><td>2.000</td><td>1.000E-05</td><td>0.100</td></tr>"""
    d = parse(page)
    c = d["conjunctions"][0]
    assert d["meta"]["primaries"] == 16902 and c["norad1"] == 25544 and c["status1"] == "+" and c["range_km"] == 0.432


def test_manoeuvre_detection():
    from sos.jobs import is_manoeuvre
    # LEO object decaying as its own drag term predicts: not a manoeuvre
    n, a0 = 15.5, 6790.0
    ndot = 5e-5                                      # MEAN_MOTION_DOT (ndot/2), rev/day^2
    exp = -2 / 3 * a0 * (2 * ndot * 1.0) / n
    assert is_manoeuvre("LEO", 1.0, a0, a0 + exp, 51.6, 51.6, ndot, n)[0] is False
    # a 2 km reboost on top of the decay: manoeuvre
    assert is_manoeuvre("LEO", 1.0, a0, a0 + exp + 2.0, 51.6, 51.6, ndot, n)[0] is True
    # a GEO inclination change of 0.05 deg: manoeuvre; 1 km in GEO: within noise
    assert is_manoeuvre("GEO", 2.0, 42164.0, 42164.0, 0.05, 0.10, 0.0, 1.0027)[0] is True
    assert is_manoeuvre("GEO", 2.0, 42164.0, 42165.0, 0.05, 0.05, 0.0, 1.0027)[0] is False
