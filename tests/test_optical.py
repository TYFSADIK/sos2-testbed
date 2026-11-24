import math

import numpy as np

from sensor.optical import Camera, correlate, detect, render, solve, star_field, streaks, subtract_stars


def test_frame_solves_and_measures_a_streak():
    cam = Camera(width=1200, height=1200)
    ra0, dec0 = math.radians(150.0), math.radians(-6.0)
    stars = star_field(ra0, dec0, 1.2, seed=3)
    rate = -15.04 * math.pi / 180 / 3600          # GEO against the stars, rad/s in RA
    sat = {"ra": ra0 + math.radians(0.2), "dec": dec0 + math.radians(0.1), "mag": 11.0, "rate_ra": rate, "rate_dec": 0.0}
    img, truth = render(cam, ra0 + math.radians(0.05), dec0 - math.radians(0.05), 0.002, stars, [sat], seed=1)
    dets, _ = detect(img)
    plate = solve(dets, *stars, ra0, dec0, cam)
    assert plate is not None and plate.rms_arcsec < 3.0
    clean, _ = subtract_stars(img, plate, *stars, dets, cam)
    sx, sy = plate.to_px(stars[0], stars[1])
    found = streaks(detect(clean)[0], math.atan2(truth[0]["dy"], truth[0]["dx"]), expected_len_px=math.hypot(truth[0]["dx"], truth[0]["dy"]),
                    star_xy=np.column_stack([sx, sy]), star_mag=stars[2])
    assert len(found) == 1
    ra, dec = plate.to_radec(found[0]["x"], found[0]["y"])
    mid_ra = sat["ra"] + rate * cam.exposure_s / 2
    err = math.hypot((ra - mid_ra) * math.cos(dec), dec - sat["dec"]) * 206265
    assert err < 3.0
    c = correlate([(ra, dec)], [{"norad": 1, "ra": mid_ra, "dec": sat["dec"]}, {"norad": 2, "ra": mid_ra + 0.01, "dec": sat["dec"]}])
    assert c[0]["norad"] == 1
    assert np.isfinite(img).all()
