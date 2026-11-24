# 1. Physical setup: an optical sensor for the GEO belt

Surveillance of Space 2 is about sensors that detect, track and characterise
objects in deep space. The hobby-scale version is a camera and a long lens on a
tracking mount: geostationary satellites are bright enough (magnitude 10 to 13)
to record in two seconds from a suburban back yard, and each frame measures
their positions to about a kilometre, independently of the public catalogue.

The software that turns frames into measurements is in `sensor/optical.py`;
`tools/optical_rehearsal.py` rehearses a full night with synthetic frames. This
page is the hardware.

## 1.1 Shopping list

| Item | Example | About (CAD) | Why |
| --- | --- | --- | --- |
| Lens | Samyang 135 mm f/2 or a used 200 mm f/2.8 | 500 | 3 to 4 arcsec per pixel, a 3 degree field: several GEO satellites per frame |
| Astronomy camera | ZWO ASI533MM (3008 x 3008, 3.76 um, mono, low read noise) | 900 | The detector the rehearsal models |
| Lens adapter | M42 or EF to the camera, with the right back-focus | 60 | |
| Tracking mount | Sky-Watcher Star Adventurer GTi or any equatorial mount with guiding | 600 | Tracks the stars, so GEO satellites draw short streaks |
| GNSS time receiver (optional) | USB GNSS puck with PPS, or the mount's GPS | 40 | Exposure timestamps within a millisecond (a GEO satellite moves 3 m per ms) |
| Laptop or Raspberry Pi 5 | runs capture (INDI/KStars or N.I.N.A.) | - | |
| Dew heater strap, power bank | | 60 | Canadian winters |

About $2,100 new; less with a used lens and mount.

## 1.2 Set up

1. Polar-align the mount (within 5 arcmin is plenty for 2 s exposures).
2. Focus on a bright star until it is smallest (a Bahtinov mask helps).
3. Time: run `chrony` with the GNSS receiver as a reference, or at least NTP.
   Write the mid-exposure time into each frame; a 0.1 s timing error is a
   300 m along-track error at GEO.
4. From Toronto (43.7 N) the belt crosses the southern sky at about 38 degrees
   elevation; point at the slot you want (the planner in the rehearsal picks
   the Canadian slots from 72 W to 119 W).

## 1.3 Observe

* Exposure 2 s, gain at unity, tracking at the sidereal rate.
* Take 10 to 20 frames per field, a few minutes apart: each satellite gives a
  short arc instead of a single point, which is what an orbit update needs.
* Avoid the two weeks around the equinoxes near local midnight: GEO satellites
  pass through Earth's shadow and disappear (eclipse season).

## 1.4 Reduce

For real frames the pipeline needs a real star catalogue in place of the
rehearsal's synthetic sky:

* plate-solve each frame with `astrometry.net` (index files 4107 to 4110 for a
  3 degree field) to get a first pointing, or use the mount's reported
  pointing;
* use Gaia DR3 stars from the field (via `astroquery`) as the catalogue the
  pipeline matches and subtracts.

Then `sensor/optical.py` does what the rehearsal measured: cubic plate
solution, star subtraction, streak measurement, correlation with the
catalogue, and flags for uncorrelated tracks and large residuals.

## 1.5 What to expect

From the rehearsal (30 fields, synthetic frames with lens distortion, seeing
and catalogue error): about 1 arcsec median position error (about 0.2 km at
GEO range), everything brighter than magnitude 14 detected in 2 s, debris
fainter than that not detected. Real skies add clouds, the Moon and light
pollution; expect fewer satellites per night, not worse positions.

## Safety and etiquette

Observing is passive and legal. Do not publish observations of classified
objects' positions; the public catalogue already omits them for a reason.
