# 4. What the rehearsal found

Every finding below came from running the code against real public data or the
rehearsal, not from reading documentation.

## Data

| What happened | Cause | What changed |
| --- | --- | --- |
| The first download of the active catalogue came back as 141 bytes | CelesTrak updates each group every 2 hours and answers a repeat download with a notice; a reachability probe had counted as the download | the ingest recognises the notice, logs it as "not-updated" and never stores it; the first snapshot keeps it as a logged example |
| No full public catalogue to download | the complete catalogue needs a Space-Track account | public groups (active, GEO protected zone, three debris clouds, recent launches, analyst objects): about 21,000 objects against SOCRATES's 33,000 |
| Catalogue numbers reused across sources | temporary numbers (e.g. 100001) appear in more than one source | matches also require the same international designator (26 Starlink records rejected) |
| SOCRATES's "top 1,000" lists were almost all docked spacecraft | ISS modules and visiting vehicles "conjunct" at 2 m every orbit | validation uses the list sorted by time instead, and co-located pairs (under 10 m/s) are reported separately on both sides |
| The SOCRATES parser read only half the conjunctions | alternate result rows carry a shading attribute the row pattern did not allow | pattern fixed; until then the missing half looked like 476 "extra" conjunctions of mine, and a wrong hypothesis (one entry per pair per week) was tested and rejected first |

## Screening

| What happened | Cause | What changed |
| --- | --- | --- |
| Same inputs as SOCRATES: times of closest approach within 0.8 ms, miss distances within 0.5 m (958 conjunctions) | | the 0.5 m is SOCRATES's own rounding (kilometres to three decimals) |
| Same hour, screened independently: all 917 SOCRATES conjunctions found, 5 that SOCRATES did not list | the 5 are real under SGP4 with identical element sets; the cause is not determined (SOCRATES uses STK/CAT, whose pre-filters may drop a few) | reported, not claimed as SOCRATES errors |
| That validation could not see over-reporting | its catalogue contained only objects SOCRATES listed in that hour | a per-hour comparison by object family against the 7-day run was added |
| The 7-day screen reported 808,531 conjunctions, about 4 times SOCRATES's rate | almost all Starlink against Starlink: 4,903 an hour against the 182 SOCRATES lists. Starlink against other objects matches (432 against 498 an hour) and so do the rest once catalogue size is allowed for (146 against 260, with 21,000 objects against 33,000) | first hypothesis: public elements for Starlink are too inaccurate (median 13.6 km off after a day) and SOCRATES screens operator orbits instead. Tested by preferring operator orbits (SupGP): Starlink against Starlink went up, not down (820,360 between two operator orbits), and all 189 Starlink pairs SOCRATES listed are between satellites that have operator orbits. Hypothesis rejected. The pairs SOCRATES does list, I match exactly; why it lists so few Starlink pairs is not documented and is **unresolved** |
| Operator orbits propagated five to six days ahead gave Starlink pairs 50 to 150 m apart at 15 km/s, and probabilities near 1 in 2,000 | the uncertainty given to operator orbits (fit RMS plus 0.3 km per day) was an assumption, and far too small beyond the span the orbits were fitted to | operator orbits keep their fit RMS for one day after epoch; after that the measured public-data model applies |
| No manoeuvres were detected on real data | three snapshots over two hours: CelesTrak had published new elements for almost no object in between | the detector is unit-tested; real history needs live mode running for days |
| The Grafana screenshot timed out | the map panel plotted every conjunction | the map shows the 5,000 with the highest probability |

## Uncertainty

| What happened | Cause | What changed |
| --- | --- | --- |
| Before the full catalogue, LEO error looked like thousands of kilometres | the only LEO matches were Starlinks launched in the last 30 days, still raising their orbits | robust per-bin sigma (median absolute error); the share over 50 km is reported separately as the manoeuvring population |
| The ISS's public elements drift from NASA's ephemeris faster than linearly (under 1 km for 3 days, 2.9 km at 5, 12.7 km at 7, 85 km at 14) | first guessed as a planned reboost; NASA's orbit-averaged semi-major axis shows none: a small mismatch in predicted decay accumulates along-track | in-track error fitted with a quadratic term |

## Optical sensor

| What happened | Cause | What changed |
| --- | --- | --- |
| Satellite positions were off by exactly 15 arcsec | truth taken at the start of the exposure, the streak measured at its middle | everything at mid-exposure |
| Then 0.075 arcsec: too good to be true | the synthetic sky had no lens distortion, seeing or catalogue error | all three added; a cubic plate model; about 1 arcsec now |
| Satellites on top of stars were missed | blends | stars subtracted using the plate solution and a fitted zero point |
| False streaks after subtraction | saturated bright stars and close double stars leave residuals | expected streak length from the tracking mode, a bright-star mask, and a residual-to-starlight test: 0 false streaks |
| An injected manoeuvre slid one satellite onto another | the test, not the pipeline | injected tests only where the sky is clear |
| One frame of 30 solved badly (197 of 1,541 stars, 3.5 arcsec) and put every satellite in it 100 arcsec off, causing false manoeuvre flags | the solver locked onto a partial match | a quality gate retries and otherwise rejects the frame; why this frame fails is open |
| Two satellites sharing a slot 11 px apart, and magnitude-15 debris | below the optics' resolution and a 2 s exposure's depth | limits, reported |

## Platform

| What happened | Cause | What changed |
| --- | --- | --- |
| Two services crashed creating the schema at the same moment | both ran `CREATE ... IF NOT EXISTS` concurrently | an advisory lock |
| The API showed 0 objects | it cached an empty catalogue before the first ingest | an empty catalogue is retried every 10 s |
| Grafana's map said "API key required" | the default basemap provider | OpenStreetMap tiles |
| The image build took 15+ minutes | the network delivered about 160 KB/s | `BASE=` builds on an existing Python image |
