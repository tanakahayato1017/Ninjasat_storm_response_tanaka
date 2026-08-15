# XRISM GNSS/TLE altitude series — pending release

Same status as [`../ninjasat/`](../ninjasat/): these files carry the
absolute GNSS-derived orbital altitude (`average_altitude_gps`) alongside
the SGP4-GNSS difference, and are pending confirmation from the XRISM
Science Team (JAXA/ISAS) before release.

The verified local staging copies are present for migration testing, but they
must be excluded from a public archive until release is approved. The monthly
GNSS inputs from which they can be regenerated are documented separately in
[`../xrism_gnss/README.md`](../xrism_gnss/README.md).

`src/06_xrism_detection_scaling.py` consumes these files. Its published output
(detection significances, 2-sigma upper limit,
own-May-response scaling predictions — matching Table 1 of the manuscript)
is already available in
[`../../figure_data/xrism_detection_limit.csv`](../../figure_data/xrism_detection_limit.csv).
