# NinjaSat GNSS/TLE altitude series — pending release

Files here (`gps_tle_average_altitude_fitting_all_severe_*_monthly_ext.csv`)
carry, alongside the SGP4-GNSS altitude difference (`altitude_diff`), the
**absolute** GNSS-derived orbital altitude (`gps_average_altitude`) at
~3-h cadence. Because the absolute altitude time series is close enough to
a direct restatement of the satellite's trajectory, this folder is treated
as pending the NinjaSat team's release approval, on the same terms as
[`../GNSS_data/`](../GNSS_data/) and [`../energy/`](../energy/).

Scripts that need this folder (e.g. `figures/fig1_lomb_scargle.py`,
`src/06_xrism_detection_scaling.py`) will run once these files are added.
The final derived/summary outputs of those scripts (rates, significances,
upper limits — not absolute altitude) are already published in
[`../../figure_data/`](../../figure_data/), matching the values reported
in the manuscript.
