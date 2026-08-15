# TLE catalogues — local staging only

The two catalogues in this folder were retrieved from
[Space-Track.org](https://www.space-track.org/) for the NinjaSat and XRISM
orbit reconstructions:

- `ninjasat_tle_per_day.csv`
- `xrism_tle_each_day.csv`

They are required locally by `src/00_sgp4_gnss_altitude.py`. They are **not to
be redistributed in the public release** pending review of the Space-Track
user agreement and confirmation of an allowed redistribution route. Exclude
both CSV files when building a public archive unless that review explicitly
clears them.

The script uses one zero-based catalogue row per month. The verified mappings
and corresponding catalogue epochs are documented directly in the script.
