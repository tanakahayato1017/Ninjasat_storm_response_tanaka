# XRISM raw monthly GNSS data — pending release

This local staging folder contains the monthly XRISM orbit products read by
`src/00_sgp4_gnss_altitude.py`:

- `xrism_gps_2405.csv`
- `xrism_gps_2410.csv`

They include time-tagged Cartesian position, geodetic position, velocity, and
osculating-orbit fields supplied by the XRISM Science Team (JAXA/ISAS). Public
redistribution is pending team approval. Do not include either CSV in the
public Zenodo archive until approval is recorded.

The verified orbit-averaged residual products needed for the manuscript's
reported results remain under `data/raw/xrism/` and `data/derived/xrism/`.
