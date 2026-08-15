# TU Delft thermospheric-density products

This folder vendors the six May/October 2024 files consumed by
`src/08_swarm_density_validation.py` (30.9 MB total): Swarm A/B POD density
v02 and GRACE-FO C accelerometer density v02c.

- Source: http://thermosphere.tudelft.nl/
- Retrieved: 2026-07-02
- Files: `SA_DNS_POD_2024_{05,10}_v02.zip`,
  `SB_DNS_POD_2024_{05,10}_v02.zip`, and
  `GC_DNS_ACC_2024_{05,10}_v02c.zip`

These are public external products. The generated three-hour table is kept
in this same external-data folder because it retains the source products'
public `alt_mean` field; it is not a NinjaSat trajectory.
