# OMNI2 hourly data

`omni2_2024.dat` is the public 2024 low-resolution hourly OMNI2 product
used for Dst, PC, AE, AL, and AU cross-checks.

`omni2_2405_with_utc.csv` and `omni2_2410_with_utc.csv` are the hourly
May and October analysis tables used by `src/07_joule_heating_budget.py`.
They retain the signed `Dst_index_nT` field and an explicit `UTC_Time`
column, matching the inputs used to generate manuscript Fig. 7.

- Source: https://omniweb.gsfc.nasa.gov/
- Retrieved: 2026-07-06
- Coverage: 2024-01-01 00:00 through 2024-12-31 23:00 UTC

The consuming scripts document the zero-based columns and fill values used
from the standard 55-field text format.
