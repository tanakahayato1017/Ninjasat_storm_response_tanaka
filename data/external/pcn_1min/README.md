# One-minute Polar Cap North index

`PCN_2405.csv` and `PCN_2410.csv` contain the one-minute Polar Cap North
(PCN) index used to generate manuscript Fig. 7 and its quoted
Joule-heating values.

- Source: DTU Space, https://pcindex.org
- Coverage: May and October 2024
- Columns: UTC timestamp and PCN in mV/m
- Fill convention: `999` denotes missing data and is converted to NaN
  before the inclusive +/-30-minute averages are calculated.
