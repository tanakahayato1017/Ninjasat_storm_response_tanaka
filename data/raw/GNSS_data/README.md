# NinjaSat raw GNSS telemetry — pending release

The raw onboard GNSS position fixes used to derive the orbital altitude
series in this study are collected and managed by the NinjaSat team at
RIKEN (Tamagawa–Enoto group). Public release of this folder is pending the
team's confirmation of the specific date range that can be shared.

The local repository staging area contains the three monthly inputs used by
the migrated pipeline:

- `conbined_gps_2405_ext.csv`
- `conbined_gps_2407_ext.csv`
- `conbined_gps_2410_ext.csv`

The misspelling `conbined` is retained because it is the filename used by the
verified notebook workflow. These monthly files have no surviving standalone
producer. They were assembled from daily onboard telemetry by concatenation,
retaining records with `gpsFixQual != 0`, and filtering to
`400 < altitude < 600 km` (the operations in `csv_connect.ipynb`, cells 0 and
3). They are therefore provided as-is once release is approved.

**In the meantime**, the derived quantities computed from this data (the
altitude residual $\Delta a$ and its smoothed rate $\dot{\widetilde{\Delta a}}$
at 1-h and 3-h cadence, for the intervals analysed in the paper) are
available in [`../derived/`](../derived/) and [`../figure_data/`](../figure_data/)
— these are sufficient to reproduce every reported number, table, and
figure without needing the raw GNSS fixes.

Do not include these monthly files in a public release until the NinjaSat team
confirms the releasable date range. The raw telemetry is available from the
NinjaSat team at RIKEN on reasonable request (see the manuscript's Data
availability section).
