# Data sources

This directory holds derived/intermediate data products (`derived/`,
`figure_data/`) and, once released, raw or near-raw NinjaSat GNSS data
under `raw/` (`raw/GNSS_data/`: position fixes; `raw/energy/`: the
specific-orbital-energy series used for the SGP4-independence check in
Fig. 1 — treated the same as raw telemetry since it is close enough to a
direct restatement of orbital radius. See each folder's own README for
status).

All files under `derived/` and `figure_data/` are licensed CC BY 4.0 (see
[`LICENSE`](LICENSE)) and are sufficient, together with the scripts in
[`../src/`](../src/), to reproduce every number, table, and figure reported
in the paper.

## Primary external data sources

None of the third-party geomagnetic/orbital data below is redistributed here
in raw form except where noted; scripts in `../src/` retrieve or expect it at
the paths documented in each script's docstring.

| Data | Source | Access |
|---|---|---|
| NinjaSat GNSS telemetry | NinjaSat team, RIKEN (Tamagawa–Enoto group) | On reasonable request; see `raw/GNSS_data/README.md` |
| XRISM GNSS telemetry | XRISM Science Team, JAXA/ISAS | Local staging under `raw/xrism_gnss/`; redistribution pending team approval |
| TLE catalogue entries | Space-Track.org / 18th Space Defense Squadron | Local staging under `external/tle/`; redistribution pending user-agreement review |
| PC index (PCN, PCS) | DTU Space | Public, https://pcindex.org |
| Dst index | World Data Center for Geomagnetism, Kyoto (via NASA/OMNI2) | Public, https://wdc.kugi.kyoto-u.ac.jp |
| AE, AL auroral-electrojet indices (1-min) | World Data Center for Geomagnetism, Kyoto | Public, https://wdc.kugi.kyoto-u.ac.jp |
| GOES-18 EXIS EUV/X-ray irradiance | NOAA NCEI | Public |
| Swarm A/B accelerometer-derived density | TU Delft thermosphere data set (Siemes et al. 2016) | Public, http://thermosphere.tudelft.nl/ |
| GRACE-FO POD-derived density | TU Delft thermosphere data set (van den IJssel et al. 2020) | Public, http://thermosphere.tudelft.nl/ |

Vendored public external products:

| Folder | Contents | Provenance |
|---|---|---|
| `external/tudelft/` | Six May/October Swarm A/B and GRACE-FO density ZIPs (v02/v02c) | `external/tudelft/README.md` |
| `external/wdc_kyoto/` | One-minute AE/AL CSVs for both storm intervals | `external/wdc_kyoto/README.md` |
| `external/omni2/` | Hourly 2024 OMNI2 text product and May/October analysis CSVs | `external/omni2/README.md` |
| `external/pcn_1min/` | One-minute PCN CSVs for May and October 2024 | `external/pcn_1min/README.md` |

Local-only, non-public staging products are documented separately under
`raw/GNSS_data/`, `raw/xrism_gnss/`, and `external/tle/`. Their README files
take precedence over the repository-wide data licence until release approval
or licence review is complete.

Retrieval dates and any version/provisional-vs-final notes are recorded in
the docstring of the script that consumes each source.

## `derived/`

Intermediate time series (altitude residual $\Delta a$, smoothed
$\widetilde{\Delta a}$ and its rate $\dot{\widetilde{\Delta a}}$, at 1-h and
3-h cadence, for NinjaSat and XRISM) that feed the correlation and
detection-significance analyses. The regenerable, non-release verification
workspace is `derived/intermediate/`, populated by
`src/00_sgp4_gnss_altitude.py` and `src/02_orbit_average_fit.py`.

The public `derived/energy/dEdt_storm_3h_redacted.csv` retains only the
rate columns needed by `src/08`; the absolute `E_J2_mean`
column was removed. Its local README documents the redaction.
`derived/correlation_inputs/` holds the n=84 correlation grids. The
three-hour TU Delft table remains under
`external/tudelft/` because it retains the external products' public mean
altitude field.

## `figure_data/`

Per-figure/table CSV outputs consumed directly by the scripts in
`../figures/`. Populated by `src/04`–`src/08`.

Note on `xrism_detection_limit.csv`: its four `mean_altitude_km` rows
(503.5/474.3 km for NinjaSat, 569.4/560.8 km for XRISM) are
monthly-mean scalars that are quoted verbatim in the manuscript's
Results section; as values already published in the co-authored paper
they are exempt from the raw-data release restriction, which covers
absolute-altitude *time series*.
