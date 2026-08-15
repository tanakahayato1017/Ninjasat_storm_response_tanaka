# NinjaSat storm-time orbital response — analysis code and data

Analysis code and derived data supporting:

> Tanaka, H. et al., "Hour-scale response of low-Earth-orbit satellite
> altitude to extreme geomagnetic storms revealed by the NinjaSat CubeSat
> during the 2024 solar maximum," *Scientific Reports* (submitted).

This repository reproduces every number, table, and figure reported in the
paper from GNSS-derived orbital altitude and public geomagnetic/solar
indices. It is a companion to, and reuses the GNSS-to-altitude and
SGP4-differencing methodology first published for, the EUV-forcing
companion study's repository:
[`Ninjasat_orbit_analysis_tanaka`](https://github.com/tanakahayato1017/Ninjasat_orbit_analysis_tanaka).

## Pipeline overview

```
NinjaSat/XRISM GNSS fixes (~10-min)         TLE catalogues (Space-Track)
              │                                       │
              └────────────────┬──────────────────────┘
                                ▼
                 src/00_sgp4_gnss_altitude.py
                 GNSS altitude + SGP4 propagation
                                │
                                ▼
                 src/02_orbit_average_fit.py
                 two-orbit NinjaSat fit / three-orbit XRISM mean,
                 Δa(t) = a_SGP4(t) − a_GNSS(t), SG smoothing
                                │
                                ▼
                 src/03_euv_baseline_regression.py
                 quiet-time EUV(140.5 nm) → Δȧ̃ regression,
                 subtracted to give the storm-time component
                                │
              ┌─────────────────┼─────────────────────┐
              ▼                 ▼                      ▼
  src/04_pcn_dst_lag_       src/05_ae_al_          src/06_xrism_detection_
  correlation.py            analysis.py            scaling.py
  (Fisher-z, 95% CI,        (1-min AE/AL,          (2σ upper limit,
   inverse-variance lag)     same procedure)         own-May scaling test)
              │                 │                      │
              └────────┬────────┴──────────┬───────────┘
                        ▼                   ▼
          src/07_joule_heating_budget.py   src/08_swarm_density_
          (Knipp et al. 2004 Q_J)           validation.py
                        │                   (independent density check)
                        └─────────┬─────────┘
                                  ▼
                    figures/fig1_lomb_scargle.py
                    (SGP4-artefact exclusion, Fig. 1 oscillation)
                                  ▼
                          figures/fig1..fig7_*.py
```

## Script ↔ manuscript section map

| Script | Manuscript section | Status |
|---|---|---|
| `src/00_sgp4_gnss_altitude.py` | Methods, NinjaSat/XRISM GNSS processing and TLE propagation | done — verified against five monthly reference files |
| `src/02_orbit_average_fit.py` | Methods, orbit averaging, Eq. 1 ($\Delta a$), and SG smoothing | done — verified for NinjaSat and XRISM |
| `src/03_euv_baseline_regression.py` | Methods, EUV baseline removal | done — published slope_difference reproduced exactly |
| `src/04_pcn_dst_lag_correlation.py` | Results §2.4–2.5, Fig. 5a–d input | done |
| `src/05_ae_al_analysis.py` | Results §2.6, Fig. 4e,f, Fig. 5e,f | done |
| `src/06_xrism_detection_scaling.py` | Results §2.7, Table 1 | done |
| `src/07_joule_heating_budget.py` | Discussion energetics, Eq. 5, Fig. 7 | done |
| `src/08_swarm_density_validation.py` | Methods (own-May-response scaling inputs) | done |
| `figures/fig1ab_method.py` | Fig. 1a,b (reimplemented; original was a manual export) | done |
| `figures/fig1_lomb_scargle.py` | Results §2.1 (Fig. 1 oscillation) | done |
| `figures/fig2b_euv_regression.py` | Fig. 2b (r = 0.767, p = 7.2e-44 reproduced) | done |
| `figures/fig3cd_euv_subtraction.py` | Fig. 3c,d | done |
| `figures/fig4ad_covariation.py` | Fig. 4a–d (+ regenerates the n=84 merged series) | done |
| `figures/fig4_5_ae_al.py` | Figs. 4e,f and 5e,f | done |
| `figures/fig5_pcn_dst_crosscorr.py` | Fig. 5a–d | done |
| `figures/fig6_xrism.py` | Fig. 6 | done |
| `figures/fig7_energetics` (in `src/07`) | Fig. 7 | done |

Not reimplemented: Fig. 2a (`fig2a_euv_timeseries.png`) and Fig. 3a,b
(`fig3a/b_storm_residual_*.pdf`) originate from the EUV companion study's
processing (Fig. 2a) and from `src/02`'s smoothing step plotted at
production time (Fig. 3a,b); Fig. 3a,b can be regenerated from the
`src/02` outputs, and Fig. 2a's provenance is the companion repository.

*(This table is updated as each script is migrated; see the repository's
commit history for progress.)*

## Setup

```bash
conda env create -f environment.yml
conda activate ninjasat-storm-response
```

## Data availability

NinjaSat and XRISM raw GNSS telemetry are staged locally pending release
approval and must not be included in a public archive until that approval;
see [`data/README.md`](data/README.md) and
[`data/raw/GNSS_data/README.md`](data/raw/GNSS_data/README.md). The local TLE
catalogues are likewise pending a Space-Track user-agreement review; see
[`data/external/tle/README.md`](data/external/tle/README.md). All
approved reference products needed to reproduce the reported results are
included under `data/derived/` and `data/figure_data/`. Regenerated
non-release intermediates are written to `data/derived/intermediate/`.

## License

Code (`src/`, `figures/`): MIT, see [`LICENSE`](LICENSE).
Data (`data/`): CC BY 4.0, see [`data/LICENSE`](data/LICENSE).

## Citation

See [`CITATION.cff`](CITATION.cff).
