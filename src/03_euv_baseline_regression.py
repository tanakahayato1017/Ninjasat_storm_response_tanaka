"""EUV baseline regression and storm-time subtraction (manuscript Methods
"EUV baseline removal", Results Sect. 2.2-2.3; produces the
`estimated_slope` and `slope_difference` series behind Fig. 3c,d and every
storm-time correlation analysis).

Three faithful stages, ported from the original analysis notebooks
(correlations_2.ipynb cell 3 and regression_test.ipynb cells 3+5; the
numerical procedure is unchanged, including its quirks, so that the
outputs reproduce the published series exactly):

Stage 1 - storm-month EUV merge. For each storm month (2024-05, 2024-10),
  average the 1-min GOES-18 EXIS EUV irradiance channels in a +/-1.57-h
  box centred on each altitude-rate epoch of the Savitzky-Golay-smoothed
  NinjaSat series, drop epochs whose box contains any NaN column, and
  attach the box means to the altitude series with a nearest-timestamp
  asof merge. Output: the storm-month merged series
  ({month}_irr_1405_correlation_thesis.csv equivalent).

Stage 2 - quiet-month regression. For each geomagnetically quiet month
  July-September 2024 and each of the seven EXIS channels, scan integer
  row shifts -35..+2 of the altitude-rate `slope` against the channel
  irradiance, pick the maximum-correlation shift, apply it, and fit a
  linear regression slope/intercept by np.polyfit. Output:
  regression_results.csv (Month, Value_Index, Slope, Intercept,
  Peak_Correlation, Peak_Shift_index, Peak_Shift_Hours).

Stage 3 - storm application. Average the 140.5-nm (irr_1405) regression
  coefficients over 2024-07..09, compute
  estimated_slope = irr_1405 * mean_slope + mean_intercept, delay it by
  the fixed 16-sample shift used in the published analysis (~50 h,
  the rounded mean quiet-month peak shift; hard-coded exactly as in the
  original), and form slope_difference = slope - estimated_slope.
  Output: regression_results_subtract_slope_{month}_thesis.csv
  equivalent - the direct input of the Fig. 4 covariation panels and of
  src/04's lag correlations.

Provenance note. Stage 3 reads the coefficients from
  data/derived/correlation_inputs/regression_results_published.csv - the
  archived coefficient table from the original analysis run - so that the
  published slope_difference series is reproduced EXACTLY. Stage 2
  regenerates the same table from the quiet-month inputs vendored here
  and reports any per-coefficient drift against the published copy: the
  quiet-month merged CSVs were re-generated after the archived
  coefficient table was produced, which shifts the recomputed
  coefficients by <~1% and the final slope_difference by <4e-6 km/sample
  (<0.05% of signal) - negligible for every reported number, but not
  bit-identical. All peak-lag shift indices are identical between the
  two vintages.

Inputs:
  data/external/goes_euv/goes_euv_{2405,2410}_combined.csv
      1-min GOES-18 EXIS EUV irradiance (public, NOAA NCEI).
  data/derived/ninjasat/gps_tle_average_altitude_fitting_all_smoothed_
      {2405,2410}_ext_thesis.csv   (SG-smoothed altitude series, src/02)
  data/derived/correlation_inputs/quiet_euv/{2407,2408,2409}_irr_{...}
      _correlation.csv             (quiet-month merged series; the same
      merge as Stage 1 applied to the quiet months - also published in
      the EUV companion study's repository)

Outputs:
  data/derived/correlation_inputs/{2405,2410}_irr_1405_correlation_thesis.csv
  data/figure_data/regression_results.csv
  data/derived/correlation_inputs/regression_subtract_slope_{2405,2410}_thesis.csv

Run (from repository root):
  python src/03_euv_baseline_regression.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
GOES = ROOT / "data" / "external" / "goes_euv" / "goes_euv_{}_combined.csv"
SMOOTHED = (ROOT / "data" / "derived" / "ninjasat" /
            "gps_tle_average_altitude_fitting_all_smoothed_{}_ext_thesis.csv")
QUIET_DIR = ROOT / "data" / "derived" / "correlation_inputs" / "quiet_euv"
OUT_MERGE = (ROOT / "data" / "derived" / "correlation_inputs" /
             "{}_irr_1405_correlation_thesis.csv")
OUT_REG = ROOT / "data" / "figure_data" / "regression_results.csv"
PUBLISHED_REG = (ROOT / "data" / "derived" / "correlation_inputs" /
                 "regression_results_published.csv")
OUT_SUB = (ROOT / "data" / "derived" / "correlation_inputs" /
           "regression_subtract_slope_{}_thesis.csv")

STORM_MONTHS = [2405, 2410]
QUIET_MONTHS = [2407, 2408, 2409]
EUV_CHANNELS = ["irr_256", "irr_284", "irr_304",
                "irr_1175", "irr_1216", "irr_1335", "irr_1405"]
BOX_HALF_H = 1.57
SHIFT_SCAN = range(-35, 3)
FIXED_STORM_SHIFT = 16  # samples (~50 h); hard-coded in the published analysis


def stage1_storm_merge(month: int) -> pd.DataFrame:
    """correlations_2.ipynb cell 3: +/-1.57-h box average of the EUV
    channels onto the smoothed altitude-rate epochs, asof-merged."""
    euv = pd.read_csv(str(GOES).format(month))
    euv["Time"] = pd.to_datetime(euv["Time"], errors="coerce")
    alt = pd.read_csv(str(SMOOTHED).format(month))
    alt["timestamp"] = pd.to_datetime(alt["timestamp"], format="mixed",
                                      errors="coerce").dt.tz_localize(None)

    aggregated = []
    for t in alt["timestamp"]:
        win = euv[(euv["Time"] >= t - pd.Timedelta(hours=BOX_HALF_H))
                  & (euv["Time"] <= t + pd.Timedelta(hours=BOX_HALF_H))]
        mean_values = win.mean(numeric_only=True)
        mean_values["timestamp"] = t
        if not mean_values.isnull().values.any():
            aggregated.append(mean_values)
    agg = pd.DataFrame(aggregated).rename(columns={"timestamp": "Time"})

    merged = pd.merge_asof(alt.sort_values("timestamp"),
                           agg.sort_values("Time"),
                           left_on="timestamp", right_on="Time",
                           direction="nearest")
    merged = merged.drop(columns=["Time"])
    out = Path(str(OUT_MERGE).format(month))
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out, index=False)
    print(f"[stage1 {month}] wrote {out} ({len(merged)} rows)")
    return merged


def stage2_quiet_regression() -> pd.DataFrame:
    """regression_test.ipynb cell 3: per quiet month and channel, peak-lag
    shift scan then linear regression of shifted slope on irradiance."""
    rows = []
    for value_index in EUV_CHANNELS:
        for month in QUIET_MONTHS:
            data = pd.read_csv(QUIET_DIR / f"{month}_{value_index}_correlation.csv")
            data["timestamp"] = pd.to_datetime(data["timestamp"])
            steps = [(data["timestamp"].iloc[i] - data["timestamp"].iloc[i - 1])
                     .total_seconds() / 3600 for i in range(1, 10)]
            time_step = sum(steps) / len(steps)
            data = data[data["timestamp"].dt.month == month - 2400]

            corr_rows = []
            for shift in SHIFT_SCAN:
                shifted = data["slope"].shift(shift)
                valid = pd.DataFrame({"s": shifted,
                                      "x": data[value_index]}).dropna()
                corr_rows.append({"Shift": shift,
                                  "Correlation": valid.corr().iloc[0, 1]})
            cdf = pd.DataFrame(corr_rows)
            cdf["Shift_Hours"] = cdf["Shift"] * time_step

            peak_idx = cdf["Correlation"].idxmax()
            peak_shift = int(cdf["Shift"].iloc[peak_idx])
            data["slope_shifted"] = data["slope"].shift(peak_shift)
            new_data = data[["slope_shifted", value_index]].dropna()
            slope, intercept = np.polyfit(new_data[value_index],
                                          new_data["slope_shifted"], 1)
            rows.append({
                "Month": month, "Value_Index": value_index,
                "Slope": slope, "Intercept": intercept,
                "Peak_Correlation": cdf["Correlation"].iloc[peak_idx],
                "Peak_Shift_index": peak_shift,
                "Peak_Shift_Hours": cdf["Shift_Hours"].iloc[peak_idx],
            })
            print(f"[stage2 {month} {value_index}] peak r = "
                  f"{cdf['Correlation'].iloc[peak_idx]:.4f} at shift "
                  f"{peak_shift} ({cdf['Shift_Hours'].iloc[peak_idx]:.1f} h)")

    results = pd.DataFrame(rows)
    OUT_REG.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUT_REG, index=False)
    print(f"[stage2] wrote {OUT_REG} ({len(results)} rows)")

    published = pd.read_csv(PUBLISHED_REG)
    merged = published.merge(results, on=["Month", "Value_Index"],
                             suffixes=("_pub", "_new"))
    drift = (merged["Slope_new"] - merged["Slope_pub"]).abs() / merged["Slope_pub"].abs()
    shifts_equal = (merged["Peak_Shift_index_pub"]
                    == merged["Peak_Shift_index_new"]).all()
    print(f"[stage2] drift vs published coefficients: max relative slope "
          f"drift = {drift.max():.2e}; all peak-shift indices identical = "
          f"{shifts_equal} (see docstring provenance note)")
    return results


def stage3_storm_subtract() -> None:
    """regression_test.ipynb cell 5: apply the 2024-07..09-averaged
    140.5-nm regression to the storm months with the fixed 16-sample
    delay, and form slope_difference. Uses the PUBLISHED coefficient
    table for exact reproduction (see docstring provenance note)."""
    results = pd.read_csv(PUBLISHED_REG)
    target = results[(results["Value_Index"] == "irr_1405")
                     & (results["Month"].isin(QUIET_MONTHS))]
    mean_slope = target["Slope"].mean()
    mean_intercept = target["Intercept"].mean()
    mean_shift = int(round(target["Peak_Shift_index"].mean()))
    print(f"[stage3] mean regression over {QUIET_MONTHS}: "
          f"slope = {mean_slope:.4e}, intercept = {mean_intercept:.4e}, "
          f"mean peak shift = {mean_shift} (published analysis uses the "
          f"fixed value {FIXED_STORM_SHIFT})")

    for month in STORM_MONTHS:
        merged = pd.read_csv(str(OUT_MERGE).format(month))
        merged["timestamp"] = pd.to_datetime(merged["timestamp"])
        merged["estimated_slope"] = (merged["irr_1405"] * mean_slope
                                     + mean_intercept)
        merged["estimated_slope"] = merged["estimated_slope"].shift(FIXED_STORM_SHIFT)
        merged["slope_difference"] = merged["slope"] - merged["estimated_slope"]
        out = Path(str(OUT_SUB).format(month))
        merged.to_csv(out, index=False)
        print(f"[stage3 {month}] wrote {out} ({len(merged)} rows)")


def main() -> None:
    for month in STORM_MONTHS:
        stage1_storm_merge(month)
    stage2_quiet_regression()
    stage3_storm_subtract()


if __name__ == "__main__":
    main()
