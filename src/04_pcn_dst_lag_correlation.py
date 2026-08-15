"""Lagged cross-correlation of the EUV-subtracted altitude rate against the
PCN and |Dst| indices, with Fisher-z 95% confidence intervals and the
inverse-variance-weighted combined PCN lag (manuscript Results Sect. 2.4
and 2.5, Methods "Statistical analysis", Fig. 5a-d input data).

Method (faithful port of the analysis notebook that produced the
manuscript numbers; the numerical procedure is unchanged):

1. Lag scan. For integer row shifts k = -8 .. +7, shift the
   EUV-subtracted altitude rate (`slope_difference`) by k rows relative to
   the index series and compute the Pearson correlation r(k) over the
   overlapping samples. The time value of each shift is k times the FIRST
   sampling interval of that month's series (3.167 h for both months —
   exactly as in the archived analysis whose outputs this script
   reproduces byte-for-byte; note the October series' median cadence is
   3.121 h, but the original analysis used the first interval), so
   tau = k * dt. Negative tau means the geomagnetic index leads the
   altitude response.

2. Peak and Fisher-z 95% CI. The peak lag tau_peak is the shift
   maximising r. Because the sampling distribution of r is skewed, the
   95% CI on the peak correlation uses the Fisher transformation
       z = arctanh(r),  SE(z) = 1/sqrt(n - 3),  n = 84 paired samples,
   with the interval z_peak +/- 1.960 SE(z) mapped back through tanh.

3. Lag uncertainty. Shifts whose r(k) lies inside the peak's 95% CI are
   statistically indistinguishable from the peak; the min/max of their
   time values define an asymmetric lag uncertainty
       tau_peak (+ (ci_max - tau_peak) / - (tau_peak - ci_min)),
   with symmetric half-width sigma_tau = (ci_max - ci_min)/2 used for the
   combination step.

4. Combined PCN lag. The two events' PCN lags are combined by
   inverse-variance weighting, w_i = 1/sigma_i^2:
       tau_combined = sum(w_i tau_i)/sum(w_i),
       sigma_combined = 1/sqrt(sum(w_i)).

Input (data/derived/correlation_inputs/, n = 84 rows each; the
`slope_difference` column is the EUV-baseline-subtracted altitude rate and
the index columns are public geomagnetic data — freely redistributable):
  {2405,2410}_PCN_correlation_thesis.csv          (column: PCN)
  {2405,2410}_Dst_index_nT_correlation_thesis.csv (column: Dst_index_nT = |Dst|)

Output:
  data/figure_data/{month}_{index}_correlation_results.csv
      columns: Shift, Correlation, Shift_Hours, In_95_CI
      (identical format to the archived notebook outputs, so the two can
      be diffed directly)
  data/figure_data/lag_correlation_summary.csv
      per index/month: n, peak r, r CI, peak lag, CI lag range,
      asymmetric/symmetric uncertainties; plus the combined PCN row

Reference values this must reproduce (manuscript):
  PCN May  r = 0.802 at tau = 0 h;   PCN Oct r = 0.852 at tau = -3.12 h
  Dst May  r = 0.805 at tau = +3.17 h; Dst Oct r = 0.825 at tau = +3.12 h
  Combined PCN lag tau = -2.19 +/- 2.64 h (October weight ~69%)

Run (from repository root):
  python src/04_pcn_dst_lag_correlation.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
IN_DIR = ROOT / "data" / "derived" / "correlation_inputs"
OUT_DIR = ROOT / "data" / "figure_data"

MONTHS = [2405, 2410]
INDICES = ["PCN", "Dst_index_nT"]
SHIFTS = range(-8, 8)
Z_CRIT = norm.ppf(0.975)


def lag_scan(month: int, index: str) -> dict:
    """Cells 11+13 of the original notebook: shift scan, peak, Fisher-z CI,
    and the CI-based lag uncertainty, for one (month, index) pair."""
    data = pd.read_csv(IN_DIR / f"{month}_{index}_correlation_thesis.csv")
    data["timestamp"] = pd.to_datetime(data["timestamp"])
    time_step = (data["timestamp"].iloc[1]
                 - data["timestamp"].iloc[0]).total_seconds() / 3600.0
    n = len(data)

    rows = []
    for shift in SHIFTS:
        shifted = data["slope_difference"].shift(shift)
        valid = pd.DataFrame({"s": shifted, "x": data[index]}).dropna()
        rows.append({"Shift": shift, "Correlation": valid.corr().iloc[0, 1]})
    df = pd.DataFrame(rows)
    df["Shift_Hours"] = df["Shift"] * time_step

    peak_idx = df["Correlation"].idxmax()
    peak_r = df["Correlation"].iloc[peak_idx]
    peak_h = df["Shift_Hours"].iloc[peak_idx]

    z_peak = 0.5 * np.log((1 + peak_r) / (1 - peak_r))
    se = np.sqrt(1.0 / (n - 3))
    r_lower = np.tanh(z_peak - Z_CRIT * se)
    r_upper = np.tanh(z_peak + Z_CRIT * se)
    df["In_95_CI"] = df["Correlation"].between(r_lower, r_upper)

    in_ci = df[df["In_95_CI"]]
    ci_min_h = in_ci["Shift_Hours"].min()
    ci_max_h = in_ci["Shift_Hours"].max()

    out = OUT_DIR / f"{month}_{index}_correlation_results.csv"
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"wrote {out}")

    return {
        "index": index, "month": month, "n": n, "time_step_h": time_step,
        "peak_r": peak_r, "r_lower": r_lower, "r_upper": r_upper,
        "peak_lag_h": peak_h, "ci_min_h": ci_min_h, "ci_max_h": ci_max_h,
        "sigma_minus_h": peak_h - ci_min_h, "sigma_plus_h": ci_max_h - peak_h,
        "half_width_h": (ci_max_h - ci_min_h) / 2.0,
    }


def main() -> None:
    summaries = []
    for index in INDICES:
        print("=" * 65)
        print(f"{index}: peak lag and 95%-CI-based lag uncertainty")
        print("=" * 65)
        for month in MONTHS:
            s = lag_scan(month, index)
            summaries.append(s)
            print(f"\n[{month}]  n = {s['n']}, dt = {s['time_step_h']:.3f} h")
            print(f"  peak r            : {s['peak_r']:.4f}  "
                  f"(95% CI [{s['r_lower']:.4f}, {s['r_upper']:.4f}])")
            print(f"  peak lag          : {s['peak_lag_h']:+.3f} h")
            print(f"  shifts inside CI  : {s['ci_min_h']:+.3f} .. {s['ci_max_h']:+.3f} h")
            print(f"  lag uncertainty   : {s['peak_lag_h']:+.3f} "
                  f"+{s['sigma_plus_h']:.3f} / -{s['sigma_minus_h']:.3f} h "
                  f"(symmetric +/- {s['half_width_h']:.3f} h)")
        print()

    # ---- inverse-variance-weighted combination of the two PCN lags ----
    pcn = [s for s in summaries if s["index"] == "PCN"]
    lags = np.array([s["peak_lag_h"] for s in pcn])
    sigmas = np.array([s["half_width_h"] for s in pcn])
    weights = 1.0 / sigmas ** 2
    tau_combined = float(np.sum(weights * lags) / np.sum(weights))
    sigma_combined = float(1.0 / np.sqrt(np.sum(weights)))

    print("=" * 65)
    print("PCN: inverse-variance-weighted combined lag")
    print("=" * 65)
    for s, w in zip(pcn, weights / np.sum(weights) * 100):
        print(f"  {s['month']}: tau = {s['peak_lag_h']:+.3f} h, "
              f"sigma = +/-{s['half_width_h']:.3f} h, weight = {w:.1f}%")
    print(f"\n  tau_combined   = {tau_combined:+.3f} h")
    print(f"  sigma_combined = +/-{sigma_combined:.3f} h")
    print(f"  result: tau = {tau_combined:.2f} +/- {sigma_combined:.2f} h")

    summary = pd.DataFrame(summaries)
    combined = pd.DataFrame([{
        "index": "PCN_combined", "month": 0, "n": sum(s["n"] for s in pcn),
        "peak_lag_h": tau_combined, "half_width_h": sigma_combined,
    }])
    out = OUT_DIR / "lag_correlation_summary.csv"
    pd.concat([summary, combined], ignore_index=True).to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
