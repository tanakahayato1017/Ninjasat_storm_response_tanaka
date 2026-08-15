"""Figure 2b: quiet-period (July 2024) linear regression of the
time-shifted altitude rate against the GOES-18 EXIS 140.5-nm irradiance
(manuscript Results Sect. 2.2, Methods "EUV baseline removal").

Faithful port of the original notebook cell (regression_test.ipynb
cell 1, run with month_index=2407, value_index='irr_1405' - the
submitted fig2b_euv_regression.pdf is a byte-copy of that cell's
regression_2407_irr_1405.pdf output): scan integer row shifts -35..+2 of
the altitude-rate `slope` against the irradiance, pick the
maximum-correlation shift, apply it, and scatter-plot the shifted rate
against irradiance with the least-squares line.

Expected values (manuscript): r = 0.767, p = 7.2e-44, peak shift -17
samples (~-54 h). The regression slope/intercept correspond to the
2407/irr_1405 row of data/figure_data/regression_results.csv (the
regenerated table; see src/03's provenance note on the ~1% vintage
drift against the archived published table).

Unit correction (author-approved, following the Fig. 6 precedent): the
submitted panel's y axis read "[km/s]" while the rate is km per ~3.16-h
sample step; the DISPLAYED rate and fit line (and the coefficient values
shown in the in-figure text box) are divided by the month's median
sampling step and labelled km/h. The regression itself is still fit in
the original per-step units so its coefficients remain identical to
regression_results.csv; r and p are scale-invariant and unchanged.

Input:
  data/derived/correlation_inputs/quiet_euv/2407_irr_1405_correlation.csv

Output:
  fig2b_euv_regression.pdf/png

Run (from repository root):
  python figures/fig2b_euv_regression.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.ticker import ScalarFormatter  # noqa: E402
from scipy.stats import pearsonr  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
IN = (ROOT / "data" / "derived" / "correlation_inputs" / "quiet_euv" /
      "2407_irr_1405_correlation.csv")
OUT = ROOT / "figures"

MONTH = 2407
VALUE_INDEX = "irr_1405"
SHIFT_SCAN = range(-35, 3)


def main() -> None:
    data = pd.read_csv(IN)
    data["timestamp"] = pd.to_datetime(data["timestamp"])
    steps = [(data["timestamp"].iloc[i] - data["timestamp"].iloc[i - 1])
             .total_seconds() / 3600 for i in range(1, 10)]
    time_step = sum(steps) / len(steps)
    data = data[data["timestamp"].dt.month == MONTH - 2400]

    corr = []
    for shift in SHIFT_SCAN:
        shifted = data["slope"].shift(shift)
        valid = pd.DataFrame({"s": shifted, "x": data[VALUE_INDEX]}).dropna()
        corr.append({"Shift": shift, "Correlation": valid.corr().iloc[0, 1]})
    cdf = pd.DataFrame(corr)
    peak_idx = cdf["Correlation"].idxmax()
    peak_shift = int(cdf["Shift"].iloc[peak_idx])
    print(f"peak correlation {cdf['Correlation'].iloc[peak_idx]:.4f} at "
          f"shift {peak_shift} ({peak_shift * time_step:.1f} h)")

    data["slope_shifted"] = data["slope"].shift(peak_shift)
    new_data = data[["timestamp", "slope_shifted", VALUE_INDEX]].dropna()

    reg_slope, intercept = np.polyfit(new_data[VALUE_INDEX],
                                      new_data["slope_shifted"], 1)
    r_value, p_value = pearsonr(new_data[VALUE_INDEX],
                                new_data["slope_shifted"])
    print(f"slope = {reg_slope:.4e}, intercept = {intercept:.4e} "
          f"(per-step units, = regression_results.csv), "
          f"r = {r_value:.4f}, p = {p_value:.2e}, n = {len(new_data)}")

    # Display in km/h (see docstring unit correction); the fit itself is
    # in per-step units so the coefficients match regression_results.csv.
    disp_y = new_data["slope_shifted"] / time_step
    disp_slope = reg_slope / time_step
    disp_intercept = intercept / time_step

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.scatter(new_data[VALUE_INDEX], disp_y, s=20, label="Data")
    ax.plot(new_data[VALUE_INDEX],
            disp_slope * new_data[VALUE_INDEX] + disp_intercept,
            color="red", label="Linear fitting line")
    ax.set_ylabel("Time shifted time rates of change\n"
                  r"$\widetilde{\dot{\Delta a}}(t)$ [km/h]", fontsize=25)
    ax.set_xlabel("EUV 140.5 nm irradiance [W/m$^2$]", fontsize=25)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((0, 0))
    ax.xaxis.set_major_formatter(formatter)
    ax.tick_params(axis="x", rotation=45, labelsize=20)
    ax.tick_params(axis="y", labelsize=20)
    ax.legend(fontsize=18)
    ax.grid()

    stats_text = (f"$r$ = {r_value:.3f}  ($p$ = {p_value:.2e})\n"
                  f"Slope = {disp_slope:.3e}\n"
                  f"Intercept = {disp_intercept:.3e}")
    ax.text(0.97, 0.05, stats_text, transform=ax.transAxes, fontsize=16,
            verticalalignment="bottom", horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      alpha=0.8))
    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = OUT / f"fig2b_euv_regression.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight",
                    dpi=150 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
