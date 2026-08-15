"""Figure 5a-d: lagged cross-correlation of the EUV-subtracted altitude
rate with PCN (panels a, b) and |Dst| (panels c, d) for the May and
October 2024 storms.

Reads the lag-scan outputs of src/04_pcn_dst_lag_correlation.py and
reproduces the manuscript's Fig. 5a-d panels: black correlation curve,
green segment marking shifts statistically indistinguishable from the
peak (inside the peak's Fisher-z 95% CI), red peak marker, and blue
dashed lines at the CI bounds on r.

Input:
  data/figure_data/{2405,2410}_{PCN,Dst_index_nT}_correlation_results.csv
  data/figure_data/lag_correlation_summary.csv

Output:
  fig5_{PCN,Dst_index_nT}_crosscorr_{2405,2410}.pdf/png

Run (from repository root):
  python figures/fig5_pcn_dst_crosscorr.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "figure_data"
OUT = ROOT / "figures"

MONTHS = [2405, 2410]
INDICES = ["PCN", "Dst_index_nT"]


def main() -> None:
    summary = pd.read_csv(DATA / "lag_correlation_summary.csv")

    for index in INDICES:
        for month in MONTHS:
            df = pd.read_csv(DATA / f"{month}_{index}_correlation_results.csv")
            s = summary[(summary["index"] == index)
                        & (summary["month"] == month)].iloc[0]

            in_ci = df["In_95_CI"]
            peak_idx = df["Correlation"].idxmax()
            peak_h = df["Shift_Hours"].iloc[peak_idx]
            peak_r = df["Correlation"].iloc[peak_idx]

            plt.figure(figsize=(10, 6))
            plt.plot(df["Shift_Hours"], df["Correlation"], marker="o",
                     label="Correlation", linestyle="-", color="k")
            plt.plot(df["Shift_Hours"][in_ci], df["Correlation"][in_ci],
                     marker="o", label="Inside 95% CI", linestyle="-",
                     color="green")
            plt.axhline(s["r_lower"], color="blue", linestyle="--",
                        label="95% Confidence Interval")
            plt.axhline(s["r_upper"], color="blue", linestyle="--")
            plt.scatter(peak_h, peak_r, color="red",
                        label="Peak Correlation", zorder=5)

            x_in = df["Shift_Hours"][in_ci]
            y_in = df["Correlation"][in_ci]
            plt.fill_between(x_in, y_in, color="green", alpha=0.3,
                             label="Confidence Interval Area")
            plt.plot([x_in.min(), x_in.min()],
                     [0, y_in.loc[x_in.idxmin()]],
                     color="green", linestyle="--", linewidth=1)
            plt.plot([x_in.max(), x_in.max()],
                     [0, y_in.loc[x_in.idxmax()]],
                     color="green", linestyle="--", linewidth=1)

            plt.xticks(np.arange(-20, 20, 2), fontsize=16)
            plt.yticks(fontsize=16)
            plt.axhline(0, color="gray", linestyle="--", linewidth=1)
            plt.xlabel("Time Shift (hours)", fontsize=22)
            plt.ylabel("Correlation Coefficient", fontsize=22)
            plt.legend(fontsize=14)
            plt.grid(True, linestyle="--", linewidth=0.7)
            plt.tight_layout()
            plt.xlim(-20, 20)
            plt.ylim(0, 0.9)

            for ext in ["pdf", "png"]:
                out = OUT / f"fig5_{index}_crosscorr_{month}.{ext}"
                plt.savefig(out, format=ext, bbox_inches="tight",
                            dpi=150 if ext == "png" else None)
                print(f"wrote {out}")
            plt.close()


if __name__ == "__main__":
    main()
