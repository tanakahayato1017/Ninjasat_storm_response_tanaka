"""Figure 3a/3b: storm-time smoothed altitude residual and its rate of
change for the May and October 2024 storms.

Port of the original plotting step (orbit_gps_tle_3.ipynb cell 17's plot
section, which produced the submitted
gps_tle_altitude_difference_smoothed_{2405,2410}.pdf = fig3a/b): the
Savitzky-Golay-smoothed residual Delta a-tilde (green, left axis) and
its derivative (orange, right axis), storm windows 2024-05-07..05-15 and
2024-10-08..10-15.

Unit correction (author-approved, following the Fig. 6 precedent): the
submitted panels' right axis read "[km/s]" while the plotted derivative
was savgol_filter(..., deriv=1) in km per ~3.1-3.2-h sample step. This
script divides the derivative by each month's median sampling step and
labels the axis km/h; curve shapes are unchanged (a per-month rescale).

Input:
  data/derived/ninjasat/gps_tle_average_altitude_fitting_all_smoothed_
      {2405,2410}_ext_thesis.csv   (src/02 output; N=9 smoothing)

Output:
  fig3a_storm_residual_2405.pdf/png
  fig3b_storm_residual_2410.pdf/png

Run (from repository root):
  python figures/fig3ab_storm_residual.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
IN = (ROOT / "data" / "derived" / "ninjasat" /
      "gps_tle_average_altitude_fitting_all_smoothed_{}_ext_thesis.csv")
OUT = ROOT / "figures"

WINDOWS = {2405: (pd.Timestamp("2024-05-07"), pd.Timestamp("2024-05-15")),
           2410: (pd.Timestamp("2024-10-08"), pd.Timestamp("2024-10-15"))}
PANEL = {2405: "fig3a_storm_residual_2405", 2410: "fig3b_storm_residual_2410"}


def main() -> None:
    for month, (lo, hi) in WINDOWS.items():
        df = pd.read_csv(str(IN).format(month))
        t = (pd.to_datetime(df["timestamp"], format="mixed", utc=True)
             .dt.tz_localize(None))
        df = df.assign(t=t).sort_values("t").reset_index(drop=True)
        step_h = np.median(np.diff(df["t"].to_numpy())
                           .astype("timedelta64[s]").astype(float)) / 3600.0
        w = df[(df["t"] >= lo) & (df["t"] <= hi)]
        rate = pd.to_numeric(w["slope"]) / step_h  # km/step -> km/h
        print(f"[{month}] step={step_h:.3f} h, n={len(w)}, "
              f"rate range [{rate.min():+.5f}, {rate.max():+.5f}] km/h")

        fig, ax2 = plt.subplots(figsize=(14, 7))
        color_altitude = "green"
        ax2.set_xlabel("Time [UTC]", fontsize=25)
        ax2.set_ylabel("Altitude Difference [km]", color=color_altitude,
                       fontsize=25)
        altitude_line, = ax2.plot(w["t"], w["altitude_diff_smoothed"],
                                  color=color_altitude,
                                  label="Altitude Difference [km]")
        ax2.tick_params(axis="y", labelcolor=color_altitude, labelsize=20)
        ax2.tick_params(axis="x", rotation=45, labelsize=20)

        ax1 = ax2.twinx()
        color_slope = "orange"
        ax1.set_ylabel("Time rates of change\n"
                       r"$\widetilde{\dot{\Delta a}}(t)$ [km/h]",
                       color=color_slope, fontsize=22)
        slope_line, = ax1.plot(w["t"], rate, color=color_slope,
                               label=r"$\widetilde{\dot{\Delta a}}(t)$ [km/h]")
        ax1.tick_params(axis="y", labelcolor=color_slope, labelsize=20)

        ax2.yaxis.grid(True, linestyle="--", color="gray", alpha=0.7)
        ax2.set_xlim(lo, hi)
        ax2.xaxis.set_major_locator(mdates.DayLocator(interval=1))
        ax2.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))

        lines = [altitude_line, slope_line]
        ax2.legend(lines, [ln.get_label() for ln in lines],
                   loc="upper left", fontsize=16)
        fig.tight_layout()
        fig.subplots_adjust(right=0.85, top=0.95)

        for ext in ["pdf", "png"]:
            out = OUT / f"{PANEL[month]}.{ext}"
            fig.savefig(out, format=ext, bbox_inches="tight",
                        dpi=150 if ext == "png" else None)
            print(f"  wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()
