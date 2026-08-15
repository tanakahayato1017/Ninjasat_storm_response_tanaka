"""Figure 1a/1b: representative one-month (July 2024) comparison of the
NinjaSat GNSS-derived and SGP4-propagated orbital altitudes, the raw
residual Delta a, and its Savitzky-Golay-smoothed series and derivative.

Reimplementation note. The submitted fig1a/fig1b PDFs were one-off manual
exports (slide pages titled "Figure2-2"/"Figure2-3") whose generating
code did not survive; this script recreates them from the archived data.
Faithfulness was verified against the submitted panels: fig1a's
altitude/residual curves come directly from the two-orbit-averaged
monthly series, and fig1b's smoothed/derivative curves reproduce the
submitted value ranges (smoothed 0.00..0.235 km; derivative
-0.0054..+0.0060 km/step) only with a 13-point Savitzky-Golay window --
the January-2025 processing vintage -- whereas the manuscript's Methods
section describes the N=9 window used for every storm-month analysis.
The submitted Fig. 1b therefore shows a W13-smoothed curve. This script
defaults to the published figure (window 13) and offers --window 9 for
the Methods-consistent variant; the difference is visual only (July is
an illustrative quiet month; no reported number derives from Fig. 1b).

Unit correction (author-approved, following the Fig. 6 precedent). The
submitted fig1b right axis was labelled "km/s" while the plotted
derivative was savgol_filter(..., deriv=1) in km per ~3.16-h sample
step. This script divides the derivative by the month's median sampling
step and labels the axis km/h; the curve shape is unchanged. The
13-point smoothing window of the published Fig. 1b is retained by the
authors' decision (it is the display window for this month-scale
overview; all storm-time analyses use N=9, now stated explicitly in the
manuscript's Fig. 1 legend).

Input (data/raw/ninjasat/ -- pending NinjaSat team release, staged
locally; see that folder's README):
  gps_tle_average_altitude_fitting_all_severe_2407_monthly_ext.csv

Output:
  fig1a_gnss_vs_sgp4.pdf/png
  fig1b_residual_derivative.pdf/png

Run (from repository root):
  python figures/fig1ab_method.py [--window 13]
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import savgol_filter  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
IN = (ROOT / "data" / "raw" / "ninjasat" /
      "gps_tle_average_altitude_fitting_all_severe_2407_monthly_ext.csv")
OUT = ROOT / "figures"

JULY = (pd.Timestamp("2024-07-01"), pd.Timestamp("2024-07-31"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", type=int, default=13,
                    help="Savitzky-Golay window (13 = published figure, "
                         "9 = Methods-consistent variant)")
    args = ap.parse_args()

    df = pd.read_csv(IN)
    t = (pd.to_datetime(df["timestamp"], format="mixed", utc=True)
         .dt.tz_localize(None))
    df = df.assign(t=t).sort_values("t").reset_index(drop=True)

    da = pd.to_numeric(df["altitude_diff"]).interpolate(method="linear")
    sm = savgol_filter(da.to_numpy(), args.window, 1)
    sl = savgol_filter(da.to_numpy(), args.window, 1, deriv=1)
    import numpy as np
    step_h = np.median(np.diff(df["t"].to_numpy())
                       .astype("timedelta64[s]").astype(float)) / 3600.0
    df = df.assign(sm=sm, sl=sl / step_h)  # km/step -> km/h
    w = df[(df["t"] >= JULY[0]) & (df["t"] <= JULY[1])]
    print(f"median step = {step_h:.3f} h")

    plt.rcdefaults()
    plt.rcParams.update({"font.size": 11, "axes.labelsize": 13,
                         "legend.fontsize": 9, "xtick.labelsize": 10,
                         "ytick.labelsize": 10})

    # ---- fig1a: GNSS vs SGP4 altitude + raw residual -------------------
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.plot(w["t"], w["gps_average_altitude"], color="blue", ls="--",
            lw=0.9, marker=".", ms=2.5, label="GPS Average Altitude")
    ax.plot(w["t"], w["tle_average_altitude"], color="red", ls="--",
            lw=0.9, marker=".", ms=2.5, label="TLE Average Altitude")
    ax.set_xlabel("Time [UTC]")
    ax.set_ylabel("Altitude [km]")
    ax2 = ax.twinx()
    ax2.plot(w["t"], w["altitude_diff"], color="green", lw=1.0,
             marker=".", ms=2.5, label="TLE - GPS Altitude Difference")
    ax2.set_ylabel("Altitude Difference [km]", color="green")
    ax2.tick_params(axis="y", labelcolor="green")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%d"))
    for lab in ax.get_xticklabels():
        lab.set_rotation(45)
        lab.set_horizontalalignment("right")
    ax.grid(True, ls="--", alpha=0.3)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines],
              loc="upper center", ncol=1)
    ax.set_xlim(*JULY)
    fig.tight_layout()
    for ext in ["pdf", "png"]:
        out = OUT / f"fig1a_gnss_vs_sgp4.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight",
                    dpi=150 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)

    # ---- fig1b: smoothed residual + derivative -------------------------
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.plot(w["t"], w["sm"], color="green", lw=1.2,
            label="Altitude Difference [km]")
    ax.set_xlabel("Time [UTC]")
    ax.set_ylabel("Altitude Difference [km]", color="green")
    ax.tick_params(axis="y", labelcolor="green")
    ax2 = ax.twinx()
    ax2.plot(w["t"], w["sl"], color="orange", lw=1.2,
             label=r"$\widetilde{\Delta \dot{a}}(t)$ [km/h]")
    ax2.set_ylabel(r"Time rate of change $\widetilde{\Delta \dot{a}}(t)$ [km/h]",
                   color="orange")
    ax2.tick_params(axis="y", labelcolor="orange")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=3))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    for lab in ax.get_xticklabels():
        lab.set_rotation(45)
        lab.set_horizontalalignment("right")
    ax.grid(True, ls="--", alpha=0.4)
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [ln.get_label() for ln in lines], loc="upper left")
    ax.set_xlim(*JULY)
    fig.tight_layout()
    for ext in ["pdf", "png"]:
        out = OUT / f"fig1b_residual_derivative.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight",
                    dpi=150 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)

    print(f"window={args.window}: smoothed range "
          f"[{w['sm'].min():.4f}, {w['sm'].max():.4f}] km; derivative "
          f"range [{w['sl'].min():+.5f}, {w['sl'].max():+.5f}] km/step")


if __name__ == "__main__":
    main()
