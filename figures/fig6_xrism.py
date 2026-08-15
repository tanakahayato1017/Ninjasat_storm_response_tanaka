"""Figure 6: NinjaSat vs. XRISM smoothed altitude-rate response (May and
October 2024 storms), in km/h.

Plots the pre-computed Savitzky-Golay first-derivative ("slope") column
from each satellite's altitude-fitting pipeline:
  - NinjaSat: data/derived/ninjasat/gps_tle_average_altitude_fitting_all_
    smoothed_{2405,2410}_ext_thesis.csv
  - XRISM: data/derived/xrism/orbital_altitude_average_fitting_{2405,2410}
    _3period_slope.csv

Both "slope" columns are savgol_filter(..., deriv=1): the analytic
derivative of a local polynomial fit to the altitude-difference series,
i.e. km per sample step (NinjaSat: ~3.17 h, XRISM: ~5.0 h). Because the two
satellites are sampled at different steps, this script converts each to a
common km/h rate by dividing by that satellite's median in-window sample
interval -- a per-satellite rescaling only; curve shapes are unchanged.

Output:
  fig6a_xrism_2405.pdf/png (May 2024)
  fig6b_xrism_2410.pdf/png (October 2024)

Run (from repository root):
  python figures/fig6_xrism.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"
OUT = ROOT / "figures"

NINJA = DATA / "ninjasat" / "gps_tle_average_altitude_fitting_all_smoothed_{}_ext_thesis.csv"
XRISM = DATA / "xrism" / "orbital_altitude_average_fitting_{}_3period_slope.csv"

XLIM = {"2405": (pd.Timestamp("2024-05-07"), pd.Timestamp("2024-05-15")),
        "2410": (pd.Timestamp("2024-10-08"), pd.Timestamp("2024-10-15"))}
PANEL = {"2405": "fig6a_xrism_2405", "2410": "fig6b_xrism_2410"}


def load_slope(path_tmpl: Path, tag: str) -> pd.DataFrame:
    df = pd.read_csv(str(path_tmpl).format(tag))
    t = (pd.to_datetime(df["timestamp"], format="mixed", utc=True)
         .dt.tz_localize(None))
    return (df.assign(dt=t, slope=df["slope"].astype(float))
              .sort_values("dt").reset_index(drop=True)[["dt", "slope"]])


def rate_km_per_h(df: pd.DataFrame, lo, hi) -> tuple[pd.DataFrame, float]:
    """slope [km/step] -> km/h using the median in-window sample step."""
    w = df[(df["dt"] >= lo) & (df["dt"] <= hi)].copy()
    step_h = np.median(np.diff(w["dt"].to_numpy())
                       .astype("timedelta64[s]").astype(float)) / 3600.0
    w["rate"] = w["slope"] / step_h
    return w, step_h


def main() -> None:
    plt.rcdefaults()
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })

    for tag in ["2405", "2410"]:
        lo, hi = XLIM[tag]
        nin, nin_step = rate_km_per_h(load_slope(NINJA, tag), lo, hi)
        xri, xri_step = rate_km_per_h(load_slope(XRISM, tag), lo, hi)
        print(f"[{tag}] NinjaSat: step={nin_step:.3f} h, n={len(nin)}, "
              f"peak={np.nanmax(nin['rate']):+.5f} km/h | "
              f"XRISM: step={xri_step:.3f} h, n={len(xri)}, "
              f"peak={np.nanmax(xri['rate']):+.5f} km/h")

        fig, ax = plt.subplots(figsize=(9, 3.2))
        ax.plot(xri["dt"], xri["rate"], color="tab:blue", lw=1.4, label="XRISM")
        ax.plot(nin["dt"], nin["rate"], color="tab:red", lw=1.4, label="NinjaSat")
        ax.set_xlim(lo, hi)
        ax.set_xlabel("Time [UTC]")
        ax.set_ylabel(r"$\widetilde{\Delta \dot{a}}(t)$ [km/h]")
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
        for lab in ax.get_xticklabels():
            lab.set_rotation(45)
            lab.set_horizontalalignment("right")
        ax.grid(True, linestyle="-", alpha=0.35)
        ax.legend(loc="upper right")
        fig.tight_layout()

        for ext in ["pdf", "png"]:
            out = OUT / f"{PANEL[tag]}.{ext}"
            fig.savefig(out, format=ext, bbox_inches="tight",
                        dpi=150 if ext == "png" else None)
            print(f"wrote {out}")
        plt.close(fig)


if __name__ == "__main__":
    main()
