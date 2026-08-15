"""Render manuscript Figs. 4e,f and 5e,f for the AE/AL analysis.

Inputs (repository-relative):
  data/figure_data/ae_al_timeseries_{2405,2410}.csv
  data/figure_data/ae_al_crosscorr_curve_{2405,2410}.csv

Outputs:
  figures/fig4e_ae_al_covariation_2405.pdf/png
  figures/fig4f_ae_al_covariation_2410.pdf/png
  figures/fig5e_ae_al_crosscorr_2405.pdf/png
  figures/fig5f_ae_al_crosscorr_2410.pdf/png

Run (from repository root):
  python figures/fig4_5_ae_al.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
import _style as st  # noqa: E402

DATA = ROOT / "data" / "figure_data"
FIGURES = ROOT / "figures"

MONTH_LABEL = {"2405": "May 2024", "2410": "October 2024"}

def save(fig, stem: str) -> None:

    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ["png", "pdf"]:
        path = FIGURES / f"{stem}.{ext}"
        fig.savefig(path, format=ext, bbox_inches="tight",
                    dpi=150 if ext == "png" else None)
        print(f"wrote {path}")

INDEX_STYLE_4 = {
    "AE": dict(color=st.C_PURPLE, ls="-", label="AE index"),
    "AL": dict(color=st.C_ORANGE, ls=(0, (4, 2)), label="AL index"),
}

XLIM_DATES = {
    "2405": (pd.Timestamp("2024-05-07"), pd.Timestamp("2024-05-18")),
    "2410": (pd.Timestamp("2024-10-06"), pd.Timestamp("2024-10-17")),
}

def covariation_panel(tag: str, stem: str) -> None:

    plt.rcdefaults()
    df = pd.read_csv(DATA / f"ae_al_timeseries_{tag}.csv", parse_dates=["datetime"])
    start_time, end_time = XLIM_DATES[tag]

    fig, ax = plt.subplots(figsize=(14, 8))

    ax.set_title("AE/AL indices and EUV-subtracted rate of change\n"
                 f"for {MONTH_LABEL[tag]}", fontsize=20)

    lines = []
    for idx in ["AE", "AL"]:
        s = INDEX_STYLE_4[idx]
        l, = ax.plot(df["datetime"], df[idx], label=s["label"], color=s["color"],
                     linestyle=s["ls"])
        lines.append(l)

    ax.set_xlabel("Time [UTC]", fontsize=25)

    ax.set_ylabel("AE, AL index [nT]", fontsize=25)

    ax.tick_params(axis="x", rotation=45, labelsize=20)

    ax.tick_params(axis="y", labelsize=20)

    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))

    ax.grid(linestyle="--", alpha=0.7)

    ax.axhline(0, color=st.C_GRAY, linewidth=0.8, zorder=0)

    ax2 = ax.twinx()

    # slope_difference is stored in km per ~3.1-3.2-h sample step; divide
    # by the median step so the km/h axis label is numerically correct
    # (the submitted fig4e/f carried the km/h label over per-step values).
    step_h = np.median(np.diff(pd.to_datetime(df["datetime"]).to_numpy())
                       .astype("timedelta64[s]").astype(float)) / 3600.0
    slope_line, = ax2.plot(df["datetime"], df["slope_difference"] / step_h,
                           "k",
                           label=r"$\widetilde{\dot{\Delta a}}(t) - "
                                 r"\widetilde{\dot{\Delta a}}_\mathrm{EUV}(t)$")

    ax2.set_ylabel("EUV-subtracted rate of change\n"
                   r"$\widetilde{\dot{\Delta a}}(t) - "
                   r"\widetilde{\dot{\Delta a}}_\mathrm{EUV}(t)$ [km/h]",
                   fontsize=20)

    ax2.tick_params(axis="y", labelsize=20)

    ax2.grid(False)

    all_lines = lines + [slope_line]
    labels = [l.get_label() for l in all_lines]

    ax.legend(all_lines, labels, loc="upper left", fontsize=16)

    ax.set_xlim(start_time, end_time)

    fig.tight_layout()

    fig.subplots_adjust(right=0.83)

    save(fig, stem)
    plt.close(fig)

CURVE_STYLE = {
    "AE": dict(color=st.C_PURPLE, ls="-"),
    "AL": dict(color=st.C_ORANGE, ls=(0, (4, 2))),
}

def _peak_row(curve: pd.DataFrame, idx: str) -> pd.Series:

    sub = curve[curve["index"] == idx]
    return sub.loc[sub["r"].abs().idxmax()]

def _plot_one_series(ax, curve, idx):

    s = CURVE_STYLE[idx]
    sub = curve[curve["index"] == idx].sort_values("lag_h").reset_index(drop=True)
    peak = _peak_row(curve, idx)
    lo, hi = sorted([peak["r_ci_lo"], peak["r_ci_hi"]])
    in_ci = sub["r"].between(lo, hi)

    ax.plot(sub["lag_h"], sub["r"], marker="o", label=f"{idx} Correlation",
           linestyle=s["ls"], color=s["color"])

    ax.plot(sub["lag_h"][in_ci], sub["r"][in_ci], marker="o",
           label="Inside 95% CI", linestyle="-", color="green")

    ax.axhline(lo, color="blue", linestyle="--", label="95% Confidence Interval")

    ax.axhline(hi, color="blue", linestyle="--")

    ax.scatter(peak["lag_h"], peak["r"], color="red",
              label=f"Peak Correlation ({idx})", zorder=5)

    x_in_ci = sub["lag_h"][in_ci]
    y_in_ci = sub["r"][in_ci]

    x_min, x_max = x_in_ci.min(), x_in_ci.max()

    y_min = y_in_ci.loc[x_in_ci.idxmin()]
    y_max = y_in_ci.loc[x_in_ci.idxmax()]

    ax.fill_between(x_in_ci, y_in_ci, color="green", alpha=0.3,
                    label="Confidence Interval Area")

    ax.plot([x_min, x_min], [0, y_min], color="green", linestyle="--", linewidth=1)

    ax.plot([x_max, x_max], [0, y_max], color="green", linestyle="--", linewidth=1)

def crosscorr_panel(tag: str, stem: str) -> None:
    plt.rcdefaults()
    curve = pd.read_csv(DATA / f"ae_al_crosscorr_curve_{tag}.csv")

    fig = plt.figure(figsize=(10, 6))
    ax = fig.gca()

    _plot_one_series(ax, curve, "AE")
    _plot_one_series(ax, curve, "AL")

    ax.set_xticks(np.arange(-20, 20, 2))
    ax.tick_params(axis="x", labelsize=16)

    ax.tick_params(axis="y", labelsize=16)

    ax.axhline(0, color="gray", linestyle="--", linewidth=1)

    ax.set_xlabel("Time Shift (hours)", fontsize=22)

    ax.set_ylabel("Correlation Coefficient", fontsize=22)

    handles, labels = ax.get_legend_handles_labels()
    seen = set()
    uniq_handles, uniq_labels = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen.add(l)
            uniq_handles.append(h)
            uniq_labels.append(l)
    ax.legend(uniq_handles, uniq_labels, fontsize=14)

    ax.grid(True, linestyle="--", linewidth=0.7)

    fig.tight_layout()

    ax.set_xlim(-20, 20)

    ax.set_ylim(-0.9, 0.9)

    save(fig, stem)
    plt.close(fig)

def main() -> None:
    covariation_panel("2405", "fig4e_ae_al_covariation_2405")
    covariation_panel("2410", "fig4f_ae_al_covariation_2410")
    crosscorr_panel("2405", "fig5e_ae_al_crosscorr_2405")
    crosscorr_panel("2410", "fig5f_ae_al_crosscorr_2410")

if __name__ == "__main__":
    main()
