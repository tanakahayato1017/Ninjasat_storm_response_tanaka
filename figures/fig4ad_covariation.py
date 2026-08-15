"""Figure 4a-d: visual co-variation of the PCN and |Dst| indices with the
EUV-subtracted altitude rate for the May and October 2024 storms, and the
n=84 merged series that feed the lag-correlation analysis (src/04).

Faithful port of the original notebook cell (correlations_2.ipynb cell 6):
for each (index, month) pair the geomagnetic index is averaged in a
+/-1.57-h box centred on each altitude-rate epoch, attached with a
nearest-timestamp asof merge, filtered to values < 900 (which also drops
fill values / empty boxes), Dst sign-flipped and displayed as |Dst|, and
windowed to 2024-05-07..05-18 / 2024-10-06..10-17 (n = 84 samples each).
The merged output equals the staged
data/derived/correlation_inputs/{month}_{index}_correlation_thesis.csv
files used by src/04 - regenerated here and verified against them.

Corrections applied (author-approved, following the Fig. 6 precedent):
(1) The submitted fig4a-d panels carried an in-figure title hard-coded
to "for May 2024" - including the October panels - a copy-paste slip in
the original notebook; panels are now titled with their actual month.
(2) The submitted right axis read "km/s" while the values are km per
~3.1-3.2-h sample step; the displayed rate is now divided by each
month's median sampling step and labelled km/h (curve shape unchanged;
the merged CSV columns verified against the staged n=84 files remain in
the original per-step units used by the correlation analysis).

Inputs:
  data/external/pcn_1min/PCN_{2405,2410}.csv      (1-min PCN; 999 -> NaN)
  data/external/omni2/omni2_{2405,2410}_with_utc.csv (hourly signed Dst)
  data/derived/correlation_inputs/regression_subtract_slope_{month}_thesis.csv
      (EUV-subtracted rate, src/03 output)

Outputs:
  fig4{a,b}_pcn_covariation_{2405,2410}.pdf/png
  fig4{c,d}_dst_covariation_{2405,2410}.pdf/png

Run (from repository root):
  python figures/fig4ad_covariation.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
PCN = ROOT / "data" / "external" / "pcn_1min" / "PCN_{}.csv"
OMNI = ROOT / "data" / "external" / "omni2" / "omni2_{}_with_utc.csv"
SLOPE = (ROOT / "data" / "derived" / "correlation_inputs" /
         "regression_subtract_slope_{}_thesis.csv")
REF = (ROOT / "data" / "derived" / "correlation_inputs" /
       "{}_{}_correlation_thesis.csv")
OUT = ROOT / "figures"

WINDOWS = {2405: (pd.Timestamp("2024-05-07"), pd.Timestamp("2024-05-18")),
           2410: (pd.Timestamp("2024-10-06"), pd.Timestamp("2024-10-17"))}
MONTH_NAME = {2405: "May 2024", 2410: "October 2024"}
PANEL = {("PCN", 2405): "fig4a_pcn_covariation_2405",
         ("PCN", 2410): "fig4b_pcn_covariation_2410",
         ("Dst_index_nT", 2405): "fig4c_dst_covariation_2405",
         ("Dst_index_nT", 2410): "fig4d_dst_covariation_2410"}
BOX_HALF_H = 1.57


def build_merged(index: str, month: int) -> pd.DataFrame:
    if index == "Dst_index_nT":
        d1 = pd.read_csv(str(OMNI).format(month))
    else:
        d1 = pd.read_csv(str(PCN).format(month))
        d1["UTC_Time"] = d1["UTC"]
        d1[index] = d1[index].replace(999, np.nan)
    d2 = pd.read_csv(str(SLOPE).format(month))

    d1["Time"] = pd.to_datetime(d1["UTC_Time"], errors="coerce")
    d2["timestamp"] = pd.to_datetime(d2["timestamp"], format="mixed",
                                     errors="coerce")
    d1["Time"] = d1["Time"].dt.tz_localize(None)
    d2["timestamp"] = d2["timestamp"].dt.tz_localize(None)

    aggregated = []
    for t in d2["timestamp"]:
        win = d1[(d1["Time"] >= t - pd.Timedelta(hours=BOX_HALF_H))
                 & (d1["Time"] <= t + pd.Timedelta(hours=BOX_HALF_H))]
        mean_values = win.mean(numeric_only=True)
        mean_values["timestamp"] = t
        aggregated.append(mean_values)
    agg = pd.DataFrame(aggregated).rename(columns={"timestamp": "Time"})
    agg["Time"] = pd.to_datetime(agg["Time"], errors="coerce")

    merged = pd.merge_asof(d2.sort_values("timestamp"),
                           agg.sort_values("Time"),
                           left_on="timestamp", right_on="Time",
                           direction="nearest").drop(columns=["Time"])

    merged = merged[merged[index] < 900]
    if index == "Dst_index_nT":
        merged[index] = -merged[index]
    lo, hi = WINDOWS[month]
    merged = merged[(merged["timestamp"] >= lo) & (merged["timestamp"] <= hi)]
    merged[index] = merged[index].abs()
    return merged


def verify(merged: pd.DataFrame, index: str, month: int) -> None:
    ref = pd.read_csv(str(REF).format(month, index))
    ref["timestamp"] = pd.to_datetime(ref["timestamp"])
    mm = ref.merge(merged[["timestamp", index, "slope_difference"]],
                   on="timestamp", suffixes=("_ref", "_new"))
    d_idx = np.nanmax(np.abs(mm[f"{index}_ref"] - mm[f"{index}_new"]))
    d_sd = np.nanmax(np.abs(mm["slope_difference_ref"]
                            - mm["slope_difference_new"]))
    print(f"  verify vs staged n=84 file: rows {len(mm)}/{len(ref)}, "
          f"max|d({index})| = {d_idx:.3e}, "
          f"max|d(slope_difference)| = {d_sd:.3e}")


def plot(merged: pd.DataFrame, index: str, month: int) -> None:
    lo, hi = WINDOWS[month]
    label = "|Dst| index" if index == "Dst_index_nT" else "PCN index"
    unit = "|Dst| index [nT]" if index == "Dst_index_nT" else "PCN index [mV/m]"

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_title(f"{label} and EUV-subtracted rate of change\n"
                 f"for {MONTH_NAME[month]}", fontsize=20)
    ax.plot(merged["timestamp"], merged[index],
            label=(f"{label} [nT]" if index == "Dst_index_nT" else label),
            color="blue")
    ax.set_xlabel("Time [UTC]", fontsize=25)
    ax.set_ylabel(unit, fontsize=25, color="blue")
    ax.tick_params(axis="x", rotation=45, labelsize=20)
    ax.tick_params(axis="y", labelsize=20, labelcolor="blue")
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax.grid(linestyle="--", alpha=0.7)

    step_h = np.median(np.diff(merged["timestamp"].to_numpy())
                       .astype("timedelta64[s]").astype(float)) / 3600.0
    ax2 = ax.twinx()
    slope_line, = ax2.plot(
        merged["timestamp"], merged["slope_difference"] / step_h, "k",
        label=r"$\widetilde{\dot{\Delta a}}(t) - "
              r"\widetilde{\dot{\Delta a}}_\mathrm{EUV}(t)$")
    ax2.set_ylabel("EUV-subtracted rate of change\n"
                   r"$\widetilde{\dot{\Delta a}}(t) - "
                   r"\widetilde{\dot{\Delta a}}_\mathrm{EUV}(t)$ [km/h]",
                   fontsize=20)
    ax2.tick_params(axis="y", labelsize=20)

    lines = ax.get_lines() + [slope_line]
    ax.legend(lines, [ln.get_label() for ln in lines],
              loc="upper left", fontsize=16)
    plt.xlim(lo, hi)
    fig.tight_layout()
    fig.subplots_adjust(right=0.83)

    for ext in ["pdf", "png"]:
        out = OUT / f"{PANEL[(index, month)]}.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight",
                    dpi=150 if ext == "png" else None)
        print(f"  wrote {out}")
    plt.close(fig)


def main() -> None:
    for index in ["PCN", "Dst_index_nT"]:
        for month in [2405, 2410]:
            print(f"[{index} {month}]")
            merged = build_merged(index, month)
            verify(merged, index, month)
            plot(merged, index, month)


if __name__ == "__main__":
    main()
