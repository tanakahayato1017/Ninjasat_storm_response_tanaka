"""Figure 3c,d: storm-time EUV baseline subtraction.

Plot the 140.5-nm GOES-18 irradiance together with the observed NinjaSat
altitude-change rate, the quiet-time EUV regression estimate, and their
residual.  The plotted series are consumed directly from the outputs of
``src/03_euv_baseline_regression.py``; no regression is recomputed here.

Unit correction (author-approved, following the Fig. 6 precedent): the
submitted panels' right axis read "[km/s]" while the underlying rates
are km per ~3.1-3.2-h sample step; the three rate series are divided by
each month's median sampling step and the axis (and its limits) are in
km/h.  Curve shapes are unchanged (a per-month rescale).

Inputs:
  data/derived/correlation_inputs/
      regression_subtract_slope_{2405,2410}_thesis.csv

Outputs:
  figures/fig3c_euv_subtraction_2405.pdf/png
  figures/fig3d_euv_subtraction_2410.pdf/png

Run (from repository root):
  python figures/fig3cd_euv_subtraction.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INPUT = (ROOT / "data" / "derived" / "correlation_inputs" /
         "regression_subtract_slope_{}_thesis.csv")
OUT = ROOT / "figures"

WINDOWS = {
    "2405": (pd.Timestamp("2024-05-07"), pd.Timestamp("2024-05-15")),
    "2410": (pd.Timestamp("2024-10-07"), pd.Timestamp("2024-10-15")),
}
STEMS = {
    "2405": "fig3c_euv_subtraction_2405",
    "2410": "fig3d_euv_subtraction_2410",
}


def load_window(month: str) -> tuple[pd.DataFrame, float]:
    """Load the precomputed subtraction series and select the plot window.

    Returns the windowed frame with the three rate columns converted from
    km per sample step to km/h, plus the median step [h] used."""
    import numpy as np
    data = pd.read_csv(str(INPUT).format(month))
    data["timestamp"] = pd.to_datetime(data["timestamp"], format="mixed")
    step_h = np.median(np.diff(data["timestamp"].to_numpy())
                       .astype("timedelta64[s]").astype(float)) / 3600.0
    for column in ["slope", "estimated_slope", "slope_difference"]:
        data[column] = data[column] / step_h
    lo, hi = WINDOWS[month]
    return data.loc[data["timestamp"].between(lo, hi)].copy(), step_h


def report_ranges(month: str, data: pd.DataFrame) -> None:
    """Print compact numerical provenance for every plotted ordinate."""
    ranges = []
    for column in ["irr_1405", "slope", "estimated_slope", "slope_difference"]:
        ranges.append(
            f"{column}=[{data[column].min():.12g}, {data[column].max():.12g}]"
        )
    print(f"[{month}] n={len(data)}; " + "; ".join(ranges))


def plot_panel(month: str, data: pd.DataFrame, step_h: float) -> None:
    """Render one submitted-style two-axis storm panel (rates in km/h)."""
    lo, hi = WINDOWS[month]
    fig, ax_left = plt.subplots(figsize=(15, 7))
    ax_right = ax_left.twinx()

    line_euv, = ax_left.plot(
        data["timestamp"], data["irr_1405"], color="purple", lw=1.5,
        label="EUV 140.5 nm irradiance",
    )
    line_observed, = ax_right.plot(
        data["timestamp"], data["slope"], color="orange", lw=1.5,
        label=r"$\widetilde{\Delta \dot{a}}(t)$ (observed)",
    )
    line_estimated, = ax_right.plot(
        data["timestamp"], data["estimated_slope"], color="red", lw=1.5,
        label=r"$\widetilde{\Delta \dot{a}}_{\mathrm{EUV}}(t)$ (estimated)",
    )
    line_residual, = ax_right.plot(
        data["timestamp"], data["slope_difference"], color="black", lw=1.5,
        label="Residual (observed − estimated)",
    )

    ax_left.set_xlim(lo, hi)
    ax_right.set_ylim(-0.03 / step_h, 0.05 / step_h)
    ax_left.set_xlabel("Time [UTC]", fontsize=22)
    ax_left.set_ylabel(r"EUV 140.5 nm irradiance [W/m$^2$]",
                       color="purple", fontsize=22)
    ax_right.set_ylabel(
        r"Time rates of change $\widetilde{\Delta \dot{a}}(t)$ [km/h]",
        fontsize=22,
    )
    ax_left.tick_params(axis="y", colors="purple", labelsize=18)
    ax_left.tick_params(axis="x", labelsize=18, rotation=45)
    ax_right.tick_params(axis="y", labelsize=18)
    ax_left.xaxis.set_major_locator(mdates.DayLocator(interval=2))
    ax_left.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))
    ax_left.grid(True, linestyle="--", alpha=0.6)
    ax_right.grid(False)
    ax_left.legend(
        [line_euv, line_observed, line_estimated, line_residual],
        [line.get_label() for line in
         [line_euv, line_observed, line_estimated, line_residual]],
        loc="upper left", fontsize=16,
    )
    fig.tight_layout()

    for extension in ["pdf", "png"]:
        output = OUT / f"{STEMS[month]}.{extension}"
        fig.savefig(output, format=extension, bbox_inches="tight",
                    dpi=150 if extension == "png" else None)
        print(f"wrote {output}")
    plt.close(fig)


def main() -> None:
    for month in ["2405", "2410"]:
        data, step_h = load_window(month)
        print(f"[{month}] median step = {step_h:.3f} h (rates shown in km/h)")
        report_ranges(month, data)
        plot_panel(month, data, step_h)


if __name__ == "__main__":
    main()
