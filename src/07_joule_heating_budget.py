"""Reproduce the manuscript Joule-heating budget and Fig. 7.

This is a portable port of ``solar_activity.ipynb`` cells 19 and 20, the
analysis that generated the paper's quoted values and Fig. 7 curves.  It
uses one-minute PCN only, averaged over the inclusive +/-30-minute window
centred on each hourly OMNI2 timestamp, and signed hourly Dst in the Knipp
et al. (2004) empirical relation (manuscript Eq. 5).

Inputs (repository-relative):
  data/external/omni2/omni2_{2405,2410}_with_utc.csv
  data/external/pcn_1min/PCN_{2405,2410}.csv

Outputs:
  data/figure_data/joule_heating_{2405,2410}.csv
  data/figure_data/joule_heating_{2405,2410}_integrated.csv
  data/figure_data/joule_heating_summary.csv
  figures/fig7_energetics_{2405,2410}.pdf
  figures/fig7_energetics_{2405,2410}.png

Run from any directory with:
  python src/07_joule_heating_budget.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OMNI_DIR = ROOT / "data" / "external" / "omni2"
PCN_DIR = ROOT / "data" / "external" / "pcn_1min"
FIGURE_DATA_DIR = ROOT / "data" / "figure_data"
FIGURE_DIR = ROOT / "figures"

MONTHS = {
    2405: {
        "plot_start": pd.Timestamp("2024-05-10 00:00:00"),
        "plot_end": pd.Timestamp("2024-05-12 00:00:00"),
        "storm_start": pd.Timestamp("2024-05-10 12:00:00"),
        "storm_end": pd.Timestamp("2024-05-12 00:00:00"),
        "quiet_start": pd.Timestamp("2024-05-05 00:00:00"),
        "quiet_end": pd.Timestamp("2024-05-10 00:00:00"),
        "euv_avg": 531.4,
        "euv_max": 597.4,
    },
    2410: {
        "plot_start": pd.Timestamp("2024-10-10 00:00:00"),
        "plot_end": pd.Timestamp("2024-10-12 00:00:00"),
        "storm_start": pd.Timestamp("2024-10-10 12:00:00"),
        "storm_end": pd.Timestamp("2024-10-12 00:00:00"),
        "quiet_start": pd.Timestamp("2024-10-05 00:00:00"),
        "quiet_end": pd.Timestamp("2024-10-10 00:00:00"),
        "euv_avg": 564.3,
        "euv_max": 577.1,
    },
}

EXPECTED = {
    "May storm peak": (1320.3, 0.05),
    "October storm peak": (1758.1, 0.05),
    "May cell-20 mean": (504.2, 0.05),
    "October cell-20 mean": (327.9, 0.05),
    "Closest quiet mean": (95.0, 2.0),
}


def load_hourly_budget(month: int) -> pd.DataFrame:
    """Apply the notebook's PCN averaging and Joule-heating equation."""
    df_dst = pd.read_csv(OMNI_DIR / f"omni2_{month}_with_utc.csv")
    df_pcn = pd.read_csv(PCN_DIR / f"PCN_{month}.csv")

    df_dst["Datetime"] = pd.to_datetime(df_dst["UTC_Time"])
    df_pcn["UTC"] = pd.to_datetime(df_pcn["UTC"])
    df_pcn["PCN"] = df_pcn["PCN"].replace(999, float("nan"))

    averaged_data = []
    for dst_time in df_dst["Datetime"]:
        start_t = dst_time - pd.Timedelta("30min")
        end_t = dst_time + pd.Timedelta("30min")
        data_in_window = df_pcn[
            (df_pcn["UTC"] >= start_t) & (df_pcn["UTC"] <= end_t)
        ]
        if not data_in_window.empty:
            averaged_data.append(data_in_window.mean(numeric_only=True))
        else:
            averaged_data.append(
                pd.Series([None] * len(df_pcn.columns), index=df_pcn.columns)
            )

    averaged_df = pd.DataFrame(averaged_data, columns=df_pcn.columns)
    result = pd.concat([df_dst, averaged_df], axis=1)
    result["heat_flux"] = (
        24.89 * result["PCN"]
        + 3.41 * result["PCN"] ** 2
        + 0.41 * result["Dst_index_nT"]
        + 0.0015 * result["Dst_index_nT"] ** 2
    )
    return result


def integrate_storm(budget: pd.DataFrame, start: pd.Timestamp,
                    end: pd.Timestamp) -> pd.DataFrame:
    """Port cell 20's inclusive storm selection and trapezoidal integral."""
    storm = budget[
        (budget["Datetime"] >= start) & (budget["Datetime"] <= end)
    ].copy()
    integral_heat_flux = np.zeros(len(storm))
    energy = np.zeros(len(storm))

    for i in range(1, len(storm)):
        delta_t = (
            storm["Datetime"].iloc[i] - storm["Datetime"].iloc[i - 1]
        ).total_seconds() / 3600
        integral_heat_flux[i] = (
            storm["heat_flux"].iloc[i] + storm["heat_flux"].iloc[i - 1]
        ) / 2 * delta_t
        energy[i] = energy[i - 1] + integral_heat_flux[i]

    storm["integral_heat_flux"] = integral_heat_flux
    storm["energy"] = energy
    return storm


def closest_quiet_window(budget: pd.DataFrame, end_limit: pd.Timestamp) -> dict:
    """Find the pre-storm whole-calendar-day window closest to 95 GW."""
    first_day = budget["Datetime"].min().normalize()
    best = None
    start = first_day
    while start < end_limit:
        end = start + pd.Timedelta(days=1)
        while end <= end_limit:
            values = budget.loc[
                (budget["Datetime"] >= start) & (budget["Datetime"] < end),
                "heat_flux",
            ]
            if not values.empty:
                candidate = {
                    "difference": abs(values.mean() - 95.0),
                    "start": start,
                    "end": end,
                    "mean": values.mean(),
                }
                if best is None or candidate["difference"] < best["difference"]:
                    best = candidate
            end += pd.Timedelta(days=1)
        start += pd.Timedelta(days=1)
    return best


def save_figure(budget: pd.DataFrame, month: int, config: dict) -> None:
    """Reproduce cell 19's Fig. 7 panel and also emit a PNG preview."""
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot(budget["Datetime"], budget["heat_flux"], label="Joule heating")
    ax.axhline(
        y=config["euv_avg"], color="violet", linestyle="--", linewidth=1.5,
        label="Average EUV power",
    )
    ax.axhline(
        y=config["euv_max"], color="black", linestyle="--", linewidth=1.5,
        label="Maximum EUV power",
    )
    ax.set_xlabel("Time [UTC]", fontsize=25)
    ax.set_ylabel("Power input [GW]", fontsize=25)
    ax.grid(True)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=20)
    ax.set_xlim(config["plot_start"], config["plot_end"])
    ax.legend(fontsize=18)
    fig.tight_layout()

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    stem = FIGURE_DIR / f"fig7_energetics_{month}"
    fig.savefig(stem.with_suffix(".pdf"), format="pdf", bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIGURE_DATA_DIR.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    results = {}
    closest_windows = []

    for month, config in MONTHS.items():
        budget = load_hourly_budget(month)
        budget.to_csv(
            FIGURE_DATA_DIR / f"joule_heating_{month}.csv", index=False
        )

        storm = integrate_storm(
            budget, config["storm_start"], config["storm_end"]
        )
        storm.to_csv(
            FIGURE_DATA_DIR / f"joule_heating_{month}_integrated.csv", index=False
        )
        save_figure(budget, month, config)

        plot_window = budget[
            (budget["Datetime"] >= config["plot_start"])
            & (budget["Datetime"] <= config["plot_end"])
        ]
        quiet = budget[
            (budget["Datetime"] >= config["quiet_start"])
            & (budget["Datetime"] < config["quiet_end"])
        ]
        peak = plot_window["heat_flux"].max()
        mean = storm["heat_flux"].mean()
        quiet_mean = quiet["heat_flux"].mean()
        closest = closest_quiet_window(budget, config["plot_start"])
        closest["month"] = month
        closest_windows.append(closest)
        summary_rows.append(
            {
                "month": month,
                "storm_peak_GW": peak,
                "storm_mean_GW": mean,
                "storm_integral_GWh": storm["energy"].iloc[-1],
                "quiet_mean_GW": quiet_mean,
                "quiet_start": config["quiet_start"],
                "quiet_end_exclusive": config["quiet_end"],
                "closest_quiet_mean_GW": closest["mean"],
                "closest_quiet_start": closest["start"],
                "closest_quiet_end_exclusive": closest["end"],
                "EUV_average_GW": config["euv_avg"],
                "EUV_maximum_GW": config["euv_max"],
            }
        )
        print(
            f"{month}: storm peak={peak:.6f} GW; "
            f"cell-20 mean={mean:.6f} GW; "
            f"quiet mean={quiet_mean:.6f} GW; "
            f"integral={storm['energy'].iloc[-1]:.6f} GWh"
        )

        event = "May" if month == 2405 else "October"
        results[f"{event} storm peak"] = peak
        results[f"{event} cell-20 mean"] = mean

    pd.DataFrame(summary_rows).to_csv(
        FIGURE_DATA_DIR / "joule_heating_summary.csv", index=False
    )

    closest = min(closest_windows, key=lambda item: item["difference"])
    results["Closest quiet mean"] = closest["mean"]
    print(
        "Closest pre-storm whole-calendar-day window to 95 GW: "
        f"{closest['start']} to {closest['end']} (end exclusive), "
        f"{closest['mean']:.6f} GW"
    )

    print("\nExpected / got / match")
    print(f"{'quantity':<26} {'expected':>10} {'got':>12} {'match':>7}")
    for quantity, (expected, tolerance) in EXPECTED.items():
        got = results[quantity]
        match = "yes" if abs(got - expected) <= tolerance else "no"
        print(f"{quantity:<26} {expected:>10.1f} {got:>12.6f} {match:>7}")


if __name__ == "__main__":
    main()
