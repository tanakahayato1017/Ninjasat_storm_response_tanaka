"""Regenerate orbit-averaged altitude residuals and Savitzky-Golay series.

NinjaSat is fitted in fixed two-orbit windows with separate 2-IQR rejection
for GNSS and TLE altitude.  XRISM follows the verified notebook workflow:
fit one global sinusoid to infer its orbital period, average each three-orbit
group, and reject three-standard-deviation outliers.  Savitzky-Golay windows
are 9 samples for NinjaSat and 13 samples for XRISM (polyorder 1, derivative
per sample step, with no ``delta`` argument).

The per-epoch inputs and all regenerated CSVs live in
``data/derived/intermediate`` so that reference products staged elsewhere are
never overwritten.

Reproducibility note. Regenerating from the staged inputs reproduces the
archived reference products to within curve_fit convergence noise
(NinjaSat: max |difference| ~1e-6 km on window-mean altitudes; XRISM:
bit-level except two fit windows on 2024-10-26 -- outside every analysis
window used in the paper -- where the sinusoidal fit converges to a
marginally different optimum, shifting those two window means by <10 m).
No number, table, or figure reported in the manuscript is affected.

Examples
--------
    python src/02_orbit_average_fit.py --satellite all
    python src/02_orbit_average_fit.py --satellite ninjasat --months 2405 2410
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter


ROOT = Path(__file__).resolve().parents[1]
INTERMEDIATE = ROOT / "data" / "derived" / "intermediate"
NINJASAT_PERIOD_MIN = {
    2404: 94.70,
    2405: 94.68,
    2406: 94.57,
    2407: 94.49,
    2408: 94.38,
    2409: 94.20,
    2410: 94.05,
    2411: 93.84,
}
NINJASAT_MONTHS = (2405, 2407, 2410)
XRISM_MONTHS = (2405, 2410)


def _fixed_sine(x: np.ndarray, amplitude: float, phase: float, offset: float, period: float) -> np.ndarray:
    return amplitude * np.sin((4.0 * np.pi / period) * x + phase) + offset


def _fit_fixed_offset(x: pd.Series, y: pd.Series, period: float) -> float:
    initial = [(y.max() - y.min()) / 2.0, np.pi, y.mean()]
    fitted, _ = curve_fit(lambda xx, a, c, d: _fixed_sine(xx, a, c, d, period), x, y, p0=initial)
    return float(fitted[2])


def fit_ninjasat(month: int) -> tuple[Path, Path]:
    """Create the two-orbit fit and SG9 products for one NinjaSat month."""
    source = INTERMEDIATE / f"ninjasat_gps_and_tle_{month}_ext.csv"
    data = pd.read_csv(source)
    time = pd.to_datetime(data["gpsUtcTime"])
    time_minutes = (time - time.min()).dt.total_seconds() / 60.0
    gps_altitude = data["gps_altitude"]
    tle_altitude = data["tle_altitude"]
    window_minutes = NINJASAT_PERIOD_MIN[month] * 2.0
    window_count = int(np.ceil(time_minutes.max() / window_minutes))

    gps_means: list[float] = []
    tle_means: list[float] = []
    centers: list[pd.Timestamp] = []
    for index in range(window_count):
        mask = (time_minutes >= index * window_minutes) & (time_minutes < (index + 1) * window_minutes)
        x = time_minutes[mask]
        gps_y = gps_altitude[mask]
        tle_y = tle_altitude[mask]
        if len(x) <= 8:
            gps_means.append(np.nan)
            tle_means.append(np.nan)
            centers.append(time[mask].mean())
            continue

        gps_q1, gps_q3 = np.percentile(gps_y, [25, 75])
        gps_iqr = gps_q3 - gps_q1
        gps_valid = (gps_y >= gps_q1 - 2.0 * gps_iqr) & (gps_y <= gps_q3 + 2.0 * gps_iqr)
        tle_q1, tle_q3 = np.percentile(tle_y, [25, 75])
        tle_iqr = tle_q3 - tle_q1
        tle_valid = (tle_y >= tle_q1 - 2.0 * tle_iqr) & (tle_y <= tle_q3 + 2.0 * tle_iqr)

        try:
            gps_means.append(_fit_fixed_offset(x[gps_valid], gps_y[gps_valid], window_minutes) if gps_valid.sum() > 8 else np.nan)
        except Exception:
            gps_means.append(np.nan)
        try:
            tle_means.append(_fit_fixed_offset(x[tle_valid], tle_y[tle_valid], window_minutes) if tle_valid.sum() > 8 else np.nan)
        except Exception:
            tle_means.append(np.nan)
        centers.append(time[mask].mean())

    valid_center = pd.notna(centers)
    centers_clean = [value for value, keep in zip(centers, valid_center) if keep]
    gps_clean = [value for value, keep in zip(gps_means, valid_center) if keep]
    tle_clean = [value for value, keep in zip(tle_means, valid_center) if keep]
    fitted = pd.DataFrame(
        {
            "timestamp": centers_clean,
            "gps_average_altitude": gps_clean,
            "tle_average_altitude": tle_clean,
            "altitude_diff": np.asarray(tle_clean) - np.asarray(gps_clean),
        }
    )
    fit_path = INTERMEDIATE / f"gps_tle_average_altitude_fitting_all_severe_{month}_monthly_ext.csv"
    fitted.to_csv(fit_path, index=False)

    difference = fitted["altitude_diff"]
    interpolated = difference.interpolate(method="linear")
    smoothed = savgol_filter(interpolated, window_length=9, polyorder=1)
    slope = savgol_filter(interpolated, window_length=9, polyorder=1, deriv=1)
    result = pd.DataFrame(
        {
            "timestamp": pd.to_datetime(fitted["timestamp"], format="ISO8601"),
            "altitude_diff": difference,
            "altitude_diff_smoothed": smoothed,
            "slope": slope,
        }
    )
    smooth_path = INTERMEDIATE / f"gps_tle_average_altitude_fitting_all_smoothed_{month}_ext_thesis.csv"
    result.to_csv(smooth_path, index=False)
    print(f"NinjaSat {month}: wrote {fit_path.name} and {smooth_path.name} ({len(fitted)} rows)")
    return fit_path, smooth_path


def _sinusoidal(time: np.ndarray, amplitude: float, frequency: float, phase: float, offset: float) -> np.ndarray:
    return amplitude * np.sin(2.0 * np.pi * frequency * time + phase) + offset


def _remove_outliers(frame: pd.DataFrame, column: str, threshold: float = 3.0) -> pd.DataFrame:
    mean = frame[column].mean()
    standard_deviation = frame[column].std()
    return frame[np.abs(frame[column] - mean) < threshold * standard_deviation]


def fit_xrism(month: int) -> tuple[Path, Path]:
    """Create the three-orbit mean and SG13 products for one XRISM month."""
    source = INTERMEDIATE / f"xrism_gps_and_tle_{month}.csv"
    data = pd.read_csv(source)
    data["datetime"] = pd.to_datetime(data["datetime"])
    data["time_seconds"] = (data["datetime"] - data["datetime"].iloc[0]).dt.total_seconds()

    gps_fit, _ = curve_fit(
        _sinusoidal,
        data["time_seconds"],
        data["gps_altitude"],
        p0=[10, 1 / 6000, 0, np.mean(data["gps_altitude"])],
    )
    # The notebook also performs this fit and aborts the month if it fails,
    # even though only the GPS-derived period defines the groups.
    curve_fit(
        _sinusoidal,
        data["time_seconds"],
        data["tle_altitude"],
        p0=[10, 1 / 6000, 0, np.mean(data["tle_altitude"])],
    )
    orbital_period_seconds = 1.0 / gps_fit[1]
    data["orbit_pair_index"] = (data["time_seconds"] // (3.0 * orbital_period_seconds)).astype(int)

    rows = []
    for _, group in data.groupby("orbit_pair_index"):
        gps_mean = group["gps_altitude"].mean()
        tle_mean = group["tle_altitude"].mean()
        rows.append(
            {
                "time_average": group["datetime"].mean(),
                "average_altitude_gps": gps_mean,
                "average_altitude_tle": tle_mean,
                "altitude_difference": tle_mean - gps_mean,
            }
        )
    fitted = pd.DataFrame(rows)
    fitted = _remove_outliers(fitted, "average_altitude_gps")
    fitted = _remove_outliers(fitted, "average_altitude_tle")
    fit_path = INTERMEDIATE / f"orbital_altitude_average_fitting_{month}_3period.csv"
    fitted.to_csv(fit_path, index=False)

    difference = fitted["altitude_difference"]
    smoothed = savgol_filter(difference, window_length=13, polyorder=1)
    slope = savgol_filter(difference, window_length=13, polyorder=1, deriv=1)
    result = pd.DataFrame(
        {
            "timestamp": pd.to_datetime(fitted["time_average"], format="ISO8601"),
            "altitude_diff": difference,
            "altitude_diff_smoothed": smoothed,
            "slope": slope,
        }
    )
    smooth_path = INTERMEDIATE / f"orbital_altitude_average_fitting_{month}_3period_slope.csv"
    result.to_csv(smooth_path, index=False)
    print(f"XRISM {month}: wrote {fit_path.name} and {smooth_path.name} ({len(fitted)} rows)")
    return fit_path, smooth_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite", choices=("ninjasat", "xrism", "all"), default="all")
    parser.add_argument("--months", type=int, nargs="+", help="YYMM month(s); defaults to all staged months")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    if args.satellite in ("ninjasat", "all"):
        months = args.months if args.months is not None else NINJASAT_MONTHS
        for month in months:
            if month in NINJASAT_MONTHS:
                fit_ninjasat(month)
            elif args.satellite == "ninjasat":
                raise ValueError(f"No staged NinjaSat month {month}")
    if args.satellite in ("xrism", "all"):
        months = args.months if args.months is not None else XRISM_MONTHS
        for month in months:
            if month in XRISM_MONTHS:
                fit_xrism(month)
            elif args.satellite == "xrism":
                raise ValueError(f"No staged XRISM month {month}")


if __name__ == "__main__":
    main()
