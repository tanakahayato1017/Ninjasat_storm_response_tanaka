"""Lomb-Scargle periodograms of the July 2024 (quiet-month) NinjaSat
altitude series, addressing the sub-diurnal oscillation visible in the raw
residual (manuscript Fig. 1a) and confirming it is not an SGP4-propagation
artefact.

(a) Raw (unsmoothed) altitude residual Delta a = a_SGP4 - a_GNSS
    (data/raw/ninjasat/gps_tle_average_altitude_fitting_all_severe_2407
    _monthly_ext.csv, column altitude_diff -- this file also carries the
    absolute gps_average_altitude column, so it is treated as pending
    NinjaSat team release like the rest of data/raw/, see that folder's
    README), exactly the green curve of manuscript Fig. 1a.
(b) The J2-corrected specific orbital energy E_J2 (data/raw/energy/
    specific_energy.csv, column E_J2_Jkg -- pending NinjaSat team release,
    see data/raw/energy/README.md), computed at each GNSS epoch directly
    from GNSS-derived position and velocity, with no SGP4 propagation
    involved:

    E_J2 = (1/2) v_i^2 - (mu/r) [1 - (J2/2)(R_E/r)^2 (3 sin^2(phi) - 1)]

    Both series are restricted to July 2024. Because E_J2 depends only on
    the GNSS state at each epoch, any periodicity it shares with the raw
    Delta a series cannot be an SGP4-propagation artefact.

Output:
  fig1_lomb_scargle.pdf/png

Run (from repository root):
  python figures/fig1_lomb_scargle.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from astropy.timeseries import LombScargle  # noqa: E402
from scipy.signal import detrend  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "derived"
OUT = ROOT / "figures"

DA_CSV = ROOT / "data" / "raw" / "ninjasat" / "gps_tle_average_altitude_fitting_all_severe_2407_monthly_ext.csv"
EJ2_CSV = ROOT / "data" / "raw" / "energy" / "specific_energy.csv"

JULY_START, JULY_END = "2024-07-01", "2024-08-01"  # [start, end)
PERIOD_MIN_H, PERIOD_MAX_H = 6.0, 50.0
N_PERIOD_SAMPLES = 20000
MARK_PERIOD_H = 12.0


def _load_da_july() -> tuple[np.ndarray, np.ndarray, float]:
    df = pd.read_csv(DA_CSV)
    df["timestamp"] = pd.to_datetime(df["timestamp"], format="mixed", utc=True)
    sub = df[(df["timestamp"] >= JULY_START) & (df["timestamp"] < JULY_END)].copy()
    t = (sub["timestamp"] - sub["timestamp"].iloc[0]).dt.total_seconds().to_numpy() / 3600.0
    y = sub["altitude_diff"].to_numpy()
    dt_median = np.median(np.diff(t))
    nyquist_period = 2.0 * dt_median
    print(f"Delta-a (Fig 1a green): n={len(sub)}, median dt={dt_median:.3f} h, "
          f"Nyquist period={nyquist_period:.3f} h")
    return t, y, nyquist_period


def _load_ej2_july() -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(EJ2_CSV)
    df["datetime"] = pd.to_datetime(df["datetime"], format="mixed", utc=True)
    sub = df[(df["datetime"] >= JULY_START) & (df["datetime"] < JULY_END)].copy()
    t = (sub["datetime"] - sub["datetime"].iloc[0]).dt.total_seconds().to_numpy() / 3600.0
    y = detrend(sub["E_J2_Jkg"].to_numpy(), type="linear")
    print(f"E_J2 proxy: n={len(sub)}, native cadence ~10 min "
          f"(median dt={np.median(np.diff(t)) * 60:.1f} min)")
    return t, y


def _periodogram(t: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    periods = np.linspace(PERIOD_MIN_H, PERIOD_MAX_H, N_PERIOD_SAMPLES)
    power = LombScargle(t, y).power(1.0 / periods)
    return periods, power


def main() -> None:
    plt.rcdefaults()
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })

    t_da, y_da, nyquist_h = _load_da_july()
    t_ej2, y_ej2 = _load_ej2_july()

    p_da, pow_da = _periodogram(t_da, y_da)
    p_ej2, pow_ej2 = _periodogram(t_ej2, y_ej2)

    peak_da = p_da[np.argmax(pow_da)]
    peak_ej2 = p_ej2[np.argmax(pow_ej2)]
    print(f"Delta-a dominant peak: P={peak_da:.3f} h, power={pow_da.max():.4f}")
    print(f"E_J2 (detrended) dominant peak in [{PERIOD_MIN_H},{PERIOD_MAX_H}] h: "
          f"P={peak_ej2:.3f} h, power={pow_ej2.max():.4f}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

    ax1.plot(p_da, pow_da, color="tab:blue", lw=1.3,
              label=r"Raw residual $\Delta a$ (Fig. 1a, green curve)")
    ax1.axvspan(PERIOD_MIN_H, nyquist_h, color="0.7", alpha=0.15, zorder=0)
    ax1.set_ylabel("Lomb-Scargle power")
    ax1.set_title("(a) Raw altitude residual $\\Delta a$ (July 2024, TLE-based)")
    ax1.legend(loc="upper right")

    ax2.plot(p_ej2, pow_ej2, color="tab:orange", lw=1.3,
              label=r"TLE-independent proxy $E_{J2}$ (linearly detrended)")
    ax2.set_xlabel("Period [h]")
    ax2.set_ylabel("Lomb-Scargle power")
    ax2.set_title("(b) TLE-independent specific-energy proxy $E_{J2}$ (July 2024, GNSS-only)")
    ax2.legend(loc="upper right")

    for ax in (ax1, ax2):
        ax.axvline(MARK_PERIOD_H, color="0.5", linestyle="--", linewidth=1.1, zorder=1)
        ax.grid(True, linestyle="--", alpha=0.3)
        ax.set_xlim(PERIOD_MIN_H, PERIOD_MAX_H)

    fig.tight_layout()

    for ext in ["pdf", "png"]:
        out = OUT / f"fig1_lomb_scargle.{ext}"
        fig.savefig(out, format=ext, bbox_inches="tight", dpi=150 if ext == "png" else None)
        print(f"wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
