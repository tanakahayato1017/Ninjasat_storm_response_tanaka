"""Validate the NinjaSat storm energy-loss proxy against independent
TU Delft Swarm A/B and GRACE-FO thermospheric densities.

This supplies the density-excess scaling inputs used in manuscript Methods
and Results Sect. 2.7. The numerical procedure is unchanged from the
pre-revision analysis: native densities are averaged into three-hour UTC
bins, then compared at zero lag with the raw and five-point SG-smoothed
NinjaSat energy-loss-rate proxy.

Inputs (repository-relative):
  data/external/tudelft/{SA,SB}_DNS_POD_2024_{05,10}_v02.zip
  data/external/tudelft/GC_DNS_ACC_2024_{05,10}_v02c.zip
  data/derived/energy/dEdt_storm_3h_redacted.csv
      The public rate-only copy; absolute orbital energy was removed.

Outputs:
  data/external/tudelft/tudelft_density_3h_storm.csv
  data/figure_data/swarm_validation.csv
  data/figure_data/swarm_validation_summary.csv

Run (from repository root):
  python src/08_swarm_density_validation.py
"""

import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
TUDELFT_DIR = ROOT / "data" / "external" / "tudelft"
IN_DEDT = ROOT / "data" / "derived" / "energy" / "dEdt_storm_3h_redacted.csv"
OUT_DENS = TUDELFT_DIR / "tudelft_density_3h_storm.csv"
OUT_TS = ROOT / "data" / "figure_data" / "swarm_validation.csv"
OUT_SUMMARY = ROOT / "data" / "figure_data" / "swarm_validation_summary.csv"

MONTHS = ["05", "10"]
POD_COLS = ["date", "time", "tsys", "alt_m", "lon", "lat", "lst", "arglat", "rho"]
ACC_COLS = ["date", "time", "tsys", "alt_m", "lon", "lat", "lst", "arglat",
            "rho_x", "rho_mean_orbit", "flag_x", "flag_mean"]
SATELLITES = {
    "Swarm_A": {"prefix": "SA_DNS_POD", "suffix": "v02", "kind": "pod"},
    "Swarm_B": {"prefix": "SB_DNS_POD", "suffix": "v02", "kind": "pod"},
    "GRACE-FO": {"prefix": "GC_DNS_ACC", "suffix": "v02c", "kind": "acc"},
}

EVENTS = {
    "May": (pd.Timestamp("2024-05-01"), pd.Timestamp("2024-06-01")),
    "Oct": (pd.Timestamp("2024-10-01"), pd.Timestamp("2024-11-01")),
}
STORM_WINDOWS = {
    "May": (pd.Timestamp("2024-05-08"), pd.Timestamp("2024-05-16")),
    "Oct": (pd.Timestamp("2024-10-08"), pd.Timestamp("2024-10-16")),
}
QUIET_WINDOWS = {
    "May": (pd.Timestamp("2024-05-05"), pd.Timestamp("2024-05-10")),
    "Oct": (pd.Timestamp("2024-10-05"), pd.Timestamp("2024-10-10")),
}

SG_POINTS = 5
SG_POLYORDER = 1
MIN_N = 10

def read_zip_month(path: Path, kind: str) -> pd.DataFrame:

    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        assert len(names) == 1, f"expected 1 file inside {path}, got {names}"
        with z.open(names[0]) as f:
            raw = f.read().decode("utf-8", errors="replace")
    lines = [ln for ln in raw.splitlines() if ln and not ln.startswith("#")]
    cols = POD_COLS if kind == "pod" else ACC_COLS
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=r"\s+", header=None,
                     names=cols)
    df["datetime"] = pd.to_datetime(df["date"] + " " + df["time"],
                                    utc=True, errors="coerce")
    df = df.dropna(subset=["datetime"])
    if kind == "pod":
        out = df[["datetime", "alt_m", "rho"]].copy()
    else:
        n_before = len(df)
        df = df[df["flag_x"] == 0.0].copy()
        if n_before - len(df):
            print(f"    {path.name}: dropped {n_before - len(df)}/{n_before} "
                  f"anomalous (flag_x!=0) records")
        out = df[["datetime", "alt_m", "rho_x"]].rename(columns={"rho_x": "rho"})
    return out

def sg_smooth_segments(values: np.ndarray, window_pts: int) -> np.ndarray:

    if window_pts <= 1:
        return values.copy()
    out = values.copy()
    idx = np.where(np.isfinite(values))[0]
    if len(idx) == 0:
        return out
    breaks = np.where(np.diff(idx) > 1)[0]
    seg_starts = np.concatenate([[0], breaks + 1])
    seg_ends = np.concatenate([breaks, [len(idx) - 1]])
    for s, e in zip(seg_starts, seg_ends):
        lo, hi = idx[s], idx[e] + 1
        if hi - lo >= window_pts:
            out[lo:hi] = savgol_filter(values[lo:hi], window_pts, SG_POLYORDER,
                                       mode="interp")
    return out

def corr_pair(x: np.ndarray, y: np.ndarray):
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < MIN_N:
        return np.nan, np.nan, n
    xm, ym = x[mask], y[mask]
    sp = float(np.corrcoef(rankdata(xm), rankdata(ym))[0, 1])
    pe = float(np.corrcoef(xm, ym)[0, 1])
    return sp, pe, n

def main() -> None:

    dens_bins = []
    for name, spec in SATELLITES.items():
        print(f"=== {name} ({spec['kind']}) ===")
        frames = []
        for mm in MONTHS:
            path = TUDELFT_DIR / f"{spec['prefix']}_2024_{mm}_{spec['suffix']}.zip"
            df = read_zip_month(path, spec["kind"])
            print(f"    loaded {path.name}: {len(df)} rows "
                  f"({df['datetime'].min()} .. {df['datetime'].max()})")
            frames.append(df)
        full = pd.concat(frames, ignore_index=True).sort_values("datetime")
        b = full.set_index("datetime").resample("3h", label="left",
                                                 closed="left").agg(
            rho_mean=("rho", "mean"),
            rho_median=("rho", "median"),
            alt_mean=("alt_m", "mean"),
            n_samples=("rho", "count"),
        )
        b = b[b["n_samples"] > 0].reset_index()
        b["satellite"] = name

        b["datetime"] = b["datetime"].dt.tz_localize(None)
        dens_bins.append(b)
        print(f"    3h bins: {len(b)}, alt_mean {b['alt_mean'].mean()/1000:.1f} km, "
              f"rho median {b['rho_mean'].median():.3e} kg/m3")

    dens = pd.concat(dens_bins, ignore_index=True)[
        ["datetime", "satellite", "rho_mean", "rho_median", "alt_mean",
         "n_samples"]]
    OUT_DENS.parent.mkdir(parents=True, exist_ok=True)
    dens.to_csv(OUT_DENS, index=False, encoding="utf-8")
    print(f"\nwrote {OUT_DENS} ({len(dens)} rows)")

    dedt = pd.read_csv(IN_DEDT, parse_dates=["datetime"])

    ts_rows = []
    summary_rows = []
    for event, (start, end) in EVENTS.items():
        grid = pd.date_range(start, end, freq="3h", inclusive="left")
        prox = dedt[dedt["event"] == event].set_index("datetime")
        proxy_raw = (-prox["dEdt_Jkg_per_day"]).reindex(grid).to_numpy(float)
        proxy_sg5 = sg_smooth_segments(proxy_raw, SG_POINTS)

        wide = pd.DataFrame({"datetime": grid, "event": event})
        for name in SATELLITES:
            col = f"rho_{name.replace('_', '').replace('-', '')}"
            s = dens[dens["satellite"] == name].set_index("datetime")["rho_mean"]
            wide[col] = s.reindex(grid).to_numpy(float)
        wide["proxy_raw"] = proxy_raw
        wide["proxy_sg5"] = proxy_sg5
        ts_rows.append(wide)

        ws, we = STORM_WINDOWS[event]
        qs, qe = QUIET_WINDOWS[event]
        storm_mask = (grid >= ws) & (grid < we)
        quiet_mask = (grid >= qs) & (grid < qe)

        rho_cols = [c for c in wide.columns if c.startswith("rho_")]
        for pcol in ["proxy_raw", "proxy_sg5"]:
            pv = wide[pcol].to_numpy(float)
            for rcol in rho_cols:
                rv = wide[rcol].to_numpy(float)
                sp, pe, n = corr_pair(rv[storm_mask], pv[storm_mask])
                summary_rows.append({
                    "event": event, "series": rcol, "corr_target": pcol,
                    "spearman": sp, "pearson": pe, "n": n,
                    "quiet_mean": np.nan, "storm_peak": np.nan,
                    "enhancement_ratio": np.nan,
                })
                print(f"[{event}] {rcol:12s} vs {pcol:9s} (lag 0, storm window): "
                      f"sp={sp:+.3f} pe={pe:+.3f} n={n}")

        for col in rho_cols + ["proxy_raw", "proxy_sg5"]:
            v = wide[col].to_numpy(float)
            quiet_mean = np.nanmean(v[quiet_mask])
            storm_peak = np.nanmax(v[storm_mask]) if np.any(
                np.isfinite(v[storm_mask])) else np.nan
            ratio = storm_peak / quiet_mean if (np.isfinite(quiet_mean)
                                                and quiet_mean > 0) else np.nan
            summary_rows.append({
                "event": event, "series": col, "corr_target": "enhancement",
                "spearman": np.nan, "pearson": np.nan,
                "n": int(np.isfinite(v[storm_mask]).sum()),
                "quiet_mean": quiet_mean, "storm_peak": storm_peak,
                "enhancement_ratio": ratio,
            })
            print(f"[{event}] {col:12s} enhancement: quiet={quiet_mean:.3e} "
                  f"peak={storm_peak:.3e} ratio={ratio:.1f}")

    ts = pd.concat(ts_rows, ignore_index=True)
    ts.to_csv(OUT_TS, index=False, encoding="utf-8")
    print(f"\nwrote {OUT_TS} ({len(ts)} rows)")

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUT_SUMMARY, index=False, encoding="utf-8")
    print(f"wrote {OUT_SUMMARY} ({len(summary)} rows)")

if __name__ == "__main__":
    main()
