"""Reproduce the AE/AL storm summaries and lag correlations reported in
manuscript Results Sect. 2.6 and used for Figs. 4e,f and 5e,f.

Inputs (repository-relative):
  data/external/wdc_kyoto/ae_al_1min_{2405,2410}.csv
      Vendored one-minute public WDC Kyoto AE/AL indices; see its README
      for the provisional-versus-real-time caveat.
  data/external/omni2/omni2_2024.dat
      Public hourly OMNI2 PC and Dst indices.
  data/derived/correlation_inputs/{2405,2410}_PCN_correlation_thesis.csv
      The paper's n=84 correlation grids and EUV-subtracted altitude rates.

The pending NinjaSat absolute-altitude monthly files remain staged under
data/raw/ninjasat/, but this step does not read them: the safe correlation
input CSVs already contain the required derived rate series.

Outputs:
  data/figure_data/ae_summary.csv
  data/figure_data/ae_lag_correlation.csv
  data/figure_data/ae_al_timeseries_{2405,2410}.csv
  data/figure_data/ae_al_crosscorr_curve_{2405,2410}.csv

Run (from repository root; uses vendored data and requires no network):
  python src/05_ae_al_analysis.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OMNI = ROOT / "data" / "external" / "omni2" / "omni2_2024.dat"
AE_1MIN = ROOT / "data" / "external" / "wdc_kyoto" / "ae_al_1min_{}.csv"

OUT_SUM = ROOT / "data" / "figure_data" / "ae_summary.csv"
OUT_LAG = ROOT / "data" / "figure_data" / "ae_lag_correlation.csv"
OUT_TS = ROOT / "data" / "figure_data" / "ae_al_timeseries_{}.csv"
OUT_XCORR = ROOT / "data" / "figure_data" / "ae_al_crosscorr_curve_{}.csv"

COL = {"Dst": 40, "PC": 51}
FILL = {"Dst": 99999, "PC": 999.9}

EVENTS = {
    "May": ("2405", pd.Timestamp("2024-05-10 17:05"),
            pd.Timestamp("2024-05-10"), pd.Timestamp("2024-05-15")),
    "Oct": ("2410", pd.Timestamp("2024-10-10 15:04"),
            pd.Timestamp("2024-10-10"), pd.Timestamp("2024-10-15")),
}

CORR_WIN = {"May": (pd.Timestamp("2024-05-08"), pd.Timestamp("2024-05-16")),
            "Oct": (pd.Timestamp("2024-10-08"), pd.Timestamp("2024-10-16"))}

BOX_HALF_H = 1.57

def load_omni():
    raw = np.loadtxt(OMNI)
    year = raw[:, 0].astype(int)
    doy = raw[:, 1].astype(int)
    hour = raw[:, 2].astype(int)
    t = (pd.to_datetime(year.astype(str), format="%Y")
         + pd.to_timedelta(doy - 1, unit="D")
         + pd.to_timedelta(hour, unit="h"))
    out = pd.DataFrame({"datetime": t})
    for k, c in COL.items():
        v = raw[:, c].astype(float)
        v[v >= FILL[k]] = np.nan
        out[k] = v
    return out

def load_ae_al_1min(tag):
    df = pd.read_csv(Path(str(AE_1MIN).format(tag)), parse_dates=["datetime"])
    return df.set_index("datetime")[["AE", "AL"]]

def box_average(src_times, src_values, target_times, half_window_h=BOX_HALF_H):

    src_t = np.asarray(src_times, dtype="datetime64[ns]")
    src_v = np.asarray(src_values, dtype=float)
    order = np.argsort(src_t)
    src_t, src_v = src_t[order], src_v[order]
    half = np.timedelta64(int(round(half_window_h * 3600)), "s")
    out = np.full(len(target_times), np.nan)
    for i, t in enumerate(pd.DatetimeIndex(target_times).values):
        lo = np.searchsorted(src_t, t - half, side="left")
        hi = np.searchsorted(src_t, t + half, side="right")
        if hi > lo:
            out[i] = np.nanmean(src_v[lo:hi])
    return out

def fisher_ci(r, n, z=1.96):
    if not np.isfinite(r) or n <= 3:
        return np.nan, np.nan
    zt = np.arctanh(r)
    se = 1.0 / np.sqrt(n - 3)
    return float(np.tanh(zt - z * se)), float(np.tanh(zt + z * se))

def shift_corr(y, x_grid, shifts):

    y = pd.Series(y)
    x = pd.Series(x_grid)
    best = (np.nan, np.nan, 0)
    rows = []
    for k in shifts:
        y_shift = y.shift(k)
        m = y_shift.notna() & x.notna()
        n = int(m.sum())
        if n < 20:
            rows.append((k, np.nan, n))
            continue
        r = float(np.corrcoef(y_shift[m], x[m])[0, 1])
        rows.append((k, r, n))
        if not np.isfinite(best[0]) or abs(r) > abs(best[0]):
            best = (r, k, n)
    return best, rows

def main() -> None:
    ae_al_1min = {tag: load_ae_al_1min(tag) for tag in ["2405", "2410"]}
    omni = load_omni()
    print("=" * 78)
    print("Phase 2 AE/AL analysis (1-min WDC Kyoto AE/AL; hourly OMNI2 Dst/PC/AU)")
    print("=" * 78)
    for tag, df in ae_al_1min.items():
        print(f"AE/AL 1min [{tag}]: {len(df)} rows, {df.index.min()} .. {df.index.max()}, "
              f"missing AE={df['AE'].isna().sum()}, AL={df['AL'].isna().sum()}")
    print(f"OMNI2: {len(omni)} rows, {omni.datetime.min()} .. {omni.datetime.max()}")
    for k in COL:
        print(f"   {k:4s}: {omni[k].notna().sum()} valid, "
              f"range [{omni[k].min():.1f}, {omni[k].max():.1f}]")

    rows = []
    print("\n" + "-" * 78)
    print("AE/AL storm-window summary (five days from 00 UT on the SSC date)")
    print("-" * 78)
    for ev, (tag, ssc, w0, w1) in EVENTS.items():
        df = ae_al_1min[tag]
        w = df[(df.index >= w0) & (df.index < w1)]
        ae = w["AE"].to_numpy(float)
        al = w["AL"].to_numpy(float)
        ok = np.isfinite(ae)
        integ = float(np.nansum(ae)) / 60.0
        hrs1000 = float(np.nansum(ae > 1000)) / 60.0
        hrs500 = float(np.nansum(ae > 500)) / 60.0
        print(f"  [{ev}] window {w0:%m-%d} .. {w1:%m-%d}  (n={ok.sum()} min)")
        print(f"     AE max   = {np.nanmax(ae):7.0f} nT   at "
              f"{w.index[int(np.nanargmax(ae))]:%m-%d %H:%M}UT")
        print(f"     AL min   = {np.nanmin(al):7.0f} nT")
        print(f"     int AE   = {integ:7.0f} nT h")
        print(f"     hours AE>1000 = {hrs1000:5.1f} h,  AE>500 = {hrs500:5.1f} h")
        for q, v in [("AE_max_nT", np.nanmax(ae)), ("AL_min_nT", np.nanmin(al)),
                     ("integrated_AE_nT_h", integ),
                     ("hours_AE_gt_1000", hrs1000), ("hours_AE_gt_500", hrs500)]:
            rows.append(dict(event=ev, quantity=q, value=float(v)))

    may = {r["quantity"]: r["value"] for r in rows if r["event"] == "May"}
    oct_ = {r["quantity"]: r["value"] for r in rows if r["event"] == "Oct"}
    print("\n  October / May ratios:")
    for q in ["AE_max_nT", "integrated_AE_nT_h", "hours_AE_gt_1000"]:
        ratio = oct_[q] / may[q]
        print(f"     {q:22s}: {ratio:.3f}")
        rows.append(dict(event="ratio_oct_over_may", quantity=q, value=ratio))

    THESIS = ROOT / "data" / "derived" / "correlation_inputs" / "{}_PCN_correlation_thesis.csv"
    oi = omni.set_index("datetime")

    lag_rows = []
    print("\n" + "-" * 78)
    print("Lag correlations between the EUV-subtracted rate and each index")
    print("  Sign convention: Shift_Hours = shift * cadence; shifts = -8..+7")
    print("-" * 78)

    for ev, (tag, ssc, _, _) in EVENTS.items():
        th = pd.read_csv(THESIS.as_posix().format(tag))
        th["dt"] = pd.to_datetime(th["timestamp"])
        th = th.sort_values("dt").reset_index(drop=True)
        y = th["slope_difference"].to_numpy(float)
        cad = np.median(np.diff(th["dt"].to_numpy())
                        .astype("timedelta64[s]").astype(float)) / 3600.0
        shifts = range(-8, 8)
        print(f"\n  [{ev}] n={len(th)}, cadence={cad:.3f} h")

        ae_al = ae_al_1min[tag]
        grids = {
            "PCN_paper": th["PCN"].to_numpy(float),
            "AE": box_average(ae_al.index, ae_al["AE"], th["dt"]),
            "AL": box_average(ae_al.index, ae_al["AL"], th["dt"]),
            "PC": box_average(oi.index, oi["PC"], th["dt"]),
            "Dst": box_average(oi.index, oi["Dst"], th["dt"]),
        }

        xcorr_rows = []
        for name, x_grid in grids.items():
            (r, L, n), curve = shift_corr(y, x_grid, shifts)
            lag_h = L * cad if np.isfinite(L) else np.nan
            lo, hi = fisher_ci(r, n)
            flag = "  <- reproduction of the paper PC analysis" if name == "PCN_paper" else ""
            print(f"     {name:10s}: peak |r| = {r:+.3f} at shift = {L!s:>3} "
                  f"(tau = {lag_h:+.2f} h, n={n}, "
                  f"95% CI on r = [{lo:+.3f}, {hi:+.3f}]){flag}")
            lag_rows.append(dict(event=ev, index=name, peak_r=r, peak_lag_h=lag_h,
                                 n=n, r_ci_lo=lo, r_ci_hi=hi, cadence_h=cad))
            if name in ("AE", "AL"):
                for k, r_k, n_k in curve:
                    lo_k, hi_k = fisher_ci(r_k, n_k)
                    xcorr_rows.append(dict(lag_h=float(k) * cad, index=name, r=r_k,
                                            n=n_k, r_ci_lo=lo_k, r_ci_hi=hi_k))

        out_xcorr = Path(str(OUT_XCORR).format(tag))
        pd.DataFrame(xcorr_rows).to_csv(out_xcorr, index=False, encoding="utf-8")
        print(f"     wrote {out_xcorr}  ({len(xcorr_rows)} rows: "
              f"{len(list(shifts))} lags x 2 indices)")

        ts_out = pd.DataFrame({"datetime": th["dt"],
                                "slope_difference": th["slope_difference"],
                                "AE": grids["AE"], "AL": grids["AL"]})
        out_ts = Path(str(OUT_TS).format(tag))
        out_ts.parent.mkdir(parents=True, exist_ok=True)
        ts_out.to_csv(out_ts, index=False, encoding="utf-8")
        print(f"     wrote {out_ts}  (n={len(ts_out)} rows, full analysis window "
              f"{th['dt'].min():%m-%d %H:%M} .. {th['dt'].max():%m-%d %H:%M})")

    OUT_SUM.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_SUM, index=False, encoding="utf-8")
    pd.DataFrame(lag_rows).to_csv(OUT_LAG, index=False, encoding="utf-8")
    print(f"\nwrote {OUT_SUM}")
    print(f"wrote {OUT_LAG}")

if __name__ == "__main__":
    main()
