"""XRISM/NinjaSat storm-response detection significance and the October
non-detection's 2-sigma upper limit and own-May-response scaling
prediction (manuscript Results Sect. 2.7, Table 1).

Addresses the question of whether the October XRISM non-detection is
physically informative or merely a sensitivity limit (Reviewer 2, Comment
R2-1): quantifies each satellite's detection significance in both storms,
and predicts the October response from each satellite's own May response,
scaled by the independently measured Swarm/GRACE-FO density ratio and each
satellite's own altitude change between the two events (which cancels the
poorly known ballistic coefficient, since each satellite is compared only
to itself).

Method (identical to the manuscript's fiducial pipeline):
  Delta a = a_SGP4 - a_GNSS  [km]
    -> Savitzky-Golay (polyorder 1, window 9 points) = Delta a-tilde
    -> centred difference                            = Delta a-dot [km/h]
  Quiet window (SSC - 5 d .. SSC): mean and standard deviation of
  Delta a-dot. Storm window (SSC .. SSC + 4 d): peak excess above the
  quiet mean, in units of the quiet-window standard deviation.

Input (data/raw/, pending NinjaSat/XRISM team release -- see the READMEs
in data/raw/ninjasat/ and data/raw/xrism/):
  data/raw/ninjasat/gps_tle_average_altitude_fitting_all_severe_{2405,2410}
    _monthly_ext.csv
  data/raw/xrism/orbital_altitude_average_fitting_{2405,2410}_3period.csv

Output:
  data/figure_data/xrism_detection_limit.csv
  (a copy of this file, matching Table 1 of the manuscript, is already
  published so the reported results can be inspected before the raw
  inputs above are released)

Run (from repository root):
  python src/06_xrism_detection_scaling.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "figure_data" / "xrism_detection_limit.csv"

NINJA = ROOT / "data" / "raw" / "ninjasat" / "gps_tle_average_altitude_fitting_all_severe_{}_monthly_ext.csv"
XRISM = ROOT / "data" / "raw" / "xrism" / "orbital_altitude_average_fitting_{}_3period.csv"

SSC = {"May": pd.Timestamp("2024-05-10 17:05"),
       "Oct": pd.Timestamp("2024-10-10 15:04")}
TAG = {"May": "2405", "Oct": "2410"}

QUIET_DAYS = 5.0    # SSC - 5 d .. SSC
STORM_DAYS = 4.0    # SSC .. SSC + 4 d
SG_WIN, SG_POLY = 9, 1

# Swarm A/B + GRACE-FO quiet-mean and storm-peak densities [kg/m^3]
# (see src/08_swarm_density_validation.py, which produces these numbers
# from the TU Delft thermospheric density product)
SWARM = {
    "SwarmA":  {"May": (9.978807270645832e-13, 5.044782215833333e-12),
                "Oct": (2.796988387460417e-12, 5.155603875e-12)},
    "SwarmB":  {"May": (6.623776158666666e-13, 3.771749067777777e-12),
                "Oct": (1.7768463510291667e-12, 3.4913856766666664e-12)},
    "GRACEFO": {"May": (9.70338008272454e-13, 4.393251107685185e-12),
                "Oct": (2.3539567048918984e-12, 4.3416376565740736e-12)},
}


def load(path_tmpl: Path, tag: str) -> pd.DataFrame:
    df = pd.read_csv(str(path_tmpl).format(tag))
    tcol = "timestamp" if "timestamp" in df.columns else "time_average"
    acol = ("altitude_diff" if "altitude_diff" in df.columns
            else "altitude_difference")
    gcol = ("gps_average_altitude" if "gps_average_altitude" in df.columns
            else "average_altitude_gps")
    t = pd.to_datetime(df[tcol], utc=True, format="mixed").dt.tz_localize(None)
    return (df.assign(dt=t, da=df[acol].astype(float), alt=df[gcol].astype(float))
              .sort_values("dt").reset_index(drop=True)[["dt", "da", "alt"]])


def rate_km_per_h(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """SG-smoothed Delta a-tilde and its centred difference Delta a-dot [km/h]."""
    sg = savgol_filter(df["da"].to_numpy(float), SG_WIN, SG_POLY, mode="interp")
    tsec = (df["dt"] - df["dt"].iloc[0]).dt.total_seconds().to_numpy()
    d = np.full(len(sg), np.nan)
    d[1:-1] = (sg[2:] - sg[:-2]) / (tsec[2:] - tsec[:-2]) * 3600.0
    return sg, d


def main() -> None:
    rows = []

    def add(**kw):
        rows.append(kw)

    print("=" * 78)
    print("XRISM / NinjaSat detection significance and the October upper limit")
    print("=" * 78)

    peaks = {}
    for sat, path in [("NinjaSat", NINJA), ("XRISM", XRISM)]:
        for ev in ["May", "Oct"]:
            df = load(path, TAG[ev])
            sg, dot = rate_km_per_h(df)
            df = df.assign(sg=sg, dot=dot)

            ssc = SSC[ev]
            q = df[(df["dt"] >= ssc - pd.Timedelta(days=QUIET_DAYS)) & (df["dt"] < ssc)]
            s = df[(df["dt"] >= ssc) & (df["dt"] <= ssc + pd.Timedelta(days=STORM_DAYS))]

            cad = np.median(np.diff(df["dt"].to_numpy())
                            .astype("timedelta64[s]").astype(float)) / 3600.0
            qm = np.nanmean(q["dot"])
            qs = np.nanstd(q["dot"], ddof=1)
            pk = np.nanmax(s["dot"])
            excess = pk - qm
            nsig = excess / qs
            ul2_exc = 2 * qs  # 2-sigma upper limit on the excess
            alt = df["alt"].mean()

            peaks[(sat, ev)] = dict(peak=pk, quiet_mean=qm, quiet_sd=qs,
                                    excess=excess, nsig=nsig, ul2_exc=ul2_exc,
                                    alt=alt)

            print(f"\n[{sat} {ev}]  n={len(df)}, cadence={cad:.2f} h, "
                  f"mean alt={alt:.1f} km")
            print(f"   quiet ({QUIET_DAYS:.0f} d): mean={qm:+.5f}  sd={qs:.5f} km/h  (n={len(q)})")
            print(f"   storm peak = {pk:+.5f} km/h  ->  excess = {excess:+.5f} "
                  f"= {nsig:.1f} sigma")
            print(f"   2-sigma upper limit on excess = {ul2_exc:.5f} km/h")

            add(satellite=sat, event=ev, quantity="mean_altitude_km", value=alt)
            add(satellite=sat, event=ev, quantity="cadence_h", value=cad)
            add(satellite=sat, event=ev, quantity="quiet_mean_dot_km_per_h", value=qm)
            add(satellite=sat, event=ev, quantity="quiet_sd_dot_km_per_h", value=qs)
            add(satellite=sat, event=ev, quantity="storm_peak_dot_km_per_h", value=pk)
            add(satellite=sat, event=ev, quantity="storm_excess_km_per_h", value=excess)
            add(satellite=sat, event=ev, quantity="significance_sigma", value=nsig)
            add(satellite=sat, event=ev, quantity="upper_limit_2sigma_km_per_h",
                value=ul2_exc)

    # -------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("Swarm / GRACE-FO absolute density-excess ratio (Oct / May)")
    print("=" * 78)
    ratios = []
    for sat, d in SWARM.items():
        dmay = d["May"][1] - d["May"][0]
        doct = d["Oct"][1] - d["Oct"][0]
        r = doct / dmay
        ratios.append(r)
        print(f"  {sat:8s} d_rho May={dmay:.3e}  Oct={doct:.3e}  Oct/May={r:.3f}")
        add(satellite=sat, event="both", quantity="drho_ratio_oct_over_may", value=r)
    rho_ratio = float(np.mean(ratios))
    print(f"\n  mean Oct/May absolute density-excess ratio = {rho_ratio:.3f}")
    add(satellite="Swarm+GRACEFO", event="both",
        quantity="mean_drho_ratio_oct_over_may", value=rho_ratio)

    # -------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("October response predicted by scaling each satellite's own May response")
    print("=" * 78)
    print("  Referencing each satellite to itself cancels its ballistic coefficient")
    print("  C_D A/m, so the prediction needs only the Oct/May density ratio and")
    print("  each satellite's own altitude change between the two events.\n")

    for sat in ["NinjaSat", "XRISM"]:
        may, oct_ = peaks[(sat, "May")], peaks[(sat, "Oct")]
        dh = may["alt"] - oct_["alt"]  # positive = lower in October
        print(f"  [{sat}] altitude May {may['alt']:.1f} -> Oct {oct_['alt']:.1f} km "
              f"(diff {dh:+.1f} km)")
        for H in [50.0, 60.0, 70.0]:
            boost = float(np.exp(dh / H))
            pred = may["excess"] * rho_ratio * boost
            obs = oct_["excess"]
            print(f"     H={H:4.0f} km: altitude factor x{boost:.2f}  "
                  f"-> predicted Oct excess = {pred:+.5f} km/h,  "
                  f"observed = {obs:+.5f} ({obs / pred:.2f} x prediction)")
            add(satellite=sat, event="Oct",
                quantity=f"predicted_excess_km_per_h_H{H:.0f}", value=pred)
            add(satellite=sat, event="Oct",
                quantity=f"obs_over_pred_H{H:.0f}", value=obs / pred)
        print()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, index=False, encoding="utf-8")
    print(f"wrote {OUT} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
