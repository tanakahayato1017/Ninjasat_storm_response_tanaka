"""Fetch one-minute AE and AL indices from WDC Kyoto and validate their
hourly means against OMNI2.

Each day prefers WDC Kyoto provisional data and falls back to the QUICKLK
real-time product when provisional data are unavailable. The vendored files
used by the offline manuscript analysis retain that provisional-versus-
real-time caveat; rerunning this fetcher later may therefore yield revised
values.

Input (repository-relative):
  data/external/omni2/omni2_2024.dat

Outputs:
  data/external/wdc_kyoto/ae_al_1min_{2405,2410}.csv

Run (from repository root; network access required):
  python src/05a_ae_al_fetch.py
"""

from pathlib import Path
from urllib.error import HTTPError
from urllib.request import urlopen

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "data" / "external" / "wdc_kyoto"
OMNI = ROOT / "data" / "external" / "omni2" / "omni2_2024.dat"

BASE_PROV = "https://wdc.kugi.kyoto-u.ac.jp/ae_provisional/{ym}/{p}{ymd}.for.request"
BASE_RT = "https://wdc.kugi.kyoto-u.ac.jp/ae_realtime/data_dir/{yyyy}/{mm}/{dd}/{p}{ymd}"

EVENTS = {
    "2405": (pd.Timestamp("2024-05-06"), pd.Timestamp("2024-05-20")),
    "2410": (pd.Timestamp("2024-10-05"), pd.Timestamp("2024-10-19")),
}
FILL = {9999, 99999}

def _fetch(url: str) -> str | None:
    try:
        with urlopen(url, timeout=30) as r:
            return r.read().decode("ascii", errors="replace")
    except HTTPError as e:
        if e.code == 404:
            return None
        raise

def _fetch_day(param: str, day: pd.Timestamp) -> str:

    ymd = f"{day:%y%m%d}"
    ym = f"{day:%Y%m}"
    url_prov = BASE_PROV.format(ym=ym, p=param, ymd=ymd)
    txt = _fetch(url_prov)
    if txt is not None:
        return txt, "provisional"
    url_rt = BASE_RT.format(yyyy=f"{day:%Y}", mm=f"{day:%m}", dd=f"{day:%d}",
                             p=param, ymd=ymd)
    txt = _fetch(url_rt)
    if txt is not None:
        return txt, "real-time"
    raise RuntimeError(f"no data for {param} {day:%Y-%m-%d} "
                        f"(tried {url_prov} and {url_rt})")

def _parse_day(text: str, day: pd.Timestamp) -> pd.DataFrame:

    rows = []
    for line in text.splitlines():
        toks = line.split()
        if len(toks) < 3 or not toks[1][:6].isdigit():
            continue
        code = toks[1]
        hh = int(code[7:9])
        vals = toks[3:63]
        if len(vals) != 60:
            raise ValueError(f"expected 60 minute values, got {len(vals)}: {line!r}")
        for mm, v in enumerate(vals):
            iv = int(v)
            rows.append((day + pd.Timedelta(hours=hh, minutes=mm),
                         np.nan if iv in FILL else float(iv)))
    return pd.DataFrame(rows, columns=["datetime", "value"])

def fetch_event(tag: str, t0: pd.Timestamp, t1: pd.Timestamp) -> pd.DataFrame:
    days = pd.date_range(t0, t1, freq="D", inclusive="left")
    out = {"AE": [], "AL": []}
    sources = {"AE": set(), "AL": set()}
    for day in days:
        for param, key in [("ae", "AE"), ("al", "AL")]:
            text, src = _fetch_day(param, day)
            out[key].append(_parse_day(text, day))
            sources[key].add(src)
        print(f"  [{tag}] {day:%Y-%m-%d}  OK")
    df = pd.merge(
        pd.concat(out["AE"], ignore_index=True).rename(columns={"value": "AE"}),
        pd.concat(out["AL"], ignore_index=True).rename(columns={"value": "AL"}),
        on="datetime", how="outer",
    ).sort_values("datetime").reset_index(drop=True)
    print(f"  [{tag}] sources used: AE={sorted(sources['AE'])} "
          f"AL={sorted(sources['AL'])}")
    print(f"  [{tag}] n={len(df)}, missing AE={df['AE'].isna().sum()}, "
          f"missing AL={df['AL'].isna().sum()}")
    return df

def load_omni_hourly() -> pd.DataFrame:
    raw = np.loadtxt(OMNI)
    t = (pd.to_datetime(raw[:, 0].astype(int).astype(str), format="%Y")
         + pd.to_timedelta(raw[:, 1].astype(int) - 1, unit="D")
         + pd.to_timedelta(raw[:, 2].astype(int), unit="h"))
    ae = raw[:, 41]
    al = raw[:, 52]
    ae[ae >= 9999] = np.nan
    al[al >= 99999] = np.nan
    return pd.DataFrame({"datetime": t, "AE_omni2": ae, "AL_omni2": al})

def validate_against_omni2(df_1min: pd.DataFrame, tag: str) -> None:

    hourly = (df_1min.set_index("datetime")[["AE", "AL"]]
              .resample("1h").mean().reset_index())
    omni = load_omni_hourly()
    m = hourly.merge(omni, on="datetime", how="inner")
    for col, ocol in [("AE", "AE_omni2"), ("AL", "AL_omni2")]:
        d = (m[col] - m[ocol]).dropna()
        if len(d) == 0:
            print(f"  [{tag}] {col}: no overlapping hours with OMNI2")
            continue
        print(f"  [{tag}] {col} vs OMNI2: n={len(d)}, "
              f"mean diff={d.mean():+.1f} nT, "
              f"|diff| median={d.abs().median():.1f} nT, "
              f"max|diff|={d.abs().max():.1f} nT, "
              f"corr={m[[col, ocol]].corr().iloc[0, 1]:.4f}")

def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for tag, (t0, t1) in EVENTS.items():
        print("=" * 78)
        print(f"{tag}: fetching {t0:%Y-%m-%d} .. {t1:%Y-%m-%d} (exclusive)")
        print("=" * 78)
        df = fetch_event(tag, t0, t1)
        out_path = OUT_DIR / f"ae_al_1min_{tag}.csv"
        df.to_csv(out_path, index=False, encoding="utf-8")
        print(f"  wrote {out_path}")
        validate_against_omni2(df, tag)

if __name__ == "__main__":
    main()
