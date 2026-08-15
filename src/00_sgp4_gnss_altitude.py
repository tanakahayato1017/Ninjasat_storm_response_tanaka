"""Regenerate the per-epoch GNSS and SGP4 altitude intermediates.

The script ports the published NinjaSat monthly workflow and the surviving
XRISM batch workflow from ``orbit_gps_tle_debug.ipynb`` and ``xrism.ipynb``.
The TLE catalogue row is fixed once per month, exactly as in those notebooks.

Inputs (relative to the repository root)
----------------------------------------
NinjaSat raw fixes:
    data/raw/GNSS_data/conbined_gps_{2405,2407,2410}_ext.csv
XRISM raw fixes:
    data/raw/xrism_gnss/xrism_gps_{2405,2410}.csv
TLE catalogues (local only pending redistribution review):
    data/external/tle/ninjasat_tle_per_day.csv
    data/external/tle/xrism_tle_each_day.csv

Outputs
-------
All regenerated files are written beneath ``data/derived/intermediate``.

Notes
-----
NinjaSat uses geodetic fixes transformed from EPSG:4326 to EPSG:4978 and
subtracts the WGS-84 equatorial radius (6378.137 km) from both vector norms.
The verified XRISM reference files were produced by the notebook's batch cell:
they use the supplied Cartesian GNSS coordinates, propagate at 60-second
epochs, merge to the nearest raw fix, and subtract 6371 km.  This historical
XRISM convention is retained because numerical fidelity to the published
intermediates is the governing requirement.

Examples
--------
Run every staged month::

    python src/00_sgp4_gnss_altitude.py --satellite all

Run one month::

    python src/00_sgp4_gnss_altitude.py --satellite ninjasat --months 2405
"""

from __future__ import annotations

import argparse
from datetime import timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from pyproj import Transformer
from sgp4.api import Satrec, jday


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "derived" / "intermediate"
EARTH_EQUATORIAL_RADIUS_KM = 6378.137
XRISM_LEGACY_RADIUS_KM = 6371.0

# Verified against the pre-migration monthly files.  The epochs are copied
# from the catalogue's EPOCH column; all are pandas' zero-based row indices.
NINJASAT_TLE_ROWS = {
    2405: (170, "2024-05-16 02:49 UTC"),
    2407: (230, "2024-07-15 04:22 UTC"),
    2410: (321, "2024-10-16 04:10 UTC"),
}
XRISM_TLE_ROWS = {
    2405: (234, "2024-05-01 13:01 UTC"),
    2410: (386, "2024-10-01 11:46 UTC"),
}


def _sgp4_position(satellite: Satrec, timestamp: pd.Timestamp) -> tuple[float, float, float]:
    """Return the TEME position in kilometres at integer-second precision."""
    dt = timestamp.to_pydatetime()
    jd, fr = jday(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    error, position, _ = satellite.sgp4(jd, fr)
    if error != 0:
        raise RuntimeError(f"SGP4 failed at {timestamp} with error code {error}")
    return position


def _eci_to_ecef(position: tuple[float, float, float], timestamp: pd.Timestamp) -> tuple[float, float, float]:
    """Apply the notebook's Greenwich-sidereal-time rotation."""
    dt = timestamp.to_pydatetime()
    jd, fr = jday(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second)
    centuries = (jd - 2451545.0) / 36525.0
    theta_deg = (
        280.46061837
        + 360.98564736629 * (jd + fr - 2451545.0)
        + 0.000387933 * centuries**2
        - centuries**3 / 38710000.0
    )
    theta = np.radians(theta_deg % 360.0)
    x_eci, y_eci, z_eci = position
    return (
        x_eci * np.cos(theta) + y_eci * np.sin(theta),
        -x_eci * np.sin(theta) + y_eci * np.cos(theta),
        z_eci,
    )


def process_ninjasat(month: int) -> Path:
    """Regenerate one NinjaSat monthly per-epoch intermediate."""
    row, epoch = NINJASAT_TLE_ROWS[month]
    raw_path = ROOT / "data" / "raw" / "GNSS_data" / f"conbined_gps_{month}_ext.csv"
    tle_path = ROOT / "data" / "external" / "tle" / "ninjasat_tle_per_day.csv"
    raw = pd.read_csv(raw_path)
    tle = pd.read_csv(tle_path)

    # Keep pyproj's authority-axis order: EPSG:4326 expects latitude then
    # longitude, matching the original notebook call.
    to_ecef = Transformer.from_crs("epsg:4326", "epsg:4978")
    x, y, z = to_ecef.transform(
        raw["gpsLatitude"], raw["gpsLongitude"], raw["gpsAltitudeMeters"]
    )
    x, y, z = np.asarray(x) / 1000.0, np.asarray(y) / 1000.0, np.asarray(z) / 1000.0
    gps_altitude = np.sqrt(x**2 + y**2 + z**2) - EARTH_EQUATORIAL_RADIUS_KM
    times = pd.to_datetime(raw["timestamp"], unit="s", utc=True)

    satellite = Satrec.twoline2rv(tle.loc[row, "TLE_LINE1"], tle.loc[row, "TLE_LINE2"])
    to_geodetic = Transformer.from_crs("EPSG:4978", "EPSG:4326")
    propagated = []
    for timestamp in times:
        position = _sgp4_position(satellite, timestamp)
        ecef = _eci_to_ecef(position, timestamp)
        lat, lon, alt = to_geodetic.transform(*(value * 1000.0 for value in ecef))
        propagated.append(
            (
                timestamp.strftime("%Y-%m-%dT%H:%M:%SZ"),
                np.linalg.norm(position) - EARTH_EQUATORIAL_RADIUS_KM,
                *ecef,
                lat,
                lon,
                alt,
            )
        )

    output = pd.DataFrame(
        {
            "gpsUtcTime": times,
            "gps_altitude": gps_altitude,
            "gps_x": x,
            "gps_y": y,
            "gps_z": z,
            "gps_lat": raw["gpsLatitude"],
            "gps_lon": raw["gpsLongitude"],
            "gps_alt": raw["gpsAltitudeMeters"],
        }
    )
    propagated_df = pd.DataFrame(
        propagated,
        columns=["time", "tle_altitude", "tle_x", "tle_y", "tle_z", "tle_lat", "tle_lon", "tle_alt"],
    )
    output = pd.concat([output, propagated_df], axis=1)
    output["tle_gps_diff"] = output["tle_altitude"] - output["gps_altitude"]
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"ninjasat_gps_and_tle_{month}_ext.csv"
    output.to_csv(path, index=False)
    print(f"NinjaSat {month}: TLE row {row} ({epoch}); wrote {path} ({len(output)} rows)")
    return path


def process_xrism(month: int) -> Path:
    """Regenerate one XRISM monthly per-epoch merged intermediate."""
    row, epoch = XRISM_TLE_ROWS[month]
    raw_path = ROOT / "data" / "raw" / "xrism_gnss" / f"xrism_gps_{month}.csv"
    tle_path = ROOT / "data" / "external" / "tle" / "xrism_tle_each_day.csv"
    raw = pd.read_csv(raw_path)
    tle = pd.read_csv(tle_path)
    raw["datetime"] = pd.to_datetime(raw["TIME (UTC)"], errors="coerce")

    satellite = Satrec.twoline2rv(tle.loc[row, "TLE_LINE1"], tle.loc[row, "TLE_LINE2"])
    current = raw["datetime"].min()
    end = raw["datetime"].max()
    propagated = []
    while current <= end:
        position = _sgp4_position(satellite, current)
        propagated.append(
            {
                "datetime": current,
                "x": position[0],
                "y": position[1],
                "z": position[2],
                "tle_altitude": np.linalg.norm(position) - XRISM_LEGACY_RADIUS_KM,
            }
        )
        current += timedelta(seconds=60)
    propagated_df = pd.DataFrame(propagated)

    raw = raw.rename(columns={"X (km)": "gps_x", "Y (km)": "gps_y", "Z (km)": "gps_z"})
    raw["gps_altitude"] = (
        np.sqrt(raw["gps_x"] ** 2 + raw["gps_y"] ** 2 + raw["gps_z"] ** 2)
        - XRISM_LEGACY_RADIUS_KM
    )
    output = pd.merge_asof(
        propagated_df.sort_values("datetime"),
        raw.sort_values("datetime"),
        on="datetime",
        direction="nearest",
        tolerance=pd.Timedelta(seconds=60),
    )
    output["position_error"] = np.sqrt(
        (output["x"] - output["gps_x"]) ** 2
        + (output["y"] - output["gps_y"]) ** 2
        + (output["z"] - output["gps_z"]) ** 2
    )
    output["altitude_difference"] = output["gps_altitude"] - output["tle_altitude"]
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"xrism_gps_and_tle_{month}.csv"
    output.to_csv(path, index=False)
    print(f"XRISM {month}: TLE row {row} ({epoch}); wrote {path} ({len(output)} rows)")
    return path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--satellite", choices=("ninjasat", "xrism", "all"), default="all")
    parser.add_argument("--months", type=int, nargs="+", help="YYMM month(s); defaults to all staged months")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.satellite in ("ninjasat", "all"):
        months = args.months if args.months is not None else sorted(NINJASAT_TLE_ROWS)
        for month in months:
            if month in NINJASAT_TLE_ROWS:
                process_ninjasat(month)
            elif args.satellite == "ninjasat":
                raise ValueError(f"No NinjaSat TLE mapping for {month}")
    if args.satellite in ("xrism", "all"):
        months = args.months if args.months is not None else sorted(XRISM_TLE_ROWS)
        for month in months:
            if month in XRISM_TLE_ROWS:
                process_xrism(month)
            elif args.satellite == "xrism":
                raise ValueError(f"No XRISM TLE mapping for {month}")


if __name__ == "__main__":
    main()
