#!/usr/bin/env python3
"""Build annual daily ERA5-Land files for the phenology pipeline.

The four volumetric soil-water layers are read from the hourly ERA5-Land
archive and averaged over each UTC calendar day. Surface solar radiation,
evaporation, and potential evaporation are read from the Earth Data Hub daily
UTC product, where ERA5-Land accumulations have already been assigned to the
correct day. Hourly data are processed lazily and are never stored locally.

Native ERA5-Land soil layers
----------------------------
swvl1 : 0-7 cm [m3 m-3]
swvl2 : 7-28 cm [m3 m-3]
swvl3 : 28-100 cm [m3 m-3]
swvl4 : 100-289 cm [m3 m-3]

Output
------
One compressed NetCDF per calendar year in ``era5land_daily_yearly``.

Required packages
-----------------
    pip install xarray dask pandas numpy h5netcdf zarr
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
HOURLY_URL = (
    "https://data.earthdatahub.destine.eu/era5/"
    "reanalysis-era5-land-no-antartica-v0.zarr"
)
DAILY_URL = (
    "https://data.earthdatahub.destine.eu/era5/"
    "era5-land-daily-utc-v1.zarr"
)

START_YEAR = 1980
END_YEAR = 2025

# Germany and Austria, with a small buffer around the study locations.
LON_MIN = 5.0
LON_MAX = 17.5
LAT_MIN = 46.0
LAT_MAX = 55.8

SOIL_VARIABLES = ["swvl1", "swvl2", "swvl3", "swvl4"]
DAILY_VARIABLES = ["ssrd", "e", "pev"]
ALL_VARIABLES = DAILY_VARIABLES + SOIL_VARIABLES

# The environment override permits the same code to run on a cluster.
PROJECT_DIR = Path(
    os.environ.get("PHENOLOGY_PROJECT_DIR", Path(__file__).resolve().parents[3])
).expanduser().resolve()
OUTPUT_DIR = PROJECT_DIR / "Data" / "Era5" / "era5land_daily_yearly"
OUTPUT_TEMPLATE = "era5land_daily_{year}.nc"
OVERWRITE = False


def open_zarr_dataset(url: str) -> xr.Dataset:
    """Open a remote Earth Data Hub Zarr dataset lazily."""
    return xr.open_dataset(
        url,
        engine="zarr",
        chunks={},
        storage_options={"client_kwargs": {"trust_env": True}},
    )


def normalize_time_coordinate(ds: xr.Dataset) -> xr.Dataset:
    """Rename the supported time coordinate to ``time``."""
    for candidate in ("time", "valid_time"):
        if candidate in ds.coords or candidate in ds.dims:
            if candidate != "time":
                ds = ds.rename({candidate: "time"})
            return ds
    raise KeyError(
        "No supported time coordinate found. Expected 'time' or 'valid_time'. "
        f"Available coordinates: {list(ds.coords)}"
    )


def subset_region(ds: xr.Dataset) -> xr.Dataset:
    """Subset the configured longitude/latitude box robustly."""
    required_coords = {"longitude", "latitude"}
    missing = required_coords.difference(ds.coords)
    if missing:
        raise KeyError(f"Dataset is missing spatial coordinates: {sorted(missing)}")

    lat_first = float(ds.latitude.values[0])
    lat_last = float(ds.latitude.values[-1])
    lat_slice = (
        slice(LAT_MAX, LAT_MIN)
        if lat_first > lat_last
        else slice(LAT_MIN, LAT_MAX)
    )

    return ds.sel(
        longitude=slice(LON_MIN, LON_MAX),
        latitude=lat_slice,
    )


def require_variables(ds: xr.Dataset, variables: list[str], source: str) -> None:
    """Fail before downloading data if a required variable is absent."""
    missing = [name for name in variables if name not in ds.data_vars]
    if missing:
        raise KeyError(
            f"Variables missing from {source}: {missing}. Available variables: "
            f"{sorted(ds.data_vars)}"
        )


def expected_dates(year: int) -> pd.DatetimeIndex:
    return pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")


def require_complete_time_axis(
    time_values: np.ndarray,
    expected: pd.DatetimeIndex,
    label: str,
) -> None:
    """Check that timestamps are complete, ordered, and unique."""
    observed = pd.DatetimeIndex(time_values)
    if not observed.is_unique:
        raise ValueError(f"{label} contains duplicate timestamps.")
    if not observed.equals(expected):
        missing = expected.difference(observed)
        extra = observed.difference(expected)
        raise ValueError(
            f"{label} has an incomplete time axis. Missing={len(missing)}, "
            f"extra={len(extra)}. First missing values: {list(missing[:5])}"
        )


def build_year(
    hourly_source: xr.Dataset,
    daily_source: xr.Dataset,
    year: int,
) -> xr.Dataset:
    """Create one complete daily dataset without materialising hourly data."""
    dates = expected_dates(year)
    hourly_index = pd.date_range(
        f"{year}-01-01 00:00:00",
        f"{year}-12-31 23:00:00",
        freq="h",
    )

    hourly = hourly_source[SOIL_VARIABLES].sel(
        time=slice(hourly_index[0], hourly_index[-1])
    )
    require_complete_time_axis(
        hourly.time.values,
        hourly_index,
        f"Hourly ERA5-Land data for {year}",
    )

    # Soil water is an instantaneous state. The daily value is therefore the
    # arithmetic mean of all 24 hourly UTC values. skipna=False prevents a day
    # with missing hourly values from being reported as a complete daily mean.
    soil_daily = hourly.resample(time="1D").mean(skipna=False)
    require_complete_time_axis(
        soil_daily.time.values,
        dates,
        f"Daily soil-water data for {year}",
    )

    accumulated_daily = daily_source[DAILY_VARIABLES].sel(
        time=slice(dates[0], dates[-1])
    )
    require_complete_time_axis(
        accumulated_daily.time.values,
        dates,
        f"Daily accumulated ERA5-Land data for {year}",
    )

    result = xr.merge(
        [accumulated_daily, soil_daily],
        join="exact",
        compat="override",
    )[ALL_VARIABLES]

    result = result.astype({name: "float32" for name in ALL_VARIABLES})

    for layer, depth in {
        1: "0-7 cm",
        2: "7-28 cm",
        3: "28-100 cm",
        4: "100-289 cm",
    }.items():
        name = f"swvl{layer}"
        result[name].attrs.update(
            {
                "long_name": f"daily mean volumetric soil water, {depth}",
                "units": "m3 m-3",
                "daily_statistic": "mean of 24 hourly UTC values",
            }
        )

    result["ssrd"].attrs.update(
        {
            "long_name": "daily surface solar radiation downwards",
            "units": result["ssrd"].attrs.get("units", "J m-2"),
            "daily_statistic": "UTC daily accumulation",
        }
    )
    result["e"].attrs.update(
        {
            "long_name": "daily evaporation",
            "units": result["e"].attrs.get("units", "m of water equivalent"),
            "daily_statistic": "UTC daily accumulation",
        }
    )
    result["pev"].attrs.update(
        {
            "long_name": "daily potential evaporation",
            "units": result["pev"].attrs.get("units", "m"),
            "daily_statistic": "UTC daily accumulation",
        }
    )

    result.attrs.update(
        {
            "title": "Daily ERA5-Land variables for winter-wheat phenology",
            "time_zone": "UTC",
            "soil_water_source": HOURLY_URL,
            "accumulated_variable_source": DAILY_URL,
            "spatial_subset": (
                f"lon {LON_MIN}-{LON_MAX}; lat {LAT_MIN}-{LAT_MAX}"
            ),
            "processing_note": (
                "swvl1-swvl4 are daily means of hourly states; ssrd, e, and "
                "pev are daily UTC accumulations from the daily product"
            ),
        }
    )
    return result


def verify_output(path: Path, year: int) -> None:
    """Verify variables, dates, and physically possible soil-water ranges."""
    with xr.open_dataset(path, engine="h5netcdf") as ds:
        require_variables(ds, ALL_VARIABLES, str(path))
        require_complete_time_axis(
            ds.time.values,
            expected_dates(year),
            f"Written output {path.name}",
        )

        for name in SOIL_VARIABLES:
            minimum = float(ds[name].min(skipna=True).values)
            maximum = float(ds[name].max(skipna=True).values)
            if not np.isfinite(minimum) or not np.isfinite(maximum):
                raise ValueError(f"{path.name}: {name} contains no finite values.")
            if minimum < 0 or maximum > 1:
                raise ValueError(
                    f"{path.name}: {name} is outside [0, 1] m3 m-3 "
                    f"(min={minimum:.4f}, max={maximum:.4f})."
                )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    xr.set_options(keep_attrs=True)

    hourly_source = normalize_time_coordinate(open_zarr_dataset(HOURLY_URL))
    daily_source = normalize_time_coordinate(open_zarr_dataset(DAILY_URL))

    try:
        require_variables(hourly_source, SOIL_VARIABLES, "hourly source")
        require_variables(daily_source, DAILY_VARIABLES, "daily source")

        hourly_source = subset_region(hourly_source)
        daily_source = subset_region(daily_source)

        encoding = {
            name: {"zlib": True, "complevel": 4, "dtype": "float32"}
            for name in ALL_VARIABLES
        }

        for year in range(START_YEAR, END_YEAR + 1):
            output_file = OUTPUT_DIR / OUTPUT_TEMPLATE.format(year=year)

            if output_file.exists() and not OVERWRITE:
                verify_output(output_file, year)
                print(f"Verified existing file: {output_file}")
                continue
  
            print(f"Processing {year}...")
            annual = build_year(hourly_source, daily_source, year)
            temporary_file = output_file.with_suffix(".nc.tmp")
            annual.to_netcdf(
                temporary_file,
                engine="h5netcdf",
                encoding=encoding,
            )
            annual.close()

            verify_output(temporary_file, year)
            temporary_file.replace(output_file)
            print(f"Saved and verified: {output_file}")
    finally:
        hourly_source.close()
        daily_source.close()


if __name__ == "__main__":
    main()
