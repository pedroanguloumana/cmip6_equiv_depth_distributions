"""Code for handing NOAA OLR Observations.
"""

import xarray as xr


def load_noaa_olr(lat_range=None, lon_range=None, time_range=None):
    """Load NOAA daily OLR, optionally subset in space and time."""

    path = (
        "/work/bb1153/b382284/OBS/NOAA_OLR/OLR-Daily_v02r00_s19790101_e20251231.nc"
    )

    olr = xr.open_dataset(path)["olr"]

    if lat_range is not None:
        olr = olr.sel(lat=slice(*lat_range))

    if lon_range is not None:
        olr = olr.sel(lon=slice(*lon_range))

    if time_range is not None:
        olr = olr.sel(time=slice(*time_range))

    return olr