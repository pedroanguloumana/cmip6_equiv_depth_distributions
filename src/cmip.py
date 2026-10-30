"""Load daily historical/AMIP rlut on the exact NOAA OLR 1-degree grid.

Configure data/cmip6_rlut_noaa1deg on each machine. The archive has one file and
one ensemble member per model; ambiguous matches fail rather than mixing members.
"""

import re

import numpy as np
import xarray as xr

from .paths import CMIP6_NOAA_GRID_DIR, PROJECT_ROOT


_FILENAME = re.compile(
    r"rlut\.day\.(?P<model>[^.]+)\.(?P<experiment>historical|amip)\."
    r"(?P<member>r\d+i\d+p\d+f\d+)\.(?P<grid>[^.]+)\.(?P<version>v\d+)\."
    r"noaa_1deg\.(?P<start_year>\d{4})-(?P<end_year>\d{4})\.nc"
)


def list_cmip_simulations(experiments=("historical", "amip")):
    """Return filename identities/paths, one per model/experiment, without field I/O.

    Each experiment is selected independently. Multiple members, source grids,
    versions or periods for one model/experiment are ambiguous and rejected.
    """
    if not CMIP6_NOAA_GRID_DIR.is_symlink() or not CMIP6_NOAA_GRID_DIR.is_dir():
        raise FileNotFoundError("Configure data/cmip6_rlut_noaa1deg as the NOAA-grid output archive symlink.")
    simulations = []
    seen = set()
    for experiment in experiments:
        if experiment not in ("historical", "amip"):
            raise ValueError(f"Unknown experiment: {experiment!r}")
        files = sorted((CMIP6_NOAA_GRID_DIR / experiment).glob("*.nc"))
        if not files:
            raise FileNotFoundError(f"No completed NOAA-grid rlut files for {experiment}.")
        for path in files:
            match = _FILENAME.fullmatch(path.name)
            if match is None or match["experiment"] != experiment:
                raise ValueError(f"Unexpected NOAA-grid CMIP rlut filename: {path}")
            key = (match["model"], experiment)
            if key in seen:
                raise ValueError(f"Multiple files found for {key}; select one member/grid/version/period explicitly.")
            seen.add(key)
            simulations.append(dict(path=path, source_id=match["model"], experiment_id=experiment,
                                    member_id=match["member"], source_grid_label=match["grid"],
                                    version=match["version"]))
    return simulations


def list_cmip_models(experiment="historical"):
    """Return sorted model IDs for one explicitly selected experiment."""
    return [item["source_id"] for item in list_cmip_simulations((experiment,))]


def inspect_cmip_rlut(model, *, experiment="historical"):
    """Read one model's source identity, calendar and coverage from its NetCDF header.

    Returns a dictionary including a Path, source_id, member_id, source_grid_label,
    analysis_grid, calendar, start and end. source_grid_label describes the
    original CMIP grid; grid_label=gr describes the new NOAA-grid product.
    """
    selections = [item for item in list_cmip_simulations((experiment,)) if item["source_id"] == model]
    if not selections:
        raise FileNotFoundError(f"No unique NOAA-grid rlut file for {model}/{experiment} in {CMIP6_NOAA_GRID_DIR}")
    selection = selections[0]
    path = selection["path"]
    # Keep CF dates even for modern Gregorian records; NumPy decoding would
    # relabel a source 'standard' calendar as proleptic_gregorian.
    with xr.open_dataset(path, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
        required = ("source_id", "experiment_id", "variant_label", "grid_label", "source_grid_label",
                    "analysis_grid", "table_id", "frequency", "regridding_fingerprint")
        missing = [key for key in required if not dataset.attrs.get(key)]
        if missing:
            raise ValueError(f"Missing CMIP identity metadata {missing} in {path.name}")
        if (dataset.attrs["source_id"] != model or dataset.attrs["experiment_id"] != experiment
                or dataset.attrs["variant_label"] != selection["member_id"]
                or dataset.attrs["source_grid_label"] != selection["source_grid_label"]
                or dataset.attrs["analysis_grid"] != "noaa_1deg" or dataset.attrs["grid_label"] != "gr"
                or dataset.attrs["table_id"] != "day" or dataset.attrs["frequency"] != "day"):
            raise ValueError(f"Filename and NOAA-grid historical/AMIP source metadata disagree in {path.name}")
        if dataset.attrs.get("member_id", dataset.attrs["variant_label"]) != dataset.attrs["variant_label"]:
            raise ValueError(f"Conflicting member_id and variant_label in {path.name}")
        field = dataset["rlut"]
        if (set(field.dims) != {"time", "lat", "lon"}
                or field.attrs.get("standard_name") != "toa_outgoing_longwave_flux"
                or field.attrs.get("units") != "W m-2"
                or field.attrs.get("positive", "up") != "up"):
            raise ValueError(f"Expected daily upward TOA longwave flux in W m-2 on (time, lat, lon): {path.name}")
        if not field.sizes["time"]:
            raise ValueError(f"Empty time coordinate in {path.name}")
        return dict(
            path=path, source_id=model, experiment_id=experiment,
            member_id=str(dataset.attrs["variant_label"]),
            source_grid_label=str(dataset.attrs["source_grid_label"]),
            analysis_grid="noaa_1deg", calendar=str(field.time.dt.calendar),
            version=selection["version"],
            regridding_fingerprint=str(dataset.attrs["regridding_fingerprint"]),
            start=str(field.time.values[0]), end=str(field.time.values[-1]),
        )


def load_cmip_rlut(model, *, experiment="historical", year_range=(1990, 2010), lat_range=(-20, 20), member_id=None):
    """Return one model's loaded rlut DataArray, in W m-2 on (time, lat, lon).

    Select inclusive complete calendar years and a latitude band in degrees north.
    Sort latitude and wrap/sort longitude to [0, 360), without interpolating space,
    changing flux sign, or converting calendars. Remove latitude-coordinate
    roundoff at 10 decimal places so mirrored endpoint latitudes match exactly.
    Daily spacing and complete years
    are required. Native cftime coordinates, units and source identity are retained.
    member_id optionally asserts the particular member expected by the caller.
    """
    info = inspect_cmip_rlut(model, experiment=experiment)
    if member_id is not None and member_id != info["member_id"]:
        raise ValueError(f"{model}: requested {member_id}, but the unique file contains {info['member_id']}")
    start_year, end_year = year_range
    if not (isinstance(start_year, int) and isinstance(end_year, int) and start_year <= end_year):
        raise ValueError("year_range must be an increasing pair of integer years.")
    if len(lat_range) != 2 or lat_range[0] >= lat_range[1]:
        raise ValueError("lat_range must be an increasing pair of latitude bounds.")

    with xr.open_dataset(info["path"], decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as dataset:
        # Year-only indexing works for both Gregorian and native CF calendars,
        # including a 360-day year whose December ends on day 30.
        field = dataset["rlut"].sel(time=slice(str(start_year), str(end_year))).sortby("lat")
        # Remapping produced coordinates such as -19.499999999999996 and
        # +19.500000000000004. Their mirrors otherwise fall microscopically
        # outside the interpolation domain and create NaNs at one band edge.
        # This changes coordinate roundoff only, not the data or grid resolution.
        lat_attrs = dict(field.lat.attrs)
        field = field.assign_coords(lat=np.round(field.lat.values, 10))
        field.lat.attrs = lat_attrs
        field = field.sel(lat=slice(*lat_range))
        if not field.sizes["time"] or not field.sizes["lat"]:
            raise ValueError(f"{model}: requested years/latitude band contain no data.")
        lon_attrs = dict(field.lon.attrs)
        field = field.assign_coords(lon=field.lon % 360).sortby("lon")
        field.lon.attrs = lon_attrs
        field = field.transpose("time", "lat", "lon")

        first, last = field.time.isel(time=0), field.time.isel(time=-1)
        if (int(first.dt.year) != start_year or int(first.dt.month) != 1 or int(first.dt.day) != 1
                or int(last.dt.year) != end_year or int(last.dt.month) != 12
                or int(last.dt.day) != int(last.dt.days_in_month)):
            raise ValueError(f"{model}: data do not cover every requested calendar year {start_year}--{end_year}.")
        times = field.time.values
        if np.issubdtype(times.dtype, np.datetime64):
            intervals = np.diff(times) / np.timedelta64(1, "D")
        else:
            intervals = np.array([(b-a).total_seconds()/86400 for a, b in zip(times[:-1], times[1:])])
        if not np.allclose(intervals, 1, rtol=0, atol=1e-9):
            raise ValueError(f"{model}: time contains gaps, duplicates or non-daily samples.")
        if (field.sizes["lon"] != 360 or not np.allclose(np.diff(field.lon), 1)
                or not np.allclose(np.diff(field.lat), 1)):
            raise ValueError(f"{model}: expected a regular global 1-degree longitude/latitude grid.")
        field = field.load()
    if not np.isfinite(field.values).all():
        raise ValueError(f"{model}: selected rlut contains missing or non-finite data.")

    # Keep rlut's physical metadata, and add the simulation identity explicitly.
    field.attrs = dict(field.attrs, **{key: info[key] for key in (
        "source_id", "experiment_id", "member_id", "source_grid_label", "analysis_grid", "calendar",
    )})
    field.attrs["source_file"] = str(info["path"].relative_to(PROJECT_ROOT))
    return field
