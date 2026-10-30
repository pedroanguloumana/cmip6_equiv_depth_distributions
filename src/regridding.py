"""Conservative remapping of daily CMIP6 OLR onto the NOAA OLR grid."""

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import tempfile
import time

import netCDF4
import numpy as np
import xarray as xr

from .paths import DATA_DIR, NOAA_OLR_FILE, PROJECT_ROOT


def noaa_target_grid():
    """Return NOAA's global 1-degree centers and cell bounds, without OLR data.

    A dummy two-dimensional field lets CDO read the grid as a NetCDF grid file.
    Coordinates and bounds are copied, not reconstructed from nominal spacing.
    """
    if not NOAA_OLR_FILE.parent.is_symlink():
        raise FileNotFoundError("Configure the data/noaa_olr archive symlink first.")
    with xr.open_dataset(NOAA_OLR_FILE) as source:
        grid = source[["lat_bnds", "lon_bnds"]].load()
    grid.encoding.pop("unlimited_dims", None)
    if (grid.sizes["lat"] != 180 or grid.sizes["lon"] != 360
            or not np.allclose(np.diff(grid.lat), 1)
            or not np.allclose(np.diff(grid.lon), 1)):
        raise ValueError("Expected the global, regular 1-degree NOAA OLR grid.")
    grid["target"] = xr.DataArray(
        np.zeros((180, 360), dtype="float32"), dims=("lat", "lon"),
        attrs=dict(long_name="NOAA OLR target grid", units="1"),
    )
    return grid


def check_daily_coverage(times, start_year, end_year):
    """Require every day of inclusive native-calendar years, without duplicates.

    Return the calendar name. Gregorian and CF calendars are supported; a
    360-day year correctly ends on December 30. No temporal interpolation occurs.
    """
    coordinate = xr.DataArray(times, dims="time")
    if not coordinate.size:
        raise ValueError("The requested period has no time samples.")
    first, last = coordinate.isel(time=0), coordinate.isel(time=-1)
    if (int(first.dt.year) != start_year or int(first.dt.month) != 1 or int(first.dt.day) != 1
            or int(last.dt.year) != end_year or int(last.dt.month) != 12
            or int(last.dt.day) != int(last.dt.days_in_month)):
        raise ValueError(f"Incomplete requested years: expected {start_year}--{end_year}.")
    if np.issubdtype(coordinate.dtype, np.datetime64):
        steps = np.diff(coordinate.values) / np.timedelta64(1, "D")
    else:
        steps = np.array([(b-a).total_seconds()/86400 for a, b in zip(times[:-1], times[1:])])
    if not np.allclose(steps, 1, rtol=0, atol=1e-9):
        raise ValueError("Time contains gaps, overlaps, duplicates or non-daily samples.")
    return str(coordinate.dt.calendar)


def regrid_rlut(task, output_root, *, cdo="cdo"):
    """Process one explicit model/experiment/member/grid/version task.

    Inputs are a planner's JSON task and a configured output symlink under data/.
    Use first-order conservative remapping of the intensive flux, in W m-2.
    Keep the complete globe, daily timestamps and native calendar. Only rlut is
    remapped: time is not averaged, and no seasonal or wave filtering is applied.

    Each native-calendar year is staged as classic NetCDF. This bounds temporary
    storage and avoids the HDF5 read/write concurrency failures observed in the
    previous workflow. Regridding weights are computed once per simulation.
    """
    output_root = Path(output_root)
    if not output_root.is_symlink() or not output_root.is_dir():
        raise FileNotFoundError(f"Configure an output archive symlink at {output_root}.")
    if not output_root.absolute().is_relative_to(DATA_DIR.absolute()):
        raise ValueError("Output archives must be exposed through a data/ symlink.")
    start_year, end_year = task["year_range"]
    identity = task["simulation"]
    grid = noaa_target_grid()
    grid_hash = hashlib.sha256(b"".join(
        np.asarray(grid[name], dtype="float64").tobytes()
        for name in ("lat", "lon", "lat_bnds", "lon_bnds")
    )).hexdigest()
    if grid_hash != task["target_grid_sha256"]:
        raise ValueError("The NOAA target grid changed after this manifest was prepared.")

    paths, file_times, file_indices = [], [], []
    time_owners, duplicate_records = {}, []
    native_grid_hash = None
    decode_times = xr.coders.CFDatetimeCoder(use_cftime=True)
    # Check every selected file's identity and time coordinate before moving any
    # climate fields. A catalog-range match alone does not establish daily cadence.
    for entry in task["files"]:
        path = PROJECT_ROOT / entry["path"]
        if not path.is_relative_to(DATA_DIR) or ".." in Path(entry["path"]).parts:
            raise ValueError("Source paths must enter through the project's data/ directory.")
        if (path.stat().st_size, path.stat().st_mtime_ns) != (entry["size"], entry["mtime_ns"]):
            raise ValueError(f"Source file changed after catalog selection: {path}")
        with xr.open_dataset(path, decode_times=decode_times) as source:
            for key in ("source_id", "experiment_id"):
                if source.attrs.get(key) != identity[key]:
                    raise ValueError(f"Source {key} conflicts with the manifest: {path}")
            expected = {"variant_label": identity["member_id"], "grid_label": task["source_grid_label"],
                        "table_id": task["table_id"], "frequency": "day"}
            if any(source.attrs.get(key) != value for key, value in expected.items()):
                raise ValueError(f"Source member/grid/table/frequency conflicts with the manifest: {path}")
            flux = source["rlut"]
            if (flux.attrs.get("standard_name") != "toa_outgoing_longwave_flux"
                    or flux.attrs.get("units") != "W m-2"
                    or flux.attrs.get("positive", "up") != "up"
                    or "time: mean" not in flux.attrs.get("cell_methods", "")
                    or flux.ndim not in (2, 3) or "time" not in flux.dims):
                raise ValueError(f"Expected daily-mean TOA outgoing flux in W m-2: {path}")
            horizontal = [coordinate for coordinate in flux.coords.values()
                          if coordinate.attrs.get("standard_name") in ("latitude", "longitude")]
            if len(horizontal) != 2:
                raise ValueError(f"Expected explicit CF latitude/longitude coordinates: {path}")
            signature = hashlib.sha256(b"".join(
                str((coordinate.name, coordinate.dims, coordinate.shape)).encode()
                + coordinate.values.tobytes() for coordinate in sorted(horizontal, key=lambda c: c.name)
            )).hexdigest()
            if native_grid_hash is not None and signature != native_grid_hash:
                raise ValueError(f"Source grid changes between files; weights cannot be reused: {path}")
            native_grid_hash = signature
            years = source.time.dt.year.values
            indices = np.flatnonzero((years >= start_year) & (years <= end_year))
            times = source.time.values[indices]
            retained = []
            overlaps = {}
            for index, stamp in zip(indices, times):
                if stamp in time_owners:
                    previous_path, previous_index = time_owners[stamp]
                    overlaps.setdefault(previous_path, []).append((previous_index, index))
                else:
                    time_owners[stamp] = (path, index)
                    retained.append(index)
            # Some CESM2-WACCM-FV2 archive chunks overlap. Retain the first
            # occurrence only after verifying exact equality of decoded OLR in
            # bounded blocks. Conflicting records are never silently discarded.
            for previous_path, pairs in overlaps.items():
                with xr.open_dataset(previous_path, decode_times=decode_times) as previous:
                    for first in range(0, len(pairs), 31):
                        old, new = zip(*pairs[first:first+31])
                        left = previous.rlut.isel(time=list(old)).values
                        right = flux.isel(time=list(new)).values
                        if not np.array_equal(left, right, equal_nan=True):
                            raise ValueError(f"Conflicting OLR values in overlapping files: {previous_path}, {path}")
                duplicate_records.append(dict(retained_source=str(previous_path.relative_to(PROJECT_ROOT)),
                                              duplicate_source=entry["path"], verified_days=len(pairs)))
                print(f"Verified {len(pairs)} identical overlapping days: {path.name}", flush=True)
            retained = np.asarray(retained, dtype=int)
            times = source.time.values[retained]
        if len(retained):
            paths.append(path)
            file_times.append(times)
            file_indices.append(retained)
    all_times = np.asarray(sorted(time_owners))
    calendar = check_daily_coverage(np.asarray(all_times), start_year, end_year)
    if calendar != task["calendar"]:
        raise ValueError("The selected files' calendar differs from the catalog evidence.")

    filename = (f"rlut.day.{identity['source_id']}.{identity['experiment_id']}."
                f"{identity['member_id']}.{task['source_grid_label']}.{task['version']}."
                f"noaa_1deg.{start_year}-{end_year}.nc")
    output = output_root / identity["experiment_id"] / filename
    output.parent.mkdir(parents=True, exist_ok=True)
    fingerprint = hashlib.sha256(json.dumps(task, sort_keys=True).encode()).hexdigest()
    if output.exists():
        with xr.open_dataset(output, decode_times=decode_times) as existing:
            if existing.attrs.get("regridding_fingerprint") != fingerprint:
                raise FileExistsError(f"An output with different provenance already exists: {output}")
            check_daily_coverage(existing.time.values, start_year, end_year)
            np.testing.assert_array_equal(existing.lat, grid.lat)
            np.testing.assert_array_equal(existing.lon, grid.lon)
        print(f"Already complete: {output}", flush=True)
        return output

    started = time.monotonic()
    environment = dict(os.environ, OMP_NUM_THREADS="1", CDO_REMAP_NORM="fracarea")
    commands = []

    def run(arguments):
        command = [str(cdo), "-O", "-L", *map(str, arguments)]
        commands.append(command)
        subprocess.run(command, env=environment, check=True)

    with tempfile.TemporaryDirectory(prefix=filename + ".partial.", dir=output.parent) as temporary:
        temporary = Path(temporary)
        target = temporary / "target.nc"
        grid.to_netcdf(target)
        weights = temporary / "weights.nc"
        years = []
        for year in range(start_year, end_year + 1):
            native = temporary / "native.nc"
            inputs = []
            for path, times, indices in zip(paths, file_times, file_indices):
                selected = indices[xr.DataArray(times).dt.year.values == year]
                if not len(selected):
                    continue
                # CDO indices are one-based. Split retained indices into
                # contiguous runs so no verified duplicate is merged twice.
                runs = np.split(selected, np.where(np.diff(selected) != 1)[0]+1)
                for run_indices in runs:
                    inputs.extend((f"-seltimestep,{run_indices[0]+1}/{run_indices[-1]+1}", "-selname,rlut", str(path)))
            run(["-f", "nc2", "-b", "F32", "-mergetime", *inputs, native])
            if not weights.exists():
                run(["-f", "nc2", f"-gencon,{target}", "-seltimestep,1", native, weights])
            year_file = temporary / f"{year}.nc"
            run(["-f", "nc2", "-b", "F32", f"-remap,{target},{weights}", native, year_file])
            native.unlink()
            years.append(year_file)
            print(f"{identity['source_id']}/{identity['experiment_id']}: remapped {year}", flush=True)
        partial = temporary / "complete.nc"
        run(["-f", "nc4", "-z", "zip_4", "-b", "F32", "-mergetime", *years, partial])

        # Preserve the reference's exact stored coordinates and bounds. CDO can
        # regenerate regular-grid coordinates with tiny floating-point changes.
        with netCDF4.Dataset(partial, "r+") as destination:
            for name in ("lat", "lon", "lat_bnds", "lon_bnds"):
                # CDO may omit bounds for a regular grid. Add the reference
                # bounds explicitly rather than relying on generated corners.
                for dimension in grid[name].dims:
                    if dimension not in destination.dimensions:
                        destination.createDimension(dimension, grid.sizes[dimension])
                if name not in destination.variables:
                    destination.createVariable(name, "f8", grid[name].dims)
                destination[name][:] = grid[name].values
                destination[name].setncatts({key: value for key, value in grid[name].attrs.items()
                                           if key != "_FillValue"})
            destination.setncatts(dict(
                source_grid_label=task["source_grid_label"], grid_label="gr",
                analysis_grid="noaa_1deg", grid="NOAA CDR OLR global 1-degree grid",
                regridding_method="first-order conservative; CDO_REMAP_NORM=fracarea",
                regridding_fingerprint=fingerprint,
            ))
        with xr.open_dataset(partial, decode_times=decode_times) as result:
            check_daily_coverage(result.time.values, start_year, end_year)
            np.testing.assert_array_equal(result.time.values, np.asarray(all_times))
            for name in ("lat", "lon", "lat_bnds", "lon_bnds"):
                np.testing.assert_array_equal(result[name], grid[name])
            if result.rlut.attrs.get("units") != "W m-2":
                raise ValueError("Regridding unexpectedly changed the flux units.")
            # Read bounded blocks: an incomplete or missing-valued product must
            # never be renamed into place as a completed analysis input.
            for first in range(0, result.sizes["time"], 31):
                if not np.isfinite(result.rlut.isel(time=slice(first, first+31)).values).all():
                    raise ValueError("Regridded rlut contains missing/non-finite cells.")
        report = dict(task=task, fingerprint=fingerprint, commands=commands,
                      verified_duplicate_records=duplicate_records,
                      workflow_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                      completed_utc=datetime.now(timezone.utc).isoformat(),
                      elapsed_seconds=time.monotonic()-started,
                      peak_process_rss_gib=max(resource.getrusage(who).ru_maxrss
                                               for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN))
                                               / (1024**3 if sys.platform == "darwin" else 1024**2),
                      cdo_version=subprocess.check_output([str(cdo), "-V"], stderr=subprocess.STDOUT, text=True))
        sidecar = temporary / "complete.json"
        sidecar.write_text(json.dumps(report, indent=2) + "\n")
        sidecar.replace(output.with_suffix(".json"))
        partial.replace(output)
    print(f"Complete: {output}", flush=True)
    return output
