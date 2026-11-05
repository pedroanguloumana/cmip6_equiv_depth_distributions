#!/usr/bin/env python3
"""Calculate five-wave OLR variance/base-point/regression maps on Levante.

From the installed project root:
    python scripts/codex_olr_regression/run.py --task-index 0
Task 0 is NOAA; 1 CNRM-CM6-1 historical; 2 HadGEM3-GC31-LL historical.
For an independent small pilot add --start-year 1990 --end-year 1992.
Outputs: data/olr_wave_regression/<years>/<source>.<experiment>.<member>.nc
and an audit JSON, with exact source identity and processing parameters.
"""

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr
import xr_ccew

from src.cmip import inspect_cmip_rlut, load_cmip_rlut
from src.paths import DATA_DIR, NOAA_OLR_FILE, PROJECT_ROOT
from src.wave_profiles import BOX_PROFILES
from src.wave_regression import calculate_olr_wave_maps


RECORDS = (
    ("NOAA_OLR", "observations", "not_applicable"),
    ("CNRM-CM6-1", "historical", "r1i1p1f2"),
    ("HadGEM3-GC31-LL", "historical", "r1i1p1f3"),
)
WAVES = ("Kelvin", "n1ERW", "MRG", "n0EIG", "n1WIG")
LATITUDE_BOUND = 20


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--task-index", type=int, choices=range(3), required=True)
    parser.add_argument("--start-year", type=int, default=1990)
    parser.add_argument("--end-year", type=int, default=2010)
    args = parser.parse_args()
    if args.end_year < args.start_year:
        parser.error("--end-year must not precede --start-year")
    if PROJECT_ROOT != Path(__file__).resolve().parents[2]:
        raise RuntimeError("Install this checkout with pip -e . before running its worker.")
    source, experiment, member = RECORDS[args.task_index]
    output_dir = DATA_DIR / "olr_wave_regression" / f"{args.start_year}-{args.end_year}"
    output = output_dir / f"{source}.{experiment}.{member}.nc"
    if output.exists() or output.with_suffix(".audit.json").exists():
        raise FileExistsError(f"Preserving existing result: {output}")

    if source == "NOAA_OLR":
        path = NOAA_OLR_FILE
        if not path.parent.is_symlink() or not path.is_file():
            raise FileNotFoundError("Configure the data/noaa_olr archive symlink.")
        with xr.open_dataset(path) as dataset:
            field = dataset["olr"].sel(time=slice(str(args.start_year), str(args.end_year))).sortby("lat")
            field = field.sel(lat=slice(-LATITUDE_BOUND, LATITUDE_BOUND)).transpose("time", "lat", "lon").load()
        if (str(field.time.values[0])[:10] != f"{args.start_year}-01-01"
                or str(field.time.values[-1])[:10] != f"{args.end_year}-12-31"):
            raise ValueError("NOAA record does not cover the requested complete years.")
        identity = dict(source_id="NOAA_CDR", source_grid_label="NOAA_1deg", analysis_grid="noaa_1deg")
    else:
        info = inspect_cmip_rlut(source, experiment=experiment)
        path = info["path"]
        field = load_cmip_rlut(source, experiment=experiment, member_id=member,
                              year_range=(args.start_year, args.end_year),
                              lat_range=(-LATITUDE_BOUND, LATITUDE_BOUND))
        identity = {key: info[key] for key in ("source_id", "source_grid_label", "analysis_grid", "version",
                                              "regridding_fingerprint")}
    if not np.allclose(np.diff(field.lat), 1) or not np.allclose(np.diff(field.lon), 1):
        raise ValueError("Expected NOAA's regular 1-degree grid.")
    print(f"{source}/{experiment}/{member}: {dict(field.sizes)}, calendar={field.time.dt.calendar}", flush=True)
    result = calculate_olr_wave_maps(field, {wave: BOX_PROFILES[wave] for wave in WAVES},
                                     seasonal_harmonics=3, taper_fraction=0.1)
    stat = path.stat()
    code_paths = [Path(__file__), PROJECT_ROOT / "src/wave_regression.py", PROJECT_ROOT / "src/wave_profiles.py"]
    code_hashes = {str(p.relative_to(PROJECT_ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in code_paths}
    library_dir = Path(xr_ccew.__file__).parent
    library_hashes = {name: hashlib.sha256((library_dir/name).read_bytes()).hexdigest()
                      for name in ("spectral.py", "filtering.py")}
    result.attrs.update(identity, dataset=source, experiment_id=experiment, member_id=member,
                        source_file=str(path.relative_to(PROJECT_ROOT)), source_size_bytes=stat.st_size,
                        source_mtime_ns=stat.st_mtime_ns, source_units=field.attrs["units"],
                        source_attributes_json=json.dumps(field.attrs, sort_keys=True, default=str),
                        equivalent_depth_restriction="none; existing rectangular BOX_PROFILES",
                        code_sha256_json=json.dumps(code_hashes, sort_keys=True),
                        xr_ccew_sha256_json=json.dumps(library_hashes, sort_keys=True),
                        created_utc=datetime.now(timezone.utc).isoformat())
    output_dir.mkdir(parents=True, exist_ok=True)
    encoding = {name: dict(zlib=True, complevel=4) for name in result.data_vars}
    result.to_netcdf(output, encoding=encoding)
    with xr.open_dataset(output) as saved:
        xr.testing.assert_allclose(result.drop_vars(["filtered_reference", "time"]),
                                   saved.drop_vars(["filtered_reference", "time"]))
    audit = dict(output=str(output.relative_to(PROJECT_ROOT)),
                 sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                 attributes=result.attrs,
                 bases=[dict(wave=wave, latitude=float(result.base_latitude.sel(wave=wave)),
                             longitude=float(result.base_longitude.sel(wave=wave)),
                             variance_W2_m4=float(result.maximum_variance.sel(wave=wave))) for wave in WAVES])
    output.with_suffix(".audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True)+"\n")
    print(f"Saved {output.relative_to(PROJECT_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
