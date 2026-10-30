"""Exercise calendar coverage and conservative CDO remapping on synthetic OLR.

Set CDO_BIN to a working CDO executable for the integration check. No archive
data are used; all source/output paths enter through temporary data/ symlinks.
"""

import hashlib
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import xarray as xr

from src import regridding


class DailyCoverageTests(unittest.TestCase):
    def test_complete_native_calendar_years_and_bad_steps(self):
        for calendar, days in (("standard", 365), ("noleap", 365), ("360_day", 360)):
            with self.subTest(calendar=calendar):
                times = xr.date_range("1990-01-01", periods=days, freq="D", calendar=calendar, use_cftime=True)
                self.assertEqual(regridding.check_daily_coverage(times.values, 1990, 1990), calendar)
                with self.assertRaisesRegex(ValueError, "gaps, overlaps"):
                    regridding.check_daily_coverage(np.delete(times.values, 30), 1990, 1990)
                with self.assertRaisesRegex(ValueError, "Incomplete"):
                    regridding.check_daily_coverage(times.values[1:], 1990, 1990)
                with self.assertRaisesRegex(ValueError, "gaps, overlaps"):
                    regridding.check_daily_coverage(np.r_[times.values, times.values], 1990, 1990)


@unittest.skipUnless(shutil.which(os.environ.get("CDO_BIN", "cdo")), "CDO_BIN/cdo is required")
class ConservativeRegriddingTests(unittest.TestCase):
    def test_conservation_exact_grid_dates_and_restart_in_three_calendars(self):
        for calendar, days in (("standard", 730), ("noleap", 730), ("360_day", 720)):
            with self.subTest(calendar=calendar), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                data = root / "data"
                data.mkdir()
                for name in ("noaa_olr", "cmip6", "output"):
                    archive = root / (name + "_archive")
                    archive.mkdir()
                    (data / name).symlink_to(archive, target_is_directory=True)
                lon, lat = np.arange(0.5, 360, 1), np.arange(-89.5, 90, 1)
                reference = xr.Dataset(
                    {"lat_bnds": (("lat", "bnds"), np.column_stack([lat-.5, lat+.5])),
                     "lon_bnds": (("lon", "bnds"), np.column_stack([lon-.5, lon+.5]))},
                    coords=dict(lat=("lat", lat, dict(standard_name="latitude", units="degrees_north", bounds="lat_bnds")),
                                lon=("lon", lon, dict(standard_name="longitude", units="degrees_east", bounds="lon_bnds"))),
                )
                reference.to_netcdf(data / "noaa_olr/reference.nc")
                times = xr.date_range("1990-01-01", periods=days, freq="D", calendar=calendar, use_cftime=True)
                source_lat, source_lon = np.arange(-89, 90, 2), np.arange(1, 360, 2)
                values = (240 + np.arange(days, dtype="float32")[:, None, None]*.01
                          + 5*np.sin(np.deg2rad(source_lat))[None, :, None]
                          + 2*np.cos(np.deg2rad(source_lon))[None, None, :]).astype("float32")
                source = xr.Dataset(
                    {"rlut": (("time", "lat", "lon"), values,
                              dict(units="W m-2", standard_name="toa_outgoing_longwave_flux", cell_methods="time: mean"))},
                    coords=dict(time=times,
                                lat=("lat", source_lat, dict(standard_name="latitude", units="degrees_north")),
                                lon=("lon", source_lon, dict(standard_name="longitude", units="degrees_east"))),
                    attrs=dict(source_id="TestModel", experiment_id="amip", variant_label="r1i1p1f1",
                               grid_label="gn", table_id="day", frequency="day"),
                )
                paths = [data/"cmip6/first.nc", data/"cmip6/second.nc"]
                source.isel(time=slice(0, 400)).to_netcdf(paths[0])
                second = source.isel(time=slice(330, None)).copy(deep=True)
                second.to_netcdf(paths[1])
                signature = hashlib.sha256(b"".join(reference[name].values.astype("float64").tobytes()
                                                   for name in ("lat", "lon", "lat_bnds", "lon_bnds"))).hexdigest()
                task = dict(simulation=dict(source_id="TestModel", experiment_id="amip", member_id="r1i1p1f1"),
                            year_range=[1990, 1991], source_grid_label="gn", table_id="day", version="v1",
                            calendar=calendar, target_grid_sha256=signature,
                            files=[dict(path=str(path.relative_to(root)), size=path.stat().st_size,
                                        mtime_ns=path.stat().st_mtime_ns) for path in paths])
                with patch.multiple(regridding, PROJECT_ROOT=root, DATA_DIR=data, NOAA_OLR_FILE=data/"noaa_olr/reference.nc"):
                    output = regridding.regrid_rlut(task, data/"output", cdo=os.environ.get("CDO_BIN", "cdo"))
                    original_mtime = output.stat().st_mtime_ns
                    regridding.regrid_rlut(task, data/"output", cdo=os.environ.get("CDO_BIN", "cdo"))
                    self.assertEqual(output.stat().st_mtime_ns, original_mtime)
                    second.rlut.values[0, 0, 0] += 1
                    second.to_netcdf(paths[1])
                    task["files"][1].update(size=paths[1].stat().st_size, mtime_ns=paths[1].stat().st_mtime_ns)
                    with self.assertRaisesRegex(ValueError, "Conflicting OLR"):
                        regridding.regrid_rlut(task, data/"output", cdo=os.environ.get("CDO_BIN", "cdo"))
                    self.assertEqual(output.stat().st_mtime_ns, original_mtime)
                with xr.open_dataset(output, decode_times=xr.coders.CFDatetimeCoder(use_cftime=True)) as result:
                    self.assertEqual(str(result.time.dt.calendar), calendar)
                    np.testing.assert_array_equal(result.time.values, times.values)
                    for name in ("lat", "lon", "lat_bnds", "lon_bnds"):
                        np.testing.assert_array_equal(result[name], reference[name])
                    self.assertEqual(result.rlut.attrs["units"], "W m-2")
                    self.assertEqual(result.attrs["source_grid_label"], "gn")
                    self.assertEqual(result.attrs["grid_label"], "gr")
                    # Global conservative remapping must preserve each day's
                    # area mean. The odd-latitude and cosine-longitude terms
                    # integrate to zero on both complete global grids.
                    area = np.sin(np.deg2rad(reference.lat_bnds[:, 1])) - np.sin(np.deg2rad(reference.lat_bnds[:, 0]))
                    means = result.rlut.weighted(area).mean(("lat", "lon")).values
                    np.testing.assert_allclose(means, 240 + np.arange(days)*.01, rtol=0, atol=5e-5)


if __name__ == "__main__":
    unittest.main()
