"""Calendar and simulation-identity checks; run with unittest discovery."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import xarray as xr
import xr_ccew as tw

from src import cmip
from src.wave_events import FIT_OPTIONS
from src.wave_fitting import EARTH_RADIUS_M, GRAVITY_M_S2, fit_event_sinusoid, resting_wave_depth


class CalendarFittingTests(unittest.TestCase):
    def test_known_wave_depth_in_each_calendar(self):
        """Identical sampled waves must recover 25 m, even around February 30."""
        c0 = np.sqrt(25 * GRAVITY_M_S2)
        beta = 2 * 7.292e-5 / EARTH_RADIUS_M
        for calendar in ("proleptic_gregorian", "standard", "noleap", "360_day"):
            times = xr.date_range("1999-11-01", periods=241, freq="D", calendar=calendar,
                                  use_cftime=calendar != "proleptic_gregorian")
            if calendar == "360_day":
                self.assertEqual((times[119].month, times[119].day), (2, 30))
            lon = np.arange(0.5, 360, 1)
            days = np.arange(len(times)) - 119
            for wave, K in (("Kelvin", 6), ("MRG", -5), ("n0EIG", 5), ("n1ERW", -6), ("n1WIG", -3)):
                with self.subTest(calendar=calendar, wave=wave):
                    k = K / EARTH_RADIUS_M
                    if wave == "Kelvin":
                        omega = c0*k
                    elif wave in ("MRG", "n0EIG"):
                        omega = (c0*k + np.sqrt((c0*k)**2 + 4*beta*c0)) / 2
                    else:
                        roots = np.roots([1, 0, -(c0*c0*k*k + 3*beta*c0), -beta*k*c0*c0])
                        positive = [root.real for root in roots if abs(root.imag) < 1e-12 and root.real > 0]
                        omega = max(positive) if wave == "n1WIG" else min(positive)
                    frequency = omega * 86400 / (2*np.pi)
                    values = 0.8 + 5*np.cos(K*np.deg2rad(lon)[None, :] - 2*np.pi*frequency*days[:, None] - 0.4)
                    field = xr.DataArray(values, dims=("time", "lon"), coords=dict(time=times, lon=lon))
                    fit, _ = fit_event_sinusoid(field, times[119], 180.5, **FIT_OPTIONS[wave])
                    self.assertTrue(fit["fit_success"])
                    np.testing.assert_allclose([fit["K_global"], fit["frequency_cpd"]], [K, frequency], rtol=1e-6, atol=1e-8)
                    depth = resting_wave_depth(wave, fit["k_rad_m"], fit["omega_rad_s"])
                    self.assertTrue(depth["dispersion_valid"])
                    self.assertAlmostEqual(depth["h_eff_m"], 25, places=4)

    def test_calendar_year_harmonics_are_removed(self):
        for calendar, year_days in (("noleap", 365), ("360_day", 360)):
            with self.subTest(calendar=calendar):
                times = xr.date_range("2000-01-01", periods=2*year_days, freq="D", calendar=calendar, use_cftime=True)
                phase = 2*np.pi*np.arange(len(times))/year_days
                field = xr.DataArray(240 + 10*np.sin(phase) + 3*np.cos(2*phase), dims="time", coords=dict(time=times))
                result = tw.remove_harmonics_of_seasonal_cycle(field, num_harmonics=3)
                np.testing.assert_allclose(result, 240, rtol=0, atol=1e-10)

    def test_fitting_rejects_a_missing_daily_sample(self):
        times = xr.date_range("2000-01-01", periods=61, freq="D", calendar="noleap", use_cftime=True)
        field = xr.DataArray(np.ones((61, 360)), dims=("time", "lon"), coords=dict(time=times, lon=np.arange(360)))
        field = field.drop_isel(time=29)
        with self.assertRaisesRegex(ValueError, "regularly sampled"):
            fit_event_sinusoid(field, times[30], 180, **FIT_OPTIONS["Kelvin"])

    def test_fit_uses_only_supplied_event_support(self):
        times = xr.date_range("2000-01-01", periods=361, freq="D")
        lon = np.arange(0.5, 360, 1)
        c0 = np.sqrt(25 * GRAVITY_M_S2)
        K = 6
        frequency = c0 * K / EARTH_RADIUS_M * 86400 / (2*np.pi)
        days = np.arange(len(times)) - 180
        values = 5*np.cos(K*np.deg2rad(lon)[None, :] - 2*np.pi*frequency*days[:, None])
        field = xr.DataArray(values, dims=("time", "lon"), coords=dict(time=times, lon=lon))
        relative_lon = (lon - 180.5 + 180) % 360 - 180
        support = np.where((abs(days[:, None]) <= 8) & (abs(relative_lon[None, :]) <= 60))
        fit, diagnostic = fit_event_sinusoid(field, times[180], 180.5,
            k_grid=FIT_OPTIONS["Kelvin"]["k_grid"],
            frequency_grid=FIT_OPTIONS["Kelvin"]["frequency_grid"],
            support_indices=support)
        self.assertEqual(diagnostic["observed"].dims, ("sample",))
        self.assertEqual(diagnostic["observed"].size, len(support[0]))
        self.assertAlmostEqual(fit["h_eff_m"], 25, places=3)


class CMIPLoaderTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.archive = self.root / "archive"
        self.archive.mkdir()
        (self.root / "data").mkdir()
        self.link = self.root / "data/cmip6_rlut_noaa1deg"
        self.link.symlink_to(self.archive, target_is_directory=True)
        self.addCleanup(patch.stopall)
        patch.object(cmip, "CMIP6_NOAA_GRID_DIR", self.link).start()
        patch.object(cmip, "PROJECT_ROOT", self.root).start()

    def write_fixture(self, calendar="proleptic_gregorian", missing_day=False, experiment="historical"):
        days = {"proleptic_gregorian": 366, "standard": 366, "noleap": 365, "360_day": 360}[calendar]
        times = xr.date_range("2000-01-01", periods=days, freq="D", calendar=calendar,
                              use_cftime=calendar != "proleptic_gregorian")
        lon = np.arange(-179.5, 180, 1)
        values = (240 + np.arange(days)[:, None, None]*0.01
                  + np.array([1, 2])[None, :, None] + np.cos(np.deg2rad(lon))[None, None, :])
        dataset = xr.Dataset(
            {"rlut": (("time", "lat", "lon"), values, dict(units="W m-2", standard_name="toa_outgoing_longwave_flux"))},
            coords=dict(time=times, lat=[0.500000000000004, -0.499999999999996], lon=lon),
            attrs=dict(source_id="TestModel", experiment_id=experiment, variant_label="r1i1p1f2",
                       grid_label="gr", source_grid_label="gn", analysis_grid="noaa_1deg",
                       regridding_fingerprint="test-fixture", table_id="day", frequency="day"),
        )
        if missing_day:
            dataset = dataset.drop_isel(time=30)
        path = self.archive / experiment / f"rlut.day.TestModel.{experiment}.r1i1p1f2.gn.v1.noaa_1deg.2000-2000.nc"
        path.parent.mkdir(exist_ok=True)
        dataset.to_netcdf(path)
        return dataset

    def test_coordinate_normalization_and_native_calendar(self):
        for calendar in ("proleptic_gregorian", "standard", "noleap", "360_day"):
            for experiment in ("historical", "amip"):
                with self.subTest(calendar=calendar, experiment=experiment):
                    original = self.write_fixture(calendar, experiment=experiment)
                    field = cmip.load_cmip_rlut("TestModel", experiment=experiment, year_range=(2000, 2000), member_id="r1i1p1f2")
                    self.assertEqual(field.dims, ("time", "lat", "lon"))
                    self.assertEqual(str(field.time.dt.calendar), calendar)
                    self.assertEqual(field.attrs["member_id"], "r1i1p1f2")
                    self.assertEqual(field.attrs["source_grid_label"], "gn")
                    self.assertEqual(field.attrs["analysis_grid"], "noaa_1deg")
                    self.assertEqual(field.attrs["experiment_id"], experiment)
                    self.assertEqual(cmip.list_cmip_models(experiment), ["TestModel"])
                    self.assertEqual(field.attrs["source_file"].split("/")[0], "data")
                    np.testing.assert_array_equal(field.lon, np.arange(0.5, 360, 1))
                    np.testing.assert_array_equal(field.lat, [-0.5, 0.5])
                    np.testing.assert_array_equal(field.sel(lon=0.5, lat=0.5), original.rlut.isel(lat=0).sel(lon=0.5))
                    projected = tw.symmetric_antisymmetric_component(field.isel(time=slice(0, 2), lon=slice(0, 2)), which="antisymmetric")
                    self.assertTrue(np.isfinite(projected.values).all())
                    with self.assertRaisesRegex(ValueError, "requested"):
                        cmip.load_cmip_rlut("TestModel", experiment=experiment, year_range=(2000, 2000), member_id="r2i1p1f2")
                    with self.assertRaisesRegex(ValueError, "every requested calendar year"):
                        cmip.load_cmip_rlut("TestModel", experiment=experiment, year_range=(1999, 2000))
        self.assertEqual(len(cmip.list_cmip_simulations()), 2)

    def test_missing_day_fails(self):
        self.write_fixture("360_day", missing_day=True)
        with self.assertRaisesRegex(ValueError, "gaps, duplicates"):
            cmip.load_cmip_rlut("TestModel", year_range=(2000, 2000))

    def test_ambiguous_model_fails(self):
        self.write_fixture()
        (self.archive / "historical/rlut.day.TestModel.historical.r2i1p1f2.gn.v1.noaa_1deg.2000-2000.nc").touch()
        with self.assertRaisesRegex(ValueError, "Multiple files"):
            cmip.list_cmip_models()


if __name__ == "__main__":
    unittest.main()
