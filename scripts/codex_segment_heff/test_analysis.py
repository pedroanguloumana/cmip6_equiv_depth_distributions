"""Scientific checks for the separate segment/depth-band experiment.

Run from the project root:
    python -m unittest discover -s scripts/codex_segment_heff -p 'test_*.py' -v
"""

import unittest

import numpy as np
import xarray as xr
from xr_ccew.matsuno import kelvin_frequency, profile_curve_frequency

from scripts.codex_segment_heff.run import (
    BOX_PROFILES, DEPTH_EDGES_M, depth_band_weights, integrate_depth_bands,
    segment_spectra,
)


class SegmentDepthTests(unittest.TestCase):
    def field(self, nt=240, calendar="proleptic_gregorian"):
        lon = np.arange(36) * 10.0 + 0.5
        t = np.arange(nt)
        carrier = 10 * np.cos(6 * np.deg2rad(lon)[None, :] - 2 * np.pi * 0.1875 * t[:, None])
        return xr.DataArray(
            np.broadcast_to(carrier[:, None, :], (nt, 4, 36)).copy(),
            dims=("time", "lat", "lon"), attrs=dict(units="W m-2"),
            coords=dict(time=xr.date_range("1990-01-01", periods=nt, freq="D", calendar=calendar,
                                          use_cftime=calendar != "proleptic_gregorian"),
                        lat=[-1.5, -0.5, 0.5, 1.5], lon=lon))

    def test_exact_segment_length_eastward_sign_and_variance(self):
        result = segment_spectra(self.field(nt=132))
        self.assertEqual(result.sizes["segment"], 2)  # includes final exact window
        self.assertEqual(result.sizes["frequency"], 96)
        self.assertAlmostEqual(float(np.diff(result.frequency)[0]), 1 / 96)
        peak = result.symmetric.isel(segment=0).stack(cell=("frequency", "zonal_wavenumber")).idxmax()
        self.assertEqual(peak.item(), (-0.1875, -6))  # first of conjugate pair
        self.assertGreater(float(result.symmetric.sel(frequency=0.1875, zonal_wavenumber=6).mean()), 20)
        np.testing.assert_allclose(result.symmetric.sum(("frequency", "zonal_wavenumber")),
                                   result.symmetric_variance, rtol=1e-12)
        np.testing.assert_allclose(result.symmetric_variance, 50, rtol=0.03)
        self.assertEqual(float(result.antisymmetric.max()), 0)
        self.assertEqual(segment_spectra(self.field(nt=96)).sizes["segment"], 1)

    def test_native_calendar_and_non_daily_rejection(self):
        result = segment_spectra(self.field(nt=132, calendar="360_day"))
        self.assertEqual(result.segment_start.dt.calendar, "360_day")
        self.assertEqual((result.segment_end.values[0] - result.segment_start.values[0]).days, 95)
        with self.assertRaisesRegex(ValueError, "daily"):
            segment_spectra(self.field().isel(time=slice(None, None, 2)))

    def test_fractional_bands_partition_and_exclude_nyquist(self):
        f, k = np.arange(1, 49) / 96, np.arange(-15, 16)
        for name, profile in BOX_PROFILES.items():
            weights, scales = depth_band_weights(f, k, profile, DEPTH_EDGES_M)
            self.assertTrue(np.isfinite(weights).all(), name)
            self.assertTrue(np.all(weights >= 0), name)
            self.assertTrue(np.all(weights.sum(axis=0) <= 1 + 1e-12), name)
            self.assertEqual(float(weights[:, -1].sum()), 0, name)
            self.assertEqual(float(weights[:, :, k == 0].sum()), 0, name)
            merged, _ = depth_band_weights(f, k, profile, [10, 80])
            np.testing.assert_allclose(weights.sum(axis=0), merged[0], atol=1e-12)
        # At K=5 these frequencies are all inside the Kelvin box: recover the
        # analytic area between the 10 m and 80 m Kelvin curves.
        weights, kelvin_scale = depth_band_weights(f, k, BOX_PROFILES["Kelvin"], DEPTH_EDGES_M)
        expected = (kelvin_frequency(5, 80) - kelvin_frequency(5, 10)) * 96
        self.assertAlmostEqual(float(weights[:, :, k == 5].sum()), expected)
        _, er_scale = depth_band_weights(f, k, BOX_PROFILES["n1ERW"], [37.5, 42.5])
        self.assertGreater(float(er_scale[0, k == -5][0]), 15)
        self.assertLess(float(kelvin_scale[5, k == 5][0]), 5)

    def test_signed_linearity_positive_clipping_and_coarse_closure(self):
        f = np.sort(-np.fft.fftfreq(96))
        k = np.arange(-16, 16)
        # Segment deviations cancel exactly. Positive clipping still yields
        # excess, despite a perfectly flat mean spectrum with no wave peak.
        power = np.ones((2, len(f), len(k))) * np.array([0.5, 1.5])[:, None, None]
        source = xr.Dataset(
            dict(symmetric=(("segment", "frequency", "zonal_wavenumber"), power),
                 antisymmetric=(("segment", "frequency", "zonal_wavenumber"), power)),
            coords=dict(segment=[0, 1], frequency=f, zonal_wavenumber=k,
                        segment_start=("segment", np.array(["1990-01-01", "1990-02-06"], dtype="datetime64[D]")),
                        segment_end=("segment", np.array(["1990-04-06", "1990-05-12"], dtype="datetime64[D]"))))
        fine = integrate_depth_bands(source)
        coarse = integrate_depth_bands(source, depth_edges_m=[10, 30, 50, 80])
        np.testing.assert_allclose(fine.signed_excess.mean("segment"), 0, atol=1e-12)
        np.testing.assert_allclose(fine.mean_spectrum_positive_excess, 0, atol=1e-12)
        self.assertGreater(float(fine.positive_excess.sum()), 0)
        for name in ("observed_power", "signed_excess", "positive_excess", "positive_ratio_excess"):
            np.testing.assert_allclose(fine[name].sum("depth"), coarse[name].sum("depth"), atol=1e-12)

    def test_identical_depth_carriers_have_different_resolved_widths(self):
        # A shared true depth must not be assumed to produce identical output
        # PDFs: the same frequency window maps to a wider ER depth response.
        widths = {}
        for wave, k in (("Kelvin", 5), ("n1ERW", -5)):
            field = self.field(nt=731)
            profile = BOX_PROFILES[wave]
            frequency = float(profile_curve_frequency(
                profile.curve_name, k, 25,
                meridional_mode_number=profile.meridional_mode_number))
            carrier = 10 * np.cos(
                k * np.deg2rad(field.lon.values)[None, :]
                - 2 * np.pi * frequency * np.arange(731)[:, None])
            field.values[:] = carrier[:, None, :]
            result = integrate_depth_bands(segment_spectra(field)).sel(wave=wave)
            power = result.mean_spectrum_positive_excess.values
            cumulative = (power / power.sum()).cumsum()
            quartiles = result.depth.values[np.searchsorted(cumulative, [.25, .5, .75])]
            self.assertLess(abs(quartiles[1] - 25), 10)
            widths[wave] = quartiles[2] - quartiles[0]
        self.assertGreater(widths["n1ERW"], widths["Kelvin"])


if __name__ == "__main__":
    unittest.main()
