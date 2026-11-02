"""Check WK99 spectrum smoothing, direction, normalization and Gaussian draws."""

import unittest
from unittest.mock import patch

import numpy as np
import xarray as xr
import xr_ccew as tw

from src.wave_significance import draw_background, estimate_wk99_background_power


class WK99BackgroundTests(unittest.TestCase):
    def components(self, nt=731, nx=48, calendar="proleptic_gregorian"):
        time = xr.date_range("2000-01-01",periods=nt,freq="D",calendar=calendar,
                             use_cftime=calendar!="proleptic_gregorian")
        lon = np.arange(nx)*360/nx+.5
        wave = 10*np.cos(6*np.deg2rad(lon)[None,:]-2*np.pi*.125*np.arange(nt)[:,None])
        noise = np.random.default_rng(8).normal(size=(nt,nx))
        field = xr.DataArray(wave+noise,dims=("time","lon"),coords=dict(time=time,lon=lon),attrs=dict(units="W m-2"))
        return xr.Dataset(dict(symmetric=field,antisymmetric=field.copy(data=noise)))

    def test_library_smoothing_direction_and_pooled_variance(self):
        components = self.components()
        with patch.object(tw,"background_spectrum",wraps=tw.background_spectrum) as smoother:
            result = estimate_wk99_background_power(components)
        self.assertEqual([(c.kwargs["passes_wavenumber"],c.kwargs["passes_frequency"]) for c in smoother.call_args_list],
                         [(10,0),(20,0),(40,0),(0,10)])
        power = result.background_power
        nt,nx = power.shape
        self.assertTrue(np.isfinite(power.values).all())
        self.assertTrue((power.values>=0).all())
        np.testing.assert_allclose(power.values,power.values[np.ix_((-np.arange(nt))%nt,(-np.arange(nx))%nx)])
        self.assertEqual(power.values[0,0],0)
        expected = float((result.segment_power_symmetric.sum()+result.segment_power_antisymmetric.sum())/2)
        self.assertAlmostEqual(float(power.mean()),expected,places=10)
        self.assertGreater(float(result.segment_power_symmetric.sum()),5*float(result.segment_power_antisymmetric.sum()))
        # A known eastward carrier must remain eastward after the library's WK
        # frequency labels are converted to NumPy FFT coordinates.
        ordered = power.sortby("fft_frequency").sortby("fft_wavenumber")
        east = float(ordered.sel(fft_frequency=-.125,fft_wavenumber=6,method="nearest"))
        west = float(ordered.sel(fft_frequency=-.125,fft_wavenumber=-6,method="nearest"))
        self.assertGreater(east,5*west)
        raw_peak = float(result.pooled_segment_power.max())
        self.assertLess(float(result.smoothed_positive_background.max()),raw_peak/5)
        self.assertEqual(result.attrs["segment_days"],96)
        self.assertEqual(result.attrs["overlap_days"],60)

    def test_gaussian_realizations_match_target_and_native_calendar(self):
        components = self.components(nt=730,nx=49,calendar="360_day")
        result = estimate_wk99_background_power(components)
        power,template = result.background_power,components.symmetric
        rng = np.random.default_rng(3)
        variance,covariance = [],[]
        for _ in range(80):
            draw = draw_background(power,template,rng)
            variance.append(float((draw**2).mean()))
            covariance.append(float(np.mean(draw.values*np.roll(draw.values,1,axis=0))))
        np.testing.assert_allclose(np.mean(variance),float(power.mean()),rtol=.04)
        np.testing.assert_allclose(np.mean(covariance),np.fft.ifft2(power.values).real[1,0],rtol=.08)
        self.assertLess(abs(float(draw.mean())),1e-12)
        self.assertEqual(draw.time.dt.calendar,"360_day")
        xr.testing.assert_equal(draw.lon,template.lon)

    def test_non_daily_and_misaligned_inputs_fail(self):
        components = self.components()
        with self.assertRaisesRegex(ValueError,"daily"):
            estimate_wk99_background_power(components.isel(time=slice(None,None,2)))
        components["antisymmetric"][0,0] = np.nan
        with self.assertRaisesRegex(ValueError,"finite"):
            estimate_wk99_background_power(components)


if __name__=="__main__":
    unittest.main()
