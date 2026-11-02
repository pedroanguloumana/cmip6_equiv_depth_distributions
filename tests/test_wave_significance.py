"""Check null covariance, commuting projections and held-out calibration."""

import unittest
from unittest.mock import patch

import numpy as np
import xarray as xr
import xr_ccew as tw

from src import wave_events as we
from src.wave_significance import (calibrate_event_peaks, distribution_diagnostics, draw_background,
                                   draw_marginal_background, estimate_background_power, project_olr,
                                   synthetic_wave_packet)


class BackgroundTests(unittest.TestCase):
    def field(self, nt=730, nx=72, calendar="proleptic_gregorian"):
        dates = xr.date_range("2000-01-01", periods=nt, freq="D", calendar=calendar,
                              use_cftime=calendar != "proleptic_gregorian")
        return xr.DataArray(np.random.default_rng(1).normal(size=(nt, nx)), dims=("time", "lon"),
            coords=dict(time=dates, lon=np.arange(nx)*360/nx+.5), attrs=dict(units="W m-2"))

    def test_envelope_seed_controls_event_footprint(self):
        field = self.field()
        day = np.arange(field.sizes["time"]) - 365
        lon = np.deg2rad(field.lon.values)
        field.data += 20*np.exp(-.5*(day[:, None]/8)**2)*np.cos(
            4*lon[None, :] - 2*np.pi*.1*day[:, None])
        sizes = []
        for cutoff in (1.0, 1.5, 2.0):
            events, sigma, support = we.detect_events(field, detection_sigma=cutoff,
                                                       return_support=True)
            strongest = max(events, key=lambda event: event["peak_amplitude_W_m2"])
            self.assertGreater(strongest["peak_amplitude_W_m2"], cutoff*sigma)
            sizes.append(len(support[strongest["event_id"]][0]))
        self.assertGreater(sizes[0], sizes[1])
        self.assertGreater(sizes[1], sizes[2])

    def test_gaussian_draw_preserves_target_covariance_and_coordinates(self):
        field = self.field(nt=512, nx=32, calendar="360_day")
        frequency, k = np.meshgrid(np.fft.fftfreq(512), np.fft.fftfreq(32), indexing="ij")
        power = 1 + 8*np.exp(-((frequency/.06)**2+(k/.15)**2))
        power[0, 0] = 0
        rng = np.random.default_rng(5)
        measured = np.zeros_like(power)
        variance, covariance = [], []
        for _ in range(100):
            draw = draw_background(power, field, rng)
            measured += abs(np.fft.fft2(draw.values))**2/draw.size/100
            variance.append(float((draw**2).mean()))
            covariance.append(np.mean(draw.values*np.roll(draw.values, 1, axis=0)))
        np.testing.assert_allclose(np.mean(variance), power.mean(), rtol=.02)
        expected_covariance = np.fft.ifft2(power).real[1, 0]
        np.testing.assert_allclose(np.mean(covariance), expected_covariance, rtol=.08)
        for mask in (abs(frequency)<.06, abs(frequency)>.2):
            np.testing.assert_allclose(measured[mask].mean(), power[mask].mean(), rtol=.03)
        self.assertEqual(draw.time.dt.calendar, "360_day")
        xr.testing.assert_equal(draw.lon, field.lon)
        self.assertLess(abs(float(draw.mean())), 1e-12)

    def test_smoothed_power_is_real_compatible_and_spreads_a_sharp_peak(self):
        field = self.field(nt=731, nx=73)
        t = np.arange(field.sizes["time"])[:, None]
        x = np.deg2rad(field.lon.values)[None, :]
        field.data = 10*np.cos(6*x-2*np.pi*.1*t)+field.values
        power = estimate_background_power(field)
        values = power.values
        nt, nx = values.shape
        np.testing.assert_allclose(values, values[np.ix_((-np.arange(nt)) % nt, (-np.arange(nx)) % nx)])
        self.assertTrue((values >= 0).all())
        # A sharp coherent peak is spread over neighboring spectral bins; the
        # null does not simply keep the observed wave-rich Fourier amplitudes.
        raw = abs(np.fft.fft2(field.values))**2/field.size
        self.assertLess(values.max(), raw.max()/20)
        self.assertGreater(power.attrs["expected_variance_W2_m4"], 10)
        with self.assertRaisesRegex(ValueError, "Hermitian"):
            broken = values.copy()
            broken[1, 2] *= 2
            draw_background(broken, field, np.random.default_rng(0))

    def test_projection_before_processing_matches_catalog_for_all_waves(self):
        for calendar in ("proleptic_gregorian", "360_day"):
            base = self.field(calendar=calendar)
            values = np.random.default_rng(2).normal(size=(730, 4, 72))
            olr = xr.DataArray(values+240, dims=("time", "lat", "lon"),
                coords=dict(time=base.time, lat=[-1.5, -.5, .5, 1.5], lon=base.lon), attrs=dict(units="W m-2"))
            projected = project_olr(olr)
            for wave, profile in we.WAVE_PROFILES.items():
                with self.subTest(calendar=calendar, wave=wave):
                    expected = we.filter_olr_for_wave(olr, wave)
                    field = projected[profile.symmetry]
                    anomaly = tw.apply_window(tw.remove_mean_and_linear_trend(
                        tw.remove_harmonics_of_seasonal_cycle(field, num_harmonics=we.SEASONAL_HARMONICS)),
                        window="tukey", pct=we.TAPER_FRACTION)
                    actual = tw.filter_field(anomaly, profile, apply_profile_symmetry=False)
                    np.testing.assert_allclose(actual.values, expected.values, rtol=1e-8, atol=1e-11)

    def test_validation_peaks_do_not_set_the_threshold(self):
        field = self.field()

        def events(peaks):
            return [dict(event_id=index, peak_time=field.time.values[100], peak_lon_deg=20,
                         duration_days=10, span_deg=30, peak_amplitude_W_m2=value)
                    for index, value in enumerate(peaks)]

        with patch.object(we, "detect_events", side_effect=[
            (events([4, 6]), 1), (events([2, 3]), 1), (events([4, 5]), 1),
            (events([1, 6]), 1), (events([1000]), 1),
        ]):
            result = calibrate_event_peaks(field, "Kelvin", n_realizations=4, quantile=.75)
        report = result["report"]
        self.assertEqual(report["candidate_peak_threshold_W_m2"], 5)
        self.assertEqual(report["heldout_false_candidates_per_record"], [1, 1])
        self.assertEqual(report["record_max_p_resolution"], 1/3)
        self.assertEqual(result["observed_events"].exceeds_null_candidate_peak_threshold.tolist(), [False, True])
        self.assertEqual(result["observed_events"].record_max_p_mc.iloc[1], 1/3)

    def test_marginal_surrogate_preserves_distribution_and_reports_spectral_error(self):
        field = self.field(nt=512, nx=32)
        field.data = np.exp(.3*field.values)
        field.data -= field.values.mean()
        power = abs(np.fft.fft2(field.values))**2/field.size
        power[0, 0] = 0
        surrogate, report = draw_marginal_background(power,field,np.random.default_rng(6))
        np.testing.assert_array_equal(np.sort(surrogate.values.ravel()),np.sort(field.values.ravel()))
        self.assertTrue(report["converged"],report)
        amplitude = np.sqrt(power*field.size)
        actual_error = np.linalg.norm(abs(np.fft.fft2(surrogate.values-surrogate.values.mean()))-amplitude)/np.linalg.norm(amplitude)
        self.assertAlmostEqual(actual_error,report["spectral_amplitude_relative_error"],places=12)
        _, limited = draw_marginal_background(power,field,np.random.default_rng(6),max_iterations=1,spectral_tolerance=1e-12)
        self.assertFalse(limited["converged"])

    def test_packet_truth_and_group_speed_for_all_branches_and_seam(self):
        from src.wave_fitting import EARTH_RADIUS_M, resting_wave_depth

        field = self.field()
        for wave,K in (("Kelvin",4),("n1ERW",-3),("MRG",-3),("n0EIG",3)):
            with self.subTest(wave=wave):
                packet, truth = synthetic_wave_packet(field,wave,K,25,center_index=365,center_lon=359.5,
                    time_sigma_periods=.5,lon_sigma_wavelengths=.2)
                equivalent, _ = synthetic_wave_packet(field,wave,K,25,center_index=365,center_lon=-.5,
                    time_sigma_periods=.5,lon_sigma_wavelengths=.2)
                np.testing.assert_allclose(packet,equivalent,rtol=0,atol=1e-13)
                depth = resting_wave_depth(wave,K/EARTH_RADIUS_M,truth["frequency_cpd"]*2*np.pi/86400)
                self.assertAlmostEqual(depth["h_eff_m"],25,places=3)
                if wave=="Kelvin":
                    self.assertAlmostEqual(truth["group_speed_deg_day"],truth["phase_speed_deg_day"],places=6)
                if wave=="MRG":
                    self.assertGreater(truth["group_speed_deg_day"],0)
                    self.assertLess(truth["phase_speed_deg_day"],0)

    def test_distribution_diagnostics_expose_variance_mixtures(self):
        field = self.field()
        field.data[:,field.lon.values>=160] *= 4
        rows = distribution_diagnostics(field)
        global_row = next(row for row in rows if row["region"]=="global" and row["season"]=="all")
        local = next(row for row in rows if row["region"]=="Pacific_160-280E" and row["season"]=="all")
        self.assertGreater(global_row["excess_kurtosis"],1)
        self.assertLess(abs(local["excess_kurtosis"]),.15)
        self.assertEqual(len(rows),25)

    def test_explicit_wig_box_supports_filtering_and_n1_injection(self):
        from src.wave_fitting import EARTH_RADIUS_M, GRAVITY_M_S2
        from src.wave_profiles import n1WIG_BOX
        from src.wave_significance import filter_projected_olr

        field = self.field()
        # Calibration passes the box explicitly to keep its recipe independent
        # of the event-fitting worker's profile selection.
        filtered = filter_projected_olr(field,"n1WIG",profile=n1WIG_BOX)
        self.assertTrue(np.isfinite(filtered.values).all())
        packet,truth = synthetic_wave_packet(field,"n1WIG",-3,8,center_index=365,
            center_lon=180,time_sigma_periods=1,lon_sigma_wavelengths=.5,profile=n1WIG_BOX)
        self.assertTrue(np.isfinite(packet.values).all())
        k = -3/EARTH_RADIUS_M
        omega = truth["frequency_cpd"]*2*np.pi/86400
        c = np.sqrt(GRAVITY_M_S2*8)
        beta = 2*7.292e-5/EARTH_RADIUS_M
        # The 3*beta*c coefficient distinguishes n=1 from the n=0 branch.
        residual = omega**3-(k*k*c*c+3*beta*c)*omega-beta*k*c*c
        self.assertLess(abs(residual/omega**3),1e-10)
        from src.wave_fitting import resting_wave_depth
        self.assertAlmostEqual(resting_wave_depth("n1WIG",k,omega)["h_eff_m"],8,places=4)
        self.assertLess(truth["phase_speed_deg_day"],0)
        self.assertTrue(truth["carrier_inside_rectangular_bounds"])


if __name__ == "__main__":
    unittest.main()
