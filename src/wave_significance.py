"""Calibrate packet amplitudes against a stationary, correlated Gaussian null.

This is a diagnostic calibration, separate from the sinusoid fit quality flags.
The background is a smoothed spectrum of the same hemispheric projection used
by the catalogs, not independent Gaussian pixels or a phase scramble retaining
the original sharp wave peaks. Its smoothing choices require sensitivity tests.
"""

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
import xarray as xr
import xr_ccew as tw

from . import wave_events as we


def project_olr(olr):
    """Return symmetric/antisymmetric northern means in W m-2 on (time, lon).

    Symmetry projection, hemispheric averaging, seasonal/trend removal and the
    Fourier filter are linear, so this reduction can precede temporal processing.
    It preserves the covariance of the averaged signal without simulating latitude.
    """
    if set(olr.dims) != {"time", "lat", "lon"} or olr.attrs.get("units") != "W m-2":
        raise ValueError("Expected OLR in W m-2 on (time, lat, lon).")
    if not np.allclose(olr.lat.values, -olr.lat.values[::-1]):
        raise ValueError("Expected increasing, paired equatorial latitudes.")
    north = olr.sel(lat=slice(*we.AVERAGE_LAT_RANGE)).mean("lat")
    south = olr.sel(lat=slice(-we.AVERAGE_LAT_RANGE[1], -we.AVERAGE_LAT_RANGE[0])).mean("lat")
    components = xr.Dataset(dict(symmetric=(north + south)/2, antisymmetric=(north - south)/2))
    for field in components.values():
        field.attrs = dict(olr.attrs)
    return components.transpose("time", "lon").load()


def estimate_background_power(field, *, frequency_smoothing_cpd=1/96, wavenumber_smoothing=2.0):
    """Estimate a smooth 2-D null spectrum from a projected daily OLR record.

    Input: finite (time, lon), daily, regular global longitude, W m-2. Output:
    power in unshifted NumPy FFT ordering, with E[abs(FFT(q))**2]/N = power.
    Thus sum(power)/N is the expected variance of generated, untapered fields.
    Gaussian kernel widths are in cycles/day and global wavenumber, respectively.
    This is a WK-inspired smoothing model, not the original WK background recipe.
    """
    if field.dims != ("time", "lon") or field.attrs.get("units") != "W m-2":
        raise ValueError("Expected projected OLR in W m-2 on (time, lon).")
    if not np.isfinite(field.values).all() or field.sizes["time"] <= 2*we.REFERENCE_EXCLUDE_DAYS:
        raise ValueError("Expected a finite record longer than the sigma-reference exclusions.")
    longitude = field.lon.values
    if (len(longitude) < 2 or not np.allclose(np.diff(longitude), 360/len(longitude))):
        raise ValueError("Expected an increasing, regular global longitude grid.")
    times = field.time.values
    if np.issubdtype(times.dtype, np.datetime64):
        steps = np.diff(times)/np.timedelta64(1, "D")
    else:
        steps = [(b-a).total_seconds()/86400 for a, b in zip(times[:-1], times[1:])]
    if not np.allclose(steps, 1):
        raise ValueError("Expected regular daily sampling in the native calendar.")
    if frequency_smoothing_cpd <= 0 or wavenumber_smoothing <= 0:
        raise ValueError("Both spectral smoothing widths must be positive.")
    anomaly = tw.remove_mean_and_linear_trend(
        tw.remove_harmonics_of_seasonal_cycle(field, num_harmonics=we.SEASONAL_HARMONICS))
    taper = tw.apply_window(xr.ones_like(field.isel(lon=0)), window="tukey", pct=we.TAPER_FRACTION)
    tapered = tw.apply_window(anomaly, window="tukey", pct=we.TAPER_FRACTION)
    nt, nx = field.shape
    power = abs(np.fft.fft2(tapered.values))**2/(nt*nx*float((taper**2).mean()))
    # Periodic boundaries are appropriate for the discrete FFT frequency torus.
    power = gaussian_filter(power, sigma=(frequency_smoothing_cpd*nt, wavenumber_smoothing), mode="wrap")
    # Hermitian power is necessary for a real Gaussian field. Do not force every
    # realization to the same variance: that would condition the Gaussian draws.
    power = (power + power[np.ix_((-np.arange(nt)) % nt, (-np.arange(nx)) % nx)])/2
    power[0, 0] = 0
    return xr.DataArray(power, dims=("fft_frequency", "fft_wavenumber"),
        coords=dict(fft_frequency=np.fft.fftfreq(nt), fft_wavenumber=np.fft.fftfreq(nx)*nx),
        name="background_power", attrs=dict(units="(W m-2)^2", ordering="unshifted NumPy FFT; not WK plotting signs",
        frequency_smoothing_cpd=float(frequency_smoothing_cpd), wavenumber_smoothing=float(wavenumber_smoothing),
        expected_variance_W2_m4=float(power.mean()),
        null_model="stationary zero-mean Gaussian; smoothed projected-anomaly spectrum"))


def draw_background(power, template, rng):
    """Draw a real, correlated Gaussian field on template's native coordinates."""
    values = np.asarray(power)
    if values.shape != template.shape or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Power must be finite, nonnegative and match template shape.")
    nt, nx = values.shape
    partner = values[np.ix_((-np.arange(nt)) % nt, (-np.arange(nx)) % nx)]
    if not np.allclose(values, partner):
        raise ValueError("Power must have Hermitian symmetry for a real-valued draw.")
    white = rng.standard_normal(template.shape)
    colored = np.fft.ifft2(np.fft.fft2(white)*np.sqrt(values)).real
    return template.copy(data=colored).rename("gaussian_background")


def estimate_wk99_background_power(components):
    """Estimate a shared WK99-style background for projected daily OLR.

    Input: symmetric/antisymmetric fields in W m-2 on (time, lon). Return a
    Dataset containing segment spectra, the 1-2-1 background and full-record
    NumPy-ordered power compatible with ``draw_background``. This adapts WK99
    to the catalog's meridional means, rather than averaging latitude spectra.

    Use 96-day segments overlapping by 60 days, average the two component
    powers, then apply xr_ccew's 1-2-1 smoother: 10/20/40 wavenumber passes
    below .2/.3/above .3 cpd, followed by 10 frequency passes. WK99 specifies
    10--40 wavenumber passes but not the breakpoints; .2/.3 cpd follow the
    corresponding transitions in NCAR's implementation. Its additional five-
    pass rule below .1 cpd is not used here. Zero frequency is retained for
    synthesis; daily Nyquist is .5 cpd. These adaptations are saved explicitly.
    """
    symmetric, antisymmetric = xr.align(components.symmetric, components.antisymmetric, join="exact")
    for field in (symmetric, antisymmetric):
        if field.dims != ("time", "lon") or field.attrs.get("units") != "W m-2":
            raise ValueError("Expected projected symmetric/antisymmetric OLR in W m-2 on (time, lon).")
        if not np.isfinite(field.values).all() or field.sizes["time"] <= 2*we.REFERENCE_EXCLUDE_DAYS:
            raise ValueError("Expected finite projected records longer than the reference exclusions.")
    nt, nx = symmetric.shape
    if nx < 2 or not np.allclose(np.diff(symmetric.lon.values), 360/nx):
        raise ValueError("Expected an increasing, regular global longitude grid.")
    times = symmetric.time.values
    steps = (np.diff(times)/np.timedelta64(1, "D") if np.issubdtype(times.dtype, np.datetime64)
             else [(b-a).total_seconds()/86400 for a,b in zip(times[:-1],times[1:])])
    if not np.allclose(steps, 1):
        raise ValueError("WK99 calibration requires regular daily sampling in the native calendar.")

    spectra = {}
    for name,field in (("symmetric", symmetric), ("antisymmetric", antisymmetric)):
        spectra[name] = tw.segment_averaged_spectrum(field, segment_days=96, overlap_days=60,
            num_harmonics=we.SEASONAL_HARMONICS, detrend_segments=True,
            window="tukey", window_pct=we.TAPER_FRACTION, average_over_lat=False).sortby("frequency").sortby("zonal_wavenumber")
    pooled = (spectra["symmetric"]+spectra["antisymmetric"])/2
    positive = pooled.sel(frequency=slice(0, .5))
    frequency = positive.frequency.values
    blocks = []
    for lower,upper,passes in ((0.,.2,10), (.2,.3,20), (.3,np.inf,40)):
        block = positive.isel(frequency=np.flatnonzero((frequency>=lower)&(frequency<upper)))
        if block.sizes["frequency"]:
            blocks.append(tw.background_spectrum(block, passes_wavenumber=passes, passes_frequency=0))
    smoothed = tw.background_spectrum(xr.concat(blocks,dim="frequency").sortby("frequency"),
                                      passes_wavenumber=0, passes_frequency=10)

    # xr_ccew labels frequency with the opposite sign to NumPy's temporal FFT.
    # Interpolate only in frequency, reflecting k for the negative WK half-plane.
    # No interpolation in k is needed: segment and full-record grids share it.
    raw_f = np.fft.fftfreq(nt)
    raw_k = np.fft.fftfreq(nx)*nx
    k_order = [int(np.argmin(abs(smoothed.zonal_wavenumber.values-k))) for k in raw_k]
    positive_values = smoothed.values[:,k_order]
    values = np.column_stack([np.interp(abs(raw_f), frequency, positive_values[:,j]) for j in range(nx)])
    reflected_k = (-np.arange(nx)) % nx
    values[raw_f>0] = values[raw_f>0][:,reflected_k]
    # Segment powers are variance per bin, whereas the generator expects
    # E[|FFT|^2]/N. Convert via the segment frequency spacing 1/96 cpd.
    values *= 96*nx
    values = (values+values[np.ix_((-np.arange(nt))%nt,reflected_k)])/2
    values[0,0] = 0
    target_variance = float(pooled.sum())
    interpolated_variance = float(values.mean())
    if target_variance<=0 or interpolated_variance<=0:
        raise ValueError("A positive pooled background variance is required.")
    # Regridding and endpoint handling alter the integral slightly. Conserve
    # pooled segment variance in the target spectrum, not in individual draws.
    scale = target_variance/interpolated_variance
    values *= scale
    power = xr.DataArray(values,dims=("fft_frequency","fft_wavenumber"),
        coords=dict(fft_frequency=raw_f,fft_wavenumber=raw_k),name="background_power",
        attrs=dict(units="(W m-2)^2",ordering="unshifted NumPy FFT; not WK plotting signs",
                   expected_variance_W2_m4=target_variance,null_model="stationary Gaussian; shared WK99-style 1-2-1 spectrum"))
    result = xr.Dataset(dict(background_power=power,segment_power_symmetric=spectra["symmetric"],
        segment_power_antisymmetric=spectra["antisymmetric"],pooled_segment_power=pooled,
        smoothed_positive_background=smoothed))
    result.attrs = dict(background_method="wk99_121",segment_days=96,overlap_days=60,
        n_segments=int(spectra["symmetric"].attrs["n_segments"]),passes_frequency=10,
        passes_wavenumber="10 for f<0.2; 20 for 0.2<=f<0.3; 40 for f>=0.3 cpd",
        symmetry_pooling="arithmetic mean of symmetric and antisymmetric projected-field powers; common background",
        projection="power of catalog meridional means, not mean of latitude powers",
        boundaries="xr_ccew edge-padded 1-2-1 smoothing on nonnegative frequencies including zero and daily Nyquist",
        interpolation="linear in |frequency|; WK-to-NumPy sign conversion; Hermitian reflection; DC zero",
        target_variance_W2_m4=target_variance,variance_before_normalization_W2_m4=interpolated_variance,
        target_spectrum_rescaling=scale,
        reference="https://atmos.washington.edu/~breth/classes/AS552/papers/WheelerKiladis1999.pdf",
        breakpoint_reference="https://github.com/NCAR/ncl/blob/develop/ni/src/examples/gsun/diagnostics_cam.ncl")
    return result


def calibrate_event_peaks(field, wave_type, *, n_realizations, seed=0, quantile=0.99,
                          frequency_smoothing_cpd=1/96, wavenumber_smoothing=2.0):
    """Calibrate the peak envelope of detected candidates, in W m-2.

    Half the realizations set a pooled null-candidate peak quantile; independent
    remaining realizations estimate false candidates per whole source-length
    record. This is neither a pixel significance threshold nor 1% catalog FDR.
    No sinusoid fitting or existing catalog modifications occur here.
    Returns power, null-events/realizations/observed-events tables and a report.
    """
    if n_realizations < 4 or not 0 < quantile < 1:
        raise ValueError("Use at least four realizations and 0 < quantile < 1.")
    profile = we.WAVE_PROFILES[wave_type]
    power = estimate_background_power(field, frequency_smoothing_cpd=frequency_smoothing_cpd,
                                      wavenumber_smoothing=wavenumber_smoothing)
    rng = np.random.default_rng(seed)
    null_rows, summaries = [], []
    n_training = n_realizations//2
    observed = None
    observed_sigma = None
    for index in range(-1, n_realizations):
        source = field if index == -1 else draw_background(power, field, rng)
        # Match catalog preprocessing and filtering after the commuting latitude
        # projection. Filtered correlated noise is allowed to create packets.
        anomaly = tw.apply_window(tw.remove_mean_and_linear_trend(
            tw.remove_harmonics_of_seasonal_cycle(source, num_harmonics=we.SEASONAL_HARMONICS)),
            window="tukey", pct=we.TAPER_FRACTION)
        filtered = tw.filter_field(anomaly, profile, apply_profile_symmetry=False)
        events, sigma = we.detect_events(filtered)
        if index == -1:
            observed, observed_sigma = pd.DataFrame(events), sigma
            continue
        split = "calibration" if index < n_training else "validation"
        for event in events:
            null_rows.append(dict(event, realization=index, split=split, reference_sigma_W_m2=sigma))
        summaries.append(dict(realization=index, split=split, candidates=len(events), reference_sigma_W_m2=sigma,
                              maximum_candidate_peak_W_m2=max((e["peak_amplitude_W_m2"] for e in events), default=0.0)))
        print(f"{wave_type} null {index+1}/{n_realizations}: {len(events)} candidates", flush=True)
    null = pd.DataFrame(null_rows)
    if null.empty or not (null.split == "calibration").any():
        raise RuntimeError("No null calibration candidates; increase record length or inspect the null spectrum.")
    training = null.loc[null.split == "calibration", "peak_amplitude_W_m2"].to_numpy()
    threshold = float(np.quantile(training, quantile, method="higher"))
    realizations = pd.DataFrame(summaries)
    counts = null.loc[null.peak_amplitude_W_m2 >= threshold].groupby("realization").size()
    realizations["above_threshold"] = realizations.realization.map(counts).fillna(0).astype(int)
    maxima = realizations.loc[realizations.split == "calibration", "maximum_candidate_peak_W_m2"].to_numpy()
    if not observed.empty:
        peaks = observed.peak_amplitude_W_m2.to_numpy()
        observed["reference_sigma_W_m2"] = observed_sigma
        observed["exceeds_null_candidate_peak_threshold"] = peaks >= threshold
        ordered = np.sort(training)
        observed["null_candidate_peak_tail_fraction"] = (len(ordered)-np.searchsorted(ordered, peaks, side="left"))/len(ordered)
        # Whole-record maxima account for searching the entire time/longitude
        # domain. Their p-value resolution depends on independent realizations,
        # not the much larger count of correlated pixels or candidate peaks.
        observed["record_max_p_mc"] = [(1+np.count_nonzero(maxima >= peak))/(len(maxima)+1) for peak in peaks]
    heldout = realizations.loc[realizations.split == "validation"]
    report = dict(wave_type=wave_type, quantile=float(quantile), seed=seed, n_realizations=n_realizations,
        calibration_realizations=n_training, validation_realizations=len(heldout),
        calibration_candidates=len(training), candidate_peak_threshold_W_m2=threshold,
        interpretation="pooled peak quantile conditional on current candidate detector; not catalog FDR or record-wise significance",
        heldout_false_candidates_per_record_mean=float(heldout.above_threshold.mean()),
        heldout_false_candidates_per_record=heldout.above_threshold.tolist(),
        record_max_p_resolution=1/(n_training+1), observed_candidates=len(observed),
        observed_above_threshold=int(observed.exceeds_null_candidate_peak_threshold.sum()) if not observed.empty else 0,
        observed_reference_sigma_W_m2=observed_sigma, record_samples=field.sizes["time"],
        calendar=str(field.time.dt.calendar), background=dict(power.attrs),
        limitations="stationary Gaussian null; smoothing sensitivity, non-Gaussian tails, local/seasonal variance and injection recovery remain to be assessed")
    return dict(power=power, null_events=null, realizations=realizations, observed_events=observed, report=report)


def draw_marginal_background(power, template, rng, *, max_iterations=80, spectral_tolerance=0.03):
    """Generate a 2-D IAAFT surrogate with an empirical marginal distribution.

    Alternately impose Fourier amplitudes from the *smooth background* and the
    sorted values of the projected anomaly. This adapts Schreiber & Schmitz's
    IAAFT iteration to time/longitude, not independent series at each longitude.
    The empirical marginal is exact; the spectrum is approximate. Return the
    surrogate and a convergence report. A failed report must not be treated as
    a matched-spectrum null. This does not preserve local/seasonal variance or
    establish that the observed marginal consists only of background convection.
    """
    if max_iterations < 1 or not 0 < spectral_tolerance < 1:
        raise ValueError("Require positive iterations and 0 < spectral_tolerance < 1.")
    ordered = np.sort(template.values.ravel())
    if not np.isfinite(ordered).all() or ordered.std() <= 0:
        raise ValueError("Expected a finite, nonconstant empirical anomaly distribution.")
    values = np.asarray(power, dtype=float).copy()
    values[0, 0] = 0
    # Identical variance is necessary for compatible spectral/marginal constraints.
    scale = ordered.var()/values.mean()
    values *= scale
    amplitude = np.sqrt(values*template.size)
    norm = np.linalg.norm(amplitude)
    initial = draw_background(values, template, rng).values
    ranks = np.argsort(initial, axis=None)
    current = np.empty(template.size)
    current[ranks] = ordered
    current = current.reshape(template.shape)
    error = np.inf
    for iteration in range(1, max_iterations+1):
        transform = np.fft.fft2(current-current.mean())
        phase = np.divide(transform, abs(transform), out=np.ones_like(transform), where=abs(transform)>0)
        spectral_match = np.fft.ifft2(amplitude*phase).real
        current.ravel()[np.argsort(spectral_match, axis=None)] = ordered
        error = float(np.linalg.norm(abs(np.fft.fft2(current-current.mean()))-amplitude)/norm)
        if error <= spectral_tolerance:
            break
    report = dict(iterations=iteration, spectral_amplitude_relative_error=error,
                  spectral_tolerance=spectral_tolerance, converged=bool(error<=spectral_tolerance),
                  target_variance_rescaling=float(scale), marginal="exact global empirical anomaly distribution")
    return template.copy(data=current).rename("marginal_background"), report


def distribution_diagnostics(field):
    """Describe signed anomalies globally and in season/longitude subsets.

    Returns descriptive moments/tails, not iid-sample significance tests. Region
    edges are explicit here; every region retains the catalog's meridional mean.
    These diagnostics cannot distinguish wave tails from background tails alone.
    """
    from scipy.stats import kurtosis, skew

    regions = {"global": np.ones(field.sizes["lon"], dtype=bool),
               "Indian_40-100E": (field.lon.values>=40)&(field.lon.values<100),
               "Maritime_100-160E": (field.lon.values>=100)&(field.lon.values<160),
               "Pacific_160-280E": (field.lon.values>=160)&(field.lon.values<280),
               "Atlantic_Africa": (field.lon.values>=280)|(field.lon.values<40)}
    seasons = np.asarray(field.time.dt.season)
    rows = []
    for region, longitudes in regions.items():
        if not longitudes.any():
            continue
        for season in ("all", "DJF", "MAM", "JJA", "SON"):
            times = np.ones(field.sizes["time"], dtype=bool) if season=="all" else seasons==season
            values = field.values[np.ix_(times, longitudes)].ravel()
            if values.size < 8:
                continue
            mean, std = float(values.mean()), float(values.std())
            quantiles = np.quantile(values, [.001, .01, .5, .99, .999])
            rows.append(dict(region=region, season=season, samples=values.size, mean_W_m2=mean, sigma_W_m2=std,
                skewness=float(skew(values)), excess_kurtosis=float(kurtosis(values)),
                q001_W_m2=quantiles[0], q01_W_m2=quantiles[1], median_W_m2=quantiles[2],
                q99_W_m2=quantiles[3], q999_W_m2=quantiles[4],
                lower_3sigma_fraction=float(np.mean(values < mean-3*std)),
                upper_3sigma_fraction=float(np.mean(values > mean+3*std))))
    return rows


def spectrum_diagnostics(field, power):
    """Compare realized and target power in broad absolute-frequency/k bins.

    Bin ratios describe covariance fidelity; a Gaussian realization is not
    expected to equal its target periodogram bin by bin. Input uses the unshifted
    NumPy FFT ordering returned by estimate_background_power.
    """
    measured = abs(np.fft.fft2(field.values-field.values.mean()))**2/field.size
    f = abs(np.fft.fftfreq(field.sizes["time"]))
    k = abs(np.fft.fftfreq(field.sizes["lon"])*field.sizes["lon"])
    rows = []
    frequency_edges = [0, 1/90, 1/30, .1, .2, .4, .500001]
    k_edges = [0, 2, 5, 10, 20, 60, field.sizes["lon"]/2+1]
    for flo, fhi in zip(frequency_edges[:-1], frequency_edges[1:]):
        for klo, khi in zip(k_edges[:-1], k_edges[1:]):
            mask = (f[:, None]>=flo)&(f[:, None]<fhi)&(k[None, :]>=klo)&(k[None, :]<khi)
            expected = float(np.asarray(power)[mask].sum())
            if expected > 0:
                rows.append(dict(f_min_cpd=flo, f_max_cpd=fhi, abs_K_min=klo, abs_K_max=khi,
                                 bins=int(mask.sum()), measured_over_target=float(measured[mask].sum()/expected)))
    return rows


def synthetic_wave_packet(template, wave_type, K_global, h_eff_m, *, center_index, center_lon,
                          time_sigma_periods, lon_sigma_wavelengths, amplitude_W_m2=1.0, profile=None):
    """Return an amplitude-modulated carrier and its known resting-wave parameters.

    Envelope widths scale with carrier period/wavelength. Its translation uses
    the dry branch's group velocity, distinct from phase speed for dispersive
    waves. It is a controlled modulated carrier, not an exact nonlinear CCEW
    solution. Longitude is periodic and time uses native daily sample positions.
    Supply a profile to use a different candidate box or an additional branch.
    """
    from xr_ccew.matsuno import profile_curve_frequency

    if profile is None:
        profile = we.WAVE_PROFILES[wave_type]
    K = float(K_global)
    curve_options = dict(meridional_mode_number=profile.meridional_mode_number)
    frequency = float(profile_curve_frequency(profile.curve_name, K, h_eff_m, **curve_options))
    step = 1e-3
    derivative = float((profile_curve_frequency(profile.curve_name, K+step, h_eff_m, **curve_options)
                        - profile_curve_frequency(profile.curve_name, K-step, h_eff_m, **curve_options))/(2*step))
    if not np.isfinite(frequency) or frequency <= 0:
        raise ValueError("The requested carrier has no positive frequency on this branch.")
    period = 1/frequency
    time_sigma = time_sigma_periods*period
    lon_sigma = lon_sigma_wavelengths*360/abs(K)
    group_deg_day = 360*derivative
    lag = np.arange(template.sizes["time"])-center_index
    x = (template.lon.values-center_lon+180) % 360-180
    envelope_x = (x[None, :]-group_deg_day*lag[:, None]+180) % 360-180
    envelope = np.exp(-.5*(lag[:, None]/time_sigma)**2-.5*(envelope_x/lon_sigma)**2)
    values = amplitude_W_m2*envelope*np.cos(K*np.deg2rad(x)[None, :]-2*np.pi*frequency*lag[:, None]-.4)
    truth = dict(wave_type=wave_type, K_global=K, h_eff_m=h_eff_m, frequency_cpd=frequency, period_days=period,
                 time_sigma_days=time_sigma, lon_sigma_deg=lon_sigma, group_speed_deg_day=group_deg_day,
                 phase_speed_deg_day=360*frequency/K, center_index=center_index, center_lon_deg=center_lon,
                 center_time=str(template.time.values[center_index]), input_amplitude_W_m2=amplitude_W_m2,
                 carrier_inside_rectangular_bounds=bool(profile.k_min<=K<=profile.k_max
                     and profile.frequency_min<=frequency<=profile.frequency_max))
    return template.copy(data=values).rename("injected_packet"), truth


def filter_projected_olr(field, wave_type, *, profile=None):
    """Preprocess/filter projected (time, lon) OLR; optionally select a new profile.

    A supplied profile overrides the legacy catalog mask while retaining the
    same seasonal removal, detrending and taper used by the event detector.
    """
    anomaly = tw.apply_window(tw.remove_mean_and_linear_trend(
        tw.remove_harmonics_of_seasonal_cycle(field, num_harmonics=we.SEASONAL_HARMONICS)),
        window="tukey", pct=we.TAPER_FRACTION)
    if profile is None:
        profile = we.WAVE_PROFILES[wave_type]
    return tw.filter_field(anomaly, profile, apply_profile_symmetry=False)


def record_max_pvalues(null_maxima, scores):
    """Monte Carlo p-values accounting for a complete source/wave record search.

    Null maxima are from independent full-length Gaussian records passed through
    the same detector. Ties count against rejection. These p-values condition on
    the specified null spectrum; they do not measure catalog false-discovery rate.
    """
    maxima = np.sort(np.asarray(null_maxima, dtype=float))
    values = np.asarray(scores, dtype=float)
    if (maxima.ndim != 1 or not maxima.size or not np.isfinite(maxima).all()
            or (maxima < 0).any() or not np.isfinite(values).all() or (values < 0).any()):
        raise ValueError("Expected finite nonnegative scores and a nonempty vector of record maxima.")
    count_ge = len(maxima)-np.searchsorted(maxima, values, side="left")
    return (1+count_ge)/(len(maxima)+1)


def record_max_threshold(null_maxima, alpha=0.01, confidence=0.95):
    """Return a strict 'peak > threshold' rule equivalent to p_mc <= alpha.

    The order-statistic interval quantifies Monte Carlo uncertainty in the target
    (1-alpha) null-maximum quantile, not background-model uncertainty. An upper
    interval beyond the sampled maxima is explicitly unbounded (None).
    """
    from scipy.stats import binom

    maxima = np.sort(np.asarray(null_maxima, dtype=float))
    record_max_pvalues(maxima, [0.0])  # Validate the shared score definition.
    if not 0 < alpha < 1 or not 0 < confidence < 1:
        raise ValueError("Require 0 < alpha, confidence < 1.")
    n = len(maxima)
    allowed_ge = int(np.floor(alpha*(n+1)+1e-12))-1
    result = dict(alpha=float(alpha), calibration_records=n, minimum_p=1/(n+1),
                  available=allowed_ge>=0, comparison="strictly greater than threshold")
    if allowed_ge < 0:
        return dict(result, threshold_W_m2=None, quantile_CI_low_W_m2=None, quantile_CI_high_W_m2=None)
    rank = n-allowed_ge  # One-based order statistic; equality does not pass.
    lower_rank = int(binom.ppf((1-confidence)/2,n,1-alpha))
    upper_rank = int(binom.ppf(1-(1-confidence)/2,n,1-alpha))+1
    return dict(result, threshold_W_m2=float(maxima[rank-1]), order_statistic_rank=rank,
                quantile_CI_low_W_m2=float(maxima[lower_rank-1]) if lower_rank>0 else 0.0,
                quantile_CI_high_W_m2=float(maxima[upper_rank-1]) if upper_rank<=n else None,
                confidence=float(confidence), interval_upper_unbounded=upper_rank>n)


def match_packet_peaks(events, field, truth, threshold):
    """Count catalog peaks within one sigma along a known moving injected packet.

    Return candidate/significant counts and the largest matching envelope peak.
    The threshold is strict, consistent with record_max_threshold. Applying this
    also to the uninjected background exposes pre-existing chance associations.
    """
    if not events:
        return dict(nearby_candidates=0, nearby_significant=0, nearby_maximum_W_m2=0.0)
    indices = field.get_index("time").get_indexer([event["peak_time"] for event in events])
    lag = indices-truth["center_index"]
    lon = np.array([event["peak_lon_deg"] for event in events])
    x = (lon-truth["center_lon_deg"]-truth["group_speed_deg_day"]*lag+180) % 360-180
    nearby = (abs(lag)<=truth["time_sigma_days"]) & (abs(x)<=truth["lon_sigma_deg"])
    peaks = np.array([event["peak_amplitude_W_m2"] for event in events])[nearby]
    return dict(nearby_candidates=len(peaks), nearby_significant=int(np.count_nonzero(peaks>threshold)),
                nearby_maximum_W_m2=float(peaks.max()) if peaks.size else 0.0)
