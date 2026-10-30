"""Fit local longitude--time sinusoids to filtered equatorial-wave events."""

import numpy as np
import xarray as xr
from scipy.optimize import least_squares


EARTH_RADIUS_M = 6_371_008.0
GRAVITY_M_S2 = 9.80665


def resting_wave_depth(wave_type, k_rad_m, omega_rad_s):
    """Infer resting-background depth from a fitted, signed (k, omega) pair.

    Supports Kelvin, n1ERW, MRG, and n0EIG. Inputs are angular wavenumber
    (rad/m) and positive angular frequency (rad/s). Return c0_m_s,
    h_eff_m, and dispersion_valid. A pair with no positive solution on its
    assigned branch returns NaN speeds/depths and a False validity flag.
    c0 is the gravity-speed parameter, not generally the phase speed omega/k.
    """
    if wave_type not in {"Kelvin", "n1ERW", "MRG", "n0EIG"}:
        raise ValueError(f"No depth inversion defined for {wave_type!r}")

    k, omega = float(k_rad_m), float(omega_rad_s)
    c0 = np.nan
    # Match the equatorial beta and rotation rate used by xr_ccew's curves.
    beta = 2 * 7.292e-5 / EARTH_RADIUS_M
    if np.isfinite(k) and np.isfinite(omega) and omega > 0:
        if wave_type == "Kelvin":
            if k > 0:
                c0 = omega / k
        elif wave_type in {"MRG", "n0EIG"}:
            denominator = beta + k * omega
            direction_ok = k < 0 if wave_type == "MRG" else k > 0
            if direction_ok and denominator > 0:
                c0 = omega**2 / denominator
        else:
            # The n=1 Matsuno relation is quadratic in c0. The larger positive
            # root selects ER; the smaller root would give a WIG interpretation.
            coefficient = k**2 + beta * k / omega
            if k < 0 and coefficient < 0:
                discriminant = (3 * beta)**2 + 4 * coefficient * omega**2
                if discriminant >= 0:
                    c0 = (-3 * beta - np.sqrt(discriminant)) / (2 * coefficient)

    valid = bool(np.isfinite(c0) and c0 > 0)
    return {
        "c0_m_s": float(c0) if valid else np.nan,
        "h_eff_m": float(c0**2 / GRAVITY_M_S2) if valid else np.nan,
        "dispersion_valid": valid,
    }


def fit_event_sinusoid(
    field: xr.DataArray,
    center_time,
    center_lon,
    *,
    k_grid,
    frequency_grid,
    half_window_days,
    half_width_deg,
    time_sigma_days,
    lon_sigma_deg,
):
    """Fit a weighted sinusoid around one event peak.

    The model is ``q = b + A*cos(theta) + B*sin(theta)``, where
    ``theta = K*longitude_offset_rad - 2*pi*f*time_offset_days``.

    Parameters
    ----------
    field : xarray.DataArray
        Real, finite filtered OLR in W m-2, with dimensions (time, lon),
        regularly sampled daily and covering all global longitudes. NumPy
        datetime coordinates and native CF calendars are both supported.
    center_time, center_lon
        Event peak time and longitude in degrees east.
    k_grid, frequency_grid : array-like
        Ascending coarse search grids: signed global wavenumber K (cycles
        per 360 degrees) and positive frequency f (cycles/day). Their
        endpoints also bound the continuous refinement. With this phase
        convention, negative K describes westward propagation.
    half_window_days, half_width_deg
        Half-widths of the rectangular fitting window, in days and degrees.
    time_sigma_days, lon_sigma_deg
        Gaussian weighting scales around the event center, in days and
        degrees. These set relative influence within the window, not errors
        on the fitted parameters.

    Returns
    -------
    result : dict
        Fitted K, f, dimensional k (rad/m) and omega (rad/s), signed phase
        speed (m/s), amplitude/offset (W m-2), phase (radians), weighted R2,
        sampled cycle counts, and optimizer/boundary flags. Phase follows
        ``q = offset + amplitude*cos(theta - phase_rad)``.
        For compatibility, ``h_eff_m`` retains the notebook's Kelvin-only
        conversion, phase_speed**2/g. For ER, MRG, and EIG use the fitted
        k and omega with that family's dispersion relation instead.
    diagnostic : dict
        Observed and fitted local DataArrays, Gaussian weights, coarse
        search grids, and their R2 surface. Coordinates are preserved;
        the local longitude coordinate is the offset from center_lon.
    """
    # Select by sample position: converting cftime dates to Gregorian dates would
    # distort a model's calendar (and fails outright on dates such as February 30).
    center_index = field.get_index("time").get_loc(center_time)
    if isinstance(center_index, slice) and center_index.stop - center_index.start == 1:
        center_index = center_index.start
    if not isinstance(center_index, (int, np.integer)):
        raise ValueError("center_time must identify exactly one daily sample.")
    first, stop = center_index - half_window_days, center_index + half_window_days + 1
    if first < 0 or stop > field.sizes["time"]:
        raise ValueError("The full event fitting window is not available.")
    patch = field.isel(time=slice(first, stop))
    center_time = field.time.values[center_index]
    if np.issubdtype(patch.time.dtype, np.datetime64):
        lag_days = (patch.time.values - center_time) / np.timedelta64(1, "D")
    else:
        lag_days = np.array([
            (value - center_time).total_seconds() / 86400.0
            for value in patch.time.values
        ])
    if not np.allclose(lag_days, np.arange(-half_window_days, half_window_days + 1)):
        raise ValueError("The fitting window must be regularly sampled once per day.")

    # Wrap before selecting longitudes: a packet at 359 degrees should include
    # nearby observations at 0 degrees. Sorting makes the local coordinate monotonic.
    relative_lon = (patch.lon.values - center_lon + 180) % 360 - 180
    order = np.argsort(relative_lon)
    order = order[np.abs(relative_lon[order]) <= half_width_deg]
    patch = patch.isel(lon=order).assign_coords(lon=relative_lon[order])
    lon_rad, lag = np.meshgrid(np.deg2rad(patch.lon.values), lag_days)

    # Emphasize the event neighborhood while fitting both signs of the OLR anomaly.
    # Multiplying residuals by sqrt(weight) makes their squares carry weight once.
    weights = np.exp(
        -0.5 * (lag / time_sigma_days) ** 2
        - 0.5 * (np.rad2deg(lon_rad) / lon_sigma_deg) ** 2
    )
    sqrt_weights = np.sqrt(weights).ravel()
    y = patch.values.ravel()
    weighted_y = sqrt_weights * y

    # Compare the wave model with the best constant-only model in the same window.
    # This weighted variance supplies the denominator of the reported R2.
    mean = np.average(y, weights=weights.ravel())
    null_sse = np.sum(weights.ravel() * (y - mean) ** 2)
    if null_sse <= 0:
        raise ValueError("A constant field cannot identify a wave.")

    # Optimize frequency as cycles per window, avoiding very different parameter
    # scales for K and cycles/day. Convert back to cycles/day after refinement.
    time_span_days = float(np.ptp(lag_days))
    x = lon_rad.ravel()
    scaled_time = lag.ravel() / time_span_days

    def residuals(parameters):
        # For fixed K and frequency the offset, cosine and sine coefficients are
        # linear unknowns. Solve them exactly at each nonlinear trial pair.
        K, cycles_per_window = parameters
        phase = K * x - 2 * np.pi * cycles_per_window * scaled_time
        design = np.column_stack([np.ones_like(phase), np.cos(phase), np.sin(phase)])
        weighted_design = sqrt_weights[:, None] * design
        coefficients = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)[0]

        # The squared norm is SSE/null_sse, hence weighted R2 = 1 - sum(residuals**2).
        return (weighted_design @ coefficients - weighted_y) / np.sqrt(null_sse)

    # Search the whole requested rectangle before local optimization. Refining
    # three promising starts reduces dependence on a single initial guess.
    coarse_pairs = np.array([
        (K, f * time_span_days) for K in k_grid for f in frequency_grid
    ])
    coarse_sse = np.array([np.sum(residuals(pair) ** 2) for pair in coarse_pairs])
    lower = [k_grid[0], frequency_grid[0] * time_span_days]
    upper = [k_grid[-1], frequency_grid[-1] * time_span_days]
    solutions = [
        least_squares(
            residuals,
            coarse_pairs[index],
            bounds=(lower, upper),
            ftol=1e-9,
            xtol=1e-9,
            gtol=1e-9,
            max_nfev=200,
        )
        for index in np.argsort(coarse_sse)[:3]
    ]
    solution = min(solutions, key=lambda fit: np.sum(fit.fun ** 2))

    # Reconstruct the signed wave field from the continuously refined parameters.
    # The solution is not restricted to integer K or to the coarse frequency grid.
    K, cycles_per_window = solution.x
    frequency = cycles_per_window / time_span_days
    phase = K * x - 2 * np.pi * cycles_per_window * scaled_time
    design = np.column_stack([np.ones_like(phase), np.cos(phase), np.sin(phase)])
    b, A, B = np.linalg.lstsq(
        sqrt_weights[:, None] * design, weighted_y, rcond=None
    )[0]
    model = (design @ np.array([b, A, B])).reshape(patch.shape)

    # K is global cycles/360 degrees; dimensional k and omega use radians.
    # omega/k keeps the propagation direction through the sign of k.
    omega = 2 * np.pi * frequency / 86400
    k = K / EARTH_RADIUS_M
    speed = omega / k
    at_edge = (
        np.isclose(solution.x, lower, rtol=0, atol=1e-4).any()
        or np.isclose(solution.x, upper, rtol=0, atol=1e-4).any()
    )

    result = dict(
        K_global=float(K),
        frequency_cpd=float(frequency),
        k_rad_m=float(k),
        omega_rad_s=float(omega),
        phase_speed_m_s=float(speed),
        # Valid as equivalent depth for resting Kelvin waves only. Other families
        # require their branch-specific inversion before squaring the gravity speed.
        h_eff_m=float(speed ** 2 / GRAVITY_M_S2),
        period_days=float(1 / frequency),
        amplitude_W_m2=float(np.hypot(A, B)),
        offset_W_m2=float(b),
        phase_rad=float(np.arctan2(B, A)),
        weighted_r2=float(1 - np.sum(solution.fun ** 2)),
        # A cycle count is positive for westward as well as eastward waves.
        spatial_cycles=float(abs(K) * np.ptp(patch.lon.values) / 360),
        temporal_cycles=float(frequency * time_span_days),
        fit_success=bool(solution.success),
        at_search_edge=bool(at_edge),
    )
    diagnostic = dict(
        observed=patch,
        fitted=patch.copy(data=model),
        weights=weights,
        k_grid=k_grid,
        frequency_grid=frequency_grid,
        coarse_r2=(1 - coarse_sse).reshape(len(k_grid), len(frequency_grid)),
    )
    return result, diagnostic
