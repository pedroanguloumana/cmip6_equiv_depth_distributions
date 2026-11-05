"""Geographical OLR variance and regression structures of filtered waves."""

import numpy as np
import xarray as xr
import xr_ccew as tw


def regress_olr_on_filtered_index(anomaly, reference):
    """Regress OLR anomalies onto a standardized convective reference index.

    Inputs are finite daily OLR [W m-2] on (time, lat, lon) and a wave-filtered
    reference [W m-2] on the identical time coordinate. Return the zero-lag
    OLR anomaly [W m-2] for a -1 population-standard-deviation reference
    excursion, its correlation, and the reference standard deviation.
    An intercept is fitted; no wave filtering or symmetry is imposed on anomaly.
    """
    if set(anomaly.dims) != {"time", "lat", "lon"} or reference.dims != ("time",):
        raise ValueError("Expected anomaly(time, lat, lon) and reference(time).")
    anomaly, reference = xr.align(anomaly, reference, join="exact")
    if anomaly.attrs.get("units") != "W m-2" or reference.attrs.get("units") != "W m-2":
        raise ValueError("Both fields must have units W m-2.")
    values = np.asarray(anomaly.transpose("time", "lat", "lon"), dtype=float)
    index = np.asarray(reference, dtype=float)
    if len(index) < 3 or not np.isfinite(values).all() or not np.isfinite(index).all():
        raise ValueError("Regression requires at least three finite paired daily samples.")
    centered = index - index.mean()
    sigma = float(np.std(centered, ddof=0))
    if sigma <= 0:
        raise ValueError("The reference index has zero variance.")

    # A centered predictor makes the target's fitted intercept drop out of the
    # cross product. A positive convection index denotes negative filtered OLR.
    response = -np.einsum("t,tij->ij", centered, values, optimize=True) / (len(index) * sigma)
    target_sigma = np.std(values, axis=0, ddof=0)
    correlation = np.divide(response, target_sigma, out=np.full_like(response, np.nan),
                            where=target_sigma > 0)
    if np.nanmax(np.abs(correlation)) > 1 + 1e-10:
        raise RuntimeError("Regression correlation exceeds one.")
    result = xr.Dataset(
        data_vars={
            "regression": (("lat", "lon"), response),
            "correlation": (("lat", "lon"), correlation),
            "reference_sigma": sigma,
        },
        coords={"lat": anomaly.lat, "lon": anomaly.lon},
    )
    result.regression.attrs.update(units="W m-2", reference_excursion_sigma=-1.0,
                                   long_name="OLR anomaly for a one-sigma convective reference excursion")
    result.correlation.attrs.update(units="1", long_name="Correlation with negative wave-filtered OLR index")
    result.reference_sigma.attrs.update(units="W m-2", ddof=0)
    return result


def calculate_olr_wave_maps(olr, profiles, *, seasonal_harmonics=3, taper_fraction=0.1):
    """Return variance, maximum-variance bases, and zero-lag OLR regressions.

    olr is finite daily OLR [W m-2] on a global regular (time, lat, lon) grid
    with mirrored latitudes. profiles maps wave names to xr_ccew WaveProfiles.
    Remove seasonal harmonics, then the mean/linear trend; taper only the
    predictor input. Statistics exclude the tapered ends. Base points maximize
    each wave's interior population variance over the supplied tropical domain.
    For the tied hemispheric maxima of symmetry-projected filters, use the
    northern point. The mapped regression target retains its full latitude
    structure and is not wave-filtered or tapered.
    """
    if set(olr.dims) != {"time", "lat", "lon"}:
        raise ValueError("Expected daily OLR on (time, lat, lon).")
    olr = olr.transpose("time", "lat", "lon")
    nt, nlat, nx = olr.shape
    if olr.attrs.get("units") != "W m-2" or not np.isfinite(olr.values).all():
        raise ValueError("Expected finite OLR in W m-2.")
    if (nlat < 2 or not np.all(np.diff(olr.lat) > 0)
            or not np.allclose(olr.lat, -olr.lat.values[::-1], atol=1e-9)
            or nx < 4 or not np.allclose(np.diff(olr.lon), 360 / nx)):
        raise ValueError("Require mirrored increasing latitudes and regular global longitudes.")
    dates = olr.time.values
    steps = (np.diff(dates) / np.timedelta64(1, "D") if np.issubdtype(dates.dtype, np.datetime64)
             else np.array([(b-a).total_seconds()/86400 for a, b in zip(dates[:-1], dates[1:])]))
    if not np.allclose(steps, 1, rtol=0, atol=1e-9):
        raise ValueError("Require uninterrupted daily sampling in the native calendar.")
    if not 0 < taper_fraction < 1:
        raise ValueError("taper_fraction must lie strictly between zero and one.")
    edge = int(np.ceil(taper_fraction * nt / 2))
    if nt - 2 * edge < 300:
        raise ValueError("Require at least 300 interior daily samples after taper exclusion.")
    interior = slice(edge, nt-edge)
    print("  removing seasonal harmonics and mean/linear trend", flush=True)
    anomaly = tw.remove_mean_and_linear_trend(
        tw.remove_harmonics_of_seasonal_cycle(olr, num_harmonics=seasonal_harmonics)).load()
    tapered = tw.apply_window(anomaly, window="tukey", pct=taper_fraction)
    maps = []
    for wave, profile in profiles.items():
        print(f"  filtering and regressing {wave}", flush=True)
        filtered = tw.filter_field(tapered, profile).transpose("time", "lat", "lon")
        variance = filtered.isel(time=interior).var("time", ddof=0)
        search = variance.sel(lat=slice(0, None)) if profile.symmetry in ("symmetric", "antisymmetric") else variance
        iy, ix = np.unravel_index(int(np.argmax(search.values)), search.shape)
        base_lat, base_lon = float(search.lat[iy]), float(search.lon[ix])
        maximum = float(search.values[iy, ix])
        if not np.isclose(maximum, float(variance.max()), rtol=1e-10, atol=1e-10):
            raise RuntimeError("Northern symmetry tie does not reproduce the full-domain maximum.")
        reference = filtered.sel(lat=base_lat, lon=base_lon, drop=True).rename("filtered_reference")
        result = regress_olr_on_filtered_index(anomaly.isel(time=interior), reference.isel(time=interior))
        if not np.isclose(float(result.reference_sigma)**2, maximum, rtol=1e-10, atol=1e-10):
            raise RuntimeError("Reference variance does not equal the selected map maximum.")
        result["filtered_variance"] = variance
        result.filtered_variance.attrs.update(units="W2 m-4", ddof=0,
                                               long_name="Interior daily wave-filtered OLR variance")
        result["filtered_reference"] = reference
        result["base_latitude"] = xr.DataArray(base_lat, attrs={"units": "degrees_north"})
        result["base_longitude"] = xr.DataArray(base_lon, attrs={"units": "degrees_east"})
        result["maximum_variance"] = xr.DataArray(maximum, attrs={"units": "W2 m-4"})
        result["filter_k_min"] = profile.k_min
        result["filter_k_max"] = profile.k_max
        result["filter_frequency_min"] = xr.DataArray(profile.frequency_min, attrs={"units": "day-1"})
        result["filter_frequency_max"] = xr.DataArray(profile.frequency_max, attrs={"units": "day-1"})
        result = result.expand_dims(wave=[wave]).assign_coords(wave_symmetry=("wave", [profile.symmetry]))
        maps.append(result)
        print(f"    base {base_lat:.1f} N, {base_lon:.1f} E; variance {maximum:.3f} W2 m-4", flush=True)
        del filtered

    combined = xr.concat(maps, dim="wave", combine_attrs="drop_conflicts")
    combined.attrs.update(
        calendar=str(olr.time.dt.calendar), record_start=str(dates[0]), record_end=str(dates[-1]),
        input_samples=nt, statistics_samples=nt-2*edge,
        statistics_start=str(dates[edge]), statistics_end=str(dates[nt-edge-1]),
        excluded_samples_per_end=edge, seasonal_harmonics=seasonal_harmonics,
        preprocessing="seasonal harmonics removed, then full-record mean/linear trend",
        predictor_taper=f"Tukey alpha={taper_fraction}; tapered ends excluded from all statistics",
        taper_fraction=taper_fraction, latitude_min=float(olr.lat.min()), latitude_max=float(olr.lat.max()),
        base_selection="separate maximum interior filtered variance per dataset/wave over supplied domain; northern hemisphere breaks symmetry ties",
        regression_target="untapered, unfiltered OLR anomaly, without symmetry projection",
        regression_scaling="negative one population-standard-deviation filtered reference excursion; fitted intercept removed",
        lag_days=0, time_selection="all interior daily samples; no event threshold or active-variance window",
        significance="descriptive maps; no significance test or spatial multiplicity correction",
    )
    return combined
