"""Check regression scaling and variance-maximum base-point selection."""

import numpy as np
import pytest
import xarray as xr

from src.wave_profiles import BOX_PROFILES
from src.wave_regression import calculate_olr_wave_maps, regress_olr_on_filtered_index


def test_regression_recovers_signed_one_sigma_structure_with_intercept():
    time = np.arange(400)
    signal = np.cos(2*np.pi*time/20)
    unrelated = np.sin(2*np.pi*time/10)
    pattern = np.array([[2.0, -3.0], [1.0, 4.0]])
    values = 250 + signal[:, None, None]*pattern + unrelated[:, None, None]
    anomaly = xr.DataArray(values, dims=("time", "lat", "lon"),
                           coords=dict(time=time, lat=[-5, 5], lon=[0, 180]), attrs={"units": "W m-2"})
    reference = xr.DataArray(17+3*signal, dims="time", coords=dict(time=time), attrs={"units": "W m-2"})
    result = regress_olr_on_filtered_index(anomaly, reference)
    np.testing.assert_allclose(result.regression, -pattern/np.sqrt(2), atol=1e-12)
    assert float(result.reference_sigma) == pytest.approx(3/np.sqrt(2))
    assert float(result.correlation[0, 0]) < 0
    with pytest.raises(ValueError, match="align|index|coordinate"):
        regress_olr_on_filtered_index(anomaly, reference.assign_coords(time=time+1))


def test_wave_map_selects_global_variance_maximum_and_preserves_target_asymmetry():
    time = np.arange(720)
    lat = np.array([-10.5, -0.5, 0.5, 10.5])
    lon = np.arange(0.5, 360, 360/48)
    carrier = np.cos(3*np.deg2rad(lon)[None, :] - 2*np.pi*time[:, None]/6)
    # The filtered predictor is symmetric; the unfiltered regression should
    # retain the deliberately stronger northern response of the target.
    amplitude = np.array([1., .4, .8, 3.])
    field = xr.DataArray(carrier[:, None, :]*amplitude[None, :, None],
                         dims=("time", "lat", "lon"),
                         coords=dict(time=np.datetime64("2000-01-01")+time, lat=lat, lon=lon),
                         attrs={"units": "W m-2"})
    maps = calculate_olr_wave_maps(field, {"Kelvin": BOX_PROFILES["Kelvin"]},
                                    seasonal_harmonics=3, taper_fraction=0.1)
    assert maps.base_latitude.item() == 10.5
    assert maps.maximum_variance.item() == pytest.approx(float(maps.filtered_variance.max()), rel=1e-12)
    assert maps.reference_sigma.item()**2 == pytest.approx(maps.maximum_variance.item(), rel=1e-12)
    np.testing.assert_allclose(maps.filtered_variance, maps.filtered_variance.isel(lat=slice(None, None, -1)), atol=1e-12)
    north = maps.regression.sel(lat=10.5).values
    south = maps.regression.sel(lat=-10.5).values
    np.testing.assert_allclose(north, 3*south, atol=1e-10)
    assert maps.attrs["statistics_samples"] == 648
