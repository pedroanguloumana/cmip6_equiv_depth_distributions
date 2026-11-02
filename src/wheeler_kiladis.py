"""Wheeler-Kiladis space-time spectral analysis.

"""

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import xr_ccew as tw
from xr_ccew.matsuno import (
    eig0_frequency, eig_frequency, er_frequency,
    kelvin_frequency, mrg_frequency, wig_frequency,
)
def spectral_power_pipeline(
    da,
    which="symmetric",
    *,
    lat_bound=15.0,
    segment_days=96,
    overlap_days=30,
    num_harmonics=3,
    window="tukey",
    window_pct=0.10,
    max_segments=None,
    time_dim="time",
    lat_dim="lat",
    lon_dim="lon",
):
    
    band = da.where(np.abs(da[lat_dim]) <= lat_bound, drop=True).sortby(lat_dim)
    component = tw.symmetric_antisymmetric_component(band, which, lat_dim=lat_dim)
    anomalies = tw.remove_mean_and_linear_trend(component, dim=time_dim)
    anomalies = tw.remove_harmonics_of_seasonal_cycle(
        anomalies, num_harmonics=num_harmonics, time_dim=time_dim
    )

    segments = tw.segment_data(
        anomalies,
        segment_days=segment_days,
        overlap_days=overlap_days,
        time_dim=time_dim,
    )
    if max_segments is not None:
        segments = segments[:max_segments]
    if not segments:
        raise ValueError(
            f"no complete {segment_days}-day segment in "
            f"{da[time_dim].min().item()}..{da[time_dim].max().item()}"
        )

    n_time = min(seg.sizes[time_dim] for seg in segments)

    powers = []
    for segment in segments:
        segment = segment.isel({time_dim: slice(0, n_time)})
        tapered = tw.apply_window(segment, dim=time_dim, window=window, pct=window_pct)
        powers.append(
            tw.power_spectrum(tapered, time_dim=time_dim, lon_dim=lon_dim)
        )

    power = (
        xr.concat(powers, dim="segment", coords="minimal", compat="override")
        .mean("segment")
        .sum(lat_dim)
    )

    power.name = f"{which}_power"
    power.attrs.update(
        {
            "component": which,
            "lat_bound": lat_bound,
            "segment_days": segment_days,
            "overlap_days": overlap_days,
            "n_segments": len(segments),
            "segment_length": n_time,
            "num_harmonics_removed": num_harmonics,
            "window": f"{window} ({window_pct:g} taper)",
            "long_name": f"{which} space-time power",
            "units": "UnNormed"
        }
    )
    return power


def add_matsuno_curves(ax, which='symmetric', he_list=(15, 25, 50)):
    _CURVES = {
            'symmetric': {
                'Kelvin': lambda k, he: kelvin_frequency(k, he),
                'n=1 ER': lambda k, he: er_frequency(k, he, n=1),
                'n=1 WIG': lambda k, he: wig_frequency(k, he, n=1),
                'n=1 EIG': lambda k, he: eig_frequency(k, he, n=1),
            },
            'antisymmetric': {
                'MRG': lambda k, he: mrg_frequency(k, he),
                'n=0 EIG': lambda k, he: eig0_frequency(k, he),
                'n=2 WIG': lambda k, he: wig_frequency(k, he, n=2),
                'n=2 EIG': lambda k, he: eig_frequency(k, he, n=2),
            },
        }
    k = np.linspace(-15, 15, 301)
    k = k[k != 0]
    for frequency_of in _CURVES[which].values():
        for he in he_list:
            ax.plot(k, frequency_of(k, he), color='red', lw=0.8)

    ax.set_ylim(0.01, 0.5)
    ax.set_xlim(-15, 15)
    return ax