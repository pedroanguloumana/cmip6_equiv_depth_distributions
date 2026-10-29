"""Custom wave profiles with very, very wide equivalent depth ranges
"""
import xr_ccew as tw


GENERAL_KELVIN = tw.WaveProfile(
    name="KW",
    aliases=("Kelvin", "Kelvin wave"),
    k_min=1,
    k_max=14,
    frequency_min=1.0 / 30.0,
    frequency_max=1.0 / 2.5,
    equivalent_depth_min=8,
    equivalent_depth_max=90,
    symmetry="symmetric",
    direction="eastward",
    curve_name="kelvin",
    description="Convectively coupled Kelvin wave.",
    references=(REFERENCE_WK99,),
)