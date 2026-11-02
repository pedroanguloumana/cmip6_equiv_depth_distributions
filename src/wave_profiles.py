"""Project wave profiles, including rectangular candidate-filter regions.

BOX_PROFILES defines the revised k/frequency limits and OLR symmetries, without
an equivalent-depth restriction. Calibration and event fitting use this collection.
"""
import xr_ccew as tw


# Bounds: signed global zonal wavenumber and frequency in cycles/day. In
# xr_ccew, an infinite depth bound disables the dispersion-curve restriction.
# Retain curve_name/mode so theory plotting and synthetic wave packets can
# still use the correct branch. These metadata do not curve the box masks.
KELVIN_BOX = tw.WaveProfile(
    name="KELVIN_BOX", k_min=1, k_max=15,
    frequency_min=1 / 14, frequency_max=1 / 2.5,
    equivalent_depth_min=0, equivalent_depth_max=float("inf"),
    symmetry="symmetric", direction="eastward", curve_name="kelvin",
    description="Rectangular Kelvin candidate filter; no equivalent-depth restriction.",
)

n1ERW_BOX = tw.WaveProfile(
    name="n1ERW_BOX", k_min=-15, k_max=-1,
    frequency_min=1 / 90, frequency_max=1 / 10,
    equivalent_depth_min=0, equivalent_depth_max=float("inf"),
    symmetry="symmetric", direction="westward", curve_name="er",
    meridional_mode_number=1,
    description="Rectangular n=1 equatorial Rossby candidate filter; no equivalent-depth restriction.",
)

MRG_BOX = tw.WaveProfile(
    name="MRG_BOX", k_min=-10, k_max=-0.25,
    frequency_min=1 / 10, frequency_max=1 / 2.25,
    equivalent_depth_min=0, equivalent_depth_max=float("inf"),
    symmetry="antisymmetric", direction="westward", curve_name="mrg",
    description="Rectangular MRG candidate filter; no equivalent-depth restriction.",
)

n0EIG_BOX = tw.WaveProfile(
    name="n0EIG_BOX", k_min=0.25, k_max=14,
    frequency_min=1 / 6, frequency_max=1 / 2.25,
    equivalent_depth_min=0, equivalent_depth_max=float("inf"),
    symmetry="antisymmetric", direction="eastward", curve_name="eig0",
    meridional_mode_number=0,
    description="Rectangular n=0 eastward inertio-gravity candidate filter; no equivalent-depth restriction.",
)

n1WIG_BOX = tw.WaveProfile(
    name="n1WIG_BOX", k_min=-15, k_max=-1,
    frequency_min=1 / 4, frequency_max=1 / 2,
    equivalent_depth_min=0, equivalent_depth_max=float("inf"),
    symmetry="symmetric", direction="westward", curve_name="wig",
    meridional_mode_number=1,
    description="Rectangular n=1 westward inertio-gravity candidate filter; no equivalent-depth restriction.",
)


# Existing wave labels plus n1WIG, for calibration and notebook loops.
BOX_PROFILES = {
    "Kelvin": KELVIN_BOX,
    "n1ERW": n1ERW_BOX,
    "n0EIG": n0EIG_BOX,
    "MRG": MRG_BOX,
    "n1WIG": n1WIG_BOX
}
