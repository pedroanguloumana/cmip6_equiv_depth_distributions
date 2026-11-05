# Notebooks

`Explore_Catalog.ipynb` loads the complete catalog once and searches in
memory. Supporting notes are in [docs/](../docs/README.md). Older executed
snapshots are preserved in `archived/catalog-audit-2026-10-27/`.

`codex_segment_heff.ipynb` compares excess OLR power in 10--80 m dispersion
bands from individual 96-day WK segments, including Kelvin--ER comparisons
and finite-resolution diagnostics. Reproduction instructions are in
[`scripts/codex_segment_heff/README.md`](../scripts/codex_segment_heff/README.md).

`codex_wave_event_komega.ipynb` compares fitted 2-sigma event locations with
WK-style spectral contours, in five-wave panels for NOAA and Historical/AMIP
model means. It includes both fit selections, dimensional k/omega axes,
dispersion curves and an overlap summary. Reproduction instructions are in
[`scripts/codex_wave_event_komega/README.md`](../scripts/codex_wave_event_komega/README.md).

`codex_cnrm_hadgem_historical_heff_wk.ipynb` compares just the historical
CNRM-CM6-1 and HadGEM3-GC31-LL runs: five-wave 2-sigma depth distributions
in linear/logarithmic views, plus symmetric/antisymmetric WK spectra on a shared
scale. Model members, native calendars and the fit selection are explicit.

`codex_cnrm_hadgem_noaa_olr_regression.ipynb` compares NOAA and the two historical
runs in five-wave geographical panels. It maps filtered OLR variance, marks
each dataset/wave's independently selected maximum, and shows zero-lag
regressions of unfiltered OLR onto those local filtered indices. Scales are
shared across datasets within each wave. Reproduction instructions are in
[`scripts/codex_olr_regression/README.md`](../scripts/codex_olr_regression/README.md).
