# Equivalent-depth bands in individual WK segments

This additional analysis asks how excess OLR power is distributed over apparent
equivalent depths of **10--80 m**, separately for Kelvin, n=1 ER, MRG, n=0 EIG
and n=1 WIG, in NOAA and the existing 34 historical / 25 AMIP simulations.
It does not change the event-fitting analysis or its products.

## Calculation

- Daily 1990--2010, inclusive, each source's native calendar; identical NOAA
  1-degree grid. All input access uses the existing `data/` symlinks.
- Symmetric/antisymmetric projection over 15 S--15 N; **average the individual
  latitude powers**, not the latitude-averaged OLR. The current WK notebook sums
  these latitude powers; dividing by the common latitude count gives physical
  mean-variance units without changing its power/background ratios.
- Remove the full-record mean/trend and three annual harmonics. Take exact
  **96-sample** windows every 36 days (60-day overlap), remove each segment's
  mean/trend and apply a 10% Tukey taper. Correct powers by mean(taper squared).
  Fourier power closes against each tapered segment's corrected mean square.
- Use 5 m depth bands: [10,15), ..., [75,80] m. At each integer zonal K, integrate
  between the band's two branch-specific Matsuno curves, clipped to the existing
  rectangular wave box. Assume zero background wind, using `BOX_PROFILES` and
  `xr_ccew.matsuno.profile_curve_frequency` for the full dispersion relations.
  K=0 and the exact daily Nyquist are excluded; propagation direction at the
  latter cannot be identified from daily samples.
- Divide spectral-cell mass fractionally between depth bands by frequency-cell
  overlap at discrete K. This avoids empty/jagged fine bins caused by assigning
  every cell to its center depth. It **assumes uniform density within each
  Fourier frequency cell** and creates no new depth resolution.
- For each dataset, pool the full-record mean symmetric/antisymmetric spectra
  and smooth with ten 1-2-1 passes in each spectral direction. Keep this
  dataset-specific background fixed across segments. This matches the current
  WK plot's simple smoothing settings; it is not the calibration workflow's
  frequency-dependent 10/20/40-pass adaptation.

For segment `s`, depth band `b`, frequency/K cell `i`, and allocation weight `w`:

```text
observed_power[s,b] = 2 sum_i w[b,i] P[s,i]
background_power[b] = 2 sum_i w[b,i] B[i]
signed_excess[s,b] = 2 sum_i w[b,i] (P[s,i] - B[i])
positive_excess[s,b] = 2 sum_i w[b,i] max(P[s,i] - B[i], 0)
positive_ratio_excess[s,b] = sum_i w[b,i] max(P[s,i] / B[i] - 1, 0)
```

The factor two includes the real field's conjugate half-plane. The first four
quantities have units W2 m-4; the last is a dimensionless weighted cell sum and
is not a contribution to OLR variance. `effective_frequency_cells` allows an
area-normalized ratio diagnostic. The output also keeps
`mean_spectrum_positive_excess = 2 sum w max(mean_s P - B, 0)`, so the effect of
rectifying individual noisy segments can be assessed directly.

## Comparisons and limits

Compare both absolute mean excess and normalized shares conditional on a wave's
10--80 m range. These shares are **power-weighted distributions**, not wave-event
counts or independent estimates of one depth per segment. A segment can contain
several depths. The saved largest-mass bin is only a descriptive spectral maximum.

Compare Kelvin and ER independently; equality is a hypothesis, not imposed.
Their spectral coverage and resolving power differ. Near 40 m at |K|=5--15,
one 1/96-cpd frequency bin corresponds to approximately 1--4 m for Kelvin and
16--27 m for n=1 ER. Low-K ER is still less informative. Five-meter output is a
display/integration grid; inspect coarser 10--30, 30--50, 50--80 m groups and the
saved one-frequency-bin depth scales. Those scales are not confidence intervals.
The rectangular boxes also truncate the accessible K range differently by depth;
effective cell counts are saved so that this geometry is visible.

The synthetic check in `test_analysis.py` gives pure Kelvin and n=1 ER carriers
the same true depth (25 m) and |K|=5. Their resolved depth distributions still
have different peaks/widths: in the checked 96-day recipe the maximum-mass bins
are centered at 22.5 m (Kelvin) and 27.5 m (ER). Their binned interquartile
widths are 5 and 10 m. This is a resolution control, not a calibrated error bar
for noisy observations; fine cross-family width/peak comparisons need matched
response or longer-window sensitivity tests.

Positive clipping has a positive noise expectation even when P has mean B.
Inspect the signed excess, rectified full-record mean and background-relative
diagnostics before interpreting a positive-excess shape as a wave preference.
No null calibration or equivalence test is provided here. With linear signed
excess, averaging segments and integrating depth bands commute exactly; keeping
segments adds variability/intermittency information, not extra resolution.
Overlapping segments and adjacent depth bins are dependent. Future confidence
intervals should resample long time blocks and refit the background, rather than
treating segment/bin values as independent. Background wind differences can
shift apparent depths; these resting-background estimates are not intrinsic
vertical-mode depths or proof of dynamical wave identity.

The new recipe follows the intent of [Wheeler and Kiladis
(1999)](https://atmos.washington.edu/~breth/classes/AS552/papers/WheelerKiladis1999.pdf).
It deliberately uses exactly 96 samples, segment detrending and 60-day overlap.
The existing plotting helper currently requests 30-day overlap and retains the
inclusive endpoint from `segment_data` (97 samples); it is left untouched.

## Run on Levante

The project must be installed in the active environment (`src` is its verified
current installed package name). No `sys.path` modification is needed.

```bash
python -m unittest discover -s scripts/codex_segment_heff -p 'test_*.py' -v
bash scripts/codex_segment_heff/submit.sh pilot
# Inspect both pilot outputs and their logs, then:
bash scripts/codex_segment_heff/submit.sh full
python scripts/codex_segment_heff/summarize.py
python scripts/codex_segment_heff/execute_notebook.py
```

`run.py --task-index 0` selects NOAA; 1--59 follow the saved calibration cohort
task list. All five waves run together for each source. Each worker preserves
existing results and writes one compact NetCDF under
`data/wave_significance/codex_segment_heff/`, containing per-segment band powers,
dates, background, geometry, frequency-to-depth scales and source provenance.
The summary worker requires all 60 records, checks identity/recipe/geometry and
power closure, then writes two compressed CSVs and `audit.json`. The companion
`notebooks/codex_segment_heff.ipynb` uses those summaries and NOAA's NetCDF.
Figures are displayed in the notebook, with no separate figure files saved.
The Historical/AMIP two-panel figure overlays all five wave families, averaging
each model's normalized positive-excess distribution with equal model weights.
Dashed curves in matching wave colors show the normalized NOAA distribution
in both panels; solid curves show the model means.
`finalize.sh` runs the audit and notebook in a small Slurm job; submit it with
`--dependency=afterok:<array-job-id>` after the full array. The executed notebook
is saved alongside the compact data products and can be copied back for review.
If rendering fails after a successful audit, rerun `execute_notebook.py` alone;
the completed summary tables are preserved.
