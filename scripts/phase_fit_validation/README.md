# Carrier-fit validation

This folder compares two estimators of carrier zonal wavenumber and frequency
on the **same active pixels** of detected 1.5-sigma events:

- `raw`: the existing constant-amplitude, offset-plus-sinusoid component fit in
  `src.wave_fitting.fit_event_sinusoid`.
- `phase`: a fit to the unit longitudinal analytic signal, with a single phase
  plane and free constant phase. The sign of zonal wavenumber determines which
  temporal phase direction represents eastward or westward propagation.

The known carrier follows the relevant equatorial dispersion curve. Six
positive amplitude envelopes are tested: smooth, skewed, bursty, two-lobed,
and two observed component envelopes ranked by temporal roughness. Empirical
envelopes come from the completed NOAA/ACCESS-CM2 calibration's saved filtered
field and are transplanted onto a new known carrier. Every packet passes through
the same preprocessing/filter before being added to an independent WK99-style
Gaussian background. No code reads archive paths outside `data/`.

The 15 runs cover NOAA OLR, ACCESS-CM2 historical, and ACCESS-CM2 AMIP for all
five waves. Each run tests three in-box carrier/depth pairs, two target
noiseless peak strengths (1.05 and 1.5 times that source/wave's 1% cutoff),
and six independent backgrounds. A clean filtered-packet case is included for
each shape/strength. Null-record maxima are fit separately as diagnostic
controls, not as a new significance calibration. The 1% cutoff and detector
are unchanged. Component matching uses the known packet's filtered footprint;
unmatched and nonsignificant trials remain in the CSV.

From the project root on Levante, first run a small pilot:

```bash
sbatch --account=bb1153 --partition=interactive --cpus-per-task=8 --mem=14G \
  --time=00:30:00 --chdir="$PWD" \
  --output=scripts/phase_fit_validation/logs/pilot.%j.log \
  --wrap="env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  \$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python \
  scripts/phase_fit_validation/run.py --source NOAA_OLR --experiment observations \
  --wave n1ERW --pilot"
```

After checking the pilot CSV, submit all cases with
`bash scripts/phase_fit_validation/submit.sh`. Then run
`python scripts/phase_fit_validation/summarize.py` in the same environment.
Outputs `trials.<source>.<experiment>.<wave>.csv` and corresponding `shapes.*.csv`
stay in this folder, as do `summary_by_shape.csv`, `summary_by_wave.csv`,
`summary_paired_by_wave.csv`, `paired_errors.csv`, `clean_checks.csv`, and
`null_controls.csv`, and `provenance.csv`. Additional summaries by morphology,
strength and clean-packet accuracy, plus a paired source/background-block
bootstrap, are also saved here. Bootstrap intervals describe these three
sources and their reused backgrounds; they do not generalize across CMIP6.
The case-specific summaries show where the known carrier is distorted by the
wave filter before background noise is added. `raw_r2_screen.csv` checks how
the existing 0.3 R2 rule relates to known carrier error on active pixels.

An additional ER control uses global K=-1 at 90 m depth (roughly 47-day period)
with `--er-k1-control`; its CSVs live in `er_k1_control/`. This distinguishes
low wavenumber from the approximately 89-day carrier at the 90-day box limit.
The main summaries exclude this control. Once all three source runs exist,
`summarize.py` also writes two `er_k1_control_*.csv` files and a compact
`er_k1_period_comparison.csv` contrasting the 47- and 89-day cases.

The Historical WIG run was repeated after a two-event indexing issue made its
two empirical envelope labels use the same observed component. The final CSVs
use distinct observed templates. The preliminary duplicate-template CSVs are
retained in `preliminary_duplicate_template/`; an earlier hard-mask pilot is in
`preliminary_hardmask/`. Neither preliminary folder enters the summaries.

## Completed result — 2026-11-01

Slurm array `28617664` completed all 15 runs; Historical WIG was replaced by
successful targeted rerun `28618239`. ER K=-1/47-day controls `28619451`–
`28619453` also completed. The main set has 3,240 paired background trials,
270 paired clean trials and 90 paired null-record maxima. Background draws are
reused across packet cases, so trial rows are correlated. All 15 final trial
CSVs have expected row counts and no recorded fit exceptions.

For detected, significant, successfully fitted packets, median absolute global
K error / median relative frequency error (raw versus phase) are:

| Wave | Raw K | Phase K | Raw frequency | Phase frequency |
| --- | ---: | ---: | ---: | ---: |
| Kelvin | 0.124 | 0.154 | 2.0% | 2.3% |
| MRG | 0.064 | 0.072 | 0.9% | 1.2% |
| n0EIG | 0.071 | 0.090 | 1.0% | 1.4% |
| n1ERW | 0.045 | 0.048 | 4.5% | 5.7% |
| n1WIG | 0.067 | 0.085 | 1.0% | 1.4% |

The raw fit also has smaller 90th-percentile errors for all five waves in
`summary_by_wave.csv`. Source/background-block bootstrap mean absolute-K error
differences (phase minus raw) are positive for all five waves in
`paired_block_bootstrap.csv`; frequency differences are positive except that
the ER interval just includes zero. Phase-only fitting is often slightly closer
on clean packets, but normalization gives weak active pixels the same influence
as strong ones after background is added. The raw active-component fit is the
better supported default for K and frequency in these cases. This comparison
uses NOAA and one CMIP6 model in two experiments; it does not establish the
ranking for every CMIP6 model.

The ER K=-1 carrier at 89.2 days lies near the 90-day box edge. Its frequency
is already shifted about 9–10% in clean controlled-shape fits, and the observed
envelope shapes cause much larger clean errors for both methods. Moving the
same K to 47.2 days removes the controlled-shape error, while some observed
irregular envelopes still shift K. See `er_k1_period_comparison.csv`. Neither
tested estimator consistently recovers the input carrier in those ambiguous
cases. On the 90 null maxima, 85 raw component fits pass
R2 >= 0.3 although none pass the calibrated peak-significance cutoff;
component R2 is descriptive fit quality, not an independent wave test.

No NOAA or CMIP6 event catalog was recomputed by this validation.

The primary comparisons are absolute global-K error and relative frequency
error, including 90th-percentile tails, on the same significant detected
trials. Report detection/significance denominators and errors among all matched
and significant events separately. A high fit score on a selected envelope
component is not an independent detection test. An envelope that varies as
fast as its carrier, or contains overlapping carriers, may no longer have a
well-defined single filtered carrier; large errors in such cases require
inspection before attributing them solely to the estimator.
