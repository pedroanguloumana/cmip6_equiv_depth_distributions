# Agent-authored CMIP6 depth and power diagnostics

This workflow extends the NOAA event-variance diagnostic to the 59 completed
1990–2010 CMIP6 simulations (34 historical, 25 AMIP). It reads the existing
1.5σ active-component fitted catalogs and the matching WK99-style calibration
spectra. It does not refit events or rerun Gaussian calibration.

From the project root on Levante, with this checkout installed:

```bash
bash scripts/codex_cmip_heff/submit.sh pilot
# Audit task 0 and its two outputs, then:
bash scripts/codex_cmip_heff/submit.sh full
# After all 59 tasks finish:
python scripts/codex_cmip_heff/summarize.py
```

Each Slurm task replays one model's native-calendar `rlut` through the catalog
filter and event detector, verifies event IDs and peaks against all five fitted
catalogs, and records active-footprint and fitted-sinusoid energies by fitted
depth. It also maps saved observed and smooth-background spectral power to
apparent depth using the same rectangular wave-filter boxes. Both diagnostics
use 1 m bins with edges 1–200 m and retain totals outside this plot range.

The output pair for each model/experiment/member is under
`data/wave_significance/codex_cmip_heff_comparison/`:

- `*.event_energy.csv.gz`: 10 rows (five waves × quality-screened/relaxed),
  with the same energy definitions as the NOAA worker in
  `scripts/codex_event_variance/run.py`. Divide bin or selected energy by
  `denominator_energy` for the fraction of the model's total preprocessed,
  projected OLR anomaly energy.
- `*.spectral_power.csv.gz`: 10 rows (five waves × model/smooth background),
  with bin, full-box, valid-depth and 12–50 m power. Divide by
  `denominator_power` for the fraction of total projected spectral power.

These are descriptive diagnostics. Event-footprint energy includes background
energy on active pixels; spectral-cell depth is an apparent dispersion depth,
not an event fit. The Gaussian spectrum is estimated from each model's OLR.
No figures are written by this workflow.

The summary audit verifies complete model/wave/selection coverage, energy
accounting, and fitted-catalog counts. It writes three compressed model-wide
tables: `all_model_event_energy.csv.gz`, `all_model_spectral_power.csv.gz`, and
`all_model_event_counts.csv.gz`. The count table has the same 1 m bins plus
50 logarithmic bins from 0.5–10,000 m. `record_days` preserves native-calendar
length so event rates can be compared per year or per 1,000 days.
The same audit also writes `noaa_event_energy.csv.gz`,
`noaa_spectral_power.csv.gz`, and `noaa_event_counts.csv.gz` with the same
binning and identity columns, using the completed NOAA diagnostic and catalogs.

`notebooks/codex_heff_comparison_cmip6.ipynb` reads these six audited tables,
shows individual-model 12–50 m metrics, and plots NOAA against Historical and
AMIP model medians and between-simulation interquartile ranges. It also shows
1 m event-count bins from 1 to 100 m for each wave, with both standard and
relaxed fit-quality selections. The individual-model panels plot all 34
Historical and 25 AMIP simulations as faint curves, their equally weighted
ensemble means as bold curves, and NOAA as black. To render its
outputs on Levante after the audit, run
`bash scripts/codex_cmip_heff/execute_notebook.sh`; this writes the executed
notebook under the diagnostic data directory for review and transfer. It does
not save standalone figure files.
