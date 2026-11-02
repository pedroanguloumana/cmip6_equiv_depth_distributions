# Agent-authored Kelvin event spectral-location diagnostic

The completed 1.5σ active-component fitted catalogs contain event-scale
global zonal wavenumber `K_global` (cycles per circumference) and frequency
`frequency_cpd` (cycles per day). The dimensional coordinates are
`k_rad_m = K_global / 6371008` and
`omega_rad_s = 2π frequency_cpd / 86400`.

On Levante, from the project root with this checkout installed, run:

```bash
python scripts/codex_kelvin_komega/extract.py
```

The script reads only the full 1990–2010 Kelvin catalogs through `data/`,
excludes pilot files, validates their sidecars, and writes
`data/wave_significance/codex_kelvin_komega/kelvin_event_fits.csv.gz` plus
`kelvin_catalog_index.csv`. It preserves every detected candidate and its fit
quality flag; the notebook makes the plotted selection explicit. These
compact CSVs are derived inputs for
`notebooks/codex_kelvin_event_komega.ipynb`. This coordinate extraction reads
no raw OLR fields, and the notebook saves no standalone figures. The executed
notebook contains 60 single-panel
joint histograms for all successful fits and another 60 for the established
fit-quality subset, using bins centered on integer global wavenumbers 1–15
and 0.02-cycle/day frequency bins.
Each map shows the percentage of that record's selected events per bin, with
the same color scale across all records within its selection.

To assign **OLR variance** to fitted event coordinates and compare it with
Fourier variance, run a NOAA/model pilot followed by the remaining records:

```bash
bash scripts/codex_kelvin_komega/submit_variance.sh pilot
# Audit the two pilot outputs and logs, then:
bash scripts/codex_kelvin_komega/submit_variance.sh full
python scripts/codex_kelvin_komega/extract_wk.py
```

`variance.py` matches event IDs and peaks against saved Kelvin filtered fields,
uses the original OLR filter if the saved field changes IDs, and checks its
energy and full-record Fourier bins against the previous CMIP/NOAA audits.
It writes one compact `variance_bins_integer_centered/*.csv.gz` per record.
`extract_wk.py`
bins the saved 96-day symmetric segment spectra and smooth pooled background
into `wk_segment_bins_integer_centered.csv.gz`. Both scripts read derived inputs through
`data/`. `notebooks/codex_kelvin_variance_komega.ipynb` plots, for each record,
quality-screened event-footprint variance, raw segment OLR variance, and
WK-style background-normalized spectral excess as separate single-panel maps.
The event-footprint variance includes background on active pixels; the
segment spectrum uses all OLR, independent of event detection.
