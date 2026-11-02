# Scripts

Routine catalog exploration needs no batch runner. The earlier one-off
audit runner is preserved in `archived/catalog-audit-2026-10-27/`.
See [the workflow guide](../docs/catalog_workflow.md).

## CMIP6 rlut on the NOAA OLR grid, 1990–2010

This workflow replaces the archived project's manually generated CDO file lists
with `cmip_catalog` selection. It processes **historical and AMIP independently**:
all models with complete 1990–2010 daily rlut coverage in the configured catalog,
one member per model/experiment. Members are tried in numeric `(r,i,p,f)` order;
the first complete member is selected. Tables, source grids and versions remain
explicit; different members or versions are never combined to fill missing years.
Source grids are tried with native `gn` first, then other labels in ascending
natural order, and the newest feasible version within each grid. This retains
the earlier workflow's GFDL-CM4 `gr1` preference over its coarser `gr2` product.
Equal historical/AMIP member labels do not establish paired initial conditions.

`src/regridding.py` retains the earlier workflow's **first-order conservative
CDO remapping**, with fractional-area normalization and float32/NetCDF4 zip-4
output. The destination is the **complete global NOAA CDR OLR grid**, copied
directly from `data/noaa_olr`: 360 longitudes at 0.5°–359.5°E and 180 latitudes
at −89.5°–89.5°N, including exact coordinate labels and cell bounds. No spatial
averaging or wave filtering is applied here. rlut remains daily-mean outgoing
longwave flux in W m⁻², with its original sign and native calendar. Thus a
360-day simulation ends on December 30, 2010.

Required symlinks are `data/cmip6` (read-only source archive), `data/noaa_olr`
(NOAA reference), and `data/cmip6_rlut_noaa1deg` (new output archive). The current
Levante output link points to `/scratch/b/b382284/cmip6_rlut_noaa1deg`, as requested
while `/work` is full. **Scratch is temporary storage**; copy completed products
to permanent storage and update this symlink when space is available. The older
IMERG-grid files and existing event catalogs are preserved.

Install the project and `cmip-catalog-levante[netcdf]` into the chosen environment;
CDO must be available on the compute nodes. For example, on Levante:

```bash
export CDO_BIN=/sw/spack-levante/cdo-2.5.3-74wnuh/bin/cdo
python scripts/prepare_cmip6_rlut_regridding.py \
    --start-year 1990 --end-year 2010 \
    --manifest data/cmip6_regrid_jobs/noaa1deg_1990-2010_v2/manifest.json
bash scripts/submit_cmip6_rlut_regridding.sh \
    data/cmip6_regrid_jobs/noaa1deg_1990-2010_v2/manifest.json --dry-run
bash scripts/submit_cmip6_rlut_regridding.sh \
    data/cmip6_regrid_jobs/noaa1deg_1990-2010_v2/manifest.json
```

The planner only reads catalog rows, attributes and coordinates. It inspects
first/last representative headers within each exact dataset, then uses the
package's full-window matching. It saves a JSON task manifest, inspected metadata,
inspection issues and match diagnostics together. Use a new manifest path for a
new selection. `--models` and `--experiments` explicitly narrow the cohort.
`--metadata` can reuse the planner's saved `inspected_metadata.parquet`; its scope
is the members already inspected in that run. This does not refresh the catalog
or reopen their representative files. The worker still validates every input.

The launcher submits one task per model/experiment/member, at most four at once.
`ACCOUNT`, `PARTITION`, `JOB_MEM`, `JOB_CPUS`, `JOB_TIME`, `MAX_CONCURRENT`,
`PYTHON_BIN` and `CDO_BIN` are explicit machine/resource settings. `ARRAY_SPEC`
can select a subset of manifest indices for a pilot or a failed-task restart.
For example, `ARRAY_SPEC=0,25` submits only those two entries. Logs are under
`data/cmip6_regrid_jobs/logs/`.

Each worker checks every selected source file's identity, grid, unchanged
size/mtime, flux metadata and actual time coordinate. Gaps fail explicitly.
Overlapping records (present in CESM2-WACCM-FV2) are compared in bounded blocks:
only exactly identical decoded OLR records are deduplicated, keeping the first
occurrence and recording the verified day counts in the sidecar. Conflicting
values fail explicitly. It stages and remaps one native-calendar year at a time, computes
conservative weights once, then joins the remapped years. Classic-NetCDF staging
avoids simultaneous HDF5 reads and writes, which caused failures in the earlier
workflow. Temporary files are cleaned up on normal completion or Python errors;
a killed job can leave a directory marked `.partial.*`.

Completed files are written atomically under the output symlink:

```text
<historical|amip>/rlut.day.<model>.<experiment>.<member>.<source-grid>.<version>.noaa_1deg.1990-2010.nc
```

Each NetCDF retains source identity and adds `source_grid_label`, `grid_label=gr`,
`analysis_grid=noaa_1deg` and a provenance fingerprint. Its JSON sidecar records
selection evidence, input files, CDO commands/version and timing. Completion
requires exact dates, NOAA coordinates/bounds, W m⁻² units and finite cells
throughout. A matching completed product is skipped; conflicting provenance
raises an error rather than overwriting an existing result.

Run the calendar checks and synthetic conservative-remapping integration check:

```bash
CDO_BIN=cdo python -m unittest discover -s tests -p test_regridding.py -v
```

The integration test checks the global area mean, exact NOAA labels/bounds,
timestamps, units, and completed-output restart behavior in Gregorian, no-leap
and 360-day calendars, using overlapping input chunks. It also verifies rejection
of conflicting duplicate values. Without CDO it is skipped; coverage checks still run.

Levante validation on 2026-10-30 passed these checks, four real-file 1990 pilots
(ACCESS-CM2, CanESM5, HadGEM3-GC31-LL and unstructured ICON-ESM-LR), and an
independent audit of the source choices. All 90 overlapping CESM2-WACCM-FV2 days
(40 AMIP, 50 historical) contained exactly identical OLR. The corrected manifest
has 34 historical and 25 AMIP tasks; only the two GFDL-CM4 grid choices differ
from the initial selection. Array `28564656` was submitted with this manifest;
it reuses the initial batch's six completed products. Audit and pilot reports
are under `data/cmip6_regrid_jobs/validation/`.

## NOAA OLR event catalogs

`noaa_olr_events.py` processes **one wave per process**. It follows the current
`NOAA_OLR_heff.ipynb` settings: 1990–2010 daily NOAA CDR OLR, seasonal harmonics
and trend removed, Tukey taper, symmetry projection over 20°S–20°N followed by
an unweighted 0–20°N mean, Hilbert-envelope candidates, local sinusoid fits,
and resting-background depth inversion for Kelvin, n=1 ER, MRG, n=0 EIG, and n=1 WIG.
Profiles, fit grids/windows, detection thresholds, and quality screens are
explicit in the shared `src/wave_events.py` module. No plots or catalogs of diagnostic arrays
are accumulated. It still needs memory for the **full record FFT**; this is
not a streaming spectral filter.

Requirements: install this project into the existing
`cmip6_equiv_depth_distributions` environment, with NumPy, SciPy, xarray,
netCDF4 and `xr-ccew`. Configure the local symlink `data/noaa_olr` to resolve to
the NOAA archive directory. The input filename is centralized in `src/paths.py`.
The scripts never write to that archive.

### Review and launch later

These commands are examples; **no analysis or jobs were run when these scripts
were added**. From the project root on Levante, preview the submission:

```bash
bash scripts/submit_noaa_olr_events.sh --dry-run
```

For a small first batch, choose one wave, a two-year record, and three events:

```bash
WAVES="Kelvin" START_DATE=1992-01-01 END_DATE=1993-12-31 MAX_EVENTS=3 \
    bash scripts/submit_noaa_olr_events.sh
```

Then submit the five independent array tasks when ready:

```bash
bash scripts/submit_noaa_olr_events.sh
```

The launcher mirrors the Python-worker/Bash-array structure of
`/work/bb1153/b382284/dyamond/code/feature_workflow`. Each array task handles one
wave; up to four may run concurrently, subject to queue availability. Defaults
are account `bb1153`, partition `interactive`, 14 GiB, 8 allocated CPUs and one hour
per task, based on the pilot below with headroom. The CPU request accommodates
the partition's memory-per-CPU limit; numerical libraries use one thread.
Override `ACCOUNT`, `PARTITION`, `JOB_MEM`, `JOB_CPUS`,
`JOB_TIME`, or `MAX_CONCURRENT` explicitly as needed. If changing partitions,
also check their walltime/memory limits in the
[DKRZ documentation](https://docs.dkrz.de/doc/levante/running-jobs/partitions-and-limits.html).
For example, 14 GiB requires at least 16 CPUs on `shared`.

`PYTHON_BIN` defaults to
`$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python`. It is used inside
each job, so activation in the submitting shell is unnecessary. Set it explicitly
if the environment lives elsewhere. Job logs and task lists go under
`DATA_DIR/noaa_olr_event_jobs/`. Do not launch the same wave/configuration twice
concurrently: its tasks would share checkpoints.

The worker can also run directly **inside an appropriate allocation**:

```bash
python scripts/noaa_olr_events.py --wave MRG --start 1990-01-01 --end 2010-12-31
```

### Files and restart behavior

Completed CSVs are written directly under `DATA_DIR`, with names of this form:

```text
NOAA_OLR.Kelvin.19900101-20101231.k1to15_dk0p5.f0p0714286to0p4_nf17.ds1p5.window.tw10_lw90.ts5_ls45.<run-id>.csv
```

Here `ds` is the envelope activity cutoff in reference standard deviations,
`tw/lw` are the time/longitude half-window, and `ts/ls` are their Gaussian widths.
Component-supported fits use `.component.` in place of the window settings.
The digest covers all scientific settings, code fingerprints, input size/mtime,
and recorded software versions. A matching `.json` stores those settings in full,
run status, reference sigma, counts, phase timings, and peak worker RSS. Individual
rows include `fit_elapsed_seconds`; the JSON summarizes the fit-time distribution
and counts full-record candidates with complete fitting windows. These support
resource estimates from a pilot. Runtime extrapolation remains approximate because
events can require different numbers of optimizer iterations. `--max-events` produces a distinct filename
marked `pilotN`, so a pilot cannot masquerade as a full catalog.

For a resource pilot, retain the **full intended time range** and limit only the
number of fitted events. A shorter input record would underestimate the FFT's
memory requirement. `peak_rss_gib` is the process high-water mark from `getrusage`,
including preprocessing peaks; compare it with Slurm accounting when sizing the
job, since process RSS does not include all job/cgroup overhead or file cache.

### Measured resource pilot

Array `28552504` completed all four pilots on 2026-10-29 (Levante local date),
using the full 1990--2010 record and the first 25 candidates of each wave. All
100 fits completed without numerical errors. Preprocessing took approximately
47 seconds per wave. Worker peak RSS is higher than Slurm's sampled MaxRSS,
so use the worker high-water mark for sizing rather than the lower samples.

| Wave | Peak worker RSS (GiB) | Mean fit (s/event) | Full candidates | Estimated full runtime (min) |
| --- | ---: | ---: | ---: | ---: |
| Kelvin | 10.10 | 0.168 | 2704 | 9--11 |
| n1ERW | 10.10 | 1.559 | 1034 | 28--31 |
| n0EIG | 10.10 | 0.191 | 3287 | 12--14 |
| MRG | 10.19 | 0.162 | 1945 | 6--8 |

The runtime ranges extrapolate the pilot's mean and 90th-percentile fit times,
including observed checkpoint overhead. They are planning estimates, not
confidence intervals or guaranteed bounds; the pilot samples only the first
25 candidates and does not establish population-wide fit quality. The 14 GiB
request leaves about 3.8 GiB above the observed peak; one hour leaves runtime
headroom for all waves. The fitting computation uses approximately one CPU.
Detailed measurements and Slurm accounting are saved in
`DATA_DIR/noaa_olr_event_jobs/benchmarks/pilot_28552504.json`.

Every selected candidate gets a row, including poor fits, invalid dispersion
solutions, and candidates without a complete fitting window. Use
`keep_for_histogram` for the existing quality screens plus a valid positive depth.
`fit_status` and `fit_error` identify numerical failures and incomplete windows.
No fitted depth is clipped to 1–200 m. A record with no detected candidates
produces an empty CSV with the normal header.

Each completed event is first written atomically into `<catalog>.csv.parts/`.
A rerun with identical inputs, code, settings, and recorded versions reloads the
filtered field but skips existing completed event parts. The final CSV is assembled
atomically only after all candidate rows are available. Checkpoints are retained;
the final CSV's existence is the completion marker. An interrupted run therefore
keeps its completed fits, and a completed run is skipped on rerun. Use
`--overwrite` for the worker or `OVERWRITE=1` for the launcher to recompute that
exact catalog and replace its generated checkpoints. Code or input changes produce
a new run ID instead of mixing old and new results.

To load a completed catalog later:

```python
import pandas as pd

events = pd.read_csv(catalog_path, parse_dates=["peak_time"])
accepted = events.loc[events["keep_for_histogram"]]
```

## CMIP6 historical and AMIP event catalogs

The active inputs are the completed **1990–2010 NOAA-grid products** exposed by
`data/cmip6_rlut_noaa1deg`, with 34 historical and 25 AMIP simulations. The older
`data/archived/data_links/cmip6_rlut` link exposes the 2000–2010 IMERG-grid
archive for historical reference; it is not used by the event worker.

`src/cmip.py` provides:

- `list_cmip_simulations(experiments=("historical", "amip"))`: filename identities
  and paths, one per model/experiment; no climate-field reads.
- `list_cmip_models(experiment="historical")`: models in one experiment.
- `inspect_cmip_rlut(model, experiment="historical")`: verify filename/header
  identity, original grid, target grid, flux metadata and native calendar.
- `load_cmip_rlut(model, experiment="historical", year_range=(1990, 2010),
  lat_range=(-20, 20), member_id=None)`: load daily `(time, lat, lon)` rlut in
  W m-2, retaining CF dates, units and simulation identity.

The loader checks complete daily calendar years and the regular 1-degree grid,
normalizes longitude to 0–360, and rejects ambiguous multiple files for one
model/experiment. Native Gregorian/standard, proleptic Gregorian, no-leap and
360-day calendars are retained without interpolation. A 360-day record ends on
2010-12-30. `source_grid_label` identifies the original CMIP grid;
`analysis_grid=noaa_1deg` identifies the actual input grid.

The NOAA and CMIP workers share all preprocessing, filtering, Hilbert-envelope
candidate detection, local sinusoid fitting, resting-background depth inversion
and quality screens in `src/wave_events.py` and `src/wave_fitting.py`. The event
CSV schema is unchanged, including peak longitude/time, native calendar, source,
experiment, member, wave, fitted parameters, `h_eff_m`, and quality flags.
No latitude localization or additional latitude columns are introduced.

### Inspect and submit

From the installed project root on Levante:

```bash
python scripts/cmip6_rlut_events.py --list-simulations
python scripts/cmip6_rlut_events.py --list-models --experiment amip
bash scripts/submit_cmip6_rlut_events.sh --dry-run
```

A pilot keeps the full intended record, limiting only the fitted event count:

```bash
MODELS="CanESM5" EXPERIMENTS="amip" WAVES="n1ERW" MAX_EVENTS=2 \
    bash scripts/submit_cmip6_rlut_events.sh
```

The full batch is:

```bash
bash scripts/submit_cmip6_rlut_events.sh
```

This submits **295 simulation/wave tasks** (59 available model/experiment/member
pairs × five waves), with at most four running concurrently. It uses available
pairs rather than a Cartesian product that would invent unavailable experiments.
The TSV task list explicitly records model, experiment, member and wave; every
worker checks that member against the NetCDF header.

When the Slurm association limits the number of simultaneous array tasks,
`run_cmip6_event_batch.py` can run one experiment's task-list entries concurrently
inside a sufficiently large single-node allocation. For example, from that
allocation and the project root:

```bash
python scripts/run_cmip6_event_batch.py \
    --tasklist data/cmip6_event_jobs/tasklists/<task-list>.tsv \
    --experiment historical --max-parallel 16 \
    --detection-sigma 1.5 --fit-support window
```

Run `--experiment amip` separately. Each subprocess uses the standard CMIP
worker, so catalog identity, checkpoints and completed-output skipping are
unchanged. Per-task output is under
`data/cmip6_event_jobs/logs/CMIP6_pooled.<job-id>/`; the batch exits nonzero if
any task fails. Size the allocation for all concurrent workers (the measured
full-record pilot peak was about 10.2 GiB per worker).

Defaults are 1990–2010, `bb1153`/`interactive`, 14 GiB, eight allocated CPUs and two
hours per task. Full 1990–2010 NOAA pilots on the same selected grid peaked at
10.2 GiB. Numerical-library threads remain one; the CPU allocation supplies the
partition's memory allowance. Override `MODELS`, `EXPERIMENTS`, `WAVES`,
`START_YEAR`, `END_YEAR`, `MAX_EVENTS`, `ACCOUNT`, `PARTITION`, `JOB_MEM`,
`JOB_CPUS`, `JOB_TIME`, `MAX_CONCURRENT`, `PYTHON_BIN` or `OVERWRITE` explicitly.
The requested years must actually exist in each source record.

One worker can also run directly inside an appropriate allocation:

```bash
python scripts/cmip6_rlut_events.py --model CanESM5 --experiment amip \
    --member r1i1p1f1 --wave Kelvin --start-year 1990 --end-year 2010
```

### Outputs, provenance and archives

Completed catalogs are saved under `data/cmip6_event_catalogs/`, with names such as:

```text
CMIP6.CanESM5.amip.r1i1p1f1.rlut.Kelvin.19900101-20101231.k1to15_dk0p5.f0p0714286to0p4_nf17.ds1p5.window.tw10_lw90.ts5_ls45.<run-id>.csv
```

`tw/lw` are the time/longitude half-window, and `ts/ls` are the Gaussian weighting
scales. JSON sidecars store the full scientific settings, source identity,
calendar, NetCDF and regridding-provenance paths/size/mtime, code hashes, versions,
quality counts and resource measurements. Filenames include experiment and member
so different simulations cannot collide. Atomic per-event checkpoints,
completed-catalog skipping, and resume behavior follow the NOAA workflow above.
Pilot catalogs have `.pilotN.csv` names and are distinct from full populations.
Use `keep_for_histogram` to select the screened event distribution.

Read CMIP peak dates as strings together with the `calendar` column; Gregorian
parsing cannot represent February 30. Additional environmental fields must be
matched to the same model, experiment, member and native-calendar dates.

Logs, task lists and current validation reports are under
`data/cmip6_event_jobs/`. Superseded 2000–2010 CMIP event products, interrupted
checkpoints and pilot products are moved to the project's existing
`archived/event_catalogs/<timestamp>/`, exposed through `data/archived`.
Each move has a manifest preserving original paths. The 1990–2010 curved-filter
NOAA and CMIP catalogs have also been archived pending the revised box-profile
calculation. Source NetCDF archives and completed regridding products remain active.

Run the calendar, known-wave-depth and experiment-aware loader checks with:

```bash
python -m unittest discover -s tests -p test_cmip_workflow.py -v
```

### Full-record validation and resource check

Preflight job `28565802` passed six tests, including 16 known-wave/depth fits
across four calendars and loader fixtures spanning both experiments and four
calendars. All 59 completed regridding products were checked against their
simulation identities and regridding provenance. The available-pair task list
contained exactly 236 entries in the full dry run.

Pilot array `28565811` processed the complete 1990–2010 record, fitting five
candidates for each of four cases: historical ACCESS-CM2/Kelvin, AMIP
CanESM5/n1ERW, historical HadGEM3-GC31-LL/MRG and AMIP CNRM-CM6-1/n0EIG.
Peak worker RSS was 10.09–10.24 GiB. All 20 rows were saved without processing
errors; one ER optimization did not converge and was flagged/excluded by the
existing quality screen. Five-event timing extrapolations were about 9, 47, 5
and 13 minutes respectively. These small-sample estimates do not bound every
model or event; the full jobs request two hours and 14 GiB for headroom.
Detailed metrics are in
`data/cmip6_event_jobs/benchmarks/noaa1deg_1990-2010_pilot_28565811.json`.
Pilot CSVs/checkpoints are archived separately from the full populations.

### Loading NOAA and CMIP6 event catalogs

Use the same loader for observations and models:

```python
from src.wave_events import load_event_catalog

sources = [("NOAA_OLR", "observations"), ("ACCESS-CM2", "historical")]
waves = ["Kelvin", "n1ERW", "MRG", "n0EIG"]
event_catalogs = {
    (model, experiment): {
        wave: load_event_catalog(model, experiment, wave, accepted_only=True)
        for wave in waves
    }
    for model, experiment in sources
}
```

CMIP6 experiments can be `historical` or `amip`; NOAA uses `NOAA_OLR` (alias
`NOAA`) and `observations`. `load_noaa_event_catalog(wave_type, ...)` and
`load_cmip6_event_catalog(model, experiment, wave_type, ...)` are also available.
The functions return pandas DataFrames with common simulation/calendar labels,
including labels supplied when loading the older NOAA catalogs.

The default `accepted_only=False` retains every candidate. Set it to `True` to
apply the catalog's existing `keep_for_histogram` flag. `year_range=(1990, 2010)`
selects the complete catalog record, rather than clipping individual peak dates;
`None` leaves the record period unconstrained. Only completed full catalogs in
the active data directories are considered. Archived files, pilots and
checkpoints are excluded. Missing catalogs raise `FileNotFoundError`; multiple
matches require explicit `member_id`, `year_range` or `run_id` selection.

`peak_time` remains a string alongside the `calendar` column so native dates such
as February 30 are preserved. The selected CSV path and JSON provenance are
available as `frame.attrs["catalog_path"]` and
`frame.attrs["catalog_metadata"]`. Loading catalogs does not require raw OLR.

Run the catalog-loader checks with:

```bash
python -m unittest discover -s tests -p test_event_catalog_loaders.py -v
```

### Candidate significance calibration (experimental)

`calibrate_wave_significance.py` generates stationary correlated Gaussian fields
from a smooth 2-D spectrum of each source's symmetric/antisymmetric hemispheric
means. It uses the current catalog preprocessing, masks and connected-component
detector. Averaging before these linear operations is verified against the full
latitude workflow. This avoids manufacturing independent latitude noise that
would disappear incorrectly under meridional averaging.

From an allocation with the project installed and the data symlinks configured:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/calibrate_wave_significance.py \
    --model NOAA_OLR --experiment observations --realizations 20

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/calibrate_wave_significance.py \
    --model ACCESS-CM2 --experiment historical --realizations 20
```

Defaults are the full 1990–2010 native-calendar record, all five waves, seed
20261030, and the 99th percentile of **detected null-candidate peak envelopes**.
The first half of the independent records calibrate that quantile; the remaining
half estimate false candidates per record. Four realizations are a smoke test
only. Twenty demonstrate the workflow but do not establish stable 1% whole-record
significance: with ten calibration records, record-maximum Monte Carlo p-values
have minimum resolution 1/11. More realizations and uncertainty assessment are
needed before adopting a significance policy.

The null is not independent Gaussian pixels, nor a phase scramble preserving
sharp observed wave peaks. FFT power is smoothed with Gaussian kernel standard
deviations `--frequency-smoothing-cpd 0.0104166667` and
`--wavenumber-smoothing 2`. These are visible provisional choices: this is
WK-inspired smoothing, not the original WK background recipe. The draws match
the resulting covariance and variance in expectation, not the exact variance of
each individual record. They are zero-mean anomalies, not absolute OLR fluxes.

Each run writes a hashed identity/settings prefix under `data/wave_significance/`:
the background spectrum (`.background.nc`), null events (`.null_events.csv`),
per-realization counts (`.realizations.csv`), observed candidates with diagnostic
tail statistics (`.observed_events.csv`), and provenance/thresholds (`.json`).
Original catalogs, detector thresholds and fit flags are untouched. A pooled
99th-percentile peak threshold is **not** 1% catalog FDR or a 99% guarantee of
no false event in the record. The report gives held-out false counts explicitly.
Its pooled candidate-tail fraction is descriptive, not an exact p-value for
dependent, selected events. Record-maximum p-values address a different,
whole-record comparison and expose their simulation-count resolution.

Remaining validation includes smoothing sensitivity, non-Gaussian tails,
geographic/seasonal variance, and injection recovery across wave parameters.
This prototype calibrates amplitude significance, not precision of k/omega or
the physical identification of a wave family.

```bash
python -m unittest discover -s tests -p test_wave_significance.py -v
```

### Validating the background assumptions

`validate_wave_significance.py` starts the methodological checks on NOAA
1990–2010. Run `--symmetry symmetric` and `--symmetry antisymmetric` as separate
allocation tasks, using the same unique `--run-label`. Outputs are isolated in
`data/wave_significance/validation/<run-label>/<symmetry>/`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/validate_wave_significance.py \
    --symmetry symmetric --run-label noaa_null_comparison --realizations 8 --injection-realizations 2
```

The first-screen comparison includes:

- Signed-anomaly moments and tails before filtering and for each wave, globally
  and by season/longitude sector. These are descriptive diagnostics, without
  treating correlated grid points as independent samples.
- Gaussian backgrounds and 2-D iterative amplitude-adjusted Fourier surrogates
  (IAAFT, adapted from Schreiber and Schmitz, 1996). Both target the same smooth
  spectrum and expected global variance. IAAFT exactly retains the empirical
  projected-anomaly marginal distribution, which also includes real waves.
  Its *approximate* spectral match must reach the recorded tolerance (default
  3% relative Fourier-amplitude error) before its events enter calibration.
  Broad frequency/wavenumber power ratios and filtered variances also record
  remaining discrepancies; global convergence does not guarantee every wave
  band is matched equally well. Failed surrogate matches are reported/excluded.
- Candidate-peak thresholds fitted on the first half of null records and false
  counts evaluated on the other half. Monte Carlo uncertainty resamples whole
  records, not individual pixels/events. Eight records are a pilot, not a final
  significance calibration; these intervals exclude uncertainty in the fitted
  background model.
- Controlled packets at depths 8/25/90 m, two K values per wave, compact/extended
  envelopes and amplitudes 2/4 times the Gaussian background's mean filtered
  reference sigma. Envelopes translate at the dry-branch group velocity; phase
  and group speeds remain distinct. Some carriers intentionally fall outside
  the rectangular filters; retained packet energy quantifies selection loss.
  The same preprocessing/filter is applied to packets before superposition with
  held-out filtered backgrounds, using linearity. The detector recomputes sigma
  on the combined field, as it would for actual data.
- Injection matches require catalog peak locations within one temporal and
  longitudinal envelope sigma along the moving packet. Background-only matches
  are recorded alongside injected matches to expose chance associations. This
  tests *detection* selection by known depth; it does not validate fitted k/omega
  or depth estimates. Small first-screen injection replication counts must be
  expanded before interpreting recovery rates precisely.

Both null families remain stationary and lack the observed regional/seasonal
variance pattern. Their difference tests the marginal constraint, not the full
physics of convective intermittency. Spectral smoothing is held fixed in this
comparison and still needs a separate sensitivity study. The script does not
change production detection thresholds, fitted catalogs, notebooks or figures.

After both tasks complete, produce the combined report with:

```bash
python scripts/summarize_wave_significance_checks.py --run-label noaa_null_comparison
```

This writes `distribution_summary.csv`, `model_comparison.csv`,
`injection_recovery_summary.csv` and `analysis_summary.json` alongside the two
task folders. It checks that the regional/seasonal partition reproduces global
variance, calculates the kurtosis expected from a Gaussian mixture with those
means/variances, and compares absolute and sigma-normalized event thresholds.
Injection summaries distinguish in-band, daily-resolvable carriers from
intentionally excluded or aliased controls. The Gaussian-mixture calculation is
descriptive; the observed distribution can contain both background and waves.

### Full Gaussian record-search calibration

```bash
bash scripts/submit_gaussian_event_calibration.sh --dry-run
bash scripts/submit_gaussian_event_calibration.sh
# For a model, after validating the NOAA method:
MODEL=ACCESS-CM2 EXPERIMENT=historical bash scripts/submit_gaussian_event_calibration.sh
```

The default protocol is NOAA 1990–2010, all five rectangular wave profiles,
1,000 Gaussian calibration records plus 1,000 independent validation records,
and seed 20261031. The score is a **detected packet's maximum Hilbert envelope**,
in W m-2, retaining the existing duration/span/connectivity requirements. It is
not a single signed-OLR pixel test. The primary flag is a **1% complete-record
Gaussian-null false-alarm test per source/experiment/member/wave**, conditional
on the fitted stationary background spectrum. It is not a 1% catalog FDR or a
joint guarantee across every wave/model, nor proof of dynamical wave identity.

For B calibration records, store each record's largest candidate envelope M_b.
For observed score T, `p_record = (1 + count(M_b >= T)) / (B + 1)`. Saved cutoffs
use the strict comparison `T > threshold`, tested to match `p_record <= alpha`
even with ties. Equality does not pass. Tables also include 5%/10% record-level
rules, pooled candidate rarity scores, and threshold curves for other policies.
No fit-quality flags or original event catalogs are overwritten.

The launcher submits dependent preparation, four first-chunk pilots, the full
checkpointed array, reduction, held-out packet injections, and final reporting.
Complete pilot chunks are reused by the main array. Defaults are 100 records
per task, at most four simultaneous simulation tasks, 4 CPUs and 6 GiB per task.
`N_CALIBRATION`, `N_VALIDATION`, `CHUNK_SIZE`, `MAX_CONCURRENT`, `START_YEAR`,
`END_YEAR`, `MODEL`, `EXPERIMENT`, `MEMBER`, `ACCOUNT`, `PARTITION`, and `PYTHON_BIN`
are explicit launcher inputs. It records job IDs in the ensemble's
`jobs/submission.txt`; it will not silently submit the same manifest twice.
Chunks are atomic and reusable after interruption. Simulations stop if their
recipe, relevant functions or numerical environment changes. Analysis-only
changes can reuse saved scores, with separately saved analysis code provenance.

Products are under
`data/wave_significance/gaussian_calibration/<model>.<experiment>.<member>.<years>.<hash>/`:

- `manifest.json`, `code_snapshot/`, `analysis_code_snapshot/`: source identity,
  native dates/calendar, all profiles and preprocessing/detection settings,
  exact seed scheme, package versions, code and configuration hashes.
- `template.<symmetry>.nc`, `spectrum.<symmetry>.nc`: projected input coordinates
  and physical units, the unsmoothed taper-corrected periodogram and fitted
  Gaussian background power. No per-realization variance forcing is applied.
- `chunks/task_*.events.csv.gz`: **every null candidate**, with its peak time,
  longitude, duration, span, score and realization/role. These permit new
  statistics/plots without repeating the simulations.
- `null_records.csv`, `record_max_ecdf.<wave>.csv`: null maxima, counts, reference
  sigma, signed-field moments and timing for every full-record realization.
- `thresholds.csv`, `threshold_convergence.csv`: cutoffs, order-statistic Monte
  Carlo uncertainty intervals, observed accepted counts, held-out probability
  of any false event with binomial intervals, and expected false-candidate counts.
- `threshold_curve.<wave>.csv`, `null_threshold_counts.<wave>.nc`: observed
  counts and complete per-realization null counts across a dense threshold grid
  including every observed score. Mean-count intervals bootstrap whole records;
  when no validation exceedance occurs, no zero upper mean-count bound is claimed.
  The binomial upper bound on record false-alarm probability remains available.
  Estimated false-count/observed-count ratios are diagnostic, not guaranteed FDR.
- `observed_significance.<wave>.csv`: all observed candidates with adjusted
  p-values, explicit significance flags, source labels and ensemble ID. Existing
  event IDs/times/longitudes are retained for subsequent propagation fitting.
- `observed_distributions.csv`, `null_histograms.<wave>.nc`,
  `spectral_bands.<wave>.csv.gz`, `spectral_summary.<wave>.csv`: regional/seasonal
  observed diagnostics, per-realization signed/envelope histograms in reference
  sigma units, and checks of simulated versus target spectral variance.
- `observed_filtered.<wave>.nc`, `example_null.<wave>.nc`: fields for later
  Hovmöller/QA plots, including one reproducible null realization and its envelope.
- `injections.csv.gz`, `injection_recovery.csv`: held-out packet recovery at
  known depths/K, compact/extended widths, and amplitudes 2/4/8 background sigma.
  Eight held-out backgrounds are used per case; phase/group speed, filtering
  losses, aliased/out-of-band controls, and pre-existing matches are recorded.
- `report.json`: final product index, threshold summary and methodological scope.

Intervals condition on the fitted Gaussian null and do not include uncertainty
in that background model. Injection tests measure detection completeness, not
accuracy of fitted k/omega or h_eff. All numerical diagnostics are retained;
the workflow does not create figure files.

```bash
python -m unittest discover -s tests -p test_gaussian_calibration.py -v
```

## WK99-style 1-2-1 null ensemble

Run a separate NOAA calibration with the same filters, detector, Gaussian draw
algorithm, calibration/validation sizes and saved diagnostic tables:

```bash
BACKGROUND_METHOD=wk99_121 bash scripts/submit_gaussian_event_calibration.sh --dry-run
BACKGROUND_METHOD=wk99_121 bash scripts/submit_gaussian_event_calibration.sh
```

Outputs are under `data/wave_significance/gaussian_calibration_wk99_121/`,
with a separate configuration hash. Existing Gaussian-kernel ensembles stay in
`gaussian_calibration/`. Existing completed ensembles retain their saved results
and code snapshots. Recipe guards may reject continuing an older manifest after
workflow code changes; use its snapshot or start a distinct new ensemble.

The new estimator uses `xr_ccew.segment_averaged_spectrum` with 96-day segments
overlapping by 60 days, then pools the symmetric and antisymmetric powers by
their arithmetic mean. It calls `xr_ccew.background_spectrum` for 10, 20 and 40
wavenumber smoothing passes at frequencies below 0.2, from 0.2 to below 0.3,
and at/above 0.3 cycles/day, followed by 10 frequency passes. The same background
is used for both symmetry classes; it is not independently rescaled to each
component's variance. Individual Gaussian draws are not variance-normalized.

This is a **WK99-based adaptation**, not an exact reproduction of the published
OLR analysis. [WK99 sections 2 and 3b](https://atmos.washington.edu/~breth/classes/AS552/papers/WheelerKiladis1999.pdf)
specify 96-day windows, approximately two-month overlap, pooling and 10 frequency
passes, with wavenumber passes increasing from 10 to 40 in two steps. The paper
does not give frequency breakpoints; 0.2 and 0.3 cycles/day are the corresponding
transitions in [NCAR's implementation](https://github.com/NCAR/ncl/blob/develop/ni/src/examples/gsun/diagnostics_cam.ncl).
We retain 10 passes below 0.2, following the paper's stated minimum, rather than
NCAR's additional five-pass rule below 0.1 cycles/day. Other explicit adaptations:

- Use powers of the catalog's latitude-averaged fields, so the spectrum represents
  the actual quantity submitted to the event detector. WK99 sums latitude powers.
- Retain the existing latitude selection, seasonal harmonics and taper; the input
  is daily data, with Nyquist frequency 0.5 cycles/day.
- Include zero frequency in edge-padded 1-2-1 smoothing for synthesis; interpolate
  the positive-frequency background onto the full-record grid, reflecting
  wavenumber for the opposite half-plane and converting WK to NumPy frequency
  signs. Set the single DC coefficient to zero and enforce Hermitian symmetry.
- Normalize the **target spectrum** after regridding to conserve the pooled
  segment variance. The scaling factor and pre/post variance are saved. The
  shared background generally differs in variance from either component alone.

`wk99_background.nc` saves the two unsmoothed segment spectra, their pooled
spectrum, the smoothed positive-frequency background and the final synthesis
power. `background_diagnostics.json` records the segment count, bandwidth recipe,
normalization and adaptations. The usual full-record spectra, raw null candidate
scores, record maxima, threshold/false-count curves, uncertainty, held-out
validation and injection recovery are also saved. Seeds are unchanged between
background methods so a given realization uses the same underlying white noise;
the two methods are paired comparisons, while calibration and validation records
remain disjoint. Differences from the existing ensemble include segment
estimation and symmetry pooling as well as the smoothing kernel, so they should
not be attributed to the kernel alone. Neither recipe removes all wave power.

## Archived outputs and current catalog status

The curved-filter catalogs, old event job records and interrupted background
calibrations were warehoused under `archived/` when the analysis moved toward
rectangular profiles and n=1 WIG. See [the archive index](../archived/README.md)
for locations and per-machine move manifests. Active event catalog
output directories are empty pending recalculation; new calibration runs write
their explicitly recorded profile configurations under `data/wave_significance/`.
Completed regridded inputs and their `data/cmip6_regrid_jobs/` provenance remain active. Archived catalogs
are not searched by `load_event_catalog`; open their recorded CSV paths explicitly
when comparing the previous method. Historical metadata and logs retain original
paths and status strings, which do not indicate currently running jobs.

## NOAA Gaussian calibration for the revised boxes

The full calibration worker explicitly selects `BOX_PROFILES` from
`src/wave_profiles.py`. The event-fitting worker now uses the same collection.
Both observed and null fields, and all injected packets, use this same profile
selection. The complete profile source and values are
included in snapshots, recipe hashes and the manifest (`profile_collection`).

Current boxes, with positive frequency and signed global wavenumber:

| Wave | Global K | Period [days] | OLR symmetry |
| --- | --- | --- | --- |
| Kelvin | 1 to 15 | 2.5 to 14 | symmetric |
| n1ERW | -15 to -1 | 10 to 90 | symmetric |
| MRG | -10 to -0.25 | 2.25 to 10 | antisymmetric |
| n0EIG | 0.25 to 14 | 2.25 to 6 | antisymmetric |
| n1WIG | -15 to -1 | 2 to 4 | symmetric |

The fit search uses the same period bounds, and a half-wavenumber grid over
integer-|K| endpoints. For MRG and n0EIG the filter's ±0.25 bound excludes
the global Fourier K=0 component; the local fit starts at |K|=1. A fit in the
180-degree longitude window needs at least half a spatial cycle, allowing
global |K|=1. The fit still requires one temporal cycle and its other recorded
quality screens. n1WIG depth inversion uses the smaller positive n=1 Matsuno
gravity-speed root; exact Nyquist fits remain boundary cases.

The activity mask uses a configurable envelope cutoff in units of the signed
filtered field's reference sigma. The baseline is 1.5; later sensitivity runs
may select 1.0 or 2.0 with `DETECTION_SIGMA` in the Gaussian calibration and
NOAA/CMIP event-fit launchers, or `--detection-sigma` in their Python workers.
The candidate peak minimum follows that same cutoff. Each choice changes the
component footprints and candidate population, so calibrate its full-record
null threshold separately before comparing significant events. The selected
value enters the manifest and event-catalog provenance and yields a distinct
output path. The 1.0 and 2.0 calibrations have not been run.

For the 1.5-sigma CMIP6 Historical/AMIP ensemble, the repeatable
`run_cmip6_gaussian_calibration.py` batch worker reads a three-column
`model, experiment, member` TSV under `data/`. Its Slurm array index selects
one simulation. It plans the WK99-style run, prepares its background, computes
100 simulation chunks in parallel, reduces per-wave record-maximum thresholds,
runs the two injection symmetry groups, and finalizes the report. Each stage
uses the existing `gaussian_event_calibration.py` worker and its checkpoints;
rerunning a failed array index reuses complete chunks. From an allocation:

```bash
python scripts/run_cmip6_gaussian_calibration.py \
    --tasklist data/wave_significance/cmip6_calibration_jobs/simulations_1990-2010_wk99_ds1p5.tsv \
    --task-index 0 --max-parallel 32
```

The worker's explicit settings are 1,000 calibration plus 1,000 held-out
validation records per model/wave, 100 records per chunk, the NOAA seed and
five box profiles. Its CMIP6 output enters through
`data/wave_significance/gaussian_calibration_wk99_121_cmip6`; configure that
path as a symlink to storage with room for the ensembles. On Levante it points
to scratch because the home quota cannot hold all 59 runs. Full results need
transfer to permanent storage after completion. Run
`python scripts/preserve_cmip6_calibration_summaries.py --expected 59` after the
Slurm array to copy complete runs' thresholds, observed significance flags,
injection summaries and provenance to the home-backed
`data/wave_significance/cmip6_calibration_summaries/`. The raw null chunks and
spectra remain on scratch; the copier exits nonzero if fewer than 59 reports
are complete.

For example, later NOAA sensitivity runs can use the same cutoff in the
calibration and fit launchers:

```bash
DETECTION_SIGMA=1 BACKGROUND_METHOD=wk99_121 bash scripts/submit_gaussian_event_calibration.sh
DETECTION_SIGMA=1 FIT_SUPPORT=component bash scripts/submit_noaa_olr_events.sh
```

Repeat with `DETECTION_SIGMA=2` for the upper sensitivity case. Select fitted
events by applying the matching calibration run's peak-amplitude threshold
directly to each fit catalog's `peak_amplitude_W_m2` (`>` the threshold).
Candidate IDs can differ slightly between the calibration and fitting paths,
so do not join significance flags to fits by event ID. For the completed CMIP6
1.5-sigma runs, this direct peak test reproduced the reported significant
counts in all 295 model/experiment/wave catalogs. When multiple complete fit
catalogs exist for one wave, pass the intended `run_id` to
`load_event_catalog` so it does not choose one implicitly.

`FIT_SUPPORT=component` in the NOAA/CMIP event-fit launchers (or
`--fit-support component` in a worker) fits the active envelope pixels of each
detected component with equal weights. This mode has no fixed time/longitude
window or Gaussian weight scales. The default `window` mode preserves the
earlier weighted fits. The mode and activity cutoff are both recorded in the
catalog configuration and filename. In either mode, event-peak significance is
a separate full-record null test; a high component-fit R2 is not another
false-alarm guarantee.

There is no equivalent-depth cutoff. The n=1 WIG branch metadata is `wig`, mode
1; the box is unchanged by this metadata correction. Its upper frequency is the
daily Nyquist value, 0.5 cycles/day. Injection reports mark carriers at or above
Nyquist as unresolved (`carrier_below_daily_nyquist=False`), separately from
resolved in-band controls; amplitude significance does not establish propagation
direction at the sampling limit.

The completed NOAA Gaussian-field run used this WK99-style background command
(do not resubmit the existing run):

```bash
MODEL=NOAA_OLR EXPERIMENT=observations BACKGROUND_METHOD=wk99_121 \
    N_CALIBRATION=1000 N_VALIDATION=1000 CHUNK_SIZE=100 MAX_CONCURRENT=4 \
    bash scripts/submit_gaussian_event_calibration.sh
```

These are Gaussian random fields whose target spectrum is the shared, pooled
WK99-style 1-2-1 estimate documented above. This is distinct from the optional
`BACKGROUND_METHOD=gaussian`, which selects the earlier Gaussian-kernel spectrum
smoother. Both generators produce Gaussian fields.

Five profiles give 100 simulation chunks, with the first chunk per profile used
as a pilot and reused later. Each wave has 1,000 calibration and 1,000 disjoint
validation records. The default injection matrix contains 1,440 trials across
five wave families. Raw null scores, record maxima, thresholds and uncertainty,
held-out false-event counts, spectra, histograms, observed flags, recovery and
provenance remain saved as before. The primary target remains 1% probability of
any null event over a full source/wave record search, separately per wave.

The NOAA run `gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307`
completed with 100/100 chunks and 1,440 injection trials. At the primary 1%
record-level cutoff, the results are:

| Wave | Threshold (W m-2) | Significant/detected candidates | Held-out false-event probability (95% CI) |
| --- | ---: | ---: | ---: |
| Kelvin | 29.25 | 19/2897 | 1.4% (0.8–2.3%) |
| n1ERW | 30.63 | 9/1029 | 1.3% (0.7–2.2%) |
| MRG | 21.54 | 6/2964 | 1.6% (0.9–2.6%) |
| n0EIG | 19.39 | 8/3660 | 0.3% (0.1–0.9%) |
| n1WIG | 12.09 | 12/3836 | 0.8% (0.3–1.6%) |

These intervals condition on the fitted null. In controlled, resolved and
in-band injections, no 2-sigma packet passed the 1% rule; 8-sigma recovery was
50/80 Kelvin, 70/80 n1ERW, 70/96 MRG, 34/48 n0EIG, and 28/64 n1WIG. These
cases sample selected depths, wavenumbers and widths, so they are not an estimate
of population-wide completeness. Significance does not validate a depth fit.

The revised 1990–2010 NOAA fitting array `28584435` completed all five catalogs.
Every catalog matched its calibration candidate IDs, times and longitudes;
peak amplitudes differed by at most 9e-6 W m-2 from operation order. The
intersection of the 1% significance flag and existing `keep_for_histogram`
fit-quality flag is **14 Kelvin, 0 n1ERW, 4 MRG, 0 n0EIG, 0 n1WIG**. All 36
significant candidates excluded by fit quality had weighted R2 below 0.3; the
branch inversion, search-edge and cycle screens were not their limiting factor.
Changing only the record-level cutoff to 5% gives 18, 0, 5, 2 and 3 retained
fits in that wave order. These are small, conditional samples, not yet robust
five-wave depth distributions.

A read-only fit-window sensitivity check on the significant ER/EIG/WIG candidates
found that a 120-degree longitude window with half-width 5 days for EIG/WIG and
45 days for ER would pass the current R2 threshold for 5/9 n1ERW, 8/8 n0EIG,
and 12/12 n1WIG candidates. These selected-event diagnostic refits alone do
not validate completeness or depth accuracy; controlled checks follow below.

The controlled 8-sigma, resolved/in-band known-depth validation then compared
the original, shorter-window, and component-pixel fits. For n1ERW, n0EIG, n1WIG,
component fits were quality-screened and within a factor of two of known depth
in 63/80, 29/48, and 24/64 trials, counting packet detection misses. Median
relative depth errors among quality-screened, significant detections were 5%,
9%, and 12%. But the strongest candidate in 8/8, 8/8, and 7/8 reused null
backgrounds also passed the component-fit R2 screen. Those null maxima are not
necessarily significant; the eight controls do not estimate a 1% false-fit
rate. The result shows that R2 on a threshold-selected footprint is descriptive,
not independent evidence of wave identity. Products are under
`data/wave_significance/noaa_component_fit_validation_b898b1ae4307/`.

## Fitted NOAA Gaussian-background event catalogs

`fit_noaa_background_events.py` regenerates one held-out WK99-style Gaussian
record from the completed NOAA calibration's saved spectrum, template and seed.
It calls the same `run_event_catalog` detector, raw component-support sinusoid
fit, quality screens and equivalent-depth inversion as the observed NOAA fit.
The CSV uses the same columns as an observed catalog. Its row identity names
the Gaussian null and realization; its JSON sidecar records the ensemble, seed,
input files and fitting configuration. The output directory is
`<ensemble>/component_fit_catalogs/`, under `data/wave_significance/`.

From the project root on Levante:

```bash
MANIFEST=data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307/manifest.json
PYTHON=$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PYTHON" scripts/fit_noaa_background_events.py \
    --manifest "$MANIFEST" --wave Kelvin --realization 1000 --max-events 3
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 "$PYTHON" scripts/fit_noaa_background_events.py \
    --manifest "$MANIFEST" --wave Kelvin --realization 1000
"$PYTHON" scripts/run_noaa_background_fit_batch.py \
    --manifest "$MANIFEST" --start 1000 --stop 1200 --max-parallel 16
START_REALIZATION=1000 STOP_REALIZATION=2000 BLOCK_SIZE=200 \
    bash scripts/submit_noaa_background_fit.sh --dry-run
START_REALIZATION=1000 STOP_REALIZATION=2000 BLOCK_SIZE=200 \
    bash scripts/submit_noaa_background_fit.sh
"$PYTHON" scripts/audit_noaa_background_fit_catalogs.py --manifest "$MANIFEST"
```

The batch runner accepts disjoint index ranges within the 1,000 validation
realizations (1000–1999), starts separate workers for all five waves, and saves
per-worker logs under `<ensemble>/jobs/component_fits.<job-id>/`. Complete
catalogs are reused on retry. Workers remove only their own per-event
checkpoints after the final CSV and JSON have completed, keeping the permanent
output small. The Slurm launcher submits five 200-realization blocks by default,
each with 16 workers and at most four blocks running concurrently (24 CPUs,
40 GiB per block). The
calibration's original null scores and thresholds remain
untouched. A full null depth distribution requires fitting all chosen records;
peak-only calibration files contain no fitted wavenumbers or depths.
The audit checks the complete CSV schema, event counts and maxima against the
saved null calibration, valid retained depths and sidecar counts. It writes
`catalog_index.1000-2000.csv` and `catalog_summary_by_wave.1000-2000.csv` only
after all 5,000 fitted catalogs pass.
