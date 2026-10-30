# Scripts

Routine catalog exploration needs no batch runner. The earlier one-off
audit runner is preserved in `data/archive/catalog-audit-2026-10-27/`.
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
and resting-background depth inversion for Kelvin, n=1 ER, MRG, and n=0 EIG.
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

Then submit the four independent array tasks when ready:

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
NOAA_OLR.Kelvin.19900101-20101231.k3to14_dk0p5.f0p0714286to0p4_nf17.tw10_lw90.ts5_ls45.<run-id>.csv
```

Here `tw/lw` are the time/longitude half-window and `ts/ls` their Gaussian widths.
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
`data/cmip6_rlut` link still exposes the 2000–2010 IMERG-grid archive for historical
reference; it is no longer used by the event worker.

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

This submits **236 simulation/wave tasks** (59 available model/experiment/member
pairs × four waves), with at most four running concurrently. It uses available
pairs rather than a Cartesian product that would invent unavailable experiments.
The TSV task list explicitly records model, experiment, member and wave; every
worker checks that member against the NetCDF header.

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
CMIP6.CanESM5.amip.r1i1p1f1.rlut.Kelvin.19900101-20101231.k3to14_dk0p5.f0p0714286to0p4_nf17.tw10_lw90.ts5_ls45.<run-id>.csv
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
Each move has a manifest preserving original paths. Current 1990–2010 NOAA
catalogs remain active; their pilot products can be archived separately. Source
NetCDF archives and completed regridding products are preserved.

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
