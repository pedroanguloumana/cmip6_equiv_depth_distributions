# CMIP6 catalog header-check findings

Recorded from completed Slurm job **28480045**, finished **2026-10-27T13:40:39.643928+00:00**.
These are an archive snapshot, not a guarantee of current or future file availability.

## Request and selection policy

- Fields: `rlut`, `pr`, and available multi-level pressure profiles of `ua`, `va`, `ta`, `hus`.
- Experiments: `amip` and `historical`; daily frequency is selected independently of table ID.
- Entire simulation-calendar years 1980–2014; no prescribed pressure grid or interpolation.
- All ensemble members are eligible. Match identities include model, institution, experiment, member and initialization year.
- One feasible table/grid/version per field; table preference `Eday`, `day`, `CFday`, `EdayZ` is explicit.
- Horizontal coordinate fingerprints must agree. This is a strict comparison policy, not merely equal grid labels.
- The model intersection across experiments allows different eligible members in each experiment. Equal member labels do not prove paired initial conditions.

## Results

| Cohort | Simulations | Models |
| --- | ---: | ---: |
| AMIP | 55 | 19 |
| Historical | 174 | 25 |
| Models represented in both | 219 individual simulations | 18 |

There are 47 model/member-label combinations represented in both experiments. This is not a count of physically paired runs.
The notebook built 1314 variable file lists **in memory** for the two-experiment model cohort; it did not write model field arrays or text file lists.

Models represented in both:

- CanESM5
- CNRM-CM6-1
- CNRM-CM6-1-HR
- CNRM-ESM2-1
- MPI-ESM-1-2-HAM
- INM-CM4-8
- INM-CM5-0
- IPSL-CM6A-LR
- MIROC6
- MPI-ESM1-2-HR
- MPI-ESM1-2-LR
- MRI-ESM2-0
- CESM2
- CESM2-FV2
- CESM2-WACCM
- CESM2-WACCM-FV2
- GFDL-CM4
- SAM0-UNICON

## Actual pressure-coordinate findings

The chosen profiles use these observed grids:

- 8 levels (hPa): 1000, 850, 700, 500, 250, 100, 50, 10
- 19 levels (hPa): 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10, 5, 1

| Variable | Eday, 19 levels: simulations | day, 8 levels: simulations |
| --- | ---: | ---: |
| ua | 115 | 114 |
| va | 116 | 113 |
| ta | 116 | 113 |
| hus | 136 | 93 |

These counts describe all 229 qualifying simulations, not only the two-experiment model intersection.
**27 of those simulations have unequal pressure grids across the four profile variables; 202 have equal arrays.**
Equal level counts alone are insufficient. In this snapshot the eight-level grid is an exact subset of the nineteen-level grid, so a later analysis could select common levels without interpolating, but that would discard additional available levels. No such reduction was applied here.

The MPI-ESM1-2-LR pilot demonstrated `rlut`/`pr` from `day` with all four pressure profiles from `Eday` (19 levels).
Its inspected `EdayZ` profiles included 19 levels and 39 levels for `ua`, but lacked longitude coordinates: a denser vertical grid does not compensate for zonal averaging when a common latitude–longitude field is required.
Known model-level profile rows (8,873 rows labelled by `level_type`) were excluded before inspection. Daily sampling does not make model coordinates interchangeable with pressure coordinates; no pressure levels were inferred merely from a table name.

## Coverage, grid diagnostics and uncertainty

- The initial filtered catalog contained 100,436 rows. After excluding explicitly known model-level profiles, 268 simulations had preliminary joint availability and some period overlap.
- There were 2,990 exact simulation/table/variable/grid/version groups to inspect, using 1 representative file per group and reusing existing header caches. This is not a count of newly opened files.
- Representative inspection recorded no failures in this run. This does **not** establish that every archive file was checked.
- Final diagnostic counts were 25 `grid_mismatch` and 14 `incomplete_window`. These are diagnostic entries, not numbers of rejected models.
- Entries lacking a latitude/longitude fingerprint (including the sampled zonal means) were removed before final matching; those exclusions are separate from the grid-mismatch counts above.
- A strict fingerprint mismatch can reflect coordinate rounding or metadata differences as well as genuinely different grids. The rejected cases were **not** diagnosed further, so these counts do not by themselves prove that regridding is physically required. The resulting cohort is conservative under the stated policy.
- Whole-year labels let `360_day` models end on December 30. No one/two-day gap tolerance or Gregorian clamping was used. The final two-experiment cohort calendars were `{'standard': 89, 'proleptic_gregorian': 81, 'noleap': 49}`; it contains no qualifying `360_day` case under this selection.
- Full coverage is based on the union/intersection of catalog file ranges. Internal missing timesteps, missing values, every file's existence/coordinates, cell bounds and vector staggering remain unverified.
- Header evidence is propagated only within the exact dataset identity. The main audit used one sample; the earlier MPI pilot used up to two. Extra pilot cache entries must not be interpreted as additional qualifying simulations.

## Saved evidence and historical records

The normal notebook now loads one complete local catalog and reuses one explicit
header-evidence table. Ordinary searches create no caches and do not rerun a header
audit. See [the workflow guide](catalog_workflow.md).

- `data/cmip6_catalog.parquet`: the existing complete catalog, renamed without duplication.
- `data/cmip6_headers.parquet`: 17,325 file rows with saved coordinate evidence from 2,990 exact datasets.
- [Original report](../data/archive/catalog-audit-2026-10-27/catalog_runs/20261027T133421Z/report.json).
- [Machine-readable findings](../data/archive/catalog-audit-2026-10-27/catalog_runs/20261027T133421Z/header_check_findings.json).
- [Original executed notebook](../data/archive/catalog-audit-2026-10-27/catalog_runs/20261027T133421Z/Explore_Catalog.executed.ipynb).
- [Simplification comparison](../data/archive/catalog-audit-2026-10-27/simplification_comparison.json).

The old subset caches, individual header JSON files, batch runner, logs and notebook
backups are retained under `data/archive/catalog-audit-2026-10-27/`. They are historical
records, not dependencies of the active notebook. Paths inside the original JSON
reports describe the layout at that earlier run; the links above locate them now.

On 2026-10-29 the existing header JSON evidence was consolidated without opening any
NetCDF files. Every selected simulation/table/version/file/common-period combination
matched the earlier workflow. Counts remain 55 AMIP and 174 historical simulations,
with 18 models and 219 individual simulations in the two-experiment model intersection.
The header table is a static study snapshot, not a live file-validity cache. New
variables, periods or changed archive files need an explicit evidence update.

The simplified notebook exposes `complete_summary`, `models_in_both`, `cohort`,
`profile_levels`, `complete.diagnostics` and `filelists`. Detailed rendered tables
can be truncated; use these variables or the original aggregate report as appropriate.

## Lessons from the historical audit

Give the `catalog.search(...)` syntax first; distinguish candidate rows from a full-coverage, coordinate-qualified simulation match. Do not turn a syntax request into a large audit without first scoping the desired verification and ensemble policy.
A single model can contain many members: the MPI pilot covered 756 dataset groups and took about 17 minutes on a login node. For a genuinely small pilot, select one model **and one member** explicitly.
The final scheduler-owned run completed in about six minutes with four isolated worker processes and persistent caches. Avoid threads with netCDF4/HDF5; keep numerical-library thread counts bounded.
Use the dedicated `cmip6_equiv_depth_distributions` environment for this project. Mac orchestration can use system Python; it does not require the unrelated `xr_ccew` environment.
The Slurm run and its reports did not depend on an active Mac SSH/Wi-Fi connection. Check the actual report/job state rather than treating a scheduler dry-run start estimate as an actual start time.

The parallel worker pool and scheduler runner described above belong to the historical audit; the simplified notebook uses neither.
