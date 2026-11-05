# Codex handoff — 2026-11-06

## OLR variance maxima and regression maps — 2026-11-06

The user requested WK-style average OLR structures for historical CNRM-CM6-1,
HadGEM3-GC31-LL and NOAA. In response to a reference-region question, they
explicitly chose **each dataset/wave's own maximum filtered-variance grid
point** and requested variance maps with the base points marked.

`notebooks/codex_cnrm_hadgem_noaa_olr_regression.ipynb` has two five-row by
three-column figures: filtered variance and zero-lag OLR regressions. Rows are
Kelvin, n1ERW, MRG, n0EIG, n1WIG; columns are NOAA/CNRM/HadGEM. Scales are
shared across datasets within each wave. A table records the 15 bases,
variances, local index sigmas, ratios to the NOAA maximum, and boundary flags.
Historical members are r1i1p1f2/r1i1p1f3; native standard/360-day calendars
remain unchanged, with NOAA proleptic Gregorian.

`src/wave_regression.py` contains the reusable calculation;
`scripts/codex_olr_regression/` contains the archive worker, Slurm launcher
and detailed README. It uses the current rectangular BOX_PROFILES with symmetry
projection and no depth cutoff, on the global 1-degree grid over 20 S–20 N.
Remove three seasonal harmonics and mean/linear trend. Taper only the predictor
input (Tukey alpha=.1), exclude its full tapered ends from all statistics, then
choose each interior variance maximum without meridional averaging. Symmetry
creates paired hemispheric maxima; use the northern tie for a consistent phase.
The OLR regression target is **unfiltered, untapered, without symmetry
projection**. Fit an intercept and scale to a **-1 sigma local filtered-OLR
excursion**, equivalent to +1 sigma convection. All interior days enter;
this is independent of the catalog's 2-sigma envelope threshold, and does not
apply the active-variance windows of Wheeler et al. (2000). Both relevant
papers and these adaptations are documented in the notebook and README.

All three 1990–1992 pilots in Slurm array `28757954` completed (36–39 s,
exit 0). All full 1990–2010 tasks in array `28757982` completed (3m12–3m35,
exit 0). Audits passed finite fields, variance/index-sigma closure, all 15
map maxima, file hashes, units, and native calendar sample counts. Statistics
use 6,902 days for NOAA/CNRM (1991-01-20 through 2009-12-12) and 6,804 for
HadGEM (1991-01-19 through 2009-12-12), after excluding 384/378 taper days
per end. Compact NetCDFs and audit JSONs are present locally and on Levante
under `data/olr_wave_regression/1990-2010/`; full filtered latitude fields
were not saved. These products are ignored by Git.

Kelvin maxima are near 0.5 N, 84.5/88.5/86.5 E for NOAA/CNRM/HadGEM.
MRG maxima are 6.5 N 171.5 E / 4.5 N 172.5 E / 5.5 N 166.5 E.
All three ER-box maxima are at 0.5 N. **HadGEM EIG maximizes at 19.5 N,
205.5 E, the northern domain edge**: the notebook flags this as domain-limited,
not an interior tropical EIG maximum. Broad-band filter labels and imposed
predictor symmetry do not independently establish wave identity. Daily Nyquist
ambiguity for the exact 0.5-cpd WIG boundary is also noted.

Two focused synthetic tests verify signed one-sigma regression scaling with an
intercept and variance-maximum selection while preserving target asymmetry.
They use standard-library unittest (the Levante environment has no pytest).
Fresh local kernel execution passed: 14 cells, two figures / 30 panels, no
errors. Both figures were visually checked; verified rendered outputs are
retained in the canonical notebook and in the ignored executed artifact under
`data/olr_wave_regression/`. Grid validation allows only the observed sub-2e-6
degree latitude roundoff in model coordinates, with no interpolation.
Existing local NOAA-distribution and remote Compare/WK notebook edits are
unrelated and must remain untouched.


## Historical CNRM--HadGEM depth and WK comparison — 2026-11-06

The user requested a new notebook comparing only the historical runs named
"CNRM-CM1" and "HadGEM3-GC1-LL", with one depth-distribution panel per wave
family and WK spectra. These were explicitly interpreted as the available
CMIP6 IDs CNRM-CM6-1 (`r1i1p1f2`, standard calendar, 7,670 days) and
HadGEM3-GC31-LL (`r1i1p1f3`, 360-day calendar, 7,560 days).

`notebooks/codex_cnrm_hadgem_historical_heff_wk.ipynb` is complete. It uses the
existing audited 2-sigma count tables and the matched saved WK cache, selecting
exact model/experiment/member identities. No new archive processing or fitting
was needed. `DEPTH_SELECTION="relaxed"` is the visible default (all successful,
positive, dispersion-valid depth fits); it can be changed to `"screened"`.
The 1--100 m / 1 m-bin figure and 0.5--10,000 m logarithmic figure each have
one panel for each of five waves, with just two source curves. All selected
events remain in the denominator; a table shows counts, medians, p10/p90 and
display coverage. Native calendars and members are retained explicitly.

A 2x2 WK figure shows models by row and symmetry by column. Colors use a shared
observed-power / own smoothed-background scale (below 1.1 white, above 2.5
saturated), and red curves show resting-state 12/25/50 m for the five analyzed
families. Spectra use the saved exact 96-day / 60-day-overlap latitude-power
analysis over 15 S--15 N and are independent of event threshold. The differing
event/WK latitude reductions and count/power normalizations are documented.

Fresh local kernel execution passed: 12 cells, three composite figures, no
errors. All three were visually checked. The canonical notebook retains the
verified outputs; its executed artifact is also local under
`data/event_threshold_comparison/pair_comparisons/`. Helper:
`/private/tmp/render_cnrm_hadgem_notebook.py`. For the default selection,
Kelvin median depths are 9.36 / 12.56 m (CNRM / HadGEM), ER 17.44 / 17.58 m,
EIG/WIG about 6.7 m in both, and MRG 30.35 / 35.58 m with broad tails.
These are descriptive single-member comparisons under the resting-wave inversion.
Existing local NOAA and remote Compare/WK notebook edits remain unrelated.

## Five-wave event histograms in WK coordinates — 2026-11-06

The user requested (k, omega) event histograms in Wheeler--Kiladis coordinates,
with a panel for each wave family, to compare their locations with WK-style
spectral enhancements. `notebooks/codex_wave_event_komega.ipynb` is complete:
six five-wave figures (NOAA, equally weighted Historical and AMIP means for
all successful fits, then the existing quality subset). Main axes use signed
global K and positive cycles/day; secondary axes give dimensional k and omega.
Each source/wave histogram is normalized before averaging models. Colors count
events once, without amplitude weighting, on one common logarithmic scale.

`scripts/codex_wave_event_komega/extract.py` read all 300 completed 2-sigma
catalogs and 60 saved 96-day WK spectra; job `28754065` completed in 24 seconds,
exit 0, on 2026-11-05. No raw OLR filtering, catalog refitting or Gaussian-null
calculation was repeated. The audit verifies source identities/input provenance,
674,102 successful fits / 182,484 screened events, dimensional coordinate
conversions, histogram closure and the 12/25/50 m branch reference curves.
Source/plotting inputs are under
`data/event_threshold_comparison/wave_event_komega_sigma2/`:
`histograms_and_wk.nc`, `catalog_index.csv`, `audit.json`. These were copied
locally with matching SHA256s. Extraction code was deployed as `23a9751`.

K cells are integer-centered (-15.5..15.5 edges); frequency centers are 0..48/96
cycles/day, with half-width endpoint cells. Nyquist events remain counted and
their bin is hatched; WK overlays exclude exact Nyquist. Primary maps include
every successful finite-coordinate fit, without a valid-depth requirement.
The second view uses `keep_for_histogram`. Dashed rectangles show carrier-fit
search bounds; binned cells can straddle those bounds.

Black contours are the saved WK power/background ratios 1.1/1.3/1.5 in each
family's appropriate symmetry. These use the existing segment analysis's
15 S--15 N latitude-power mean and ten K/ten frequency 1-2-1 smoothing passes.
Event detection uses a northern 0--20 latitude mean after symmetry projection
over 20 S--20 N. The differing reductions and count versus variance weighting
are explicit caveats. Orange curves are resting-state 12/25/50 m theory.

The descriptive overlap table excludes zero/Nyquist cells and reports both
event overlap with WK ratio >=1.1 and a geometric area reference, weighted for
partial cells within the fit-search domain. It is not a statistical null model.
For NOAA successful fits, event overlaps are Kelvin 98.1%, ER 96.8%, MRG 17.7%,
EIG 5.2%, WIG 37.8%; corresponding available-area shares are 54.6%, 55.3%,
9.2%, 14.1%, 41.6%. These shared-input comparisons do not validate physical
wave identity or depth accuracy by themselves.

The notebook ran from a fresh local kernel: 11 cells, six multi-panel figures,
zero errors; all six were visually inspected. Its executed artifact is also
preserved locally as `codex_wave_event_komega.reviewed.executed.ipynb` under
the output directory. Helper: `/private/tmp/render_wave_event_komega.py`.

Preserve existing user edits: local `codex_heff_distribution_noaa.ipynb` has its
overview selection changed to relaxed; Levante has edits in `Compare_heff.ipynb`
and `Wheeler-Kiladis_Spectra.ipynb`. None belongs in this task's commit.

## Remake the older distribution notebooks with 2-sigma events — 2026-11-03

The user requested that the older distribution notebooks use the completed
2-sigma catalogs, and an all-wave single-panel overlay near the top of the
NOAA notebook. Scope is `codex_heff_distribution_noaa.ipynb`,
`codex_heff_comparison_cmip6.ipynb` and `Compare_heff.ipynb`. The WK-segment
method does not use an event threshold. Gaussian comparisons remain deferred.
Original notebook copies are preserved on both machines in
`data/event_threshold_comparison/notebooks_sigma15_cb7a3c9/` and in Git.

The remade notebooks are complete with verified 2-sigma inputs and outputs.
The NOAA notebook has a top all-wave histogram with an explicit screened/relaxed
selector; all five curves use the corresponding 2-sigma event populations.
The comparison draft's undefined variables/unfinished cell were replaced with
an explicit 2-sigma audit-index/table selection, retaining NOAA and ACCESS-CM2
Historical as its chosen sources. No catalog loader/source hashes were changed.

`scripts/codex_sigma2_events/remake_distributions.py` replays raw data for all
60 source records and five waves, checks exact saved 2-sigma IDs/peaks, and
computes 1 m/log depth counts plus new component-footprint/fitted energies.
The old reusable energy helper now accepts an explicit threshold; its default
remains 1.5 for reproduction. Observed spectral-power tables are threshold-
independent and reused without any background comparison or computation.

Pilots `28679407_0/2/32` (NOAA, ACCESS-CM2 Historical, HadGEM3-GC31-LL AMIP
with a 360-day calendar) completed all five waves each, exit 0, in about four
minutes with 14 GiB allocations. Code commits `850b372` and `814d03e` are
pushed and deployed. Full replay array `28679573` used two 30-record shards,
12 workers per shard (48 CPUs/170 GiB/1 h); completed pilots are reused. Audit
`28679575` completed after both shards. The shards finished in 12m34s/12m23s,
and the audit in 5 seconds, all exit 0. All 60 record replays / 300 catalogs
passed source, exact detection, count, histogram-closure and energy-bound checks.
Submission receipts and per-record logs are in
`data/cmip6_event_jobs/sigma2_component/distribution_remake/`. Do not rerun this
completed calculation. Outputs are under
`data/event_threshold_comparison/distributions_sigma2/`: 600 count rows,
600 energy rows, plus the reused actual-data spectra (295 model and 5 NOAA rows).
The six summary CSVs and audit were copied locally and their SHA256s match.

All three notebooks ran in fresh local kernels: NOAA 18 cells / 26 figures,
CMIP6 18 cells / 35 figures, Compare 7 cells / 10 figures, zero errors. Visual
checks covered the NOAA all-wave overlay, count/energy/tail panels, model spread,
observed spectral power and the selected-source comparison. Verified outputs
were promoted to the canonical notebooks; executed artifacts are also retained
locally as `distributions_sigma2/<notebook>.sigma2.executed.ipynb`. The helper
`/private/tmp/render_remade_distribution_notebooks.py` records this execution.

For screened 12--50 m footprints, NOAA's energy share of total projected OLR
energy is 2.173% Kelvin versus Historical/AMIP model medians 0.866/0.888%,
0.253% EIG versus 0.0148/0.0123%, and 0.0523% WIG versus 0.0026/0.0027%.
These remain conditional descriptive diagnostics; no Gaussian comparison or
statistical significance claim was added. Spectral powers do not depend on the
event threshold and are labeled separately in the remade notebooks.

Levante has an unrelated user edit in `notebooks/Wheeler-Kiladis_Spectra.ipynb`.
Preserve it; no resets, stashes, commits or overwrites of that file are needed.

## Two-sigma event envelopes and event-implied depths — 2026-11-02

The user requested new NOAA and CMIP6 event catalogs with a 2-sigma envelope
activity threshold to assess whether stronger components improve inferred
depths. They clarified that this means the event envelope definition and
explicitly deferred all Gaussian-background work. Run event detection and the
existing raw active-component carrier fit/inversion only. Keep every candidate
and the existing validity/quality flags; preserve the 1.5-sigma catalogs.

No detector/fitter changes were needed: the standard NOAA/CMIP workers already
accept `--detection-sigma 2 --fit-support component`. The period, cohort,
profiles, preprocessing, sigma reference, fit grids and quality screens are
unchanged. A higher threshold changes component boundaries/IDs/support, so do
not filter the old fitted catalogs or join IDs across thresholds.

Five NOAA and ten CanESM5 Historical/AMIP three-event pilots completed in arrays
`28660725` and `28660727`, with all 45 fits successful. Checks against matching
1.5-sigma catalogs verified the same reference sigma and recipe except the
threshold/pilot limit. Peak observed/model worker RSS is about 10.25 GiB.

`scripts/codex_sigma2_events/submit.py` submitted the full computation:
NOAA array `28660987` (five waves, one concurrent), CMIP6 array `28660988`
(four shards, 85/85 Historical and 60/65 AMIP tasks, up to 12 workers each,
96 CPUs/170 GiB per shard), audit `28660989` dependent on both arrays. The
full 295-task cohort is pinned to
`data/cmip6_event_jobs/tasklists/simulations_waves_20261030_205901_2794174.tsv`.
Submission/shard records and logs are in `data/cmip6_event_jobs/sigma2_component/`.
The launcher refuses duplicate submissions once its `submission.json` exists.
All five NOAA tasks and four CMIP shards completed with exit 0. The final AMIP
shard's pending time limit was reduced from 4 h to 1 h; it completed in 15m12s,
without changing worker concurrency or any scientific settings. Audit `28660989`
completed in 14 seconds, exit 0. Do not resubmit this completed calculation.

New `.ds2.component.` CSV/JSON pairs go to the standard NOAA `data/` root and
`data/cmip6_event_catalogs/`. The audit compares all 600 complete catalogs
(both thresholds, five waves, 60 records), checking recipes, identity, flags,
counts, depth units and histogram closure. It writes compact summary and
5 m depth histograms to `data/event_threshold_comparison/`; no background
calibration or significance selection is used. All 600 catalogs passed the audit,
including exact source/wave recipe equality except the two threshold fields.
The compact CSVs and audit were copied locally. The comparison notebook,
`notebooks/codex_event_threshold_comparison.ipynb`, was executed locally from a
fresh kernel: now 16 cells, 13 embedded single-panel figures, no errors. Visual
checks covered the Kelvin/ER/MRG valid-depth curves, Kelvin/MRG screened curves
and all three added Hovmoller examples. Executed artifacts are preserved under
`data/event_threshold_comparison/`; the latest is
`codex_event_threshold_comparison.hovmoller.reviewed.executed.ipynb`.

The user requested cells illustrating the new envelope threshold with an example
Hovmoller. `scripts/codex_sigma2_events/prepare_hovmoller.py` replays the complete
NOAA Kelvin record and checks all event IDs, peaks, longitudes and amplitudes
against the full 2-sigma catalog. It selects the existing quality-screened event
nearest 1 July 2000, without selecting on depth: event 1650, peak 26 June 2000
at 168.5 E, h_eff=10.662 m, carrier c=10.225 m/s, R2=0.909461. Preparation job
`28664031` completed in 61 seconds (exit 0), also verifying the full selected
support, reconstructed R2 and branch depth. The small cache
`data/event_threshold_comparison/noaa_kelvin_sigma2_example.nc` (40 days x 181
longitudes) was copied locally. Source recipe and preparation hash are retained.
The job's exact source snapshot is in
`data/cmip6_event_jobs/sigma2_component/hovmoller.prepare.py`.

The five new notebook cells show signed OLR with 1.5/2-sigma envelope cutoffs,
an envelope close-up and the fitted negative-OLR crests inside the chosen active
component. The envelope includes both OLR signs. Carrier-phase propagation,
rather than the envelope boundary's motion, supplies depth. No catalogs were
refitted or background calculations performed. Existing cached products remain
preserved; the notebook needs only compact local inputs for fresh execution.

Across the fixed cohort, candidates decrease 808,059 -> 674,102 (16.6%), while
the existing quality-screened subset decreases 295,033 -> 182,484 (38.1%).
Median component-fit R2 rises and median sampled temporal cycles fall in all
300 paired source/wave comparisons; the screened fraction falls in 299/300.
The valid-depth interquartile range narrows in 281/300. NOAA's valid-depth
medians shift by less than 0.6 m for every wave. At 2 sigma, median depths in
metres are:

| Wave | NOAA | Historical model median | AMIP model median |
| --- | ---: | ---: | ---: |
| Kelvin | 14.8 | 11.6 | 12.2 |
| n1ERW | 20.0 | 17.7 | 17.4 |
| MRG | 54.0 | 35.1 | 33.8 |
| n0EIG | 8.8 | 6.7 | 6.5 |
| n1WIG | 7.9 | 6.7 | 6.7 |

Model entries are medians of individual source event medians, not pooled-event
medians. Paired model depth-median changes are generally below 1 m except MRG
(median change +3.5 m Historical / +2.9 m AMIP). NOAA MRG remains sensitive to
quality selection: 54.0 m for all valid fits versus 19.0 m screened, with 23.3%
of valid depths above the 200 m plot limit. Histogram denominators retain tails.
Sigma is dataset- and wave-specific; a common relative threshold does not impose
equal absolute amplitudes. Shorter support, resting-background inversion and
assigned branches limit interpretation. Higher R2 and narrower histograms do
not by themselves establish more accurate physical depths.

## Individual WK-segment depth bands — 2026-11-02

The user additionally requested a two-panel Historical/AMIP figure comparing
all wave families. It now appears near the top of
`notebooks/codex_segment_heff.ipynb`: shared axes, consistent wave colors,
5 m bins over 10--80 m, equally weighted means of individually normalized
positive-segment-excess distributions (34 Historical and 25 AMIP models).
NOAA observations are overlaid in both panels as dashed curves, using the
same wave colors and normalization; solid curves show the model means.
The notebook was executed locally from the compact audited products in a fresh
kernel: 16 cells, 18 figures, zero errors; the two-panel figure was visually
checked. The previous executed product is preserved, and the additional local
execution is `data/wave_significance/codex_segment_heff/codex_segment_heff.wave_panels_with_noaa.executed.ipynb`.

The user requested an additional spectral approach over 1990--2010, retaining
individual 96-day segments and integrating excess signal near dispersion curves.
They narrowed the scientific focus to 10--80 m and cross-family comparisons,
especially whether Kelvin and ER have the same apparent-depth distribution.
The implementation uses 5 m output bands and also compares broad 10--30,
30--50 and 50--80 m groups. This is a power-weighted spectral diagnostic,
not an event-count distribution or a statistical equivalence test.

Code is separate in `scripts/codex_segment_heff/`; the companion notebook is
`notebooks/codex_segment_heff.ipynb`. Inputs use the existing NOAA/CMIP `data/`
symlinks. The recipe averages latitude powers over 15 S--15 N after symmetry
projection, removes three annual harmonics and mean/trends, uses exact 96
daily samples with 60-day overlap, and corrects a 10% Tukey taper's variance
loss. A dataset-specific fixed background pools the full-record S/A spectra
and uses ten 1-2-1 passes in K and frequency, matching the simple WK plotting
smoother rather than the calibration's 10/20/40-pass adaptation. Depth bands
use fractional frequency-cell overlap at integer K, clipped to the current
wave boxes. Exact daily Nyquist is excluded because direction is ambiguous.

Outputs retain signed P-B, pointwise positive segment excess, positive ratio
excess, and positive excess of the mean spectrum separately. Positive clipping
has noise bias. Fractional bin allocation creates no new resolution: near
40 m at |K|=5--15, one Fourier-frequency bin spans roughly 1--4 m for Kelvin
and 16--27 m for n=1 ER. Adjacent depth bands and overlapping segments are
dependent. Mean signed excess is exactly a depth rebinning of the mean WK
spectrum. See the workflow README for units, geometry and further caveats.

The user explicitly requested committing/pushing and matching both checkouts
before continuing. Local snapshot `cbe5596` was merged with existing GitHub
history in `386d795`. Commit `77745a4` preserved Levante's newer NOAA/WK
notebooks, `Compare_heff.ipynb` and `src/wheeler_kiladis.py`. The newer local
agent-authored comparison notebooks were retained. Levante's previous source
state is preserved in `stash@{0}` (`pre-sync-20261102-preserved-source`); no
archive inputs or analysis outputs were stashed. Both checkouts were clean
at that commit. Commit `6b60a07` adds the new notebook and completion scripts.

All 31 existing project tests and four new scientific checks passed locally;
the four new checks also passed in Levante's project environment. NOAA and
ACCESS-CM2 historical pilots `28655973_0/1` completed in 82/91 seconds with
about 2.3/3.1 GB maximum RSS, matching integration geometry and passing signed
power/variance checks. Full array `28656027` completed all remaining 58 records;
Slurm permitted five concurrent user jobs. Job `28656087` passed the 60-record
audit (34 Historical, 25 AMIP, NOAA) and saved 4,200 depth rows and 63,210
wave/segment rows, then failed only because the new render wrapper imported
`nbformat`, absent from this environment. The wrapper now uses the project's
verified existing `jupyter_client`/JSON approach. Render-only job `28657184`
completed in 12 seconds: 14 cells, 17 standalone embedded figures, no cell
errors. No scientific calculation or summary was rerun for that repair.

Products are under `data/wave_significance/codex_segment_heff/`: 60 NetCDFs,
`depth_summary.csv.gz`, `segment_summary.csv.gz`, `audit.json` and the executed
notebook. The two CSVs, audit, NOAA NetCDF and executed notebook were copied
locally. All 300 normalized wave distributions sum to one; the saved 14-band
segment sums close within 3e-14 W2 m-4, and the audit fingerprint matches the
current analysis source. Five workflow tests, including the equal-25-m
Kelvin/ER response control, now pass on both machines (36 local tests total).
The delivered notebook retains the verified outputs and adds descriptive findings.

NOAA's positive-segment-excess shares in 10--30 m are Kelvin 56.7% and ER 52.6%;
model-median shares are 62.3%/53.0% in Historical and 59.7%/54.4% in AMIP.
Integrated signed Kelvin excess [W2 m-4] is NOAA 12.249 versus Historical/AMIP
medians 5.446/3.536; ER is 10.069 versus 12.595/13.462. Thus broadly similar
normalized shapes can conceal strength differences. These are apparent-depth,
one-record, one-member, conditional spectral comparisons, not significance or
physical-equivalence results. Visual checks covered the Kelvin/ER shapes,
broad-band comparison, depth-resolution plot and NOAA segment heatmap.

## NOAA–CMIP6 fitted-depth and power comparison — 2026-11-02

The user requested the remaining CMIP6 calculations for comparing fitted
equivalent-depth distributions and the allocation of OLR power against NOAA.
The completed 1.5σ active-component catalogs (295 full model/wave CSVs) and
59 WK99-style Gaussian calibrations are reused. New agentic code is confined
to `scripts/codex_cmip_heff/`, and the comparison notebook is
`notebooks/codex_heff_comparison_cmip6.ipynb`; no `src/` changes were made.

`run.py` replays each model's native-calendar `rlut` through the exact catalog
filter and detector, verifies event IDs and peaks against the saved fits, and
bins event-footprint and fitted-sinusoid energies by fitted depth (1 m bins,
1–200 m). It also bins the saved observed and smoothed-background Fourier
power by apparent depth on the same actual discrete filter boxes. Both energy
diagnostics use total projected OLR variance denominators. Outputs go to
`data/wave_significance/codex_cmip_heff_comparison/` on Levante. The pilot
ACCESS-CM2 Historical task `28635118` completed all five waves and passed
energy checks. The full 58-task array is `28635206`, up to five jobs at once.
All 59 model output pairs are complete. One full-array task failed initially
(index 14, CanESM5 Historical): its full catalog was present but
the worker also found an older `pilot3` catalog. The script now excludes pilot
files. The first repair `28635451` failed before reading data because manual
submission omitted the Slurm project-root export; it cancelled the original
dependent audit and notebook. The Slurm fallback was fixed, and only task 14
was resubmitted with explicit exports as repair array `28636090`, which
completed and matched all five saved catalogs.

The model-wide summary job `28636097` completed in 14 seconds after the full
array and repair. It ran `summarize.py` to audit all 59 pairs,
recount saved fitted depths in linear/log bins, and produce six compact tables:
three model-wide and three matching NOAA tables. Two initial notebook attempts
failed only in the execution wrapper (JSON source lines and forced non-inline
Matplotlib backend); the corrected job `28636369` completed with 14 cells and
15 embedded figures, zero cell errors. Its executed output was copied to the
local `notebooks/codex_heff_comparison_cmip6.ipynb`, and the six compact CSVs
were copied to local `data/wave_significance/codex_cmip_heff_comparison/`.
The notebook uses individual-model integrated
12–50 m tables plus NOAA and model-ensemble count, event-energy and spectral
plots. Model-specific Gaussian null **spectra** are included; null *event-fit*
catalogs are not needed for this direct observed/model comparison.

The output folder has 59 `*.event_energy.csv.gz` and 59
`*.spectral_power.csv.gz` pairs, six summary CSVs, and one executed notebook.
For quality-screened 12–50 m event-footprint energy, NOAA versus
Historical/AMIP model medians (% of total projected OLR energy) are Kelvin
3.180 versus 1.559/1.441, EIG 0.441 versus 0.031/0.025, and WIG 0.080 versus
0.006/0.006. ER is broadly similar. The spectral-power comparison confirms
larger NOAA shallow Kelvin/EIG shares; WIG's NOAA observed/background spectral
ratio is only 1.014, so its model difference is not itself a distinct spectral
peak. These are descriptive one-record/one-member comparisons, not formal
significance tests. Source hashes for `src/cmip.py`, `src/wave_events.py`,
`src/wave_significance.py`, `src/wave_profiles.py`, and `src/wave_fitting.py`
matched between local and Levante before launch.

## NOAA Gaussian-background event fitting — 2026-11-01

The user requested fitted equivalent-depth catalogs for NOAA's Gaussian
background events, using the same 1.5-sigma detector, raw active-component
sinusoid fit, fit-quality screens, inversion and CSV columns as the observed
NOAA catalogs. The repeatable worker, pooled runner and Slurm launcher are
`scripts/fit_noaa_background_events.py`,
`scripts/run_noaa_background_fit_batch.py` and
`scripts/submit_noaa_background_fit.sh`; instructions are in `scripts/README.md`.
The calibration run is `data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307`.
Outputs are one same-schema CSV plus JSON provenance per wave/realization under
that run's `component_fit_catalogs/` directory. Each CSV has
null identity and `member_id=null_<realization>`. The original calibration
thresholds and score files are untouched. Only the 1,000 held-out validation
realizations (indices 1000–1999) are fitted, for 5,000 final catalogs.

The old strict calibration manifest guard rejects later changes to the detector,
simulation and preparation source code. The saved `draw`, `filter`, profile,
spectral library and recipe hashes still match. Regeneration of held-out
realization 1000 for Kelvin and MRG exactly reproduced saved candidate counts,
peak maxima and reference sigma before fitting. The five-wave three-event
pilot (jobs `28620219`–`28620223`) completed, with the saved candidate count
reproduced for every wave, three successful fits per wave, matching observed
NOAA CSV columns and no leftover parts. Earlier pilot jobs `28620207`–
`28620211` failed before reading data because JSON lists were compared directly
with Python tuples; that comparison was corrected, with no output modified.

The five full realization-1000 pilots (`28620243`–`28620247`) completed with
exit 0. They produced 4,290 Kelvin, 1,404 ER, 4,108 EIG, 3,301 MRG and
4,495 WIG rows, each exactly matching the saved null candidate count and
maximum envelope peak. Every fit succeeded; their worker RSS was about
1.25–1.26 GiB, and each run took 1.5–3.7 minutes. The five permanent CSVs
occupied 8.4 MiB; checkpoint parts were removed only after final completion.

Full five-block Slurm array **`28620332`** completed all 5,000 worker catalogs
with zero worker failures. It used 200 validation realizations and five waves
per block (1,000 catalogs/block), at most four concurrent blocks, 24 CPUs and
40 GiB per block on `interactive`, with 16 simultaneous workers. Every block
exited 0 after about 3.5–3.6 h. Block logs are
`<ensemble>/jobs/component_fits.28620332_<task>.log`, and individual logs are
`<ensemble>/jobs/component_fits.<Slurm-child-job-id>/<wave>.<realization>.log`.
The five-record partial audit `catalog_index.1000-1001.csv` and
`catalog_summary_by_wave.1000-1001.csv` passed and contain 17,598 fitted
null events. Full audit job **`28620431`** completed with exit 0 in 2m52s.
`scripts/audit_noaa_background_fit_catalogs.py` wrote
`catalog_index.1000-2000.csv` and `catalog_summary_by_wave.1000-2000.csv`
after verifying all 5,000 catalogs have the observed NOAA schema and match
saved null candidate counts/maxima and sidecar counts. Its log is
`<ensemble>/jobs/component_fit_audit.28620431.log`.
The index has 5,000 rows; the catalogs contain 17,497,857 null events, of which
17,497,746 fitted successfully and 6,875,864 passed `keep_for_histogram`.
Retained counts by wave are Kelvin 1,097,898, ER 388,352, EIG 1,880,192,
MRG 1,062,948 and WIG 2,446,474. The permanent catalog directory occupies
8.2 GiB, with no leftover `.parts` directories. NOAA's original threshold
calibration and observed catalogs remain untouched. The next scientific step
is to compare observed and null equivalent-depth distributions separately by
wave, including absolute count and normalized-shape comparisons. No such
population significance claim has yet been made. The background distribution
need not be flat and remains conditional on the stationary Gaussian spectral
null.

## Active-component catalog recomputation — 2026-11-01

The user chose the existing raw filtered-OLR sinusoid fitter with support
restricted to the detected 1.5-sigma Hilbert-envelope event component. This is
`FIT_SUPPORT=component`: the fitter uses only active component pixels with equal
weights, a constant-amplitude sinusoid plus offset, and no Gaussian window or
`time_sigma_days`. It computes K, frequency, dispersion validity and `h_eff_m`
from each fit. The completed `.window.` catalogs are preserved; new filenames
contain `.component.`.

Local and Levante hashes of both fitting modules, both workers and the batch
runner matched before submission. The NOAA and CMIP input symlinks resolved.
Five NOAA pilots (array `28619656`) and ten CanESM5 Historical/AMIP pilots
(array `28619662`) all completed with exit 0. Each wrote three rows and all 45
pilot fits reported `fit_status=ok`; their full-record detected candidate
counts were NOAA 1,029–3,836 and CanESM5 1,096–4,026, consistent with prior
1.5-sigma detections. The first three events are not a quality-screen sample.

Full NOAA five-wave array `28619703` is submitted on `interactive`, with 8 CPUs,
14 GiB and 2 h per task, maximum four concurrent. The CMIP task list is
`data/cmip6_event_jobs/tasklists/simulations_waves_20261030_205901_2794174.tsv`:
295 unique tasks, 170 Historical and 125 AMIP, with all five waves for 59
simulations. Pooled jobs `28619729` (Historical) and `28619730` (AMIP) are
submitted on `interactive`, each with 112 CPUs, 200 GiB, 8 h and up to 16
worker processes. They depended on successful completion of pilot array
`28619662`, which has occurred. Both pooled jobs were RUNNING at the last
check. Each worker writes its own normal checkpoint parts and final CSV.
Monitor with `squeue -j 28619703,28619729,28619730 -r`; NOAA logs are in
`data/noaa_olr_event_jobs/logs/NOAA_OLR.28619703_<task>.log`, and pooled CMIP
logs are `data/cmip6_event_jobs/logs/CMIP6_component_pooled.<job>.log`, with
per-task logs under `data/cmip6_event_jobs/logs/CMIP6_pooled.<job>/`.

All five full NOAA tasks completed with exit 0. Their CSVs contain exactly the
expected 2,897 Kelvin, 2,964 MRG, 3,660 EIG, 1,029 ER and 3,836 WIG rows,
and every fit succeeded. `keep_for_histogram` counts are 767, 806, 1,409, 215
and 2,036 respectively. Applying the existing NOAA 1% peak thresholds gives
significant-candidate / significant-and-keep counts: Kelvin 19/14, MRG 6/3,
EIG 8/5, ER 9/5, WIG 12/9. Significance remains a separate amplitude test;
these quality yields are not false-alarm estimates.

At the first pooled CMIP survey, 35 full catalogs had completed across both
experiments and all five waves. Every CSV row count matched its sidecar and
all but one individual fit succeeded. Per-catalog candidate ranges were
2,744–3,272 Kelvin, 2,398–3,105 MRG, 3,081–3,816 EIG, 1,108–1,214 ER and
3,136–3,890 WIG. ER `keep_for_histogram` ranged 199–287 per catalog. These
are reasonable 1.5-sigma detector populations, not yet a full-sample result.

As the remaining CMIP tasks finish, verify every CSV row count matches its
sidecar `catalog_candidates`, inspect fit success and quality by wave, and
confirm all 295 final CMIP catalogs complete. To select significant fits for
depth distributions, apply each existing 1% background-calibrated peak cutoff
to that fitted catalog's `peak_amplitude_W_m2`; do not join old significance
flags by event ID. Compare quality yields with `.window.` catalogs, especially
for ER, before interpreting depth distributions. The 1.0- and 2.0-sigma
sensitivity catalogs and their calibrations remain future work.

## Active-component carrier-fit validation — 2026-11-01

The requested comparison of raw constant-amplitude and phase-normalized fits
on **identical active event pixels** is complete. Code, run instructions and
59 saved diagnostic CSVs are in `scripts/phase_fit_validation/` locally and on
Levante. Main Slurm array `28617664` completed all 15 NOAA/ACCESS-CM2
Historical/AMIP wave runs; one Historical WIG empirical-template selection was
corrected and rerun as `28618239`. Three ER K=-1/47-day control jobs
`28619451`–`28619453` completed. All final trial files passed row-count and
fit-exception checks; selected local/remote CSV hashes match.

The main comparison uses 3,240 paired noisy packet trials: three in-box
carriers, six smooth/irregular/observed envelope forms, two peak strengths
relative to the 1% cutoff and six independent backgrounds for each of 15
source/wave cases. There are also 270 clean paired cases and 90 paired null
maxima. Backgrounds are reused across packet cases, so trial rows are not
independent. The raw active-component fit has smaller median absolute K and
relative frequency errors than phase-only fitting for all five waves; the
90th-percentile errors also favor raw. ER median errors are global K 0.045
raw versus 0.048 phase, and relative frequency 4.5% raw versus 5.7% phase.
The block bootstrap in `paired_block_bootstrap.csv` favors raw for mean K error
in all five waves, but the ER frequency interval includes zero. This supports
raw active-component fitting as the current default for K/f, with phase-only
results as a sensitivity diagnostic. It is conditional on NOAA and one CMIP6
model (two experiments), not proof of universal model performance.

For ER K=-1 at about 89 days, the box/filter edge shifts frequency even in
clean controlled envelopes; observed irregular envelopes make both estimators
substantially worse. A K=-1/47-day control removes most clean controlled-shape
error, though some observed shapes still shift K. `er_k1_period_comparison.csv`
separates those effects. Existing component R2 >= 0.3 is not a significance
test: 85/90 null maxima pass it, while none exceeds the calibrated 1% peak
cutoff. No event catalogs were refit in this validation; the existing CMIP6
catalogs still use the old fixed windows.

## CMIP6 WK99 calibration completion — 2026-11-01

All 59 requested 1.5-sigma CMIP6 WK99-style Gaussian calibrations are
**complete**: pilot array `28599568` (2/2 tasks), main array `28599571`
(57/57 tasks), and reduced-result preservation job `28599710` all exited 0.
Every run has a complete report, five-wave thresholds and 1,440 injection
trials. `data/wave_significance/cmip6_calibration_summaries/index.json` records
all 59 preserved reduced runs on home storage. Full null chunks and spectra
remain under the scratch-backed `data/wave_significance/
gaussian_calibration_wk99_121_cmip6` link. These are per-model/experiment/
member/wave 1% full-record peak tests conditional on each stationary Gaussian
background, not joint across models or a proof of wave identity.

A read-only audit (`/private/tmp/audit_cmip_significance_fit_join.py` on the Mac)
checked all 295 wave calibrations against completed fit catalogs. Reports and
threshold files are present, but a few candidate IDs differ between the two
numerically near-equivalent preprocessing paths; matching significance flags
to fits by event ID is unsafe. A second audit
(`/private/tmp/audit_cmip_threshold_on_fits.py`) applied each calibrated 1%
threshold directly to its fitted catalog's envelope peaks: **all 295 selected
counts match** the calibration's reported observed counts, with no fitted peak
closer than 0.00033 W m-2 to its threshold. Use this direct peak test for depth
selection rather than joining saved flags by ID.

Across the 59 simulations, significant candidates / significant fits passing
the existing `keep_for_histogram` screen are Kelvin 649/461, MRG 159/125,
n0EIG 525/245, n1ERW 491/10, and n1WIG 481/292. n1ERW fixed-window fits
remain severely limited. These totals pool models/experiments and are not a
single population distribution or cross-model false-alarm guarantee. Held-out
false-record fractions vary by wave/model; they condition on the fitted null.

## CMIP6 WK99 calibration launch — 2026-10-31

The user requested the 1.5-sigma Gaussian background calibration for all
59 available CMIP6 simulations (34 historical, 25 amip), matching the NOAA
WK99-style 1-2-1 method and five rectangular wave profiles. The 59-row
simulation list is `data/wave_significance/cmip6_calibration_jobs/simulations_1990-2010_wk99_ds1p5.tsv`.
The repeatable `scripts/run_cmip6_gaussian_calibration.py` runs the existing
plan, prepare, 100 simulation chunks, threshold reduction, two injection groups
and final report for one simulation per Slurm array task. It uses 1,000 null
calibration plus 1,000 held-out validation records per wave, seed 20261031 and
32 concurrent simulation chunk processes per Slurm task. Complete chunks are
reused on retry. The runner and modified calibration worker were synced to
Levante and syntax/import checks passed there.

Home quota had only about 17 GB free while 59 NOAA-sized ensembles would need
roughly 80 GB; `/work` reported no free blocks. The CMIP-only output path
`data/wave_significance/gaussian_calibration_wk99_121_cmip6` is a Levante
symlink to `/scratch/b/b382284/cmip6_equiv_depth_distributions/wave_significance/gaussian_calibration_wk99_121_cmip6`.
The completed NOAA calibration remains in its original home directory. Scratch
is temporary; copy reduced final products to permanent storage after completion.

Pilot Slurm array **`28599568`**, tasks 0 (ACCESS-CM2 historical) and 34
(ACCESS-CM2 amip), **completed with exit 0**. Both have 100/100 chunks,
complete reports and 1,440 injection trials. Their observed-candidate counts
match the completed 1.5-sigma fit catalogs for all ten wave/experiment pairs.
All event IDs, peak times and longitudes match except one AMIP n1WIG event
(ID 2440) whose nearly tied peak shifts by one day; its envelope score differs
by only 3.3e-6 W m-2. This should be handled explicitly when joining later
significance flags to fit rows. The remaining 57 simulations are
submitted as array **`28599571`** with `--dependency=afterok:28599568` and
at most five concurrent array tasks; this matches the account's interactive
five-job and 256-CPU limits (48 CPUs, 80 GiB, 2 h per task). The dependency
released: five tasks were RUNNING and 52 pending at the last check, with no
failed Slurm tasks. Across the seven started runs, two reports were complete,
one was calibrated, and four were still simulating. Do not resubmit duplicates.
If a task fails, inspect
its per-run `jobs/pooled/` logs and rerun only that array index; existing chunks
will be reused. The 1% significance target is per source/member/experiment/wave
record, conditional on the fitted stationary Gaussian null, not joint across
all CMIP models or waves.

Preservation job **`28599710`** depends on the full array finishing (even if
some tasks fail). It runs `scripts/preserve_cmip6_calibration_summaries.py`,
copying reduced tables, observed significance flags, injection summaries and
code provenance for complete runs into home-backed
`data/wave_significance/cmip6_calibration_summaries/`. It exits nonzero if fewer
than 59 reports are complete. Full null chunks and spectra remain on scratch.

## CMIP6 fit completion — 2026-10-31

The two 1.5-sigma fixed-window CMIP6 pooled fitting jobs are **complete**:
Historical `28589555` (170/170 worker tasks, 2h39m) and AMIP `28589561`
(125/125 worker tasks, 1h56m), both Slurm exit 0 with no worker failures.
All 295 full catalogs are present for 59 simulations (34 historical, 25 amip)
and five waves. A read-only audit (`/private/tmp/audit_cmip_completed.py` on the
Mac) checked every JSON/CSV pair: complete status, 1.5-sigma window settings,
exact detected-candidate/CSV-row counts, and finite positive accepted depths;
there were zero issues. The minimum per-catalog optimizer success fraction was
at least 99.8% in each experiment/wave group.

Detected-event counts by wave across all models: Kelvin 2,363–3,310; MRG
2,108–3,151; n0EIG 2,484–3,848; n1ERW 886–1,245; n1WIG 2,941–4,152.
These are raw 1.5-sigma activity components, not calibrated significant events.
Fit-quality counts for n1ERW are very low (AMIP median 2; Historical median 0),
consistent with the fixed-window NOAA ER problem. All five-wave depth
distributions still need model-specific significance calibration and careful
fit-quality interpretation. No CMIP6 Gaussian null calibration was submitted.

## CMIP6 fit submission — 2026-10-31

The user requested 1.5-sigma CMIP6 Historical and AMIP event detection,
fixed-window sinusoid fitting and `h_eff_m` catalogs now. The 59 available
model/experiment/member simulations (34 historical, 25 amip) give 295 five-wave
tasks. The remote source/worker/launcher hashes matched the local development
clone. `data/cmip6_rlut_noaa1deg` resolved and the task list contains 170
historical plus 125 amip wave tasks. No CMIP6 Gaussian significance calibration
was requested or submitted; these catalogs will have fit-quality flags but no
model-specific calibrated significance flags until that separate work runs.

A 10-task CanESM5 pilot across both experiments/five waves is job `28589300`
on `shared` (three events fitted per wave). All ten tasks passed. Detected
candidate counts were 1,096–1,191 ER, 2,843–3,166 Kelvin/MRG/EIG and
4,020–4,026 WIG, broadly comparable to NOAA at 1.5 sigma. An initial full
array `28589321` was submitted then canceled
while pending, because this account's shared association permits only two
running jobs despite the array's requested concurrency. Its TSV task list was
preserved at `data/cmip6_event_jobs/tasklists/simulations_waves_20261030_205901_2794174.tsv`.
To use the available parallelism, `scripts/run_cmip6_event_batch.py` runs up to
16 ordinary catalog workers within each single-node allocation, with separate
per-task logs and normal catalog checkpoints. Compute-partition jobs `28589421`
and `28589429` were canceled while pending due to queue priority. The active
replacements are `28589555` (Historical, 170 tasks) and `28589561` (AMIP, 125
tasks), both RUNNING on `interactive` at the last check. Each requested 8 h,
112 CPUs and 200 GiB; this fits the account's 256-CPU interactive QoS limit.
The user asked to survey first completed full catalogs for reasonable candidate
counts. At the first survey, 14 full catalogs had completed across Kelvin, MRG
and n0EIG, with every catalog row count matching its detected candidate count
and nearly all fits returning successfully. Among the first 7–8 full detections
per wave, candidate ranges were Kelvin 2,744–3,272; MRG 2,398–3,105; n0EIG
3,035–3,816; n1ERW 1,090–1,214; and n1WIG 3,136–3,890. These are broadly
comparable to NOAA's 1.5-sigma detected counts, with no empty or explosive
populations. The first five completed MRG catalogs had 1,120–1,180 fits passing
the existing fit-quality screen; this is separate from significance, which
awaits CMIP6-specific calibration. The two Slurm jobs were still running at
this check. Queue and runtimes prevent a guaranteed completion time.

## Latest continuation — 2026-10-31

The user requested that later sensitivity runs allow envelope activity cutoffs
of 1.0, 1.5 (baseline), and 2.0 times the filtered field's reference sigma.
The current code accepts those three values in the NOAA/CMIP fit workers and
Slurm launchers and in the full Gaussian calibration plan/launcher. A chosen
cutoff changes active components, candidate IDs and fit support; it is recorded
in each catalog/manifest and produces a distinct output identity. Every cutoff
requires its own full-record null calibration and matching candidate join.
Only the 1.5-sigma NOAA WK99-style calibration has been completed; 1.0 and 2.0
have not been submitted. The current default fit mode remains the old window;
`FIT_SUPPORT=component` fits only active component pixels without Gaussian
time/longitude weights. The component-fit validation described in
`scripts/README.md` is mixed, so treat its R2 as descriptive and do not equate it
with independent significance evidence.

Three remote launcher dry-runs propagated 1, 1.5 and 2 to both calibration and
NOAA fits. Their Gaussian recipe hashes differ. The detector footprint test and
focused calendar/fitting tests passed on Levante. The Gaussian integration suite
passed 4/4 tests in Slurm job `28589094` (exit 0). See the dated continuation
below for calibration and fitting results.
The older `Current run` and `Next steps` sections are historical snapshots; the
calibration and revised fixed-window NOAA catalogs are complete.

## Continuation — 2026-10-31

The NOAA WK99-style Gaussian calibration is **complete**: report status
`complete`, 100/100 chunks, 1,440 injection trials, and all Slurm stages exit 0.
Primary 1% thresholds and candidate counts below remain correct. Held-out
full-record false-event probabilities (95% binomial CIs) are Kelvin 1.4%
(0.8–2.3%), n1ERW 1.3% (0.7–2.2%), MRG 1.6% (0.9–2.6%), n0EIG 0.3%
(0.1–0.9%), n1WIG 0.8% (0.3–1.6%). Calibration thresholds stabilize within
roughly 0–0.2 W m-2 from 600 to 1,000 records, though 95% quantile intervals
at 1,000 records are wider. For controlled resolved/in-band 8-sigma packets,
recovery was 50/80 Kelvin, 70/80 n1ERW, 70/96 MRG, 34/48 n0EIG, and 28/64
n1WIG. No 2-sigma packet passed. These are selected packet cases, not
population-wide completeness or depth-fit validation.

`src/wave_events.py` now uses `BOX_PROFILES` for all five families with matching
fit period grids; local fit K starts at |K|=1 for MRG/EIG. The user chose a
half-spatial-cycle quality screen so global |K|=1 fits may pass the 180-degree
window. `src/wave_fitting.py` now inverts the smaller positive n=1 Matsuno root
for WIG. Launchers include n1WIG. These changes are uncommitted. The modified
source/test/launchers were selectively synced to Levante; the remote WK notebook
and user-added plotting helper were not touched. Focused Levante tests passed:
6/6 calendar/fitting tests (Slurm `28584175`) and 8/8 significance tests
(`28584257`). Five full-record, three-event NOAA pilots passed (`28584324`).

The full revised NOAA event-fit array **`28584435` completed all five catalogs**
with exit 0. A read-only join check (`/private/tmp/check_noaa_fit_join.py` on
the Mac) matched all candidate IDs, times and longitudes to the calibration
flags; peak differences were below 9e-6 W m-2. At the 1% record-level cutoff,
usable fits (`keep_for_histogram`) number **14 Kelvin, 0 n1ERW, 4 MRG, 0 n0EIG,
0 n1WIG**. All 36 significant candidates that fail this fit-quality flag have
weighted R2 below 0.3, while the branch, edge and cycle screens pass. Loosening
the significance rule to 5% alone yields only 18, 0, 5, 2, 3 usable fits in
that wave order.

A read-only shorter-window diagnostic (`/private/tmp/diagnose_noaa_fit_windows.py`
on the Mac) refit selected significant events from saved filtered NOAA fields.
Using a 120-degree longitude window and 10-day EIG/WIG or 90-day ER duration
raised quality passes to 8/8 n0EIG, 12/12 n1WIG and 5/9 n1ERW; median absolute
log10 depth changes against the original fits were 0.109, 0.089 and 0.064 dex.
This same-event sensitivity is **not** fit-accuracy or null-fit validation.
Before replacing catalogs, test known-depth packet injections and null controls
with the new windows. The user has been asked whether to pursue that validation
or retain the sparse current result. CMIP6 calibration and revised fits have not
been launched.

## Goal and scope

Estimate event-level distributions of resting-background equivalent depth
(`h_eff_m`) from observed NOAA OLR and compare with CMIP6 historical/AMIP.
The immediate task is to finish and assess the **NOAA Gaussian-background
significance calibration for the revised rectangular profiles**, including n=1
WIG. The full run is already submitted: monitor/review it, do not duplicate it.
Only NOAA calibration is currently running; revised event fits and CMIP6
calibrations have not been launched.

## Environment and working rules

- Local development clone: `/Users/pedro/Projects/Analysis/cmip6_equiv_depth_distributions`.
- Primary/data checkout: `/home/b/b382284/Projects/cmip6_equiv_depth_distributions`
  on `levante.dkrz.de`, account `bb1153`.
- Edit locally; run archive-dependent work on Levante via SSH/Slurm. Read
  `AGENTS.md`. All analysis data access must use `data/` links and `src/paths.py`.
- Installed project package is actually named `src`. Levante Python:
  `~/.conda/envs/cmip6_equiv_depth_distributions/bin/python`.
  Local system Python lacks NumPy/xarray; remote system `python3` is 3.6.
- Prefer `ssh -o BatchMode=yes -o ConnectTimeout=12 levante.dkrz.de ...`.
  Sandbox DNS failures have required tool escalation; they were not SSH-key failures.
- Slurm jobs run independently of the laptop, Wi-Fi and this chat.
- Preserve user edits and completed results. No new figures unless requested.
  The user explicitly requested two panels for the WK notebook.

## Current run — checked 2026-10-30 13:21 UTC

Relative to the Levante project root:

```text
data/wave_significance/gaussian_calibration_wk99_121/
NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307/
```

The two lines above form one run directory, called `$RUN` below.

- `manifest.json` status: **calibrated**; **100/100 complete chunk receipts**.
- Preparation `28574243`, pilots `28574244`, main array `28574245`, and threshold
  reduction `28574246` succeeded. The threshold reduction took 24 seconds.
- Injection array **`28574247_0` / `_1`** was RUNNING; final report **`28574248`**
  was pending its dependency. No injection completion receipts or `report.json`
  existed at the last check. Recheck live status; these facts will age.
- Five waves × 2,000 records = 10,000 filtered wave-records: 1,000 calibration
  plus 1,000 disjoint validation records per wave. Chunk size 100, concurrency 4,
  seed `20261031`; pilot chunks are reused, not additional samples.
- All **15 tests passed**, validation job `28574066` (4m22s, exit 0).
  Log: `data/wave_significance/box_workflow_validation/tests.28574066.log`.
  Coverage includes both background estimators, five-wave end-to-end reduction,
  saved flags, injection reporting, and the n=1 WIG dispersion polynomial.

Useful status commands:

```bash
ssh -o BatchMode=yes -o ConnectTimeout=12 levante.dkrz.de \
  'squeue -u b382284 -o "%.22i %.25j %.12T %.12M %R"'
ssh -o BatchMode=yes -o ConnectTimeout=12 levante.dkrz.de \
  'sacct --starttime 2026-10-30 -j 28574245,28574246,28574247,28574248 -X --format=JobID,State,Elapsed,ExitCode'
```

On Levante:

```bash
cd /home/b/b382284/Projects/cmip6_equiv_depth_distributions
RUN=data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307
cat "$RUN/jobs/submission.txt"
tail -n 30 "$RUN/jobs/inject.28574247_0.log"
tail -n 30 "$RUN/jobs/inject.28574247_1.log"
```

The launcher schedules all stages automatically. It refuses a duplicate
submission when `jobs/submission.txt` exists. Do not rerun simulation or reduction
just because the final injection/report stage is unfinished. If a stage fails,
inspect its log and reuse complete chunks rather than regenerate the ensemble.

## Scientific choices and current profiles

`src/wave_profiles.py:BOX_PROFILES` is the calibration's explicit selection.
There is **no equivalent-depth cutoff**: an infinite upper depth disables the
curve mask, while branch metadata remains available for theory/injections.

| Key | Global zonal K | Period [days] | OLR symmetry |
| --- | --- | --- | --- |
| Kelvin | 1 to 15 | 2.5 to 14 | symmetric |
| n1ERW | -15 to -1 | 10 to 90 | symmetric |
| MRG | -10 to -0.25 | 2.25 to 10 | antisymmetric |
| n0EIG | 0.25 to 14 | 2.25 to 6 | antisymmetric |
| n1WIG | -15 to -1 | 2 to 4 | symmetric |

These are the user's revised bounds. WIG metadata was corrected to
`curve_name="wig", meridional_mode_number=1`; its box was preserved.
Global Fourier K is integer, so ±0.25 excludes K=0 just as ±0.5 would.
Exact K=0 was discussed as a separate n=0 category; none was implemented.
WIG reaches daily Nyquist (0.5 cpd); injection summaries distinguish carriers
at/above Nyquist from resolved controls. Amplitude significance is not proof of
propagation direction at that limit.

Preprocessing: daily 1990–2010, latitude -20 to 20 degrees, symmetry projection
followed by unweighted northern 0–20-degree averaging; remove three annual
harmonics, mean/trend, apply Tukey `pct=0.1`. No background-wind correction.
Detection uses the **longitudinal Hilbert envelope**, not a signed OLR peak:
8-neighbor components, periodic longitude, detection/peak sigma 1.5, minimum
2-day duration and 10-degree span; sigma excludes 150 days at each record end.

The null fields are stationary, zero-mean correlated Gaussian anomalies:
`IFFT(sqrt(P_background) * FFT(white_Gaussian_noise)).real`.
Variance matches the target spectrum **in expectation**, not by rescaling every
draw. The current run uses `BACKGROUND_METHOD=wk99_121`: 96-day segments with
60-day overlap, pooled symmetric/antisymmetric powers, xr_ccew 1–2–1 smoothing
(10/20/40 K passes below 0.2/from 0.2–0.3/above 0.3 cpd, then 10 frequency passes).
It is a documented WK99-based adaptation to the catalog's latitude means,
including regridding to the full-record FFT with pooled-variance conservation.
See `scripts/README.md` and `estimate_wk99_background_power` for exact details.
`BACKGROUND_METHOD=gaussian` instead means the older Gaussian-kernel **smoother**;
both methods generate Gaussian random fields. The launcher default remains
`gaussian`, so explicitly select `wk99_121` when reproducing the current method.

Primary significance targets **1% probability of any null event in the entire
21-year source/wave search**, separately per wave. This is not a pointwise 99th
percentile, catalog FDR, or a joint guarantee across models/waves. Apply a strict
`peak > threshold`, equivalent to the saved plus-one Monte Carlo p-value <= 0.01.
The 5%/10% alternatives and full threshold curves are also saved. Thresholds are
conditional on the fitted stationary null, not proof of modal identity or fit quality.

## Outputs already available and still to assess

Primary 1% results from `thresholds.csv` (envelope peak, W m-2; not yet reviewed
alongside injection completeness):

| Wave | Threshold | Selected / candidates | Held-out probability of any false event |
| --- | ---: | ---: | ---: |
| Kelvin | 29.253845 | 19 / 2897 | 1.4% |
| n1ERW | 30.630716 | 9 / 1029 | 1.3% |
| MRG | 21.543027 | 6 / 2964 | 1.6% |
| n0EIG | 19.389234 | 8 / 3660 | 0.3% |
| n1WIG | 12.087805 | 12 / 3836 | 0.8% |

Saved CIs and convergence must be inspected before declaring calibration
adequate. The primary rule is conservative and retains only 6–19 candidates per
wave; discuss the error/completeness tradeoff with the user before changing it.
Alternative thresholds can be assessed from saved scores without a new ensemble.

`$RUN` contains: `thresholds.csv`, `threshold_convergence.csv`,
`threshold_curve.<wave>.csv`, `null_records.csv`, `null_threshold_counts.<wave>.nc`,
`observed_significance.<wave>.csv` (p-values and separate significance flags),
raw null candidates in `chunks/task_*.events.csv.gz`, histograms/spectral summaries,
`wk99_background.nc`, `background_diagnostics.json`, source/analysis snapshots,
configuration and seeds. Only representative null fields are saved; other fields
can be regenerated from spectra/seeds. Finalization will add combined injections,
`injection_recovery.csv` and `report.json` (1,440 injection trials expected).

Next steps:
1. Confirm both injection jobs and final report succeed; inspect recovery by
   depth, width, amplitude and resolved/in-band status, plus threshold/validation CIs.
2. Report the completed NOAA calibration and retained populations to the user.
   Do not silently equate significance with a valid h_eff fit.
3. Before recomputing event-fit catalogs, update the legacy fitting workflow to
   the revised boxes, fit grids and n1WIG support. **This has not been done.**
   `src/wave_events.py:WAVE_PROFILES` still contains four old curved profiles;
   its `FIT_OPTIONS` are old, and `resting_wave_depth` supports only four families.
   Add/test the WIG branch inversion; `(omega/k)^2/g` applies only to Kelvin.
   Revisit the one-spatial-cycle screen: with a roughly 180-degree fit window it
   excludes |K|<2 despite the new boxes admitting K=1. This is a scientific choice.
   Changed masks require fresh detection and fits even at coincident event peaks;
   old curved-filter catalog IDs/fits cannot simply receive the new flags.
4. CMIP6 needs its own calibration per model/experiment/member, using a common
   method/error criterion rather than NOAA's numerical thresholds. Code supports
   model inputs, but no new CMIP calibration or fitting batch has been launched.

The user intentionally does not need event latitude centroids; time/longitude
plus meridional means/cross-sections are sufficient for later environment studies.
The Gaussian assumption was accepted as a working null after earlier marginal
surrogate checks; stationarity/local variance and background-estimator sensitivity
remain limitations. Changed masks require refiltering/detection of nulls; unchanged
background spectra/seeds may still be reusable when preprocessing/source match.

## Files, data and completed cleanup

- Calibration: `scripts/gaussian_event_calibration.py`,
  `scripts/submit_gaussian_event_calibration.sh`, `src/wave_significance.py`.
  Helpers now accept explicit profiles so calibration uses boxes without
  silently changing legacy fitting defaults. Do not edit scientific inputs/code
  while jobs run: recipe/environment fingerprints reject mismatches.
- Fits/loaders: `src/wave_fitting.py`, `src/wave_events.py`, `src/cmip.py`,
  `src/noaa_olr.py`; worker/launcher scripts and details in `scripts/README.md`.
- Remote `notebooks/Wheeler-Kiladis_Spectra.ipynb` is newer than the local copy.
  User-added `src/wheeler_kiladis.py` exists remotely but is absent locally.
  Fetch latest files before editing; never overwrite them from a stale local copy.
  Notebook spectra and plotting helper are separate from the calibration estimator.
- NOAA input: `data/noaa_olr -> /work/bb1153/b382284/OBS/NOAA_OLR`, file
  `OLR-Daily_v02r00_s19790101_e20251231.nc` (select 1990–2010).
- CMIP input: `data/cmip6_rlut_noaa1deg -> /scratch/b/b382284/cmip6_rlut_noaa1deg`:
  34 historical + 25 AMIP simulations, one member each, NOAA 1-degree daily grid,
  1990–2010, native calendars. Regridding is complete; retain provenance in
  `data/cmip6_regrid_jobs/`. Scratch products still need eventual permanent storage.
- Queries/planning use `data/cmip6_catalog.parquet` and `cmip6_headers.parquet`.
  `data/catalogs -> /pool/data/Catalogs` remains only for explicit refreshes;
  normal queries do not need that link. `data/cmip6` remains the raw-archive portal.
- Cleanup is **complete**; do not repeat it. `data/archive/` was consolidated into
  `archived/catalog-audit-2026-10-27/`. `data/archived -> ../archived` is the portal.
  Batch `20261030T104134Z_pre_box_profiles` preserved 4 NOAA + 127 complete CMIP
  catalogs, 4 partial CMIP runs, all checkpoints/job records, and old null studies.
  See `archived/README.md` and `archived/move_manifests/*.json`: 284,390 files,
  1,197,039,456 bytes verified preserved on Levante. No data was deleted.
  Legacy IMERG-grid input pointer moved to `archived/data_links/cmip6_rlut`;
  `CMIP6_RLUT_DIR` now routes through `data/archived`. External targets were untouched.
  Old CSVs/sidecars/status strings retain historical paths; normal loaders exclude
  archives and currently find no revised event-fit catalogs. Archive payloads are
  ignored by Git; the index and small move manifests are versionable.

## Working-copy caution

Local `main` latest commit is `aa55e4a` (fallback multiwave catalog implementation).
Current calibration, profile, cleanup/docs and test work is **uncommitted**, with
several new files untracked; no later commit/push was requested. Preserve the dirty
worktree; do not reset/clean or blindly stage archive payloads. Both machines have
been selectively synced, not fully mirrored. Remote Git may need a module on PATH.
Local origin is HTTPS; keep SSH/Git credentials local and inspect remotes before
publication. No analysis reruns or unrelated file changes are needed merely to
continue this handoff.
