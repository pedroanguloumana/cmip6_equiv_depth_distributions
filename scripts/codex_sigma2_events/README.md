# Two-sigma envelope components and event-implied equivalent depths

This run changes only the existing envelope activity threshold from **1.5 to
2 times the filtered OLR reference standard deviation**. It uses the same
1990--2010 NOAA record, 34 Historical / 25 AMIP simulations, native calendars,
five rectangular wave profiles, preprocessing, raw component-supported sinusoid
fit, branch-specific depth inversion and existing fit-quality screens.

For filtered OLR `q`, active pixels satisfy
`abs(Hilbert_longitude(q)) > 2 * std(q_reference)`. The reference excludes 150
days at each record end. Connected active pixels define the event component;
all its pixels enter the sinusoid fit with equal weights. A higher threshold
changes boundaries, can split components, and shortens fit support, so these
are new catalogs, not a peak-amplitude subset of the 1.5-sigma fits.

There is **no Gaussian-background calculation or significance selection**.
The 1.5-sigma catalogs are preserved. The new complete CSVs and JSON sidecars
use `.ds2.component.` under `data/` (NOAA) and `data/cmip6_event_catalogs/`
(CMIP6). Fitted depths are not restricted to 10--80 m or clipped to 200 m.

## Run on Levante

The existing standard workers already support the required settings. First run
and inspect three-event pilots, with this project installed in the environment:

```bash
DETECTION_SIGMA=2 FIT_SUPPORT=component MAX_EVENTS=3 MAX_CONCURRENT=2 \
    bash scripts/submit_noaa_olr_events.sh
DETECTION_SIGMA=2 FIT_SUPPORT=component MODELS=CanESM5 MAX_EVENTS=3 MAX_CONCURRENT=2 \
    bash scripts/submit_cmip6_rlut_events.sh
python scripts/codex_sigma2_events/submit.py
```

The full launcher checks those 15 pilots, pins the same 295 model/wave tasks
as the completed 1.5-sigma calculation, and submits:

- Five NOAA waves, one concurrent task.
- Four CMIP6 shards (85/85 Historical and 60/65 AMIP tasks), up to 12 worker
  subprocesses per shard, using the existing `run_cmip6_event_batch.py`.
  Allocations are 96 CPUs / 170 GiB / 4 h each; previous observed/model workers
  peaked at 10.26 GiB. Five running Slurm jobs cover all 49 concurrent workers.
- An audit depending on successful completion of NOAA and all four shards.

Submission IDs and shard task lists are recorded in
`data/cmip6_event_jobs/sigma2_component/`. The launcher refuses duplicate
submission when its `submission.json` exists. The standard workers preserve
completed catalogs and resume their own checkpoints. No fitter source changes
are needed for this run.

## Assess the threshold change

`summarize.py` reads all 600 full catalogs: both thresholds for five waves in
NOAA and each of 59 simulations. It checks identities, source configuration,
reference sigma, CSV/sidecar counts, the unchanged fit-quality flags, depth
units and histogram closure. For each matched source/wave it requires the
recorded recipe to differ only in the two envelope threshold fields.

Results go to `data/event_threshold_comparison/`:

- `catalog_summary.csv`: source paths/run IDs, counts, fit R2, support sizes,
  search-edge and existing quality-screen fractions, depth quantiles and tails.
- `depth_histograms.csv.gz`: 5 m bins from 0--200 m for all successful,
  dispersion-valid fits and for the existing `keep_for_histogram` subset.
  Density uses **all** selected events as the denominator; the integral over
  the plotted range can be below one when events lie above 200 m. Tail counts
  are retained explicitly.
- `audit.json`: completed cohort and verification details.

Full source CSVs retain every detected candidate and its fit/validity flags.
`run_id` in the summary can be passed to the existing `load_event_catalog`
function to select the precise threshold-specific catalog; no arbitrary
newest-file choice is needed. Event IDs are local to a particular detection
threshold and must not be joined across the two populations.

Higher R2 alone does not establish a better physical depth estimate: shrinking
support can make a sinusoid easier to fit while reducing sampled cycles. Compare
R2, support/cycle counts, search boundaries and depth distributions together.
Counts and distributions here are descriptive, with no background significance
interpretation.
