# Two-sigma envelope components and event-implied equivalent depths

This run changes only the existing envelope activity threshold from **1.5 to
2 times the filtered OLR reference standard deviation**. It uses the same
1990--2010 NOAA record, 34 Historical / 25 AMIP simulations, native calendars,
five rectangular wave profiles, preprocessing, raw component-supported sinusoid
fit, branch-specific depth inversion and existing fit-quality screens.

For filtered OLR `q`, active pixels satisfy
`abs(Hilbert_longitude(q)) > 2 * std(q_reference)`. The reference excludes 150
days at each record end. Sigma is dataset- and wave-specific, so the threshold
does not impose a common absolute OLR amplitude across sources.
Eight-neighbor connected active pixels, with periodic longitude, define the
event component. The existing minimum duration of 2 days, minimum span of 10
degrees and exclusion of components touching the record endpoints are unchanged.
All component pixels enter the sinusoid fit with equal weights. A higher threshold
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
  peaked at 10.26 GiB. The plan permits up to 49 concurrent workers including
  NOAA; Slurm initially admitted two model shards plus NOAA during this run.
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

For example, inspect the exact NOAA Kelvin catalog on Levante with:

```python
import pandas as pd
from src.paths import DATA_DIR
from src.wave_events import load_event_catalog

index = pd.read_csv(DATA_DIR / "event_threshold_comparison/catalog_summary.csv")
row = index.loc[(index.model == "NOAA_OLR") & (index.wave == "Kelvin")
                & (index.detection_sigma == 2)].squeeze()
events = load_event_catalog(row.model, row.experiment, row.wave, run_id=row.run_id)
events[["event_id", "h_eff_m", "weighted_r2", "keep_for_histogram"]]
```

Higher R2 alone does not establish a better physical depth estimate: shrinking
support can make a sinusoid easier to fit while reducing sampled cycles. Compare
R2, support/cycle counts, search boundaries and depth distributions together.
Counts and distributions here are descriptive, with no background significance
interpretation.
`notebooks/codex_event_threshold_comparison.ipynb` reads the compact products
to compare counts, fit diagnostics and depth distributions at both thresholds.

The 2026-11-02 run completed all 300 new catalogs and passed the 600-catalog
audit (job `28660989`). The notebook retains 10 distribution figures, comparing
both thresholds with NOAA and equally weighted Historical/AMIP model means,
separately for all valid-depth fits and the existing quality-screened subset.
The full catalogs remain on Levante; the compact audit tables are also local.

## Example Hovmoller

On Levante, prepare the small plotting input with:

```bash
python scripts/codex_sigma2_events/prepare_hovmoller.py
```

This replays the full 1990--2010 NOAA Kelvin filter and 2-sigma detector, checks
all event IDs, peak times, longitudes and amplitudes against the saved catalog,
and caches `data/event_threshold_comparison/noaa_kelvin_sigma2_example.nc`.
It uses the existing quality-screened event nearest 1 July 2000, with event ID
breaking ties and no selection on depth. The chosen event is 1650, peaking
26 June 2000 at 168.5 degrees east, with h_eff=10.662 m and R2=0.909461.
The cache preserves all component pixels, verifies the saved fit and depth,
and records the source recipe. Existing caches are preserved rather than replaced.

The notebook adds three single-panel Hovmoller figures: signed OLR with 1.5/2
sigma envelope cutoffs, an envelope close-up with the active component, and
signed OLR with fitted negative-OLR crests inside that component. The cached
40-day, 181-longitude slice can be copied locally for fresh-kernel execution;
the archive record is needed only for preparation. No catalog refitting or
background calculation is performed. Carrier-phase propagation, rather than
the envelope boundary's motion, supplies the implied depth.

## Remake the older distribution notebooks

`remake_distributions.py` uses the audited 2-sigma catalogs to produce the old
notebooks' 1 m linear and 50-bin logarithmic distributions. It also replays raw
NOAA/CMIP6 OLR through the exact catalog filter and detector to recompute energy
on the new event footprints. Existing observed spectral powers are independent
of the event threshold and are reused. Gaussian comparisons remain deferred.

Run workers on Levante in compute allocations, with the project installed:

```bash
python -m scripts.codex_sigma2_events.remake_distributions --task-index 0
# Repeat for the 60 distinct records (indices 0..59), then:
python -m scripts.codex_sigma2_events.remake_distributions --summarize
```

Each worker needs the same memory allowance as the event-catalog workers
(14 GiB used for pilots). Cohort order is the distinct model/experiment/member
order in the 2-sigma rows of `data/event_threshold_comparison/catalog_summary.csv`.
Outputs and complete provenance receipts go to its `distributions_sigma2/`
subdirectory, preserving all 1.5-sigma products. The summary refuses incomplete
records and verifies counts, histogram closure and event-energy bounds before
writing the six notebook input tables and `audit.json`.

After successful pilots at indices 0 (NOAA), 2 (ACCESS-CM2 Historical) and 32
(HadGEM3-GC31-LL AMIP, 360-day calendar), submit the two pooled replay jobs and
their dependent audit with:

```bash
python -m scripts.codex_sigma2_events.submit_distribution_remake
```

Each shard runs up to 12 workers with 48 CPUs and 170 GiB for one hour. Complete
pilot outputs are reused. Submission IDs and per-record logs are retained under
`data/cmip6_event_jobs/sigma2_component/distribution_remake/`; a submission
receipt prevents accidental duplicate launches.

The 2026-11-03 remake is complete: pilot array `28679407`, full replay array
`28679573` and audit `28679575` all passed. All six compact input tables and
the audit are available on Levante and locally, with matching SHA256 hashes.
The following canonical notebooks now explicitly require the 2-sigma audit and
retain fresh-kernel rendered outputs:

- `notebooks/codex_heff_distribution_noaa.ipynb`: 26 figures, including a
  single-panel all-wave overview near the top. `OVERVIEW_SELECTION` chooses
  the existing screened subset or all valid-depth fits.
- `notebooks/codex_heff_comparison_cmip6.ipynb`: 35 count, individual-model,
  event-energy and threshold-independent observed-power figures.
- `notebooks/Compare_heff.ipynb`: 10 figures for the explicitly selected NOAA
  and ACCESS-CM2 Historical records; exact 2-sigma catalog run IDs are displayed.

The three original notebook versions are preserved in Git and in
`data/event_threshold_comparison/notebooks_sigma15_cb7a3c9/` on both machines.
Gaussian comparisons are omitted from the remake rather than mixing 2-sigma
observed events with 1.5-sigma nulls. No new null ensembles were calculated.
