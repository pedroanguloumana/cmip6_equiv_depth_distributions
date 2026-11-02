# Archived project material

This is the single warehouse for superseded data and analysis outputs. The local
`data/archived -> ../archived` symlink exposes it to analysis code. The former
`data/archive/` directory has been consolidated here. No scientific data or
completed products were deleted during this cleanup.

## Current cleanup: 20261030T104134Z_pre_box_profiles

The previous event catalogs used curved 1–200 m equivalent-depth filters. They
were superseded by the revised rectangular profile definitions and planned n=1
WIG analysis. Slurm jobs were cancelled before archiving. Active output locations
are available for recalculation; moving old results does not update the fitting
workers, detection settings or branch inversions for the new profiles.

| Location | Contents |
| --- | --- |
| `catalog-audit-2026-10-27/` | Previous catalog-audit caches, header JSON evidence, reports and workflow snapshots; moved from `data/archive/`. |
| `event_catalogs/20261030T104134Z_pre_box_profiles/noaa/` | Four completed NOAA catalogs, JSON sidecars and all event checkpoints. |
| `event_catalogs/20261030T104134Z_pre_box_profiles/cmip6_event_catalogs/` | 127 completed CMIP6 catalogs and four interrupted catalogs, including provenance and checkpoints. |
| `event_catalogs/20261030T104134Z_pre_box_profiles/cmip6_event_jobs/` | Previous CMIP6 fitting task lists, logs, pilot and validation records. |
| `event_catalogs/20261030T104134Z_pre_box_profiles/noaa_olr_event_jobs/` | Previous NOAA fitting job, pilot and resource records. |
| `wave_significance/20261030T104134Z_pre_box_profiles/` | Gaussian/marginal validation diagnostics, partial full Gaussian calibration, planned WK99 calibration, raw scores, cached spectra, snapshots and seeds. |
| `data_links/cmip6_rlut` | Legacy 2000–2010 IMERG-grid rlut symlink. Only the pointer moved; its external target is unchanged. |
| `move_manifests/` | Per-machine source-to-destination mappings and verification summaries for this cleanup. |
| `cleanup_snapshots/20261030T104134Z_pre_box_profiles/` | Pre-cleanup inventory and copies of the small documentation/path files changed during consolidation. |

Older notebooks and `event_catalogs/` archive batches already present here were
left in place. The local development clone contains only the small diagnostic
copies it previously held; the complete datasets and catalog/checkpoint archives
remain on Levante.

## Interpreting and retrieving archived results

Catalogs under this directory are historical method comparisons, not current
box-profile distributions. Open the CSVs explicitly when needed; the normal
catalog loaders search only active output locations. Keep each CSV with its JSON
sidecar and checkpoint directory. Original JSON contents, code snapshots, file
paths and recorded status strings are preserved. A historical `running` status
may describe a subsequently cancelled job; it is not a live queue status.

Unfiltered background spectra and seeds may be reusable when their dataset and
preprocessing settings match. Filtered-null scores, thresholds and event fits
remain specific to their original profiles and require recalculation for changed
masks. No automatic reuse or fallback into these archived outputs is enabled.

To reverse a move, consult the relevant manifest and restore the recorded entry
only while its original destination is unused; never overwrite a newer run.
Active NOAA/raw CMIP archives, the NOAA-grid CMIP files, catalog Parquet tables and
regridding provenance were kept at their original input locations.

Archive payloads are ignored by Git. This index and the small move manifests are
kept with the project documentation.
