# Local data

All analysis inputs enter through this directory. Configure archive symlinks
separately on each machine; do not copy raw climate fields into the repository.

## Active inputs

- `noaa_olr/`: NOAA CDR OLR archive.
- `cmip6/`: the CMIP6 data archive link.
- `catalogs/`: optional DKRZ catalog-source link, used only by the explicit
  catalog-refresh recipe. Normal queries and regridding planning use the local
  Parquet tables below.
- `cmip6_catalog.parquet` and `cmip6_headers.parquet`: complete catalog and saved
  header evidence used by the catalog notebook.
- `cmip6_rlut_noaa1deg/`: completed historical/AMIP rlut for **1990–2010**, one
  explicit ensemble member per model/experiment, on the exact NOAA 1-degree grid.
- `cmip6_regrid_jobs/`: planning manifests, diagnostics and provenance for those
  regridded inputs; these remain current after changes to the wave profiles.

The Levante NOAA-grid link currently points to
`/scratch/b/b382284/cmip6_rlut_noaa1deg`. These temporary products still need a
later transfer to permanent storage and a symlink update. Files retain source
member/grid/version, native calendar, cell bounds and W m-2 units. See the
[regridding instructions](../scripts/README.md#cmip6-rlut-on-the-noaa-olr-grid-19902010).

## Event outputs

Future CMIP6 catalogs go in `cmip6_event_catalogs/`, NOAA catalogs are written as
`NOAA_OLR.<wave>....csv` in this directory, and background calibration outputs go
in `wave_significance/`. Event job records go in `cmip6_event_jobs/` and
`noaa_olr_event_jobs/`. These output locations were cleared of superseded results
by moving their contents to the archive. The NOAA Gaussian calibration now uses
`BOX_PROFILES`, including n=1 WIG, and writes fresh calibration runs in
`wave_significance/`. The revised 1.5-sigma NOAA and CMIP6 event-fit catalogs
have completed. On Levante, the larger CMIP6 WK99 calibration ensembles use
`wave_significance/gaussian_calibration_wk99_121_cmip6`, a local symlink to
scratch storage. Scratch is temporary; preserve reduced calibration products
in permanent storage after the ensemble finishes.

## Archive

`archived -> ../archived` exposes the single canonical archive while keeping
analysis data access under `data/`. See [the archive index](../archived/README.md)
for old catalogs, checkpoints, job records, null diagnostics and move manifests.

The former `data/archive/catalog-audit-2026-10-27/` is now
`archived/catalog-audit-2026-10-27/` at the project root. The older 2000–2010
IMERG-grid input pointer is `data/archived/data_links/cmip6_rlut`; its external
files were not moved. Archived results require explicit paths and are never
silently substituted for current catalogs by the loaders.
