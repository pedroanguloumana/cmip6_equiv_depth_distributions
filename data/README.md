# Local data

The active metadata tables are `cmip6_catalog.parquet` and
`cmip6_headers.parquet`. `cmip6` and `catalogs` are local archive symlinks.
Older audit material is under `archive/`. No climate field arrays are
copied here. See [the workflow guide](../docs/catalog_workflow.md).

`noaa_olr` exposes the NOAA CDR archive. `cmip6_rlut` exposes the local historical
daily rlut files already regridded to `imerg_1deg_trop`. These symlinks are configured
separately on each machine; analysis code accesses them through `src/paths.py`.
The current CMIP directory contains 15 models for 2000--2010, including Gregorian,
no-leap and 360-day calendars. Derived model catalogs go in `cmip6_event_catalogs/`
and job manifests/logs/validation reports in `cmip6_event_jobs/`.

`cmip6_rlut_noaa1deg` is the separate output symlink for the new **1990–2010
historical and AMIP** regridding workflow. On Levante it currently points to
`/scratch/b/b382284/cmip6_rlut_noaa1deg` while `/work` is full. These temporary
products must later be copied to permanent storage and the symlink updated.
Files retain one explicit member, source grid/version, calendar and W m-2 units,
and use the exact global NOAA OLR grid, including cell bounds. Planning manifests,
metadata evidence, diagnostics and logs are in `cmip6_regrid_jobs/`.
See [the script instructions](../scripts/README.md#cmip6-rlut-on-the-noaa-olr-grid-19902010).

The active CMIP6 event workflow reads `cmip6_rlut_noaa1deg` for both historical
and AMIP simulations over 1990–2010. `cmip6_rlut` remains the older 2000–2010
input archive. Superseded event products and pilot CSVs/checkpoints are moved
through `data/archived` to the project's `archived/event_catalogs/` folder, with
move manifests. Current NOAA event catalogs remain active.
