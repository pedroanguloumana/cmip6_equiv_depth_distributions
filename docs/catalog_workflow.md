# Catalog workflow

## Normal use: load once, search in memory

```python
from cmip_catalog import load_catalog
from src.paths import CMIP6_CATALOG_FILE

catalog = load_catalog("CMIP6", source=CMIP6_CATALOG_FILE)
candidates = catalog.search(
    variable_id=["rlut", "pr", "ua", "va", "ta", "hus"],
    experiment_id=["amip", "historical"],
    frequency="daily",
)
```

Load the complete catalog once per kernel. Keep `catalog` in memory and
change `.search(...)` calls as needed. No query is written to disk.
`search()` returns catalog rows, not a confirmed complete simulation.
`.match(...)` checks that all requested fields coexist in one simulation
with common coverage and compatible grid evidence. It is needed only when
constructing that joint cohort or file lists.

## The two active data tables

| File | Purpose |
| --- | --- |
| `data/cmip6_catalog.parquet` | One complete local catalog: 7,109,975 file rows in this snapshot; the existing approximately 101 MB file was renamed, not duplicated |
| `data/cmip6_headers.parquet` | One table of coordinate evidence reconstructed from the earlier header audit; no climate field arrays |

The header table is a dated study input for the six variables, two
experiments and 1980–2014 selection. It records actual pressure values,
horizontal coordinate fingerprints, calendars, sample paths and propagated
evidence. It is **not** a general statement about every catalog entry and
is not automatically refreshed when files change. Ordinary catalog
exploration does not need it. The notebook joins it only for the stronger
pressure-qualified selection. New variables, periods or archive changes
may require additional metadata; missing evidence must remain unknown.

The previous audit's reports, notebook snapshots, runner and redundant
caches are preserved under `data/archive/catalog-audit-2026-10-27/`.
They are historical records, not inputs to normal notebook execution.

## Why the archive symlinks exist

`data/cmip6` points to the CMIP6 archive, and `data/catalogs` points to the
DKRZ catalog directory. These are pointers, not copied data. The project's
`AGENTS.md` explicitly requires archive access through `data/` symlinks;
the library itself does not require them. They keep machine-local roots
out of the analysis code. Shared path names are in `src/paths.py`.

`data/catalogs` is used only when explicitly refreshing the complete
catalog. `data/cmip6` gives the selected archive paths a project-local
entry point. The current mapping derives the archive root from that link;
if moving to a mirror whose catalog uses a different URI prefix, supply
the correct explicit `remap_paths` mapping before using the saved evidence.

## Refresh the complete catalog explicitly

```python
from cmip_catalog import load_catalog
from src.paths import DATA_DIR, CMIP6_CATALOG_FILE

catalog = load_catalog(
    "CMIP6", source=DATA_DIR / "catalogs" / "dkrz_cmip6_disk.json",
)
catalog.df.to_parquet(CMIP6_CATALOG_FILE, index=False, compression="zstd")
```

This deliberately replaces the single local copy. Do not run it for every
query. It neither reads NetCDF field arrays nor refreshes the header table.
Keep the saved evidence's date/limitations in mind after an archive refresh.

## Scientific rules retained

- Daily frequency comes from metadata, independently of table ID.
- Pressure and hybrid model coordinates are different kinds; level count
  alone cannot make them interchangeable.
- Each complete match uses one model, experiment, member and compatible
  horizontal grid, with full catalog-range coverage of 1980–2014.
- Table preference is `Eday`, `day`, `CFday`, `EdayZ`, subject to coordinate
  checks and a feasible joint period; there is no maximum-level optimization.
- The four pressure profiles retain their own actual arrays. No vertical
  interpolation or common-level reduction is applied.
- Models must have a qualifying simulation in both experiments. Members
  can differ between experiments; equal labels do not prove physical pairs.

See [the recorded findings](cmip6_header_check.md) for counts, pressure
arrays, rejected alternatives and the limitations of representative evidence.

## What changed in the simplification

The full catalog has one descriptive filename. Per-query Parquet caches,
thousands of individual header JSON files, parallel header scans, progress
reporting and the Slurm runner have been removed from the active workflow.
Previously saved header evidence was consolidated without opening NetCDF
files again. The old material was archived, not discarded. The notebook
keeps the ordinary search and the stronger scientific match visible.

## Validation on 2026-10-29

All eight code cells ran in order in a fresh Python/IPython process on Levante, with embedded outputs saved. The complete 7,109,975-row catalog was loaded once. The selected simulation/table/version/file/common-period combinations matched the original audit exactly: 55 AMIP simulations, 174 historical simulations, and 18 models represented in both (219 individual simulations). The 27 cases with different pressure grids and 1,314 in-memory variable file lists were also reproduced. No NetCDF files were reopened, no per-query caches were created, and no scheduler job was submitted.
