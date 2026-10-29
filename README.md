# CMIP6 equivalent-depth distributions

Start with [`notebooks/Explore_Catalog.ipynb`](notebooks/Explore_Catalog.ipynb).
It loads the complete CMIP6 catalog once, searches in memory, and reuses
saved coordinate evidence for the daily AMIP/historical pressure-profile
selection. Routine exploration creates no per-query caches and runs no
archive header scan or scheduler job.

- [Catalog workflow, data layout and refresh instructions](docs/catalog_workflow.md)
- [Header-check findings and scientific limitations](docs/cmip6_header_check.md)

The existing Levante environment is `cmip6_equiv_depth_distributions`.
The project's installed package is currently named `src`; its shared paths
live in `src/paths.py`. The separate catalog library is `cmip_catalog`.
