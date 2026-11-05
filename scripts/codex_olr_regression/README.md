# OLR variance and regression maps

This workflow compares NOAA CDR OLR, CNRM-CM6-1 historical `r1i1p1f2`, and
HadGEM3-GC31-LL historical `r1i1p1f3`, using daily 1990–2010 OLR on the NOAA
1° grid. Configure the existing `data/noaa_olr` and `data/cmip6_rlut_noaa1deg`
symlinks on Levante; code never reads an alternative archive root.

The five rectangular `BOX_PROFILES` match the event analysis, including their
symmetric/antisymmetric projections and **no equivalent-depth restriction**.
The geographical domain is 20°S–20°N at all longitudes. Three seasonal harmonics
and the full-record mean/linear trend are removed. Only the predictor input is
Tukey-tapered (alpha=0.1); both tapered ends are excluded from every statistic.
Dates and the excluded sample count are recorded separately for each native
calendar. No calendar interpolation is performed.

For each dataset/wave, the base is the grid point with the largest interior
filtered variance (population variance, units W² m⁻⁴). A symmetry-projected
field has tied northern/southern variance maxima; the northern point fixes
the phase convention. There is no averaging over a reference box or latitude
before choosing the maximum.

An ordinary least-squares regression with intercept relates the **unfiltered,
untapered OLR anomaly field**, without symmetry projection, to that point's
filtered OLR index. Maps show the fitted anomaly for a **−1σ** index excursion,
in W m⁻². Negative map values represent enhanced convection. The reference
index is centered and standardized using the same interior days. All interior
days enter; this is independent of the catalog's 2σ envelope detection rule.

This follows the filtered-index structure diagnostic in
[Wheeler, Kiladis and Webster (2000), section 3](https://webster.eas.gatech.edu/Papers/Webster2000a.pdf),
alongside the filtered-variance maps of
[Wheeler and Kiladis (1999), section 4](https://wind.mit.edu/~emanuel/tropical/wheeler_kiladis_1999.pdf).
It is an adaptation: current rectangular masks replace their depth-bounded
filters; all interior days replace their active-variance windows; and −1σ
replaces their −2σ plotting amplitude. No statistical-significance mask is
claimed. A strong band-limited pattern alone does not establish wave identity.

Install this checkout (`pip install -e .`) and `xr_ccew`. Run a small pilot
from the project root on a compute node:

```bash
python scripts/codex_olr_regression/run.py --task-index 0 --start-year 1990 --end-year 1992
```

Then submit all three full records:

```bash
bash scripts/codex_olr_regression/submit.sh
```

Tasks 0/1/2 select NOAA/CNRM/HadGEM. The launcher requests two CPUs and 24 GiB
per dataset for up to 30 minutes. It accepts `--start-year` and `--end-year`
to launch a separate pilot. Workers refuse to overwrite completed products.
Inspect `squeue` and the existing output/logs before resubmitting.

Outputs under `data/olr_wave_regression/<start>-<end>/` are one compact NetCDF
per dataset (variance, regressions, correlations, base coordinates, complete
filtered reference series, native dates, and filter bounds) and an audit JSON
with source metadata, code hashes, source size/mtime and product SHA256.
Raw filtered latitude fields are not saved. Copy these compact files to the
local `data/` tree to run `notebooks/codex_cnrm_hadgem_noaa_olr_regression.ipynb`.
The notebook embeds its figures and does not save standalone figure files.
Plotting needs Matplotlib, pandas and Cartopy with Natural Earth coastlines
at a scale of 1:110 million.
Its product-grid check permits the observed sub-2e-6-degree coordinate roundoff;
no extra interpolation or regridding is performed. The notebook flags HadGEM's
EIG maximum at 19.5°N as a domain-edge base point.
