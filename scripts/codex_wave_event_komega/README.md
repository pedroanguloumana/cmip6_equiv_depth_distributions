# Five-wave event locations in WK coordinates

Run on Levante from the project root in the installed project environment:

```bash
python scripts/codex_wave_event_komega/extract.py
```

The extractor reads the 300 audited 2-sigma active-component catalogs, covering
NOAA, 34 Historical and 25 AMIP simulations over 1990--2010. No raw-data filtering,
event fitting or Gaussian-realization calculation is repeated. Existing outputs
are preserved rather than overwritten.

Outputs are in `data/event_threshold_comparison/wave_event_komega_sigma2/`:

- `histograms_and_wk.nc`: event-count histograms, subset totals, matched observed
  WK powers/backgrounds, fit bounds and resting-state 12/25/50 m dispersion curves.
- `catalog_index.csv`: exact source/member/wave catalogs, run IDs and SHA256s.
- `audit.json`: identities, bin closure, dimensional-unit and curve checks.

The primary selection is every successful fit with finite k and omega. It does
not select on inferred depth or branch-inversion validity. The second selection
is the existing `keep_for_histogram` flag. Counts use one event per bin without
amplitude weights. K-bin edges are integer-centered from -15.5 to 15.5; frequency
cells match the 1/96 cycles/day WK grid, with half-width endpoint bins at zero
and daily Nyquist. All selected events remain in the denominator. Exact Nyquist
events are counted separately and the top bin is marked as directionally ambiguous.

WK contours reuse the `codex_segment_heff` spectra: 96 daily samples, 60-day
overlap, symmetry-projected latitude powers averaged over 15 S--15 N, with a
shared S/A background smoothed by ten 1-2-1 passes in K and frequency. This is
the project's WK-style plotting adaptation. Event detection instead uses the
northern latitude-mean projection over 0--20 degrees after symmetry projection
over 20 S--20 N. The diagnostics share each source record but use different
latitude reductions. Event counts and spectral variance need not have identical
shapes, and their agreement is a descriptive comparison.

`notebooks/codex_wave_event_komega.ipynb` presents five-panel figures for NOAA
and equally weighted Historical/AMIP ensembles, with the two fit selections.
Main axes use signed global K and positive frequency; secondary axes show
dimensional k and omega. Dashed rectangles identify each family's fit-search
domain. Contours show the saved WK power/background ratio, and colored line
styles identify 12/25/50 m reference curves. Fitted coordinates are continuous;
the histogram binning does not imply integer-only fits or 96-day event durations.
