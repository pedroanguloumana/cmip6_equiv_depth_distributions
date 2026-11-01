# Agent-authored NOAA event variance diagnostic

This diagnostic replays the saved 1.5σ NOAA Gaussian validation records and
redetects their active event supports. It joins each support by `event_id` to
the **existing** fitted component catalog. It does not recalibrate the null or
refit a sinusoid. The observed replay reads the original NOAA file through
`data/noaa_olr` and uses the exact fitted-catalog filter path. Run from the
project root in its installed Python environment on Levante:

```bash
python scripts/codex_event_variance/run.py --observed
python scripts/codex_event_variance/run.py --start 1000 --stop 1025
```

Repeat the second command over disjoint ranges through 2000. Outputs are
compressed CSVs under the saved calibration's `codex_event_variance/` directory
through `data/`. Each row represents one wave, realization, and selection
(`screened` or `relaxed`). Columns contain the denominator, total filtered
energy, assigned energies, and 1 m bin energies from 1 to 200 m.

After the observed and `1000:1001` pilots have passed, submit the remaining
held-out records with `bash scripts/codex_event_variance/submit.sh`. It launches
40 disjoint Slurm tasks (at most eight concurrent) and writes logs under
`codex_event_variance/jobs/`.

`support_*` is the filtered OLR sum of squares on active component pixels
relative to that record's filtered-field mean. `fitted_*` is the fitted
sinusoid's explained sum of squares, using its saved R² and each component's
local sum of squares. Both use the sum of squares of the two preprocessed,
projected OLR components as denominator. They are localization diagnostics,
not a unique decomposition of raw latitude–longitude OLR into physical waves.
