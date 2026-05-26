# examples/

Replication do-files that mirror the three Python notebooks shipped with `ValidMLInference`.

## Planned do-files

- `remote_work.do` — mirrors `../../ValidMLInference-main/remote_work.ipynb`. Reproduces Table 1 of BCHS 2025 (Effect of Remote Work on Log Salary).
- `topic_model_example.do` — mirrors `topic_model_example.ipynb`. Reproduces Table 2 of BCHS 2025 (CEO time allocation and firm performance).
- `synthetic_example.do` — mirrors `synthetic_example.ipynb`. Compares OLS vs. bias-corrected estimators on simulated data with known mislabeling probability.

## Conventions

- Each do-file is self-contained: starts with `clear all`, sets the adopath, loads the shipped data, runs the analysis, prints/exports results.
- Tables are written to a `output/` subfolder alongside the do-file when applicable.
- No external dependencies beyond Stata 17 + the installed `ValidMLInference-stata` package.
