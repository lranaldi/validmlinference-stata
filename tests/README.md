# tests/

Numerical parity harness: the Stata port is run on a fixed set of cases and
compared to **pre-saved Python reference values** held as per-case CSV
bundles under `fixtures/`.

## Layout

```
tests/
├── refresh_fixtures.py        Python: regenerates fixtures from the upstream package
├── run_all_tests.do           Stata:  loops over fixtures, runs each command, PASS/FAIL
├── test_load.do               Smoke test: library builds, all 7 stubs/commands load
├── test_ols.do                Standalone parity test for `ols` (Milestone 2)
├── test_ols_bca.do            Standalone parity test for `ols_bca` (Milestone 3)
├── test_ols_bcm.do            Standalone parity test for `ols_bcm` (Milestone 3)
├── test_ols_bca_topic.do      Standalone parity test for `ols_bca_topic` (Milestone 4)
├── test_ols_bcm_topic.do      Standalone parity test for `ols_bcm_topic` (Milestone 4)
└── fixtures/
    └── <case>/
        ├── meta.txt           KEY=VALUE per line (command, options, tolerances, dims)
        ├── input.csv          observation-level data with header row
        ├── ref_coef.csv       name, coef, se  (one row per coefficient)
        ├── ref_vcov.csv       row, col, vcov  (long format; d^2 rows)
        └── W.csv, S.csv, B.csv  matrix data, no header (topic-model cases only)
```

## Cases currently fixturized (Milestone 5)

| Case | Command | n | d | Notes |
|---|---|---:|---:|---|
| `ols_synthetic` | `ols` | 200 | 4 | 3 standard-normal regressors + intercept |
| `ols_remote_work` | `ols` | 16,315 | 2 | log_salary on wfh_rwo, real shipped data |
| `ols_bca_synthetic` | `ols_bca` | 500 | 3 | binary x1, fpr=0.05, m=2000 |
| `ols_bcm_synthetic` | `ols_bcm` | 500 | 3 | same data, multiplicative correction |
| `ols_bca_topic_synthetic` | `ols_bca_topic` | 300 | 6 | r=3, c=5, v=8, q=2, k=0.5 |
| `ols_bcm_topic_small_k` | `ols_bcm_topic` | 300 | 6 | k=0.005 → rho(Γ)≈0.40 (rho<1 branch) |
| `ols_bcm_topic_large_k` | `ols_bcm_topic` | 300 | 6 | k=20  → rho(Γ)≈1.6e+03 (BCA fallback) |

MLE cases (`one_step`, `one_step_gmm`) will be added in Milestones 6–7 once those commands exist.

## How to refresh

```bash
# 1. Regenerate the Python references (needs the `vmli` conda env)
/c/Users/loren/.conda/envs/vmli/python.exe ValidMLInference-stata/tests/refresh_fixtures.py

# 2. Run the Stata harness
# (from Stata, with current directory at the project root)
do ValidMLInference-stata/tests/run_all_tests.do
```

`run_all_tests.do` walks up from the current working directory to find the
project root, so it works whether you launch it from the project root, from
`tests/`, or from `src/`.

## Conventions

### Float precision

Two precision leaks were diagnosed and closed during Milestone 5; both fixes
are in the harness now but worth knowing about:

1. `refresh_fixtures.py` calls `jax.config.update("jax_enable_x64", True)`
   **before** importing the upstream `ValidMLInference` package. Without this,
   JAX silently downcasts to float32 inside `_reorder_intercept_first` and the
   Python reference loses ~7 sig digits.
2. `run_all_tests.do` uses `import delimited ..., asdouble` for every
   fixture CSV. Without `asdouble`, Stata stores variables as `float` and the
   subsequent `mkmat` / `st_data` calls see truncated values.

See `Notes/porting_decisions.md` for the full diagnosis.

### Coefficient ordering

Both Python and Stata standardize `e(b)` to "intercept first, then
alphabetical" (Python `_standardize_coefficient_order`; Stata
`vmli_alpha_perm` in `vmli_core.mata`). The harness aligns refs to Stata
coefficient names rather than relying on positional ordering, so any name
drift surfaces as a Mata error rather than a silent miscomparison.

### Tolerances

Default per-fixture (set in `meta.txt`):

- Closed-form (`ols`, `ols_bca`, `ols_bcm`, `*_topic`): `tol_coef=1e-8`, `tol_vcov=1e-7`.
- MLE (`one_step`, `one_step_gmm`): `tol_coef=1e-3` (will be added in
  Milestones 6–7).

Tighter / looser per-case tolerances can be overridden in the meta. Document
the reason in the `note=` field of the same meta file.
