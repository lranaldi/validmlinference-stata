# src/ado/

One `.ado` file per public Stata command. Each `.ado` is a thin wrapper that:

1. Parses the command syntax (`syntax depvar varlist [if] [in], ...`).
2. Builds the design matrix.
3. Calls the corresponding Mata routine in `../mata/`.
4. Posts results via `ereturn post b V`, with `e(cmd)`, `e(N)`, and option-specific scalars.

## Files

- `ols.ado` — implemented (v0.1.0)
- `ols_bca.ado` — implemented (v0.2.0; option `generated`, not `generated_var`)
- `ols_bcm.ado` — implemented (v0.2.0)
- `ols_bca_topic.ado` — implemented (v0.3.0; options `wmatrix`/`smatrix`/`bmatrix`)
- `ols_bcm_topic.ado` — implemented (v0.3.0)
- `one_step.ado` — stub only (Milestone 6)
- `one_step_gmm.ado` — stub only (Milestone 7)

## Conventions

- Coefficient ordering posted in `e(b)`: intercept first, then non-intercept names sorted alphabetically.
- Every estimation command sets `e(cmd)` to its own name (lowercase) so `replay` works.
- Options that need to round-trip (`fpr`, `m`, `generated_var`, `k`, `homoskedastic`, `nguess`, `maxiter`, `seed`) are echoed in `e()`.
