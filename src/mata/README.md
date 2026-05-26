# src/mata/

Mata library files containing the numerical routines. The `.ado` wrappers in `../ado/` call into these functions.

## Current files

- `vmli_core.mata` — single Mata file holding everything ported to date (v0.3.0):
  - `vmli_version()`, `vmli_alpha_perm()` — version + intercept-first/alphabetical permutation.
  - `vmli_ols_fit()`, `vmli_ols_run()` — HC0 OLS.
  - `vmli_bc_prepare()`, `vmli_ols_bca_fit()`, `vmli_ols_bcm_fit()` and `_run` drivers — binary-label bias correction.
  - `vmli_bc_topic_prepare()`, `vmli_ols_bca_topic_fit()`, `vmli_ols_bcm_topic_fit()` and `_run` drivers — topic-model bias correction.
  - Private helpers `_vmli_post_results()`, `_vmli_find_col()`, `_vmli_topic_names()`, `_vmli_topic_read_inputs()`.
- `_build.do` — drops, recompiles, and re-indexes `lvmli.mlib`. Adds both `vmli_*()` and `_vmli_*()` patterns.

Future MLE machinery (`one_step`, `one_step_gmm`) will likely land in a separate `vmli_mle.mata` to keep `vmli_core.mata` from becoming unwieldy.

## Conventions

- All public functions are prefixed `vmli_` (ValidMLInference) to avoid namespace collisions.
- Functions return either a single matrix or a struct with named members; document each in the file header.
- Use `invsym()` for symmetric matrix inverses, `lusolve()` for general linear solves.
- HC0 variance: `V = invsym(sXX) * Omega * invsym(sXX) / n^2` — matches the Python `_ols_core` formula exactly.
- Build the library with `mata mlib create lvmli, replace` and `mata mlib index`.
