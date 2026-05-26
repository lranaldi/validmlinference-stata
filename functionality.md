# ValidMLInference (Stata) — API reference

Stata-side API reference, mirroring `../ValidMLInference-main/functionality.md` but adapted to Stata syntax. All seven commands are implemented (v0.5.0, 2026-05-22).

| Command | Method | Variance | Parity status |
|---|---|---|---|
| [`ols`](#ols) | OLS | HC0 | ≤ 1e-8 / 1e-7 |
| [`ols_bca`](#ols_bca) | OLS + additive correction (binary label) | HC0 + fpr term | ≤ 1e-8 / 1e-7 |
| [`ols_bcm`](#ols_bcm) | OLS + multiplicative correction (binary label) | HC0 + fpr term | ≤ 1e-8 / 1e-7 |
| [`ols_bca_topic`](#ols_bca_topic) | OLS + additive correction (topic shares) | HC0 (uncorrected) | ≤ 1e-8 / 1e-7 |
| [`ols_bcm_topic`](#ols_bcm_topic) | OLS + multiplicative correction (topic shares) | HC0 (uncorrected) | ≤ 1e-8 / 1e-7 |
| [`one_step`](#one_step) | Joint MLE (normal / Laplace / Student-t outcome) | OIM (numerical Hessian) | ≤ 1e-3 / 1e-2 |
| [`one_step_gmm`](#one_step_gmm) | Joint MLE (Gaussian mixture outcome) | OIM (numerical Hessian) | ≤ 1e-3 / 1e-2 |

Parity status is reported as `max |db|` / `max |dV|` against the upstream Python implementation on the shipped fixtures.

---

## `ols`

**Status:** implemented (v0.1.0, 2026-05-18).

**Syntax:**
```
ols depvar varlist [if] [in] [, noconstant]
```

**Description:** OLS with heteroskedasticity-consistent **HC0** standard errors (no degrees-of-freedom adjustment). Mirrors Python `ValidMLInference.ols`. Coefficient ordering in `e(b)`: intercept first, then non-intercept names sorted alphabetically — matching Python's `_standardize_coefficient_order`.

**Note for Stata users.** `regress y x, robust` returns HC1 (the `(n)/(n-k)` inflated variant). `ols` returns HC0, so standard errors here are slightly smaller by a factor of √((n-k)/n). The choice matches the upstream Python package and is required for parity-test fixtures to align.

**Options:**
- `noconstant` — fit without an intercept.

**Posted estimates:** `e(b)`, `e(V)`, `e(N)`, `e(cmd)`, `e(cmdline)`, `e(depvar)`, `e(title)`, `e(vce)`, `e(vcetype)`, `e(properties)`.

**Example:**
```stata
use "ValidMLInference-stata/data/remote_work_data.dta", clear
gen log_salary = log(salary)
ols log_salary wfh_rwo
```

---

## `ols_bca`

**Status:** implemented (v0.2.0, 2026-05-18).

**Syntax:**
```
ols_bca depvar varlist [if] [in], fpr(real) m(real) [generated(varname) noconstant]
```

**Description:** Additive bias correction for OLS when one covariate is a binary AI/ML-generated label. Requires an external estimate of the classifier's false-positive rate `fpr` and the sample size `m` used to estimate it. Standard errors are inflated by an `fpr(1-fpr)/m` term to reflect uncertainty in `fpr`. Mirrors Python `ValidMLInference.ols_bca`.

**Options:**
- `fpr(real)` — estimated false-positive rate of the classifier. Must lie in `[0, 1]`.
- `m(real)` — size of the external sample used to estimate `fpr`. Set very large (e.g. `m(1e10)`) when `fpr` is treated as known.
- `generated(varname)` — name of the AI/ML-generated covariate. Defaults to the first variable in the varlist. (Renamed from Python's `generated_var` because Stata option names disallow underscores.)
- `noconstant` — omit the constant term.

**Posted estimates:** `e(b)`, `e(V)`, `e(N)`, `e(fpr)`, `e(m)`, `e(generated_var)`, `e(cmd)`, `e(cmdline)`, `e(depvar)`, `e(title)`, `e(vce)`, `e(vcetype)`, `e(properties)`.

**Example:**
```stata
use "ValidMLInference-stata/data/remote_work_data.dta", clear
gen log_salary = log(salary)
ols_bca log_salary wfh_rwo, fpr(0.009) m(1000) generated(wfh_rwo)
```

---

## `ols_bcm`

**Status:** implemented (v0.2.0, 2026-05-18).

**Syntax:**
```
ols_bcm depvar varlist [if] [in], fpr(real) m(real) [generated(varname) noconstant]
```

**Description:** Multiplicative bias correction for OLS when one covariate is a binary AI/ML-generated label. Options identical to `ols_bca`. Recommended over `ols_bca` for binary imputed labels per BCHS 2025. Mirrors Python `ValidMLInference.ols_bcm`.

**Posted estimates:** identical to `ols_bca` with `e(cmd) == "ols_bcm"`.

**Example:**
```stata
use "ValidMLInference-stata/data/remote_work_data.dta", clear
gen log_salary = log(salary)
ols_bcm log_salary wfh_rwo, fpr(0.009) m(1000) generated(wfh_rwo)
```

---

## `ols_bca_topic`

**Status:** implemented (v0.3.0, 2026-05-18).

**Syntax:**
```
ols_bca_topic depvar [Q_varlist] [if] [in], wmatrix(matname) smatrix(matname) bmatrix(matname) k(real) [noconstant]
```

**Description:** Additive bias correction for regressions with topic-model-generated regressors. The design matrix is `Xhat = (Theta, Q)` where `Theta = W * S'` are the document-topic shares. The correction is `b_corr = (I + Gamma) * b0` with `Gamma = k * sqrt(n) * (Xhat'Xhat)^{-1} * A`, where the topic-share block of `A` is

```
Omega = S * (B B')^{-1} * B * (B' ⊙ (B' * mean(W))) * (B B')^{-1} * S' - Theta'Theta / n.
```

Mirrors Python `ValidMLInference.ols_bca_topic`.

**Options:**
- `wmatrix(matname)` — n × c document-by-component matrix (row-stochastic).
- `smatrix(matname)` — r × c topic-by-component matrix (row-stochastic).
- `bmatrix(matname)` — c × v component-by-vocabulary matrix (row-stochastic).
- `k(real)` — bias-correction scaling parameter (typically `O(1)`).
- `noconstant` — omit the constant term.

The Stata option names are `wmatrix`/`smatrix`/`bmatrix` (not single letters) because Stata's `syntax` parser would otherwise collide on the `b()` shortcut for matrices and the internal `tempname b`.

**Coefficient ordering:** the r topic-share columns are named `topic_1`, ..., `topic_r`. Final `e(b)` ordering: `_cons` first, then non-intercept names sorted alphabetically (so `_cons q1 q2 topic_1 topic_2 topic_3` for two Q covariates).

**Posted estimates:** `e(b)`, `e(V)`, `e(N)`, `e(k)`, `e(W_matrix)`, `e(S_matrix)`, `e(B_matrix)`, `e(cmd)`, `e(cmdline)`, `e(depvar)`, `e(title)`, `e(vce)`, `e(vcetype)`, `e(properties)`.

**Example:**
```stata
matrix W = ...   // n x c
matrix S = ...   // r x c
matrix B = ...   // c x v
ols_bca_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0.5)
```

---

## `ols_bcm_topic`

**Status:** implemented (v0.3.0, 2026-05-18).

**Syntax:** identical to `ols_bca_topic`:
```
ols_bcm_topic depvar [Q_varlist] [if] [in], wmatrix(matname) smatrix(matname) bmatrix(matname) k(real) [noconstant]
```

**Description:** Multiplicative bias correction for regressions with topic-model-generated regressors. Uses the same `Gamma` as `ols_bca_topic`, but inverts `(I - Gamma)`:

- If the spectral radius `rho(Gamma) < 1`, the correction is `b_corr = (I - Gamma)^{-1} * b0` (computed via `lusolve`).
- Otherwise (numerically unstable), the routine falls back to the additive (BCA) formula `b_corr = (I + Gamma) * b0` and prints a warning.

Mirrors Python `ValidMLInference.ols_bcm_topic`.

**Options:** same as `ols_bca_topic`.

**Posted estimates:** same as `ols_bca_topic` with `e(cmd) == "ols_bcm_topic"`.

**Example:**
```stata
ols_bcm_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0.5)
```

---

## `one_step`

**Status:** implemented (v0.6.0, 2026-05-26).

**Syntax:**
```
one_step depvar varlist [if] [in] [, generated(varname) homoskedastic dist(name) df(#) noconstant]
```

**Description:** Maximum-likelihood joint estimation treating one binary covariate as an AI/ML-generated noisy label for an unobserved latent treatment. The latent treatment is integrated out via a 2 × 2 misclassification table; the residual distribution is selected by `dist()` (normal, Laplace, or Student-t) with separately-parameterized scales in the two latent classes (heteroskedastic by default). Use this when no external estimate of the false-positive rate is available. Mirrors Python `ValidMLInference.one_step`.

The optimizer is Mata's `optimize()` with BFGS on the negative log-likelihood. The likelihood is computed entirely in log-space (logsumexp throughout) to stay finite when the optimizer probes parameters that would otherwise drive `exp(log sigma)` past the IEEE double overflow threshold. Variance is the top-`d` block of `pinv(H)`, where `H` is the numerical Hessian at the optimum.

**Options:**
- `generated(varname)` — name of the binary AI/ML-generated covariate. Defaults to the first variable in the varlist. Must take values 0 and 1.
- `homoskedastic` — constrain `sigma_0 = sigma_1`. Default is heteroskedastic.
- `dist(name)` — residual distribution: `normal` (default), `laplace`, or `t`. The scale parameter `sigma` is reinterpreted per distribution (Laplace `b`, Student-`t` scale); the theta layout is unchanged.
- `df(#)` — degrees of freedom for `dist(t)`. Required and must be `> 0` when `dist(t)` is set; rejected with any other distribution. Fixed (not estimated).
- `noconstant` — omit the constant term.

**Posted estimates:** `e(b)`, `e(V)`, `e(N)`, `e(ll)`, `e(converged)`, `e(iterations)`, `e(homoskedastic)`, `e(distcode)`, `e(df)` (only when `dist(t)`), `e(generated_var)`, `e(cmd)`, `e(cmdline)`, `e(depvar)`, `e(dist)`, `e(title)`, `e(vce)`, `e(vcetype)`, `e(properties)`.

**Departures from Python.** Python's `one_step(..., distribution=...)` accepts an arbitrary JAX-traceable PDF callable; Stata exposes the three most common families as a closed menu (no clean Stata analogue for a JAX-traceable callable). Starting values use the Gaussian pdf for the latent-X* imputation step regardless of `dist()`, matching Python's `get_starting_values_unlabeled_jax_with_treatment_idx` (hard-coded Gaussian `pdf_func`).

**Parity numbers (synthetic fixtures, vs Python upstream):**
- Gaussian: max |db| ≈ 1.2e-6, max |dV| ≈ 1.1e-7.
- Laplace: max |db| ≈ 2.9e-4, max |dV| ≈ 2.1e-3.
- Student-t(df=5): max |db| ≈ 2.7e-6, max |dV| ≈ 6.6e-9.

**Example:**
```stata
use "ValidMLInference-stata/data/remote_work_data.dta", clear
gen log_salary = log(salary)
one_step log_salary wfh_rwo, generated(wfh_rwo)
one_step log_salary wfh_rwo, generated(wfh_rwo) dist(laplace)
one_step log_salary wfh_rwo, generated(wfh_rwo) dist(t) df(5)
```

---

## `one_step_gmm`

**Status:** implemented (v0.5.0, 2026-05-22). Stata-side command name shortened from the Python `one_step_gaussian_mixture`.

**Syntax:**
```
one_step_gmm depvar varlist [if] [in] [, generated(varname) k(integer) nguess(integer) maxiter(integer) seed(integer) noconstant]
```

**Description:** Joint MLE using a `k`-component Gaussian mixture for the outcome residual in each latent-treatment class. The mixture can accommodate skewness, fat tails, or bimodality while still identifying the misclassification probabilities. Mirrors Python `ValidMLInference.one_step_gaussian_mixture`.

Multistart: the optimizer runs from `nguess` starts. Start 0 uses the heteroskedastic OLS-based starts directly; subsequent starts perturb them by Gaussian noise (scaled 0.5 on coefficient slots, 0.3 on misclassification logits, 1.0 elsewhere — matching Python's `_one_step_gaussian_mixture_core`). The best attempt by negative log-likelihood wins. The likelihood is computed in log-space throughout (logsumexp for both the inner component mixtures and the outer 2-class mixture), with Python's `+1e-12` numerical stabilizer applied in log-domain.

**Options:**
- `generated(varname)` — name of the binary AI/ML-generated covariate. Defaults to the first variable in the varlist.
- `k(integer)` — number of mixture components (default 2; must be ≥ 2).
- `nguess(integer)` — number of multistart attempts (default 10; must be ≥ 1).
- `maxiter(integer)` — maximum BFGS iterations per attempt (default 100).
- `seed(integer)` — seed for the multistart perturbations (default 0). Stata's RNG is unrelated to JAX's PRNGKey, so per-attempt parity for `i > 0` is not expected; the `i = 0` attempt is deterministic and identical across implementations.
- `homoskedastic` — **rejected at runtime.** Python's `get_starting_values_unlabeled_gaussian_mixture(homosked=True)` returns a theta of length `d + 1 + 3k` but `unpack_theta(homosked=True)` reads `d - 1 + 5k`; supporting it in Stata would require deviating from upstream layout. See `Notes/porting_decisions.md`.
- `noconstant` — omit the constant term.

**Posted estimates:** `e(b)`, `e(V)`, `e(N)`, `e(k)`, `e(nguess)`, `e(maxiter)`, `e(seed)`, `e(ll)`, `e(converged)`, `e(iterations)`, `e(best_idx)`, `e(n_finished)`, `e(homoskedastic)`, `e(generated_var)`, `e(cmd)`, `e(cmdline)`, `e(depvar)`, `e(title)`, `e(vce)`, `e(vcetype)`, `e(properties)`.

**Reduced-dimension optimization at symmetric starts.** For `k = 2`, the OLS-based starts place four mixture-structure parameters at exactly zero, making the two components mathematically identical. The likelihood is then analytically flat in those four directions, and Mata's numerical-derivative BFGS would abort. `one_step_gmm` detects this case and runs the optimizer (and the Hessian step) on the active subspace. This matches Python's effective behavior — jaxopt's exact-zero autodiff gradient on the same subspace keeps those parameters at zero throughout. The coefficient-block variance is unchanged by the dimension reduction because the Hessian is block-diagonal between active and symmetric slots at that point. See `Notes/porting_decisions.md`.

**Example:**
```stata
one_step_gmm y x1 x2, generated(x1) k(2) nguess(10) seed(0)
```

---

## Data loaders

| Python | Stata equivalent |
|---|---|
| `remote_work_data()` | `use "data/remote_work_data.dta", clear` |
| `topic_model_data()` | `do "data/load_topic_model_data.do"` (TBD) |
| `fed_sentiment_data()` | `use "data/fed_sentiment_meetings.dta", clear` + matrix loaders for `B` and `kappa` |

---

## Notes for users coming from the Python package

- **Formula syntax is dropped.** Use Stata varlists (e.g. `i.soc_2021_2`) instead of `C(soc_2021_2)`.
- **`intercept` flag → `noconstant` option.** Stata convention is "constant by default; opt out".
- **Returns live in `e(...)`.** Use `matrix list e(b)` and `matrix list e(V)`, or `estimates store`/`esttab`/`coefplot` as with any native estimator.
- **Coefficient ordering** is standardized to intercept-first, then alphabetical — matching the Python package, not Stata's default ordering.
