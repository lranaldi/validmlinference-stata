"""
refresh_fixtures.py -- generate Python reference fixtures for the Stata port

Runs the upstream ValidMLInference package on a fixed set of deterministic
cases, then writes per-case CSV bundles into ./fixtures/<case>/. The Stata
side reads the bundles back via run_all_tests.do and checks numerical parity.

Layout per case:
    fixtures/<case>/meta.txt        # command, options, tolerance, dims
    fixtures/<case>/input.csv       # observation-level data (header row)
    fixtures/<case>/W.csv,S.csv,B.csv  # (topic cases only) matrices, no header
    fixtures/<case>/ref_coef.csv    # name,coef,se
    fixtures/<case>/ref_vcov.csv    # square matrix, header row + index col

Stata-side conventions: the Python intercept name 'Intercept' is rewritten to
'_cons' here so the fixture row labels match what `ereturn post` writes into
e(b).

Run:
    .conda/envs/vmli/python.exe ValidMLInference-stata/tests/refresh_fixtures.py

Pinned against upstream ValidMLInference 1.4.0.
"""

from __future__ import annotations

import os
import sys
import json
from pathlib import Path

# Enable JAX 64-bit. Upstream ValidMLInference (1.4.0) routes some final-step
# reorderings through jnp.take in _reorder_intercept_first; without this flag,
# JAX silently downcasts to float32 and the public coefficients lose ~7 sig
# digits. We want the fixtures to reflect the *intended* float64 math, not the
# inadvertent precision loss.
import jax
jax.config.update("jax_enable_x64", True)

import numpy as np
import pandas as pd
import jax.numpy as jnp
import jax.scipy.stats as jss

from ValidMLInference import (
    ols,
    ols_bca,
    ols_bcm,
    ols_bca_topic,
    ols_bcm_topic,
    one_step,
    one_step_gaussian_mixture,
    remote_work_data,
)


# ---------------------------------------------------------------------------
# Distribution callables for the `distribution=` argument of one_step /
# one_step_gaussian_mixture. Signature must be pdf(y, loc, scale) returning
# a JAX array, and must be hashable (functions are) because upstream uses
# @partial(jit, static_argnames=('distribution',)).
# ---------------------------------------------------------------------------
def laplace_pdf(y, loc, scale):
    return jss.laplace.pdf(y, loc=loc, scale=scale)


def _make_t_pdf(df):
    def t_pdf(y, loc, scale):
        return jss.t.pdf(y, df=df, loc=loc, scale=scale)
    return t_pdf


t5_pdf = _make_t_pdf(5.0)

try:
    from importlib.metadata import version as _pkg_version
    VMLI_VERSION = _pkg_version("ValidMLInference")
except Exception:
    VMLI_VERSION = "1.4.0"

THIS = Path(__file__).resolve().parent
OUT_ROOT = THIS / "fixtures"
OUT_ROOT.mkdir(exist_ok=True)


def _rename_intercept(names):
    return ["_cons" if n in ("Intercept", "(Intercept)") else n for n in names]


def _write_case(case, meta, input_df, ref_coef, ref_vcov, matrices=None):
    out = OUT_ROOT / case
    out.mkdir(exist_ok=True)

    # Stable meta: KEY=VALUE per line. Stata reads with `file read`.
    with open(out / "meta.txt", "w") as f:
        for k, v in meta.items():
            f.write(f"{k}={v}\n")

    fmt = "%.17g"
    input_df.to_csv(out / "input.csv", index=False, float_format=fmt)
    ref_coef.to_csv(out / "ref_coef.csv", index=False, float_format=fmt)

    # ref_vcov in long format so the Stata side can import without colliding
    # with leading-underscore names like "_cons".
    long_rows = []
    for i, ni in enumerate(ref_vcov.index):
        for j, nj in enumerate(ref_vcov.columns):
            long_rows.append({"row": ni, "col": nj,
                              "vcov": ref_vcov.iat[i, j]})
    pd.DataFrame(long_rows, columns=["row", "col", "vcov"]).to_csv(
        out / "ref_vcov.csv", index=False, float_format=fmt)

    if matrices:
        for name, M in matrices.items():
            pd.DataFrame(M).to_csv(out / f"{name}.csv", index=False,
                                   header=False, float_format=fmt)


def _result_to_tables(result):
    """Convert a RegressionResult into (ref_coef DF, ref_vcov DF) with
    intercept renamed to _cons."""
    names = _rename_intercept(list(result.names))
    b = np.asarray(result.coef).ravel()
    V = np.asarray(result.vcov)
    se = np.sqrt(np.diag(V))

    ref_coef = pd.DataFrame({"name": names, "coef": b, "se": se})
    ref_vcov = pd.DataFrame(V, index=names, columns=names)
    return ref_coef, ref_vcov


# ---------------------------------------------------------------------------
# Case builders
# ---------------------------------------------------------------------------

def case_ols_synthetic():
    rng = np.random.default_rng(20260518)
    n = 200
    X = rng.standard_normal((n, 3))
    eps = rng.standard_normal(n)
    beta = np.array([1.0, -0.5, 0.25])
    Y = X @ beta + 0.4 + eps

    names = ["x1", "x2", "x3"]
    res = ols(Y=Y, X=X, names=names, intercept=True)
    ref_coef, ref_vcov = _result_to_tables(res)

    df = pd.DataFrame({"y": Y, "x1": X[:, 0], "x2": X[:, 1], "x3": X[:, 2]})
    meta = {
        "case": "ols_synthetic",
        "command": "ols",
        "depvar": "y",
        "covars": "x1 x2 x3",
        "options": "",
        "tol_coef": 1e-8,
        "tol_vcov": 1e-7,
        "n": n,
        "k": len(names) + 1,
        "vmli_version": "1.4.0",
        "note": "Pure OLS on 3 standard-normal regressors + intercept; deterministic seed 20260518.",
    }
    _write_case("ols_synthetic", meta, df, ref_coef, ref_vcov)


def case_ols_remote_work():
    df_full = remote_work_data().copy()
    df_full = df_full[df_full["salary"] > 0].reset_index(drop=True)
    df_full["log_salary"] = np.log(df_full["salary"])
    Y = df_full["log_salary"].values
    X = df_full[["wfh_rwo"]].values.astype(float)

    res = ols(Y=Y, X=X, names=["wfh_rwo"], intercept=True)
    ref_coef, ref_vcov = _result_to_tables(res)

    df = pd.DataFrame({"log_salary": Y, "wfh_rwo": X[:, 0]})
    meta = {
        "case": "ols_remote_work",
        "command": "ols",
        "depvar": "log_salary",
        "covars": "wfh_rwo",
        "options": "",
        "tol_coef": 1e-8,
        "tol_vcov": 1e-7,
        "n": len(df),
        "k": 2,
        "vmli_version": "1.4.0",
        "note": "Univariate OLS on the shipped remote_work_data after dropping non-positive salaries.",
    }
    _write_case("ols_remote_work", meta, df, ref_coef, ref_vcov)


def _synthetic_binary_data(seed=20260518, n=500):
    rng = np.random.default_rng(seed)
    x1 = (rng.uniform(size=n) < 0.4).astype(float)  # binary label
    x2 = rng.standard_normal(n)
    eps = rng.standard_normal(n) * 0.5
    Y = 0.7 * x1 - 0.3 * x2 + 0.2 + eps
    return Y, x1, x2


def case_ols_bca_synthetic():
    Y, x1, x2 = _synthetic_binary_data()
    Xhat = np.column_stack([x1, x2])
    res = ols_bca(
        Y=Y,
        Xhat=Xhat,
        names=["x1", "x2"],
        fpr=0.05,
        m=2000,
        generated_var="x1",
        intercept=True,
    )
    ref_coef, ref_vcov = _result_to_tables(res)
    df = pd.DataFrame({"y": Y, "x1": x1, "x2": x2})
    meta = {
        "case": "ols_bca_synthetic",
        "command": "ols_bca",
        "depvar": "y",
        "covars": "x1 x2",
        "options": "fpr(0.05) m(2000) generated(x1)",
        "tol_coef": 1e-8,
        "tol_vcov": 1e-7,
        "n": len(Y),
        "k": 3,
        "vmli_version": "1.4.0",
        "note": "Additive bias correction on a synthetic binary regressor; fpr=0.05, m=2000.",
    }
    _write_case("ols_bca_synthetic", meta, df, ref_coef, ref_vcov)


def case_ols_bcm_synthetic():
    Y, x1, x2 = _synthetic_binary_data()
    Xhat = np.column_stack([x1, x2])
    res = ols_bcm(
        Y=Y,
        Xhat=Xhat,
        names=["x1", "x2"],
        fpr=0.05,
        m=2000,
        generated_var="x1",
        intercept=True,
    )
    ref_coef, ref_vcov = _result_to_tables(res)
    df = pd.DataFrame({"y": Y, "x1": x1, "x2": x2})
    meta = {
        "case": "ols_bcm_synthetic",
        "command": "ols_bcm",
        "depvar": "y",
        "covars": "x1 x2",
        "options": "fpr(0.05) m(2000) generated(x1)",
        "tol_coef": 1e-8,
        "tol_vcov": 1e-7,
        "n": len(Y),
        "k": 3,
        "vmli_version": "1.4.0",
        "note": "Multiplicative bias correction on a synthetic binary regressor; fpr=0.05, m=2000.",
    }
    _write_case("ols_bcm_synthetic", meta, df, ref_coef, ref_vcov)


def _synthetic_topic_data(seed=20260518, n=300, r=3, c=5, v=8, q=2):
    rng = np.random.default_rng(seed)
    W = rng.uniform(size=(n, c))
    W = W / W.sum(axis=1, keepdims=True)
    S = rng.uniform(size=(r, c))
    S = S / S.sum(axis=1, keepdims=True)
    B = rng.uniform(size=(c, v))
    B = B / B.sum(axis=1, keepdims=True)
    Theta = W @ S.T
    Q = rng.standard_normal((n, q))
    bt = np.array([1.0, 0.5, -0.4])
    bq = np.array([0.2, -0.3])
    eps = rng.standard_normal(n)
    Y = Theta @ bt + Q @ bq + 0.7 + eps
    return Y, Q, W, S, B


def _topic_case(case_name, fit_fn, k_val, note):
    Y, Q, W, S, B = _synthetic_topic_data()
    names = [f"Q_{i+1}" for i in range(Q.shape[1])]
    # Pass user-friendly names; topic_1..topic_r are filled in by the impl.
    res = fit_fn(Y=Y, Q=Q, W=W, S=S, B=B, k=k_val, intercept=True, names=None)
    # Result names from upstream: topic_1..topic_r, Q_1..Q_q, Intercept
    # We rename "Q_1"->"q1", "Q_2"->"q2" to match what the Stata test will
    # pass in via varlist.
    out_names = []
    for n in res.names:
        if n.startswith("Q_"):
            out_names.append("q" + n[2:])
        else:
            out_names.append(n)
    res.names = out_names
    ref_coef, ref_vcov = _result_to_tables(res)

    # input_df only holds y, q1, q2 (W/S/B go to separate matrix CSVs)
    n = len(Y)
    df_cols = {"y": Y}
    for i in range(Q.shape[1]):
        df_cols[f"q{i+1}"] = Q[:, i]
    df = pd.DataFrame(df_cols)

    cmd = {
        "bca": "ols_bca_topic",
        "bcm": "ols_bcm_topic",
    }["bcm" if "bcm" in case_name else "bca"]

    meta = {
        "case": case_name,
        "command": cmd,
        "depvar": "y",
        "covars": "q1 q2",
        "options": f"wmatrix(W) smatrix(S) bmatrix(B) k({k_val})",
        "tol_coef": 1e-8,
        "tol_vcov": 1e-7,
        "n": n,
        "k": 6,
        "r": S.shape[0],
        "c": S.shape[1],
        "v": B.shape[1],
        "q": Q.shape[1],
        "kval": k_val,
        "vmli_version": "1.4.0",
        "note": note,
    }
    _write_case(case_name, meta, df, ref_coef, ref_vcov,
                matrices={"W": W, "S": S, "B": B})


def case_ols_bca_topic():
    _topic_case(
        "ols_bca_topic_synthetic",
        ols_bca_topic,
        k_val=0.5,
        note="Additive topic-model correction; small synthetic recipe (n=300, r=3, c=5, v=8, q=2); k=0.5.",
    )


def case_ols_bcm_topic_small_k():
    _topic_case(
        "ols_bcm_topic_small_k",
        ols_bcm_topic,
        k_val=0.005,
        note="Multiplicative topic-model correction with small k=0.005 -> rho<1 branch.",
    )


def case_ols_bcm_topic_large_k():
    _topic_case(
        "ols_bcm_topic_large_k",
        ols_bcm_topic,
        k_val=20.0,
        note="Multiplicative topic-model correction with large k=20 -> rho>=1 BCA fallback.",
    )


def case_one_step_synthetic():
    rng = np.random.default_rng(20260518)
    n = 2000

    # Latent true binary X*
    x_star = (rng.uniform(size=n) < 0.4).astype(float)

    # Observed/noisy x1: flip ~10% of labels (symmetric misclassification)
    flip = (rng.uniform(size=n) < 0.10)
    x1   = np.where(flip, 1.0 - x_star, x_star).astype(float)

    # Auxiliary continuous regressor
    x2 = rng.standard_normal(n)

    # Heteroskedastic Gaussian errors keyed off the latent treatment
    sigma0, sigma1 = 0.5, 0.8
    eps = np.where(
        x_star == 1,
        rng.standard_normal(n) * sigma1,
        rng.standard_normal(n) * sigma0,
    )
    Y = 0.2 + 0.7 * x_star - 0.3 * x2 + eps

    Xhat = np.column_stack([x1, x2])
    res = one_step(
        Y=Y,
        Xhat=Xhat,
        names=["x1", "x2"],
        generated_var="x1",
        homoskedastic=False,
        intercept=True,
    )
    ref_coef, ref_vcov = _result_to_tables(res)

    df = pd.DataFrame({"y": Y, "x1": x1, "x2": x2})
    meta = {
        "case": "one_step_synthetic",
        "command": "one_step",
        "depvar": "y",
        "covars": "x1 x2",
        "options": "generated(x1)",
        "tol_coef": 1e-3,
        "tol_vcov": 1e-2,
        "n": n,
        "k": 3,
        "vmli_version": "1.4.0",
        "note": "One-step MLE on a binary AI/ML-generated regressor; "
                "true fpr=fnr=0.10, heteroskedastic Gaussian errors "
                "(sigma0=0.5, sigma1=0.8). MLE tolerance reflects optimizer "
                "differences between JAX-LBFGS and Mata-BFGS.",
    }
    _write_case("one_step_synthetic", meta, df, ref_coef, ref_vcov)


def case_one_step_gmm_synthetic():
    """One-step GMM (k=2) MLE on a binary AI-generated regressor with
    bimodal residuals.

    Design notes:
      * nguess is pinned to 1 for the parity fixture. JAX's PRNGKey stream
        and Mata's `rseed`+`rnormal` stream are unrelated, so the i>0
        multistart attempts can't be compared draw-by-draw. The i=0
        attempt -- OLS-based starts, no perturbation -- is identical across
        implementations, and for a well-identified mixture problem it
        reaches the global optimum on its own. This fixture exercises that
        deterministic path; multistart correctness is checked separately.
      * Residuals are a centered 2-component Gaussian mixture in each
        latent-treatment class, so k=2 is genuinely identified.
    """
    rng = np.random.default_rng(20260522)
    n = 2000

    x_star = (rng.uniform(size=n) < 0.4).astype(float)
    flip   = (rng.uniform(size=n) < 0.10)
    x1     = np.where(flip, 1.0 - x_star, x_star).astype(float)
    x2     = rng.standard_normal(n)

    # Centered bimodal residuals (weighted mean 0 in each class) so the
    # mixture is identified relative to the regression intercept.
    eps = np.empty(n)
    z   = rng.standard_normal(n)
    u   = rng.uniform(size=n)

    # Class 1 (x_star == 1): 0.6 * N(-0.4, 0.3) + 0.4 * N(+0.6, 0.5)
    in1 = (x_star == 1)
    c1  = u[in1] < 0.6
    eps[in1] = np.where(c1, -0.4 + 0.3 * z[in1], 0.6 + 0.5 * z[in1])

    # Class 0 (x_star == 0): 0.7 * N(-0.3, 0.2) + 0.3 * N(+0.7, 0.4)
    in0 = (x_star == 0)
    c0  = u[in0] < 0.7
    eps[in0] = np.where(c0, -0.3 + 0.2 * z[in0], 0.7 + 0.4 * z[in0])

    Y = 0.2 + 0.7 * x_star - 0.3 * x2 + eps

    Xhat = np.column_stack([x1, x2])
    res = one_step_gaussian_mixture(
        Y=Y,
        Xhat=Xhat,
        names=["x1", "x2"],
        generated_var="x1",
        k=2,
        homosked=False,
        nguess=1,
        maxiter=500,
        seed=0,
        intercept=True,
    )
    ref_coef, ref_vcov = _result_to_tables(res)

    df = pd.DataFrame({"y": Y, "x1": x1, "x2": x2})
    meta = {
        "case": "one_step_gmm_synthetic",
        "command": "one_step_gmm",
        "depvar": "y",
        "covars": "x1 x2",
        "options": "generated(x1) k(2) nguess(1) maxiter(500) seed(0)",
        "tol_coef": 1e-3,
        "tol_vcov": 1e-2,
        "n": n,
        "k": 3,
        "vmli_version": "1.4.0",
        "note": "One-step Gaussian-mixture MLE with k=2 components, "
                "heteroskedastic across the latent treatment. "
                "Bimodal centered residuals in each class. "
                "nguess=1 to fix the deterministic i=0 path for "
                "parity (multistart RNG differs between JAX and Mata).",
    }
    _write_case("one_step_gmm_synthetic", meta, df, ref_coef, ref_vcov)


# ---------------------------------------------------------------------------
# Non-Gaussian one_step / one_step_gmm fixtures
# ---------------------------------------------------------------------------

def _gen_one_step_data(seed, n, dist, df=None, b0=0.5, b1=0.8):
    """Synthetic one_step data with class-heteroskedastic residuals drawn
    from `dist` (normal | laplace | t). x_star is the latent treatment,
    x1 is the observed label with ~10% symmetric flips."""
    rng = np.random.default_rng(seed)
    x_star = (rng.uniform(size=n) < 0.4).astype(float)
    flip   = (rng.uniform(size=n) < 0.10)
    x1     = np.where(flip, 1.0 - x_star, x_star).astype(float)
    x2     = rng.standard_normal(n)

    scale_per_obs = np.where(x_star == 1, b1, b0)
    if dist == "normal":
        eps = rng.standard_normal(n) * scale_per_obs
    elif dist == "laplace":
        eps = rng.laplace(loc=0.0, scale=1.0, size=n) * scale_per_obs
    elif dist == "t":
        eps = rng.standard_t(df, size=n) * scale_per_obs
    else:
        raise ValueError(dist)
    Y = 0.2 + 0.7 * x_star - 0.3 * x2 + eps
    return Y, x1, x2


def _case_one_step_dist(name, dist, jax_pdf, df=None, opt_extra=""):
    Y, x1, x2 = _gen_one_step_data(seed=20260518, n=2000, dist=dist, df=df)
    Xhat = np.column_stack([x1, x2])
    res = one_step(
        Y=Y,
        Xhat=Xhat,
        names=["x1", "x2"],
        generated_var="x1",
        homoskedastic=False,
        intercept=True,
        distribution=jax_pdf,
    )
    ref_coef, ref_vcov = _result_to_tables(res)
    df_in = pd.DataFrame({"y": Y, "x1": x1, "x2": x2})

    options = f"generated(x1) dist({dist})"
    if dist == "t":
        options += f" df({int(df) if float(df).is_integer() else df})"
    if opt_extra:
        options += " " + opt_extra

    meta = {
        "case": name,
        "command": "one_step",
        "depvar": "y",
        "covars": "x1 x2",
        "options": options,
        "tol_coef": 1e-3,
        "tol_vcov": 1e-2,
        "n": len(Y),
        "k": 3,
        "vmli_version": "1.4.0",
        "note": f"One-step MLE with dist={dist}"
                + (f" df={df}" if dist == "t" else "")
                + "; true fpr=fnr=0.10, heteroskedastic residuals "
                  "(b0=0.5, b1=0.8). Tolerance accounts for "
                  "JAX-LBFGS vs Mata-BFGS optimizer differences.",
    }
    _write_case(name, meta, df_in, ref_coef, ref_vcov)


def case_one_step_synthetic_laplace():
    _case_one_step_dist("one_step_synthetic_laplace", "laplace", laplace_pdf)


def case_one_step_synthetic_t5():
    _case_one_step_dist("one_step_synthetic_t5", "t", t5_pdf, df=5.0)


CASES = [
    case_ols_synthetic,
    case_ols_remote_work,
    case_ols_bca_synthetic,
    case_ols_bcm_synthetic,
    case_ols_bca_topic,
    case_ols_bcm_topic_small_k,
    case_ols_bcm_topic_large_k,
    case_one_step_synthetic,
    case_one_step_gmm_synthetic,
    case_one_step_synthetic_laplace,
    case_one_step_synthetic_t5,
]


def main():
    print(f"ValidMLInference version: {VMLI_VERSION}")
    print(f"Writing fixtures into: {OUT_ROOT}")
    for fn in CASES:
        name = fn.__name__.replace("case_", "")
        print(f"  -> {name}")
        fn()
    print("done.")


if __name__ == "__main__":
    main()
