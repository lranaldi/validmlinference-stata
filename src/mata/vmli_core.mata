*! version 0.5.0  22may2026  Lorenzo Ranaldi
*! ValidMLInference (Stata port) -- Mata core
*!
*! Provides:
*!   - vmli_version()              package version string (smoke test)
*!   - vmli_alpha_perm()           intercept-first / alphabetical permutation
*!   - vmli_ols_fit()              OLS with HC0 variance
*!   - vmli_ols_run()              driver for ols.ado
*!   - vmli_bc_prepare()           shared b0/V0/Gamma setup for binary-label BC
*!   - vmli_ols_bca_fit()          additive bias correction (binary label)
*!   - vmli_ols_bcm_fit()          multiplicative bias correction (binary label)
*!   - vmli_ols_bca_run()          driver for ols_bca.ado
*!   - vmli_ols_bcm_run()          driver for ols_bcm.ado
*!   - vmli_bc_topic_prepare()     shared setup for topic-model BC
*!   - vmli_ols_bca_topic_fit()    additive bias correction (topic shares)
*!   - vmli_ols_bcm_topic_fit()    multiplicative bias correction (topic shares)
*!   - vmli_ols_bca_topic_run()    driver for ols_bca_topic.ado
*!   - vmli_ols_bcm_topic_run()    driver for ols_bcm_topic.ado
*!   - vmli_normal_pdf()           Gaussian pdf, vectorized
*!   - vmli_subset_std()           population std on a 0/1-masked subset
*!   - _vmli_log_pdf()             log-pdf dispatch: normal | laplace | t
*!   - vmli_theta_to_pars()        unpack one_step parameter vector
*!   - vmli_one_step_starts()      starting values for one_step MLE
*!   - vmli_one_step_obj()         optimize() d0 evaluator (neg log-likelihood)
*!   - vmli_one_step_fit()         driver: starts -> optimize -> Hessian -> V
*!   - vmli_one_step_run()         Stata-facing driver for one_step.ado
*!   - vmli_log_softmax_ref()      log-softmax with implicit 0 reference
*!   - vmli_centered_means()       cumsum->append-0->subtract weighted mean
*!   - vmli_log_mixture_pdf()      log-density of a Gaussian mixture
*!   - vmli_gmm_starts()           starting values for one_step_gmm MLE
*!   - vmli_gmm_obj()              optimize() d0 evaluator for mixture MLE
*!   - vmli_gmm_fit()              multistart driver for one_step_gmm
*!   - vmli_gmm_run()              Stata-facing driver for one_step_gmm.ado
*!   - _vmli_gmm_symmetric_slots() private: detect flat-gradient slots
*!   - _vmli_gmm_obj_reduced()     private: reduced-dim d0 evaluator wrapper
*!   - _vmli_post_results()        private: alpha-perm + push to Stata matrices
*!   - _vmli_find_col()            private: locate a name in a string rowvector

version 17

mata:

mata set matastrict on

// ---------------------------------------------------------------------------
// vmli_version
//   Returns the package version. A successful call confirms lvmli.mlib
//   compiled and is on Mata's library path.
// ---------------------------------------------------------------------------
string scalar vmli_version()
{
    return("0.5.0")
}

// ---------------------------------------------------------------------------
// vmli_alpha_perm(names)
//   Permutation putting the intercept (if any) first, then sorting the
//   remaining names alphabetically. Mirrors Python's
//   _standardize_coefficient_order. Recognized intercept labels:
//   "_cons", "Intercept", "(Intercept)".
// ---------------------------------------------------------------------------
real rowvector vmli_alpha_perm(string rowvector names)
{
    real scalar      d, i, cons_idx
    real rowvector   all_idx, non_cons_idx
    real colvector   sort_idx

    d       = cols(names)
    all_idx = 1..d

    cons_idx = 0
    for (i = 1; i <= d; i++) {
        if (names[i] == "_cons" |
            names[i] == "Intercept" |
            names[i] == "(Intercept)") {
            cons_idx = i
            break
        }
    }

    if (cons_idx == 0) {
        sort_idx = order(names', 1)
        return(sort_idx')
    }

    non_cons_idx = select(all_idx, all_idx :!= cons_idx)
    sort_idx     = order(names[non_cons_idx]', 1)
    return((cons_idx, non_cons_idx[sort_idx']))
}

// ---------------------------------------------------------------------------
// vmli_ols_result struct
//   Container for OLS results so estimators that build on OLS
//   (ols_bca, ols_bcm, topic variants) can share the same shape.
// ---------------------------------------------------------------------------
struct vmli_ols_result {
    real colvector  b
    real matrix     V
    real matrix     XX_inv     // (X'X)^{-1}: needed by bias-correction stages
    real scalar     n
}

// ---------------------------------------------------------------------------
// vmli_ols_fit(Y, X)
//   OLS coefficients and HC0 (no df adjustment) robust variance:
//      b      = (X'X)^{-1} X'Y
//      Omega  = sum_i u_i^2 x_i x_i'   (=  X' diag(u^2) X)
//      V      = (X'X)^{-1} Omega (X'X)^{-1}
//   Matches Python _ols_core exactly.
// ---------------------------------------------------------------------------
struct vmli_ols_result scalar vmli_ols_fit(real colvector Y, real matrix X)
{
    struct vmli_ols_result scalar  r
    real matrix    Omega
    real colvector u

    r.n      = rows(X)
    r.XX_inv = invsym(quadcross(X, X))
    r.b      = r.XX_inv * quadcross(X, Y)

    u        = Y - X * r.b
    Omega    = quadcross(X, u :^ 2, X)
    r.V      = r.XX_inv * Omega * r.XX_inv

    return(r)
}

// ---------------------------------------------------------------------------
// vmli_ols_run(yname, xnames, tousename, addconst, bname, Vname)
//   Stata-facing driver. Reads the data via st_data(), fits OLS with HC0
//   variance, reorders coefficients (intercept first, then alphabetical),
//   and writes results into the Stata matrices `bname` and `Vname` with
//   column/row stripes set. Sample size goes into r(N).
// ---------------------------------------------------------------------------
void vmli_ols_run(
    string scalar  yname,
    string scalar  xnames,
    string scalar  tousename,
    real scalar    addconst,
    string scalar  bname,
    string scalar  Vname)
{
    real colvector  Y
    real matrix     X
    string rowvector names
    struct vmli_ols_result scalar r

    Y     = st_data(., yname, tousename)
    names = tokens(xnames)
    X     = st_data(., names, tousename)

    if (addconst) {
        X     = X, J(rows(X), 1, 1)
        names = names, "_cons"
    }

    r = vmli_ols_fit(Y, X)

    _vmli_post_results(r.b, r.V, names, r.n, bname, Vname)
}

// ===========================================================================
// Closed-form bias-correction estimators (Milestone 3)
// ===========================================================================

// ---------------------------------------------------------------------------
// vmli_bc_setup struct
//   Shared state for closed-form bias-correction estimators (ols_bca,
//   ols_bcm). Holds the uncorrected OLS fit plus the Gamma matrix that
//   the correction stages multiply through.
// ---------------------------------------------------------------------------
struct vmli_bc_setup {
    real colvector  b0
    real matrix     V0
    real matrix     Gamma
    real scalar     n
    real scalar     d
}

// ---------------------------------------------------------------------------
// vmli_bc_prepare(Y, X, target_idx)
//   Runs OLS (HC0) and builds Gamma = n * (X'X)^{-1} * A, where A is the
//   d x d matrix with 1 in position (target_idx, target_idx) and 0
//   elsewhere. Matches Python's setup in _ols_bca_core / _ols_bcm_core:
//     sXX   = (1/n) X'X
//     Gamma = solve(sXX, A) = n * (X'X)^{-1} * A.
// ---------------------------------------------------------------------------
struct vmli_bc_setup scalar vmli_bc_prepare(real colvector Y,
                                            real matrix X,
                                            real scalar target_idx)
{
    struct vmli_bc_setup scalar  s
    struct vmli_ols_result scalar r
    real matrix A

    r    = vmli_ols_fit(Y, X)
    s.b0 = r.b
    s.V0 = r.V
    s.n  = r.n
    s.d  = cols(X)

    A = J(s.d, s.d, 0)
    A[target_idx, target_idx] = 1
    s.Gamma = s.n * r.XX_inv * A

    return(s)
}

// ---------------------------------------------------------------------------
// vmli_bc_result struct
//   Output of a bias-corrected estimator.
// ---------------------------------------------------------------------------
struct vmli_bc_result {
    real colvector  b
    real matrix     V
}

// ---------------------------------------------------------------------------
// vmli_ols_bca_fit(Y, X, target_idx, fpr, m)
//   Additive bias correction (Python _ols_bca_core):
//     b_corr = b0 + fpr * Gamma * b0
//     M      = I + fpr * Gamma
//     V_corr = M V0 M'
//            + (fpr (1-fpr) / m) * Gamma (V0 + b_corr b_corr') Gamma'
// ---------------------------------------------------------------------------
struct vmli_bc_result scalar vmli_ols_bca_fit(real colvector Y,
                                              real matrix    X,
                                              real scalar    target_idx,
                                              real scalar    fpr,
                                              real scalar    m)
{
    struct vmli_bc_setup  scalar s
    struct vmli_bc_result scalar bc
    real matrix I_d, M, extra

    s = vmli_bc_prepare(Y, X, target_idx)

    bc.b  = s.b0 + fpr * (s.Gamma * s.b0)
    I_d   = I(s.d)
    M     = I_d + fpr * s.Gamma
    extra = (fpr * (1 - fpr) / m) *
            (s.Gamma * (s.V0 + bc.b * bc.b') * s.Gamma')
    bc.V  = M * s.V0 * M' + extra

    return(bc)
}

// ---------------------------------------------------------------------------
// vmli_ols_bcm_fit(Y, X, target_idx, fpr, m)
//   Multiplicative bias correction (Python _ols_bcm_core):
//     M      = I - fpr * Gamma
//     b_corr = M^{-1} b0
//     V_corr = M^{-1} V0 (M^{-1})'
//            + (fpr (1-fpr) / m) * Gamma (V0 + b_corr b_corr') Gamma'
// ---------------------------------------------------------------------------
struct vmli_bc_result scalar vmli_ols_bcm_fit(real colvector Y,
                                              real matrix    X,
                                              real scalar    target_idx,
                                              real scalar    fpr,
                                              real scalar    m)
{
    struct vmli_bc_setup  scalar s
    struct vmli_bc_result scalar bc
    real matrix I_d, M, M_inv, extra

    s     = vmli_bc_prepare(Y, X, target_idx)

    I_d   = I(s.d)
    M     = I_d - fpr * s.Gamma
    M_inv = luinv(M)
    bc.b  = M_inv * s.b0
    extra = (fpr * (1 - fpr) / m) *
            (s.Gamma * (s.V0 + bc.b * bc.b') * s.Gamma')
    bc.V  = M_inv * s.V0 * M_inv' + extra

    return(bc)
}

// ---------------------------------------------------------------------------
// _vmli_find_col(names, target)
//   Private helper: locate `target` in the string rowvector `names`.
//   Aborts with an error if not found, signaling rc=111 to Stata.
// ---------------------------------------------------------------------------
real scalar _vmli_find_col(string rowvector names, string scalar target)
{
    real scalar i

    for (i = 1; i <= cols(names); i++) {
        if (names[i] == target) return(i)
    }

    errprintf("variable %s not found in design matrix\n", target)
    exit(111)
}

// ---------------------------------------------------------------------------
// _vmli_post_results(b, V, names, n, bname, Vname)
//   Private helper used by every estimator driver: apply the intercept-first
//   / alphabetical permutation to (b, V, names), then push the matrices to
//   Stata with column/row stripes set. Also writes r(N).
// ---------------------------------------------------------------------------
void _vmli_post_results(real colvector  b,
                        real matrix     V,
                        string rowvector names,
                        real scalar     n,
                        string scalar   bname,
                        string scalar   Vname)
{
    real rowvector perm
    string matrix  cstripe

    perm  = vmli_alpha_perm(names)
    b     = b[perm]
    V     = V[perm, perm]
    names = names[perm]

    st_matrix(bname, b')
    st_matrix(Vname, V)

    cstripe = J(cols(names), 1, ""), names'
    st_matrixcolstripe(bname, cstripe)
    st_matrixcolstripe(Vname, cstripe)
    st_matrixrowstripe(Vname, cstripe)

    st_numscalar("r(N)", n)
}

// ---------------------------------------------------------------------------
// vmli_ols_bca_run(yname, xnames, tousename, addconst, fpr, m, gen_var, bname, Vname)
//   Stata-facing driver for ols_bca: reads data, builds the design matrix
//   (appending _cons if requested), locates the target column by name,
//   computes the additive correction, and posts results.
// ---------------------------------------------------------------------------
void vmli_ols_bca_run(string scalar  yname,
                      string scalar  xnames,
                      string scalar  tousename,
                      real scalar    addconst,
                      real scalar    fpr,
                      real scalar    m,
                      string scalar  gen_var,
                      string scalar  bname,
                      string scalar  Vname)
{
    real colvector  Y
    real matrix     X
    string rowvector names
    real scalar     target_idx
    struct vmli_bc_result scalar bc

    Y     = st_data(., yname, tousename)
    names = tokens(xnames)
    X     = st_data(., names, tousename)

    if (addconst) {
        X     = X, J(rows(X), 1, 1)
        names = names, "_cons"
    }

    target_idx = _vmli_find_col(names, gen_var)
    bc         = vmli_ols_bca_fit(Y, X, target_idx, fpr, m)

    _vmli_post_results(bc.b, bc.V, names, rows(X), bname, Vname)
}

// ---------------------------------------------------------------------------
// vmli_ols_bcm_run -- as above, but multiplicative correction.
// ---------------------------------------------------------------------------
void vmli_ols_bcm_run(string scalar  yname,
                      string scalar  xnames,
                      string scalar  tousename,
                      real scalar    addconst,
                      real scalar    fpr,
                      real scalar    m,
                      string scalar  gen_var,
                      string scalar  bname,
                      string scalar  Vname)
{
    real colvector  Y
    real matrix     X
    string rowvector names
    real scalar     target_idx
    struct vmli_bc_result scalar bc

    Y     = st_data(., yname, tousename)
    names = tokens(xnames)
    X     = st_data(., names, tousename)

    if (addconst) {
        X     = X, J(rows(X), 1, 1)
        names = names, "_cons"
    }

    target_idx = _vmli_find_col(names, gen_var)
    bc         = vmli_ols_bcm_fit(Y, X, target_idx, fpr, m)

    _vmli_post_results(bc.b, bc.V, names, rows(X), bname, Vname)
}

// ===========================================================================
// Closed-form bias corrections for topic-model regressors (Milestone 4)
// ===========================================================================

// ---------------------------------------------------------------------------
// vmli_bc_topic_setup struct
//   Holds the uncorrected OLS fit on Xhat = [Theta, Q] plus the Gamma
//   matrix derived from the topic-model bias term Omega.
//   Mirrors Python ols_bc_topic_internal.
// ---------------------------------------------------------------------------
struct vmli_bc_topic_setup {
    real colvector  b0
    real matrix     V0
    real matrix     Gamma
    real scalar     n
    real scalar     d        // total columns of Xhat = r + cols(Q)
    real scalar     r        // number of topics
}

// ---------------------------------------------------------------------------
// vmli_bc_topic_prepare(Y, Q, W, S, B, k)
//   Builds Theta = W S' (n x r topic-share matrix), stacks with Q to form
//   Xhat, runs HC0 OLS, then constructs Omega and Gamma:
//
//     mW    = colMeans(W)              (c x 1)
//     M     = B' .* (B' mW)             (v x c)   broadcasts (v,1) over cols
//     Omega = S inv(BB') B M inv(BB') S' - Theta'Theta / n   (r x r)
//     A     = (d x d) with Omega in the top-left r x r block, 0 elsewhere
//     Gamma = (k / sqrt(n)) * solve(sXX, A) = k * sqrt(n) * (X'X)^{-1} A
//
//   Shape contract:
//     Y : n x 1     Q : n x q  (may be q=0)
//     W : n x c     S : r x c     B : c x v
// ---------------------------------------------------------------------------
struct vmli_bc_topic_setup scalar vmli_bc_topic_prepare(
    real colvector Y,
    real matrix    Q,
    real matrix    W,
    real matrix    S,
    real matrix    B,
    real scalar    k)
{
    struct vmli_bc_topic_setup scalar  s
    struct vmli_ols_result     scalar  r_ols
    real matrix    Theta, Xhat, A, Bt, BBt, BBt_inv, Mmat, Omega
    real colvector mW

    Theta = W * S'
    s.r   = rows(S)

    if (cols(Q) > 0) {
        Xhat = (Theta, Q)
    }
    else {
        Xhat = Theta
    }
    s.d = cols(Xhat)

    r_ols = vmli_ols_fit(Y, Xhat)
    s.b0  = r_ols.b
    s.V0  = r_ols.V
    s.n   = r_ols.n

    mW      = mean(W)'                              // c x 1
    Bt      = B'                                    // v x c
    BBt     = B * Bt                                // c x c
    BBt_inv = invsym(BBt)
    Mmat    = Bt :* (Bt * mW)                       // v x c (broadcast)
    Omega   = S * BBt_inv * B * Mmat * BBt_inv * S' ///
              - quadcross(Theta, Theta) / s.n

    A = J(s.d, s.d, 0)
    A[|1, 1 \ s.r, s.r|] = Omega

    s.Gamma = (k * sqrt(s.n)) * r_ols.XX_inv * A

    return(s)
}

// ---------------------------------------------------------------------------
// vmli_ols_bca_topic_fit
//   Additive topic-model bias correction (Python ols_bca_topic):
//     b_corr = (I + Gamma) * b0
//     V_corr = V0           (no extra SE inflation in the Python source)
// ---------------------------------------------------------------------------
struct vmli_bc_result scalar vmli_ols_bca_topic_fit(
    real colvector Y,
    real matrix    Q,
    real matrix    W,
    real matrix    S,
    real matrix    B,
    real scalar    k)
{
    struct vmli_bc_topic_setup scalar s
    struct vmli_bc_result      scalar bc

    s    = vmli_bc_topic_prepare(Y, Q, W, S, B, k)
    bc.b = (I(s.d) + s.Gamma) * s.b0
    bc.V = s.V0

    return(bc)
}

// ---------------------------------------------------------------------------
// vmli_ols_bcm_topic_fit
//   Multiplicative topic-model bias correction (Python ols_bcm_topic):
//     If rho(Gamma) < 1:  b_corr = (I - Gamma)^{-1} b0
//     Else:                b_corr = (I + Gamma) b0       [BCA fallback]
//     V_corr = V0
// ---------------------------------------------------------------------------
struct vmli_bc_result scalar vmli_ols_bcm_topic_fit(
    real colvector Y,
    real matrix    Q,
    real matrix    W,
    real matrix    S,
    real matrix    B,
    real scalar    k)
{
    struct vmli_bc_topic_setup scalar s
    struct vmli_bc_result      scalar bc
    numeric rowvector eig
    real scalar       rho

    s   = vmli_bc_topic_prepare(Y, Q, W, S, B, k)
    eig = eigenvalues(s.Gamma)
    rho = max(abs(eig))

    if (rho < 1) {
        bc.b = lusolve(I(s.d) - s.Gamma, s.b0)
    }
    else {
        bc.b = (I(s.d) + s.Gamma) * s.b0
    }
    bc.V = s.V0

    return(bc)
}

// ---------------------------------------------------------------------------
// _vmli_topic_names(r)
//   Generate "topic_1", "topic_2", ..., "topic_r" as a string rowvector.
// ---------------------------------------------------------------------------
string rowvector _vmli_topic_names(real scalar r)
{
    string rowvector names
    real scalar i

    names = J(1, r, "")
    for (i = 1; i <= r; i++) {
        names[i] = "topic_" + strofreal(i)
    }
    return(names)
}

// ---------------------------------------------------------------------------
// _vmli_topic_read_inputs
//   Shared input reader for both topic-model drivers. Reads Y, Q (possibly
//   with intercept appended), W/S/B from Stata matrices. Validates shapes.
//   Returns Y, Q, W, S, B by reference (via Mata's struct/pointer convention
//   would be cleanest, but here we just return a tuple via out-args).
//
//   We use a small struct to bundle the inputs since Mata can't return
//   multiple typed values cleanly.
// ---------------------------------------------------------------------------
struct vmli_topic_inputs {
    real colvector   Y
    real matrix      Q
    real matrix      W
    real matrix      S
    real matrix      B
    string rowvector q_names
}

struct vmli_topic_inputs scalar _vmli_topic_read_inputs(
    string scalar yname,
    string scalar qnames,
    string scalar tousename,
    string scalar w_matname,
    string scalar s_matname,
    string scalar b_matname,
    real scalar   addconst)
{
    struct vmli_topic_inputs scalar in
    real scalar n

    in.Y = st_data(., yname, tousename)
    n    = rows(in.Y)

    if (qnames != "") {
        in.q_names = tokens(qnames)
        in.Q       = st_data(., in.q_names, tousename)
    }
    else {
        in.q_names = J(1, 0, "")
        in.Q       = J(n, 0, .)
    }

    if (addconst) {
        in.Q       = in.Q, J(n, 1, 1)
        in.q_names = in.q_names, "_cons"
    }

    in.W = st_matrix(w_matname)
    in.S = st_matrix(s_matname)
    in.B = st_matrix(b_matname)

    if (rows(in.W) != n) {
        errprintf("matrix %s must have %g rows (matching Y); got %g\n",
                  w_matname, n, rows(in.W))
        exit(503)
    }
    if (cols(in.S) != cols(in.W)) {
        errprintf("matrix %s must have %g cols (matching W); got %g\n",
                  s_matname, cols(in.W), cols(in.S))
        exit(503)
    }
    if (rows(in.B) != cols(in.W)) {
        errprintf("matrix %s must have %g rows (matching cols(W)); got %g\n",
                  b_matname, cols(in.W), rows(in.B))
        exit(503)
    }

    return(in)
}

// ---------------------------------------------------------------------------
// vmli_ols_bca_topic_run / vmli_ols_bcm_topic_run
//   Stata-facing drivers.
// ---------------------------------------------------------------------------
void vmli_ols_bca_topic_run(
    string scalar yname,
    string scalar qnames,
    string scalar tousename,
    string scalar w_matname,
    string scalar s_matname,
    string scalar b_matname,
    real scalar   k,
    real scalar   addconst,
    string scalar out_b,
    string scalar out_V)
{
    struct vmli_topic_inputs scalar in
    struct vmli_bc_result    scalar bc
    string rowvector names

    in    = _vmli_topic_read_inputs(yname, qnames, tousename,
                                    w_matname, s_matname, b_matname, addconst)
    bc    = vmli_ols_bca_topic_fit(in.Y, in.Q, in.W, in.S, in.B, k)
    names = _vmli_topic_names(rows(in.S)), in.q_names

    _vmli_post_results(bc.b, bc.V, names, rows(in.Y), out_b, out_V)
}

void vmli_ols_bcm_topic_run(
    string scalar yname,
    string scalar qnames,
    string scalar tousename,
    string scalar w_matname,
    string scalar s_matname,
    string scalar b_matname,
    real scalar   k,
    real scalar   addconst,
    string scalar out_b,
    string scalar out_V)
{
    struct vmli_topic_inputs scalar in
    struct vmli_bc_result    scalar bc
    string rowvector names

    in    = _vmli_topic_read_inputs(yname, qnames, tousename,
                                    w_matname, s_matname, b_matname, addconst)
    bc    = vmli_ols_bcm_topic_fit(in.Y, in.Q, in.W, in.S, in.B, k)
    names = _vmli_topic_names(rows(in.S)), in.q_names

    _vmli_post_results(bc.b, bc.V, names, rows(in.Y), out_b, out_V)
}

// ===========================================================================
// One-step joint MLE for a binary AI/ML-generated regressor (Milestone 6)
// Mirrors Python `_one_step_core_with_treatment_idx`.
//
// Likelihood: a Gaussian outcome model where the treatment column at
// `target_idx` is observed with classification noise. The latent treatment X*
// is integrated out via a 2x2 misclassification table parameterized by
// (w00, w01, w10, w11), with heteroskedastic residual variance (sigma0,
// sigma1) optionally constrained to sigma0 == sigma1 (homoskedastic).
//
// Parameter vector theta of length (d + 5) [hetero] or (d + 4) [homo]:
//   theta[1..d]       = b               regression coefficients
//   theta[d+1..d+3]   = (v1, v2, v3)    log-odds parameterization of (w00, w01, w10) vs w11:
//                                       w_ij = exp(v) / (1 + exp(v1) + exp(v2) + exp(v3))
//                                       w11  = 1   / (1 + exp(v1) + exp(v2) + exp(v3))
//   theta[d+4]        = log(sigma0)
//   theta[d+5]        = log(sigma1)     (omitted when homoskedastic)
// ===========================================================================

// ---------------------------------------------------------------------------
// vmli_normal_pdf(y, loc, scale)
//   Gaussian density. `y` and `loc` are column vectors of equal length;
//   `scale` is a positive scalar. Returns a column vector.
// ---------------------------------------------------------------------------
real colvector vmli_normal_pdf(real colvector y,
                               real colvector loc,
                               real scalar    scale)
{
    real colvector z
    z = (y :- loc) :/ scale
    return( exp(-0.5 :* z:^2) :/ (sqrt(2 * pi()) * scale) )
}

// ---------------------------------------------------------------------------
// vmli_subset_std(u, mask)
//   Population standard deviation of `u` on the subset where mask == 1.
//   Divisor is sum(mask), matching Python `subset_std` (jnp.std with
//   default ddof=0). Returns . if the subset is empty.
// ---------------------------------------------------------------------------
real scalar vmli_subset_std(real colvector u, real colvector mask)
{
    real scalar n_mask, mean_val, var

    n_mask = sum(mask)
    if (n_mask == 0) return(.)

    mean_val = sum(u :* mask) / n_mask
    var      = sum(mask :* (u :- mean_val):^2) / n_mask
    return(sqrt(var))
}

// ---------------------------------------------------------------------------
// _vmli_log_pdf(Y, loc, log_scale, dist_code, df)
//   Vectorized log-pdf for the parametric residual menu used by one_step
//   and one_step_gmm.
//
//     dist_code = 1   Normal     log N(y; loc, sigma=exp(log_scale))
//     dist_code = 2   Laplace    log Laplace(y; loc, b=exp(log_scale))
//     dist_code = 3   Student-t  log t_df(y; loc, sigma=exp(log_scale))
//
//   `Y` and `loc` are column vectors of equal length; `log_scale` is a
//   real scalar; `df` is the Student-t degrees of freedom (ignored when
//   dist_code != 3). All computations are in log-space so the value stays
//   finite at any log_scale (no overflow from forming sigma directly).
//
//   For Student-t the constant
//       c_t = lngamma((df+1)/2) - lngamma(df/2) - 0.5 log(df pi)
//   is computed once per call (df is the same for the whole sample).
// ---------------------------------------------------------------------------
real colvector _vmli_log_pdf(real colvector Y,
                             real colvector loc,
                             real scalar    log_scale,
                             real scalar    dist_code,
                             real scalar    df)
{
    real colvector z, lp
    real scalar    inv_scale, log_sqrt2pi, c_t, half_dfp1

    inv_scale = exp(-log_scale)
    z         = (Y :- loc) :* inv_scale

    if (dist_code == 1) {
        log_sqrt2pi = 0.5 * log(2 * pi())
        lp = (-log_sqrt2pi - log_scale) :- 0.5 :* z:^2
    }
    else if (dist_code == 2) {
        lp = (-log(2) - log_scale) :- abs(z)
    }
    else if (dist_code == 3) {
        c_t       = lngamma((df + 1) / 2) - lngamma(df / 2) - 0.5 * log(df * pi())
        half_dfp1 = 0.5 * (df + 1)
        lp = (c_t - log_scale) :- half_dfp1 :* log(1 :+ z:^2 :/ df)
    }
    else {
        errprintf("vmli: unknown distribution code %g\n", dist_code)
        exit(459)
    }
    return(lp)
}

// ---------------------------------------------------------------------------
// vmli_one_step_pars struct
//   Unpacked parameter container, returned by vmli_theta_to_pars().
// ---------------------------------------------------------------------------
struct vmli_one_step_pars {
    real colvector  b
    real scalar     w00, w01, w10, w11
    real scalar     sigma0, sigma1
}

// ---------------------------------------------------------------------------
// vmli_theta_to_pars(theta, d, homoskedastic)
//   Unpack the optimization parameter vector. Mirrors Python
//   `theta_to_pars_jax`. Probabilities are recovered via softmax-with-
//   reference: e^v[i] / (1 + sum e^v[j]) for i=00,01,10 and 1/(1+sum) for 11.
// ---------------------------------------------------------------------------
struct vmli_one_step_pars scalar vmli_theta_to_pars(real colvector theta,
                                                   real scalar    d,
                                                   real scalar    homoskedastic)
{
    struct vmli_one_step_pars scalar p
    real scalar e00, e01, e10, s

    p.b  = theta[1..d]

    e00 = exp(theta[d+1])
    e01 = exp(theta[d+2])
    e10 = exp(theta[d+3])
    s   = 1 + e00 + e01 + e10
    p.w00 = e00 / s
    p.w01 = e01 / s
    p.w10 = e10 / s
    p.w11 = 1   / s

    p.sigma0 = exp(theta[d+4])
    p.sigma1 = homoskedastic ? p.sigma0 : exp(theta[d+5])

    return(p)
}

// ---------------------------------------------------------------------------
// vmli_one_step_starts(Y, Xhat, target_idx, homoskedastic)
//   Closed-form starting values. Mirrors Python
//   `get_starting_values_unlabeled_jax_with_treatment_idx`:
//     1. Fit OLS, take residual std (population).
//     2. Impute latent X* by comparing pdf(y; mu, sigma) to pdf(y; mu +/- b_t,
//        sigma) for each observation.
//     3. Form joint frequencies of (observed, imputed), floor at 0.001,
//        renormalize, take log-odds vs w11.
//     4. Compute residual stds on the imputed-0 and imputed-1 subsets;
//        fill NaNs by copying across.
// ---------------------------------------------------------------------------
real colvector vmli_one_step_starts(real colvector Y,
                                    real matrix    Xhat,
                                    real scalar    target_idx,
                                    real scalar    homoskedastic)
{
    real colvector b, u, mu, ind, X_imp, mask0, mask1, theta0, cond1, cond2
    real colvector pdf_ref, pdf_alt1, pdf_alt2
    real scalar    sigma, te, n
    real scalar    f00, f01, f10, f11, sw
    real scalar    w00, w01, w10, w11, v1, v2, v3
    real scalar    sigma0, sigma1, p_val, sigma_comb

    n = rows(Y)

    b = invsym(quadcross(Xhat, Xhat)) * quadcross(Xhat, Y)
    u = Y - Xhat * b

    // Population std (ddof=0), to match jnp.std default
    sigma = sqrt( sum( (u :- mean(u)):^2 ) / n )

    mu  = Xhat * b
    te  = b[target_idx]
    ind = Xhat[., target_idx]

    pdf_ref  = vmli_normal_pdf(Y, mu, sigma)
    pdf_alt1 = vmli_normal_pdf(Y, mu :- te, sigma)
    pdf_alt2 = vmli_normal_pdf(Y, mu :+ te, sigma)

    // For obs with ind==1: X*=1 if pdf(Y;mu,sigma) > pdf(Y;mu-te,sigma)
    // For obs with ind==0: X*=1 if pdf(Y;mu+te,sigma) > pdf(Y;mu,sigma)
    cond1 = (pdf_ref  :> pdf_alt1)
    cond2 = (pdf_alt2 :> pdf_ref)
    X_imp = (ind :== 1) :* cond1 :+ (ind :!= 1) :* cond2

    f00 = sum( (ind :== 0) :* (X_imp :== 0) ) / n
    f01 = sum( (ind :== 0) :* (X_imp :== 1) ) / n
    f10 = sum( (ind :== 1) :* (X_imp :== 0) ) / n
    f11 = sum( (ind :== 1) :* (X_imp :== 1) ) / n

    w00 = max((f00, 0.001))
    w01 = max((f01, 0.001))
    w10 = max((f10, 0.001))
    w11 = max((f11, 0.001))
    sw  = w00 + w01 + w10 + w11
    w00 = w00 / sw
    w01 = w01 / sw
    w10 = w10 / sw
    w11 = w11 / sw

    v1 = log(w00 / w11)
    v2 = log(w01 / w11)
    v3 = log(w10 / w11)

    mask0  = (X_imp :== 0)
    mask1  = (X_imp :== 1)
    sigma0 = vmli_subset_std(u, mask0)
    sigma1 = vmli_subset_std(u, mask1)
    if (sigma0 == . | sigma0 == 0) sigma0 = sigma1
    if (sigma1 == . | sigma1 == 0) sigma1 = sigma0

    if (homoskedastic) {
        p_val      = mean(X_imp)
        sigma_comb = sigma1 * p_val + sigma0 * (1 - p_val)
        theta0     = b \ v1 \ v2 \ v3 \ log(sigma_comb)
    }
    else {
        theta0 = b \ v1 \ v2 \ v3 \ log(sigma0) \ log(sigma1)
    }
    return(theta0)
}

// ---------------------------------------------------------------------------
// vmli_one_step_obj(todo, p, Y, Xhat, target_idx, homoskedastic,
//                   dist_code, df, fv, g, H)
//   Mata `optimize()` d0 evaluator. Computes the negative log-likelihood
//   of the joint model in log-space so it stays finite at any theta:
//
//     log w_ij computed via softmax-with-reference (log_den = logsumexp);
//     log f(y; mu, sigma) computed by _vmli_log_pdf (normal | laplace | t)
//     directly from log_scale, never from the pdf value itself; and the
//     per-observation mixture combined with a logsumexp(la, lb) over the
//     two latent-class contributions.
//
//   For dist_code != 1 the parameter `sigma` is reinterpreted as the
//   distribution's scale parameter (Laplace b, Student-t scale). The
//   theta layout is unchanged.
// ---------------------------------------------------------------------------
void vmli_one_step_obj(real scalar    todo,
                       real rowvector p,
                       real colvector Y,
                       real matrix    Xhat,
                       real scalar    target_idx,
                       real scalar    homoskedastic,
                       real scalar    dist_code,
                       real scalar    df,
                       real scalar    fv,
                       real rowvector g,
                       real matrix    H)
{
    real colvector theta, b, mu, ind, m0, m1
    real colvector la1, lb1, la0, lb0, mx1, mx0, lp1, lp0, log_p
    real scalar    d, te, v1, v2, v3
    real scalar    log_sigma0, log_sigma1
    real scalar    M, log_den
    real scalar    log_w00, log_w01, log_w10, log_w11

    theta = p'
    d     = cols(Xhat)
    b     = theta[1..d]

    v1 = theta[d+1]
    v2 = theta[d+2]
    v3 = theta[d+3]

    // logsumexp(0, v1, v2, v3) — the normalizer of the softmax-with-reference
    M       = max((0, v1, v2, v3))
    log_den = M + log(exp(0 - M) + exp(v1 - M) + exp(v2 - M) + exp(v3 - M))

    log_w00 = v1 - log_den
    log_w01 = v2 - log_den
    log_w10 = v3 - log_den
    log_w11 = 0  - log_den

    log_sigma0 = theta[d+4]
    log_sigma1 = homoskedastic ? log_sigma0 : theta[d+5]

    mu  = Xhat * b
    te  = b[target_idx]
    ind = Xhat[., target_idx]

    m0 = mu :- te
    m1 = mu :+ te

    // ind == 1:  log P(Y | ind=1) = logsumexp(
    //                                log_w11 + log f(Y; mu, sigma1),
    //                                log_w10 + log f(Y; m0, sigma0))
    // ind == 0:  log P(Y | ind=0) = logsumexp(
    //                                log_w01 + log f(Y; m1, sigma1),
    //                                log_w00 + log f(Y; mu, sigma0))
    la1 = log_w11 :+ _vmli_log_pdf(Y, mu, log_sigma1, dist_code, df)
    lb1 = log_w10 :+ _vmli_log_pdf(Y, m0, log_sigma0, dist_code, df)
    la0 = log_w01 :+ _vmli_log_pdf(Y, m1, log_sigma1, dist_code, df)
    lb0 = log_w00 :+ _vmli_log_pdf(Y, mu, log_sigma0, dist_code, df)

    // Element-wise logsumexp(la, lb) = max + log(exp(la-max) + exp(lb-max))
    mx1 = (la1 :> lb1) :* la1 :+ (la1 :<= lb1) :* lb1
    lp1 = mx1 :+ log(exp(la1 :- mx1) :+ exp(lb1 :- mx1))

    mx0 = (la0 :> lb0) :* la0 :+ (la0 :<= lb0) :* lb0
    lp0 = mx0 :+ log(exp(la0 :- mx0) :+ exp(lb0 :- mx0))

    log_p = (ind :== 1) :* lp1 :+ (ind :!= 1) :* lp0
    fv    = -sum(log_p)
}

// ---------------------------------------------------------------------------
// vmli_one_step_result struct
//   Container for one_step output: coefficient block and its variance.
// ---------------------------------------------------------------------------
struct vmli_one_step_result {
    real colvector  b
    real matrix     V
    real scalar     converged
    real scalar     iterations
    real scalar     loglik
}

// ---------------------------------------------------------------------------
// vmli_one_step_fit(Y, Xhat, target_idx, homoskedastic, dist_code, df)
//   Run starts -> optimize() (BFGS, d0 evaluator) -> numerical Hessian ->
//   V = pinv(H)[1..d, 1..d]. Mirrors Python `_one_step_jax_core_with_treatment_idx`.
//   Starting values use the Gaussian pdf for the latent-X* imputation step
//   regardless of dist_code -- this matches Python, whose
//   get_starting_values_unlabeled_jax_with_treatment_idx hard-codes a
//   Gaussian pdf_func in the start phase even when distribution= is passed.
// ---------------------------------------------------------------------------
struct vmli_one_step_result scalar vmli_one_step_fit(real colvector Y,
                                                    real matrix    Xhat,
                                                    real scalar    target_idx,
                                                    real scalar    homoskedastic,
                                                    real scalar    dist_code,
                                                    real scalar    df)
{
    struct vmli_one_step_result scalar res
    transmorphic S
    real rowvector theta0, theta_hat
    real matrix    H, V_full
    real scalar    d

    d      = cols(Xhat)
    theta0 = vmli_one_step_starts(Y, Xhat, target_idx, homoskedastic)'

    S = optimize_init()
    optimize_init_evaluator(S, &vmli_one_step_obj())
    optimize_init_evaluatortype(S, "d0")
    optimize_init_argument(S, 1, Y)
    optimize_init_argument(S, 2, Xhat)
    optimize_init_argument(S, 3, target_idx)
    optimize_init_argument(S, 4, homoskedastic)
    optimize_init_argument(S, 5, dist_code)
    optimize_init_argument(S, 6, df)
    optimize_init_which(S, "min")
    optimize_init_technique(S, "bfgs")
    optimize_init_conv_maxiter(S, 500)
    optimize_init_conv_ptol(S, 1e-12)
    optimize_init_conv_vtol(S, 1e-12)
    optimize_init_tracelevel(S, "none")
    optimize_init_params(S, theta0)

    theta_hat = optimize(S)
    H         = optimize_result_Hessian(S)
    V_full    = pinv(H)

    res.b          = theta_hat[1..d]'
    res.V          = V_full[1..d, 1..d]
    res.converged  = optimize_result_converged(S)
    res.iterations = optimize_result_iterations(S)
    res.loglik     = -optimize_result_value(S)

    return(res)
}

// ---------------------------------------------------------------------------
// vmli_one_step_run(yname, xnames, tousename, addconst, gen_var, homosked,
//                   dist_code, df, bname, Vname)
//   Stata-facing driver. Reads Y/Xhat from the active sample, validates
//   the treatment column is binary 0/1, locates target_idx, runs the MLE,
//   and posts the d-block of coefficients and variance.
// ---------------------------------------------------------------------------
void vmli_one_step_run(string scalar yname,
                       string scalar xnames,
                       string scalar tousename,
                       real scalar   addconst,
                       string scalar gen_var,
                       real scalar   homoskedastic,
                       real scalar   dist_code,
                       real scalar   df,
                       string scalar bname,
                       string scalar Vname)
{
    real colvector  Y, tcol, uvals
    real matrix     X
    string rowvector names
    real scalar     target_idx
    struct vmli_one_step_result scalar res

    Y     = st_data(., yname, tousename)
    names = tokens(xnames)
    X     = st_data(., names, tousename)

    if (addconst) {
        X     = X, J(rows(X), 1, 1)
        names = names, "_cons"
    }

    target_idx = _vmli_find_col(names, gen_var)

    tcol  = X[., target_idx]
    uvals = uniqrows(tcol)
    if (rows(uvals) != 2) {
        errprintf("treatment variable '%s' must be binary 0/1 (%g distinct values found)\n",
                  gen_var, rows(uvals))
        exit(459)
    }
    if (uvals[1] != 0 | uvals[2] != 1) {
        errprintf("treatment variable '%s' must take values 0 and 1\n", gen_var)
        exit(459)
    }

    res = vmli_one_step_fit(Y, X, target_idx, homoskedastic, dist_code, df)

    _vmli_post_results(res.b, res.V, names, rows(X), bname, Vname)

    st_numscalar("r(converged)",  res.converged)
    st_numscalar("r(iterations)", res.iterations)
    st_numscalar("r(ll)",         res.loglik)
}

// ===========================================================================
// One-step joint MLE with Gaussian-mixture residuals (Milestone 7)
// Mirrors Python `_one_step_gaussian_mixture_core`.
//
// Outcome model: same misclassification structure as `one_step`, but the
// per-class residual is a k-component Gaussian mixture rather than a single
// Gaussian. The mixture in each latent class has its own component weights,
// mean offsets, and (heteroskedastic) component stds.
//
// Parameter vector theta layout (heteroskedastic only, length d - 1 + 6k):
//   theta[1..d]              b              regression coefficients
//   theta[d+1..d+3]          v_misclass     logits of (w00, w01, w10) vs w11
//   theta[d+4..d+2+k]        v0             logits of cluster-0 component
//                                           weights (length k-1; last comp ref)
//   theta[d+3+k..d+1+2k]     v1             logits of cluster-1 component weights
//   theta[d+2+2k..d+3k]      m0p            cluster-0 mean increments (length k-1)
//   theta[d+1+3k..d-1+4k]    m1p            cluster-1 mean increments
//   theta[d+4k..d-1+5k]      logs0          log of cluster-0 component stds (length k)
//   theta[d+5k..d-1+6k]      logs1          log of cluster-1 component stds
//
// Means are recovered as mu_raw = [cumsum(m_incr), 0] then centered to a
// weighted mean of zero so the mixture is identified relative to the
// regression-mean component.
//
// Variance: numerical Hessian at the best multistart optimum, V = pinv(H)[1..d, 1..d].
//
// Upstream caveat: the Python `get_starting_values_unlabeled_gaussian_mixture`
// homoskedastic branch returns a theta with length d + 1 + 3k, but
// `unpack_theta(homosked=True)` reads d - 1 + 5k slots. We only port the
// heteroskedastic branch in M7; the homoskedastic option is rejected with
// an explanatory error. See Notes/porting_decisions.md.
// ===========================================================================

// ---------------------------------------------------------------------------
// vmli_log_softmax_ref(v)
//   log-softmax with an implicit reference of 0 appended at the end.
//   Input  v: real rowvector, length m.
//   Output:   real rowvector, length m+1, where
//             out[j] = v[j] - logsumexp(0, v[1..m]) for j <= m,
//             out[m+1] = -logsumexp(0, v[1..m]).
//   Numerically stable via max-shift.
// ---------------------------------------------------------------------------
real rowvector vmli_log_softmax_ref(real rowvector v)
{
    real rowvector all_v, log_w
    real scalar    M, log_den

    all_v   = v, 0
    M       = max(all_v)
    log_den = M + log(sum(exp(all_v :- M)))
    log_w   = all_v :- log_den
    return(log_w)
}

// ---------------------------------------------------------------------------
// vmli_centered_means(m_incr, log_w)
//   Build the centered mean vector for a k-component mixture from k-1
//   increment parameters. Mirrors Python's
//     mu_raw = jnp.concatenate([jnp.cumsum(m_incr), jnp.zeros(1)])
//     mu     = mu_raw - jnp.dot(omega, mu_raw)
//   so that sum_j w_j mu_j = 0 (the mixture is centered at zero, leaving
//   the regression mu as the location parameter).
// ---------------------------------------------------------------------------
real rowvector vmli_centered_means(real rowvector m_incr, real rowvector log_w)
{
    real rowvector mu_raw, w, mu
    real scalar    shift

    mu_raw = runningsum(m_incr), 0
    w      = exp(log_w)
    shift  = sum(w :* mu_raw)
    mu     = mu_raw :- shift
    return(mu)
}

// ---------------------------------------------------------------------------
// vmli_log_mixture_pdf(Y, center, log_w, mu_off, log_sigma)
//   For each observation i, returns
//     log_mix[i] = log( sum_j exp(log_w[j]) N(Y[i]; center[i] + mu_off[j], sigma[j]) )
//   computed entirely in log-space.
//
//   Y         : n x 1
//   center    : n x 1   (e.g., the regression mean or a shifted version)
//   log_w     : 1 x k   (log component weights)
//   mu_off    : 1 x k   (component mean offsets relative to center)
//   log_sigma : 1 x k   (log component stds)
// ---------------------------------------------------------------------------
real colvector vmli_log_mixture_pdf(real colvector Y,
                                    real colvector center,
                                    real rowvector log_w,
                                    real rowvector mu_off,
                                    real rowvector log_sigma)
{
    real matrix    log_terms
    real colvector log_mix, max_t, resid, z
    real rowvector inv_sigma
    real scalar    k, log_sqrt2pi, j

    k           = cols(log_w)
    log_sqrt2pi = 0.5 * log(2 * pi())
    inv_sigma   = exp(-log_sigma)
    resid       = Y - center

    log_terms = J(rows(Y), k, .)
    for (j = 1; j <= k; j++) {
        z              = (resid :- mu_off[j]) :* inv_sigma[j]
        log_terms[., j] = (log_w[j] - log_sqrt2pi - log_sigma[j]) :- 0.5 :* z:^2
    }

    max_t   = rowmax(log_terms)
    log_mix = max_t :+ log(rowsum(exp(log_terms :- max_t)))
    return(log_mix)
}

// ---------------------------------------------------------------------------
// vmli_gmm_starts(Y, Xhat, target_idx, k, homoskedastic)
//   Heteroskedastic-only starting values, mirroring Python
//   `get_starting_values_unlabeled_gaussian_mixture(homosked=False)`:
//     b, sigma  : OLS fit + population residual std (ddof=0).
//     X_imputed : per-obs latent-treatment imputation comparing pdfs.
//     w_ij      : floor 1e-3, renormalize, log-odds vs w11.
//     sigma_g   : residual std on imputed-class subset (fill NaN by copying).
//     v0_block  : 0.01 * arange(k-1)      (cluster-0 component logits)
//     v1_block  : 0.005 * arange(k-1)     (cluster-1 component logits)
//     m0_block  : 0.1 * arange(k-1)       (cluster-0 mean increments)
//     m1_block  : 0.15 * arange(k-1)      (cluster-1 mean increments)
//     logs0/1   : log(sigma_g) repeated k times
// ---------------------------------------------------------------------------
real colvector vmli_gmm_starts(real colvector Y,
                               real matrix    Xhat,
                               real scalar    target_idx,
                               real scalar    k,
                               real scalar    homoskedastic)
{
    real colvector b, u, mu, ind, X_imp, mask0, mask1, theta0
    real colvector pdf_ref, pdf_alt1, pdf_alt2, cond1, cond2
    real scalar    sigma, te, n
    real scalar    f00, f01, f10, f11, sw
    real scalar    w00, w01, w10, w11, v1, v2, v3
    real scalar    sigma0, sigma1, log_sigma0, log_sigma1
    real rowvector v0_blk, v1_blk, m0_blk, m1_blk, logs0_blk, logs1_blk, idx

    if (homoskedastic) {
        errprintf("homoskedastic Gaussian mixture not currently supported ")
        errprintf("(upstream Python theta-length mismatch); use heteroskedastic.\n")
        exit(459)
    }

    n = rows(Y)

    b     = invsym(quadcross(Xhat, Xhat)) * quadcross(Xhat, Y)
    u     = Y - Xhat * b
    sigma = sqrt(sum((u :- mean(u)):^2) / n)

    mu  = Xhat * b
    te  = b[target_idx]
    ind = Xhat[., target_idx]

    pdf_ref  = vmli_normal_pdf(Y, mu, sigma)
    pdf_alt1 = vmli_normal_pdf(Y, mu :- te, sigma)
    pdf_alt2 = vmli_normal_pdf(Y, mu :+ te, sigma)

    cond1 = (pdf_ref  :> pdf_alt1)
    cond2 = (pdf_alt2 :> pdf_ref)
    X_imp = (ind :== 1) :* cond1 :+ (ind :!= 1) :* cond2

    f00 = sum((ind :== 0) :* (X_imp :== 0)) / n
    f01 = sum((ind :== 0) :* (X_imp :== 1)) / n
    f10 = sum((ind :== 1) :* (X_imp :== 0)) / n
    f11 = sum((ind :== 1) :* (X_imp :== 1)) / n

    w00 = max((f00, 0.001))
    w01 = max((f01, 0.001))
    w10 = max((f10, 0.001))
    w11 = max((f11, 0.001))
    sw  = w00 + w01 + w10 + w11
    w00 = w00 / sw
    w01 = w01 / sw
    w10 = w10 / sw
    w11 = w11 / sw

    v1 = log(w00 / w11)
    v2 = log(w01 / w11)
    v3 = log(w10 / w11)

    mask0  = (X_imp :== 0)
    mask1  = (X_imp :== 1)
    sigma0 = vmli_subset_std(u, mask0)
    sigma1 = vmli_subset_std(u, mask1)
    if (sigma0 == . | sigma0 == 0) sigma0 = sigma1
    if (sigma1 == . | sigma1 == 0) sigma1 = sigma0
    if (sigma0 == . | sigma0 == 0) sigma0 = sigma
    if (sigma1 == . | sigma1 == 0) sigma1 = sigma
    log_sigma0 = log(sigma0)
    log_sigma1 = log(sigma1)

    idx = 0..(k-2)

    v0_blk    = 0.01  :* idx
    v1_blk    = 0.005 :* idx
    m0_blk    = 0.1   :* idx
    m1_blk    = 0.15  :* idx
    logs0_blk = J(1, k, log_sigma0)
    logs1_blk = J(1, k, log_sigma1)

    theta0 = b           \
             v1 \ v2 \ v3 \
             v0_blk'     \
             v1_blk'     \
             m0_blk'     \
             m1_blk'     \
             logs0_blk'  \
             logs1_blk'
    return(theta0)
}

// ---------------------------------------------------------------------------
// vmli_gmm_obj(todo, p, Y, Xhat, target_idx, k, homoskedastic, fv, g, H)
//   Mata optimize() d0 evaluator: negative log-likelihood of the Gaussian-
//   mixture one-step model, computed entirely in log-space (logsumexp for
//   the softmax denominators AND the inner/outer mixtures). Python's
//   stabilizer of +1e-12 inside log() is matched in log-domain via
//   logaddexp(log_P, log(1e-12)).
// ---------------------------------------------------------------------------
void vmli_gmm_obj(real scalar    todo,
                  real rowvector p,
                  real colvector Y,
                  real matrix    Xhat,
                  real scalar    target_idx,
                  real scalar    k,
                  real scalar    homoskedastic,
                  real scalar    fv,
                  real rowvector g,
                  real matrix    H)
{
    real colvector theta, b, mu, ind
    real colvector mu_minus_te, mu_plus_te
    real colvector log_mix1_h1, log_mix0_h1, log_mix1_h0, log_mix0_h0
    real colvector la, lb, mx_lp, lp_h1, lp_h0, log_p, mx_eps
    real rowvector vm, v0, v1, m0p, m1p, log_w0, log_w1
    real rowvector mu0, mu1, log_s0, log_s1
    real scalar    d, te, M, log_den, te_target
    real scalar    log_w00, log_w01, log_w10, log_w11
    real scalar    log_eps, i

    theta = p'
    d     = cols(Xhat)
    b     = theta[1..d]

    i  = d + 1
    vm = theta[i..i+2]'
    i  = i + 3

    M       = max((0, vm))
    log_den = M + log(exp(0 - M) + sum(exp(vm :- M)))
    log_w00 = vm[1] - log_den
    log_w01 = vm[2] - log_den
    log_w10 = vm[3] - log_den
    log_w11 = 0     - log_den

    v0 = theta[i..i+k-2]'
    i  = i + (k - 1)
    log_w0 = vmli_log_softmax_ref(v0)

    v1 = theta[i..i+k-2]'
    i  = i + (k - 1)
    log_w1 = vmli_log_softmax_ref(v1)

    m0p = theta[i..i+k-2]'
    i   = i + (k - 1)
    mu0 = vmli_centered_means(m0p, log_w0)

    m1p = theta[i..i+k-2]'
    i   = i + (k - 1)
    mu1 = vmli_centered_means(m1p, log_w1)

    log_s0 = theta[i..i+k-1]'
    i      = i + k

    if (homoskedastic) {
        log_s1 = log_s0
    }
    else {
        log_s1 = theta[i..i+k-1]'
        i      = i + k
    }

    mu        = Xhat * b
    te_target = b[target_idx]
    ind       = Xhat[., target_idx]

    mu_minus_te = mu :- te_target
    mu_plus_te  = mu :+ te_target

    // ind == 1: outer mixture over (latent X*=1, sigma1) vs (latent X*=0, sigma0)
    log_mix1_h1 = vmli_log_mixture_pdf(Y, mu,          log_w1, mu1, log_s1)
    log_mix0_h1 = vmli_log_mixture_pdf(Y, mu_minus_te, log_w0, mu0, log_s0)

    la    = log_w11 :+ log_mix1_h1
    lb    = log_w10 :+ log_mix0_h1
    mx_lp = (la :> lb) :* la :+ (la :<= lb) :* lb
    lp_h1 = mx_lp :+ log(exp(la :- mx_lp) :+ exp(lb :- mx_lp))

    // ind == 0: outer mixture over (latent X*=1, sigma1) vs (latent X*=0, sigma0)
    log_mix1_h0 = vmli_log_mixture_pdf(Y, mu_plus_te, log_w1, mu1, log_s1)
    log_mix0_h0 = vmli_log_mixture_pdf(Y, mu,         log_w0, mu0, log_s0)

    la    = log_w01 :+ log_mix1_h0
    lb    = log_w00 :+ log_mix0_h0
    mx_lp = (la :> lb) :* la :+ (la :<= lb) :* lb
    lp_h0 = mx_lp :+ log(exp(la :- mx_lp) :+ exp(lb :- mx_lp))

    log_p = (ind :== 1) :* lp_h1 :+ (ind :!= 1) :* lp_h0

    // Python adds +1e-12 inside the log; in log-space that is
    //   log(exp(log_p) + 1e-12) = logaddexp(log_p, log(1e-12)).
    log_eps = log(1e-12)
    mx_eps  = (log_p :> log_eps) :* log_p :+ (log_p :<= log_eps) :* log_eps
    log_p   = mx_eps :+ log(exp(log_p :- mx_eps) :+ exp(log_eps :- mx_eps))

    fv = -sum(log_p)
}

// ---------------------------------------------------------------------------
// _vmli_gmm_symmetric_slots(theta, d, k)
//   Return a row vector of slot indices (in theta) that are *exactly* zero
//   at the OLS-based starts and lie along an exact symmetry direction of
//   the mixture likelihood. At these slots the objective is analytically
//   flat -- f(theta + h e_j) == f(theta - h e_j) for any h -- so the
//   numerical gradient is identically zero and Mata's BFGS cannot
//   compute a meaningful FD derivative there (it aborts with rc=5,
//   "could not calculate numerical derivatives").
//
//   For k = 2 every entry of (v0, v1, m0p, m1p) is symmetric at the
//   starts; for k >= 3 only the first entry of each block is symmetric
//   (the remaining entries are 0.01*j, 0.005*j, 0.1*j, 0.15*j -- non-zero
//   for j >= 1 in 0-indexed Python). Returns J(1, 0, .) when nothing is
//   symmetric (i.e., a perturbed restart has broken every symmetry).
//
//   The detection is done numerically: any slot whose absolute value is
//   below `1e-12` is treated as on the symmetric subspace.
// ---------------------------------------------------------------------------
real rowvector _vmli_gmm_symmetric_slots(real rowvector theta,
                                         real scalar    d,
                                         real scalar    k)
{
    real rowvector idx, candidates
    real scalar    j, base

    // Slot positions of the first element of each (v0, v1, m0p, m1p) block.
    // theta[d+1..d+3] = vm; v0 starts at d+4.
    base       = d + 4
    candidates = base, base + (k-1), base + 2*(k-1), base + 3*(k-1)

    idx = J(1, 0, .)
    for (j = 1; j <= 4; j++) {
        if (abs(theta[candidates[j]]) < 1e-12) {
            idx = idx, candidates[j]
        }
    }

    // For k = 2 each block has exactly one entry (which is `base + i*0`
    // == base for i=0, base+0 for the next, etc.) so candidates collapse;
    // remove duplicates while preserving order.
    if (cols(idx) > 1) {
        real rowvector uniq
        real scalar    i, m, seen
        uniq = J(1, 0, .)
        m    = cols(idx)
        for (i = 1; i <= m; i++) {
            seen = 0
            for (j = 1; j <= cols(uniq); j++) {
                if (uniq[j] == idx[i]) seen = 1
            }
            if (!seen) uniq = uniq, idx[i]
        }
        idx = uniq
    }
    return(idx)
}

// ---------------------------------------------------------------------------
// _vmli_gmm_obj_reduced(...)
//   Wrapper d0 evaluator over a *reduced* parameter vector that omits the
//   slots listed in `fix_idx` (an external Mata global). The fixed slots
//   take their value from `fix_val` (also external). This lets us run
//   Mata's BFGS on the active dimensions only when the full-dimensional
//   problem has analytically flat directions at the start.
// ---------------------------------------------------------------------------
void _vmli_gmm_obj_reduced(real scalar    todo,
                           real rowvector p_reduced,
                           real colvector Y,
                           real matrix    Xhat,
                           real scalar    target_idx,
                           real scalar    k,
                           real scalar    homoskedastic,
                           real scalar    fv,
                           real rowvector g,
                           real matrix    H)
{
    external real rowvector _vmli_gmm_fix_idx
    external real rowvector _vmli_gmm_fix_val
    external real scalar    _vmli_gmm_full_len

    real rowvector p_full, active_idx, all_idx, in_fix
    real scalar    j, m, n_fix, n_full

    n_full = _vmli_gmm_full_len
    n_fix  = cols(_vmli_gmm_fix_idx)

    // active_idx = setdiff(1..n_full, _vmli_gmm_fix_idx)
    all_idx    = 1..n_full
    in_fix     = J(1, n_full, 0)
    for (j = 1; j <= n_fix; j++) {
        in_fix[_vmli_gmm_fix_idx[j]] = 1
    }
    active_idx = select(all_idx, in_fix :== 0)

    p_full = J(1, n_full, 0)
    m      = cols(active_idx)
    for (j = 1; j <= m; j++) {
        p_full[active_idx[j]] = p_reduced[j]
    }
    for (j = 1; j <= n_fix; j++) {
        p_full[_vmli_gmm_fix_idx[j]] = _vmli_gmm_fix_val[j]
    }

    vmli_gmm_obj(todo, p_full, Y, Xhat, target_idx, k, homoskedastic, fv, g, H)
}

// ---------------------------------------------------------------------------
// vmli_gmm_result struct
//   Container for one_step_gmm output: coefficient block, variance, and
//   multistart diagnostics.
// ---------------------------------------------------------------------------
struct vmli_gmm_result {
    real colvector  b
    real matrix     V
    real scalar     converged
    real scalar     iterations
    real scalar     loglik
    real scalar     best_idx       // index of the multistart attempt that won
    real scalar     n_finished     // number of attempts that returned a finite value
}

// ---------------------------------------------------------------------------
// vmli_gmm_fit(Y, Xhat, target_idx, k, homoskedastic, nguess, maxiter, seed)
//   Multistart fit. The i=0 attempt uses the OLS-derived starting values
//   directly; subsequent attempts add Gaussian perturbations with the same
//   schedule Python uses:
//     noise_scale = 0.05 + 0.02 * (i / nguess)
//     scale on coefficient slots: 0.5
//     scale on misclass-logit slots: 0.3
//   The best attempt by minimum nll (and finite) wins. Numerical Hessian is
//   evaluated at the winning theta; V = pinv(H)[1..d, 1..d].
//
//   RNG note: Mata uses Stata's RNG (set via rseed(seed)), not JAX's PRNGKey.
//   For parity, the i=0 attempt is identical across implementations, so
//   when the OLS-based starts reach the global minimum (the typical case
//   for an identified mixture problem), the best-of-nguess agrees up to
//   numerical-Hessian precision regardless of RNG differences.
// ---------------------------------------------------------------------------
struct vmli_gmm_result scalar vmli_gmm_fit(real colvector Y,
                                           real matrix    Xhat,
                                           real scalar    target_idx,
                                           real scalar    k,
                                           real scalar    homoskedastic,
                                           real scalar    nguess,
                                           real scalar    maxiter,
                                           real scalar    seed)
{
    struct vmli_gmm_result scalar res
    transmorphic   S
    real rowvector theta0, theta_try, theta_hat, best_theta, noise, scale_vec
    real rowvector sym, active, all_idx, in_fix, theta_try_red, theta_hat_red
    real matrix    H, V_full
    real scalar    d, p_len, i, j, rc, best_loss, this_loss, noise_scale
    real scalar    best_idx, n_finished, conv_best, iter_best, use_reduced

    external real rowvector _vmli_gmm_fix_idx
    external real rowvector _vmli_gmm_fix_val
    external real scalar    _vmli_gmm_full_len

    d         = cols(Xhat)
    theta0    = vmli_gmm_starts(Y, Xhat, target_idx, k, homoskedastic)'
    p_len     = cols(theta0)

    // Scale multiplier for the noise vector: 0.5 on b slots, 0.3 on the
    // three misclassification logits, 1.0 elsewhere.
    scale_vec           = J(1, p_len, 1)
    scale_vec[1..d]     = J(1, d, 0.5)
    scale_vec[d+1..d+3] = J(1, 3, 0.3)

    rseed(seed)

    best_loss  = .
    best_theta = theta0
    best_idx   = 0
    n_finished = 0
    conv_best  = 0
    iter_best  = 0

    for (i = 0; i < nguess; i++) {

        if (i == 0) {
            theta_try = theta0
        }
        else {
            noise_scale = 0.05 + 0.02 * (i / nguess)
            noise       = rnormal(1, p_len, 0, 1) :* (noise_scale :* scale_vec)
            theta_try   = theta0 :+ noise
        }

        // Detect exact symmetry slots in theta_try. If any exist, optimize
        // over the reduced parameter space with those slots fixed at 0
        // (their value at the symmetric subspace). This matches what
        // JAX-LBFGS effectively does via exact-zero autodiff gradients;
        // Mata's numerical FD would abort otherwise.
        sym         = _vmli_gmm_symmetric_slots(theta_try, d, k)
        use_reduced = (cols(sym) > 0)

        if (use_reduced) {
            _vmli_gmm_fix_idx  = sym
            _vmli_gmm_fix_val  = J(1, cols(sym), 0)
            _vmli_gmm_full_len = p_len

            all_idx = 1..p_len
            in_fix  = J(1, p_len, 0)
            for (j = 1; j <= cols(sym); j++) {
                in_fix[sym[j]] = 1
            }
            active        = select(all_idx, in_fix :== 0)
            theta_try_red = theta_try[active]

            S = optimize_init()
            optimize_init_evaluator(S, &_vmli_gmm_obj_reduced())
            optimize_init_evaluatortype(S, "d0")
            optimize_init_argument(S, 1, Y)
            optimize_init_argument(S, 2, Xhat)
            optimize_init_argument(S, 3, target_idx)
            optimize_init_argument(S, 4, k)
            optimize_init_argument(S, 5, homoskedastic)
            optimize_init_which(S, "min")
            optimize_init_technique(S, "bfgs")
            optimize_init_conv_maxiter(S, maxiter)
            optimize_init_conv_ptol(S, 1e-12)
            optimize_init_conv_vtol(S, 1e-12)
            optimize_init_tracelevel(S, "none")
            optimize_init_params(S, theta_try_red)

            rc = _optimize(S)
            if (rc != 0) continue

            theta_hat_red = optimize_result_params(S)
            this_loss     = optimize_result_value(S)

            // Map the reduced theta back to full length.
            theta_hat = J(1, p_len, 0)
            for (j = 1; j <= cols(active); j++) {
                theta_hat[active[j]] = theta_hat_red[j]
            }
        }
        else {
            S = optimize_init()
            optimize_init_evaluator(S, &vmli_gmm_obj())
            optimize_init_evaluatortype(S, "d0")
            optimize_init_argument(S, 1, Y)
            optimize_init_argument(S, 2, Xhat)
            optimize_init_argument(S, 3, target_idx)
            optimize_init_argument(S, 4, k)
            optimize_init_argument(S, 5, homoskedastic)
            optimize_init_which(S, "min")
            optimize_init_technique(S, "bfgs")
            optimize_init_conv_maxiter(S, maxiter)
            optimize_init_conv_ptol(S, 1e-12)
            optimize_init_conv_vtol(S, 1e-12)
            optimize_init_tracelevel(S, "none")
            optimize_init_params(S, theta_try)

            rc = _optimize(S)
            if (rc != 0) continue

            theta_hat = optimize_result_params(S)
            this_loss = optimize_result_value(S)
        }

        if (this_loss >= .) continue
        n_finished = n_finished + 1

        if (best_loss >= . | this_loss < best_loss) {
            best_loss  = this_loss
            best_theta = theta_hat
            best_idx   = i
            conv_best  = optimize_result_converged(S)
            iter_best  = optimize_result_iterations(S)
        }
    }

    if (best_loss >= .) {
        errprintf("one_step_gmm: all %g multistart attempts failed to return a finite value\n", nguess)
        exit(498)
    }

    // Final Hessian at the best theta. If best_theta still lies on a
    // symmetric subspace (flat gradient in some slots), evaluating the
    // full-dimensional Hessian via Mata FD will fail with rc=5 in those
    // directions. We restrict the Hessian to the active subspace; the
    // block on the coefficient slots (1..d) is the same in either case,
    // because the Hessian is block-diagonal between active and symmetric
    // slots at the symmetric point (off-diagonals are derivatives of an
    // identically-zero gradient).
    sym = _vmli_gmm_symmetric_slots(best_theta, d, k)

    if (cols(sym) > 0) {
        _vmli_gmm_fix_idx  = sym
        _vmli_gmm_fix_val  = J(1, cols(sym), 0)
        _vmli_gmm_full_len = p_len

        all_idx = 1..p_len
        in_fix  = J(1, p_len, 0)
        for (j = 1; j <= cols(sym); j++) {
            in_fix[sym[j]] = 1
        }
        active = select(all_idx, in_fix :== 0)

        S = optimize_init()
        optimize_init_evaluator(S, &_vmli_gmm_obj_reduced())
        optimize_init_evaluatortype(S, "d0")
        optimize_init_argument(S, 1, Y)
        optimize_init_argument(S, 2, Xhat)
        optimize_init_argument(S, 3, target_idx)
        optimize_init_argument(S, 4, k)
        optimize_init_argument(S, 5, homoskedastic)
        optimize_init_which(S, "min")
        optimize_init_technique(S, "bfgs")
        optimize_init_conv_maxiter(S, 0)
        optimize_init_tracelevel(S, "none")
        optimize_init_params(S, best_theta[active])
        (void) _optimize(S)

        H = optimize_result_Hessian(S)
        // Coefficient slots 1..d are never fixed (fix_idx lives in the
        // mixture-structure block, slots d+4 onward), so they sit at
        // positions 1..d inside `active`. Hence pinv(H_reduced)[1..d, 1..d]
        // is the coefficient variance.
        V_full = pinv(H)
    }
    else {
        S = optimize_init()
        optimize_init_evaluator(S, &vmli_gmm_obj())
        optimize_init_evaluatortype(S, "d0")
        optimize_init_argument(S, 1, Y)
        optimize_init_argument(S, 2, Xhat)
        optimize_init_argument(S, 3, target_idx)
        optimize_init_argument(S, 4, k)
        optimize_init_argument(S, 5, homoskedastic)
        optimize_init_which(S, "min")
        optimize_init_technique(S, "bfgs")
        optimize_init_conv_maxiter(S, 0)
        optimize_init_tracelevel(S, "none")
        optimize_init_params(S, best_theta)
        (void) _optimize(S)

        H      = optimize_result_Hessian(S)
        V_full = pinv(H)
    }

    res.b          = best_theta[1..d]'
    res.V          = V_full[1..d, 1..d]
    res.converged  = conv_best
    res.iterations = iter_best
    res.loglik     = -best_loss
    res.best_idx   = best_idx
    res.n_finished = n_finished

    return(res)
}

// ---------------------------------------------------------------------------
// vmli_gmm_run(yname, xnames, tousename, addconst, gen_var, k, homosked,
//              nguess, maxiter, seed, bname, Vname)
//   Stata-facing driver. Reads Y/Xhat from the active sample, validates
//   the treatment column is binary 0/1, runs the multistart MLE, and
//   posts the d-block of coefficients and variance.
// ---------------------------------------------------------------------------
void vmli_gmm_run(string scalar yname,
                  string scalar xnames,
                  string scalar tousename,
                  real scalar   addconst,
                  string scalar gen_var,
                  real scalar   k,
                  real scalar   homoskedastic,
                  real scalar   nguess,
                  real scalar   maxiter,
                  real scalar   seed,
                  string scalar bname,
                  string scalar Vname)
{
    real colvector  Y, tcol, uvals
    real matrix     X
    string rowvector names
    real scalar     target_idx
    struct vmli_gmm_result scalar res

    Y     = st_data(., yname, tousename)
    names = tokens(xnames)
    X     = st_data(., names, tousename)

    if (addconst) {
        X     = X, J(rows(X), 1, 1)
        names = names, "_cons"
    }

    target_idx = _vmli_find_col(names, gen_var)

    tcol  = X[., target_idx]
    uvals = uniqrows(tcol)
    if (rows(uvals) != 2) {
        errprintf("treatment variable '%s' must be binary 0/1 (%g distinct values found)\n",
                  gen_var, rows(uvals))
        exit(459)
    }
    if (uvals[1] != 0 | uvals[2] != 1) {
        errprintf("treatment variable '%s' must take values 0 and 1\n", gen_var)
        exit(459)
    }

    res = vmli_gmm_fit(Y, X, target_idx, k, homoskedastic, nguess, maxiter, seed)

    _vmli_post_results(res.b, res.V, names, rows(X), bname, Vname)

    st_numscalar("r(converged)",  res.converged)
    st_numscalar("r(iterations)", res.iterations)
    st_numscalar("r(ll)",         res.loglik)
    st_numscalar("r(best_idx)",   res.best_idx)
    st_numscalar("r(n_finished)", res.n_finished)
}

end
