*! test_ols_bca_topic.do -- Milestone 4 parity check for ols_bca_topic
*!
*! Generates synthetic topic-model data (small n, small dims) inside
*! Mata, posts y/q variables to Stata, posts W/S/B as Stata matrices,
*! then exercises ols_bca_topic with:
*!   1. k=0 sanity: should reduce to ols on (Theta, Q).
*!   2. Formula check: hand-computed Mata reference matches .ado output.
*!   3. Names: e(b) columns are intercept-first then alphabetical.

version 17
set more off
clear all

* --- Locate project root ---
local pkg "ValidMLInference-stata"
local marker "`pkg'/src/ado/ols_bca_topic.ado"
local found 0
foreach up in "." ".." "../.." {
    cap confirm file "`up'/`marker'"
    if !_rc {
        qui cd "`up'"
        local found 1
        continue, break
    }
}
if !`found' {
    di as error "test_ols_bca_topic.do: cannot locate `marker' from cwd"
    exit 601
}
local proj "`c(pwd)'"

qui cd "`pkg'/src/mata"
qui do _build.do
qui cd "`proj'"

qui adopath ++ "`pkg'/src/ado"
qui adopath ++ "`pkg'/src/mata"
qui mata: mata mlib index

cap program drop ols
cap program drop ols_bca_topic

* ---------------------------------------------------------------------------
* Synthetic topic-model data, small dims
*   n=300 documents, r=3 topics, c=5 mixture dim, v=8 vocab, q=2 covariates
* ---------------------------------------------------------------------------
clear
set obs 300
gen byte _id = _n

mata:
    rseed(20260518)
    nn = 300
    rr = 3
    cc = 5
    vv = 8
    qq = 2
    W = uniform(nn, cc)
    W = W :/ rowsum(W)
    S = uniform(rr, cc)
    S = S :/ rowsum(S)
    B = uniform(cc, vv)
    B = B :/ rowsum(B)
    Theta = W * S'
    Q = rnormal(nn, qq, 0, 1)
    bt = (1.0 \ 0.5 \ -0.4)
    bq = (0.2 \ -0.3)
    cval = 0.7
    eps = rnormal(nn, 1, 0, 1)
    Y = (Theta * bt + Q * bq + eps) :+ cval
    (void) st_addvar("double", ("y", "q1", "q2"))
    st_store(., "y",  Y)
    st_store(., "q1", Q[, 1])
    st_store(., "q2", Q[, 2])
    st_matrix("W", W)
    st_matrix("S", S)
    st_matrix("B", B)
end

* ===========================================================================
* Check 1: k=0 reduces to ols
* ===========================================================================
di as txt _newline(2) "============== Check 1: k=0 reduces to ols on (Theta, Q) =============="

* Build the topic-share variables Theta = W S' so we can run ols on the same
* design matrix that ols_bca_topic would form internally.
mata:
    Wm = st_matrix("W")
    Sm = st_matrix("S")
    Th = Wm * Sm'
    (void) st_addvar("double", ("theta_1", "theta_2", "theta_3"))
    st_store(., "theta_1", Th[, 1])
    st_store(., "theta_2", Th[, 2])
    st_store(., "theta_3", Th[, 3])
end

qui ols y theta_1 theta_2 theta_3 q1 q2
matrix b_ols = e(b)
matrix V_ols = e(V)

qui ols_bca_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0)
matrix b_k0 = e(b)
matrix V_k0 = e(V)

matrix diff_b = b_ols - b_k0
matrix diff_V = V_ols - V_k0
mata:
    db = st_matrix("diff_b")
    dV = st_matrix("diff_V")
    st_numscalar("max_db", max(abs(db)))
    st_numscalar("max_dV", max(abs(dV)))
end

local tol = 1e-12
if (max_db < `tol' & max_dV < `tol') {
    di as result "  Check 1 PASS  (max |db|=" %10.2e max_db ", max |dV|=" %10.2e max_dV ")"
}
else {
    di as error "  Check 1 FAIL  (max |db|=" %10.2e max_db ", max |dV|=" %10.2e max_dV ")"
    exit 9
}

* ===========================================================================
* Check 2: hand-built reference matches .ado on a nontrivial k
* ===========================================================================
di as txt _newline(2) "============== Check 2: formula parity (k=0.5) =============="

local kval = 0.5
qui ols_bca_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(`kval')
matrix b_ado = e(b)
matrix V_ado = e(V)

mata:
    Y  = st_data(., "y")
    Q1 = st_data(., ("q1", "q2"))
    Q1 = Q1, J(rows(Y), 1, 1)         // append intercept
    W  = st_matrix("W")
    S  = st_matrix("S")
    B  = st_matrix("B")
    k  = `kval'

    Theta = W * S'
    Xhat  = (Theta, Q1)
    n     = rows(Xhat)
    d     = cols(Xhat)
    r     = rows(S)

    XX_inv = invsym(quadcross(Xhat, Xhat))
    b0 = XX_inv * quadcross(Xhat, Y)
    u  = Y - Xhat * b0
    V0 = XX_inv * quadcross(Xhat, u :^ 2, Xhat) * XX_inv

    mW   = mean(W)'
    Bt   = B'
    BBt  = B * Bt
    BBti = invsym(BBt)
    Mmat = Bt :* (Bt * mW)
    Omega = S * BBti * B * Mmat * BBti * S' - quadcross(Theta, Theta) / n

    A = J(d, d, 0)
    A[|1, 1 \ r, r|] = Omega
    Gamma = (k * sqrt(n)) * XX_inv * A

    b_corr = (I(d) + Gamma) * b0
    V_corr = V0

    // names = topic_1, topic_2, topic_3, q1, q2, _cons; perm puts _cons first
    // then alphabetizes: _cons, q1, q2, topic_1, topic_2, topic_3
    // current order: 1=topic_1, 2=topic_2, 3=topic_3, 4=q1, 5=q2, 6=_cons
    // perm = (6, 4, 5, 1, 2, 3)
    perm = (6, 4, 5, 1, 2, 3)
    b_ref = b_corr[perm]
    V_ref = V_corr[perm, perm]

    b_ado = st_matrix("b_ado")'
    V_ado = st_matrix("V_ado")

    st_numscalar("max_db", max(abs(b_ado - b_ref)))
    st_numscalar("max_dV", max(abs(vec(V_ado - V_ref))))
end

local tol = 1e-10
if (max_db < `tol' & max_dV < `tol') {
    di as result "  Check 2 PASS  (max |db|=" %10.2e max_db ", max |dV|=" %10.2e max_dV ")"
}
else {
    di as error "  Check 2 FAIL  (max |db|=" %10.2e max_db ", max |dV|=" %10.2e max_dV ")"
    exit 9
}

* ===========================================================================
* Check 3: column names ordering
* ===========================================================================
di as txt _newline(2) "============== Check 3: e(b) column ordering =============="

qui ols_bca_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0.5)
local cnames : colnames e(b)
local expected "_cons q1 q2 topic_1 topic_2 topic_3"
if ("`cnames'" == "`expected'") {
    di as result "  Check 3 PASS  (`cnames')"
}
else {
    di as error "  Check 3 FAIL  got `cnames', expected `expected'"
    exit 9
}

di as result _newline "test_ols_bca_topic: all 3 checks PASS."
