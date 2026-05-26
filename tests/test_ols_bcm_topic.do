*! test_ols_bcm_topic.do -- Milestone 4 parity check for ols_bcm_topic
*!
*! Four checks:
*!   1. k=0 sanity: reduces to ols.
*!   2. Formula check (rho<1 branch): (I-Gamma)^{-1} b0 matches a hand-built
*!      Mata reference.
*!   3. Formula check (rho>=1 branch): falls back to the BCA formula
*!      (I+Gamma) b0; should equal ols_bca_topic output with the same inputs.
*!   4. Column ordering: intercept-first, then alphabetical.

version 17
set more off
clear all

* --- Locate project root ---
local pkg "ValidMLInference-stata"
local marker "`pkg'/src/ado/ols_bcm_topic.ado"
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
    di as error "test_ols_bcm_topic.do: cannot locate `marker' from cwd"
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
cap program drop ols_bcm_topic

* ---------------------------------------------------------------------------
* Synthetic topic-model data (same recipe as test_ols_bca_topic)
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
    Y = (Theta * bt + Q * bq + rnormal(nn, 1, 0, 1)) :+ 0.7
    (void) st_addvar("double", ("y", "q1", "q2", "theta_1", "theta_2", "theta_3"))
    st_store(., "y",       Y)
    st_store(., "q1",      Q[, 1])
    st_store(., "q2",      Q[, 2])
    st_store(., "theta_1", Theta[, 1])
    st_store(., "theta_2", Theta[, 2])
    st_store(., "theta_3", Theta[, 3])
    st_matrix("W", W)
    st_matrix("S", S)
    st_matrix("B", B)
end

* ===========================================================================
* Check 1: k=0 reduces to ols
* ===========================================================================
di as txt _newline(2) "============== Check 1: k=0 reduces to ols on (Theta, Q) =============="

qui ols y theta_1 theta_2 theta_3 q1 q2
matrix b_ols = e(b)
matrix V_ols = e(V)

qui ols_bcm_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0)
matrix b_k0 = e(b)
matrix V_k0 = e(V)

matrix diff_b = b_ols - b_k0
matrix diff_V = V_ols - V_k0
mata:
    st_numscalar("max_db", max(abs(st_matrix("diff_b"))))
    st_numscalar("max_dV", max(abs(st_matrix("diff_V"))))
end

local tol = 1e-12
if (max_db < `tol' & max_dV < `tol') {
    di as result "  Check 1 PASS"
}
else {
    di as error "  Check 1 FAIL  (max |db|=" %10.2e max_db ")"
    exit 9
}

* ===========================================================================
* Check 2: rho < 1 branch matches (I - Gamma)^{-1} b0
* ===========================================================================
di as txt _newline(2) "============== Check 2: rho<1 branch (small k) =============="

* With this synthetic recipe rho_raw(XX_inv*A) ~ 4.6, so we need
* k * sqrt(n) * rho_raw < 1 -> k < ~0.013. k=0.005 keeps rho safely below 1.
local kval = 0.005
qui ols_bcm_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(`kval')
matrix b_ado = e(b)
matrix V_ado = e(V)

mata:
    Y = st_data(., "y")
    Q1 = st_data(., ("q1", "q2"))
    Q1 = Q1, J(rows(Y), 1, 1)
    W = st_matrix("W")
    S = st_matrix("S")
    B = st_matrix("B")
    k = `kval'
    Theta = W * S'
    Xhat = (Theta, Q1)
    n = rows(Xhat)
    d = cols(Xhat)
    r = rows(S)
    XX_inv = invsym(quadcross(Xhat, Xhat))
    b0 = XX_inv * quadcross(Xhat, Y)
    u  = Y - Xhat * b0
    V0 = XX_inv * quadcross(Xhat, u :^ 2, Xhat) * XX_inv
    mW = mean(W)'
    Bt = B'
    BBti = invsym(B * Bt)
    Mmat = Bt :* (Bt * mW)
    Omega = S * BBti * B * Mmat * BBti * S' - quadcross(Theta, Theta) / n
    A = J(d, d, 0)
    A[|1, 1 \ r, r|] = Omega
    Gamma = (k * sqrt(n)) * XX_inv * A
    rho = max(abs(eigenvalues(Gamma)))
    st_numscalar("rho_small", rho)
    b_corr = lusolve(I(d) - Gamma, b0)
    V_corr = V0
    perm = (6, 4, 5, 1, 2, 3)
    b_ref = b_corr[perm]
    V_ref = V_corr[perm, perm]
    b_ado = st_matrix("b_ado")'
    V_ado = st_matrix("V_ado")
    st_numscalar("max_db", max(abs(b_ado - b_ref)))
    st_numscalar("max_dV", max(abs(vec(V_ado - V_ref))))
end

di as txt "  rho(Gamma) at k=`kval' = " %8.4f rho_small
local tol = 1e-10
if (rho_small < 1 & max_db < `tol' & max_dV < `tol') {
    di as result "  Check 2 PASS  (rho<1 branch, parity within `tol')"
}
else {
    di as error  "  Check 2 FAIL  (rho=" %6.4f rho_small ", max |db|=" %10.2e max_db ")"
    exit 9
}

* ===========================================================================
* Check 3: rho >= 1 branch falls back to BCA formula
*   We pick a large k that pushes rho(Gamma) >= 1, then verify that
*   ols_bcm_topic output equals ols_bca_topic output at the same k.
* ===========================================================================
di as txt _newline(2) "============== Check 3: rho>=1 branch falls back to BCA =============="

local kval = 20
qui ols_bca_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(`kval')
matrix b_bca = e(b)

qui ols_bcm_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(`kval')
matrix b_bcm = e(b)

mata:
    Y = st_data(., "y")
    Q1 = st_data(., ("q1", "q2"))
    Q1 = Q1, J(rows(Y), 1, 1)
    W = st_matrix("W")
    S = st_matrix("S")
    B = st_matrix("B")
    Theta = W * S'
    Xhat = (Theta, Q1)
    n = rows(Xhat)
    d = cols(Xhat)
    r = rows(S)
    XX_inv = invsym(quadcross(Xhat, Xhat))
    mW = mean(W)'
    Bt = B'
    BBti = invsym(B * Bt)
    Mmat = Bt :* (Bt * mW)
    Omega = S * BBti * B * Mmat * BBti * S' - quadcross(Theta, Theta) / n
    A = J(d, d, 0)
    A[|1, 1 \ r, r|] = Omega
    Gamma = (`kval' * sqrt(n)) * XX_inv * A
    rho = max(abs(eigenvalues(Gamma)))
    st_numscalar("rho_k20", rho)
    b_bca = st_matrix("b_bca")
    b_bcm = st_matrix("b_bcm")
    st_numscalar("max_db", max(abs(b_bca - b_bcm)))
end

di as txt "  rho(Gamma) at k=`kval' = " %8.4f rho_k20
local tol = 1e-12
if (rho_k20 >= 1 & max_db < `tol') {
    di as result "  Check 3 PASS  (BCM == BCA when rho>=1)"
}
else {
    di as error  "  Check 3 FAIL  (rho=" %6.4f rho_k20 ", max |b_bca - b_bcm|=" %10.2e max_db ")"
    exit 9
}

* ===========================================================================
* Check 4: column ordering
* ===========================================================================
di as txt _newline(2) "============== Check 4: e(b) column ordering =============="

qui ols_bcm_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0.5)
local cnames : colnames e(b)
local expected "_cons q1 q2 topic_1 topic_2 topic_3"
if ("`cnames'" == "`expected'") {
    di as result "  Check 4 PASS  (`cnames')"
}
else {
    di as error "  Check 4 FAIL  got `cnames', expected `expected'"
    exit 9
}

di as result _newline "test_ols_bcm_topic: all 4 checks PASS."
