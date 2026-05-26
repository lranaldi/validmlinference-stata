*! test_ols_bcm.do -- Milestone 3 parity check for ols_bcm
*!
*! Three checks, parallel to test_ols_bca.do:
*!   1. fpr=0 sanity: ols_bcm should reproduce ols exactly.
*!   2. Formula check: a hand-built Mata reference reproduces the
*!      .ado output to machine precision on a nontrivial (fpr, m).
*!   3. generated() option: switching the target column changes the
*!      result and the command still runs.

version 17
set more off
clear all

* --- Locate project root ---
local pkg "ValidMLInference-stata"
local marker "`pkg'/src/ado/ols_bcm.ado"
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
    di as error "test_ols_bcm.do: cannot locate `marker' from cwd"
    exit 601
}
local proj "`c(pwd)'"

qui cd "`pkg'/src/mata"
do _build.do
qui cd "`proj'"

adopath ++ "`pkg'/src/ado"
adopath ++ "`pkg'/src/mata"
mata: mata mlib index

cap program drop ols
cap program drop ols_bcm

* --- Synthetic data with a noisy binary covariate ---
clear
set seed 20260518
set obs 2000
gen byte x_true = uniform() > 0.4
gen byte x1 = x_true
replace x1 = 1 if x_true == 0 & uniform() < 0.10
replace x1 = 0 if x_true == 1 & uniform() < 0.05
gen x2 = rnormal()
gen x3 = rnormal()
gen y  = 1 + 0.5*x_true + 2*x2 - 0.3*x3 + rnormal() * (1 + 0.5*abs(x2))

* ===========================================================================
* Check 1: fpr=0 reduces to ols
* ===========================================================================
di as txt _newline(2) "============== Check 1: fpr=0 reduces to ols =============="

qui ols y x1 x2 x3
matrix b_ols  = e(b)
matrix V_ols  = e(V)

qui ols_bcm y x1 x2 x3, fpr(0) m(100) generated(x1)
matrix b_bcm0 = e(b)
matrix V_bcm0 = e(V)

matrix diff_b = b_ols - b_bcm0
matrix diff_V = V_ols - V_bcm0

mata:
    db = st_matrix("diff_b")
    dV = st_matrix("diff_V")
    max_db = max(abs(db))
    max_dV = max(abs(dV))
    printf("  max |delta b|  = %12.4e\n", max_db)
    printf("  max |delta V|  = %12.4e\n", max_dV)
    st_numscalar("max_db", max_db)
    st_numscalar("max_dV", max_dV)
end

local tol = 1e-12
if (max_db < `tol' & max_dV < `tol') {
    di as result "  Check 1 PASS"
}
else {
    di as error  "  Check 1 FAIL"
    exit 9
}

* ===========================================================================
* Check 2: hand-built Mata reference matches .ado output
* ===========================================================================
di as txt _newline(2) "============== Check 2: formula parity (fpr=0.07, m=500) =============="

local fpr  = 0.07
local m    = 500
qui ols_bcm y x1 x2 x3, fpr(`fpr') m(`m') generated(x1)
matrix b_ado = e(b)
matrix V_ado = e(V)

mata:
    Y = st_data(., "y")
    X = st_data(., ("x1", "x2", "x3"))
    X = X, J(rows(X), 1, 1)
    n = rows(X)

    XX_inv = invsym(quadcross(X, X))
    b0 = XX_inv * quadcross(X, Y)
    u  = Y - X * b0
    V0 = XX_inv * quadcross(X, u :^ 2, X) * XX_inv

    d = cols(X)
    A = J(d, d, 0)
    A[1, 1] = 1
    Gamma = n * XX_inv * A

    fpr = `fpr'
    m   = `m'

    // BCM correction (formula straight from Python _ols_bcm_core)
    Id    = I(d)
    M     = Id - fpr * Gamma
    M_inv = luinv(M)
    b_corr = M_inv * b0
    extra  = (fpr * (1 - fpr) / m) * (Gamma * (V0 + b_corr * b_corr') * Gamma')
    V_corr = M_inv * V0 * M_inv' + extra

    perm = (4, 1, 2, 3)
    b_ref = b_corr[perm]
    V_ref = V_corr[perm, perm]

    b_ado = st_matrix("b_ado")'
    V_ado = st_matrix("V_ado")

    max_db = max(abs(b_ado - b_ref))
    max_dV = max(abs(vec(V_ado - V_ref)))

    printf("  max |delta b|  = %12.4e\n", max_db)
    printf("  max |delta V|  = %12.4e\n", max_dV)

    st_numscalar("max_db", max_db)
    st_numscalar("max_dV", max_dV)
end

local tol = 1e-10
if (max_db < `tol' & max_dV < `tol') {
    di as result "  Check 2 PASS (within `tol' of hand-built reference)"
}
else {
    di as error  "  Check 2 FAIL"
    exit 9
}

* ===========================================================================
* Check 3: generated() option swaps the target column
* ===========================================================================
di as txt _newline(2) "============== Check 3: generated() targets a different column =============="

qui ols_bcm y x1 x2 x3, fpr(0.05) m(500)
matrix b_def = e(b)

qui ols_bcm y x1 x2 x3, fpr(0.05) m(500) generated(x2)
matrix b_x2  = e(b)

matrix diff_b = b_def - b_x2
mata:
    db = st_matrix("diff_b")
    max_db = max(abs(db))
    printf("  max |b_default - b_generated(x2)| = %12.4e\n", max_db)
    st_numscalar("max_db", max_db)
end

if (max_db > 1e-6) {
    di as result "  Check 3 PASS (target column swap changes the estimate, as expected)"
}
else {
    di as error  "  Check 3 FAIL"
    exit 9
}

di as result _newline "test_ols_bcm: all 3 checks PASS."
