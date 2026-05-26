*! test_ols_bca.do -- Milestone 3 parity check for ols_bca
*!
*! Three checks:
*!   1. fpr=0 sanity: ols_bca should reproduce ols exactly.
*!   2. Formula check: a hand-built Mata reference reproduces the
*!      .ado output to machine precision on a nontrivial (fpr, m).
*!   3. generated() option: switching the target column changes the
*!      result and the command still runs without error.

version 17
set more off
clear all

* --- Locate project root ---
local pkg "ValidMLInference-stata"
local marker "`pkg'/src/ado/ols_bca.ado"
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
    di as error "test_ols_bca.do: cannot locate `marker' from cwd"
    exit 601
}
local proj "`c(pwd)'"

* --- Build mlib + set adopath ---
qui cd "`pkg'/src/mata"
do _build.do
qui cd "`proj'"

adopath ++ "`pkg'/src/ado"
adopath ++ "`pkg'/src/mata"
mata: mata mlib index

cap program drop ols
cap program drop ols_bca

* --- Synthetic data with a noisy binary covariate ---
clear
set seed 20260518
set obs 2000
gen byte x_true = uniform() > 0.4
gen byte x1 = x_true
* Add classification noise: ~10% false positives, ~5% false negatives
replace x1 = 1 if x_true == 0 & uniform() < 0.10
replace x1 = 0 if x_true == 1 & uniform() < 0.05
gen x2 = rnormal()
gen x3 = rnormal()
gen y  = 1 + 0.5*x_true + 2*x2 - 0.3*x3 + rnormal() * (1 + 0.5*abs(x2))

* ===========================================================================
* Check 1: fpr=0 sanity
* ===========================================================================
di as txt _newline(2) "============== Check 1: fpr=0 reduces to ols =============="

qui ols y x1 x2 x3
matrix b_ols  = e(b)
matrix V_ols  = e(V)

qui ols_bca y x1 x2 x3, fpr(0) m(100) generated(x1)
matrix b_bca0 = e(b)
matrix V_bca0 = e(V)

matrix diff_b = b_ols - b_bca0
matrix diff_V = V_ols - V_bca0

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
    di as error  "  Check 1 FAIL (delta b = " %10.2e max_db ", delta V = " %10.2e max_dV ")"
    exit 9
}

* ===========================================================================
* Check 2: hand-built Mata reference matches ols_bca on a nontrivial case
* ===========================================================================
di as txt _newline(2) "============== Check 2: formula parity (fpr=0.07, m=500) =============="

local fpr  = 0.07
local m    = 500
qui ols_bca y x1 x2 x3, fpr(`fpr') m(`m') generated(x1)
matrix b_ado = e(b)
matrix V_ado = e(V)

mata:
    // Read data from current sample
    Y = st_data(., "y")
    X = st_data(., ("x1", "x2", "x3"))
    X = X, J(rows(X), 1, 1)            // append intercept at end
    n = rows(X)

    // HC0 OLS
    XX_inv = invsym(quadcross(X, X))
    b0 = XX_inv * quadcross(X, Y)
    u  = Y - X * b0
    V0 = XX_inv * quadcross(X, u :^ 2, X) * XX_inv

    // Gamma at target_idx = 1 (= x1)
    d = cols(X)
    A = J(d, d, 0)
    A[1, 1] = 1
    Gamma = n * XX_inv * A

    fpr = `fpr'
    m   = `m'

    // BCA correction (formula straight from Python _ols_bca_core)
    b_corr = b0 + fpr * (Gamma * b0)
    Id     = I(d)
    M      = Id + fpr * Gamma
    extra  = (fpr * (1 - fpr) / m) * (Gamma * (V0 + b_corr * b_corr') * Gamma')
    V_corr = M * V0 * M' + extra

    // Reorder: intercept-first then alphabetical
    //   design names = (x1, x2, x3, _cons); perm = (_cons, x1, x2, x3) = (4,1,2,3)
    perm = (4, 1, 2, 3)
    b_ref = b_corr[perm]
    V_ref = V_corr[perm, perm]

    st_matrix("b_ref", b_ref')
    st_matrix("V_ref", V_ref)

    // Compare to .ado output
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

* Run with default target (x1) and with target=x2; the b/V should differ.
qui ols_bca y x1 x2 x3, fpr(0.05) m(500)
matrix b_def = e(b)

qui ols_bca y x1 x2 x3, fpr(0.05) m(500) generated(x2)
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
    di as error  "  Check 3 FAIL (target swap had no effect)"
    exit 9
}

di as result _newline "test_ols_bca: all 3 checks PASS."
