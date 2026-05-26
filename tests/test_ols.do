*! test_ols.do -- Milestone 2 parity check for ols
*!
*! Generates a synthetic regression dataset, runs both our `ols` and
*! Stata's native `regress ..., vce(hc0)`, and compares them. The two
*! must agree to within floating-point noise on every coefficient and
*! standard error because both implement the same HC0 estimator on the
*! same data.

version 17
set more off
clear all

* ---------------------------------------------------------------------------
* Locate project root
* ---------------------------------------------------------------------------
local pkg "ValidMLInference-stata"
local marker "`pkg'/src/ado/ols.ado"
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
    di as error "test_ols.do: cannot locate `marker' from cwd"
    exit 601
}
local proj "`c(pwd)'"

* ---------------------------------------------------------------------------
* (Re)build the Mata library so we test the current source
* ---------------------------------------------------------------------------
qui cd "`pkg'/src/mata"
do _build.do
qui cd "`proj'"

adopath ++ "`pkg'/src/ado"
adopath ++ "`pkg'/src/mata"
mata: mata mlib index

* Clear any prior in-memory definition of `ols` so we pick up the new ado.
cap program drop ols

* ---------------------------------------------------------------------------
* Synthetic data: y = 1 + 0.5*x1 + 2*x2 - 0.3*x3 + heteroskedastic noise
* (heteroskedasticity is what makes HC0 differ meaningfully from OLS-classic)
* ---------------------------------------------------------------------------
clear
set seed 20260518
set obs 1000
gen x1 = rnormal()
gen x2 = rnormal()
gen x3 = rnormal()
gen e  = rnormal() * (1 + 0.5 * abs(x2))
gen y  = 1 + 0.5*x1 + 2*x2 - 0.3*x3 + e

* ---------------------------------------------------------------------------
* Run our ols
* ---------------------------------------------------------------------------
di as txt _newline(2) "============== ols y x1 x2 x3 =============="
ols y x1 x2 x3

matrix b_ols = e(b)
matrix V_ols = e(V)

* ---------------------------------------------------------------------------
* Run native regress, robust (HC1), then scale variance down to HC0.
*   Stata's regress lacks a direct HC0 option, but
*       V_HC1 = (n/(n-k)) * V_HC0   =>   V_HC0 = (df_r / n) * V_HC1
* ---------------------------------------------------------------------------
di as txt _newline(2) "============== regress y x1 x2 x3, robust (HC1) =============="
regress y x1 x2 x3, robust

local n_reg  = e(N)
local df_r   = e(df_r)
local scale  = `df_r' / `n_reg'

matrix b_reg = e(b)
matrix V_reg = `scale' * e(V)    // converted to HC0 for parity comparison

* ---------------------------------------------------------------------------
* Compare coefficient by coefficient (lookup by name, so ordering is moot)
* ---------------------------------------------------------------------------
di as txt _newline "============== Parity vs. regress, robust (rescaled to HC0) =============="
local maxdiff_b  = 0
local maxdiff_se = 0
foreach v in _cons x1 x2 x3 {
    local b_o  = b_ols[1, "`v'"]
    local b_r  = b_reg[1, "`v'"]
    local se_o = sqrt(V_ols["`v'", "`v'"])
    local se_r = sqrt(V_reg["`v'", "`v'"])

    local d_b  = abs(`b_o'  - `b_r')
    local d_se = abs(`se_o' - `se_r')

    if (`d_b'  > `maxdiff_b')  local maxdiff_b  = `d_b'
    if (`d_se' > `maxdiff_se') local maxdiff_se = `d_se'

    di as txt "  `v':" _col(12)            ///
       "b_ols = " %14.10f `b_o'            ///
       "  b_reg = " %14.10f `b_r'          ///
       "  |delta b| = "  %10.2e `d_b'      ///
       "  |delta se| = " %10.2e `d_se'
}

di as txt _newline "max |delta b|  = " %10.2e `maxdiff_b'
di as txt          "max |delta se| = " %10.2e `maxdiff_se'

local tol = 1e-10
if (`maxdiff_b' < `tol' & `maxdiff_se' < `tol') {
    di as result _newline "test_ols: PASS (within `tol' of regress, vce(hc0))."
}
else {
    di as error _newline "test_ols: FAIL"
    di as error "  max |delta b|  = " %10.2e `maxdiff_b'  " (tol `tol')"
    di as error "  max |delta se| = " %10.2e `maxdiff_se' " (tol `tol')"
    exit 9
}

* ---------------------------------------------------------------------------
* Sanity: check that e(b) column order is intercept-first then alphabetical
* ---------------------------------------------------------------------------
qui ols y x1 x2 x3
local expected "_cons x1 x2 x3"
local cnames : colnames e(b)
if ("`cnames'" == "`expected'") {
    di as result "test_ols: e(b) column order OK (`cnames')"
}
else {
    di as error "test_ols: e(b) column order is `cnames', expected `expected'"
    exit 9
}

* ---------------------------------------------------------------------------
* Sanity: noconstant
* ---------------------------------------------------------------------------
di as txt _newline "============== ols y x1 x2 x3, noconstant =============="
ols y x1 x2 x3, noconstant
local cnames_nc : colnames e(b)
if (strpos("`cnames_nc'", "_cons") == 0) {
    di as result "test_ols: noconstant correctly drops intercept (`cnames_nc')"
}
else {
    di as error "test_ols: noconstant did NOT drop _cons (`cnames_nc')"
    exit 9
}
