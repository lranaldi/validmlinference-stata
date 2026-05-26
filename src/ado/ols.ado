*! version 0.1.0  18may2026  Lorenzo Ranaldi
*! ValidMLInference (Stata port) -- ols
*!
*! OLS with HC0 (no df adjustment) robust standard errors. Mirrors the
*! Python ValidMLInference.ols function exactly. Differs from Stata's
*! `regress ..., robust`, which applies HC1's (n)/(n-k) inflation.
*!
*! Coefficient ordering posted in e(b): intercept first, then non-intercept
*! names alphabetically -- matching Python's _standardize_coefficient_order
*! so parity-test fixtures align.

program define ols, eclass
    version 17

    syntax varlist(numeric min=2) [if] [in] [, noCONStant]

    gettoken depvar indepvars : varlist

    marksample touse
    markout `touse' `depvar' `indepvars'

    qui count if `touse'
    if (r(N) == 0) error 2000
    local N = r(N)

    local addconst = ("`constant'" != "noconstant")

    tempname b V
    mata: vmli_ols_run("`depvar'", "`indepvars'", "`touse'", ///
                       `addconst', "`b'", "`V'")

    ereturn post `b' `V', esample(`touse') depname(`depvar')
    ereturn scalar N         = `N'
    ereturn local  cmd       "ols"
    ereturn local  cmdline   `"ols `0'"'
    ereturn local  title     "OLS with HC0 standard errors"
    ereturn local  vce       "robust (HC0)"
    ereturn local  vcetype   "Robust"
    ereturn local  depvar    "`depvar'"
    ereturn local  properties "b V"

    di _newline as txt "`e(title)' (ValidMLInference)"
    di as txt "Observations: " as result %9.0gc e(N)
    di
    ereturn display
end
