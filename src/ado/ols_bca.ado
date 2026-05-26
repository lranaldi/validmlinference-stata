*! version 0.2.0  18may2026  Lorenzo Ranaldi
*! ValidMLInference (Stata port) -- ols_bca
*!
*! Additive bias correction for OLS when one covariate is an AI/ML-generated
*! binary label. Requires an estimate of the classifier's false-positive
*! rate (fpr) and the sample size used to estimate it (m). Standard errors
*! account for uncertainty in the fpr estimate.
*!
*! Coefficient ordering in e(b): intercept first, then non-intercept names
*! alphabetically (matches Python's _standardize_coefficient_order).

program define ols_bca, eclass
    version 17

    syntax varlist(numeric min=2) [if] [in], ///
        FPR(real) M(real)                    ///
        [GENerated(varname) noCONStant]

    gettoken depvar indepvars : varlist

    marksample touse
    markout `touse' `depvar' `indepvars'

    qui count if `touse'
    if (r(N) == 0) error 2000
    local N = r(N)

    * --- option validation ---
    if (`fpr' < 0 | `fpr' > 1) {
        di as error "fpr must be in [0, 1]; got `fpr'"
        exit 198
    }
    if (`m' <= 0) {
        di as error "m must be positive; got `m'"
        exit 198
    }

    local addconst = ("`constant'" != "noconstant")

    * --- generated variable ---
    *   Default: first variable in the varlist.
    if ("`generated'" == "") {
        gettoken target : indepvars
    }
    else {
        local target "`generated'"
        local found 0
        foreach v of local indepvars {
            if ("`v'" == "`target'") local found 1
        }
        if (!`found') {
            di as error "generated(`target') not found in varlist (`indepvars')"
            exit 111
        }
    }

    tempname b V
    mata: vmli_ols_bca_run("`depvar'", "`indepvars'", "`touse'",   ///
                           `addconst', `fpr', `m', "`target'",     ///
                           "`b'", "`V'")

    ereturn post `b' `V', esample(`touse') depname(`depvar')
    ereturn scalar N             = `N'
    ereturn scalar fpr           = `fpr'
    ereturn scalar m             = `m'
    ereturn local  generated_var "`target'"
    ereturn local  cmd           "ols_bca"
    ereturn local  cmdline       `"ols_bca `0'"'
    ereturn local  title         "OLS with additive bias correction"
    ereturn local  vce           "robust (HC0) + bias-correction"
    ereturn local  vcetype       "Robust"
    ereturn local  depvar        "`depvar'"
    ereturn local  properties    "b V"

    di _newline as txt "`e(title)' (ValidMLInference)"
    di as txt "Observations: " as result %9.0gc e(N)              ///
       _col(40) as txt "fpr: " as result %7.4f e(fpr)
    di as txt "Generated variable: " as result "`target'"         ///
       _col(40) as txt "m:   " as result %9.0g e(m)
    di
    ereturn display
end
