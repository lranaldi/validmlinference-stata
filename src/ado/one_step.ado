*! version 0.6.0  26may2026  Lorenzo Ranaldi
*! ValidMLInference (Stata port) -- one_step
*!
*! Joint MLE for an outcome model with a noisy/AI-generated binary
*! regressor. The latent treatment is integrated out via a 2x2
*! misclassification table; residual variance can be heteroskedastic in
*! the latent treatment (default) or homoskedastic. The residual
*! distribution is chosen from a closed menu via dist().
*!
*! Mirrors Python ValidMLInference.one_step. The Python `distribution=`
*! callback is exposed in Stata as a closed menu: dist(normal|laplace|t).
*! For dist(t) the user-supplied degrees of freedom df(#) is fixed (not
*! estimated), matching Python's pattern (user picks the callable).
*!
*! Coefficient ordering posted in e(b): intercept first, then non-intercept
*! names alphabetically (matches Python `_standardize_coefficient_order`).

program define one_step, eclass
    version 17

    syntax varlist(numeric min=2) [if] [in], ///
        [GENerated(varname)                  ///
         HOMOSKedastic                       ///
         DIST(string)                        ///
         DF(real 0)                          ///
         noCONStant]

    gettoken depvar indepvars : varlist

    marksample touse
    markout `touse' `depvar' `indepvars'

    qui count if `touse'
    if (r(N) == 0) error 2000
    local N = r(N)

    local addconst = ("`constant'" != "noconstant")
    local homosk   = ("`homoskedastic'" != "")

    * --- distribution menu ---
    if ("`dist'" == "") local dist "normal"
    local dist = strlower("`dist'")
    if      ("`dist'" == "normal")  local distcode = 1
    else if ("`dist'" == "laplace") local distcode = 2
    else if ("`dist'" == "t")       local distcode = 3
    else {
        di as error `"dist(`dist') not recognized; allowed: normal, laplace, t"'
        exit 198
    }

    if (`distcode' == 3) {
        if (`df' <= 0) {
            di as error "dist(t) requires df(#) with df > 0"
            exit 198
        }
    }
    else {
        if (`df' != 0) {
            di as error "df() is only used with dist(t); drop df() or set dist(t)"
            exit 198
        }
    }

    * --- generated/treatment variable: default to first covariate ---
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
    mata: vmli_one_step_run("`depvar'", "`indepvars'", "`touse'",   ///
                            `addconst', "`target'", `homosk',       ///
                            `distcode', `df',                       ///
                            "`b'", "`V'")

    * --- capture the auxiliary scalars Mata stashed in r() before ereturn ---
    local converged  = r(converged)
    local iterations = r(iterations)
    local ll         = r(ll)

    local distlabel = strproper("`dist'")
    if ("`dist'" == "t") {
        local distlabel = "Student-t(df=`df')"
        local title "One-step joint MLE (`distlabel')"
    }
    else {
        local title "One-step joint MLE (`distlabel')"
    }

    ereturn post `b' `V', esample(`touse') depname(`depvar')
    ereturn scalar N             = `N'
    ereturn scalar converged     = `converged'
    ereturn scalar iterations    = `iterations'
    ereturn scalar ll            = `ll'
    ereturn scalar homoskedastic = `homosk'
    ereturn scalar distcode      = `distcode'
    if (`distcode' == 3) ereturn scalar df = `df'
    ereturn local  dist          "`dist'"
    ereturn local  generated_var "`target'"
    ereturn local  cmd           "one_step"
    ereturn local  cmdline       `"one_step `0'"'
    ereturn local  title         "`title'"
    ereturn local  vce           "OIM (numerical Hessian)"
    ereturn local  vcetype       "OIM"
    ereturn local  depvar        "`depvar'"
    ereturn local  properties    "b V"

    di _newline as txt "`e(title)' (ValidMLInference)"
    di as txt "Observations: "        as result %9.0gc e(N)               ///
       _col(40) as txt "Iterations: " as result %4.0f  e(iterations)
    di as txt "Treatment:    "        as result "`target'"                ///
       _col(40) as txt "Log-lik:    " as result %10.4f e(ll)
    di as txt "Hetero. var:  "        as result cond(`homosk', "no", "yes") ///
       _col(40) as txt "Residual:   " as result "`distlabel'"
    if (!e(converged)) {
        di as error  "WARNING: optimizer did not converge."
    }
    di
    ereturn display
end
