*! version 0.5.0  22may2026  Lorenzo Ranaldi
*! ValidMLInference (Stata port) -- one_step_gmm
*!
*! Joint MLE for a Gaussian-mixture outcome model where one regressor is a
*! noisy/AI-generated binary label. The per-class residual is a k-component
*! Gaussian mixture (heteroskedastic across the latent treatment).
*!
*! Mirrors Python ValidMLInference.one_step_gaussian_mixture, with the
*! Stata-friendly command name. The Python `distribution=` callback is not
*! exposed (upstream's GMM public interface does not accept it -- the
*! component density is hard-coded to Gaussian via `mixture_pdf`). The
*! `homoskedastic` option is rejected: the upstream Python homoskedastic
*! branch has a theta-length mismatch (starts has d+1+3k slots but
*! unpack_theta reads d-1+5k). See Notes/porting_decisions.md.
*!
*! Coefficient ordering posted in e(b): intercept first, then non-intercept
*! names alphabetically (matches Python `_standardize_coefficient_order`).

program define one_step_gmm, eclass
    version 17

    syntax varlist(numeric min=2) [if] [in], ///
        [GENerated(varname)                  ///
         K(integer 2)                        ///
         HOMOSKedastic                       ///
         NGuess(integer 10)                  ///
         MAXiter(integer 100)                ///
         SEED(integer 0)                     ///
         noCONStant]

    gettoken depvar indepvars : varlist

    marksample touse
    markout `touse' `depvar' `indepvars'

    qui count if `touse'
    if (r(N) == 0) error 2000
    local N = r(N)

    if (`k' < 2) {
        di as error "k(`k') must be >= 2"
        exit 198
    }
    if (`nguess' < 1) {
        di as error "nguess(`nguess') must be >= 1"
        exit 198
    }
    if (`maxiter' < 1) {
        di as error "maxiter(`maxiter') must be >= 1"
        exit 198
    }

    local addconst = ("`constant'" != "noconstant")
    local homosk   = ("`homoskedastic'" != "")

    if (`homosk') {
        di as error "homoskedastic Gaussian mixture not currently supported"
        di as error "(upstream Python theta-length mismatch). Use the default heteroskedastic option."
        exit 198
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
    mata: vmli_gmm_run("`depvar'", "`indepvars'", "`touse'",        ///
                      `addconst', "`target'", `k', `homosk',        ///
                      `nguess', `maxiter', `seed',                  ///
                      "`b'", "`V'")

    local converged  = r(converged)
    local iterations = r(iterations)
    local ll         = r(ll)
    local best_idx   = r(best_idx)
    local n_finished = r(n_finished)

    ereturn post `b' `V', esample(`touse') depname(`depvar')
    ereturn scalar N             = `N'
    ereturn scalar k             = `k'
    ereturn scalar nguess        = `nguess'
    ereturn scalar maxiter       = `maxiter'
    ereturn scalar seed          = `seed'
    ereturn scalar converged     = `converged'
    ereturn scalar iterations    = `iterations'
    ereturn scalar ll            = `ll'
    ereturn scalar best_idx      = `best_idx'
    ereturn scalar n_finished    = `n_finished'
    ereturn scalar homoskedastic = `homosk'
    ereturn local  generated_var "`target'"
    ereturn local  cmd           "one_step_gmm"
    ereturn local  cmdline       `"one_step_gmm `0'"'
    ereturn local  title         "One-step joint MLE (Gaussian mixture, k=`k')"
    ereturn local  vce           "OIM (numerical Hessian)"
    ereturn local  vcetype       "OIM"
    ereturn local  depvar        "`depvar'"
    ereturn local  properties    "b V"

    di _newline as txt "`e(title)' (ValidMLInference)"
    di as txt "Observations: "        as result %9.0gc e(N)               ///
       _col(40) as txt "Components: " as result %4.0f  e(k)
    di as txt "Treatment:    "        as result "`target'"                ///
       _col(40) as txt "Log-lik:    " as result %10.4f e(ll)
    di as txt "Multistart:   "        as result %3.0f `n_finished'        ///
                                      as txt   "/"                        ///
                                      as result %-3.0f `nguess'           ///
       _col(40) as txt "Best try:   " as result %4.0f  e(best_idx)        ///
       _col(56) as txt "Iter (best):" as result %5.0f  e(iterations)
    if (!e(converged)) {
        di as error  "WARNING: best optimizer attempt did not converge."
    }
    di
    ereturn display
end
