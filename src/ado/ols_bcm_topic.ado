*! version 0.3.0  18may2026  Lorenzo Ranaldi
*! ValidMLInference (Stata port) -- ols_bcm_topic
*!
*! Multiplicative bias correction for OLS with topic-model-generated
*! regressors. Falls back to the additive (BCA) formula whenever the
*! spectral radius of Gamma is at least 1; otherwise inverts (I - Gamma).
*!
*! Same inputs as ols_bca_topic. Mirrors Python ValidMLInference.ols_bcm_topic.

program define ols_bcm_topic, eclass
    version 17

    syntax varlist(numeric min=1) [if] [in],     ///
        Wmatrix(name) Smatrix(name) Bmatrix(name) ///
        K(real)                                   ///
        [noCONStant]

    gettoken depvar covars : varlist

    marksample touse
    markout `touse' `depvar' `covars'

    qui count if `touse'
    if (r(N) == 0) error 2000
    local N = r(N)

    foreach mat in wmatrix smatrix bmatrix {
        cap confirm matrix ``mat''
        if _rc {
            di as error "matrix ``mat'' (option `mat'()) not found"
            exit 198
        }
    }

    local addconst = ("`constant'" != "noconstant")

    tempname b V
    mata: vmli_ols_bcm_topic_run("`depvar'", "`covars'", "`touse'",     ///
                                 "`wmatrix'", "`smatrix'", "`bmatrix'", ///
                                 `k', `addconst', "`b'", "`V'")

    ereturn post `b' `V', esample(`touse') depname(`depvar')
    ereturn scalar N        = `N'
    ereturn scalar k        = `k'
    ereturn local  cmd      "ols_bcm_topic"
    ereturn local  cmdline  `"ols_bcm_topic `0'"'
    ereturn local  title    "OLS with multiplicative topic-model bias correction"
    ereturn local  vce      "robust (HC0)"
    ereturn local  vcetype  "Robust"
    ereturn local  depvar   "`depvar'"
    ereturn local  W_matrix "`wmatrix'"
    ereturn local  S_matrix "`smatrix'"
    ereturn local  B_matrix "`bmatrix'"
    ereturn local  properties "b V"

    di _newline as txt "`e(title)' (ValidMLInference)"
    di as txt "Observations: " as result %9.0gc e(N) ///
       _col(40) as txt "k: " as result %9.4f e(k)
    di
    ereturn display
end
