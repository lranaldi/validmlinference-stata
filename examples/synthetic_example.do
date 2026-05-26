// examples/synthetic_example.do
//
// Simulation example from Battaglia, Christensen, Hansen & Sacher (2025).
// Mirrors `ValidMLInference-main/synthetic_example.ipynb`. Compares four
// estimators on simulated data where a binary covariate is observed with
// classification noise:
//     1. OLS on the noisy label (biased)
//     2. ols_bca (additive bias correction)
//     3. ols_bcm (multiplicative bias correction)
//     4. one_step (joint MLE)
//
// The classification structure: for each obs, draw u ~ U(0,1) and set
//   - u in [0, fpr]:           X=1, Xhat=0  (false negative)
//   - u in (fpr, 2*fpr]:       X=0, Xhat=1  (false positive)
//   - u in (2*fpr, p+fpr]:     X=1, Xhat=1  (true positive)
//   - u in (p+fpr, 1]:         X=0, Xhat=0  (true negative)
// then Y = b0 + b1*X + (s1*X + s0*(1-X)) * eps,  eps ~ N(0,1).
//
// Run from the project root. Optionally override the number of simulations:
//     do "ValidMLInference-stata/examples/synthetic_example.do"          // nsim=1000
//     do "ValidMLInference-stata/examples/synthetic_example.do" 50       // nsim=50 (smoke test)

clear all
version 17
set more off

args nsim_arg
local nsim = cond("`nsim_arg'" == "", 1000, real("`nsim_arg'"))

local proj "`c(pwd)'"
capture confirm file "`proj'/ValidMLInference-stata/src/ado/ols.ado"
if _rc {
    local proj "C:/Users/loren/Dropbox/Replication_AIValidInference"
}
adopath ++ "`proj'/ValidMLInference-stata/src/ado"
adopath ++ "`proj'/ValidMLInference-stata/src/mata"
adopath ++ "`proj'/ValidMLInference-stata/src/help"
quietly mata: mata mlib index

// ---------------------------------------------------------------------------
// Simulation parameters (BCHS 2025 simulation design)
// ---------------------------------------------------------------------------

local n    16000
local m    1000
local p    = 0.05
local kappa_str = 1.0
local fpr  = `kappa_str' / sqrt(`n')

local b0   = 10.0
local b1   = 1.0
local s0   = 0.3
local s1   = 0.5

di _newline(2) as text "{hline 78}"
di as text "Synthetic example (BCHS 2025 simulation design)"
di as text "{hline 78}"
di as text "nsim          = " as result `nsim'
di as text "n (training)  = " as result `n'
di as text "m (test)      = " as result `m'
di as text "p = P(X=1)    = " as result %5.3f `p'
di as text "fpr           = " as result %7.5f `fpr'
di as text "(b0, b1)      = (" as result %5.2f `b0' as text ", " as result %5.2f `b1' as text ")"
di as text "(s0, s1)      = (" as result %5.2f `s0' as text ", " as result %5.2f `s1' as text ")"

tempfile simres
postfile sim int(simid method) double(b0_hat b1_hat se0_hat se1_hat) ///
    using `"`simres'"', replace

set seed 20260522

local started = c(current_time)
di _newline as text "Started simulations at " as result "`started'"

forvalues i = 1/`nsim' {
    quietly {
        clear
        set obs `=`n' + `m''
        generate double u    = runiform()
        generate byte   X    = (u <= `fpr') | (u > 2 * `fpr' & u <= `p' + `fpr')
        generate byte   Xhat = (u >  `fpr' & u <= `p' + `fpr')
        generate double eps  = rnormal()
        generate double Y    = `b0' + `b1' * X + (`s1' * X + `s0' * (1 - X)) * eps
        generate byte   train = _n <= `n'

        // fpr_hat from the test sample, matching the Python convention
        // (mean of Xhat * (1 - X), which equals P(Xhat=1, X=0) on the test split).
        generate double fp = Xhat * (1 - X) if !train
        summarize fp if !train, meanonly
        local fpr_hat = r(mean)
        drop fp
    }

    // -- Method 1: OLS on noisy label
    capture quietly ols Y Xhat if train
    if !_rc {
        matrix b = e(b)
        matrix V = e(V)
        post sim (`i') (1) (b[1, 1]) (b[1, 2]) (sqrt(V[1, 1])) (sqrt(V[2, 2]))
    }

    // -- Method 2: additive correction
    capture quietly ols_bca Y Xhat if train, generated(Xhat) fpr(`fpr_hat') m(`m')
    if !_rc {
        matrix b = e(b)
        matrix V = e(V)
        post sim (`i') (2) (b[1, 1]) (b[1, 2]) (sqrt(V[1, 1])) (sqrt(V[2, 2]))
    }

    // -- Method 3: multiplicative correction
    capture quietly ols_bcm Y Xhat if train, generated(Xhat) fpr(`fpr_hat') m(`m')
    if !_rc {
        matrix b = e(b)
        matrix V = e(V)
        post sim (`i') (3) (b[1, 1]) (b[1, 2]) (sqrt(V[1, 1])) (sqrt(V[2, 2]))
    }

    // -- Method 4: one_step (single-Gaussian MLE)
    capture quietly one_step Y Xhat if train, generated(Xhat)
    if !_rc {
        matrix b = e(b)
        matrix V = e(V)
        post sim (`i') (4) (b[1, 1]) (b[1, 2]) (sqrt(V[1, 1])) (sqrt(V[2, 2]))
    }

    if mod(`i', 100) == 0 | `i' == `nsim' {
        di as text "  Done " as result `i' as text "/" as result `nsim' ///
            as text "  (elapsed start " as result "`started'" as text ", now " ///
            as result c(current_time) as text ")"
    }
}

postclose sim

// ---------------------------------------------------------------------------
// Tabulate coverage of 95% CI for the slope and bias of each estimator
// ---------------------------------------------------------------------------

use `"`simres'"', clear

label define methodlbl 1 "OLS" 2 "ols_bca" 3 "ols_bcm" 4 "one_step"
label values method methodlbl

generate byte covered_b1 = abs(b1_hat - `b1') <= 1.96 * se1_hat

di _newline(2) as text "{hline 78}"
di as text "Coverage of 95% CI for b1 (true b1 = " as result %5.3f `b1' as text ")"
di as text "{hline 78}"
table method, statistic(mean covered_b1) nformat(%6.3f)

di _newline as text "{hline 78}"
di as text "Average estimate and standard error across simulations"
di as text "{hline 78}"
di as text "      method        avg b0        avg b1       avg se0       avg se1   [2.5%, 97.5%] of b1"
di as text "{hline 95}"
forvalues k = 1/4 {
    quietly summarize b0_hat  if method == `k', meanonly
    local ab0 = r(mean)
    quietly summarize b1_hat  if method == `k', meanonly
    local ab1 = r(mean)
    quietly summarize se0_hat if method == `k', meanonly
    local as0 = r(mean)
    quietly summarize se1_hat if method == `k', meanonly
    local as1 = r(mean)
    quietly _pctile b1_hat if method == `k', p(2.5 97.5)
    local q_lo = r(r1)
    local q_hi = r(r2)
    local mname : label methodlbl `k'
    di as text %12s "`mname'" "  " ///
       as result %12.5f `ab0' "  " %12.5f `ab1' "  " ///
       %12.5f `as0' "  " %12.5f `as1' as text "  [" ///
       as result %6.3f `q_lo' as text ", " as result %6.3f `q_hi' as text "]"
}
