*! test_pkg_install.do -- sanity-check `net install` from a local path.
*!
*! Simulates a clean install: redirects PLUS (and PERSONAL) to a tempdir
*! so the dev `adopath` is not consulted, runs `net install` from the
*! local `ValidMLInference-stata/` source, then verifies that:
*!   (a) every public command is discoverable on the adopath,
*!   (b) every help file renders,
*!   (c) a minimal one_step call works end-to-end (the .mlib loads).
*!
*! Run from the project root:
*!     do "ValidMLInference-stata/tests/test_pkg_install.do"

clear all
version 17
set more off

* ---------------------------------------------------------------------------
* Locate the package source. `from()` for `net install` must point at the
* directory containing stata.toc + validmlinference.pkg. Try a few likely
* locations so the do-file works whether launched from the project root,
* from ValidMLInference-stata/, or from tests/.
* ---------------------------------------------------------------------------
local pkg_src ""
foreach cand in                            ///
        "`c(pwd)'/ValidMLInference-stata"  ///
        "`c(pwd)'"                         ///
        "`c(pwd)'/.."                      ///
        "C:/Users/loren/Dropbox/Replication_AIValidInference/ValidMLInference-stata" {
    capture confirm file "`cand'/validmlinference.pkg"
    if !_rc {
        local pkg_src "`cand'"
        continue, break
    }
}
if "`pkg_src'" == "" {
    di as error "test_pkg_install.do: cannot find validmlinference.pkg"
    di as error "  cwd was: `c(pwd)'"
    exit 601
}

* ---------------------------------------------------------------------------
* Redirect PLUS / PERSONAL to a fresh temp directory so this install does not
* pollute the user's existing ado tree (and so adopath does not silently
* resolve commands from the dev src/ado).
* ---------------------------------------------------------------------------
tempfile junk
local tmp = subinstr("`junk'", "\", "/", .)
local tmp = subinstr("`tmp'", ".tmp", "_pkgtest", .)
capture mkdir "`tmp'"
capture mkdir "`tmp'/plus"
capture mkdir "`tmp'/personal"

local saved_plus     "`c(sysdir_plus)'"
local saved_personal "`c(sysdir_personal)'"

sysdir set PLUS     "`tmp'/plus/"
sysdir set PERSONAL "`tmp'/personal/"

* Remove the dev `src/ado` / `src/mata` / `src/help` from adopath if they're
* there (they would mask the fresh install otherwise).
foreach sub in src/ado src/mata src/help {
    capture adopath - "`pkg_src'/`sub'"
}

di as txt _newline "=== ValidMLInference net-install sanity check ==="
di as txt "package source: `pkg_src'"
di as txt "temp PLUS:      `c(sysdir_plus)'"
di as txt "temp PERSONAL:  `c(sysdir_personal)'"

* ---------------------------------------------------------------------------
* Run `net install`. The `replace` is harmless on a fresh PLUS but lets the
* script be rerun without manual cleanup.
* ---------------------------------------------------------------------------
di _newline as txt ">>> net install validmlinference"
capture noisily net install validmlinference, from("`pkg_src'") replace
if _rc {
    di as error "net install failed (rc=" _rc ")"
    sysdir set PLUS     "`saved_plus'"
    sysdir set PERSONAL "`saved_personal'"
    exit _rc
}

* Force Stata to refresh its Mata library index in case lvmli.mlib was
* added after the session started.
quietly mata: mata mlib index

* ---------------------------------------------------------------------------
* (a) every command discoverable
* ---------------------------------------------------------------------------
di _newline as txt ">>> which <each command>"
local n_cmd_ok = 0
foreach cmd in ols ols_bca ols_bcm ols_bca_topic ols_bcm_topic one_step one_step_gmm {
    capture which `cmd'
    if _rc {
        di as error "  FAIL: `cmd' not on adopath after install"
    }
    else {
        di as result "  PASS: `cmd'"
        local ++n_cmd_ok
    }
}

* ---------------------------------------------------------------------------
* (b) every help file renders (no rc)
* ---------------------------------------------------------------------------
di _newline as txt ">>> help <each command>  (capture only -- not displayed)"
local n_help_ok = 0
foreach cmd in ols ols_bca ols_bcm ols_bca_topic ols_bcm_topic one_step one_step_gmm {
    capture quietly help `cmd'
    if _rc {
        di as error "  FAIL: help `cmd' (rc=" _rc ")"
    }
    else {
        di as result "  PASS: help `cmd'"
        local ++n_help_ok
    }
}

* ---------------------------------------------------------------------------
* (c) minimal end-to-end run -- exercises the .mlib too
* ---------------------------------------------------------------------------
di _newline as txt ">>> end-to-end smoke test: one_step on a tiny synthetic sample"

set seed 20260526
quietly {
    clear
    set obs 1000
    generate double xstar = runiform() < 0.4
    generate double flip  = runiform() < 0.10
    generate double x1    = cond(flip, 1 - xstar, xstar)
    generate double x2    = rnormal()
    generate double eps   = rnormal() * cond(xstar, 0.8, 0.5)
    generate double y     = 0.2 + 0.7*xstar - 0.3*x2 + eps
}

local smoke_ok = 0
capture noisily one_step y x1 x2, generated(x1)
if !_rc {
    di as result "  PASS: one_step ran (b[x1] = " %6.3f _b[x1] ")"
    local smoke_ok = 1
}
else {
    di as error  "  FAIL: one_step (rc=" _rc ")"
}

* Verify the dist() option parser is wired (do not run the MLE -- model fit on
* arbitrary data is not part of a packaging sanity check; parity-test fixtures
* under tests/fixtures/ exercise the laplace/t paths properly).
capture noisily one_step y x1 x2, generated(x1) dist(bogus)
if _rc == 198 {
    di as result "  PASS: dist() parser rejects unknown values (rc=198)"
}
else {
    di as error  "  FAIL: dist(bogus) should error with rc=198 (got " _rc ")"
    local smoke_ok = 0
}

* ---------------------------------------------------------------------------
* Restore session sysdir so we don't leak the tempdir into the user's state.
* ---------------------------------------------------------------------------
sysdir set PLUS     "`saved_plus'"
sysdir set PERSONAL "`saved_personal'"

di _newline _newline as txt "=== Summary ==="
di as txt    "  commands discovered: " as result %2.0f `n_cmd_ok'  as txt " / 7"
di as txt    "  help pages rendered: " as result %2.0f `n_help_ok' as txt " / 7"
di as txt    "  smoke test:          " as result cond(`smoke_ok', "PASS", "FAIL")

if `n_cmd_ok' == 7 & `n_help_ok' == 7 & `smoke_ok' {
    di _newline as result "ALL CHECKS PASS -- package ready for distribution."
}
else {
    di _newline as error "Sanity check FAILED -- see messages above."
    exit 9
}
