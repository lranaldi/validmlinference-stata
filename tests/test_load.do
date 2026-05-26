*! test_load.do -- Milestone 1 smoke test: confirm the package loads
*!
*! Builds the Mata library, sets the adopath, calls vmli_version(), and
*! invokes every .ado stub to confirm each one is found and exits 198
*! ("not yet implemented") as expected.
*!
*! Usage:  set the project root with PROJECT_ROOT or run from the project
*!         root containing the ValidMLInference-stata/ folder.

version 17
set more off
clear all

* ---------------------------------------------------------------------------
* 1. Locate the project root
*
*    Marker file: ValidMLInference-stata/src/ado/ols.ado .
*    Try cwd first; if not found, try walking up to two levels (covers the
*    common case of being launched from .../ValidMLInference-stata/tests/).
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
    di as error "test_load.do: cannot find `marker' from cwd = `c(pwd)'"
    di as error "Run from the project root or from anywhere inside it."
    exit 601
}

local proj "`c(pwd)'"
di as txt "Project root: `proj'"

* ---------------------------------------------------------------------------
* 2. Build the Mata library
* ---------------------------------------------------------------------------
di as txt _newline "--- Building Mata library ---"
qui cd "`pkg'/src/mata"
do _build.do
qui cd "`proj'"

* ---------------------------------------------------------------------------
* 3. Configure ado/mlib search paths
* ---------------------------------------------------------------------------
adopath ++ "`pkg'/src/ado"
adopath ++ "`pkg'/src/mata"
mata: mata mlib index

* ---------------------------------------------------------------------------
* 4. Smoke test 1 -- Mata library loaded
* ---------------------------------------------------------------------------
di as txt _newline "--- Mata library check ---"
mata: vmli_version()

* ---------------------------------------------------------------------------
* 5. Smoke test 2 -- every planned command is discoverable on the adopath
* ---------------------------------------------------------------------------
di as txt _newline "--- Command discoverability ---"
local all_cmds ols ols_bca ols_bcm ols_bca_topic ols_bcm_topic one_step one_step_gmm
local fail 0
foreach cmd of local all_cmds {
    cap qui which `cmd'
    if (_rc == 0) {
        di as result "  OK   `cmd' is on the adopath"
    }
    else {
        di as error  "  FAIL `cmd' not found (rc=`=_rc')"
        local ++fail
    }
}

* ---------------------------------------------------------------------------
* 6. Summary
* ---------------------------------------------------------------------------
di as txt _newline "--- Summary ---"
if (`fail' == 0) {
    di as result "test_load: PASS (Mata library loaded, all commands discoverable, stubs OK)."
}
else {
    di as error  "test_load: FAIL (`fail' issue(s))."
    exit 9
}
