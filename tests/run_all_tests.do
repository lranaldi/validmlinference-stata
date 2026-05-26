*! run_all_tests.do -- Master parity harness (Milestone 5)
*!
*! Discovers every subfolder in tests/fixtures/, reads its meta.txt, runs the
*! matching Stata command, and compares e(b)/e(V) to the Python reference
*! stored in ref_coef.csv / ref_vcov.csv. Tolerances come from the fixture
*! meta (tol_coef / tol_vcov).
*!
*! Comparison is by name: both Python and Stata standardize coefficient order
*! to "_cons first, then alphabetical", so the row orderings should match,
*! but we still align by name to surface any mismatch loudly.

version 17
set more off
clear all

* ---------------------------------------------------------------------------
* Locate project root and build library
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
    di as error "run_all_tests.do: cannot locate `marker' from cwd"
    exit 601
}
local proj "`c(pwd)'"

qui cd "`pkg'/src/mata"
qui do _build.do
qui cd "`proj'"

qui adopath ++ "`pkg'/src/ado"
qui adopath ++ "`pkg'/src/mata"
qui mata: mata mlib index

foreach prog in ols ols_bca ols_bcm ols_bca_topic ols_bcm_topic one_step one_step_gmm {
    cap program drop `prog'
}

* ---------------------------------------------------------------------------
* Mata helper -- comparison kernel. Defined up front so it's in scope when
* the loop below calls it. Reads e(b)/e(V) and aligns to Mata globals
* g_ref_names / g_ref_coef / g_vcov_rownames / g_ref_vcov by name.
* ---------------------------------------------------------------------------
mata:
void _vmli_build_vcov_from_long()
{
    string colvector row_s, col_s, uniq_names
    real colvector vcov_v, found
    real scalar n_rows, d, i, idx_r, idx_c

    external string colvector g_vcov_rownames
    external real matrix      g_ref_vcov

    row_s  = st_sdata(., "row")
    col_s  = st_sdata(., "col")
    vcov_v = st_data(., "vcov")
    n_rows = length(row_s)

    uniq_names = J(0, 1, "")
    for (i = 1; i <= n_rows; i++) {
        if (length(selectindex(uniq_names :== row_s[i])) == 0) {
            uniq_names = uniq_names \ row_s[i]
        }
    }
    d = length(uniq_names)

    g_vcov_rownames = uniq_names
    g_ref_vcov = J(d, d, .)

    for (i = 1; i <= n_rows; i++) {
        found = selectindex(uniq_names :== row_s[i])
        idx_r = found[1]
        found = selectindex(uniq_names :== col_s[i])
        idx_c = found[1]
        g_ref_vcov[idx_r, idx_c] = vcov_v[i]
    }
}

void _vmli_test_compare()
{
    real rowvector b_sta, b_ref_aligned
    real matrix V_sta, V_ref_aligned
    string matrix stripe
    string colvector sta_names
    real scalar i, j, ii, jj, d
    real colvector found

    external string colvector g_ref_names
    external real colvector   g_ref_coef
    external string colvector g_vcov_rownames
    external real matrix      g_ref_vcov

    b_sta = st_matrix("e(b)")
    V_sta = st_matrix("e(V)")
    stripe = st_matrixcolstripe("e(b)")
    sta_names = stripe[, 2]
    d = length(sta_names)

    b_ref_aligned = J(1, d, .)
    for (i = 1; i <= d; i++) {
        found = selectindex(g_ref_names :== sta_names[i])
        if (length(found) != 1) {
            errprintf("name '%s' missing or duplicated in ref_coef\n",
                      sta_names[i])
            exit(111)
        }
        b_ref_aligned[i] = g_ref_coef[found[1]]
    }

    V_ref_aligned = J(d, d, .)
    for (i = 1; i <= d; i++) {
        found = selectindex(g_vcov_rownames :== sta_names[i])
        if (length(found) != 1) {
            errprintf("name '%s' missing in ref_vcov rows\n", sta_names[i])
            exit(111)
        }
        ii = found[1]
        for (j = 1; j <= d; j++) {
            found = selectindex(g_vcov_rownames :== sta_names[j])
            if (length(found) != 1) {
                errprintf("name '%s' missing in ref_vcov cols\n",
                          sta_names[j])
                exit(111)
            }
            jj = found[1]
            V_ref_aligned[i, j] = g_ref_vcov[ii, jj]
        }
    }

    st_numscalar("max_db", max(abs(b_sta - b_ref_aligned)))
    st_numscalar("max_dV", max(abs(vec(V_sta - V_ref_aligned))))
}
end

* ---------------------------------------------------------------------------
* Discover cases
* ---------------------------------------------------------------------------
local fixroot "`pkg'/tests/fixtures"
local cases : dir "`fixroot'" dirs "*"

if `"`cases'"' == "" {
    di as error "no fixtures found under `fixroot' -- run refresh_fixtures.py first"
    exit 601
}

di as txt _newline "=== ValidMLInference Stata parity harness ==="
di as txt "Project root: `proj'"
di as txt "Fixtures:     `fixroot'"

local n_total = 0
local n_pass  = 0
local n_fail  = 0
local failed  ""

foreach case of local cases {
    local case_dir "`fixroot'/`case'"
    local n_total = `n_total' + 1

    di _newline as txt "------ `case' ------"

    * --- read meta.txt -----------------------------------------------------
    cap confirm file "`case_dir'/meta.txt"
    if _rc {
        di as error "  meta.txt missing -> FAIL"
        local n_fail = `n_fail' + 1
        local failed `"`failed' `case'"'
        continue
    }

    foreach key in case command depvar covars options tol_coef tol_vcov n note {
        local m_`key' ""
    }

    qui file open _mf using "`case_dir'/meta.txt", read text
    file read _mf line
    while r(eof) == 0 {
        local eq = strpos(`"`macval(line)'"', "=")
        if `eq' > 0 {
            local key = trim(substr(`"`macval(line)'"', 1, `eq' - 1))
            local val = trim(substr(`"`macval(line)'"', `eq' + 1, .))
            local m_`key' `"`macval(val)'"'
        }
        file read _mf line
    }
    file close _mf

    local cmd      "`m_command'"
    local depvar   "`m_depvar'"
    local covars   "`m_covars'"
    local opts     `"`m_options'"'
    local tol_b    = real("`m_tol_coef'")
    local tol_v    = real("`m_tol_vcov'")
    if missing(`tol_b') local tol_b = 1e-8
    if missing(`tol_v') local tol_v = 1e-7

    * --- read reference coef/vcov into Mata globals ------------------------
    cap confirm file "`case_dir'/ref_coef.csv"
    if _rc {
        di as error "  ref_coef.csv missing -> FAIL"
        local n_fail = `n_fail' + 1
        local failed `"`failed' `case'"'
        continue
    }
    qui import delimited "`case_dir'/ref_coef.csv", ///
        delimiter(",") varnames(1) stringcols(1) asdouble clear
    mata: g_ref_names = st_sdata(., "name")
    mata: g_ref_coef  = st_data(., "coef")

    qui import delimited "`case_dir'/ref_vcov.csv", ///
        delimiter(",") varnames(1) stringcols(1 2) asdouble clear
    mata: _vmli_build_vcov_from_long()

    * --- load W/S/B matrices (topic cases) ---------------------------------
    foreach M in W S B {
        cap confirm file "`case_dir'/`M'.csv"
        if !_rc {
            qui import delimited "`case_dir'/`M'.csv", ///
                delimiter(",") varnames(nonames) asdouble clear
            qui mkmat _all, matrix(`M')
        }
    }

    * --- load actual input data -------------------------------------------
    cap confirm file "`case_dir'/input.csv"
    if _rc {
        di as error "  input.csv missing -> FAIL"
        local n_fail = `n_fail' + 1
        local failed `"`failed' `case'"'
        continue
    }
    qui import delimited "`case_dir'/input.csv", ///
        delimiter(",") varnames(1) asdouble clear

    * --- assemble and run the Stata command --------------------------------
    local stcmd `"`cmd' `depvar' `covars'"'
    if `"`opts'"' != "" {
        local stcmd `"`stcmd', `opts'"'
    }
    cap noi qui `stcmd'
    if _rc {
        di as error "  command failed (rc=" _rc "): `stcmd' -> FAIL"
        local n_fail = `n_fail' + 1
        local failed `"`failed' `case'"'
        continue
    }

    * --- compare e(b)/e(V) to refs, aligned by name -----------------------
    mata: _vmli_test_compare()

    local pass_b = (max_db < `tol_b')
    local pass_v = (max_dV < `tol_v')

    if (`pass_b' & `pass_v') {
        di as result "  PASS  (max|db|=" %10.2e max_db ", max|dV|=" %10.2e max_dV ")"
        local n_pass = `n_pass' + 1
    }
    else {
        di as error "  FAIL  (max|db|=" %10.2e max_db " tol=" %8.1e `tol_b' ///
                   ", max|dV|=" %10.2e max_dV " tol=" %8.1e `tol_v' ")"
        local n_fail = `n_fail' + 1
        local failed `"`failed' `case'"'
    }
}

di _newline _newline as txt "=== Summary ==="
di as txt    "  total: `n_total'"
di as result "  pass : `n_pass'"
if `n_fail' > 0 {
    di as error "  fail : `n_fail' [`failed' ]"
    exit 9
}
else {
    di as result _newline "ALL CASES PASS within declared tolerances."
}
