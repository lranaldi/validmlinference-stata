// load.do — read the topic-model CSV matrices into Stata matrices.
//
// Usage from a do-file:
//     do "<path>/ValidMLInference-stata/data/topic_model/load.do" "<path>/ValidMLInference-stata/data/topic_model"
//
// The single argument is the directory containing the CSV matrices.
//
// After this script runs, the following Stata matrices exist:
//     W_full   916 x 2  document-topic shares (full sample)
//     W_samp   916 x 2  document-topic shares (10% subsample)
//     B_full   2 x V    topic-word distributions (full sample)
//     B_samp   2 x V    topic-word distributions (10% subsample)
//     S        1 x 2    selection matrix picking the first (leadership) topic
//     LDA      n_lda x 2  C_i for full and subsample (used to build kappa)
//
// The current dataset in memory is preserved.

args here

if `"`here'"' == "" {
    di as error "load.do: pass the data directory as the first argument."
    exit 198
}
// Normalise to a trailing slash.
if substr(`"`here'"', -1, 1) != "/" & substr(`"`here'"', -1, 1) != "\" {
    local here `"`here'/"'
}

quietly {
    preserve

    foreach m in W_full W_samp B_full B_samp LDA {
        capture matrix drop `m'
    }

    foreach pair in "theta_full W_full" "theta_samp W_samp" ///
                    "beta_full B_full"  "beta_samp B_samp"  ///
                    "lda_data LDA" {
        local file  : word 1 of `pair'
        local mat   : word 2 of `pair'
        import delimited `"`here'`file'.csv"', asdouble varnames(nonames) clear
        mkmat _all, matrix(`mat')
    }

    matrix S = (1, 0)

    restore
}
