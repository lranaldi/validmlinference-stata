{smcl}
{* *! version 0.5.0  22may2026  Lorenzo Ranaldi}{...}
{viewerjumpto "Syntax" "ols_bca_topic##syntax"}{...}
{viewerjumpto "Description" "ols_bca_topic##description"}{...}
{viewerjumpto "Options" "ols_bca_topic##options"}{...}
{viewerjumpto "Stored results" "ols_bca_topic##results"}{...}
{viewerjumpto "Examples" "ols_bca_topic##examples"}{...}
{viewerjumpto "References" "ols_bca_topic##references"}{...}
{viewerjumpto "Author" "ols_bca_topic##author"}{...}

{title:Title}

{phang}
{bf:ols_bca_topic} {hline 2} Additive bias correction for topic-model-generated regressors


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:ols_bca_topic} {it:depvar} [{it:Q_varlist}] {ifin}{cmd:,}
{opt wmatrix(matname)} {opt smatrix(matname)} {opt bmatrix(matname)}
{opt k(real)}
[{it:options}]

{synoptset 24 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opt wmatrix(matname)}}n {c x} c document-by-component matrix (row-stochastic).{p_end}
{synopt :{opt smatrix(matname)}}r {c x} c topic-by-component matrix (row-stochastic).{p_end}
{synopt :{opt bmatrix(matname)}}c {c x} v component-by-vocabulary matrix (row-stochastic).{p_end}
{synopt :{opth k(real)}}bias-correction scaling parameter (typically of order 1).{p_end}

{syntab:Optional}
{synopt :{opt nocons:tant}}omit the intercept term.{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ols_bca_topic} fits a regression in which the regressors are document-
level topic shares estimated from a topic model (plus optional additional
covariates {it:Q_varlist}), and applies an {it:additive} bias correction to
account for the topic-model estimation error in those shares.

{pstd}
The design matrix is {it:Xhat = (Theta, Q)} with topic shares
{it:Theta = W S'} (n {c x} r). The correction multiplies the uncorrected OLS
coefficient vector by {it:(I + Gamma)}, where Gamma is built from the topic-
model bias term

{p 8 8 2}{it:Omega = S (B B')^(-1) B (B'  *  (B' mean(W))) (B B')^(-1) S' - Theta'Theta / n,}{p_end}

{pstd}
yielding {it:Gamma = k * sqrt(n) * (Xhat'Xhat)^(-1) * A}, where {it:A} has
{it:Omega} in its top-left r {c x} r block. The variance posted in {cmd:e(V)}
is the uncorrected HC0 variance — matching the upstream Python implementation,
which does not inflate SEs for the topic-model correction.

{pstd}
The topic-share columns are named {cmd:topic_1}, ..., {cmd:topic_r} in
{cmd:e(b)} (intercept first, then non-intercept names alphabetically). See
{help ols_bcm_topic:ols_bcm_topic} for the multiplicative analogue.


{marker options}{...}
{title:Options}

{phang}
{opt wmatrix(matname)} is required. n {c x} c document-by-component matrix
from the topic model. Each row should sum to 1.

{phang}
{opt smatrix(matname)} is required. r {c x} c topic-by-component matrix.

{phang}
{opt bmatrix(matname)} is required. c {c x} v component-by-vocabulary matrix.

{phang}
{opth k(real)} is required. Bias-correction scaling parameter. In BCHS 2025
this is typically chosen of order 1.

{phang}
{opt noconstant} omits the intercept term.

{pstd}
{it:Why not just `b()`, `s()`, `w()` for the matrices?}
{cmd:syntax} would clash with internal {cmd:tempname b}. The fully spelled-out
forms {cmd:wmatrix()}, {cmd:smatrix()}, {cmd:bmatrix()} keep the parse
unambiguous.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ols_bca_topic} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}bias-correction scaling{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ols_bca_topic}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(W_matrix)}}name of the W matrix passed in{p_end}
{synopt:{cmd:e(S_matrix)}}name of the S matrix passed in{p_end}
{synopt:{cmd:e(B_matrix)}}name of the B matrix passed in{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(vce)}}{cmd:robust}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust (HC0)}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}bias-corrected coefficient vector{p_end}
{synopt:{cmd:e(V)}}HC0 variance from the uncorrected OLS fit{p_end}


{marker examples}{...}
{title:Examples}

{phang}Small synthetic example (n=300, r=3 topics, c=5 components, v=8 vocab terms):{p_end}

{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set seed 20260518}{p_end}
{phang}{cmd:. set obs 300}{p_end}
{phang}{cmd:. * Build W, S, B in Mata to keep the example self-contained}{p_end}
{phang}{cmd:. mata:}{p_end}
{phang}{cmd:.   W = runiform(300, 5) ; W = W :/ rowsum(W) ; st_matrix("W", W)}{p_end}
{phang}{cmd:.   S = runiform(3, 5)   ; S = S :/ rowsum(S) ; st_matrix("S", S)}{p_end}
{phang}{cmd:.   B = runiform(5, 8)   ; B = B :/ rowsum(B) ; st_matrix("B", B)}{p_end}
{phang}{cmd:. end}{p_end}
{phang}{cmd:. generate double q1 = rnormal()}{p_end}
{phang}{cmd:. generate double q2 = rnormal()}{p_end}
{phang}{cmd:. generate double y  = rnormal()  // fill in your DGP here}{p_end}
{phang}{cmd:. ols_bca_topic y q1 q2, wmatrix(W) smatrix(S) bmatrix(B) k(0.5)}{p_end}


{marker references}{...}
{title:References}

{phang}
Battaglia, L., Christensen, T., Hansen, S., and Sacher, S. 2025.
Inference for Regression with Variables Generated by AI or Machine Learning.
{browse "https://arxiv.org/abs/2402.15585":arXiv:2402.15585}.

{phang}
Christensen, T., and Hansen, S. 2025.
Performing Valid Inference with AI/ML-Generated Covariates: A Guide for
Empirical Practice.


{marker author}{...}
{title:Author}

{pstd}Lorenzo Ranaldi.{p_end}
{pstd}Stata port of the Python package {bf:ValidMLInference} by
Kurczynski, Christensen, Hansen, Battaglia, and Sacher.{p_end}


{title:Also see}

{p 4 14 2}
Help: {help ols}, {help ols_bca}, {help ols_bcm},
{help ols_bcm_topic}, {help one_step}, {help one_step_gmm}
{p_end}
