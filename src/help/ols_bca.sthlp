{smcl}
{* *! version 0.5.0  22may2026  Lorenzo Ranaldi}{...}
{viewerjumpto "Syntax" "ols_bca##syntax"}{...}
{viewerjumpto "Description" "ols_bca##description"}{...}
{viewerjumpto "Options" "ols_bca##options"}{...}
{viewerjumpto "Stored results" "ols_bca##results"}{...}
{viewerjumpto "Examples" "ols_bca##examples"}{...}
{viewerjumpto "References" "ols_bca##references"}{...}
{viewerjumpto "Author" "ols_bca##author"}{...}

{title:Title}

{phang}
{bf:ols_bca} {hline 2} Additive bias correction for a binary AI/ML-generated regressor


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:ols_bca} {it:depvar} {it:indepvars} {ifin}{cmd:,}
{opt fpr(real)} {opt m(real)}
[{it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt :{opth fpr(real)}}estimated false-positive rate of the classifier (in [0, 1]).{p_end}
{synopt :{opth m(real)}}size of the external sample used to estimate {opt fpr}.{p_end}

{syntab:Optional}
{synopt :{opth gen:erated(varname)}}name of the AI/ML-generated covariate (default: first variable in {it:indepvars}).{p_end}
{synopt :{opt nocons:tant}}omit the intercept term.{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ols_bca} estimates a linear regression in which one covariate is a
binary AI/ML-generated label observed with classification noise, and applies
an {it:additive} bias correction to recover an unbiased coefficient on the
underlying latent variable.

{pstd}
The correction needs an external estimate of the classifier's false-positive
rate {opt fpr} (the probability that the label is 1 when the latent is 0),
together with the sample size {opt m} used to obtain that estimate.
{opt m} enters the standard-error formula as an
{it:fpr (1 - fpr) / m} inflation term that propagates uncertainty in the
{opt fpr} estimate into inference on the regression coefficient. When
{opt fpr} is treated as known exactly, set {opt m} to a very large number
(e.g. {cmd:m(1e10)}).

{pstd}
Let {it:b0} and {it:V0} be the uncorrected OLS coefficients and HC0 variance,
let {it:Gamma = n * (X'X)^(-1) * A} where {it:A} has 1 in the (target, target)
entry and 0 elsewhere, and let {it:M = I + fpr * Gamma}. Then

{p 8 8 2}{it:b_corr  = b0 + fpr * Gamma * b0}{p_end}
{p 8 8 2}{it:V_corr  = M V0 M' + (fpr (1 - fpr) / m) * Gamma (V0 + b_corr b_corr') Gamma'.}{p_end}

{pstd}
{cmd:ols_bca} mirrors Python's {bf:ValidMLInference.ols_bca}. See
{help ols_bcm:ols_bcm} for the multiplicative analogue, which is the
recommended estimator for binary imputed labels per BCHS 2025.


{marker options}{...}
{title:Options}

{phang}
{opth fpr(real)} is required. Estimated false-positive rate of the
classifier. Must lie in [0, 1].

{phang}
{opth m(real)} is required. Size of the external sample used to estimate
{opt fpr}. Pass a large number (e.g. 1e10) when {opt fpr} is treated as
known exactly.

{phang}
{opth generated(varname)} names the AI/ML-generated covariate to which the
correction is applied. Defaults to the first variable in {it:indepvars}.
(The Python keyword {cmd:generated_var} is renamed to {cmd:generated} here
because Stata option names cannot contain underscores.)

{phang}
{opt noconstant} omits the intercept term.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ols_bca} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(fpr)}}false-positive rate used in the correction{p_end}
{synopt:{cmd:e(m)}}external sample size{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ols_bca}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(generated_var)}}name of the generated covariate{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(vce)}}{cmd:robust}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:Robust (HC0)}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}bias-corrected coefficient vector{p_end}
{synopt:{cmd:e(V)}}corrected variance, including the {it:fpr (1 - fpr) / m} inflation{p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}


{marker examples}{...}
{title:Examples}

{phang}Synthetic example with a 10% false-positive rate:{p_end}

{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set seed 20260518}{p_end}
{phang}{cmd:. set obs 2000}{p_end}
{phang}{cmd:. generate double xstar = runiform() < 0.4}{p_end}
{phang}{cmd:. generate double flip  = runiform() < 0.10}{p_end}
{phang}{cmd:. generate double x1    = cond(flip, 1 - xstar, xstar)}{p_end}
{phang}{cmd:. generate double x2    = rnormal()}{p_end}
{phang}{cmd:. generate double y     = 0.2 + 0.7*xstar - 0.3*x2 + rnormal()*0.5}{p_end}
{phang}{cmd:. ols_bca y x1 x2, fpr(0.10) m(2000) generated(x1)}{p_end}


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
Help: {help ols}, {help ols_bcm}, {help ols_bca_topic},
{help ols_bcm_topic}, {help one_step}, {help one_step_gmm}
{p_end}
