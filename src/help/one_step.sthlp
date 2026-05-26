{smcl}
{* *! version 0.6.0  26may2026  Lorenzo Ranaldi}{...}
{viewerjumpto "Syntax" "one_step##syntax"}{...}
{viewerjumpto "Description" "one_step##description"}{...}
{viewerjumpto "Options" "one_step##options"}{...}
{viewerjumpto "Stored results" "one_step##results"}{...}
{viewerjumpto "Examples" "one_step##examples"}{...}
{viewerjumpto "References" "one_step##references"}{...}
{viewerjumpto "Author" "one_step##author"}{...}

{title:Title}

{phang}
{bf:one_step} {hline 2} Joint MLE for a regression with a binary AI/ML-generated covariate


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:one_step} {it:depvar} {it:indepvars} {ifin}
[{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opth gen:erated(varname)}}name of the binary AI/ML-generated covariate (default: first variable in {it:indepvars}).{p_end}
{synopt :{opt homosk:edastic}}assume a common error scale across the latent treatment.{p_end}
{synopt :{opt dist(name)}}residual distribution: {cmd:normal} (default), {cmd:laplace}, or {cmd:t}.{p_end}
{synopt :{opt df(#)}}degrees of freedom for {cmd:dist(t)} (required when {cmd:dist(t)}; must be > 0).{p_end}
{synopt :{opt nocons:tant}}omit the intercept term.{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:one_step} fits the joint regression / measurement-error model from
BCHS 2025 by maximum likelihood. One regressor is a binary AI/ML-generated
label observed with classification noise; the latent treatment is integrated
out via a 2 {c x} 2 misclassification table. Conditional on the latent
treatment the residual follows the distribution selected by {cmd:dist()}
({cmd:normal}, {cmd:laplace}, or Student-{cmd:t}; default {cmd:normal}),
with separately-parameterized scales in the two latent classes
(default: heteroskedastic).

{pstd}
This is the right estimator when no external estimate of the false-positive
rate is available: the model identifies the misclassification probabilities
jointly with the regression coefficients, given enough variation in the
outcome.

{pstd}
The optimization runs in Mata via {help mf_optimize:optimize()} using BFGS
on the negative log-likelihood, which is computed entirely in log-space
(logsumexp throughout) to stay numerically well-behaved even when the
optimizer probes parameters that would otherwise drive
{it:exp(log sigma)} past the IEEE double overflow threshold. Variance is the
top-d block of {it:pinv(H)}, where H is the numerical Hessian at the
optimum.


{marker options}{...}
{title:Options}

{phang}
{opth generated(varname)} names the binary AI/ML-generated covariate.
Defaults to the first variable in {it:indepvars}. The selected variable must
take values 0 and 1 only.

{phang}
{opt homoskedastic} constrains {it:sigma_0 = sigma_1}. By default
{cmd:one_step} estimates the two latent-class residual scales separately.

{phang}
{opt dist(name)} selects the residual distribution. Allowed values:
{cmd:normal} (default), {cmd:laplace}, or {cmd:t}. The scale parameter
{it:sigma} is reinterpreted per distribution (Laplace {it:b}, Student-{it:t}
scale); the theta layout is unchanged. Mirrors Python's
{cmd:one_step(..., distribution=...)} as a closed menu instead of a
JAX-traceable callable.

{phang}
{opt df(#)} sets the degrees of freedom for {cmd:dist(t)}. Required and
must be {it:> 0} when {cmd:dist(t)} is specified; an error is raised if
supplied with any other distribution. The {it:df} is fixed (not estimated),
matching the Python pattern.

{phang}
{opt noconstant} omits the intercept term.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:one_step} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(ll)}}log-likelihood at the optimum{p_end}
{synopt:{cmd:e(converged)}}1 if the optimizer converged, 0 otherwise{p_end}
{synopt:{cmd:e(iterations)}}number of BFGS iterations{p_end}
{synopt:{cmd:e(homoskedastic)}}1 if {opt homoskedastic} was specified{p_end}
{synopt:{cmd:e(distcode)}}internal distribution code (1=normal, 2=laplace, 3=t){p_end}
{synopt:{cmd:e(df)}}degrees of freedom (posted only when {cmd:dist(t)}){p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:one_step}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(dist)}}residual distribution name (normal/laplace/t){p_end}
{synopt:{cmd:e(generated_var)}}name of the binary generated covariate{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(vce)}}{cmd:OIM (numerical Hessian)}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:OIM}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}regression coefficients (intercept first, then alphabetical){p_end}
{synopt:{cmd:e(V)}}top-d block of {it:pinv(numerical Hessian)} at the optimum{p_end}


{marker examples}{...}
{title:Examples}

{phang}Synthetic example mirroring the M6 parity fixture (true fpr=fnr=0.10, heteroskedastic Gaussian errors):{p_end}

{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set seed 20260518}{p_end}
{phang}{cmd:. set obs 2000}{p_end}
{phang}{cmd:. generate double xstar = runiform() < 0.4}{p_end}
{phang}{cmd:. generate double flip  = runiform() < 0.10}{p_end}
{phang}{cmd:. generate double x1    = cond(flip, 1 - xstar, xstar)}{p_end}
{phang}{cmd:. generate double x2    = rnormal()}{p_end}
{phang}{cmd:. generate double eps   = rnormal() * cond(xstar, 0.8, 0.5)}{p_end}
{phang}{cmd:. generate double y     = 0.2 + 0.7*xstar - 0.3*x2 + eps}{p_end}
{phang}{cmd:. one_step y x1 x2, generated(x1)}{p_end}

{phang}Same model with Laplace residuals:{p_end}
{phang}{cmd:. one_step y x1 x2, generated(x1) dist(laplace)}{p_end}

{phang}Student-{it:t} residuals with five degrees of freedom:{p_end}
{phang}{cmd:. one_step y x1 x2, generated(x1) dist(t) df(5)}{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
{it:Python parity.}  Mata's BFGS and Python's jaxopt LBFGS converge to
nearly identical optima on the synthetic test cases. Gaussian fixture:
max |db| {c 0210} 1e-6, max |dV| {c 0210} 1e-7. Laplace fixture:
max |db| {c 0210} 3e-4, max |dV| {c 0210} 2e-3. Student-{it:t}{c 0210}5 fixture:
max |db| {c 0210} 3e-6, max |dV| {c 0210} 7e-9. All inside the declared 1e-3
/ 1e-2 MLE tolerances.

{pstd}
{it:Distribution menu vs callback.}  Python's
{cmd:one_step(..., distribution=...)} accepts an arbitrary JAX-traceable
PDF callable; Stata exposes the three most common cases as a closed menu
(no clean Stata analogue for a JAX-traceable callable). Starting values use
the Gaussian pdf for the latent-X* imputation step regardless of
{cmd:dist()}, which matches the Python implementation's
{cmd:get_starting_values_unlabeled_jax_with_treatment_idx} (hard-coded
Gaussian {it:pdf_func}).

{pstd}
{it:Numerical caveats.}  The Laplace log-density is non-smooth at
{it:y = loc}; in continuous samples no observation lies exactly on the
kink, but the numerical Hessian can be slightly noisy near the optimum.
The Student-{it:t} log-density is smooth and converges as cleanly as the
Gaussian case.


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
Help: {help one_step_gmm}, {help ols_bca}, {help ols_bcm}
{p_end}
