{smcl}
{* *! version 0.5.0  22may2026  Lorenzo Ranaldi}{...}
{viewerjumpto "Syntax" "one_step_gmm##syntax"}{...}
{viewerjumpto "Description" "one_step_gmm##description"}{...}
{viewerjumpto "Options" "one_step_gmm##options"}{...}
{viewerjumpto "Stored results" "one_step_gmm##results"}{...}
{viewerjumpto "Examples" "one_step_gmm##examples"}{...}
{viewerjumpto "References" "one_step_gmm##references"}{...}
{viewerjumpto "Author" "one_step_gmm##author"}{...}

{title:Title}

{phang}
{bf:one_step_gmm} {hline 2} Joint MLE with a Gaussian-mixture error term


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmd:one_step_gmm} {it:depvar} {it:indepvars} {ifin}
[{cmd:,} {it:options}]

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opth gen:erated(varname)}}name of the binary AI/ML-generated covariate (default: first variable in {it:indepvars}).{p_end}
{synopt :{opth k(integer)}}number of mixture components (default 2; must be >= 2).{p_end}
{synopt :{opth ng:uess(integer)}}number of multistart attempts (default 10; must be >= 1).{p_end}
{synopt :{opth max:iter(integer)}}maximum BFGS iterations per attempt (default 100).{p_end}
{synopt :{opth seed(integer)}}seed for the multistart perturbations (default 0).{p_end}
{synopt :{opt homosk:edastic}}{it:not currently supported} — see {it:Remarks}.{p_end}
{synopt :{opt nocons:tant}}omit the intercept term.{p_end}
{synoptline}


{marker description}{...}
{title:Description}

{pstd}
{cmd:one_step_gmm} is the Gaussian-mixture analogue of {help one_step:one_step}.
Outcome residuals are modeled as a {opt k}-component Gaussian mixture in each
latent-treatment class — i.e., the model can accommodate skewness, fat tails,
or bimodality in the residuals while still identifying the misclassification
probabilities.

{pstd}
The optimizer is run from {opt nguess} starting points: the first uses the
OLS-based starts directly, subsequent attempts perturb the starts by Gaussian
noise (scaled 0.5 on coefficient slots, 0.3 on the misclassification logits,
1.0 elsewhere — matching Python's
{bf:_one_step_gaussian_mixture_core}). The best attempt by negative log-
likelihood wins.

{pstd}
The likelihood is computed in log-space throughout (logsumexp for both the
inner component mixtures and the outer 2-class mixture, plus Python's +1e-12
stabilizer in log-domain).


{marker options}{...}
{title:Options}

{phang}
{opth generated(varname)} names the binary generated covariate. Defaults to
the first variable in {it:indepvars}. Must take values 0 and 1.

{phang}
{opth k(integer)} sets the number of mixture components per latent class.
Default is 2; must be >= 2.

{phang}
{opth nguess(integer)} sets the number of multistart attempts. Default is 10.
The first attempt always uses the OLS-based starts unchanged; subsequent
attempts add noise per Python's schedule.

{phang}
{opth maxiter(integer)} sets the maximum number of BFGS iterations per
multistart attempt. Default is 100.

{phang}
{opth seed(integer)} seeds Stata's RNG for the multistart perturbations.
Default 0. (Stata's RNG is unrelated to JAX's PRNGKey, so per-attempt parity
with the Python implementation is not expected for {it:i > 0}; the {it:i = 0}
attempt is deterministic and identical across implementations.)

{phang}
{opt homoskedastic} is currently rejected at runtime. The upstream Python
{cmd:get_starting_values_unlabeled_gaussian_mixture(homosked=True)} returns
a theta of length {it:d + 1 + 3k} but {cmd:unpack_theta(homosked=True)} reads
{it:d - 1 + 5k}; supporting it in Stata would require a Stata-side canonical
layout that deviates from upstream. See {cmd:Notes/porting_decisions.md} in
the repository.

{phang}
{opt noconstant} omits the intercept term.


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:one_step_gmm} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{p2col 5 24 28 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of mixture components used{p_end}
{synopt:{cmd:e(nguess)}}multistart attempts requested{p_end}
{synopt:{cmd:e(maxiter)}}BFGS iteration limit per attempt{p_end}
{synopt:{cmd:e(seed)}}seed used for multistart perturbations{p_end}
{synopt:{cmd:e(ll)}}log-likelihood at the winning optimum{p_end}
{synopt:{cmd:e(converged)}}1 if the winning attempt converged, 0 otherwise{p_end}
{synopt:{cmd:e(iterations)}}BFGS iterations used by the winning attempt{p_end}
{synopt:{cmd:e(best_idx)}}index (0-based) of the winning multistart attempt{p_end}
{synopt:{cmd:e(n_finished)}}number of attempts that returned a finite value{p_end}
{synopt:{cmd:e(homoskedastic)}}always 0 (option not supported){p_end}

{p2col 5 24 28 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:one_step_gmm}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(generated_var)}}name of the binary generated covariate{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(vce)}}{cmd:OIM (numerical Hessian)}{p_end}
{synopt:{cmd:e(vcetype)}}{cmd:OIM}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 24 28 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}regression coefficients (intercept first, then alphabetical){p_end}
{synopt:{cmd:e(V)}}top-d block of {it:pinv(numerical Hessian)} at the winning optimum{p_end}


{marker examples}{...}
{title:Examples}

{phang}Synthetic example mirroring the M7 parity fixture (k=2 mixture, bimodal centered residuals):{p_end}

{phang}{cmd:. clear}{p_end}
{phang}{cmd:. set seed 20260522}{p_end}
{phang}{cmd:. set obs 2000}{p_end}
{phang}{cmd:. generate double xstar = runiform() < 0.4}{p_end}
{phang}{cmd:. generate double flip  = runiform() < 0.10}{p_end}
{phang}{cmd:. generate double x1    = cond(flip, 1 - xstar, xstar)}{p_end}
{phang}{cmd:. generate double x2    = rnormal()}{p_end}
{phang}{cmd:. // bimodal residuals: see Notes for the full DGP}{p_end}
{phang}{cmd:. generate double eps   = rnormal()  // simplified for help; see Notes}{p_end}
{phang}{cmd:. generate double y     = 0.2 + 0.7*xstar - 0.3*x2 + eps}{p_end}
{phang}{cmd:. one_step_gmm y x1 x2, generated(x1) k(2) nguess(10) seed(0)}{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
{it:Reduced-dimension optimization at symmetric starts.}  For {opt k} = 2,
the OLS-based starts place four mixture-structure parameters at exactly
zero, which makes the two components mathematically identical. The
likelihood is then analytically flat in those four directions, and Mata's
numerical-derivative pass would abort. {cmd:one_step_gmm} detects this case
and runs the optimizer (and the Hessian step) on the active subspace. This
matches Python's effective behavior, where jaxopt's exact-zero autodiff
gradient on the same subspace keeps those parameters at zero throughout.
The coefficient-block variance is unchanged by the dimension reduction
because the Hessian is block-diagonal between active and symmetric slots at
that point.

{pstd}
{it:Multistart RNG.}  Mata uses Stata's RNG (seeded by {opt seed}); the
Python implementation uses JAX's PRNGKey. The two RNGs are unrelated, so
the {it:i > 0} perturbations differ between implementations and the
{it:best-of-nguess} optimum may land on a different local mode in each. The
deterministic {it:i = 0} attempt is identical across implementations and is
what the shipped parity fixture tests.


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
Help: {help one_step}, {help ols_bca}, {help ols_bcm}
{p_end}
