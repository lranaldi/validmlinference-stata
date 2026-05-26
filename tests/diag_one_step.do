* Diagnostic for Milestone 6 one_step parity failure.

clear all
set more off

local proj "C:/Users/loren/Dropbox/Replication_AIValidInference"
qui cd "`proj'/ValidMLInference-stata/src/mata"
qui do _build.do
qui cd "`proj'"
qui adopath ++ "`proj'/ValidMLInference-stata/src/ado"
qui adopath ++ "`proj'/ValidMLInference-stata/src/mata"
qui mata: mata mlib index

import delimited using "`proj'/ValidMLInference-stata/tests/fixtures/one_step_synthetic/input.csv", asdouble clear

di "n = " _N

mata:
real colvector Y
real matrix    Xhat
real rowvector th0
real scalar    fv, target_idx, homosked, d, i
struct vmli_one_step_pars scalar pars

Y    = st_data(., "y")
Xhat = st_data(., ("x1", "x2"))
Xhat = Xhat, J(rows(Xhat), 1, 1)
d    = cols(Xhat)
target_idx = 1
homosked   = 0

th0 = vmli_one_step_starts(Y, Xhat, target_idx, homosked)'

printf("\n--- starting values (theta0) ---\n")
for (i = 1; i <= cols(th0); i++) printf("  th0[%g] = %14.6g\n", i, th0[i])

pars = vmli_theta_to_pars(th0', d, homosked)
printf("\n--- unpacked pars at theta0 ---\n")
for (i = 1; i <= d; i++) printf("  b[%g] = %14.6g\n", i, pars.b[i])
printf("  (w00, w01, w10, w11) = (%g, %g, %g, %g)\n", pars.w00, pars.w01, pars.w10, pars.w11)
printf("  (sigma0, sigma1)     = (%g, %g)\n", pars.sigma0, pars.sigma1)

fv = .
vmli_one_step_obj(0, th0, Y, Xhat, target_idx, homosked, fv, J(0,0,.), J(0,0,.))
printf("\nobjective at theta0 = %20.10f\n\n", fv)

printf("--- optimize() with tracelevel=params ---\n")

transmorphic S
S = optimize_init()
optimize_init_evaluator(S, &vmli_one_step_obj())
optimize_init_evaluatortype(S, "d0")
optimize_init_argument(S, 1, Y)
optimize_init_argument(S, 2, Xhat)
optimize_init_argument(S, 3, target_idx)
optimize_init_argument(S, 4, homosked)
optimize_init_which(S, "min")
optimize_init_technique(S, "bfgs")
optimize_init_conv_maxiter(S, 500)
optimize_init_conv_ptol(S, 1e-12)
optimize_init_conv_vtol(S, 1e-12)
optimize_init_tracelevel(S, "params")
optimize_init_params(S, th0)

real rowvector th_hat
real scalar rc
rc = _optimize(S)
printf("\nreturn code from _optimize = %g\n", rc)
if (rc == 0) {
    th_hat = optimize_result_params(S)
    printf("converged   = %g\n", optimize_result_converged(S))
    printf("iterations  = %g\n", optimize_result_iterations(S))
    printf("value       = %g\n", optimize_result_value(S))
}

end
