# Root finding for halofit's R_nl and the reionization redshift (2026-10-27)

**Status 2026-10-27:** implemented (branch `perf/root-finders`, on `perf/spline-integration-weights`). Results
move within the tolerances of the solvers they replace. Evidence tags: [M] measured here, [C] stated by the
code.

## Problem

Two bisections dominate stages that CLASS++ spends much longer in than CAMB
(`benchmarks/camb_anatomy/FINDINGS.md`, item F) [M, macOS `sample`, 1 thread]:

* **The reionization redshift.** With `tau_reio` as input, `thermodynamics_reionization` bisected z_reio until
  the bracket in tau was narrower than `reionization_optical_depth_tol` (1e-4) [C]. Every step rebuilt the
  whole reionization history, an adaptive sampling with a background interpolation at every trial step. About
  17 histories per run: 12 of the ~20 ms of background + thermodynamics.
* **halofit's R_nl**, where sigma(R) = 1, at every time step: a bisection in log R to |sigma - 1| <=
  `halofit_tol_sigma` (1e-6), twenty-odd integrals over k, then two more for the derivatives [C].

## Approach

* **Reionization:** regula falsi with the Illinois modification on tau(z_reio) - tau_reio. It keeps the
  bracket, so it cannot fail where bisection would not, and converges superlinearly on the smooth, growing
  tau(z). It stops when the tau of the last sampled history is within the tolerance (relative 1e-4), and that
  history is the one kept, as before. tau = 0 at z_reio = 0 stays a bracket value, not a sampled one.
* **halofit:** Newton's method on ln sigma^2 as a function of ln R. Its derivative, -sum2/sum1, is the integral
  that halofit needs anyway, and its integrand is the sigma^2 integrand times 2 (kR)^2. So one pass over k with
  one exp gives sigma^2 and both derivative integrals, and the final step already has all three. A step that
  would leave the bracket bisects instead. The stopping rule is unchanged.

## Verification

* Against the branch below, 8 codecmp LCDM draws, l_max_scalars 3000 [M]: halofit P_NL moves by <= 7e-6
  (z = 0-3, k = 0.01-10/Mpc), within `halofit_tol_sigma`. The lensed C_l move by <= 0.019 CV (TT, l = 2-3000).
  That is the reionization tolerance: 1e-4 in tau is 1.1e-5 in the amplitude e^(-2 tau).
* Which is closer: both solvers run to a tolerance of 1e-9 agree to <= 2.7e-4 CV (the halofit change). At the
  default tolerance, the Illinois history is <= 0.006 CV from its converged one, against <= 0.014 for
  bisection [M].
* Stages, 1 thread, master precision, l_max_scalars 3000, best of 3, 4 draws [M]: background + thermodynamics
  29.5 -> 17 ms, non-linear 13 -> 11 ms.

## Not done here

What is left in both stages is table interpolation, pow and exp [M, `sample`]: the background at every trial
step of the reionization sampling, and the linear P(k) interpolated onto halofit's integration grid at every
time step. Both codes run RECFAST in this benchmark (`benchmarks/codecmp/models.py`), where CAMB takes 2.5 ms
for background + thermodynamics [M].
