# Spline integrals as precomputed weights (2026-10-27)

**Status 2026-10-27:** implemented (branch `perf/spline-integration-weights`). Results equal to rounding.
Evidence tags: [M] measured here, [C] stated by the code.

## Problem

Two of the stages that CLASS++ spends much longer in than CAMB (`benchmarks/camb_anatomy/FINDINGS.md`, item F)
repeat the same spline integral over a fixed grid [M, macOS `sample`, 1 thread, l_max_scalars 3000]:

* **C_l** (`SpectraModule::spectra_compute_cl`): for every l and every spectrum type, `array_spline` over the
  ~2300 k values of the transfer grid, then `array_integrate_all_trapzd_or_spline`. The spline is about two
  thirds of the stage.
* **halofit** (`NonlinearModule::nonlinear_halofit_integrate`): a natural spline and its integral over the same
  k grid for every trial radius of the bisection for R_nl (twenty-odd per time step), and for the two
  derivative integrals.

## Approach

The spline integral is linear in the integrand y: the second derivatives m solve A m = B y, with A and B set
by the grid alone, and the integral is t.y - c.m (t the trapezoidal weights, c the h^3/24 corrections). So it
equals w.y with w = t - B^T z, where A^T z = c: one transposed tridiagonal solve per grid [C].

`array_spline_integration_weights(x, n, index_start_spline, spline_mode, w)` (`tools/arrays.cpp`) builds w
for both end conditions in use (`_SPLINE_EST_DERIV_`, `_SPLINE_NATURAL_`) and for a trapezoidal part below
`index_start_spline`, matching `array_spline` followed by `array_integrate_all_trapzd_or_spline`.

* The spectra module builds the weights of the transfer k grid (and of the full-Limber grid of C_l^phiphi)
  once per mode in `spectra_cls`. In the closed case the trapezoidal part below the flat approximation stays.
  Each C_l is then a dot product, and the integrand table loses its second-derivative columns.
* halofit builds the weights of its integration grid once per time step; each integral is a sum inside the
  integrand loop.

## Verification

* `test-spline-integral`: on uneven random grids of 3 to 2263 points, both end conditions and three starts of
  the spline part, the weights reproduce the spline integral to 2.7e-16 of sum |w y| [M].
* Against master, largest difference relative to the largest value of each spectrum [M]: lensed LCDM to l =
  3000, closed and open universes, tensors, number counts and lensing in two bins, CDM isocurvature: <= 2e-14
  (lensed BB; everything else <= 7e-15). On master with #477, <= 4e-13: the acoustic spline's envelope fit
  passes the rounding on to TT, EE and TE. halofit P_NL(k, z) at z = 0-4, LCDM, w0wa with Pk_equal, a 0.3 eV
  neutrino: <= 5e-15.
* Stages, 1 thread, master precision, l_max_scalars 3000, 4 draws, best of 3 [M]: non-linear 25 -> 13.5 ms,
  C_l 17.5 -> 4.5 ms, run 0.765 -> 0.741 s (-3.1%).

## Not done here

* halofit still bisects for R_nl: twenty-odd integrals per time step, each with an exp per k. Newton's method
  on ln sigma(ln R) needs the derivative integral, whose integrand is the sigma integrand times 2 (kR)^2 (one
  exp for both), and should converge in ~4 steps. That moves P_NL within `halofit_tol_sigma`, so it is a
  separate change.
* HMcode's sigma(R) integrals (`nonlinear_sigmas`) use the same spline integral and could take the same
  weights.
* Found on the way, on master: halofit with `pk_eq = yes`, lensing and `z_max_pk = 4` fails in
  `array_search_bracket` (tau below the Pk_equal table) [M].
