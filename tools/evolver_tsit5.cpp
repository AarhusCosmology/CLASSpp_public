#include "evolver_tsit5.h"

#include "evolver_erk_impl.h"

/**
 * Tsitouras 5(4) explicit adaptive Runge-Kutta integrator with 4th-order dense
 * output (Ch. Tsitouras, Comput. Math. Appl. 62 (2011) 770-775). Conforms to the
 * canonical CLASS evolver signature (shared with evolver_ndf15 / evolver_rk).
 *
 * Same shape as Dormand-Prince 5(4): 7 stages, FSAL, so 6 derivative evaluations
 * per accepted step and a 4th-order continuous extension over the existing
 * stages. It differs only in the coefficients, which were optimised for a
 * smaller principal 5th-order truncation-error norm under the first-column
 * simplifying assumption alone.
 *
 * The integrator lives in evolver_erk.cpp; this file only binds the tableau, so
 * tsit5-vs-rkdp45 is a comparison of tableaux and not of two implementations.
 *
 * minimum_variation, evaluate_timescale and timestep_over_timescale are part of
 * the shared signature but unused here: this solver performs its own embedded
 * error control.
 */
void evolver_tsit5(EvolverDerivs derivs,
                   double x_ini,
                   double x_end,
                   double* y,
                   int y_size,
                   void* parameters_and_workspace_for_derivs,
                   const EvolverOptions& options) {
  EvolverOptionsCheck(options,
                      "evolver_tsit5",
                      x_ini,
                      {EvolverFeature::ErkController,
                       EvolverFeature::Stats,
                       EvolverFeature::Histograms,
                       EvolverFeature::AbsTol});

  evolver_erk_run<ErkTsitouras54>(derivs,
                                  x_ini,
                                  x_end,
                                  y,
                                  y_size,
                                  parameters_and_workspace_for_derivs,
                                  options);
}
