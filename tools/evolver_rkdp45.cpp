#include "evolver_rkdp45.h"

#include "evolver_erk_impl.h"

/**
 * Dormand-Prince 5(4) explicit adaptive Runge-Kutta integrator with 4th-order
 * dense output. Ported from branch PhD2024-EBH. Conforms to the canonical CLASS
 * evolver signature (shared with evolver_ndf15 / evolver_rk).
 *
 * The integrator itself now lives in evolver_erk.cpp, which drives any embedded
 * explicit pair from a tableau; this file only binds the Dormand-Prince tableau
 * to the CLASS signature. evolver_tsit5 binds a different one to the SAME core,
 * which is what makes a comparison between them a comparison of tableaux.
 *
 * minimum_variation, evaluate_timescale and timestep_over_timescale are part of
 * the shared signature but unused here: this solver performs its own embedded
 * error control.
 */
void evolver_rkdp45(EvolverDerivs derivs,
                    double x_ini,
                    double x_end,
                    double* y,
                    int y_size,
                    void* parameters_and_workspace_for_derivs,
                    const EvolverOptions& options) {
  EvolverOptionsCheck(options,
                      "evolver_rkdp45",
                      x_ini,
                      {EvolverFeature::ErkController,
                       EvolverFeature::Stats,
                       EvolverFeature::Histograms,
                       EvolverFeature::AbsTol});

  evolver_erk_run<ErkDormandPrince45>(derivs,
                                      x_ini,
                                      x_end,
                                      y,
                                      y_size,
                                      parameters_and_workspace_for_derivs,
                                      options);
}
