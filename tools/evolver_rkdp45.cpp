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
void evolver_rkdp45(
    void (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double minimum_variation,
    void (*evaluate_timescale)(double x, void* parameters_and_workspace, double* timescale),
    double timestep_over_timescale,
    double* x_sampling,
    int x_size,
    void (*output)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace),
    void (*print_variables)(double x, double y[], double dy[], void* parameters_and_workspace),
    /* Part of the shared evolver signature; unused here. See evolver_rkdp45.h. */
    void (* /*derivs_diagonal*/)(
        double x, double* y, double* diag, void* parameters_and_workspace)) {
  (void) minimum_variation;
  (void) evaluate_timescale;
  (void) timestep_over_timescale;

  evolver_erk_run<ErkDormandPrince45>(derivs,
                                      x_ini,
                                      x_end,
                                      y,
                                      used_in_output,
                                      y_size,
                                      parameters_and_workspace_for_derivs,
                                      tolerance,
                                      x_sampling,
                                      x_size,
                                      output,
                                      print_variables);
}
