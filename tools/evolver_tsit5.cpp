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
void evolver_tsit5(
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
    /* Part of the shared evolver signature; unused here. See evolver_tsit5.h. */
    void (* /*derivs_diagonal*/)(
        double x, double* y, double* diag, void* parameters_and_workspace)) {
  (void) minimum_variation;
  (void) evaluate_timescale;
  (void) timestep_over_timescale;

  evolver_erk_run<ErkTsitouras54>(derivs,
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
