#ifndef __EVOTSIT5__
#define __EVOTSIT5__

#include "common.h"

/**************************************************************/

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
    /* Jacobian DIAGONAL callback, part of the shared CLASS evolver signature and
       unused here -- as minimum_variation, evaluate_timescale and
       timestep_over_timescale already are. Only evolver_etd consumes it. May be
       null. */
    void (*derivs_diagonal)(double x, double* y, double* diag, void* parameters_and_workspace));

/**************************************************************/

#endif
