#include "evolver_rkck.h"

#include <vector>

int evolver_rk(
    int (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double minimum_variation,
    int (*evaluate_timescale)(double x, void* parameters_and_workspace, double* timescale),
    double timestep_over_timescale,
    double* x_sampling,
    int x_size,
    int (*output)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace),
    int (*print_variables)(double x, double y[], double dy[], void* parameters_and_workspace)) {
  int next_index_x;
  double x1, x2 = 0., timestep, timescale;
  struct generic_integrator_workspace gi;
  std::vector<double> dy(y_size);
  short call_output;

  class_test(x_ini > x_sampling[x_size - 1],
             "called with x=%e, last x_sampling=%e",
             x_ini,
             x_sampling[x_size - 1]);

  next_index_x = 0;

  while (x_sampling[next_index_x] < x_ini)
    next_index_x++;

  initialize_generic_integrator(y_size, &gi);

  x1 = x_ini;

  call_output = false;

  while ((x1 < x_end) && (next_index_x < x_size)) {
    (*evaluate_timescale)(x1, parameters_and_workspace_for_derivs, &timescale);

    timestep = timestep_over_timescale * timescale;

    class_test(fabs(timestep / x1) < minimum_variation,
               "integration step =%e < machine precision : leads either to numerical error or "
               "infinite loop",
               fabs(timestep / x1));

    if (x1 + 2. * timestep < x_sampling[next_index_x]) {
      x2 = x1 + timestep;
    }
    else {
      x2          = x_sampling[next_index_x];
      call_output = true;
    }

    if (x2 > x_end) {
      x2          = x_end;
      call_output = false;
    }

    if (print_variables != nullptr) {
      if (x1 == x_ini) {
        (*derivs)(x1, y, dy.data(), parameters_and_workspace_for_derivs);
      }

      (*print_variables)(x1, y, dy.data(), parameters_and_workspace_for_derivs);
    }

    generic_integrator(derivs,
                       x1,
                       x2,
                       y,
                       parameters_and_workspace_for_derivs,
                       tolerance,
                       x1 * minimum_variation,
                       &gi);

    if (call_output) {
      (*derivs)(x2, y, dy.data(), parameters_and_workspace_for_derivs);

      (*output)(x2, y, dy.data(), next_index_x, parameters_and_workspace_for_derivs);

      call_output = false;

      next_index_x++;
    }

    x1 = x2;
  }

  /* a last call is compulsory to ensure that all quantitites in
     y,dy,parameters_and_workspace_for_derivs are updated to the last
     point in the covered range */
  (*derivs)(x1, y, dy.data(), parameters_and_workspace_for_derivs);

  if (print_variables != nullptr)
    (*print_variables)(x1, y, dy.data(), parameters_and_workspace_for_derivs);

  cleanup_generic_integrator(&gi);

  return _SUCCESS_;
}
