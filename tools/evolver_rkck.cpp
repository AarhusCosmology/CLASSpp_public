#include "evolver_rkck.h"

#include <limits>
#include <vector>

void evolver_rk(EvolverDerivs derivs,
                double x_ini,
                double x_end,
                double* y,
                int y_size,
                void* parameters_and_workspace_for_derivs,
                const EvolverOptions& options) {
  EvolverOptionsCheck(options, "evolver_rk", x_ini, {EvolverFeature::Timescale});

  const double tolerance                    = options.rtol;
  const double* x_sampling                  = options.x_sampling;
  const int x_size                          = options.x_sampling_size;
  const int* used_in_output                 = options.used_in_output;
  const EvolverOutput output                = options.output;
  const EvolverPrint print_variables        = options.print_variables;
  const EvolverTimescale evaluate_timescale = options.evaluate_timescale;
  const double timestep_over_timescale      = options.timestep_over_timescale;
  /* Every caller filled this with ppr->smallest_allowed_variation, which is
     never read from the input file and so is always DBL_EPSILON. Taken directly
     rather than passed through sixteen arguments. */
  const double minimum_variation = std::numeric_limits<double>::epsilon();

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
}
