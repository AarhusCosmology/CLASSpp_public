/** @file spline_integral_test.cpp
 * The spline integration helpers of tools/arrays.cpp against an exact integral.
 *
 * On one interval, the integral of a cubic spline is
 *   h/2 (y_i + y_{i+1}) - h^3/24 (y''_i + y''_{i+1}),
 * the trapezoidal rule minus a curvature correction. CLASS++ had a + sign on
 * the correction in all three helpers, which doubles the trapezoidal error
 * instead of removing it; class_public fixed it after v3.2 (credits: C.
 * Radermacher). Here the spline result must beat the trapezoidal rule by a
 * large factor on a smooth function.
 *
 * array_spline_integration_weights must reproduce array_spline followed by
 * array_integrate_all_trapzd_or_spline as a dot product, on uneven grids, for
 * any start of the spline part, both end conditions, down to three points.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "arrays.h"

int main() {
  /* integral of sin(x) from 0 to pi = 2, on 21 points */
  const int n        = 21;
  const double exact = 2.;
  enum { ix, iy, iddy, iint, ncol };
  std::vector<double> x(n);
  std::vector<double> table(n * ncol);
  for (int i = 0; i < n; i++) {
    x[i]                 = _PI_ * i / (n - 1.);
    table[i * ncol + ix] = x[i];
    table[i * ncol + iy] = sin(x[i]);
  }
  array_spline(table.data(), ncol, n, ix, iy, iddy, _SPLINE_EST_DERIV_);

  double trapezoid = 0.;
  for (int i = 0; i < n - 1; i++) {
    trapezoid += 0.5 * (table[i * ncol + iy] + table[(i + 1) * ncol + iy]) * (x[i + 1] - x[i]);
  }

  double all_spline;
  double trapzd_or_spline;
  array_integrate_all_spline(table.data(), ncol, n, ix, iy, iddy, &all_spline);
  array_integrate_all_trapzd_or_spline(table.data(), ncol, n, 0, ix, iy, iddy, &trapzd_or_spline);
  array_integrate_spline_table_line_to_line(x.data(), n, table.data(), ncol, iy, iddy, iint);
  const double line_to_line = table[(n - 1) * ncol + iint];

  const double err_trapezoid = fabs(trapezoid - exact);
  printf("trapezoid %.3e, all_spline %.3e, trapzd_or_spline %.3e, line_to_line %.3e\n",
         err_trapezoid,
         fabs(all_spline - exact),
         fabs(trapzd_or_spline - exact),
         fabs(line_to_line - exact));

  /* the spline integrals must be at least 100 times better than the trapezoidal rule */
  assert(fabs(all_spline - exact) < 1e-2 * err_trapezoid);
  assert(fabs(trapzd_or_spline - exact) < 1e-2 * err_trapezoid);
  assert(fabs(line_to_line - exact) < 1e-2 * err_trapezoid);

  /* integration weights against the spline integral they stand for */
  std::mt19937 rng(20261027);
  std::uniform_real_distribution<double> uniform(0., 1.);
  double worst = 0.;
  for (const short mode : {_SPLINE_EST_DERIV_, _SPLINE_NATURAL_}) {
    for (const int m : {3, 4, 5, 10, 200, 2263}) {
      for (const double start_fraction : {0., 0.37, 1.}) {
        enum { jx, jy, jddy, mcol };
        std::vector<double> grid(m);
        std::vector<double> data(m * mcol);
        double xi = 1e-4;
        for (int i = 0; i < m; i++) {
          xi                  += 1e-3 * (0.2 + uniform(rng)) * (1. + 0.01 * i);
          grid[i]              = xi;
          data[i * mcol + jx]  = xi;
          data[i * mcol + jy]  = cos(40. * xi) * exp(-xi) + 0.3 * (uniform(rng) - 0.5);
        }
        const int start = std::min(static_cast<int>(start_fraction * (m - 1)), m - 1);
        array_spline(data.data(), mcol, m, jx, jy, jddy, mode);
        double direct;
        array_integrate_all_trapzd_or_spline(data.data(), mcol, m, start, jx, jy, jddy, &direct);

        std::vector<double> w(m);
        array_spline_integration_weights(grid.data(), m, start, mode, w.data());
        double dot   = 0.;
        double scale = 0.;
        for (int i = 0; i < m; i++) {
          dot   += w[i] * data[i * mcol + jy];
          scale += fabs(w[i] * data[i * mcol + jy]);
        }
        worst = std::max(worst, fabs(dot - direct) / scale);
      }
    }
  }
  printf("integration weights against the spline integral: worst %.2e of sum |w y|\n", worst);
  assert(worst < 1e-13);

  printf("spline integrals OK\n");
  return 0;
}
