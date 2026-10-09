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
 */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "arrays.h"

int main() {
  /* integral of sin(x) from 0 to pi = 2, on 21 points */
  const int n        = 21;
  const double exact = 2.;
  enum { ix, iy, iddy, iint, ncol };
  std::vector<double> x(n), table(n * ncol);
  for (int i = 0; i < n; i++) {
    x[i]                 = _PI_ * i / (n - 1.);
    table[i * ncol + ix] = x[i];
    table[i * ncol + iy] = sin(x[i]);
  }
  array_spline(table.data(), ncol, n, ix, iy, iddy, _SPLINE_EST_DERIV_);

  double trapezoid = 0.;
  for (int i = 0; i < n - 1; i++)
    trapezoid += 0.5 * (table[i * ncol + iy] + table[(i + 1) * ncol + iy]) * (x[i + 1] - x[i]);

  double all_spline, trapzd_or_spline;
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
  printf("spline integrals OK\n");
  return 0;
}
