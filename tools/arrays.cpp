/**
 * module with tools for manipulating arrays
 * Julien Lesgourgues, 18.04.2010
 */

#include "arrays.h"

#include <algorithm>

void array_derive_spline(const double* x_array,
                         int n_lines,
                         double* array,
                         const double* array_splined,
                         int n_columns,
                         int index_y,
                         int index_dydx) {
  int i;

  double h;

  class_test(index_dydx == index_y,
             "Output column %d must differ from input columns %d",
             index_dydx,
             index_y);

  class_test(n_lines < 2, "no possible derivation with less than two lines");

  for (i = 0; i < n_lines - 1; i++) {
    h = x_array[i + 1] - x_array[i];
    class_test(h == 0, "h=0, stop to avoid division by zero");

    array[i * n_columns + index_dydx] =
        (array[(i + 1) * n_columns + index_y] - array[i * n_columns + index_y]) / h -
        h / 6. *
            (array_splined[(i + 1) * n_columns + index_y] +
             2. * array_splined[i * n_columns + index_y]);
  }

  h = x_array[n_lines - 1] - x_array[n_lines - 2];

  array[(n_lines - 1) * n_columns + index_dydx] =
      (array[(n_lines - 1) * n_columns + index_y] - array[(n_lines - 2) * n_columns + index_y]) /
          h +
      h / 6. *
          (2. * array_splined[(n_lines - 1) * n_columns + index_y] +
           array_splined[(n_lines - 2) * n_columns + index_y]);
}

void array_derive_spline_table_line_to_line(const double* x_array,
                                            int n_lines,
                                            double* array,
                                            int n_columns,
                                            int index_y,
                                            int index_ddy,
                                            int index_dy) {
  int i;

  double h;

  class_test(index_ddy == index_y,
             "Output column %d must differ from input columns %d",
             index_ddy,
             index_y);

  class_test(index_ddy == index_dy,
             "Output column %d must differ from input columns %d",
             index_ddy,
             index_dy);

  class_test(n_lines < 2, "no possible derivation with less than two lines");

  for (i = 0; i < n_lines - 1; i++) {
    h = x_array[i + 1] - x_array[i];
    class_test(h == 0, "h=0, stop to avoid division by zero")

        array[i * n_columns + index_dy] =
            (array[(i + 1) * n_columns + index_y] - array[i * n_columns + index_y]) / h -
            h / 6. *
                (array[(i + 1) * n_columns + index_ddy] + 2. * array[i * n_columns + index_ddy]);
  }

  h = x_array[n_lines - 1] - x_array[n_lines - 2];

  array[(n_lines - 1) * n_columns + index_dy] =
      (array[(n_lines - 1) * n_columns + index_y] - array[(n_lines - 2) * n_columns + index_y]) /
          h +
      h / 6. *
          (2. * array[(n_lines - 1) * n_columns + index_ddy] +
           array[(n_lines - 2) * n_columns + index_ddy]);
}

void array_integrate_spline_table_line_to_line(const double* x_array,
                                               int n_lines,
                                               double* array,
                                               int n_columns,
                                               int index_y,
                                               int index_ddy,
                                               int index_inty) {
  int i;

  double h;

  *(array + 0 * n_columns + index_inty) = 0.;

  for (i = 0; i < n_lines - 1; i++) {
    h = (x_array[i + 1] - x_array[i]);

    *(array + (i + 1) * n_columns + index_inty) =
        *(array + i * n_columns + index_inty) +
        (array[i * n_columns + index_y] + array[(i + 1) * n_columns + index_y]) * h / 2. +
        (array[i * n_columns + index_ddy] + array[(i + 1) * n_columns + index_ddy]) * h * h * h /
            24.;
  }
}

void array_spline(double* array,
                  int n_columns,
                  int n_lines,
                  int index_x, /** from 0 to (n_columns-1) */
                  int index_y,
                  int index_ddydx2,
                  short spline_mode) {
  int i, k;
  double p, qn, sig, un;
  std::vector<double> vec_u(n_lines - 1);
  double* u = vec_u.data();
  double dy_first;
  double dy_last;

  class_test(n_lines < 3, "n_lines=%d, while routine needs n_lines >= 3");

  if (spline_mode == _SPLINE_NATURAL_) {
    *(array + 0 * n_columns + index_ddydx2) = u[0] = 0.0;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_first = ((*(array + 2 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) *
                      (*(array + 2 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) *
                      (*(array + 1 * n_columns + index_y) - *(array + 0 * n_columns + index_y)) -
                  (*(array + 1 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) *
                      (*(array + 1 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) *
                      (*(array + 2 * n_columns + index_y) - *(array + 0 * n_columns + index_y))) /
                 ((*(array + 2 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) *
                  (*(array + 1 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) *
                  (*(array + 2 * n_columns + index_x) - *(array + 1 * n_columns + index_x)));

      *(array + 0 * n_columns + index_ddydx2) = -0.5;

      u[0] = (3. / (*(array + 1 * n_columns + index_x) - *(array + 0 * n_columns + index_x))) *
             ((*(array + 1 * n_columns + index_y) - *(array + 0 * n_columns + index_y)) /
                  (*(array + 1 * n_columns + index_x) - *(array + 0 * n_columns + index_x)) -
              dy_first);
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  for (i = 1; i < n_lines - 1; i++) {
    sig = (*(array + i * n_columns + index_x) - *(array + (i - 1) * n_columns + index_x)) /
          (*(array + (i + 1) * n_columns + index_x) - *(array + (i - 1) * n_columns + index_x));

    p = sig * *(array + (i - 1) * n_columns + index_ddydx2) + 2.0;
    *(array + i * n_columns + index_ddydx2) = (sig - 1.0) / p;
    u[i] = (*(array + (i + 1) * n_columns + index_y) - *(array + i * n_columns + index_y)) /
               (*(array + (i + 1) * n_columns + index_x) - *(array + i * n_columns + index_x)) -
           (*(array + i * n_columns + index_y) - *(array + (i - 1) * n_columns + index_y)) /
               (*(array + i * n_columns + index_x) - *(array + (i - 1) * n_columns + index_x));
    u[i] =
        (6.0 * u[i] /
             (*(array + (i + 1) * n_columns + index_x) - *(array + (i - 1) * n_columns + index_x)) -
         sig * u[i - 1]) /
        p;
  }

  if (spline_mode == _SPLINE_NATURAL_) {
    qn = 0.;
    un = 0.;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_last = ((*(array + (n_lines - 3) * n_columns + index_x) -
                  *(array + (n_lines - 1) * n_columns + index_x)) *
                     (*(array + (n_lines - 3) * n_columns + index_x) -
                      *(array + (n_lines - 1) * n_columns + index_x)) *
                     (*(array + (n_lines - 2) * n_columns + index_y) -
                      *(array + (n_lines - 1) * n_columns + index_y)) -
                 (*(array + (n_lines - 2) * n_columns + index_x) -
                  *(array + (n_lines - 1) * n_columns + index_x)) *
                     (*(array + (n_lines - 2) * n_columns + index_x) -
                      *(array + (n_lines - 1) * n_columns + index_x)) *
                     (*(array + (n_lines - 3) * n_columns + index_y) -
                      *(array + (n_lines - 1) * n_columns + index_y))) /
                ((*(array + (n_lines - 3) * n_columns + index_x) -
                  *(array + (n_lines - 1) * n_columns + index_x)) *
                 (*(array + (n_lines - 2) * n_columns + index_x) -
                  *(array + (n_lines - 1) * n_columns + index_x)) *
                 (*(array + (n_lines - 3) * n_columns + index_x) -
                  *(array + (n_lines - 2) * n_columns + index_x)));

      qn = 0.5;
      un = (3. / (*(array + (n_lines - 1) * n_columns + index_x) -
                  *(array + (n_lines - 2) * n_columns + index_x))) *
           (dy_last - (*(array + (n_lines - 1) * n_columns + index_y) -
                       *(array + (n_lines - 2) * n_columns + index_y)) /
                          (*(array + (n_lines - 1) * n_columns + index_x) -
                           *(array + (n_lines - 2) * n_columns + index_x)));
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  *(array + (n_lines - 1) * n_columns + index_ddydx2) =
      (un - qn * u[n_lines - 2]) / (qn * *(array + (n_lines - 2) * n_columns + index_ddydx2) + 1.0);

  for (k = n_lines - 2; k >= 0; k--)
    *(array + k * n_columns + index_ddydx2) = *(array + k * n_columns + index_ddydx2) *
                                                  *(array + (k + 1) * n_columns + index_ddydx2) +
                                              u[k];
}

void array_spline_table_line_to_line(const double* x, /* vector of size x_size */
                                     int n_lines,
                                     double* array,
                                     int n_columns,
                                     int index_y,
                                     int index_ddydx2,
                                     short spline_mode) {
  int i, k;
  double p, qn, sig, un;
  std::vector<double> vec_u(n_lines - 1);
  double* u = vec_u.data();
  double dy_first;
  double dy_last;

  if (spline_mode == _SPLINE_NATURAL_) {
    *(array + 0 * n_columns + index_ddydx2) = u[0] = 0.0;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_first = ((x[2] - x[0]) * (x[2] - x[0]) *
                      (*(array + 1 * n_columns + index_y) - *(array + 0 * n_columns + index_y)) -
                  (x[1] - x[0]) * (x[1] - x[0]) *
                      (*(array + 2 * n_columns + index_y) - *(array + 0 * n_columns + index_y))) /
                 ((x[2] - x[0]) * (x[1] - x[0]) * (x[2] - x[1]));
      *(array + 0 * n_columns + index_ddydx2) = -0.5;
      u[0] = (3. / (x[1] - x[0])) *
             ((*(array + 1 * n_columns + index_y) - *(array + 0 * n_columns + index_y)) /
                  (x[1] - x[0]) -
              dy_first);
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  for (i = 1; i < n_lines - 1; i++) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);

    p = sig * *(array + (i - 1) * n_columns + index_ddydx2) + 2.0;
    *(array + i * n_columns + index_ddydx2) = (sig - 1.0) / p;
    u[i] = (*(array + (i + 1) * n_columns + index_y) - *(array + i * n_columns + index_y)) /
               (x[i + 1] - x[i]) -
           (*(array + i * n_columns + index_y) - *(array + (i - 1) * n_columns + index_y)) /
               (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  if (spline_mode == _SPLINE_NATURAL_) {
    qn = 0.;
    un = 0.;
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      dy_last = ((x[n_lines - 3] - x[n_lines - 1]) * (x[n_lines - 3] - x[n_lines - 1]) *
                     (*(array + (n_lines - 2) * n_columns + index_y) -
                      *(array + (n_lines - 1) * n_columns + index_y)) -
                 (x[n_lines - 2] - x[n_lines - 1]) * (x[n_lines - 2] - x[n_lines - 1]) *
                     (*(array + (n_lines - 3) * n_columns + index_y) -
                      *(array + (n_lines - 1) * n_columns + index_y))) /
                ((x[n_lines - 3] - x[n_lines - 1]) * (x[n_lines - 2] - x[n_lines - 1]) *
                 (x[n_lines - 3] - x[n_lines - 2]));
      qn      = 0.5;
      un      = (3. / (x[n_lines - 1] - x[n_lines - 2])) *
                (dy_last - (*(array + (n_lines - 1) * n_columns + index_y) -
                            *(array + (n_lines - 2) * n_columns + index_y)) /
                               (x[n_lines - 1] - x[n_lines - 2]));
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  *(array + (n_lines - 1) * n_columns + index_ddydx2) =
      (un - qn * u[n_lines - 2]) / (qn * *(array + (n_lines - 2) * n_columns + index_ddydx2) + 1.0);

  for (k = n_lines - 2; k >= 0; k--)
    *(array + k * n_columns + index_ddydx2) = *(array + k * n_columns + index_ddydx2) *
                                                  *(array + (k + 1) * n_columns + index_ddydx2) +
                                              u[k];
}

void array_spline_table_lines(const double* x, /* vector of size x_size */
                              int x_size,
                              const double* y_array, /* array of size x_size*y_size with elements
						  y_array[index_x*y_size+index_y] */
                              int y_size,
                              double* ddy_array, /* array of size x_size*y_size */
                              short spline_mode) {
  std::vector<double> vec_u((x_size - 1) * y_size);
  std::vector<double> vec_p(y_size);
  std::vector<double> vec_qn(y_size);
  std::vector<double> vec_un(y_size);
  double* u  = vec_u.data();
  double* p  = vec_p.data();
  double* qn = vec_qn.data();
  double* un = vec_un.data();
  double sig;
  int index_x;
  int index_y;
  double dy_first;
  double dy_last;

  if (x_size == 2)
    spline_mode =
        _SPLINE_NATURAL_;  // in the case of only 2 x-values, only the natural spline method is appropriate, for _SPLINE_EST_DERIV_ at least 3 x-values are needed.

  index_x = 0;

  if (spline_mode == _SPLINE_NATURAL_) {
    for (index_y = 0; index_y < y_size; index_y++) {
      ddy_array[index_x * y_size + index_y] = u[index_x * y_size + index_y] = 0.0;
    }
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      for (index_y = 0; index_y < y_size; index_y++) {
        dy_first = ((x[2] - x[0]) * (x[2] - x[0]) *
                        (y_array[1 * y_size + index_y] - y_array[0 * y_size + index_y]) -
                    (x[1] - x[0]) * (x[1] - x[0]) *
                        (y_array[2 * y_size + index_y] - y_array[0 * y_size + index_y])) /
                   ((x[2] - x[0]) * (x[1] - x[0]) * (x[2] - x[1]));

        ddy_array[index_x * y_size + index_y] = -0.5;

        u[index_x * y_size + index_y] =
            (3. / (x[1] - x[0])) *
            ((y_array[1 * y_size + index_y] - y_array[0 * y_size + index_y]) / (x[1] - x[0]) -
             dy_first);
      }
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  for (index_x = 1; index_x < x_size - 1; index_x++) {
    sig = (x[index_x] - x[index_x - 1]) / (x[index_x + 1] - x[index_x - 1]);

    for (index_y = 0; index_y < y_size; index_y++) {
      p[index_y] = sig * ddy_array[(index_x - 1) * y_size + index_y] + 2.0;

      ddy_array[index_x * y_size + index_y] = (sig - 1.0) / p[index_y];

      u[index_x * y_size + index_y] =
          (y_array[(index_x + 1) * y_size + index_y] - y_array[index_x * y_size + index_y]) /
              (x[index_x + 1] - x[index_x]) -
          (y_array[index_x * y_size + index_y] - y_array[(index_x - 1) * y_size + index_y]) /
              (x[index_x] - x[index_x - 1]);

      u[index_x * y_size + index_y] = (6.0 * u[index_x * y_size + index_y] /
                                           (x[index_x + 1] - x[index_x - 1]) -
                                       sig * u[(index_x - 1) * y_size + index_y]) /
                                      p[index_y];
    }
  }

  if (spline_mode == _SPLINE_NATURAL_) {
    for (index_y = 0; index_y < y_size; index_y++) {
      qn[index_y] = un[index_y] = 0.0;
    }
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      for (index_y = 0; index_y < y_size; index_y++) {
        dy_last = ((x[x_size - 3] - x[x_size - 1]) * (x[x_size - 3] - x[x_size - 1]) *
                       (y_array[(x_size - 2) * y_size + index_y] -
                        y_array[(x_size - 1) * y_size + index_y]) -
                   (x[x_size - 2] - x[x_size - 1]) * (x[x_size - 2] - x[x_size - 1]) *
                       (y_array[(x_size - 3) * y_size + index_y] -
                        y_array[(x_size - 1) * y_size + index_y])) /
                  ((x[x_size - 3] - x[x_size - 1]) * (x[x_size - 2] - x[x_size - 1]) *
                   (x[x_size - 3] - x[x_size - 2]));

        qn[index_y] = 0.5;

        un[index_y] = (3. / (x[x_size - 1] - x[x_size - 2])) *
                      (dy_last - (y_array[(x_size - 1) * y_size + index_y] -
                                  y_array[(x_size - 2) * y_size + index_y]) /
                                     (x[x_size - 1] - x[x_size - 2]));
      }
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  index_x = x_size - 1;

  for (index_y = 0; index_y < y_size; index_y++) {
    ddy_array[index_x * y_size + index_y] =
        (un[index_y] - qn[index_y] * u[(index_x - 1) * y_size + index_y]) /
        (qn[index_y] * ddy_array[(index_x - 1) * y_size + index_y] + 1.0);
  }

  for (index_x = x_size - 2; index_x >= 0; index_x--) {
    for (index_y = 0; index_y < y_size; index_y++) {
      ddy_array[index_x * y_size + index_y] = ddy_array[index_x * y_size + index_y] *
                                                  ddy_array[(index_x + 1) * y_size + index_y] +
                                              u[index_x * y_size + index_y];
    }
  }
}

void array_spline_table_columns(const double* x, /* vector of size x_size */
                                int x_size,
                                const double* y_array, /* array of size x_size*y_size with elements
					  y_array[index_y*x_size+index_x] */
                                int y_size,
                                double* ddy_array, /* array of size x_size*y_size */
                                short spline_mode) {
  std::vector<double> vec_u((x_size - 1) * y_size);
  std::vector<double> vec_p(y_size);
  std::vector<double> vec_qn(y_size);
  std::vector<double> vec_un(y_size);
  double* u  = vec_u.data();
  double* p  = vec_p.data();
  double* qn = vec_qn.data();
  double* un = vec_un.data();
  double sig;
  int index_x;
  int index_y;
  double dy_first;
  double dy_last;

  if (x_size == 2)
    spline_mode =
        _SPLINE_NATURAL_;  // in the case of only 2 x-values, only the natural spline method is appropriate, for _SPLINE_EST_DERIV_ at least 3 x-values are needed.

  index_x = 0;

  if (spline_mode == _SPLINE_NATURAL_) {
    for (index_y = 0; index_y < y_size; index_y++) {
      ddy_array[index_y * x_size + index_x] = 0.0;
      u[index_x * y_size + index_y]         = 0.0;
    }
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      class_test(x[2] - x[0] == 0., "x[2]=%g, x[0]=%g, stop to avoid seg fault", x[2], x[0]);
      class_test(x[1] - x[0] == 0., "x[1]=%g, x[0]=%g, stop to avoid seg fault", x[1], x[0]);
      class_test(x[2] - x[1] == 0., "x[2]=%g, x[1]=%g, stop to avoid seg fault", x[2], x[1]);

      for (index_y = 0; index_y < y_size; index_y++) {
        dy_first = ((x[2] - x[0]) * (x[2] - x[0]) *
                        (y_array[index_y * x_size + 1] - y_array[index_y * x_size + 0]) -
                    (x[1] - x[0]) * (x[1] - x[0]) *
                        (y_array[index_y * x_size + 2] - y_array[index_y * x_size + 0])) /
                   ((x[2] - x[0]) * (x[1] - x[0]) * (x[2] - x[1]));

        ddy_array[index_y * x_size + index_x] = -0.5;

        u[index_x * y_size + index_y] =
            (3. / (x[1] - x[0])) *
            ((y_array[index_y * x_size + 1] - y_array[index_y * x_size + 0]) / (x[1] - x[0]) -
             dy_first);
      }
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  for (index_x = 1; index_x < x_size - 1; index_x++) {
    sig = (x[index_x] - x[index_x - 1]) / (x[index_x + 1] - x[index_x - 1]);

    for (index_y = 0; index_y < y_size; index_y++) {
      p[index_y] = sig * ddy_array[index_y * x_size + (index_x - 1)] + 2.0;

      ddy_array[index_y * x_size + index_x] = (sig - 1.0) / p[index_y];

      u[index_x * y_size + index_y] =
          (y_array[index_y * x_size + (index_x + 1)] - y_array[index_y * x_size + index_x]) /
              (x[index_x + 1] - x[index_x]) -
          (y_array[index_y * x_size + index_x] - y_array[index_y * x_size + (index_x - 1)]) /
              (x[index_x] - x[index_x - 1]);

      u[index_x * y_size + index_y] = (6.0 * u[index_x * y_size + index_y] /
                                           (x[index_x + 1] - x[index_x - 1]) -
                                       sig * u[(index_x - 1) * y_size + index_y]) /
                                      p[index_y];
    }
  }

  if (spline_mode == _SPLINE_NATURAL_) {
    for (index_y = 0; index_y < y_size; index_y++) {
      qn[index_y] = un[index_y] = 0.0;
    }
  }
  else {
    if (spline_mode == _SPLINE_EST_DERIV_) {
      for (index_y = 0; index_y < y_size; index_y++) {
        dy_last = ((x[x_size - 3] - x[x_size - 1]) * (x[x_size - 3] - x[x_size - 1]) *
                       (y_array[index_y * x_size + (x_size - 2)] -
                        y_array[index_y * x_size + (x_size - 1)]) -
                   (x[x_size - 2] - x[x_size - 1]) * (x[x_size - 2] - x[x_size - 1]) *
                       (y_array[index_y * x_size + (x_size - 3)] -
                        y_array[index_y * x_size + (x_size - 1)])) /
                  ((x[x_size - 3] - x[x_size - 1]) * (x[x_size - 2] - x[x_size - 1]) *
                   (x[x_size - 3] - x[x_size - 2]));

        qn[index_y] = 0.5;

        un[index_y] = (3. / (x[x_size - 1] - x[x_size - 2])) *
                      (dy_last - (y_array[index_y * x_size + (x_size - 1)] -
                                  y_array[index_y * x_size + (x_size - 2)]) /
                                     (x[x_size - 1] - x[x_size - 2]));
      }
    }
    else {
      class_stop("Spline mode not identified: %d", spline_mode);
    }
  }

  index_x = x_size - 1;

  for (index_y = 0; index_y < y_size; index_y++) {
    ddy_array[index_y * x_size + index_x] =
        (un[index_y] - qn[index_y] * u[(index_x - 1) * y_size + index_y]) /
        (qn[index_y] * ddy_array[index_y * x_size + (index_x - 1)] + 1.0);
  }

  for (index_x = x_size - 2; index_x >= 0; index_x--) {
    for (index_y = 0; index_y < y_size; index_y++) {
      ddy_array[index_y * x_size + index_x] = ddy_array[index_y * x_size + index_x] *
                                                  ddy_array[index_y * x_size + (index_x + 1)] +
                                              u[index_x * y_size + index_y];
    }
  }
}

void array_spline_table_columns2(const double* x, /* vector of size x_size */
                                 int x_size,
                                 const double* y_array, /* array of size x_size*y_size with elements
					  y_array[index_y*x_size+index_x] */
                                 int y_size,
                                 double* ddy_array, /* array of size x_size*y_size */
                                 short spline_mode) {
  std::vector<double> vec_u((x_size - 1) * y_size);
  std::vector<double> vec_p(y_size);
  std::vector<double> vec_qn(y_size);
  std::vector<double> vec_un(y_size);
  double* u  = vec_u.data();
  double* p  = vec_p.data();
  double* qn = vec_qn.data();
  double* un = vec_un.data();
  double sig;
  int index_x;
  int index_y;
  double dy_first;
  double dy_last;

  if (x_size == 2)
    spline_mode =
        _SPLINE_NATURAL_;  // in the case of only 2 x-values, only the natural spline method is appropriate, for _SPLINE_EST_DERIV_ 3 x-values are needed.

#pragma omp parallel shared(x,               \
                                x_size,      \
                                y_array,     \
                                y_size,      \
                                ddy_array,   \
                                spline_mode, \
                                p,           \
                                qn,          \
                                un,          \
                                u) private(index_y, index_x, sig, dy_first, dy_last)
  {
#pragma omp for schedule(dynamic)

    for (index_y = 0; index_y < y_size; index_y++) {
      if (spline_mode == _SPLINE_NATURAL_) {
        ddy_array[index_y * x_size + 0] = 0.0;
        u[0 * y_size + index_y]         = 0.0;
      }
      else {
        dy_first = ((x[2] - x[0]) * (x[2] - x[0]) *
                        (y_array[index_y * x_size + 1] - y_array[index_y * x_size + 0]) -
                    (x[1] - x[0]) * (x[1] - x[0]) *
                        (y_array[index_y * x_size + 2] - y_array[index_y * x_size + 0])) /
                   ((x[2] - x[0]) * (x[1] - x[0]) * (x[2] - x[1]));

        ddy_array[index_y * x_size + 0] = -0.5;

        u[0 * y_size + index_y] = (3. / (x[1] - x[0])) *
                                  ((y_array[index_y * x_size + 1] - y_array[index_y * x_size + 0]) /
                                       (x[1] - x[0]) -
                                   dy_first);
      }

      for (index_x = 1; index_x < x_size - 1; index_x++) {
        sig = (x[index_x] - x[index_x - 1]) / (x[index_x + 1] - x[index_x - 1]);

        p[index_y] = sig * ddy_array[index_y * x_size + (index_x - 1)] + 2.0;

        ddy_array[index_y * x_size + index_x] = (sig - 1.0) / p[index_y];

        u[index_x * y_size + index_y] =
            (y_array[index_y * x_size + (index_x + 1)] - y_array[index_y * x_size + index_x]) /
                (x[index_x + 1] - x[index_x]) -
            (y_array[index_y * x_size + index_x] - y_array[index_y * x_size + (index_x - 1)]) /
                (x[index_x] - x[index_x - 1]);

        u[index_x * y_size + index_y] = (6.0 * u[index_x * y_size + index_y] /
                                             (x[index_x + 1] - x[index_x - 1]) -
                                         sig * u[(index_x - 1) * y_size + index_y]) /
                                        p[index_y];
      }

      if (spline_mode == _SPLINE_NATURAL_) {
        qn[index_y] = un[index_y] = 0.0;
      }
      else {
        dy_last = ((x[x_size - 3] - x[x_size - 1]) * (x[x_size - 3] - x[x_size - 1]) *
                       (y_array[index_y * x_size + (x_size - 2)] -
                        y_array[index_y * x_size + (x_size - 1)]) -
                   (x[x_size - 2] - x[x_size - 1]) * (x[x_size - 2] - x[x_size - 1]) *
                       (y_array[index_y * x_size + (x_size - 3)] -
                        y_array[index_y * x_size + (x_size - 1)])) /
                  ((x[x_size - 3] - x[x_size - 1]) * (x[x_size - 2] - x[x_size - 1]) *
                   (x[x_size - 3] - x[x_size - 2]));

        qn[index_y] = 0.5;

        un[index_y] = (3. / (x[x_size - 1] - x[x_size - 2])) *
                      (dy_last - (y_array[index_y * x_size + (x_size - 1)] -
                                  y_array[index_y * x_size + (x_size - 2)]) /
                                     (x[x_size - 1] - x[x_size - 2]));
      }

      index_x = x_size - 1;

      ddy_array[index_y * x_size + index_x] =
          (un[index_y] - qn[index_y] * u[(index_x - 1) * y_size + index_y]) /
          (qn[index_y] * ddy_array[index_y * x_size + (index_x - 1)] + 1.0);

      for (index_x = x_size - 2; index_x >= 0; index_x--) {
        ddy_array[index_y * x_size + index_x] = ddy_array[index_y * x_size + index_x] *
                                                    ddy_array[index_y * x_size + (index_x + 1)] +
                                                u[index_x * y_size + index_y];
      }
    }
  }
}

void array_integrate_all_spline(const double* array,
                                int n_columns,
                                int n_lines,
                                int index_x, /** from 0 to (n_columns-1) */
                                int index_y,
                                int index_ddy,
                                double* result) {
  int i;
  double h;

  *result = 0;

  for (i = 0; i < n_lines - 1; i++) {
    h = (array[(i + 1) * n_columns + index_x] - array[i * n_columns + index_x]);

    *result += (array[i * n_columns + index_y] + array[(i + 1) * n_columns + index_y]) * h / 2. +
               (array[i * n_columns + index_ddy] + array[(i + 1) * n_columns + index_ddy]) * h * h *
                   h / 24.;
  }
}

void array_integrate_all_trapzd_or_spline(const double* array,
                                          int n_columns,
                                          int n_lines,
                                          int index_start_spline,
                                          int index_x, /** from 0 to (n_columns-1) */
                                          int index_y,
                                          int index_ddy,
                                          double* result) {
  int i;
  double h;

  class_test((index_start_spline < 0) || (index_start_spline >= n_lines),
             "index_start_spline outside of range");

  *result = 0;

  /* trapezoidal integration till given index */

  for (i = 0; i < index_start_spline; i++) {
    h = (array[(i + 1) * n_columns + index_x] - array[i * n_columns + index_x]);

    *result += (array[i * n_columns + index_y] + array[(i + 1) * n_columns + index_y]) * h / 2.;
  }

  /* then, spline integration */

  for (i = index_start_spline; i < n_lines - 1; i++) {
    h = (array[(i + 1) * n_columns + index_x] - array[i * n_columns + index_x]);

    *result += (array[i * n_columns + index_y] + array[(i + 1) * n_columns + index_y]) * h / 2. +
               (array[i * n_columns + index_ddy] + array[(i + 1) * n_columns + index_ddy]) * h * h *
                   h / 24.;
  }
}

/**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
void array_interpolate(const double* array,
                       int n_columns,
                       int n_lines,
                       int index_x, /** from 0 to (n_columns-1) */
                       double x,
                       int* last_index,
                       double* result,
                       int result_size) /** from 1 to n_columns */
{
  int inf, sup, mid, i;
  double weight;

  inf = 0;
  sup = n_lines - 1;

  if (*(array + inf * n_columns + index_x) < *(array + sup * n_columns + index_x)) {
    class_test(x < *(array + inf * n_columns + index_x),
               "x=%e < x_min=%e",
               x,
               *(array + inf * n_columns + index_x));

    class_test(x > *(array + sup * n_columns + index_x),
               "x=%e > x_max=%e",
               x,
               *(array + sup * n_columns + index_x));

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x < *(array + mid * n_columns + index_x)) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(x < *(array + sup * n_columns + index_x),
               "x=%e < x_min=%e",
               x,
               *(array + sup * n_columns + index_x));

    class_test(x > *(array + inf * n_columns + index_x),
               "x=%e > x_max=%e",
               x,
               *(array + inf * n_columns + index_x));

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x > *(array + mid * n_columns + index_x)) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  *last_index = inf;

  weight = (x - *(array + inf * n_columns + index_x)) /
           (*(array + sup * n_columns + index_x) - *(array + inf * n_columns + index_x));

  for (i = 0; i < result_size; i++)
    *(result + i) = *(array + inf * n_columns + i) * (1. - weight) +
                    weight * *(array + sup * n_columns + i);

  *(result + index_x) = x;
}

/**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
void array_interpolate_spline(const double* x_array,
                              int n_lines,
                              const double* array,
                              const double* array_splined,
                              int n_columns,
                              double x,
                              int* last_index,
                              double* result,
                              int result_size) /** from 1 to n_columns */
{
  int inf, sup, mid, i;
  double h, a, b;

  inf = 0;
  sup = n_lines - 1;

  if (x_array[inf] < x_array[sup]) {
    class_test(x < x_array[inf], "x=%e < x_min=%e", x, x_array[inf]);

    class_test(x > x_array[sup], "x=%e > x_max=%e", x, x_array[sup]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x < x_array[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(x < x_array[sup], "x=%e < x_min=%e", x, x_array[sup]);

    class_test(x > x_array[inf], "x=%e > x_max=%e", x, x_array[inf]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x > x_array[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x - x_array[inf]) / h;
  a = 1 - b;

  for (i = 0; i < result_size; i++)
    *(result + i) = a * *(array + inf * n_columns + i) + b * *(array + sup * n_columns + i) +
                    ((a * a * a - a) * *(array_splined + inf * n_columns + i) +
                     (b * b * b - b) * *(array_splined + sup * n_columns + i)) *
                        h * h / 6.;
}

/**
  * Get the y[i] for which y[i]>c
  *
  * Called by nonlinear_HMcode()
  */
void array_search_bisect(int n_lines, const double* array, double c, int* last_index) {
  int inf, sup, mid;

  inf = 0;
  sup = n_lines - 1;

  if (array[inf] < array[sup]) {
    class_test(c < array[inf], "c=%e < y_min=%e", c, array[inf]);

    class_test(c > array[sup], "c=%e > y_max=%e", c, array[sup]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (c < array[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(c < array[sup], "x=%e < x_min=%e", c, array[sup]);

    class_test(c > array[inf], "x=%e > x_max=%e", c, array[inf]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (c > array[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  *last_index = inf;
}

/**
  * interpolate to get y_i(x), when x and y_i are in different arrays
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
void array_interpolate_linear(const double* x_array,
                              int n_lines,
                              const double* array,
                              int n_columns,
                              double x,
                              int* last_index,
                              double* result,
                              int result_size) /** from 1 to n_columns */
{
  int inf, sup, mid, i;
  double h, a, b;

  inf = 0;
  sup = n_lines - 1;

  if (x_array[inf] < x_array[sup]) {
    class_test(x < x_array[inf], "x=%e < x_min=%e", x, x_array[inf]);

    class_test(x > x_array[sup], "x=%e > x_max=%e", x, x_array[sup]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x < x_array[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(x < x_array[sup], "x=%e < x_min=%e", x, x_array[sup]);

    class_test(x > x_array[inf], "x=%e > x_max=%e", x, x_array[inf]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x > x_array[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x - x_array[inf]) / h;
  a = 1 - b;

  for (i = 0; i < result_size; i++)
    *(result + i) = a * *(array + inf * n_columns + i) + b * *(array + sup * n_columns + i);
}

/**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
void array_interpolate_growing_closeby(const double* array,
                                       int n_columns,
                                       int n_lines,
                                       int index_x, /** from 0 to (n_columns-1) */
                                       double x,
                                       int* last_index,
                                       double* result,
                                       int result_size) /** from 1 to n_columns */
{
  int inf, sup, i;
  double weight;

  inf = *last_index;

  while (x < *(array + inf * n_columns + index_x)) {
    inf--;
    class_test(inf < 0, "x=%e < x_min=%e", x, array[index_x]);
  }
  sup = inf + 1;
  while (x > *(array + sup * n_columns + index_x)) {
    sup++;
    class_test(sup > (n_lines - 1),
               "x=%e > x_max=%e",
               x,
               array[(n_lines - 1) * n_columns + index_x]);
  }
  inf = sup - 1;

  *last_index = inf;

  weight = (x - *(array + inf * n_columns + index_x)) /
           (*(array + sup * n_columns + index_x) - *(array + inf * n_columns + index_x));

  for (i = 0; i < result_size; i++)
    *(result + i) = *(array + inf * n_columns + i) * (1. - weight) +
                    weight * *(array + sup * n_columns + i);

  *(result + index_x) = x;
}

/**
  * interpolate to get y(x), when x and y are two columns of the same array, x is arranged in growing order, and the point x is presumably close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
void array_interpolate_one_growing_closeby(const double* array,
                                           int n_columns,
                                           int n_lines,
                                           int index_x, /** from 0 to (n_columns-1) */
                                           double x,
                                           int* last_index,
                                           int index_y,
                                           double* result) {
  int inf, sup;
  double weight;

  inf = *last_index;

  while (x < *(array + inf * n_columns + index_x)) {
    inf--;
    class_test(inf < 0, "x=%e < x_min=%e", x, array[index_x]);
  }
  sup = inf + 1;
  while (x > *(array + sup * n_columns + index_x)) {
    sup++;
    class_test(sup > (n_lines - 1),
               "x=%e > x_max=%e",
               x,
               array[(n_lines - 1) * n_columns + index_x]);
  }
  inf = sup - 1;

  *last_index = inf;

  weight = (x - *(array + inf * n_columns + index_x)) /
           (*(array + sup * n_columns + index_x) - *(array + inf * n_columns + index_x));

  *result = *(array + inf * n_columns + index_y) * (1. - weight) +
            *(array + sup * n_columns + index_y) * weight;
}

/**
  * interpolate to get y_i(x), when x and y_i are all columns of the same array, x is arranged in growing order, and the point x is presumably very close to the previous point x from the last call of this function.
  *
  * Called by background_at_eta(); background_eta_of_z(); background_solve(); thermodynamics_at_z().
  */
void array_interpolate_spline_growing_closeby(const double* x_array,
                                              int n_lines,
                                              const double* array,
                                              const double* array_splined,
                                              int n_columns,
                                              double x,
                                              int* last_index,
                                              double* result,
                                              int result_size) /** from 1 to n_columns */
{
  int inf, sup, i;
  double h, a, b;

  /*
  if (*last_index < 0) {
    class_stop( "problem with last_index =%d < 0", *last_index);
  }
  if (*last_index > (n_lines-1)) {
    class_stop( "problem with last_index =%d > %d", *last_index,n_lines-1);
  }
  */

  inf = *last_index;
  class_test(inf < 0 || inf > (n_lines - 1),
             "*lastindex=%d out of range [0:%d]\n",
             inf,
             n_lines - 1);
  while (x < x_array[inf]) {
    inf--;
    class_test(inf < 0, "x=%e < x_min=%e", x, x_array[0]);
  }
  sup = inf + 1;
  while (x > x_array[sup]) {
    sup++;
    class_test(sup > (n_lines - 1), "x=%e > x_max=%e", x, x_array[n_lines - 1]);
  }
  inf = sup - 1;

  *last_index = inf;

  h = x_array[sup] - x_array[inf];
  b = (x - x_array[inf]) / h;
  a = 1 - b;

  for (i = 0; i < result_size; i++)
    *(result + i) = a * *(array + inf * n_columns + i) + b * *(array + sup * n_columns + i) +
                    ((a * a * a - a) * *(array_splined + inf * n_columns + i) +
                     (b * b * b - b) * *(array_splined + sup * n_columns + i)) *
                        h * h / 6.;
}

/**
 * interpolate linearily to get y_i(x), when x and y_i are in two different arrays
 *
 * Called by transfer_interpolate_sources(); transfer_functions_at_k(); perturb_sources_at_eta().
 */
void array_interpolate_two(const double* array_x,
                           int n_columns_x,
                           int index_x, /** from 0 to (n_columns_x-1) */
                           const double* array_y,
                           int n_columns_y,
                           int n_lines, /** must be the same for array_x and array_y */
                           double x,
                           double* result,
                           int result_size) /** from 1 to n_columns_y */
{
  int inf, sup, mid, i;
  double weight;

  inf = 0;
  sup = n_lines - 1;

  if (array_x[inf * n_columns_x + index_x] < array_x[sup * n_columns_x + index_x]) {
    class_test(x < array_x[inf * n_columns_x + index_x],
               "x=%e < x_min=%e",
               x,
               array_x[inf * n_columns_x + index_x]);

    class_test(x > array_x[sup * n_columns_x + index_x],
               "x=%e > x_max=%e",
               x,
               array_x[sup * n_columns_x + index_x]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x < array_x[mid * n_columns_x + index_x]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(x < *(array_x + sup * n_columns_x + index_x),
               "x=%e < x_min=%e",
               x,
               *(array_x + sup * n_columns_x + index_x));

    class_test(x > *(array_x + inf * n_columns_x + index_x),
               "x=%e > x_max=%e",
               x,
               *(array_x + inf * n_columns_x + index_x));

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x > *(array_x + mid * n_columns_x + index_x)) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  weight = (x - *(array_x + inf * n_columns_x + index_x)) /
           (*(array_x + sup * n_columns_x + index_x) - *(array_x + inf * n_columns_x + index_x));

  for (i = 0; i < result_size; i++)
    *(result + i) = *(array_y + i * n_lines + inf) * (1. - weight) +
                    weight * *(array_y + i * n_lines + sup);
}

/**
 * Same as array_interpolate_two, but with order of indices exchanged in array_y
 */
void array_interpolate_two_bis(const double* array_x,
                               int n_columns_x,
                               int index_x, /** from 0 to (n_columns_x-1) */
                               const double* array_y,
                               int n_columns_y,
                               int n_lines, /** must be the same for array_x and array_y */
                               double x,
                               double* result,
                               int result_size) /** from 1 to n_columns_y */
{
  int inf, sup, mid, i;
  double weight;

  inf = 0;
  sup = n_lines - 1;

  if (array_x[inf * n_columns_x + index_x] < array_x[sup * n_columns_x + index_x]) {
    class_test(x < array_x[inf * n_columns_x + index_x],
               "x=%e < x_min=%e",
               x,
               array_x[inf * n_columns_x + index_x]);

    class_test(x > array_x[sup * n_columns_x + index_x],
               "x=%e > x_max=%e",
               x,
               array_x[sup * n_columns_x + index_x]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x < array_x[mid * n_columns_x + index_x]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(x < *(array_x + sup * n_columns_x + index_x),
               "x=%e < x_min=%e",
               x,
               *(array_x + sup * n_columns_x + index_x));

    class_test(x > *(array_x + inf * n_columns_x + index_x),
               "x=%e > x_max=%e",
               x,
               *(array_x + inf * n_columns_x + index_x));

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x > *(array_x + mid * n_columns_x + index_x)) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  weight = (x - *(array_x + inf * n_columns_x + index_x)) /
           (*(array_x + sup * n_columns_x + index_x) - *(array_x + inf * n_columns_x + index_x));

  for (i = 0; i < result_size; i++)
    *(result + i) = *(array_y + inf * n_columns_y + i) * (1. - weight) +
                    weight * *(array_y + sup * n_columns_y + i);
}

/**
 * interpolate linearily to get y_i(x), when x and y_i are in two different arrays
 *
 * Called by transfer_interpolate_sources(); transfer_functions_at_k(); perturb_sources_at_eta().
 */
void array_interpolate_two_arrays_one_column(
    const double* array_x, /* assumed to be a vector (i.e. one column array) */
    const double* array_y,
    int n_columns_y,
    int index_y, /* between 0 and (n_columns_y-1) */
    int n_lines, /** must be the same for array_x and array_y */
    double x,
    double* result) {
  int inf, sup, mid;
  double weight;
  double epsilon = 1e-9;

  inf = 0;
  sup = n_lines - 1;

  if (array_x[inf] < array_x[sup]) {
    class_test(x < array_x[inf] - epsilon, "x=%e < x_min=%e", x, array_x[inf]);

    class_test(x > array_x[sup] + epsilon, "x=%e > x_max=%e", x, array_x[sup]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x < array_x[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  else {
    class_test(x < array_x[sup] - epsilon, "x=%e < x_min=%e", x, array_x[sup]);

    class_test(x > array_x[inf] + epsilon, "x=%e > x_max=%e", x, array_x[inf]);

    while (sup - inf > 1) {
      mid = (int) (0.5 * (inf + sup));
      if (x > array_x[mid]) {
        sup = mid;
      }
      else {
        inf = mid;
      }
    }
  }

  weight = (x - array_x[inf]) / (array_x[sup] - array_x[inf]);

  *result = array_y[index_y * n_lines + inf] * (1. - weight) +
            weight * array_y[index_y * n_lines + sup];
}

/**
 * cubic interpolation of array with equally space abscisses
 */

void array_interpolate_cubic_equal(
    double x0, double dx, const double* yarray, int Nx, double x, double* result) {
  int i;
  double frac;

  class_test((dx > 0 && (x < x0 || x > x0 + dx * (Nx - 1))),
             "x=%e out of range [%e %e]",
             x,
             x0,
             x0 + dx * (Nx - 1));

  class_test((dx < 0 && (x > x0 || x < x0 + dx * (Nx - 1))),
             "x=%e out of range [%e %e]",
             x,
             x0 + dx * (Nx - 1),
             x0);

  i = (int) floor((x - x0) / dx);
  if (i < 1)
    i = 1;
  if (i > Nx - 3)
    i = Nx - 3;
  frac    = (x - x0) / dx - i;
  yarray += i - 1;

  *result = -yarray[0] * frac * (1. - frac) * (2. - frac) / 6. +
            yarray[1] * (1. + frac) * (1. - frac) * (2. - frac) / 2. +
            yarray[2] * (1. + frac) * frac * (2. - frac) / 2. +
            yarray[3] * (1. + frac) * frac * (frac - 1.) / 6.;
}

void array_interpolate_parabola(double x1,
                                double x2,
                                double x3,
                                double x,
                                double y1,
                                double y2,
                                double y3,
                                double* y,
                                double* dy,
                                double* ddy) {
  double a, b, c;

  /*
    a x_i**2 + b x_i + c = y_i

    a (x1**2-x2**2) + b (x1-x2) = y1-y2
    a (x3**2-x2**2) + b (x3-x2) = y3-y2

    a (x1**2-x2**2)(x3**2-x2**2) + b (x1-x2)(x3**2-x2**2) = (y1-y2)(x3**2-x2**2)
    a (x3**2-x2**2)(x1**2-x2**2) + b (x3-x2)(x1**2-x2**2) = (y3-y2)(x1**2-x2**2)

    b = [(y1-y2)(x3**2-x2**2) - (y3-y2)(x1**2-x2**2)]/(x1-x2)(x3-x2)(x3-x1)

  */

  b = ((y1 - y2) * (x3 - x2) * (x3 + x2) - (y3 - y2) * (x1 - x2) * (x1 + x2)) / (x1 - x2) /
      (x3 - x2) / (x3 - x1);

  a = (y1 - y2 - b * (x1 - x2)) / (x1 - x2) / (x1 + x2);

  c = y2 - b * x2 - a * x2 * x2;

  *y   = a * x * x + b * x + c;
  *dy  = 2. * a * x + b;
  *ddy = 2. * a;
}

void array_smooth(double* array,
                  int n_columns,
                  int n_lines,
                  int index, /** from 0 to (n_columns-1) */
                  int radius) {
  std::vector<double> vec_smooth(n_lines);
  double* smooth = vec_smooth.data();
  int i, j, jmin, jmax;
  double weigth;

  for (i = 0; i < n_lines; i++) {
    smooth[i] = 0.;
    weigth    = 0.;
    jmin      = std::max(i - radius, 0);
    jmax      = std::min(i + radius, n_lines - 1);
    for (j = jmin; j <= jmax; j++) {
      smooth[i] += array[j * n_columns + index];
      weigth    += 1.;
    }
    smooth[i] /= weigth;
  }

  for (i = 0; i < n_lines; i++)
    array[i * n_columns + index] = smooth[i];
}

/**
 * Compute quadrature weights for the trapezoidal integration method, xhen x is in gorwing order.
 *
 * @param x                     Input: Grid points on which f() is known.
 * @param n                     Input: number of grid points.
 * @param w_trapz               Output: Weights of the trapezoidal method.
 * @return the error status
 */

/**
 * Compute quadrature weights for the trapezoidal integration method, when x is in decreasing order.
 *
 * @param x                     Input: Grid points on which f() is known.
 * @param n                     Input: number of grid points.
 * @param w_trapz               Output: Weights of the trapezoidal method.
 * @return the error status
 */

void array_trapezoidal_mweights(const double* x, int n, double* w_trapz) {
  int i;

  /* Case with just one point. */
  if (n == 1) {
    w_trapz[0] = 1.0;
  }
  else if (n > 1) {
    //Set edgeweights:
    w_trapz[0]     = 0.5 * (x[0] - x[1]);
    w_trapz[n - 1] = 0.5 * (x[n - 2] - x[n - 1]);
    //Set inner weights:
    for (i = 1; i < (n - 1); i++) {
      w_trapz[i] = 0.5 * (x[i - 1] - x[i + 1]);
    }
  }
}

/**
 * Compute integral of function using trapezoidal method.
 *
 * @param integrand             Input: The function we are integrating.
 * @param n                     Input: Compute integral on grid [0;n-1].
 * @param w_trapz               Input: Weights of the trapezoidal method.
 * @param I                     Output: The integral.
 * @return the error status
 */

void array_trapezoidal_integral(const double* integrand, int n, const double* w_trapz, double* I) {
  int i;
  double res = 0.0;
  for (i = 0; i < n; i++) {
    res += integrand[i] * w_trapz[i];
  }
  *I = res;
}

/**
 * Compute convolution integral of product of two functions using trapezoidal method.
 *
 * @param integrand1            Input: Function 1.
 * @param integrand2            Input: Function 2.
 * @param n                     Input: Compute integral on grid [0;n-1].
 * @param w_trapz               Input: Weights of the trapezoidal method.
 * @param I                     Output: The integral.
 * @return the error status
 */

void array_trapezoidal_convolution(
    const double* integrand1, const double* integrand2, int n, const double* w_trapz, double* I) {
  int i;
  double res = 0.0;
  for (i = 0; i < n; i++) {
    res += integrand1[i] * integrand2[i] * w_trapz[i];
  }
  *I = res;
}
