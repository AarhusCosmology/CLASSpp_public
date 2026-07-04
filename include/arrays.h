/**
 * definitions for module thermodynamics.c
 */

#ifndef __ARRAYS__
#define __ARRAYS__

#include <vector>

#include "common.h"

#ifdef _WIN32
#define _restrict
#endif

#define _SPLINE_NATURAL_ 0   /**< natural spline: ddy0=ddyn=0 */
#define _SPLINE_EST_DERIV_ 1 /**< spline with estimation of first derivative on both edges */

void array_derive_spline(const double* x_array,
                         int n_lines,
                         double* array,
                         const double* array_splined,
                         int n_columns,
                         int index_y,
                         int index_dydx);

void array_derive_spline_table_line_to_line(const double* x_array,
                                            int n_lines,
                                            double* array,
                                            int n_columns,
                                            int index_y,
                                            int index_ddy,
                                            int index_dy);

void array_spline(double* array,
                  int n_columns,
                  int n_lines,
                  int index_x, /** from 0 to (n_columns-1) */
                  int index_y,
                  int index_ddydx2,
                  short spline_mode);

void array_spline_table_line_to_line(const double* x, /* vector of size x_size */
                                     int x_size,
                                     double* array,
                                     int n_columns,
                                     int index_y,
                                     int index_ddydx2,
                                     short spline_mode);

void array_spline_table_columns(const double* x,
                                int x_size,
                                const double* y_array,
                                int y_size,
                                double* ddy_array,
                                short spline_mode);

void array_spline_table_columns2(const double* x,
                                 int x_size,
                                 const double* y_array,
                                 int y_size,
                                 double* ddy_array,
                                 short spline_mode);

void array_spline_table_lines(const double* x,
                              int x_size,
                              const double* y_array,
                              int y_size,
                              double* ddy_array,
                              short spline_mode);

void array_integrate_all_spline(const double* array,
                                int n_columns,
                                int n_lines,
                                int index_x,
                                int index_y,
                                int index_ddy,
                                double* result);

void array_integrate_all_trapzd_or_spline(const double* array,
                                          int n_columns,
                                          int n_lines,
                                          int index_start_spline,
                                          int index_x, /** from 0 to (n_columns-1) */
                                          int index_y,
                                          int index_ddy,
                                          double* result);

void array_integrate_spline_table_line_to_line(const double* x_array,
                                               int n_lines,
                                               double* array,
                                               int n_columns,
                                               int index_y,
                                               int index_ddy,
                                               int index_inty);

void array_interpolate_spline(const double* x_array,
                              int n_lines,
                              const double* array,
                              const double* array_splined,
                              int n_columns,
                              double x,
                              int* last_index,
                              double* result,
                              int result_size); /** from 1 to n_columns */

void array_search_bisect(int n_lines, const double* array, double c, int* last_index);

void array_interpolate_linear(const double* x_array,
                              int n_lines,
                              const double* array,
                              int n_columns,
                              double x,
                              int* last_index,
                              double* result,
                              int result_size); /** from 1 to n_columns */

void array_interpolate_one_growing_closeby(const double* array,
                                           int n_columns,
                                           int n_lines,
                                           int index_x, /** from 0 to (n_columns-1) */
                                           double x,
                                           int* last_index,
                                           int index_y,
                                           double* result);

void array_interpolate_spline_growing_closeby(const double* x_array,
                                              int n_lines,
                                              const double* array,
                                              const double* array_splined,
                                              int n_columns,
                                              double x,
                                              int* last_index,
                                              double* result,
                                              int result_size); /** from 1 to n_columns */

void array_interpolate_two(const double* array_x,
                           int n_columns_x,
                           int index_x, /** from 0 to (n_columns_x-1) */
                           const double* array_y,
                           int n_lines, /** must be the same for array_x and array_y */
                           double x,
                           double* result,
                           int result_size);

void array_interpolate_two_bis(const double* array_x,
                               int n_columns_x,
                               int index_x, /** from 0 to (n_columns_x-1) */
                               const double* array_y,
                               int n_columns_y,
                               int n_lines, /** must be the same for array_x and array_y */
                               double x,
                               double* result,
                               int result_size); /** from 1 to n_columns_y */

void array_interpolate_two_arrays_one_column(
    const double* array_x, /* assumed to be a vector (i.e. one column array) */
    const double* array_y,
    int index_y,
    int n_lines, /** must be the same for array_x and array_y */
    double x,
    double* result);

void array_interpolate_cubic_equal(
    double x0, double dx, const double* yarray, int Nx, double x, double* result);

void array_interpolate_parabola(double x1,
                                double x2,
                                double x3,
                                double x,
                                double y1,
                                double y2,
                                double y3,
                                double* y,
                                double* dy,
                                double* ddy);

void array_smooth(double* array,
                  int n_columns,
                  int n_lines,
                  int index, /** from 0 to (n_columns-1) */
                  int radius);

void array_trapezoidal_mweights(const double* x, int n, double* w_trapz);

void array_trapezoidal_integral(const double* integrand, int n, const double* w_trapz, double* I);

void array_trapezoidal_convolution(
    const double* integrand1, const double* integrand2, int n, const double* w_trapz, double* I);

#endif
