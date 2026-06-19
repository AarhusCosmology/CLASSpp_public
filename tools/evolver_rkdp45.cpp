#include "evolver_rkdp45.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

/**
 * Dormand-Prince 4(5) explicit adaptive Runge-Kutta integrator with 4th-order
 * dense output. Ported from branch PhD2024-EBH. Conforms to the canonical CLASS
 * evolver signature (shared with evolver_ndf15 / evolver_rk).
 *
 * minimum_variation, evaluate_timescale and timestep_over_timescale are part of
 * the shared signature but unused here: this solver performs its own embedded
 * error control.
 */
int evolver_rkdp45(
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
  (void) minimum_variation;
  (void) evaluate_timescale;
  (void) timestep_over_timescale;

  /* Contract (shared with evolver_rk): x_sampling is a non-null array of
     output points, and tolerance > 0 (threshold below divides by it). */
  class_test(x_sampling == nullptr, "rkdp45 requires a non-null x_sampling array");
  class_test(tolerance <= 0., "rkdp45 requires tolerance > 0 (got %e)", tolerance);

  const int neq          = y_size;
  const double rtol      = tolerance;
  const double abstol    = 1e-15; /* matches ndf15 */
  const double threshold = abstol / rtol;
  const double pow_grow  = 0.2;
  const int s            = 7; /* DP45 stages */

  /* Butcher tableau (Dormand-Prince 4(5)) */
  const double ci[s] = {0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0};
  const double bi[s] =
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
  const double bi_diff[s] = {71.0 / 57600.0,
                             0.0,
                             -71.0 / 16695.0,
                             71.0 / 1920.0,
                             -17253.0 / 339200.0,
                             22.0 / 525.0,
                             -1.0 / 40.0};
  /* zero-init is load-bearing: only the strict lower triangle is set below;
     the unused upper triangle/diagonal must read as 0. */
  double ai[s][s] = {{0.0}};
  ai[1][0]        = 0.2;
  ai[2][0]        = 3.0 / 40.0;
  ai[2][1]        = 9.0 / 40.0;
  ai[3][0]        = 44.0 / 45.0;
  ai[3][1]        = -56.0 / 15.0;
  ai[3][2]        = 32.0 / 9.0;
  ai[4][0]        = 19372.0 / 6561.0;
  ai[4][1]        = -25360.0 / 2187.0;
  ai[4][2]        = 64448.0 / 6561.0;
  ai[4][3]        = -212.0 / 729.0;
  ai[5][0]        = 9017.0 / 3168.0;
  ai[5][1]        = -355.0 / 33.0;
  ai[5][2]        = 46732.0 / 5247.0;
  ai[5][3]        = 49.0 / 176.0;
  ai[5][4]        = -5103.0 / 18656.0;
  ai[6][0]        = 35.0 / 384.0;
  ai[6][1]        = 0.0;
  ai[6][2]        = 500.0 / 1113.0;
  ai[6][3]        = 125.0 / 192.0;
  ai[6][4]        = -2187.0 / 6784.0;
  ai[6][5]        = 11.0 / 84.0;

  /* 4th-order dense-output interpolation coefficients */
  const double ixx[5][3] = {{1500.0 / 371.0, -1000.0 / 159.0, 1000.0 / 371.0},
                            {-125.0 / 32.0, 125.0 / 12.0, -375.0 / 64.0},
                            {9477.0 / 3392.0, -729.0 / 106.0, 25515.0 / 6784.0},
                            {-11.0 / 7.0, 11.0 / 3.0, -55.0 / 28.0},
                            {1.5, -4.0, 2.5}};
  const double i01 = -183.0 / 64.0, i02 = 37.0 / 12.0, i03 = -145.0 / 128.0;

  std::vector<double> ynew(neq), ytemp(neq), err(neq), yinterp(neq), dyinterp(neq);
  std::vector<double> ki(s * neq);
  double bi_vec_y[s], bi_vec_dy[s];

  double t = x_ini;

  /* initialise k0 = f(t, y) */
  (*derivs)(t, y, ki.data(), parameters_and_workspace_for_derivs);

  const double hmax = fabs(x_end - x_ini) / 10.0;
  double absh;
  if (x_size > 1)
    absh = std::min(hmax, fabs(x_sampling[1] - x_sampling[0]));
  else
    absh = hmax;
  if (absh == 0.0)
    absh = hmax;

  /* initial step from derivative scale */
  double rh = 0.0;
  for (int k = 0; k < neq; k++) {
    double maxtmp = std::max(fabs(y[k]), threshold);
    rh            = std::max(rh, fabs(ki[k]) / maxtmp);
  }
  rh /= 0.8 * pow(rtol, pow_grow);
  if (absh * rh > 1.0)
    absh = 1.0 / rh;

  const int tdir = (x_end > x_ini) ? 1 : -1;
  double hnew    = absh * tdir;
  int nofailed   = _TRUE_;

  int idx = 0;
  if ((t - x_end) * tdir < 0.0)
    for (idx = 0; (idx < x_size) && ((x_sampling[idx] - t) * tdir < 0.0); idx++)
      ;

  while ((t - x_end) * tdir < 0.0) {
    double h          = hnew;
    const double hmin = 100.0 * DBL_MIN * fabs(t);
    class_test(fabs(h) < hmin,
               "rkdp45: step size %e fell below minimum %e at x=%e",
               fabs(h),
               hmin,
               t);
    if (fabs(h) > 0.9 * fabs(x_end - t))
      h = x_end - t;

    /* stage 0 (FSAL: ki[0..neq) already holds k0) */
    for (int k = 0; k < neq; k++) {
      ynew[k] = y[k] + h * bi[0] * ki[k];
      err[k]  = h * bi_diff[0] * ki[k];
    }
    for (int i = 1; i < s; i++) {
      for (int k = 0; k < neq; k++)
        ytemp[k] = y[k];
      for (int j = 0; j < i; j++)
        for (int k = 0; k < neq; k++)
          ytemp[k] += h * ai[i][j] * ki[j * neq + k];
      (*derivs)(t + ci[i] * h,
                ytemp.data(),
                ki.data() + i * neq,
                parameters_and_workspace_for_derivs);
      for (int k = 0; k < neq; k++) {
        ynew[k] += h * bi[i] * ki[i * neq + k];
        err[k]  += h * bi_diff[i] * ki[i * neq + k];
      }
    }

    double errmax = 0.0;
    for (int k = 0; k < neq; k++) {
      double errtemp = fabs(err[k] / std::max(threshold, fabs(ynew[k])));
      if (errtemp > errmax)
        errmax = errtemp;
    }

    if (errmax > rtol) {
      /* step rejected. Standard ode45 control: first failure shrinks by the
       * error-proportional factor and marks nofailed; consecutive failures
       * halve. (This corrects a transcription bug in the PhD2024-EBH branch
       * where nofailed was never set, making the halving path dead code.) */
      if (nofailed == _TRUE_) {
        nofailed = _FALSE_;
        hnew = tdir * std::max(hmin, fabs(h) * std::max(0.1, 0.8 * pow(rtol / errmax, pow_grow)));
      }
      else {
        hnew = tdir * std::max(hmin, 0.5 * fabs(h));
      }
      continue;
    }

    /* step accepted */
    if (print_variables != nullptr)
      (*print_variables)(t + h,
                         ynew.data(),
                         ki.data() + 6 * neq,
                         parameters_and_workspace_for_derivs);

    nofailed = _TRUE_;
    hnew     = tdir * std::max(hmin, fabs(h) * std::max(0.1, 0.8 * pow(rtol / errmax, pow_grow)));
    const double tnew = t + h;

    /* emit output at all sampling points within (t, tnew] */
    for (; (idx < x_size) && ((tnew - x_sampling[idx]) * tdir >= 0.0); idx++) {
      if (tnew == x_sampling[idx]) {
        (*output)(tnew, ynew.data(), ki.data() + 6 * neq, idx, parameters_and_workspace_for_derivs);
      }
      else {
        const double ti  = x_sampling[idx];
        const double ss1 = (ti - t) / h, ss2 = ss1 * ss1, ss3 = ss2 * ss1, ss4 = ss2 * ss2;
        bi_vec_y[0]  = ss1 + i01 * ss2 + i02 * ss3 + i03 * ss4;
        bi_vec_dy[0] = 1.0 + i01 * 2.0 * ss1 + i02 * 3.0 * ss2 + i03 * 4.0 * ss3;
        for (int i = 2; i < 7; i++) {
          bi_vec_y[i]  = ixx[i - 2][0] * ss2 + ixx[i - 2][1] * ss3 + ixx[i - 2][2] * ss4;
          bi_vec_dy[i] = ixx[i - 2][0] * 2 * ss1 + ixx[i - 2][1] * 3 * ss2 +
                         ixx[i - 2][2] * 4 * ss3;
        }
        for (int k = 0; k < neq; k++) {
          if (used_in_output[k] == _TRUE_) {
            yinterp[k]  = y[k];
            dyinterp[k] = 0.0;
            for (int i = 0; i < 7; i++) {
              if (i != 1) {
                yinterp[k]  += h * bi_vec_y[i] * ki[i * neq + k];
                dyinterp[k] += bi_vec_dy[i] * ki[i * neq + k];
              }
            }
          }
        }
        (*output)(ti, yinterp.data(), dyinterp.data(), idx, parameters_and_workspace_for_derivs);
      }
    }

    for (int k = 0; k < neq; k++) {
      y[k]  = ynew[k];
      ki[k] = ki[6 * neq + k]; /* FSAL: last stage becomes next k0 */
    }
    t = tnew;
  }

  return _SUCCESS_;
}
