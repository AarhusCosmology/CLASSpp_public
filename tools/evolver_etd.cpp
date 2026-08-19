#include "evolver_etd.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include "nonfinite.h"

namespace {

/** phi_k(z) = sum_{j>=0} z^j / (j+k)!, i.e.
 *      phi0 = e^z,  phi1 = (e^z-1)/z,  phi2 = (e^z-1-z)/z^2,
 *      phi3 = (e^z-1-z-z^2/2)/z^3,
 *  evaluated together because ETDRK4 needs all of them at the same argument.
 *
 *  The closed forms cancel catastrophically as z -> 0: phi3's numerator is the
 *  difference of quantities agreeing to z^3/6, so it loses ~3 digits per decade of
 *  |z|. Below |z| = 0.5 the series is used instead, truncated after z^8; the first
 *  dropped term for phi3 is z^9/12! = 4e-12 at the switch point, i.e. 2e-11 relative
 *  -- far below any tolerance this evolver runs at. Above 0.5 the direct forms lose
 *  at most ~1.8 digits, leaving ~14.
 *
 *  Large NEGATIVE z is the operating regime (a relaxation) and is benign: e^z -> 0 so
 *  every phi_k -> -1/z or smaller. Large POSITIVE z would overflow exp; it means a
 *  GROWING diagonal, where the step is wrong regardless, so z is clamped and the
 *  error controller is left to reject rather than being handed an inf it cannot
 *  compare (`errtemp > errmax` is false for NaN, which would silently accept).
 */
inline void PhiFuncs(double z, double& e0, double& p1, double& p2, double& p3) {
  constexpr double kSeries = 0.5;
  constexpr double kZMax   = 500.0;  // exp(500) is finite; exp(710) is not
  // Exactly zero is the common case, not a corner: only the interacting species report
  // a diagonal at all, so in a standard CLASS run nearly every component takes this
  // branch on every stage and every dense-output sample.
  if (z == 0.) {
    e0 = 1.;
    p1 = 1.;
    p2 = 0.5;
    p3 = 1.0 / 6.0;
    return;
  }
  if (z > kZMax)
    z = kZMax;
  if (std::fabs(z) < kSeries) {
    // phi_k = sum_j z^j/(j+k)!. The coefficient at term j is 1/(j+k)!, so advancing
    // j divides it by (j+k) -- which is why the three divisors differ.
    p1        = 0.;
    p2        = 0.;
    p3        = 0.;
    double f1 = 1.0, f2 = 0.5, f3 = 1.0 / 6.0, zp = 1.;  // j = 0: 1/1!, 1/2!, 1/3!
    for (int j = 0; j <= 8; ++j) {
      if (j > 0) {
        zp *= z;
        f1 /= (double) (j + 1);
        f2 /= (double) (j + 2);
        f3 /= (double) (j + 3);
      }
      p1 += f1 * zp;
      p2 += f2 * zp;
      p3 += f3 * zp;
    }
    e0 = std::exp(z);
    return;
  }
  e0              = std::exp(z);
  const double e  = std::expm1(z);
  const double z2 = z * z;
  p1              = e / z;
  p2              = (e - z) / z2;
  p3              = (e - z - 0.5 * z2) / (z2 * z);
}

}  // namespace

void evolver_etd(
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
    void (*derivs_diagonal)(double x, double* y, double* diag, void* parameters_and_workspace)) {
  (void) minimum_variation;
  (void) evaluate_timescale;
  (void) timestep_over_timescale;

  class_test(x_sampling == nullptr, "etd requires a non-null x_sampling array");
  class_test(tolerance <= 0., "etd requires tolerance > 0 (got %e)", tolerance);

  const int neq          = y_size;
  const double rtol      = tolerance;
  const double abstol    = 1e-15; /* matches ndf15 / rkdp45 */
  const double threshold = abstol / rtol;

  /* The accepted solution is order 4 and the embedded companion below is order 3, so
     the measured difference is an O(h^4) local quantity and 1/4 is its exponent.

     This is LOCAL EXTRAPOLATION, as in DP45: what is controlled is the error of the
     solution that is thrown away, and what is advanced is the better one. It costs 4
     RHS per attempted step, against ~11 for the step doubling it replaced, and is
     more accurate at the same tolerance -- the median |dP/P| against rkdp45 improved
     3x at Gamma=1e5 and 16x at Gamma=1e6. */
  const double pow_grow = 0.25;

  /* Eight scratch blocks: the stage remainders N1..N4, the stage states U2..U4, tmp. */
  std::vector<double> f0(neq), ynew(neq), diag(neq), err(neq);
  std::vector<double> yinterp(neq), dyinterp(neq), fnew(neq), scratch((size_t) 8 * neq);

  double t = x_ini;
  (*derivs)(t, y, f0.data(), parameters_and_workspace_for_derivs);

  const double hmax = std::fabs(x_end - x_ini) / 10.0;
  double absh = (x_size > 1) ? std::min(hmax, std::fabs(x_sampling[1] - x_sampling[0])) : hmax;
  if (absh == 0.0)
    absh = hmax;

  double rh = 0.0;
  for (int k = 0; k < neq; k++)
    rh = std::max(rh, std::fabs(f0[k]) / std::max(std::fabs(y[k]), threshold));
  rh /= 0.8 * pow(rtol, pow_grow);
  if (absh * rh > 1.0)
    absh = 1.0 / rh;

  const int tdir = (x_end > x_ini) ? 1 : -1;
  double hnew    = absh * tdir;
  bool nofailed  = true;

  int idx = 0;
  if ((t - x_end) * tdir < 0.0)
    for (idx = 0; (idx < x_size) && ((x_sampling[idx] - t) * tdir < 0.0); idx++)
      ;

  while ((t - x_end) * tdir < 0.0) {
    double h          = hnew;
    const double hmin = 100.0 * DBL_MIN * std::fabs(t);
    class_test(std::fabs(h) < hmin,
               "etd: step size %e fell below minimum %e at x=%e",
               std::fabs(h),
               hmin,
               t);
    if (std::fabs(h) > 0.9 * std::fabs(x_end - t))
      h = x_end - t;
    /* NO-PROGRESS GUARD, after the clamp because the clamp is what can produce the
       offending h: hmin above is 100*DBL_MIN*|t|, an absolute floor many orders below
       the SPACING of t, so `t + h == t` can hold while that test passes and the loop
       spins forever at one t. Same failure mode, and same fix, as the ndf15 guard;
       see nonfinite.h for why the test must live in a non-fast-math TU. Cannot fire on
       a legitimate final step: the loop runs only while t != x_end, so h = x_end - t
       moves t by construction. */
    class_test(StepMakesNoProgress(t, h),
               "etd: step size %e underflows the resolution of x=%e (t+h==t), so the "
               "integration cannot advance. The step controller is being driven to zero "
               "by something upstream.",
               h,
               t);

    /* Diagonal at the CURRENT state, re-linearised every attempted step. It is
       frozen only across the step, which is what makes this a Rosenbrock rather
       than a semilinear ETD method. */
    if (derivs_diagonal != nullptr)
      (*derivs_diagonal)(t, y, diag.data(), parameters_and_workspace_for_derivs);
    else
      std::fill(diag.begin(), diag.end(), 0.);

    /* One ETDRK4 step of size hh from state `yin` (whose RHS is `fin`) into `yout`.
       The diagonal is held fixed across the step. The stage remainders N1..N4 are left
       in scratch for the embedded estimator and the continuous extension below. */
    auto etdrk4 = [&](const double* yin, const double* fin, double tin, double hh, double* yout) {
      double* N1  = scratch.data() + 0 * neq;
      double* N2  = scratch.data() + 1 * neq;
      double* N3  = scratch.data() + 2 * neq;
      double* N4  = scratch.data() + 3 * neq;
      double* U2  = scratch.data() + 4 * neq;
      double* U3  = scratch.data() + 5 * neq;
      double* U4  = scratch.data() + 6 * neq;
      double* tmp = scratch.data() + 7 * neq;

      // N(u) = f(u) - L u with L = diag(J). The whole scheme is exact whenever N is
      // constant, which is the late-time collision: a fixed relaxation rate onto a
      // slowly-drifting attractor.
      for (int k = 0; k < neq; k++) {
        double e0, p1, p2, p3;
        PhiFuncs(0.5 * hh * diag[k], e0, p1, p2, p3);
        N1[k] = fin[k] - diag[k] * yin[k];
        U2[k] = e0 * yin[k] + 0.5 * hh * p1 * N1[k];
      }
      (*derivs)(tin + 0.5 * hh, U2, tmp, parameters_and_workspace_for_derivs);
      for (int k = 0; k < neq; k++) {
        double e0, p1, p2, p3;
        PhiFuncs(0.5 * hh * diag[k], e0, p1, p2, p3);
        N2[k] = tmp[k] - diag[k] * U2[k];
        U3[k] = e0 * yin[k] + 0.5 * hh * p1 * N2[k];
      }
      (*derivs)(tin + 0.5 * hh, U3, tmp, parameters_and_workspace_for_derivs);
      for (int k = 0; k < neq; k++) {
        double e0, p1, p2, p3;
        PhiFuncs(0.5 * hh * diag[k], e0, p1, p2, p3);
        N3[k] = tmp[k] - diag[k] * U3[k];
        U4[k] = e0 * U2[k] + 0.5 * hh * p1 * (2. * N3[k] - N1[k]);
      }
      (*derivs)(tin + hh, U4, tmp, parameters_and_workspace_for_derivs);
      for (int k = 0; k < neq; k++) {
        double e0, p1, p2, p3;
        PhiFuncs(hh * diag[k], e0, p1, p2, p3);
        N4[k]           = tmp[k] - diag[k] * U4[k];
        const double c1 = p1 - 3. * p2 + 4. * p3;
        const double c2 = 2. * p2 - 4. * p3;
        const double c4 = -p2 + 4. * p3;
        yout[k]         = e0 * yin[k] + hh * (c1 * N1[k] + c2 * (N2[k] + N3[k]) + c4 * N4[k]);
      }
    };

    etdrk4(y, f0.data(), t, h, ynew.data());

    /* Endpoint derivative for the embedded estimate and, if accepted, for the next
       step's f0. Let N5 = f(t+h, ynew) - L ynew. The third-order companion replaces
       the Cox-Matthews c4(z) * N4 contribution by c4(z) * N5:

         y_emb = y_high + h c4(z) (N5 - N4).

       At L=0 this is the classical RK4 tableau with a fifth endpoint-derivative stage
       and weights [1/6, 1/3, 1/3, 0, 1/6]. For a constant remainder N, N4=N5, so the
       stiff scalar relaxation stays exact and does not limit the step. */
    (*derivs)(t + h, ynew.data(), fnew.data(), parameters_and_workspace_for_derivs);
    {
      const double* N4 = scratch.data() + 3 * neq;
      for (int k = 0; k < neq; k++) {
        double e0, p1, p2, p3;
        PhiFuncs(h * diag[k], e0, p1, p2, p3);
        const double c4 = -p2 + 4. * p3;
        const double N5 = fnew[k] - diag[k] * ynew[k];
        err[k]          = h * c4 * (N4[k] - N5);
      }
    }

    /* Same accumulate-inside-the-reduction discipline as rkdp45: `errtemp > errmax`
       is FALSE for NaN, so a NaN component would leave errmax at 0 and the step
       would be ACCEPTED, completing the run silently with a poisoned state. */
    double errmax = 0.0;
    int nonfinite = 0;
    for (int k = 0; k < neq; k++) {
      const double errtemp  = std::fabs(err[k] / std::max(threshold, std::fabs(ynew[k])));
      nonfinite            |= (int) (IsNonFinite(ynew[k]) || IsNonFinite(err[k]));
      if (errtemp > errmax)
        errmax = errtemp;
    }

    if (nonfinite != 0) {
      int kbad = 0;
      for (int k = 0; k < neq; k++) {
        if (IsNonFinite(ynew[k]) || IsNonFinite(err[k])) {
          kbad = k;
          break;
        }
      }
      class_stop(
          "etd: non-finite state at x=%e (h=%e): y[%d]=%e, err[%d]=%e, diag[%d]=%e. "
          "A non-finite DIAGONAL points at the species' BackgroundDerivsDiagonal "
          "rather than at this evolver.",
          t,
          h,
          kbad,
          ynew[kbad],
          kbad,
          err[kbad],
          kbad,
          diag[kbad]);
    }

    if (errmax > rtol) {
      if (nofailed) {
        nofailed = false;
        hnew     = tdir *
                   std::max(hmin, std::fabs(h) * std::max(0.1, 0.8 * pow(rtol / errmax, pow_grow)));
      }
      else {
        hnew = tdir * std::max(hmin, 0.5 * std::fabs(h));
      }
      continue;
    }

    /* Step accepted. fnew was already computed for the embedded estimate and becomes
       both this step's reported derivative and the next step's f0. */
    const double tnew = t + h;

    if (print_variables != nullptr)
      (*print_variables)(tnew, ynew.data(), fnew.data(), parameters_and_workspace_for_derivs);

    /* Dense output. The two endpoint cases are exact and cost nothing; the interior
       one is the exponential continuous extension described below. */
    for (; (idx < x_size) && ((tnew - x_sampling[idx]) * tdir >= 0.0); idx++) {
      if (x_sampling[idx] == t) {
        (*output)(t, y, f0.data(), idx, parameters_and_workspace_for_derivs);
      }
      else if (tnew == x_sampling[idx]) {
        (*output)(tnew, ynew.data(), fnew.data(), idx, parameters_and_workspace_for_derivs);
      }
      else {
        /* Interpolate the REMAINDER N quadratically in s through its three stage
           samples -- N1 at s=0, (N2+N3)/2 at s=h/2, N4 at s=h -- and then integrate
           the variation-of-constants formula against that quadratic EXACTLY. With
           P(s/h) = A + B*th + C*th^2 and th = (ti - t)/h,

             y(t + th*h) = e^{th z} y_n + th*h [ A phi1 + B th phi2 + 2 C th^2 phi3 ],

           every phi evaluated at th*z. The three integrals used are the standard
           phi identities int_0^1 e^{w(1-s)} s^{k-1}/(k-1)! ds = phi_k(w).

           At th=1 the weights collapse to Cox-Matthews' own c1, c2, c4, so the curve
           meets the accepted step exactly rather than merely closely; at th=0 it
           returns y_n. For constant N it is the exact solution at every th, like the
           step itself.

           This is what a cubic Hermite through the step's endpoints cannot do: when
           the step deliberately jumps a relaxation layer of width 1/|diag|, no
           polynomial through the two ends knows the layer is there, and the SAMPLED
           table -- which in the perturbations is the entire deliverable -- is then far
           worse than the trajectory that produced it. Here the layer is carried by
           e^{th z} and the phi's, exactly as in the step. It costs no RHS evaluations:
           the stages are still in `scratch`, and dy follows from dy = L y + N. */
        const double ti  = x_sampling[idx];
        const double th  = (ti - t) / h;
        const double th2 = th * th;
        const double* N1 = scratch.data() + 0 * neq;
        const double* N2 = scratch.data() + 1 * neq;
        const double* N3 = scratch.data() + 2 * neq;
        const double* N4 = scratch.data() + 3 * neq;
        for (int k = 0; k < neq; k++) {
          if (used_in_output[k] == 0)
            continue;
          double e0, p1, p2, p3;
          PhiFuncs(th * h * diag[k], e0, p1, p2, p3);
          const double mid = N2[k] + N3[k];
          const double A   = N1[k];
          const double B   = -3. * N1[k] + 2. * mid - N4[k];
          const double C   = 2. * N1[k] - 2. * mid + 2. * N4[k];
          yinterp[k]       = e0 * y[k] + th * h * (A * p1 + B * th * p2 + 2. * C * th2 * p3);
          dyinterp[k]      = diag[k] * yinterp[k] + (A + B * th + C * th2);
        }
        (*output)(ti, yinterp.data(), dyinterp.data(), idx, parameters_and_workspace_for_derivs);
      }
    }

    nofailed = true;
    /* The growth factor is capped, which rkdp45 does not need to do: there errmax is
       never exactly zero, but here it is whenever the remainder N is constant over the
       step -- the case this scheme integrates exactly -- and rtol/0 would hand the next
       step an infinity. It survives only because the x_end clamp catches it, and a
       single jump to the end of the interval is right for a genuinely linear problem
       but is repeatedly rejected and halved for a nearly-linear one. */
    constexpr double kMaxGrowth = 10.0;
    hnew = tdir *
           std::max(hmin,
                    std::fabs(h) *
                        std::max(0.1, std::min(kMaxGrowth, 0.8 * pow(rtol / errmax, pow_grow))));
    for (int k = 0; k < neq; k++) {
      y[k]  = ynew[k];
      f0[k] = fnew[k];
    }
    t = tnew;
  }
}
