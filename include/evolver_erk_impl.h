#ifndef __EVOERKIMPL__
#define __EVOERKIMPL__

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <vector>

#include "common.h"  // class_test, class_stop
#include "evolver_erk.h"
#include "nonfinite.h"

/**
 * The one driver behind every embedded explicit RK evolver in CLASS.
 *
 * This is an IMPLEMENTATION header, included by exactly the thin per-method
 * wrappers (evolver_rkdp45.cpp, evolver_tsit5.cpp) and by the unit test. It is a
 * template rather than an ordinary function so each tableau's coefficients stay
 * compile-time constants inside the stage loops -- see the note in
 * evolver_erk.h on why that is load-bearing rather than cosmetic.
 */
template <class Tab>
void evolver_erk_run(
    void (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double* x_sampling,
    int x_size,
    void (*output)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace),
    void (*print_variables)(double x, double y[], double dy[], void* parameters_and_workspace)) {
  /* Every step-size rule below reads cfg.safety and cfg.fac_min, including the
     `legacy` ones, which used to hard-code 0.8 and 0.1. At the defaults those ARE
     0.8 and 0.1, so `legacy` is unchanged; the point is that a knob the input
     documents as the safety factor should not be silently ignored by the shipped
     controller. */

  /* Contract (shared with evolver_rk): x_sampling is a non-null array of
     output points, and tolerance > 0 (threshold below divides by it). */
  class_test(x_sampling == nullptr, "%s requires a non-null x_sampling array", Tab::kName);
  class_test(tolerance <= 0., "%s requires tolerance > 0 (got %e)", Tab::kName, tolerance);

  const ErkControllerConfig cfg = evolver_erk_config(); /* one snapshot per call */
  const bool track              = evolver_erk_stats_enabled();

  const int neq          = y_size;
  const double rtol      = tolerance;
  const double abstol    = 1e-15; /* matches ndf15 */
  const double threshold = abstol / rtol;
  const double pow_grow  = 1.0 / Tab::kOrder;
  constexpr int s        = Tab::kStages;
  constexpr int deg      = Tab::kDenseDeg;

  /* The derivative at the right end of an accepted step. FSAL pairs already have
     it as the last stage; otherwise it costs one extra evaluation, which also
     supplies the next step's k0. */
  constexpr int idx_end = Tab::kFsal ? (s - 1) : s;
  constexpr int k_rows  = Tab::kFsal ? s : (s + 1);

  std::vector<double> ynew(neq), ytemp(neq), err(neq), yinterp(neq), dyinterp(neq);
  std::vector<double> ki(static_cast<size_t>(k_rows) * neq);
  double bi_vec_y[s], bi_vec_dy[s], thpow[deg + 1];

  /* History for the quintic Hermite interpolant: the value and the derivative at
     the END of the step before this one. Allocated only when that interpolant is
     selected -- it is two more vectors of the state size per workspace. */
  const bool want_hermite = (cfg.interpolant == ErkInterpolant::hermite3);
  std::vector<double> y_prev, f_prev;
  double t_prev  = 0.0;
  bool have_prev = false;
  if (want_hermite) {
    y_prev.resize(neq);
    f_prev.resize(neq);
  }

  double t = x_ini;

  /* initialise k0 = f(t, y) */
  (*derivs)(t, y, ki.data(), parameters_and_workspace_for_derivs);
  if (track)
    erk_detail::CountDerivs();

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
  rh /= cfg.safety * pow(rtol, pow_grow);
  if (absh * rh > 1.0)
    absh = 1.0 / rh;

  const int tdir = (x_end > x_ini) ? 1 : -1;
  double hnew    = absh * tdir;
  bool nofailed  = true;

  /* PI memory: the previous ACCEPTED step's error ratio. Hairer's DOPRI5 seeds
     this at 1e-4 so the first steps are not throttled by an absent history. */
  double ratio_old      = 1e-4;
  const double alpha_pi = pow_grow - 0.75 * cfg.pi_beta;

  int idx = 0;
  if ((t - x_end) * tdir < 0.0)
    for (idx = 0; (idx < x_size) && ((x_sampling[idx] - t) * tdir < 0.0); idx++)
      ;

  const bool step_to_output = (cfg.output_stepping == ErkOutputStepping::exact);

  while ((t - x_end) * tdir < 0.0) {
    double h          = hnew;
    const double hmin = 100.0 * DBL_MIN * fabs(t);
    class_test(fabs(h) < hmin,
               "%s: step size %e fell below minimum %e at x=%e",
               Tab::kName,
               fabs(h),
               hmin,
               t);
    if (fabs(h) > 0.9 * fabs(x_end - t))
      h = x_end - t;

    /* In `exact` output stepping the step is truncated so it ends on the next
       sampling point. `clamped` records that the length was imposed rather than
       chosen: the controller must then NOT read its own proposal back from this
       step, because a truncated step has a smaller error for reasons that say
       nothing about how large the next step may be. */
    bool clamped     = false;
    double t_clamped = 0.0;
    if (step_to_output && idx < x_size) {
      const double dt_out = (x_sampling[idx] - t) * tdir;
      /* `<=` and not `<`: when the step already ends on the point, t + h still
         only lands within a rounding of it, which is enough to miss the equality
         test at emission and push an exactly-known state through the interpolant. */
      if (dt_out > 0.0 && dt_out <= fabs(h)) {
        h         = tdir * dt_out;
        clamped   = true;
        t_clamped = x_sampling[idx];
      }
    }

    /* stage 0 (FSAL: ki[0..neq) already holds k0) */
    for (int k = 0; k < neq; k++) {
      ynew[k] = y[k] + h * Tab::b[0] * ki[k];
      err[k]  = h * Tab::e[0] * ki[k];
    }
    for (int i = 1; i < s; i++) {
      for (int k = 0; k < neq; k++)
        ytemp[k] = y[k];
      for (int j = 0; j < i; j++) {
        const double aij = Tab::A[i * s + j];
        for (int k = 0; k < neq; k++)
          ytemp[k] += h * aij * ki[j * neq + k];
      }
      (*derivs)(t + Tab::c[i] * h,
                ytemp.data(),
                ki.data() + i * neq,
                parameters_and_workspace_for_derivs);
      if (track)
        erk_detail::CountDerivs();
      for (int k = 0; k < neq; k++) {
        ynew[k] += h * Tab::b[i] * ki[i * neq + k];
        err[k]  += h * Tab::e[i] * ki[i * neq + k];
      }
    }

    /* The non-finite flag is accumulated INSIDE the error reduction, not tested on
       errmax afterwards, because the reduction cannot carry a NaN out:
       `errtemp > errmax` is FALSE for NaN, so a NaN component leaves errmax at 0 and
       the step is ACCEPTED -- h then grows and the run completes, silently, with a
       poisoned state. Inf behaves differently (inf > x is true), which is why the
       observed decaying-NCDM failure hung on rejection instead of returning garbage.
       Both modes have to be caught, so the test is per component. */
    double errmax = 0.0;
    int nonfinite = 0;
    for (int k = 0; k < neq; k++) {
      double errtemp  = fabs(err[k] / std::max(threshold, fabs(ynew[k])));
      nonfinite      |= (int) (IsNonFinite(ynew[k]) || IsNonFinite(err[k]));
      if (errtemp > errmax)
        errmax = errtemp;
    }
    /* `errmax` stays the quantity reported in diagnostics; `errnorm` is what the
       controller acts on. They coincide unless the RMS norm is selected, which
       is a second reduction rather than a modification of the first so that the
       max-norm path keeps its exact historical arithmetic. */
    double errnorm = errmax;
    if (cfg.norm == ErkErrorNorm::rms_norm) {
      double acc = 0.0;
      for (int k = 0; k < neq; k++) {
        const double e_k  = err[k] / std::max(threshold, fabs(ynew[k]));
        acc              += e_k * e_k;
      }
      errnorm = sqrt(acc / neq);
    }

    /* A non-finite state must ABORT, not be handed back to the step controller. The
       O(neq) search for the offending index runs only on the failure path.

       On the Inf path the run does not stop, it GRINDS: errmax > rtol is taken, the
       step is rejected, h shrinks, and the only exit is class_test(fabs(h) < hmin)
       with hmin = 100*DBL_MIN*|x|. h converges ONTO that value without going below
       it, so the guard never fires. Measured on the decaying-NCDM sector at Gamma=1e8
       (dr_N_q=30, k=0.5): NaN at x=523.76, then h pinned at 1.1654e-303 with 99.4% of
       steps rejected, repeating identically -- two multi-hour runs burned on what
       should have been an immediate error.

       IsNonFinite is a bit test rather than std::isnan because classpp is built with
       -ffast-math, under which std::isnan may be folded to false (see nonfinite.h). */
    if (nonfinite != 0) {
      int kbad = 0;
      for (int k = 0; k < neq; k++) {
        if (IsNonFinite(ynew[k]) || IsNonFinite(err[k])) {
          kbad = k;
          break;
        }
      }
      class_stop(
          "%s: non-finite state at x=%e (h=%e, errmax=%e): y[%d]=%e, err[%d]=%e. The "
          "integration cannot recover -- rejecting and shrinking h only converges onto the "
          "step-size floor. Something upstream produced NaN/Inf.",
          Tab::kName,
          t,
          h,
          errmax,
          kbad,
          ynew[kbad],
          kbad,
          err[kbad]);
    }

    if (errnorm > rtol) {
      if (track)
        erk_detail::CountRejected(t, errnorm / rtol);
      if (cfg.kind == ErkControllerKind::legacy) {
        /* Standard ode45 control: first failure shrinks by the error-proportional
           factor and marks nofailed; consecutive failures halve. (This corrects a
           transcription bug in the PhD2024-EBH branch where nofailed was never
           set, making the halving path dead code.) */
        if (nofailed) {
          nofailed = false;
          hnew = tdir * std::max(hmin,
                                 fabs(h) * std::max(cfg.fac_min,
                                                    cfg.safety * pow(rtol / errnorm, pow_grow)));
        }
        else {
          hnew = tdir * std::max(hmin, 0.5 * fabs(h));
        }
      }
      else {
        double fac = cfg.safety * pow(rtol / errnorm, pow_grow);
        fac        = std::max(cfg.fac_min, std::min(fac, 1.0));
        if (!nofailed)
          fac = std::min(fac, 0.5); /* repeated failures: force a real cut */
        nofailed = false;
        hnew     = tdir * std::max(hmin, fabs(h) * fac);
      }
      continue;
    }

    /* step accepted */
    if (track)
      erk_detail::CountAccepted(t, errnorm / rtol);

    if (!Tab::kFsal) {
      for (int k = 0; k < neq; k++)
        ytemp[k] = ynew[k];
      (*derivs)(t + h,
                ytemp.data(),
                ki.data() + idx_end * neq,
                parameters_and_workspace_for_derivs);
      if (track)
        erk_detail::CountDerivs();
    }

    if (print_variables != nullptr)
      (*print_variables)(t + h,
                         ynew.data(),
                         ki.data() + idx_end * neq,
                         parameters_and_workspace_for_derivs);

    if (clamped) {
      /* keep the controller's own proposal; this step's length was imposed */
      hnew = tdir * std::max(hmin, fabs(hnew));
    }
    else if (cfg.kind == ErkControllerKind::legacy) {
      hnew = tdir *
             std::max(hmin,
                      fabs(h) * std::max(cfg.fac_min, cfg.safety * pow(rtol / errnorm, pow_grow)));
    }
    else {
      const double ratio = errnorm / rtol;
      double fac;
      if (cfg.kind == ErkControllerKind::pi && ratio > 0.0)
        fac = cfg.safety * pow(ratio, -alpha_pi) * pow(ratio_old, cfg.pi_beta);
      else
        fac = cfg.safety * pow(rtol / errnorm, pow_grow);
      fac = std::max(cfg.fac_min, std::min(fac, cfg.fac_max));
      /* Never grow immediately after a rejection: the step that just failed is
         evidence that the local scale shrank. Measured on LCDM + Mnu this is
         worth little on its own -- with the growth cap it moves the rejection
         rate only from 13.3% to 12.1% -- so it is the PI memory term above, not
         this, that takes rejections down to ~6%. */
      if (!nofailed)
        fac = std::min(fac, 1.0);
      hnew      = tdir * std::max(hmin, fabs(h) * fac);
      ratio_old = std::max(ratio, 1e-4);
    }
    nofailed = true;
    /* A truncated step ends ON the sampling point by construction. Recomputing
       t + h instead would land within a rounding of it, which is enough to miss
       the equality test below and send an exactly-known state through the
       interpolant. */
    const double tnew = clamped ? t_clamped : (t + h);

    /* emit output at all sampling points within (t, tnew] */
    for (; (idx < x_size) && ((tnew - x_sampling[idx]) * tdir >= 0.0); idx++) {
      if (tnew == x_sampling[idx]) {
        if (track)
          erk_detail::CountExact();
        (*output)(tnew,
                  ynew.data(),
                  ki.data() + idx_end * neq,
                  idx,
                  parameters_and_workspace_for_derivs);
      }
      else {
        if (track)
          erk_detail::CountDense();
        const double ti = x_sampling[idx];
        const double th = (ti - t) / h;
        /* theta^m, m = 0..deg. Squaring for m = 2 and m = 4 is not an
           optimisation: it is the association the hand-written Dormand-Prince
           interpolant used, and keeping it makes this driver reproduce the
           previous rkdp45 output bit for bit. */
        thpow[0] = 1.0;
        thpow[1] = th;
        if (deg >= 2)
          thpow[2] = th * th;
        if (deg >= 3)
          thpow[3] = thpow[2] * th;
        if (deg >= 4)
          thpow[4] = thpow[2] * thpow[2];
        for (int m = 5; m <= deg; m++)
          thpow[m] = thpow[m - 1] * th;
        for (int i = 0; i < s; i++) {
          double vy = 0.0, vd = 0.0;
          for (int m = 0; m < deg; m++) {
            const double r  = Tab::R[i * deg + m];
            vy             += r * thpow[m + 1];
            vd             += r * (m + 1) * thpow[m];
          }
          bi_vec_y[i]  = vy;
          bi_vec_dy[i] = vd;
        }
        for (int k = 0; k < neq; k++) {
          if (used_in_output[k]) {
            yinterp[k]  = y[k];
            dyinterp[k] = 0.0;
            for (int i = 0; i < s; i++) {
              yinterp[k]  += h * bi_vec_y[i] * ki[i * neq + k];
              dyinterp[k] += bi_vec_dy[i] * ki[i * neq + k];
            }
          }
        }

        if (want_hermite && have_prev) {
          /* Quintic Hermite through (t_prev, t, tnew), value and derivative at
             each. The six basis functions depend only on the node positions and
             on ti, so they are built once and then applied to every component --
             the same cost per component as the tableau extension it replaces.

             h_i(x)    = [1 - 2 (x - x_i) L_i'(x_i)] L_i(x)^2
             hbar_i(x) = (x - x_i) L_i(x)^2
             with L_i the Lagrange basis of the three nodes. */
          const double x0 = t_prev, x1 = t, x2 = tnew;
          const double d01 = x0 - x1, d02 = x0 - x2, d12 = x1 - x2;
          /* Coincident nodes would make this singular; they cannot occur (each
             node is a strictly later accepted step end) but a zero-length step
             would, so fall back rather than divide by it. */
          if (d01 != 0.0 && d02 != 0.0 && d12 != 0.0) {
            const double u0 = ti - x0, u1 = ti - x1, u2 = ti - x2;
            const double L0 = (u1 * u2) / (d01 * d02);
            const double L1 = (u0 * u2) / (-d01 * d12);
            const double L2 = (u0 * u1) / (d02 * d12);
            /* L_i'(x_i) = sum_{j != i} 1/(x_i - x_j) */
            const double Lp0 = 1.0 / d01 + 1.0 / d02;
            const double Lp1 = -1.0 / d01 + 1.0 / d12;
            const double Lp2 = -1.0 / d02 - 1.0 / d12;
            /* derivative of L_i at ti, for the interpolant's own derivative */
            const double dL0 = (u1 + u2) / (d01 * d02);
            const double dL1 = (u0 + u2) / (-d01 * d12);
            const double dL2 = (u0 + u1) / (d02 * d12);

            const double L0s = L0 * L0, L1s = L1 * L1, L2s = L2 * L2;
            const double w0 = (1.0 - 2.0 * u0 * Lp0) * L0s;
            const double w1 = (1.0 - 2.0 * u1 * Lp1) * L1s;
            const double w2 = (1.0 - 2.0 * u2 * Lp2) * L2s;
            const double v0 = u0 * L0s, v1 = u1 * L1s, v2 = u2 * L2s;
            const double dw0 = -2.0 * Lp0 * L0s + (1.0 - 2.0 * u0 * Lp0) * 2.0 * L0 * dL0;
            const double dw1 = -2.0 * Lp1 * L1s + (1.0 - 2.0 * u1 * Lp1) * 2.0 * L1 * dL1;
            const double dw2 = -2.0 * Lp2 * L2s + (1.0 - 2.0 * u2 * Lp2) * 2.0 * L2 * dL2;
            const double dv0 = L0s + u0 * 2.0 * L0 * dL0;
            const double dv1 = L1s + u1 * 2.0 * L1 * dL1;
            const double dv2 = L2s + u2 * 2.0 * L2 * dL2;

            const double* f1 = ki.data();                 /* f at t     (FSAL k0) */
            const double* f2 = ki.data() + idx_end * neq; /* f at tnew            */
            for (int k = 0; k < neq; k++) {
              if (used_in_output[k]) {
                yinterp[k]  = w0 * y_prev[k] + w1 * y[k] + w2 * ynew[k] + v0 * f_prev[k] +
                              v1 * f1[k] + v2 * f2[k];
                dyinterp[k] = dw0 * y_prev[k] + dw1 * y[k] + dw2 * ynew[k] + dv0 * f_prev[k] +
                              dv1 * f1[k] + dv2 * f2[k];
              }
            }
          }
        }
        (*output)(ti, yinterp.data(), dyinterp.data(), idx, parameters_and_workspace_for_derivs);
      }
    }

    if (want_hermite) {
      /* This step's LEFT end becomes the history for the next step. It must be
         copied before y and ki are overwritten below -- and only for the
         components the interpolant can ever be asked about, since this copy runs
         on every accepted step while the interpolant runs on roughly half of
         them. */
      for (int k = 0; k < neq; k++) {
        if (used_in_output[k]) {
          y_prev[k] = y[k];
          f_prev[k] = ki[k];
        }
      }
      t_prev    = t;
      have_prev = true;
    }
    for (int k = 0; k < neq; k++) {
      y[k]  = ynew[k];
      ki[k] = ki[idx_end * neq + k]; /* FSAL: last stage becomes next k0 */
    }
    t = tnew;
  }
}

#endif
