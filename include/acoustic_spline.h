#ifndef ACOUSTIC_SPLINE_H
#define ACOUSTIC_SPLINE_H

#include <vector>

/**
 * Interpolation in l of a CMB power spectrum known at sampled multipoles, built for its acoustic
 * oscillation.
 *
 * C_l is sampled every l_linstep multipoles and interpolated in between. A cubic spline treats the
 * acoustic oscillation as an error, which grows as (h omega)^4. This interpolates
 * u_l = l(l+1) C_l / S(l) instead, with
 *
 *  - S(l) a smooth positive envelope: exp of a Chebyshev polynomial in log l, fitted by least squares
 *    to the log of the running maximum of |l(l+1) C_l| over one period. It takes out the damping,
 *    so that what is left oscillates with a nearly constant amplitude.
 *  - a spline in the space {1, l, cos(omega l), sin(omega l)} on each interval (a spline under
 *    imaginary tension), C2 at the nodes, with omega = 2 theta_s, the frequency of the acoustic
 *    peaks in l. It is exact for a sinusoid of that frequency plus a linear trend, and becomes the
 *    cubic spline as omega -> 0. Its end conditions keep u'' + omega^2 u linear over the last three
 *    nodes, which every function of the space satisfies, so the ends are exact for it too.
 *
 * It passes through every sample, so only the values between them change. Nothing is fitted
 * across the oscillation: there is no template, and a model whose peaks move is handled like any
 * other. docs/superpowers/specs/2026-10-27-acoustic-l-interpolation-design.md
 */
class AcousticSpline {
 public:
  /** l[0..n-1]: the sampled multipoles, increasing; cl: C_l at them, cl[i * stride]; omega: 2 theta_s. */
  AcousticSpline(const double* l, int n, const double* cl, int stride, double omega);

  /** C_l at l, for l in [l.front(), l.back()]. */
  double operator()(double l) const;

  /** Whether this interpolation applies: false for a spectrum that is zero everywhere, and for
      samples too sparse for the frequency (omega h near pi makes the spline singular). The
      caller keeps its cubic spline then. */
  static bool Applies(const double* l, int n, const double* cl, int stride, double omega);

 private:
  double Envelope(double l) const;

  std::vector<double> l_;         /* sampled multipoles */
  std::vector<double> u_;         /* l(l+1) C_l / S(l) at them */
  std::vector<double> m_;         /* u'' at them, from the spline equations */
  std::vector<double> chebyshev_; /* coefficients of log S in t(log l) */
  double log_l_lo_ = 0.;
  double log_l_hi_ = 0.;
  double omega_    = 0.;
};

#endif  // ACOUSTIC_SPLINE_H
