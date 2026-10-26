/** @file acoustic_spline.cpp Interpolation in l of CMB power spectra, built for their acoustic
 * oscillation: see include/acoustic_spline.h and
 * docs/superpowers/specs/2026-10-27-acoustic-l-interpolation-design.md
 */

#include "acoustic_spline.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

#include "constants.h"

namespace {

/* Envelope: degree of the Chebyshev polynomial in log l, and how far below the largest |D| the
   running maximum may fall before it is floored (a window where the spectrum vanishes). */
constexpr int kEnvelopeDegree   = 6;
constexpr double kEnvelopeFloor = 1e-30;
/* omega h must stay away from pi, where the spline equations are singular */
constexpr double kMaxPhaseStep = 2.0;

double Chebyshev(int k, double t) {
  double t_km1 = 1.;
  double t_k   = t;
  if (k == 0) {
    return 1.;
  }
  for (int j = 1; j < k; j++) {
    const double t_kp1 = 2. * t * t_k - t_km1;
    t_km1              = t_k;
    t_k                = t_kp1;
  }
  return t_k;
}

}  // namespace

bool AcousticSpline::Applies(const double* l, int n, const double* cl, int stride, double omega) {
  if ((n < kEnvelopeDegree + 2) || !(omega > 0.)) {
    return false;
  }
  bool nonzero = false;
  for (int i = 0; i < n; i++) {
    nonzero = nonzero || (cl[static_cast<size_t>(i) * stride] != 0.);
    if ((i > 0) && (omega * (l[i] - l[i - 1]) > kMaxPhaseStep)) {
      return false;
    }
  }
  return nonzero;
}

AcousticSpline::AcousticSpline(const double* l, int n, const double* cl, int stride, double omega)
    : omega_(omega) {
  if (!Applies(l, n, cl, stride, omega)) {
    throw std::logic_error("AcousticSpline built where it does not apply (see Applies)");
  }
  l_.resize(n);
  std::vector<double> D(n);
  double D_max = 0.;
  for (int i = 0; i < n; i++) {
    l_[i] = l[i];
    D[i]  = cl[static_cast<size_t>(i) * stride] * l_[i] * (l_[i] + 1.);
    D_max = std::max(D_max, std::fabs(D[i]));
  }

  /* envelope: log of the running maximum of |D| over one period, fitted by a Chebyshev polynomial
     in log l (normal equations; the basis keeps them well conditioned at this degree) */
  const double half_period = _PI_ / omega_;
  log_l_lo_                = std::log(l_.front());
  log_l_hi_                = std::log(l_.back());
  const int degree         = kEnvelopeDegree;
  std::array<std::array<double, kEnvelopeDegree + 2>, kEnvelopeDegree + 1> normal{};
  for (int i = 0; i < n; i++) {
    double running_max = 0.;
    for (int j = 0; j < n; j++) {
      if (std::fabs(l_[j] - l_[i]) <= half_period) {
        running_max = std::max(running_max, std::fabs(D[j]));
      }
    }
    const double log_env = std::log(std::max(running_max, kEnvelopeFloor * D_max));
    const double t       = 2. * (std::log(l_[i]) - log_l_lo_) / (log_l_hi_ - log_l_lo_) - 1.;
    std::array<double, kEnvelopeDegree + 1> basis{};
    for (int k = 0; k <= degree; k++) {
      basis[k] = Chebyshev(k, t);
    }
    for (int r = 0; r <= degree; r++) {
      for (int c = 0; c <= degree; c++) {
        normal[r][c] += basis[r] * basis[c];
      }
      normal[r][degree + 1] += basis[r] * log_env;
    }
  }
  /* Gaussian elimination with partial pivoting on the (degree+1) x (degree+2) augmented system */
  for (int c = 0; c <= degree; c++) {
    int pivot = c;
    for (int r = c + 1; r <= degree; r++) {
      if (std::fabs(normal[r][c]) > std::fabs(normal[pivot][c])) {
        pivot = r;
      }
    }
    std::swap(normal[c], normal[pivot]);
    for (int r = c + 1; r <= degree; r++) {
      const double f = normal[r][c] / normal[c][c];
      for (int k = c; k <= degree + 1; k++) {
        normal[r][k] -= f * normal[c][k];
      }
    }
  }
  chebyshev_.assign(degree + 1, 0.);
  for (int r = degree; r >= 0; r--) {
    double s = normal[r][degree + 1];
    for (int k = r + 1; k <= degree; k++) {
      s -= normal[r][k] * chebyshev_[k];
    }
    chebyshev_[r] = s / normal[r][r];
  }

  u_.resize(n);
  for (int i = 0; i < n; i++) {
    u_[i] = D[i] / Envelope(l_[i]);
  }

  /* The spline: on [l_i, l_i+1], with a = l_i+1 - l, b = l - l_i,
       s = -(m_i sin(w a) + m_i+1 sin(w b)) / (w^2 sin(w h)) + (u_i + m_i/w^2) a/h + (u_i+1 + m_i+1/w^2) b/h,
     so s(l_i) = u_i and s''(l_i) = m_i. Continuity of s' at the interior nodes gives
       alpha_i-1 m_i-1 + (beta_i-1 + beta_i) m_i + alpha_i m_i+1 = (u_i+1 - u_i)/h_i - (u_i - u_i-1)/h_i-1,
     with alpha -> h/6 and beta -> h/3 as w -> 0 (the cubic spline). */
  const double w2 = omega_ * omega_;
  std::vector<double> alpha(n - 1);
  std::vector<double> beta(n - 1);
  for (int i = 0; i < n - 1; i++) {
    const double h  = l_[i + 1] - l_[i];
    const double wh = omega_ * h;
    alpha[i]        = (omega_ / std::sin(wh) - 1. / h) / w2;
    beta[i]         = (1. / h - omega_ * std::cos(wh) / std::sin(wh)) / w2;
  }
  /* End conditions: g = m + w^2 u collinear over the first and over the last three nodes (exact
     for every function of the space), each eliminated into its neighbouring row. */
  auto collinear = [&](int i0, int i1, int i2, double& c0, double& c1, double& c2, double& rhs) {
    c0  = l_[i2] - l_[i1];
    c1  = -(l_[i2] - l_[i0]);
    c2  = l_[i1] - l_[i0];
    rhs = -w2 * (c0 * u_[i0] + c1 * u_[i1] + c2 * u_[i2]);
  };
  double c0 = 0.;
  double c1 = 0.;
  double c2 = 0.;
  double r0 = 0.;
  double d0 = 0.;
  double d1 = 0.;
  double d2 = 0.;
  double rn = 0.;
  collinear(0, 1, 2, c0, c1, c2, r0);
  collinear(n - 1, n - 2, n - 3, d0, d1, d2, rn);

  /* tridiagonal system for m_1 .. m_n-2: sub, diag, super, rhs */
  const int nn = n - 2;
  std::vector<double> sub(nn);
  std::vector<double> diag(nn);
  std::vector<double> sup(nn);
  std::vector<double> rhs(nn);
  for (int k = 0; k < nn; k++) {
    const int i = k + 1;
    sub[k]      = alpha[i - 1];
    diag[k]     = beta[i - 1] + beta[i];
    sup[k]      = alpha[i];
    rhs[k] = (u_[i + 1] - u_[i]) / (l_[i + 1] - l_[i]) - (u_[i] - u_[i - 1]) / (l_[i] - l_[i - 1]);
  }
  /* m_0 = (r0 - c1 m_1 - c2 m_2) / c0 into row 1; m_n-1 = (rn - d1 m_n-2 - d2 m_n-3) / d0 into row n-2 */
  diag[0] -= sub[0] * c1 / c0;
  if (nn > 1) {
    sup[0] -= sub[0] * c2 / c0;
  }
  rhs[0]       -= sub[0] * r0 / c0;
  diag[nn - 1] -= sup[nn - 1] * d1 / d0;
  if (nn > 1) {
    sub[nn - 1] -= sup[nn - 1] * d2 / d0;
  }
  rhs[nn - 1] -= sup[nn - 1] * rn / d0;
  /* Thomas algorithm (the interior rows are diagonally dominant for w h < kMaxPhaseStep) */
  for (int k = 1; k < nn; k++) {
    const double f  = sub[k] / diag[k - 1];
    diag[k]        -= f * sup[k - 1];
    rhs[k]         -= f * rhs[k - 1];
  }
  m_.assign(n, 0.);
  m_[nn] = rhs[nn - 1] / diag[nn - 1];
  for (int k = nn - 2; k >= 0; k--) {
    m_[k + 1] = (rhs[k] - sup[k] * m_[k + 2]) / diag[k];
  }
  m_[0]     = (r0 - c1 * m_[1] - c2 * m_[2]) / c0;
  m_[n - 1] = (rn - d1 * m_[n - 2] - d2 * m_[n - 3]) / d0;
}

double AcousticSpline::Envelope(double l) const {
  const double t = 2. * (std::log(l) - log_l_lo_) / (log_l_hi_ - log_l_lo_) - 1.;
  /* Clenshaw recurrence for sum_k c_k T_k(t) */
  double b1 = 0.;
  double b2 = 0.;
  for (int k = static_cast<int>(chebyshev_.size()) - 1; k >= 1; k--) {
    const double b0 = 2. * t * b1 - b2 + chebyshev_[k];
    b2              = b1;
    b1              = b0;
  }
  return std::exp(t * b1 - b2 + chebyshev_[0]);
}

double AcousticSpline::operator()(double l) const {
  const int n     = static_cast<int>(l_.size());
  int i           = static_cast<int>(std::upper_bound(l_.begin(), l_.end(), l) - l_.begin()) - 1;
  i               = std::clamp(i, 0, n - 2);
  const double h  = l_[i + 1] - l_[i];
  const double a  = l_[i + 1] - l;
  const double b  = l - l_[i];
  const double w2 = omega_ * omega_;
  const double s  = -(m_[i] * std::sin(omega_ * a) + m_[i + 1] * std::sin(omega_ * b)) /
                        (w2 * std::sin(omega_ * h)) +
                    (u_[i] + m_[i] / w2) * a / h + (u_[i + 1] + m_[i + 1] / w2) * b / h;
  return Envelope(l) * s / (l * (l + 1.));
}
