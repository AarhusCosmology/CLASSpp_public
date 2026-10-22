#include "asymptotic_extrapolation.h"

#include <cmath>
#include <cstddef>

AsymptoticExtrapolation::AsymptoticExtrapolation(
    const double* k, const double* source, int n, int stride, int power)
    : power_(power) {
  if (n < 1) {
    return;
  }
  auto y_at = [&](int i) {
    return source[static_cast<std::ptrdiff_t>(i) * stride] * std::pow(k[i], power);
  };
  const double k_max = k[n - 1];
  ln_k_max_          = std::log(k_max);
  inv_k_max_         = 1. / k_max;
  y_max_             = y_at(n - 1);

  /* least squares for y - y_max = A u + C v, u = ln k - ln k_max, v = 1/k - 1/k_max */
  double suu = 0.;
  double suv = 0.;
  double svv = 0.;
  double suy = 0.;
  double svy = 0.;
  int points = 0;
  for (int i = n - 2; (i >= 0) && (k[i] >= k_max / kFitRange); --i) {
    const double u  = std::log(k[i]) - ln_k_max_;
    const double v  = 1. / k[i] - inv_k_max_;
    const double y  = y_at(i) - y_max_;
    suu            += u * u;
    suv            += u * v;
    svv            += v * v;
    suy            += u * y;
    svy            += v * y;
    ++points;
  }
  const double det = suu * svv - suv * suv;
  if ((points >= 2) && (std::abs(det) > 1e-12 * suu * svv)) {
    slope_ln_k_  = (suy * svv - svy * suv) / det;
    slope_inv_k_ = (svy * suu - suy * suv) / det;
  }
  else if (points >= 1) {
    slope_ln_k_ = suy / suu;
  }
}

double AsymptoticExtrapolation::operator()(double k) const {
  const double y = y_max_ + slope_ln_k_ * (std::log(k) - ln_k_max_) +
                   slope_inv_k_ * (1. / k - inv_k_max_);
  return y / std::pow(k, power_);
}
