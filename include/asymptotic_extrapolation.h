#pragma once

/**
 * Continuation of a perturbation source past the last wavenumber of its table.
 *
 * Deep inside the horizon in matter domination, the density contrast grows as ln k with a
 * 1/k correction, and the potentials follow it through the Poisson equation
 * (k^2 phi ~ delta). So s(k) k^power = A ln k + B + C/k, with power 0 for a density and 2
 * for a potential. A and C are fitted by least squares to the table points with k in
 * [k_max/3, k_max), and B so that the continuation passes through the last point: it joins
 * the table without a step. Fitted at k_max = 1.8/Mpc, it reproduces solved modes to 0.2% up
 * to 40/Mpc (docs/superpowers/specs/2026-10-21-lensing-source-extrapolation-design.md).
 *
 * With one point in the fitting range it fits the ln k slope only; with none it is constant.
 */
class AsymptoticExtrapolation {
 public:
  AsymptoticExtrapolation() = default;

  /**
   * @param k       increasing wavenumbers of the table, k[0..n-1]
   * @param source  the source at k[i] is source[i * stride]
   * @param n       number of table points
   * @param stride  distance between consecutive k in source
   * @param power   0 for a density, 2 for a potential
   */
  AsymptoticExtrapolation(const double* k, const double* source, int n, int stride, int power);

  /** The continued source at k (meant for k >= the last k of the table). */
  [[nodiscard]] double operator()(double k) const;

 private:
  double ln_k_max_    = 0.;
  double inv_k_max_   = 0.;
  double y_max_       = 0.;
  double slope_ln_k_  = 0.;
  double slope_inv_k_ = 0.;
  int power_          = 0;
};
