#include "scalar_field_potential.h"

#include <cmath>
#include <utility>

ScalarFieldPotential DefaultScalarFieldPotential() {
  ScalarFieldPotential p;
  p.V = [](double phi, const std::vector<double>& params) {
    const double lambda = params[0], alpha = params[1], A = params[2], B = params[3];
    return std::exp(-lambda * phi) * (std::pow(phi - B, alpha) + A);
  };
  p.dV = [](double phi, const std::vector<double>& params) {
    const double lambda = params[0], alpha = params[1], A = params[2], B = params[3];
    const double Ve  = std::exp(-lambda * phi);
    const double Vp  = std::pow(phi - B, alpha) + A;
    const double dVp = alpha * std::pow(phi - B, alpha - 1);
    return -lambda * Ve * Vp + Ve * dVp;
  };
  p.ddV = [](double phi, const std::vector<double>& params) {
    const double lambda = params[0], alpha = params[1], A = params[2], B = params[3];
    const double Ve   = std::exp(-lambda * phi);
    const double Vp   = std::pow(phi - B, alpha) + A;
    const double dVe  = -lambda * Ve;
    const double dVp  = alpha * std::pow(phi - B, alpha - 1);
    const double ddVe = lambda * lambda * Ve;
    const double ddVp = alpha * (alpha - 1.) * std::pow(phi - B, alpha - 2);
    return ddVe * Vp + 2 * dVe * dVp + Ve * ddVp;
  };
  // Reproduces the historical ComputeShootingGuess logic exactly (H0 unused):
  //   tuning_index == 0 -> {sqrt(3/Omega), -0.5*sqrt(3)*Omega^-1.5}
  //   otherwise         -> {params[tuning_index], 1}
  p.shooting_guess = [](double omega0_scf,
                        double /*H0*/,
                        const std::vector<double>& params,
                        int tuning_index) -> std::pair<double, double> {
    if (tuning_index == 0)
      return {std::sqrt(3.0 / omega0_scf), -0.5 * std::sqrt(3.0) * std::pow(omega0_scf, -1.5)};
    return {params[tuning_index], 1.0};
  };
  return p;
}

ScalarFieldPotential AxionScalarFieldPotential() {
  ScalarFieldPotential p;
  p.V = [](double phi, const std::vector<double>& q) {
    const double m = q[0], f = q[1], n = q[2];
    const double u = 1. - std::cos(phi / f);
    return m * m * f * f * std::pow(u, n);
  };
  p.dV = [](double phi, const std::vector<double>& q) {
    const double m = q[0], f = q[1], n = q[2];
    const double u = 1. - std::cos(phi / f);
    return m * m * f * n * std::pow(u, n - 1.) * std::sin(phi / f);
  };
  p.ddV = [](double phi, const std::vector<double>& q) {
    const double m = q[0], f = q[1], n = q[2];
    const double u = 1. - std::cos(phi / f);
    // ddV = m^2 n u^(n-1) [ (n-1)(2-u) + (1-u) ]: sin^2 = u(2-u) removes the
    // u^(n-2) singularity, so this is regular at u -> 0 for all n >= 1.
    return m * m * n * std::pow(u, n - 1.) * ((n - 1.) * (2. - u) + (1. - u));
  };
  p.shooting_guess = [](double omega,
                        double H0,
                        const std::vector<double>& q,
                        int /*tuning_index*/) -> std::pair<double, double> {
    const double f = q[1], n = q[2], theta = q[3];
    const double u0 = 1. - std::cos(theta);
    const double m  = std::sqrt(3. * H0 * H0 * omega / (f * f * std::pow(u0, n)));
    return {m, m / (2. * omega)};  // dm/dOmega for m = sqrt(c*Omega)
  };
  return p;
}
