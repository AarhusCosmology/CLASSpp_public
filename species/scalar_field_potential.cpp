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
