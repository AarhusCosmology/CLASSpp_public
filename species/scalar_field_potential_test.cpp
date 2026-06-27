#include "scalar_field_potential.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
  const std::vector<double> p    = {1.22, 2.0, 0.5, 0.1};  // lambda, alpha, A, B
  const double phi               = 0.7;
  const ScalarFieldPotential pot = DefaultScalarFieldPotential();

  // Default bundle reproduces V = exp(-lambda*phi)*((phi-B)^alpha + A).
  const double lambda = p[0], alpha = p[1], A = p[2], B = p[3];
  const double Ve          = std::exp(-lambda * phi);
  const double V_expected  = Ve * (std::pow(phi - B, alpha) + A);
  const double dV_expected = -lambda * Ve * (std::pow(phi - B, alpha) + A) +
                             Ve * alpha * std::pow(phi - B, alpha - 1);
  assert(std::fabs(pot.V(phi, p) - V_expected) < 1e-12);
  assert(std::fabs(pot.dV(phi, p) - dV_expected) < 1e-12);

  // An injected pure-1EXP bundle (params = [V0, lambda]) is honored.
  ScalarFieldPotential exp1{[](double f, const std::vector<double>& q) {
                              return q[0] * std::exp(-q[1] * f);
                            },
                            [](double f, const std::vector<double>& q) {
                              return -q[1] * q[0] * std::exp(-q[1] * f);
                            },
                            [](double f, const std::vector<double>& q) {
                              return q[1] * q[1] * q[0] * std::exp(-q[1] * f);
                            }};
  const std::vector<double> q = {3.0, 1.22};
  assert(std::fabs(exp1.V(phi, q) - 3.0 * std::exp(-1.22 * phi)) < 1e-12);

  std::printf("scalar_field_potential tests passed\n");
  return 0;
}
