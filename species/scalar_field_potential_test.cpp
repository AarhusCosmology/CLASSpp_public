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

  // ── Axion bundle: V = m^2 f^2 (1 - cos(phi/f))^n, params = [m, f, n, Theta_ini]
  {
    const ScalarFieldPotential ax = AxionScalarFieldPotential();
    for (double n : {1.0, 2.0, 2.5, 3.0}) {
      const std::vector<double> ap = {2.0e-3, 0.5, n, 2.0};  // m, f, n, Theta
      const double m = ap[0], f = ap[1];
      // V matches the closed form away from the minimum.
      for (double ph : {0.05, 0.4, 1.0, 1.4}) {
        const double u     = 1.0 - std::cos(ph / f);
        const double V_exp = m * m * f * f * std::pow(u, n);
        assert(std::fabs(ax.V(ph, ap) - V_exp) < 1e-12 * std::fabs(V_exp) + 1e-300);
        // dV, ddV against central finite differences.
        const double h      = 1e-6;
        const double dV_fd  = (ax.V(ph + h, ap) - ax.V(ph - h, ap)) / (2 * h);
        const double ddV_fd = (ax.V(ph + h, ap) - 2 * ax.V(ph, ap) + ax.V(ph - h, ap)) / (h * h);
        assert(std::fabs(ax.dV(ph, ap) - dV_fd) < 1e-5 * std::fabs(dV_fd) + 1e-12);
        assert(std::fabs(ax.ddV(ph, ap) - ddV_fd) < 1e-3 * std::fabs(ddV_fd) + 1e-8);
      }
      // Regularity at the minimum: no NaN/inf, and for n=1 ddV(0) = m^2 exactly.
      assert(std::isfinite(ax.V(0.0, ap)));
      assert(std::isfinite(ax.dV(0.0, ap)));
      assert(std::isfinite(ax.ddV(0.0, ap)));
    }
    const std::vector<double> a1 = {2.0e-3, 0.5, 1.0, 2.0};
    assert(std::fabs(ax.ddV(0.0, a1) - a1[0] * a1[0]) < 1e-15);

    // Frozen-field shooting guess: V(Theta*f) = 3 H0^2 Omega at the guessed m.
    const double H0 = 2.2e-4, omega = 0.05;
    const auto [mg, dmdo]  = ax.shooting_guess(omega, H0, a1, 0);
    std::vector<double> ag = a1;
    ag[0]                  = mg;
    assert(std::fabs(ax.V(a1[3] * a1[1], ag) - 3 * H0 * H0 * omega) < 1e-10 * 3 * H0 * H0 * omega);
    assert(std::fabs(dmdo - mg / (2 * omega)) < 1e-12 * mg / omega);
  }

  std::printf("scalar_field_potential tests passed\n");
  return 0;
}
