#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "background.h"
#include "scalar_field.h"

int main() {
  background pba{};
  const double beta                = -2.0;
  const std::vector<double> params = {1.22, 0.0, 0.0, 0.0};  // pure exp: V = exp(-1.22 phi)
  ScalarFieldSpecies scf(pba,
                         /*omega0=*/0.7,
                         params,
                         /*tuning=*/0,
                         /*attractor=*/true,
                         /*phi_ini=*/1.,
                         /*phi_prime_ini=*/1.,
                         DefaultScalarFieldPotential(),
                         beta);
  assert(std::fabs(scf.beta() - beta) < 1e-15);

  // Register the background/integration indices into fresh counters.
  int index_bg = 0, index_bi = 0;
  scf.RegisterBackgroundIndices(index_bg);
  scf.RegisterIntegrationIndices(index_bi);
  std::vector<double> pvecback(index_bg, 0.), pvecback_B(index_bi, 0.);

  const double a = 0.5, phi = 0.3, phi_prime = 0.4;
  pvecback_B[scf.bi_phi_index()]       = phi;
  pvecback_B[scf.bi_phi_prime_index()] = phi_prime;
  scf.ComputeBackground(a, pvecback_B.data(), pvecback.data());

  const double V            = std::exp(-1.22 * phi);
  const double rho_expected = ((1. - 2. * beta) * phi_prime * phi_prime / (2. * a * a) + V) / 3.;
  // rho is stored at the scalar field's background rho slot; Rho() reads it.
  assert(std::fabs(scf.Rho(pvecback.data()) - rho_expected) < 1e-12);

  std::printf("scalar_field beta background test passed\n");
  return 0;
}
