#include <cassert>
#include <cmath>
#include <cstdio>

#include "common.h"

int main() {
  const double h = 0.67556, T = 2.7255;

  // Round-trip identity in both directions.
  const double og = Omega0gFromTcmb(T, h);
  assert(std::fabs(TcmbFromOmega0g(og, h) - T) < 1e-10);
  const double og2 = 6.0e-5;
  assert(std::fabs(Omega0gFromTcmb(TcmbFromOmega0g(og2, h), h) - og2) < 1e-15);

  // Sanity: photon Omega0 at default cosmology is ~5e-5.
  assert(og > 4.0e-5 && og < 7.0e-5);

  std::printf("photons formula tests passed (Omega0_g=%.6e)\n", og);
  return 0;
}
