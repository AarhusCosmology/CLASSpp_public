// The Wigner d^l_{mm'}(mu) of the lensing module, computed for kChunk nodes at once and stored
// [l][node], against exact values: d^l_{mm'}(0) = delta_{mm'}, d^l_{mm'}(pi) = (-1)^(l-m)
// delta_{m,-m'}, and d^l_00(mu) = P_l(mu). The nodes sit at different positions of the chunk, so
// a slip in the indexing of the table shows up as well as one in a recurrence.
// docs/superpowers/specs/2026-10-26-lensing-wigner-recursions-design.md
#include "wigner_d.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr int kChunk = wigner::kChunk;
constexpr int kLmax  = 3000;

using Recurrence = void (*)(const std::array<double, kChunk>&, int, double*);

struct Case {
  const char* name;
  Recurrence d;
  int m;
  int mp;
};

double At(const std::vector<double>& d, int l, int j) {
  return d[static_cast<size_t>(l) * kChunk + j];
}

}  // namespace

int main() {
  /* node 3 at mu = 1, node 20 at mu = -1, the others spread in between */
  std::array<double, kChunk> mu{};
  for (int j = 0; j < kChunk; j++) {
    mu[j] = cos(0.05 + 3.0 * j / kChunk);
  }
  mu[3]  = 1.;
  mu[20] = -1.;

  const std::array<Case, 12> cases{{{"d00", wigner::D00, 0, 0},
                                    {"d11", wigner::D11, 1, 1},
                                    {"d1m1", wigner::D1m1, 1, -1},
                                    {"d2m2", wigner::D2m2, 2, -2},
                                    {"d22", wigner::D22, 2, 2},
                                    {"d20", wigner::D20, 2, 0},
                                    {"d31", wigner::D31, 3, 1},
                                    {"d3m1", wigner::D3m1, 3, -1},
                                    {"d3m3", wigner::D3m3, 3, -3},
                                    {"d40", wigner::D40, 4, 0},
                                    {"d4m2", wigner::D4m2, 4, -2},
                                    {"d4m4", wigner::D4m4, 4, -4}}};

  std::vector<double> d(static_cast<size_t>(kLmax + 1) * kChunk);
  for (const Case& c : cases) {
    c.d(mu, kLmax, d.data());
    double worst = 0.;
    for (int l = std::max(std::abs(c.m), std::abs(c.mp)); l <= kLmax; l++) {
      const double at_0  = (c.m == c.mp) ? 1. : 0.;
      const double at_pi = (c.m == -c.mp) ? (((l - c.m) % 2 == 0) ? 1. : -1.) : 0.;
      worst              = std::max(worst, std::fabs(At(d, l, 3) - at_0));
      worst              = std::max(worst, std::fabs(At(d, l, 20) - at_pi));
    }
    if (worst > 1e-9) {
      std::printf("%s: |d - exact| = %.2e at mu = +-1\n", c.name, worst);
      assert(false);
    }
  }

  /* d^l_00 = P_l at every node */
  wigner::D00(mu, kLmax, d.data());
  for (int j = 0; j < kChunk; j++) {
    double p_lm1 = 1.;
    double p_l   = mu[j];
    for (int l = 1; l < kLmax; l++) {
      const double p_lp1 = ((2 * l + 1) * mu[j] * p_l - l * p_lm1) / (l + 1);
      p_lm1              = p_l;
      p_l                = p_lp1;
      assert(std::fabs(At(d, l + 1, j) - p_l) < 1e-9);
    }
  }

  std::printf("lensing Wigner d-function tests passed\n");
  return 0;
}
