#include "energy_deposition.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <initializer_list>

namespace {

/* The Galli et al. 2013 table nodes we assert against, straight from
   arXiv:1306.0563 Table 1 as distributed with class_public
   (external/heating/Galli_et_al_2013.dat): x, heat, lya, ion_H, ion_He, lowE. */
struct Node {
  double x, heat, lya, ion_H, ion_He, lowE;
};

// clang-format off
const Node kNodes[] = {
    {0.000000, 0.151800, 0.323526, 0.350798, 0.024367, 0.142202},
    {0.001000, 0.210027, 0.291280, 0.341822, 0.023134, 0.128324},
    {0.010000, 0.338316, 0.238304, 0.301893, 0.021302, 0.103617},
    {0.100000, 0.654816, 0.121705, 0.175739, 0.012852, 0.053388},
    {0.500000, 0.923644, 0.026168, 0.043901, 0.004914, 0.011310},
    {0.990000, 0.995299, 0.000290, 0.001700, 0.002924, 0.000118},
};
// clang-format on

bool close(double a, double b, double tol) {
  return std::fabs(a - b) <= tol;
}

void test_galli_reproduces_its_own_nodes() {
  for (const Node& n : kNodes) {
    EnergyDeposition dep = energy_deposition_fractions(deposition_Galli_2013, n.x);
    assert(close(dep.heat, n.heat, 1e-12));
    assert(close(dep.lya, n.lya, 1e-12));
    assert(close(dep.ion_H, n.ion_H, 1e-12));
    assert(close(dep.ion_He, n.ion_He, 1e-12));
    assert(close(dep.lowE, n.lowE, 1e-12));
  }
}

/* Every channel is a fraction of the deposited energy, so each lies in [0,1].
   Their sum is close to but not exactly 1: the published rows are computed per
   channel and sum to between 0.993 and 1.019 across the table, so this bounds the
   interpolant to that range rather than enforcing energy conservation. The
   interesting part is the lower bound -- a spline across the shoulder near x = 1
   dips below zero, and the implementation must not let that through. */
void test_galli_fractions_stay_physical() {
  for (int i = 0; i <= 2000; i++) {
    double x             = 1.2 * i / 2000.;
    EnergyDeposition dep = energy_deposition_fractions(deposition_Galli_2013, x);
    double channels[]    = {dep.heat, dep.ion_H, dep.ion_He, dep.lya, dep.lowE};
    double sum           = 0.;
    for (double c : channels) {
      assert(c >= 0.);
      assert(c <= 1.);
      sum += c;
    }
    assert(sum >= 0.99);
    assert(sum <= 1.02);
  }
}

/* Above the last tabulated ionization fraction everything is heat, and the model
   must stay defined for the x > 1 values reached while helium is still ionized
   (x = x_H + f_He x_He can reach 1 + 2 f_He). */
void test_galli_saturates_when_fully_ionized() {
  for (double x : {1.2, 1.5, 2.0}) {
    EnergyDeposition dep = energy_deposition_fractions(deposition_Galli_2013, x);
    assert(close(dep.heat, 1., 1e-12));
    assert(close(dep.ion_H, 0., 1e-12));
    assert(close(dep.ion_He, 0., 1e-12));
    assert(close(dep.lya, 0., 1e-12));
    assert(close(dep.lowE, 0., 1e-12));
  }
}

/* The legacy model must reproduce the two analytic fits the thermodynamics module
   carried inline before this change -- including the fact that the Lyman-alpha
   channel reused the ionization coefficient, which is what RECFAST's
   (1-C)/L_alpha term assumed. Compared at a tight relative tolerance rather than
   exactly: classpp is built with -ffast-math and this test executable is not, so
   std::pow can differ in the last ULP across that boundary. */
void test_legacy_reproduces_the_inline_fits() {
  for (int i = 1; i < 100; i++) {
    double x             = i / 100.;
    EnergyDeposition dep = energy_deposition_fractions(deposition_legacy, x);

    double expected_heat =
        std::fmin(0.996857 * (1. - std::pow(1. - std::pow(x, 0.300134), 1.51035)), 1.0);
    double expected_ion = 0.369202 * std::pow(1. - std::pow(x, 0.463929), 1.70237);

    assert(close(dep.heat, expected_heat, 1e-12 * expected_heat));
    assert(close(dep.ion_H, expected_ion, 1e-12 * expected_ion));
    assert(close(dep.lya, expected_ion, 1e-12 * expected_ion));
    assert(dep.ion_He == 0.);
    assert(dep.lowE == 0.);
  }
}

void test_legacy_saturates_when_fully_ionized() {
  for (double x : {1.0, 1.1}) {
    EnergyDeposition dep = energy_deposition_fractions(deposition_legacy, x);
    assert(dep.heat == 1.);
    assert(dep.ion_H == 0.);
    assert(dep.lya == 0.);
    assert(dep.ion_He == 0.);
  }
}

/* The two models disagree, and the disagreement is the reason for the change:
   the legacy fit sends nothing into heat at vanishing ionization, while the
   published table sends 15%. Pin the gap so a silent swap of one for the other
   cannot pass unnoticed. */
void test_the_models_actually_differ_at_low_ionization() {
  EnergyDeposition table = energy_deposition_fractions(deposition_Galli_2013, 1e-4);
  EnergyDeposition fit   = energy_deposition_fractions(deposition_legacy, 1e-4);

  assert(close(table.heat, 0.15188, 1e-5));
  assert(table.heat > 1.5 * fit.heat);
  assert(table.lya > 0.3);
  assert(close(fit.lya, fit.ion_H, 1e-12));
}

}  // namespace

int main() {
  test_galli_reproduces_its_own_nodes();
  test_galli_fractions_stay_physical();
  test_galli_saturates_when_fully_ionized();
  test_legacy_reproduces_the_inline_fits();
  test_legacy_saturates_when_fully_ionized();
  test_the_models_actually_differ_at_low_ionization();

  std::printf("energy deposition tests passed\n");
  return 0;
}
