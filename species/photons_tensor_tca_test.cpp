/* The leading-order tensor tight-coupling closure, checked against the hierarchy
   it closes rather than against transcribed numbers.

   PhotonsSpecies::PerturbTensorDerivs integrates

     F_0' = -4/3 theta - kappa' (F_0 + sqrt6 P_2) + sqrt6 hdot
     G_0' = -k G_1     - kappa' (G_0 - sqrt6 P_2)

   With F_0, G_0 = O(1/kappa') the O(1) balance as kappa' -> infinity is

     kappa' (F_0 + sqrt6 P_2) = sqrt6 hdot        (temperature, against the GW source)
     kappa' (G_0 - sqrt6 P_2) = 0                 (polarisation)

   and P_2 is the hierarchy's own multipole combination evaluated on the closed
   state. Those three statements have a unique root, so a closure that drifts
   from the hierarchy -- or a call site that grows its own constants -- fails
   here. Before #398 the three sites disagreed by a factor of 7.2. */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>

#include "constants.h"
#include "photons.h"

namespace {

/* Arbitrary, nonzero, and not 1: a dropped or spurious factor of gwdot or
   dkappa has to show up rather than cancel. */
constexpr double kGwdot  = 3.7e-5;
constexpr double kDkappa = 2.4e3;

bool Close(double a, double b) {
  return std::fabs(a - b) <= 1e-12 * std::max({1.0, std::fabs(a), std::fabs(b)});
}

/* Vanishes against the natural scale of the balance it came from, sqrt6*hdot --
   comparing a residual to a bare 0 would only test the tolerance. */
bool Vanishes(double residual) {
  return std::fabs(residual) <= 1e-12 * _SQRT6_ * std::fabs(kGwdot);
}

/* The temperature equation's stiff term balances the gravitational-wave source. */
void test_closure_balances_the_temperature_equation() {
  const PhotonsSpecies::TensorTcaClosure c = PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa);
  assert(Vanishes(kDkappa * (c.F0 + _SQRT6_ * c.P2) - _SQRT6_ * kGwdot));
}

/* The polarisation equation has no source, so its stiff term balances to zero. */
void test_closure_balances_the_polarisation_equation() {
  const PhotonsSpecies::TensorTcaClosure c = PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa);
  assert(Vanishes(kDkappa * (c.G0 - _SQRT6_ * c.P2)));
}

/* P_2 must be the hierarchy's own combination evaluated on the closed state,
   not an independently written constant. This is the assertion the three call
   sites failed: substituting the closure into the combination gave +0.13608
   gwdot/kappa' while the source branch used +0.97980. */
void test_P2_is_the_hierarchy_combination_of_the_closed_state() {
  const PhotonsSpecies::TensorTcaClosure c = PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa);
  const double from_multipoles             = PhotonsSpecies::TensorP2(c.F0, 0., 0., c.G0, 0., 0.);
  assert(Close(from_multipoles, c.P2));
}

/* The polarisation-to-temperature ratio implied by the balance. The old closure
   satisfied this one, which is why the error survived: it is necessary, not
   sufficient. It is a leading-order statement only -- see the first-order tests. */
void test_polarisation_monopole_is_minus_a_quarter_of_the_temperature_monopole() {
  const PhotonsSpecies::TensorTcaClosure c = PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa);
  assert(Close(c.G0 / c.F0, -0.25));
}

/* Both old branches had P_2 positive for positive hdot. The sign is the whole
   difference between adding to and subtracting from the B-mode source. */
void test_P2_opposes_the_gravitational_wave_source() {
  const PhotonsSpecies::TensorTcaClosure c = PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa);
  assert(c.P2 < 0.);
  assert(c.F0 > 0.);
}

/* Every component scales as hdot/kappa' and nothing else -- a stray additive
   term or a wrong power of kappa' shows up here and nowhere else. */
void test_closure_scales_as_gwdot_over_dkappa() {
  const PhotonsSpecies::TensorTcaClosure c = PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa);
  const PhotonsSpecies::TensorTcaClosure scaled = PhotonsSpecies::TensorTightCoupling(3. * kGwdot,
                                                                                      2. * kDkappa);
  assert(Close(scaled.F0, 1.5 * c.F0));
  assert(Close(scaled.G0, 1.5 * c.G0));
  assert(Close(scaled.P2, 1.5 * c.P2));
}

/* The combination is the one written in PerturbTensorDerivs. Pinning one
   off-diagonal coefficient keeps a future edit from quietly renormalising it. */
void test_TensorP2_uses_the_hierarchy_coefficients() {
  /* Only the shear slot populated: P_2 = -(2/7) sigma / sqrt6. */
  const double sigma = 1.3e-4;
  assert(Close(PhotonsSpecies::TensorP2(0., sigma, 0., 0., 0., 0.), -(2. / 7.) * sigma / _SQRT6_));
  /* Only the polarisation quadrupole: P_2 = -(6/7) G_2 / sqrt6. */
  const double g2 = -8.1e-5;
  assert(Close(PhotonsSpecies::TensorP2(0., 0., 0., 0., g2, 0.), -(6. / 7.) * g2 / _SQRT6_));
}

/* ---------------------------------------------------------------------------
   First order in tau_c. F_0 and G_0 vary on the timescale over which tau_c
   itself varies, and through recombination that is not slow: the dropped terms
   are kappa'^-1 F_0' and kappa'^-1 G_0', both O(tau_c) relative to the leading
   closure, whereas the shear and G_2 the multipole combination also omits are
   O(tau_c^2) relative. Keeping the derivative terms, the balance becomes

     kappa' (F_0 + sqrt6 P_2) = sqrt6 hdot - F_0'
     kappa' (G_0 - sqrt6 P_2) = -G_0'

   evaluated on the leading solution, F_0' = (4 sqrt6/3) D and G_0' = -(sqrt6/3) D
   with D = d/dtau (hdot/kappa'). --------------------------------------------- */

constexpr double kDrift = 5.1e-2; /* D, in the same units as gwdot */

/* Zero drift must reproduce the leading closure exactly -- the default argument
   is what every call site that has no derivative to offer relies on. */
void test_zero_drift_reproduces_the_leading_closure() {
  const PhotonsSpecies::TensorTcaClosure lead = PhotonsSpecies::TensorTightCoupling(kGwdot,
                                                                                    kDkappa);
  const PhotonsSpecies::TensorTcaClosure first =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, 0.);
  assert(Close(first.F0, lead.F0));
  assert(Close(first.G0, lead.G0));
  assert(Close(first.P2, lead.P2));
}

void test_first_order_balances_the_temperature_equation() {
  const PhotonsSpecies::TensorTcaClosure c =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, kDrift);
  const double dF0 = 4. * _SQRT6_ / 3. * kDrift; /* derivative of the leading F_0 */
  assert(Vanishes(kDkappa * (c.F0 + _SQRT6_ * c.P2) - (_SQRT6_ * kGwdot - dF0)));
}

void test_first_order_balances_the_polarisation_equation() {
  const PhotonsSpecies::TensorTcaClosure c =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, kDrift);
  const double dG0 = -_SQRT6_ / 3. * kDrift; /* derivative of the leading G_0 */
  assert(Vanishes(kDkappa * (c.G0 - _SQRT6_ * c.P2) + dG0));
}

/* P_2 stays the hierarchy's own combination at first order too. */
void test_first_order_P2_is_the_hierarchy_combination() {
  const PhotonsSpecies::TensorTcaClosure c =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, kDrift);
  assert(Close(PhotonsSpecies::TensorP2(c.F0, 0., 0., c.G0, 0., 0.), c.P2));
}

/* The -1/4 ratio is leading order only. If the first-order closure still
   satisfied it, the drift terms would have entered F_0 and G_0 in proportion --
   they do not (-22/9 against +13/9), and that asymmetry is the correction. */
void test_first_order_breaks_the_minus_one_quarter_ratio() {
  const PhotonsSpecies::TensorTcaClosure c =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, kDrift);
  assert(!Close(c.G0 / c.F0, -0.25));
}

/* The drift enters linearly, on top of the leading solution. */
void test_drift_enters_linearly() {
  const PhotonsSpecies::TensorTcaClosure lead = PhotonsSpecies::TensorTightCoupling(kGwdot,
                                                                                    kDkappa);
  const PhotonsSpecies::TensorTcaClosure a =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, kDrift);
  const PhotonsSpecies::TensorTcaClosure b =
      PhotonsSpecies::TensorTightCoupling(kGwdot, kDkappa, 2. * kDrift);
  assert(Close(b.F0 - lead.F0, 2. * (a.F0 - lead.F0)));
  assert(Close(b.G0 - lead.G0, 2. * (a.G0 - lead.G0)));
  assert(Close(b.P2 - lead.P2, 2. * (a.P2 - lead.P2)));
}

/* D itself: d/dtau of (hdot/kappa'), from the quantities a call site has. */
void test_drift_is_the_derivative_of_gwdot_over_dkappa() {
  const double gw_second = -7.3e-3, ddkappa = -1.9e5;
  const double expected = gw_second / kDkappa - kGwdot * ddkappa / (kDkappa * kDkappa);
  assert(Close(PhotonsSpecies::TensorTcaDrift(kGwdot, kDkappa, gw_second, ddkappa), expected));
}

}  // namespace

int main() {
  test_closure_balances_the_temperature_equation();
  test_closure_balances_the_polarisation_equation();
  test_P2_is_the_hierarchy_combination_of_the_closed_state();
  test_polarisation_monopole_is_minus_a_quarter_of_the_temperature_monopole();
  test_P2_opposes_the_gravitational_wave_source();
  test_closure_scales_as_gwdot_over_dkappa();
  test_TensorP2_uses_the_hierarchy_coefficients();

  test_zero_drift_reproduces_the_leading_closure();
  test_first_order_balances_the_temperature_equation();
  test_first_order_balances_the_polarisation_equation();
  test_first_order_P2_is_the_hierarchy_combination();
  test_first_order_breaks_the_minus_one_quarter_ratio();
  test_drift_enters_linearly();
  test_drift_is_the_derivative_of_gwdot_over_dkappa();

  std::printf("tensor tight-coupling closure tests passed\n");
  return 0;
}
