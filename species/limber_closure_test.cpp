/**
 * Unit tests for the Limber closure of a continuously-sourced hierarchy, shared
 * by DarkRadiationSpecies (dcdm_dr, dncdm_dr) and WdmDecayProductSpecies
 * (dcdm_wdm's massive daughter).
 *
 * See docs/superpowers/specs/2026-09-07-dr-limber-closure-design.md for the
 * derivation and .../2026-09-09-massive-daughter-limber-closure-design.md for
 * the massive case. These tests pin the two numbers the closure rests on:
 *
 *   - limber_closure::Order(l): the exact flat moment ratio M_0 / M_{-1} of j_l,
 *     which is where the Bessel weight in the shell superposition sits.
 *   - limber_closure::Cot(L, k, K): cot_K evaluated at the Limber distance for
 *     that order, in the units of perturb_workspace::cotKgen.
 */

#include "limber_closure.h"

#include <cassert>
#include <cmath>
#include <cstdio>

namespace {

// x_eff(l) = 2(l+1)/l [Gamma((l+1)/2)/Gamma(l/2)]^2, computed independently in
// python3 from math.lgamma. These are the values tabulated in issue #418.
struct OrderCase {
  int l;
  double x_eff;
};
constexpr OrderCase kOrderTable[] = {
    {5, 5.432488724203},
    {10, 10.464385502859},
    {17, 17.478593179881},
    {25, 25.485301167508},
    {40, 40.490742503069},
    {100, 100.496268772145},
};

void test_order_matches_closed_form() {
  for (const auto& c : kOrderTable) {
    const double L = limber_closure::Order(c.l);
    assert(std::abs(L - c.x_eff) / c.x_eff < 1.e-12);
  }
}

/* The whole point of the closure is that L is Limber's l+1/2, not MB's k*tau.
   Pin the approach to it: L = l + 1/2 - 3/(8l) + 3/(16 l^2) + O(l^-3), so the
   residual after removing the 1/l term must scale as 1/l^2 with coefficient
   3/16 -- a plain "L ~ l+1/2" assertion would pass for the wrong function. */
void test_order_asymptote() {
  for (const auto& c : kOrderTable) {
    const double l     = c.l;
    const double L     = limber_closure::Order(c.l);
    const double resid = L - (l + 0.5 - 3. / (8. * l));
    assert(resid > 0.);  // approaches l+1/2 from below
    assert(std::abs(resid * l * l - 3. / 16.) < 0.5 / l);
  }
  // ... and the closure coefficient (2l+1)/L is therefore ~2: linear
  // extrapolation of the ladder, versus MB's ~ -F_{L-1}.
  for (const auto& c : kOrderTable) {
    const double cl = (2. * c.l + 1.) / limber_closure::Order(c.l);
    assert(cl > 2.);
    assert(cl < 2. + 1. / (static_cast<double>(c.l) * c.l));
  }
}

void test_flat_cot_is_reciprocal_order() {
  const double L = limber_closure::Order(17);
  for (double k : {1.e-4, 1.e-2, 0.1, 1., 10.})
    assert(limber_closure::Cot(L, k, 0.) == 1. / L);
}

/* Independent evaluation: solve the curved-Limber turning point the way
   TransferModule::transfer_limber does (asin / asinh on the hyperspherical
   argument), then evaluate cot_K there directly. The closed form under test
   must reproduce it. */
double CotAtTurningPointClosed(double L, double k, double K) {
  const double sqrtK = std::sqrt(K);
  const double q     = std::sqrt(k * k + K);
  const double chi   = std::asin(sqrtK * L / q) / sqrtK;
  return sqrtK / k / std::tan(sqrtK * chi);
}

double CotAtTurningPointOpen(double L, double k, double K) {
  const double sqrtK = std::sqrt(-K);
  const double q     = std::sqrt(k * k + K);
  const double chi   = std::asinh(sqrtK * L / q) / sqrtK;
  return sqrtK / k / std::tanh(sqrtK * chi);
}

void test_curved_cot_matches_turning_point() {
  for (int l : {4, 17, 40}) {
    const double L = limber_closure::Order(l);

    for (double K : {1.e-8, 1.e-7, 1.e-6}) {
      const double k   = 0.1;
      const double ref = CotAtTurningPointClosed(L, k, K);
      const double got = limber_closure::Cot(L, k, K);
      assert(std::abs(got - ref) / std::abs(ref) < 1.e-12);
    }

    for (double K : {-1.e-8, -1.e-6, -1.e-4}) {
      const double k   = 0.1;
      const double ref = CotAtTurningPointOpen(L, k, K);
      const double got = limber_closure::Cot(L, k, K);
      assert(std::abs(got - ref) / std::abs(ref) < 1.e-12);
    }
  }
}

/* Curvature enters as a correction of order K L^2 / k^2, with the sign that
   makes an open universe free-stream faster. Both must vanish smoothly as
   K -> 0, or the closure would jump at the flat limit. */
void test_curvature_sign_and_continuity() {
  const double L    = limber_closure::Order(17);
  const double k    = 0.1;
  const double flat = 1. / L;

  assert(limber_closure::Cot(L, k, 1.e-6) < flat);
  assert(limber_closure::Cot(L, k, -1.e-6) > flat);

  for (double K : {1.e-14, -1.e-14})
    assert(std::abs(limber_closure::Cot(L, k, K) - flat) / flat < 1.e-9);
}

/* In a closed universe a mode carries no multipole above nu = q/sqrt(K): the
   turning point does not exist and s_l clamps to zero there. The closure must
   clamp the same way rather than take the square root of a negative number. */
void test_closed_space_clamps_above_nu() {
  const double L = limber_closure::Order(17);
  const double k = 0.1, K = 1.e-4;  // sqrt(K) L / q = 1.74 > 1: no turning point
  assert(std::isfinite(limber_closure::Cot(L, k, K)));
  assert(limber_closure::Cot(L, k, K) == 0.);
}

/* The massive-daughter case (WdmDecayProductSpecies). Both closures are written
   as dF_L/dtau = (qk/eps) F_{L-1} - D F_L, and inverting that against the exact
   hierarchy gives the recurrence each one asserts,

       F_{L+1} = -F_{L-1} + D F_L (2L+1) / [(L+1)(qk/eps)],

   so the coefficient of F_L is (2L+1) C with C = D eps/((L+1) q k) ... i.e. D
   carries the transport rate and C does not. MB takes C = eps/(q k tau), whose
   eps/q cancels the transport rate and leaves D = (1+L) k cotKgen -- no mass
   dependence at all. The Limber argument carries no eps/q, so the transport
   rate SURVIVES: D = (1+L)(qk/eps) Cot(L). This asserts exactly that, by
   reconstructing F_{L+1} from each D and checking the coefficient of F_L. */
double ImpliedNextMultipoleCoefficient(double D, int L, double qk_div_epsilon) {
  return D * (2. * L + 1.) / ((L + 1.) * qk_div_epsilon);
}

void test_massive_daughter_damping_carries_the_transport_rate() {
  const int L    = 10;
  const double k = 0.1, tau = 1.4e4;
  const double cotKgen = 1. / (k * tau);
  const double Lord    = limber_closure::Order(L);

  // A daughter at 2% of c: eps/q = 50, which is exactly the factor at stake.
  for (double q_over_eps : {1.0, 0.02, 0.001}) {
    const double qk_div_epsilon = k * q_over_eps;

    // MB: D = (1+L) k cotKgen. Its implied coefficient must be the reciprocal of
    // the distance THIS momentum has streamed, eps/(q k tau) -- so it scales with
    // 1/q_over_eps and MB is not mass-blind after all, it is chi-blind.
    const double d_mb = (1. + L) * k * cotKgen;
    assert(std::abs(ImpliedNextMultipoleCoefficient(d_mb, L, qk_div_epsilon) -
                    (2. * L + 1.) * cotKgen / q_over_eps) < 1.e-12);

    // Limber: D = (1+L)(qk/eps) Cot(L). Its implied coefficient is (2L+1)/L_ord
    // for every momentum -- the Limber argument is a multipole order, not a
    // distance, so it must NOT depend on q/eps.
    const double d_lim = (1. + L) * qk_div_epsilon * limber_closure::Cot(Lord, k, 0.);
    assert(std::abs(ImpliedNextMultipoleCoefficient(d_lim, L, qk_div_epsilon) -
                    (2. * L + 1.) / Lord) < 1.e-12);

    // And the two differ by exactly the ratio of the arguments: dropping the
    // transport rate (D = (1+L) k Cot) puts the q/eps back where Limber has none,
    // which is the measured no-op control in the spec.
    const double d_naive = (1. + L) * k * limber_closure::Cot(Lord, k, 0.);
    assert(std::abs(ImpliedNextMultipoleCoefficient(d_naive, L, qk_div_epsilon) -
                    (2. * L + 1.) / (Lord * q_over_eps)) < 1.e-12);
  }
}

}  // namespace

int main() {
  test_order_matches_closed_form();
  test_order_asymptote();
  test_flat_cot_is_reciprocal_order();
  test_curved_cot_matches_turning_point();
  test_curvature_sign_and_continuity();
  test_closed_space_clamps_above_nu();
  test_massive_daughter_damping_carries_the_transport_rate();
  printf("limber_closure_test: all tests passed\n");
  return 0;
}
