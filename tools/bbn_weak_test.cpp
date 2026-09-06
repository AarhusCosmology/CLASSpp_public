/** The n <-> p weak rates.
 *
 *  Two properties pin these down almost completely, and both are exact rather
 *  than approximate, which makes them strong tests:
 *
 *    - as T -> 0 every occupation vanishes, only the decay channel survives, and
 *      the rate must equal 1/tau_n by construction;
 *    - at a COMMON temperature the three channel pairs each carry a factor
 *      exp(-Delta m / T) between them, so the total ratio is exactly that,
 *      independent of the phase-space integrals, the Coulomb factor, or the
 *      quadrature.
 *
 *  The second is what catches a swapped blocking factor or a channel integrated
 *  over the wrong range: those change the integrals but would have to conspire
 *  across all three pairs to preserve the ratio.
 */

#include "bbn_weak.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "bbn_nuclides.h"
#include "bbn_rates.h"

namespace {

bool Close(double a, double b, double rtol) {
  return std::fabs(a - b) <= rtol * std::fabs(b);
}

constexpr double kTauN = 878.4;

void test_free_decay_normalization() {
  const BbnWeakRates weak(kTauN);
  double np = 0.0, pn = 0.0;
  /* Cold enough that every occupation has underflowed away. */
  weak.Rates(MeVFromT9(1.0e-4), MeVFromT9(1.0e-4), &np, &pn);
  assert(Close(np, 1.0 / kTauN, 1.0e-10));
  assert(pn < 1.0e-30 * np);

  /* The normalization must track tau_n linearly. */
  const BbnWeakRates slow(2.0 * kTauN);
  double np2 = 0.0, pn2 = 0.0;
  slow.Rates(MeVFromT9(1.0e-4), MeVFromT9(1.0e-4), &np2, &pn2);
  assert(Close(np2, 0.5 * np, 1.0e-10));
}

/* Detailed balance at a common temperature. Exact, and independent of every
   detail of the phase-space integrals. */
void test_detailed_balance_at_common_temperature() {
  const BbnWeakRates weak(kTauN);
  for (double T9 : {100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1}) {
    const double T = MeVFromT9(T9);
    double np = 0.0, pn = 0.0;
    weak.Rates(T, T, &np, &pn);
    assert(Close(pn / np, std::exp(-kDeltaMnpMeV / T), 1.0e-12));
  }
}

/* The two temperatures must be genuinely separate. After e+- annihilation the
   neutrinos are colder than the photons, and conflating them is a silent bug
   that leaves every single-temperature test above still passing. */
void test_two_temperatures_are_distinct() {
  const BbnWeakRates weak(kTauN);
  const double T    = MeVFromT9(3.0);
  const double T_nu = T / 1.40102; /* post-annihilation ratio */

  double same_np = 0.0, same_pn = 0.0;
  double split_np = 0.0, split_pn = 0.0;
  weak.Rates(T, T, &same_np, &same_pn);
  weak.Rates(T, T_nu, &split_np, &split_pn);

  assert(!Close(split_np, same_np, 1.0e-3));
  /* Colder neutrinos mean fewer of them to capture, so both directions slow. */
  assert(split_np < same_np);
  assert(split_pn < same_pn);
  /* And with two temperatures the exact exp(-dm/T) relation no longer holds. */
  assert(!Close(split_pn / split_np, std::exp(-kDeltaMnpMeV / T), 1.0e-3));
}

/* Far above the mass difference the two directions equalise. Asserted as a
   monotone approach rather than a fixed tolerance: at T9 = 300 the ratio is
   still exp(-1.293/25.9) = 0.951, so a "within 1%" bound there would be simply
   wrong rather than demanding. */
void test_high_temperature_limit() {
  const BbnWeakRates weak(kTauN);
  double previous = 0.0;
  for (double T9 : {30.0, 100.0, 300.0, 1000.0, 3000.0}) {
    const double T = MeVFromT9(T9);
    double np = 0.0, pn = 0.0;
    weak.Rates(T, T, &np, &pn);
    const double ratio = pn / np;
    assert(ratio > previous); /* monotonically rising towards 1 */
    assert(ratio < 1.0);      /* and never crossing it: p -> n stays the uphill way */
    previous = ratio;
  }
  assert(previous > 0.995);
}

/* Freeze-out happens where the rates cross H ~ 1 s^-1, around T9 ~ 10. This is a
   coarse check that the overall scale is right, not just the ratios: a rate
   wrong by a large factor would move freeze-out and hence Y_p, while leaving
   every ratio test above intact. */
void test_rate_scale_at_freeze_out() {
  const BbnWeakRates weak(kTauN);
  double np = 0.0, pn = 0.0;
  weak.Rates(MeVFromT9(10.0), MeVFromT9(10.0), &np, &pn);
  assert(np > 0.1 && np < 10.0);
}

/* Convergence is algebraic here, not spectral -- see the note on the constructor.
   The bound is the measured 256-node figure, and the 64-node control below shows
   it is a real constraint rather than one any node count would meet. */
void test_quadrature_is_converged() {
  const BbnWeakRates coarse(kTauN, 256);
  const BbnWeakRates fine(kTauN, 2048);
  for (double T9 : {300.0, 100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1}) {
    const double T = MeVFromT9(T9);
    double cn = 0.0, cp = 0.0, fn = 0.0, fp = 0.0;
    coarse.Rates(T, T, &cn, &cp);
    fine.Rates(T, T, &fn, &fp);
    assert(Close(cn, fn, 1.0e-7));
    assert(Close(cp, fp, 1.0e-7));
  }

  const BbnWeakRates too_coarse(kTauN, 64);
  const double T = MeVFromT9(3.0);
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
  too_coarse.Rates(T, T, &a, &b);
  fine.Rates(T, T, &c, &d);
  assert(!Close(a, c, 1.0e-7));
}

/** The splined table, checked directly against the integrals it replaces along a
 *  representative trajectory (T_nu falling from T_gamma to the post-annihilation
 *  ratio). Checked at points BETWEEN the nodes, which is the only place a spline
 *  can be wrong. */
void test_weak_table_reproduces_the_integrals() {
  const BbnWeakRates exact(kTauN);
  const double T9_initial = 100.0;
  const double tau_end    = std::log(T9_initial / 0.01);
  const int nodes         = 512;

  std::vector<double> tau(nodes), T_gamma(nodes), T_nu(nodes);
  for (int i = 0; i < nodes; ++i) {
    tau[i]     = tau_end * i / (nodes - 1.0);
    T_gamma[i] = MeVFromT9(T9_initial * std::exp(-tau[i]));
    /* A stand-in for the real trajectory: T_nu/T_gamma slides from 1 to
       (4/11)^(1/3) across the annihilation, which is the behaviour the table has
       to follow. Its exact shape does not matter to this test. */
    const double x = 1.0 / (1.0 + std::exp((tau[i] - 4.5) * 2.0));
    T_nu[i]        = T_gamma[i] * (x + (1.0 - x) / 1.40102);
  }
  const BbnWeakTable table(exact, tau, T_gamma, T_nu);

  double worst_np = 0.0;
  for (int i = 0; i + 1 < nodes; ++i) {
    const double mid   = 0.5 * (tau[i] + tau[i + 1]);
    const double T     = MeVFromT9(T9_initial * std::exp(-mid));
    const double frac  = 1.0 / (1.0 + std::exp((mid - 4.5) * 2.0));
    const double T_nu_ = T * (frac + (1.0 - frac) / 1.40102);

    double enp = 0.0, epn = 0.0, tnp = 0.0, tpn = 0.0;
    exact.Rates(T, T_nu_, &enp, &epn);
    table.Rates(mid, &tnp, &tpn);
    worst_np = std::max(worst_np, std::fabs(tnp / enp - 1.0));
  }
  assert(worst_np < 1.0e-6);

  /* Outside the tabulated range the endpoints are held rather than extrapolated;
     the evolver can ask marginally past tau_end on its final step. */
  double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
  table.Rates(tau_end, &a, &b);
  table.Rates(tau_end + 1.0, &c, &d);
  assert(Close(a, c, 1.0e-12) && Close(b, d, 1.0e-12));
}

/** Nonsense construction arguments must abort, not produce plausible numbers.
 *
 *  A non-positive lifetime divides into the normalization and would hand back
 *  rates of the WRONG SIGN -- which propagates into an n/p ratio and a Y_p that
 *  look like results. A negative node count reaches std::vector as an enormous
 *  size_t. Neither is caught by any other test here, because every other test
 *  constructs the object correctly. */
void test_construction_guards() {
  for (double bad_tau : {0.0, -1.0, -878.4}) {
    bool threw = false;
    try {
      BbnWeakRates bad(bad_tau);
    }
    catch (const std::exception&) {
      threw = true;
    }
    assert(threw);
  }

  for (int bad_nodes : {-16, 0, 1, 7}) {
    bool threw = false;
    try {
      BbnWeakRates bad(kTauN, bad_nodes);
    }
    catch (const std::exception&) {
      threw = true;
    }
    assert(threw);
  }

  /* The smallest accepted node count must still construct. */
  BbnWeakRates minimal(kTauN, 8);
  double np = 0.0, pn = 0.0;
  minimal.Rates(MeVFromT9(1.0), MeVFromT9(1.0), &np, &pn);
  assert(np > 0.0 && pn > 0.0);
}

}  // namespace

int main() {
  test_free_decay_normalization();
  test_detailed_balance_at_common_temperature();
  test_two_temperatures_are_distinct();
  test_high_temperature_limit();
  test_rate_scale_at_freeze_out();
  test_quadrature_is_converged();
  test_weak_table_reproduces_the_integrals();
  test_construction_guards();
  std::printf("bbn weak-rate tests passed\n");
  return 0;
}
