/** Reaction-rate data file, its guards, and detailed balance.
 *
 *  The point of the detailed-balance check is that it is INDEPENDENT: the
 *  reverse rate the solver uses is built from per-nuclide thermal scales, while
 *  the expectation here groups the masses into a single ratio. Both come from the
 *  same physics, but a transposed index or a dropped degeneracy in one would not
 *  appear in the other.
 */

#include "bbn_rates.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <string>

namespace {

std::string RatesFile() {
  return std::string(__CLASSDIR__) + "/bbn/rates_reaclib.dat";
}
std::string TempFile(const char* name) {
  return std::string("/tmp/bbn_rates_test_") + name;
}

bool Close(double a, double b, double rtol) {
  return std::fabs(a - b) <= rtol * std::fabs(b);
}

/* (1 MeV)^3 in cm^-3, recomputed here rather than shared, so that a wrong
   conversion constant in bbn_rates.cpp cannot agree with itself. */
const double kMeV3ToInvCm3 = std::pow(1.0 / 1.973269804e-11, 3.0);

/** Reverse/forward ratio, grouping masses instead of forming a thermal scale per
 *  nuclide. Independent arithmetic, same physics. */
double ExpectedReverseOverForward(const BbnReaction& r, double T9) {
  const double T = MeVFromT9(T9);

  double ratio      = ReactionPrefactor(r.reactants, r.n_reactants) /
                      ReactionPrefactor(r.products, r.n_products);
  double mass_ratio = 1.0;
  for (int i = 0; i < r.n_reactants; ++i) {
    ratio      *= kNuclides[r.reactants[i]].g;
    mass_ratio *= NuclearMassMeV(kNuclides[r.reactants[i]]);
  }
  for (int i = 0; i < r.n_products; ++i) {
    ratio      /= kNuclides[r.products[i]].g;
    mass_ratio /= NuclearMassMeV(kNuclides[r.products[i]]);
  }

  /* One leftover factor of (m T / 2 pi)^{3/2} per unbalanced massive particle. */
  const int surplus  = r.n_reactants - r.n_products;
  ratio             *= std::pow(mass_ratio, 1.5);
  ratio             *= std::pow(T / (2.0 * M_PI), 1.5 * surplus);
  ratio             *= std::pow(kMeV3ToInvCm3, surplus);

  return ratio * std::exp(-ReactionQValueMeV(r) / T);
}

void test_network_topology_is_conserved() {
  for (int r = 0; r < kBbnNumReactions; ++r) {
    const BbnReaction& rxn = kReactions[r];
    int A_in = 0, Z_in = 0, A_out = 0, Z_out = 0;
    for (int i = 0; i < rxn.n_reactants; ++i) {
      A_in += kNuclides[rxn.reactants[i]].A;
      Z_in += kNuclides[rxn.reactants[i]].Z;
    }
    for (int i = 0; i < rxn.n_products; ++i) {
      A_out += kNuclides[rxn.products[i]].A;
      Z_out += kNuclides[rxn.products[i]].Z;
    }
    assert(A_in == A_out);
    assert(Z_in == Z_out);
    /* Every reaction in this network is exothermic in the written direction. */
    assert(ReactionQValueMeV(rxn) > 0.0);
  }
  /* Spot-check two Q-values against their measured values, so a corrupted mass
     table cannot pass merely by being self-consistent. */
  assert(Close(ReactionQValueMeV(kReactions[kRxnNPG]), 2.2246, 1.0e-4));
  assert(Close(ReactionQValueMeV(kReactions[kRxnTDN]), 17.589, 1.0e-4));
}

void test_identical_particle_prefactors() {
  /* d + d must carry 1/2; li7 + p -> he4 + he4 must carry it on the product side. */
  assert(ReactionPrefactor(kReactions[kRxnDDN].reactants, 2) == 0.5);
  assert(ReactionPrefactor(kReactions[kRxnDDN].products, 2) == 1.0);
  assert(ReactionPrefactor(kReactions[kRxnLi7PA].reactants, 2) == 1.0);
  assert(ReactionPrefactor(kReactions[kRxnLi7PA].products, 2) == 0.5);
}

void test_detailed_balance() {
  const BbnRates rates(RatesFile());
  for (int r = 0; r < kBbnNumReactions; ++r) {
    for (double T9 : {0.05, 0.2, 1.0, 3.0, 10.0}) {
      const double got  = rates.ReverseOverForward(r, T9);
      const double want = ExpectedReverseOverForward(kReactions[r], T9);
      assert(Close(got, want, 1.0e-10));
      /* And the composed reverse rate is exactly forward times that ratio. */
      assert(Close(rates.Reverse(r, T9), rates.Forward(r, T9) * got, 1.0e-12));
    }
  }
}

void test_rates_are_positive_and_finite() {
  const BbnRates rates(RatesFile());
  for (int r = 0; r < kBbnNumReactions; ++r) {
    for (double T9 = 0.01; T9 <= 100.0; T9 *= 1.4) {
      const double f = rates.Forward(r, T9);
      assert(std::isfinite(f) && f >= 0.0);
      assert(std::isfinite(rates.Reverse(r, T9)));
    }
  }
}

/* The fit ceiling: above kReaclibT9Max the fit argument is frozen, so the forward
   rate stops changing -- while the reverse rate keeps responding to temperature
   through the Q-value, because that is thermodynamics rather than a fit. Without
   the ceiling the a5*T9^{5/3} term reaches exp(207) at T9 = 100. */
void test_fit_ceiling() {
  const BbnRates rates(RatesFile());
  for (int r = 0; r < kBbnNumReactions; ++r) {
    const double at_max = rates.Forward(r, kReaclibT9Max);
    assert(Close(rates.Forward(r, 50.0), at_max, 1.0e-14));
    assert(Close(rates.Forward(r, 100.0), at_max, 1.0e-14));
    assert(rates.Forward(r, 1.0) != at_max);
  }
  /* The reverse rate must still vary above the ceiling, or equilibrium would be
     frozen at the wrong composition. */
  assert(!Close(rates.ReverseOverForward(kRxnNPG, 100.0),
                rates.ReverseOverForward(kRxnNPG, 50.0),
                1.0e-6));
}

/* Absence must be fatal, never silently zero: a reaction with no coefficients
   would evaluate to a perfectly valid-looking rate of exactly zero. */
void test_guards_fail_loudly() {
  bool threw = false;
  try {
    BbnRates missing("/nonexistent/rates.dat");
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);

  /* A file naming a reaction the code does not know. */
  const std::string unknown = TempFile("unknown.dat");
  {
    std::ofstream out(unknown);
    out << "not_a_reaction  xx  1 0 0 0 0 0 0\n";
  }
  threw = false;
  try {
    BbnRates bad(unknown);
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);

  /* A file that is well-formed but incomplete. */
  const std::string incomplete = TempFile("incomplete.dat");
  {
    std::ofstream out(incomplete);
    out << kReactions[kRxnNPG].name << "  xx  1 0 0 0 0 0 0\n";
  }
  threw = false;
  try {
    BbnRates partial(incomplete);
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);

  /* A malformed coefficient row. */
  const std::string malformed = TempFile("malformed.dat");
  {
    std::ofstream out(malformed);
    out << kReactions[kRxnNPG].name << "  xx  1 2 3\n";
  }
  threw = false;
  try {
    BbnRates bad_row(malformed);
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);
}

}  // namespace

int main() {
  test_network_topology_is_conserved();
  test_identical_particle_prefactors();
  test_detailed_balance();
  test_rates_are_positive_and_finite();
  test_fit_ceiling();
  test_guards_fail_loudly();
  std::printf("bbn rates tests passed\n");
  return 0;
}
