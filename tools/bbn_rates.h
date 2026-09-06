#pragma once

/** @file bbn_rates.h Forward rates from a REACLIB data file; reverse rates by detailed balance.
 *
 *  Design: docs/superpowers/specs/2026-09-06-bbn-nuclear-network-design.md section 4
 *
 *  Forward rates are DATA, read from bbn/rates_reaclib.dat. Reverse rates are not:
 *  they are computed here from the nuclide degeneracies and masses in
 *  bbn_nuclides.h. There is deliberately no way to supply a reverse rate, so a
 *  wrong rate file produces a wrong abundance -- which a test can catch -- rather
 *  than a violation of detailed balance, which shows up only as a slow drift away
 *  from nuclear statistical equilibrium and looks like a solver problem.
 */

#include <string>
#include <vector>

#include "bbn_reactions.h"

/** Rate units, stated once because getting them wrong is the classic BBN bug:
 *
 *    two-body forward       <sigma v>            [cm^3 s^-1]   (already /N_A)
 *    two-body reverse       <sigma v>            [cm^3 s^-1]
 *    radiative-capture reverse (photodisintegration)           [s^-1]
 *
 *  The asymmetry is real, not an inconsistency: a radiative capture a+b -> c has
 *  two massive particles on one side and one on the other, so detailed balance
 *  carries a leftover thermal factor with dimensions of number density. */
class BbnRates {
 public:
  /** Loads and validates the whole network. Throws (severe) if the file is
   *  missing, malformed, names a reaction the code does not know, or omits one it
   *  does. There are no defaults: a missing rate must fail loudly rather than
   *  quietly produce a plausible-looking abundance. */
  explicit BbnRates(const std::string& rates_file);

  /** Forward rate at T9 = T / 10^9 K. */
  double Forward(int reaction, double T9) const;

  /** Reverse rate at T9, by detailed balance from the forward rate. */
  double Reverse(int reaction, double T9) const;

  /** Reverse/forward ratio -- the pure detailed-balance factor, independent of
   *  any rate data. Exposed so bbn_rates_test can check it against an
   *  independently computed Saha equilibrium without going through a rate file. */
  double ReverseOverForward(int reaction, double T9) const;

  const std::string& source_file() const {
    return source_file_;
  }

  /** Provenance label(s) carried from the data file, for the verbose printout. */
  const std::string& provenance() const {
    return provenance_;
  }

 private:
  /** One REACLIB seven-parameter set. Several may apply to one reaction; they sum. */
  struct FitSet {
    double a[7];
  };

  /** Detailed balance factorizes as C * T^(3/2 * surplus) * exp(-Q/T), where the
   *  constant collects the degeneracies, the masses and the unit conversion. None
   *  of that depends on temperature, so it is formed once in the constructor
   *  rather than rebuilt per nuclide per reaction on every right-hand side --
   *  which is where a pow() call per nuclide used to come from. */
  struct DetailedBalance {
    double constant; /**< C, including the identical-particle prefactors */
    double q_mev;    /**< Q-value, positive for exothermic */
    int surplus;     /**< massive reactants minus massive products */
  };

  std::vector<std::vector<FitSet>> sets_; /**< indexed by BbnReactionIndex */
  std::vector<DetailedBalance> balance_;  /**< indexed by BbnReactionIndex */
  std::string source_file_;
  std::string provenance_;
};

/** Temperature conversions. T9 is the natural variable for reaction rates
 *  (REACLIB is written in it); MeV is the natural variable for everything else. */
/** Upper end of the REACLIB fits' validity. Above this the fit argument is held
 *  fixed; see the comment in BbnRates::Forward for why that is safe and why
 *  simply starting the integration lower is not. */
inline constexpr double kReaclibT9Max = 10.0;

inline constexpr double kMeVPerT9 = 8.617333262e-2; /* k_B * 10^9 K, in MeV */
inline constexpr double kT9PerMeV = 1.0 / kMeVPerT9;

inline constexpr double MeVFromT9(double T9) {
  return T9 * kMeVPerT9;
}
inline constexpr double T9FromMeV(double T) {
  return T * kT9PerMeV;
}
