#include "bbn_rates.h"

#include <cmath>
#include <fstream>
#include <set>
#include <sstream>

#include "errors.h"

namespace {

/** hbar*c in MeV cm. Converts natural-unit number densities (MeV^3) to cm^-3,
 *  which is what the REACLIB rate units expect. */
constexpr double kHbarCMeVcm = 1.973269804e-11;

/** (1 MeV)^3 expressed in cm^-3. */
const double kMeV3ToInvCm3 = 1.0 / (kHbarCMeVcm * kHbarCMeVcm * kHbarCMeVcm);

constexpr double kAvogadro = 6.02214076e23;

/** Temperature-independent part of the thermal number-density scale: the full
 *  scale is this times T^{3/2}.
 *
 *  It is the factor in the Maxwell-Boltzmann equilibrium density
 *  n^eq = g (m T / 2 pi hbar^2)^{3/2} exp((mu - m)/T), and it is where the
 *  NUCLEAR rather than atomic mass is required: the atomic electrons take no part
 *  in the reaction, and unlike in the Q-value these factors do not cancel between
 *  the two sides when the two sides have different particle counts.
 *
 *  Split off the temperature so that the whole detailed-balance constant can be
 *  formed once in the constructor. Re-forming it per nuclide per reaction on
 *  every right-hand side was a pow() call per nuclide, for a number that never
 *  changes. */
double ThermalScaleConstant(const BbnNuclide& nuc) {
  const double theta = NuclearMassMeV(nuc) / (2.0 * M_PI);
  return std::pow(theta, 1.5) * kMeV3ToInvCm3;
}

}  // namespace

BbnRates::BbnRates(const std::string& rates_file) : source_file_(rates_file) {
  sets_.resize(kBbnNumReactions);

  std::ifstream file(rates_file);
  class_test_severe(!file.is_open(),
                    "could not open BBN reaction-rate file '%s'",
                    rates_file.c_str());

  std::set<std::string> labels;
  std::string line;
  int line_number = 0;

  while (std::getline(file, line)) {
    ++line_number;
    const std::size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] == '#') {
      continue;
    }

    std::istringstream in(line);
    std::string name, label;
    FitSet set;
    in >> name >> label;
    for (double& coefficient : set.a) {
      in >> coefficient;
    }
    class_test_severe(in.fail(),
                      "malformed line %d of BBN rate file '%s': expected "
                      "<reaction> <label> a0..a6, got '%s'",
                      line_number,
                      rates_file.c_str(),
                      line.c_str());

    /* Name matching is the whole safety story of the data file: it may supply
       numbers for reactions the code knows about, and nothing else. An unknown
       name is far more likely to be a typo silently disabling a real rate than a
       deliberate extension, so it is fatal rather than ignored. */
    int reaction = -1;
    for (int r = 0; r < kBbnNumReactions; ++r) {
      if (name == kReactions[r].name) {
        reaction = r;
        break;
      }
    }
    class_test_severe(reaction < 0,
                      "BBN rate file '%s' line %d names reaction '%s', which is not in "
                      "the network defined by tools/bbn_reactions.h",
                      rates_file.c_str(),
                      line_number,
                      name.c_str());

    sets_[reaction].push_back(set);
    labels.insert(label);
  }

  /* Absence must fail, not skip: a reaction with no coefficients would evaluate
     to exactly zero, which is a physically meaningful (and wrong) rate rather
     than an obvious error. */
  for (int r = 0; r < kBbnNumReactions; ++r) {
    class_test_severe(sets_[r].empty(),
                      "BBN rate file '%s' supplies no rate for reaction '%s'",
                      rates_file.c_str(),
                      kReactions[r].name);
  }

  for (const std::string& label : labels) {
    provenance_ += (provenance_.empty() ? "" : " ") + label;
  }

  /* Detailed balance, precomputed. At equilibrium the two directions have equal
     rates per volume, so
         pref_f lambda_f prod_R n_i^eq  =  pref_r lambda_r prod_P n_j^eq
     and every chemical potential cancels because the reaction conserves them.
     What survives is degeneracies, thermal scales and the Q-value -- and only the
     last two carry any temperature dependence. */
  balance_.resize(kBbnNumReactions);
  for (int r = 0; r < kBbnNumReactions; ++r) {
    const BbnReaction& rxn = kReactions[r];
    DetailedBalance& db    = balance_[r];

    db.constant = ReactionPrefactor(rxn.reactants, rxn.n_reactants) /
                  ReactionPrefactor(rxn.products, rxn.n_products);
    for (int i = 0; i < rxn.n_reactants; ++i) {
      const BbnNuclide& nuc  = kNuclides[rxn.reactants[i]];
      db.constant           *= nuc.g * ThermalScaleConstant(nuc);
    }
    for (int i = 0; i < rxn.n_products; ++i) {
      const BbnNuclide& nuc  = kNuclides[rxn.products[i]];
      db.constant           /= nuc.g * ThermalScaleConstant(nuc);
    }
    db.q_mev   = ReactionQValueMeV(rxn);
    db.surplus = rxn.n_reactants - rxn.n_products;
  }
}

double BbnRates::Forward(int reaction, double T9) const {
  /* REACLIB fits are valid over roughly T9 in [0.01, 10] and MUST NOT be
     extrapolated above that: the a4*T9 and a5*T9^{5/3} terms are fitting
     conveniences, not physics, and at T9 = 100 they produce exp(207) -- a rate
     wrong by ninety orders of magnitude, which presents as the integrator
     failing rather than as a bad rate.
 
     Holding the fit argument at its ceiling is safe because of what the rate
     does up there. Above T9 ~ 10 every nuclide heavier than a nucleon sits in
     nuclear statistical equilibrium, where the ABUNDANCES are fixed by detailed
     balance and the rate magnitude only decides how quickly equilibrium is
     restored. So the ceiling is applied to the fit argument alone: the
     Q-value and thermal factors in ReverseOverForward always see the true
     temperature, and the equilibrium the network relaxes to is therefore still
     exactly right. bbn_solver_test checks this by varying T9_initial.
 
     The integration cannot simply start below the ceiling instead: at T9 = 10
     the weak rates are only ~1.6 times H, so the n/p ratio is already leaving
     equilibrium and an equilibrium initial condition would be wrong. */
  const double T9_fit = T9 < kReaclibT9Max ? T9 : kReaclibT9Max;

  double rate            = 0.0;
  const double cube_root = std::cbrt(T9_fit);
  const double log_t9    = std::log(T9_fit);

  for (const FitSet& set : sets_[reaction]) {
    /* The JINA REACLIB seven-parameter form. Verified against pynucastro's own
       evaluator over T9 in [0.01, 10] by scripts/gen_bbn_rates.py at generation
       time, so a transposed coefficient cannot reach this file unnoticed. */
    const double exponent  = set.a[0] + set.a[1] / T9_fit + set.a[2] / cube_root +
                             set.a[3] * cube_root + set.a[4] * T9_fit +
                             set.a[5] * T9_fit * cube_root * cube_root + set.a[6] * log_t9;
    rate                  += std::exp(exponent);
  }

  /* The file carries N_A<sigma v>; the network works in number densities. */
  return rate / kAvogadro;
}

double BbnRates::ReverseOverForward(int reaction, double T9) const {
  const DetailedBalance& db = balance_[reaction];
  const double T_mev        = MeVFromT9(T9);

  /* Every reaction in this network has surplus 0 (two bodies on each side) or 1
     (radiative capture), so the temperature power is either absent or a single
     T^{3/2}; the general pow() remains only as a fallback. */
  double thermal = 1.0;
  if (db.surplus == 1) {
    thermal = T_mev * std::sqrt(T_mev);
  }
  else if (db.surplus != 0) {
    thermal = std::pow(T_mev, 1.5 * db.surplus);
  }

  /* Exothermic Q > 0 suppresses the reverse direction. This underflows to zero
     at low temperature for the large-Q reactions, which is correct and is why it
     is written as a bare exp rather than guarded: zero is a finite number, and
     the -ffast-math builds only forbid materialising infinities. */
  return db.constant * thermal * std::exp(-db.q_mev / T_mev);
}

double BbnRates::Reverse(int reaction, double T9) const {
  return Forward(reaction, T9) * ReverseOverForward(reaction, T9);
}
