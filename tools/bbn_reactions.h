#pragma once

/** @file bbn_reactions.h Reaction topology for the primordial nucleosynthesis network.
 *
 *  Design: docs/superpowers/specs/2026-09-06-bbn-nuclear-network-design.md section 4
 *
 *  Topology lives here, in code; forward rate COEFFICIENTS live in the data file
 *  bbn/rates_reaclib.dat and are matched to these entries by name. The split is
 *  deliberate: a rate file can supply wrong numbers, but it cannot invent a
 *  reaction, change what a reaction does, or contradict a Q-value, because none
 *  of those are expressible in it.
 *
 *  Multiplicity is expressed by REPEATING an index: d + d appears as {kNucD, kNucD}.
 *  That keeps the identical-particle factors derivable (see ReactionPrefactor)
 *  rather than tabulated, so they cannot be forgotten for a new reaction.
 *
 *  Photons are ABSENT from both sides. They carry no chemical potential and no
 *  number density that enters the rate equations, so a radiative capture a+b -> c
 *  is stored with two reactants and one product; the resulting asymmetry in the
 *  massive-particle count is exactly what turns its reverse rate into a
 *  photodisintegration rate with units of s^-1. See bbn_rates.cpp.
 */

#include "bbn_nuclides.h"

inline constexpr int kBbnMaxReactionSide = 2;

struct BbnReaction {
  const char* name; /**< must match a name in the rate file */
  int n_reactants;
  int n_products;
  int reactants[kBbnMaxReactionSide]; /**< nuclide indices, repeated for multiplicity */
  int products[kBbnMaxReactionSide];
};

enum BbnReactionIndex {
  kRxnNPG = 0, /* n(p,g)d        -- the deuterium bottleneck; sets when BBN starts */
  kRxnDPG,     /* d(p,g)he3      -- the dominant D destruction channel             */
  kRxnDDN,     /* d(d,n)he3                                                        */
  kRxnDDP,     /* d(d,p)t                                                          */
  kRxnHe3NP,   /* he3(n,p)t                                                        */
  kRxnTDN,     /* t(d,n)he4      -- the main route into He4                        */
  kRxnHe3DP,   /* he3(d,p)he4                                                      */
  kRxnHe3AG,   /* he3(he4,g)be7  -- the dominant mass-7 production channel         */
  kRxnTAG,     /* t(he4,g)li7                                                      */
  kRxnBe7NP,   /* be7(n,p)li7                                                      */
  kRxnLi7PA,   /* li7(p,he4)he4  -- the main mass-7 destruction channel            */
  kRxnADG,     /* he4(d,g)li6                                                      */
  kRxnLi6PA,   /* li6(p,he4)he3                                                    */
  kRxnTPG,     /* t(p,g)he4                                                        */
  kBbnNumReactions
};

/** The fourteen strong reactions of the Kawano key network.
 *
 *  Free neutron decay and the two-directional weak n <-> p rates are NOT here.
 *  They are not rate-file data: they come from phase-space integrals normalized to
 *  the tau_n input (bbn_weak.h), and they are the only reactions in the problem
 *  that do not conserve nucleon species counts through a nuclear matrix element. */
inline constexpr BbnReaction kReactions[kBbnNumReactions] = {
    /*  name              nR nP  reactants              products              */
    {"n(p,g)d", 2, 1, {kNucN, kNucP}, {kNucD, -1}},
    {"d(p,g)he3", 2, 1, {kNucD, kNucP}, {kNucHe3, -1}},
    {"d(d,n)he3", 2, 2, {kNucD, kNucD}, {kNucN, kNucHe3}},
    {"d(d,p)t", 2, 2, {kNucD, kNucD}, {kNucP, kNucT}},
    {"he3(n,p)t", 2, 2, {kNucHe3, kNucN}, {kNucP, kNucT}},
    {"t(d,n)he4", 2, 2, {kNucT, kNucD}, {kNucN, kNucHe4}},
    {"he3(d,p)he4", 2, 2, {kNucHe3, kNucD}, {kNucP, kNucHe4}},
    {"he3(he4,g)be7", 2, 1, {kNucHe3, kNucHe4}, {kNucBe7, -1}},
    {"t(he4,g)li7", 2, 1, {kNucT, kNucHe4}, {kNucLi7, -1}},
    {"be7(n,p)li7", 2, 2, {kNucBe7, kNucN}, {kNucP, kNucLi7}},
    {"li7(p,he4)he4", 2, 2, {kNucLi7, kNucP}, {kNucHe4, kNucHe4}},
    {"he4(d,g)li6", 2, 1, {kNucHe4, kNucD}, {kNucLi6, -1}},
    {"li6(p,he4)he3", 2, 2, {kNucLi6, kNucP}, {kNucHe4, kNucHe3}},
    {"t(p,g)he4", 2, 1, {kNucT, kNucP}, {kNucHe4, -1}},
};

/** Identical-particle factor 1/prod(nu_i!) for one side of a reaction.
 *
 *  With at most two particles a side, the only case that bites is a pair of
 *  identical ones, giving 1/2. Derived from the index list rather than tabulated
 *  so that adding d+d -> ... or ... -> he4+he4 cannot silently omit it. */
inline constexpr double ReactionPrefactor(const int* side, int n) {
  return (n == 2 && side[0] == side[1]) ? 0.5 : 1.0;
}

/** Q-value in MeV, computed from the atomic masses in kNuclides.
 *
 *  Positive means exothermic. Never stored: a tabulated Q could drift out of
 *  agreement with the masses that detailed balance uses, and that disagreement
 *  would show up only as a slow unphysical drift away from equilibrium. Atomic
 *  masses are correct here because every reaction conserves Z, so the electron
 *  masses cancel between the two sides. */
inline constexpr double ReactionQValueMeV(const BbnReaction& r) {
  double q = 0.0;
  for (int i = 0; i < r.n_reactants; ++i) {
    q += kNuclides[r.reactants[i]].atomic_mass_mev;
  }
  for (int i = 0; i < r.n_products; ++i) {
    q -= kNuclides[r.products[i]].atomic_mass_mev;
  }
  return q;
}
