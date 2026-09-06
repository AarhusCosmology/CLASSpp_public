#pragma once

/** @file bbn_nuclides.h Nuclide data for the primordial nucleosynthesis network.
 *
 *  Design: docs/superpowers/specs/2026-09-06-bbn-nuclear-network-design.md
 *
 *  This table is the code-side half of the "topology in code, rates in a data
 *  file" split (spec section 4). It exists so that Q-values and detailed-balance
 *  reverse rates are computed from ONE set of masses and degeneracies, and can
 *  therefore never disagree with each other.
 *
 *  Masses are ATOMIC masses (nucleus + Z electrons), the convention nuclear mass
 *  excesses are tabulated in. That is the right choice for Q-values, because Z is
 *  conserved by every reaction here and the electron masses cancel exactly.
 *  Detailed balance needs NUCLEAR masses instead, so NuclearMassMeV() subtracts
 *  them back off; see the comment there.
 */

#include <cstddef>

/** Electron mass, MeV. Also the scale at which e+- annihilation reheats the
 *  photons, so it appears in the plasma sector too. */
inline constexpr double kElectronMassMeV = 0.51099895000;

struct BbnNuclide {
  const char* name;
  int A;                  /**< mass number */
  int Z;                  /**< charge number */
  double g;               /**< spin degeneracy 2J+1 */
  double atomic_mass_mev; /**< atomic mass (nucleus + Z electrons), MeV */
};

/** Index into kNuclides. The ORDER IS THE STATE-VECTOR ORDER of the solver, so
 *  it is part of the interface: bbn_solver.cpp indexes abundances by these. */
enum BbnNuclideIndex {
  kNucN = 0,
  kNucP,
  kNucD,
  kNucT,
  kNucHe3,
  kNucHe4,
  kNucLi6,
  kNucLi7,
  kNucBe7,
  kBbnNumNuclides
};

/** The nine-nuclide Kawano key network. Degeneracies are 2J+1 for the ground
 *  state; every nuclide here has a ground state far below the ~MeV temperatures
 *  of interest, so partition functions are 1 and detailed balance needs no
 *  partition-function ratio.
 *
 *  Values from AME2020 via pynucastro 3.0.0 (see scripts/gen_bbn_rates.py). */
inline constexpr BbnNuclide kNuclides[kBbnNumNuclides] = {
    /* name    A  Z  g    atomic mass [MeV] */
    {"n", 1, 0, 2.0, 939.56542052},
    {"p", 1, 1, 2.0, 938.78307348},
    {"d", 2, 1, 3.0, 1876.12392774},
    {"t", 3, 1, 2.0, 2809.43211816},
    {"he3", 3, 2, 2.0, 2809.41352614},
    {"he4", 4, 2, 1.0, 3728.40132555},
    {"li6", 6, 3, 3.0, 5603.05149492},
    {"li7", 7, 3, 4.0, 6535.36582194},
    {"be7", 7, 4, 4.0, 6536.22771694},
};

/** Nuclear (not atomic) mass, MeV.
 *
 *  Detailed balance weights each species by its own m^{3/2} through the thermal
 *  factor (m T / 2 pi hbar^2)^{3/2}, and that factor belongs to the NUCLEUS: the
 *  atomic electrons are not participating in the reaction and are not present at
 *  BBN temperatures in bound form anyway. Unlike in the Q-value, these do not
 *  cancel between the two sides -- a radiative capture has two nuclei on one side
 *  and one on the other -- so using atomic masses here would be a ~0.05% error per
 *  nuclide in the reverse rates. Electron binding energy (eV scale) is ignored. */
inline constexpr double NuclearMassMeV(const BbnNuclide& nuc) {
  return nuc.atomic_mass_mev - nuc.Z * kElectronMassMeV;
}

/** Neutron-proton mass difference, MeV. Sets the n/p equilibrium ratio and hence,
 *  through freeze-out, essentially all of Y_p. Atomic masses would give the
 *  neutron-to-hydrogen difference, which is wrong by an electron mass here, so
 *  this is deliberately built from the nuclear masses. */
inline constexpr double kDeltaMnpMeV = NuclearMassMeV(kNuclides[kNucN]) -
                                       NuclearMassMeV(kNuclides[kNucP]);
