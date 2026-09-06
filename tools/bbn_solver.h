#pragma once

/** @file bbn_solver.h Coupled stiff solve for the primordial light-element abundances.
 *
 *  Design: docs/superpowers/specs/2026-09-06-bbn-nuclear-network-design.md section 3
 *
 *  One coupled system, as Kawano/PArthENoPE/AlterBBN do, rather than a plasma
 *  stage splined into a network stage. Twelve equations: the scale factor, the
 *  baryon-to-photon ratio, time, and nine abundances, with the photon temperature
 *  as the INDEPENDENT variable (it is monotonic, and it makes the integration
 *  range exactly known instead of estimated).
 *
 *  This is a cosmology-free library: it takes numbers, not a CLASS module, so
 *  bbn_solver_test can drive it without constructing a Cosmology.
 */

#include <functional>
#include <string>
#include <vector>

#include "bbn_nuclides.h"

struct BbnInput {
  double omega_b    = 0.02242; /**< Omega_b h^2 */
  double delta_neff = 0.0;     /**< N_eff - 3.046, the sBBN table's axis convention */
  double tau_n      = 878.4;   /**< free neutron lifetime, s */
  double T_cmb      = 2.7255;  /**< photon temperature today, K */

  double T9_initial = 100.0;  /**< start of the integration, 10^9 K */
  double T9_final   = 0.01;   /**< end of the integration, 10^9 K */
  double tolerance  = 1.0e-8; /**< relative tolerance handed to the evolver */

  std::string rates_file; /**< path to the REACLIB data file; required */

  /** Nodes for the splined weak-rate table, or 0 to evaluate the phase-space
   *  integrals directly at every step.
   *
   *  The integrals are ~87% of a right-hand-side evaluation, so tabulating them
   *  is where nearly all of the runtime goes. 0 is the reference path, kept so
   *  that bbn_weak_test and bbn_solver_test can measure what the table costs in
   *  accuracy rather than assume it is free. */
  int weak_table_points = 512;

  /** Optional extra energy density in MeV^4 as a function of the scale factor.
   *
   *  The seam of spec section 8. Left null, the expansion is Standard-Model plus
   *  whatever delta_neff asks for. Filled, it lets a decaying species or any
   *  other exotic component drive BBN directly, which a single delta_neff cannot
   *  express. Nothing in the prototype fills it; the interface exists so that it
   *  does not have to be reopened later. */
  std::function<double(double a)> rho_extra;

  /** If non-empty, write the abundance history here (one row per sampled T9). */
  std::string history_file;
  int history_rows = 400;
};

struct BbnResult {
  /** The two candidate helium conventions, both reported.
   *
   *  They differ through the He4 binding energy by roughly half a percent, and
   *  spec section 6.1 deliberately declines to assert which one CLASS's YHe
   *  means -- the comparison against sBBN_2017.dat decides it. Reporting both
   *  keeps that an empirical question rather than a silent assumption. */
  double Yp_nucleon = 0.0; /**< 4 * Y_He4, the nucleon-number fraction */
  double Yp_mass    = 0.0; /**< rho_He4 / rho_b, the true mass fraction */

  double DoverH   = 0.0;
  double He3overH = 0.0; /**< includes tritium, which beta-decays to He3 */
  double Li7overH = 0.0; /**< includes Be7, which electron-captures to Li7 */
  double Li6overH = 0.0;

  double eta_final     = 0.0; /**< baryon-to-photon ratio today */
  double eta10         = 0.0; /**< 10^10 * eta_final */
  double aT_ratio      = 0.0; /**< (aT)_final / (aT)_initial; (11/4)^(1/3) if instantaneous */
  double neutrino_neff = 0.0; /**< total effective species used in the expansion */
  double time_final    = 0.0; /**< s */

  std::vector<double> abundances; /**< final Y_i = n_i/n_b, indexed by BbnNuclideIndex */
};

/** Solve. Throws on a malformed rate file or a failed integration. */
BbnResult SolveBbn(const BbnInput& input);
