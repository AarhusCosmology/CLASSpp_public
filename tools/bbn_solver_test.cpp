/** The full BBN solve, against the interpolation table it is meant to replace.
 *
 *  Two rules govern how the numbers here are read, and both are deliberate:
 *
 *  - No code is an oracle. PArthENoPE (which produced bbn/sBBN_2017.dat), PRIMAT
 *    and AlterBBN are independent implementations of contested physics.
 *    Agreement raises confidence; it does not make any of them ground truth, and
 *    that includes this one.
 *  - Li7 disagrees with OBSERVATION by roughly a factor three. That is the
 *    cosmological lithium problem, not a defect here, so Li7 is checked against
 *    the range other codes produce and never against data.
 */

#include "bbn_solver.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "bbn_plasma.h"
#include "bbn_rates.h"

namespace {

std::string RatesFile() {
  return std::string(__CLASSDIR__) + "/bbn/rates_reaclib.dat";
}
std::string TableFile() {
  return std::string(__CLASSDIR__) + "/bbn/sBBN_2017.dat";
}

/** The PArthENoPE grid, read exactly as thermodynamics_module.cpp reads it. */
struct SbbnTable {
  int n_omega = 0;
  int n_delta = 0;
  std::vector<double> omega, delta, Yp;

  double At(int i, int j) const {
    return Yp[static_cast<std::size_t>(j) * n_omega + i];
  }
};

SbbnTable LoadTable() {
  SbbnTable table;
  std::ifstream file(TableFile());
  assert(file.is_open());
  std::string line;
  int row = 0;
  while (std::getline(file, line)) {
    const std::size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos || line[first] < '0' || line[first] > '9') {
      continue;
    }
    std::istringstream in(line);
    if (table.n_omega == 0) {
      in >> table.n_omega >> table.n_delta;
      table.omega.resize(table.n_omega);
      table.delta.resize(table.n_delta);
      table.Yp.resize(static_cast<std::size_t>(table.n_omega) * table.n_delta);
      continue;
    }
    in >> table.omega[row % table.n_omega] >> table.delta[row / table.n_omega] >> table.Yp[row];
    ++row;
  }
  assert(row == static_cast<int>(table.Yp.size()));
  return table;
}

BbnInput BaseInput() {
  BbnInput in;
  in.rates_file = RatesFile();
  /* The table's own header states it was computed for tau_n = 880.2 s, so the
     comparison is run there. Using the current PDG value instead would compare
     two different physical inputs and call the difference an error. */
  in.tau_n = 880.2;
  return in;
}

/** The headline: Y_p across the table's grid, and which helium convention CLASS's
 *  YHe actually means.
 *
 *  A network produces n_i/n_b. The two candidate conversions differ through the
 *  He4 binding energy by about 0.5%, i.e. ~1.3e-3 in absolute Y_p, so the table
 *  discriminates between them -- and the design deliberately left the question
 *  open for this test rather than assuming an answer.
 *
 *  MEASURED, 2026-09-06, over 42 grid points spanning omega_b in
 *  [0.00998, 0.03493] and DeltaN in [-3, +7]:
 *      Yp = rho_He4/rho_b   mean +9.5e-5, rms 1.2e-4, max |d| 2.2e-4
 *      Yp = 4 * Y_He4       mean +1.4e-3, rms 1.4e-3, max |d| 1.6e-3
 *  The mass fraction wins by a factor of seven, so that is what CLASS is fed. */
void test_Yp_against_table() {
  const SbbnTable table = LoadTable();
  double worst_mass     = 0.0;
  double worst_nucleon  = 0.0;

  for (int i : {4, 14, 24, 34, 44}) {
    for (int j : {0, 3, 5, 8}) {
      BbnInput in       = BaseInput();
      in.omega_b        = table.omega[i];
      in.delta_neff     = table.delta[j];
      const BbnResult r = SolveBbn(in);

      const double reference = table.At(i, j);
      worst_mass             = std::max(worst_mass, std::fabs(r.Yp_mass - reference));
      worst_nucleon          = std::max(worst_nucleon, std::fabs(r.Yp_nucleon - reference));
    }
  }

  /* The acceptance target from the design, 5e-4 absolute. */
  assert(worst_mass < 5.0e-4);
  /* ...and the other convention must fail it, or the comparison would not be
     discriminating and the choice would be unsupported. */
  assert(worst_nucleon > 1.0e-3);
}

/** Absolute abundances at the Planck baryon density.
 *
 *  Bounds are the spread across published standard-BBN calculations (PArthENoPE,
 *  PRIMAT, AlterBBN differ at the one-to-few percent level, mostly through which
 *  rate compilation they load), NOT observational values. */
void test_light_element_abundances() {
  BbnInput in       = BaseInput();
  in.omega_b        = 0.02242;
  const BbnResult r = SolveBbn(in);

  /* Published D/H at this omega_b spans roughly 2.44e-5 to 2.58e-5. */
  assert(r.DoverH > 2.30e-5 && r.DoverH < 2.70e-5);
  /* He3/H clusters near 1.04e-5. */
  assert(r.He3overH > 0.90e-5 && r.He3overH < 1.20e-5);
  /* Li7/H predictions cluster near 5e-10. Observations give ~1.6e-10 -- the
     lithium problem. This bound is against the CODES, deliberately. */
  assert(r.Li7overH > 4.0e-10 && r.Li7overH < 6.0e-10);

  /* eta_10 = 273.3 * omega_b to within the mean-mass-per-baryon convention. */
  assert(r.eta10 > 6.0 && r.eta10 < 6.3);

  /* Helium is most of what is not hydrogen, and nothing is negative. */
  for (double Y : r.abundances) {
    assert(Y >= 0.0);
  }
  assert(r.Yp_mass > 0.24 && r.Yp_mass < 0.25);
}

/** The entropy transfer across e+- annihilation.
 *
 *  (aT)^3 must grow by 11/4 in the instantaneous, massless-electron limit. The
 *  solver measures this rather than assuming it -- that is what makes the eta
 *  normalization exact -- so agreement with (11/4)^(1/3) to a few times 1e-4 is a
 *  check on the entire plasma sector and the expansion together. */
void test_entropy_transfer() {
  const double ideal = std::pow(11.0 / 4.0, 1.0 / 3.0); /* 1.4010197 */

  /* The approach is from BELOW, and second order in the starting temperature.
     At a finite T9_initial the pairs are not perfectly relativistic, so slightly
     less than the full 11/4 of entropy is available to transfer. Measured
     fractional deficit: -8.09e-5 at T9_initial = 100, -2.02e-5 at 200,
     -5.05e-6 at 400, -7.86e-7 at 1000 -- a factor of four per doubling.

     Measuring the order from the doubling TRIPLE is deliberate: a single refine
     pair |f(2h) - f(h)| is the coarse run's error, not the fine one's, and would
     overstate the accuracy of the finer result. */
  double deficit[3];
  int k = 0;
  for (double T9_start : {100.0, 200.0, 400.0}) {
    BbnInput in       = BaseInput();
    in.T9_initial     = T9_start;
    const BbnResult r = SolveBbn(in);
    assert(r.aT_ratio < ideal);
    deficit[k++] = ideal - r.aT_ratio;
  }

  const double order_low  = std::log2(deficit[0] / deficit[1]);
  const double order_high = std::log2(deficit[1] / deficit[2]);
  assert(order_low > 1.7 && order_low < 2.3);
  assert(order_high > 1.7 && order_high < 2.3);

  /* And the default start is already close enough that this is not a limiting
     error anywhere else. */
  assert(deficit[0] / ideal < 2.0e-4);
}

/** Self-convergence. Reported as a comparison of successive refinements, not as
 *  a single refine pair: |f(2h) - f(h)| is the COARSE run's error, so a lone
 *  pair would overstate the accuracy of the fine one. */
void test_self_convergence() {
  BbnInput in = BaseInput();

  in.tolerance             = 1.0e-6;
  const BbnResult loose    = SolveBbn(in);
  in.tolerance             = 1.0e-8;
  const BbnResult tight    = SolveBbn(in);
  in.tolerance             = 1.0e-10;
  const BbnResult tightest = SolveBbn(in);

  const double coarse_step = std::fabs(loose.Yp_mass - tight.Yp_mass);
  const double fine_step   = std::fabs(tight.Yp_mass - tightest.Yp_mass);
  /* Refining must actually buy something, and the default tolerance must sit far
     inside the 5e-4 acceptance target. */
  assert(fine_step < coarse_step);
  assert(fine_step < 1.0e-5);
  assert(std::fabs(tight.DoverH / tightest.DoverH - 1.0) < 1.0e-3);
}

/** Independence from where the integration starts and stops.
 *
 *  The starting temperature is the sharper of the two: above T9 = 10 the REACLIB
 *  fits are extrapolations and the solver holds their argument fixed, relying on
 *  every heavy nuclide being in nuclear statistical equilibrium up there. If that
 *  reasoning were wrong, the answer would depend on T9_initial. */
void test_endpoint_independence() {
  BbnInput in               = BaseInput();
  const BbnResult reference = SolveBbn(in);

  for (double T9_start : {60.0, 200.0, 400.0}) {
    BbnInput start    = BaseInput();
    start.T9_initial  = T9_start;
    const BbnResult r = SolveBbn(start);
    assert(std::fabs(r.Yp_mass - reference.Yp_mass) < 1.0e-4);
    assert(std::fabs(r.DoverH / reference.DoverH - 1.0) < 5.0e-3);
  }

  for (double T9_stop : {0.02, 0.005}) {
    BbnInput stop     = BaseInput();
    stop.T9_final     = T9_stop;
    const BbnResult r = SolveBbn(stop);
    assert(std::fabs(r.Yp_mass - reference.Yp_mass) < 1.0e-5);
    assert(std::fabs(r.DoverH / reference.DoverH - 1.0) < 1.0e-3);
  }
}

/** The neutron lifetime is the single weak-sector input, and Y_p must respond to
 *  it: a longer-lived neutron leaves more neutrons at freeze-out and more helium.
 *  The published sensitivity is dY_p/dtau_n ~ 2e-4 per second. */
void test_neutron_lifetime_sensitivity() {
  BbnInput in                 = BaseInput();
  in.tau_n                    = 870.0;
  const BbnResult short_lived = SolveBbn(in);
  in.tau_n                    = 890.0;
  const BbnResult long_lived  = SolveBbn(in);

  assert(long_lived.Yp_mass > short_lived.Yp_mass);
  const double slope = (long_lived.Yp_mass - short_lived.Yp_mass) / 20.0;
  assert(slope > 1.0e-4 && slope < 3.0e-4);
}

/** The rho_extra seam must actually reach the expansion rate. Extra radiation
 *  makes the universe expand faster, freezing out more neutrons and raising Y_p,
 *  exactly as a positive DeltaN_eff does. */
void test_rho_extra_seam() {
  BbnInput plain          = BaseInput();
  const BbnResult without = SolveBbn(plain);

  BbnInput seeded = BaseInput();
  /* One extra neutrino-like species, expressed through the seam rather than
     through delta_neff, must reproduce delta_neff = 1 closely. */
  const double T_nu_initial = MeVFromT9(plain.T9_initial);
  seeded.rho_extra          = [T_nu_initial](double a) {
    return BbnPlasma::RhoNeutrino(T_nu_initial / a, 1.0);
  };
  const BbnResult via_seam = SolveBbn(seeded);

  BbnInput counted          = BaseInput();
  counted.delta_neff        = 1.0;
  const BbnResult via_delta = SolveBbn(counted);

  assert(via_seam.Yp_mass > without.Yp_mass);
  assert(std::fabs(via_seam.Yp_mass - via_delta.Yp_mass) < 1.0e-6);
}

/** The splined weak rates must not change the answer.
 *
 *  They are ~87% of a right-hand-side evaluation, so tabulating them is where
 *  nearly all of the runtime went; the whole saving is worthless if it moves the
 *  result. The reference path (weak_table_points = 0) evaluates the phase-space
 *  integrals directly at every step and is kept precisely so this comparison can
 *  be made rather than assumed.
 *
 *  MEASURED 2026-09-06: at the default 512 nodes the table shifts Y_p by 1.8e-10
 *  and D/H by 2.1e-9, while running 11x faster. Beyond ~512 nodes the difference
 *  stops falling because it is then the evolver's own rtol noise, not the
 *  interpolation. */
void test_weak_table_matches_exact_integrals() {
  BbnInput exact            = BaseInput();
  exact.weak_table_points   = 0;
  const BbnResult reference = SolveBbn(exact);

  BbnInput tabulated      = BaseInput();
  const BbnResult splined = SolveBbn(tabulated);

  assert(std::fabs(splined.Yp_mass - reference.Yp_mass) < 1.0e-8);
  assert(std::fabs(splined.DoverH / reference.DoverH - 1.0) < 1.0e-7);
  assert(std::fabs(splined.Li7overH / reference.Li7overH - 1.0) < 1.0e-5);

  /* A deliberately starved table must be measurably worse, or the bound above
     would be satisfied by any node count and would test nothing. */
  BbnInput starved          = BaseInput();
  starved.weak_table_points = 8;
  const BbnResult coarse    = SolveBbn(starved);
  assert(std::fabs(coarse.Yp_mass - reference.Yp_mass) > 1.0e-8);
}

/** Malformed inputs must fail rather than return an abundance.
 *
 *  SolveBbn is a library entry point, reachable without going through
 *  InputModule's parsing, so it cannot rely on that layer having checked
 *  anything. weak_table_points = 1 is the sharp one: the sampling grid divides by
 *  (points - 1). */
void test_input_guards() {
  const auto rejects = [](const BbnInput& in) {
    try {
      SolveBbn(in);
    }
    catch (const std::exception&) {
      return true;
    }
    return false;
  };

  BbnInput no_rates = BaseInput();
  no_rates.rates_file.clear();
  assert(rejects(no_rates));

  BbnInput inverted   = BaseInput();
  inverted.T9_initial = 0.01;
  inverted.T9_final   = 100.0;
  assert(rejects(inverted));

  for (double bad_tau : {0.0, -878.4}) {
    BbnInput in = BaseInput();
    in.tau_n    = bad_tau;
    assert(rejects(in));
  }

  BbnInput bad_omega = BaseInput();
  bad_omega.omega_b  = 0.0;
  assert(rejects(bad_omega));

  BbnInput bad_tol  = BaseInput();
  bad_tol.tolerance = 0.0;
  assert(rejects(bad_tol));

  for (int bad_points : {1, 2, 3, -8}) {
    BbnInput in          = BaseInput();
    in.weak_table_points = bad_points;
    assert(rejects(in));
  }

  /* 0 means "no table, exact integrals" and must remain legal. */
  BbnInput exact          = BaseInput();
  exact.weak_table_points = 0;
  assert(!rejects(exact));
}

}  // namespace

int main() {
  test_Yp_against_table();
  test_light_element_abundances();
  test_entropy_transfer();
  test_self_convergence();
  test_endpoint_independence();
  test_neutron_lifetime_sensitivity();
  test_rho_extra_seam();
  test_weak_table_matches_exact_integrals();
  test_input_guards();
  std::printf("bbn solver tests passed\n");
  return 0;
}
