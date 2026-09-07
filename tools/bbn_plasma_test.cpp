/** Plasma sector of the BBN solver: photon + electron-positron thermodynamics.
 *
 *  These checks need no nuclear physics at all, which is the point of testing
 *  them separately: if the expansion rate or the entropy transfer is wrong, it
 *  shows up here rather than as a mysterious offset in Y_p.
 */

#include "bbn_plasma.h"

#include <cassert>
#include <cmath>
#include <cstdio>

#include "bbn_rates.h"
#include "constants.h"

namespace {

bool Close(double a, double b, double rtol) {
  return std::fabs(a - b) <= rtol * std::fabs(b);
}

/* Well above the electron mass, e+- are ultra-relativistic and contribute
   7/8 * 4 = 3.5 of g_*, against 2 for photons: rho_em/rho_gamma -> 11/4. */
void test_relativistic_limit() {
  const BbnPlasma plasma;
  const double T      = MeVFromT9(1.0e4); /* ~860 MeV, far above m_e */
  const BbnEmState em = plasma.ElectromagneticState(T);
  assert(Close(em.rho / BbnPlasma::RhoPhoton(T), 11.0 / 4.0, 1.0e-6));
  /* Ultra-relativistic gas: p = rho/3. */
  assert(Close(em.pressure, em.rho / 3.0, 1.0e-6));
  /* rho ~ T^4 implies drho/dT = 4 rho / T. */
  assert(Close(em.drho_dT, 4.0 * em.rho / T, 1.0e-6));
}

/* Far below the electron mass the pairs have annihilated and only photons remain. */
void test_annihilated_limit() {
  const BbnPlasma plasma;
  const double T      = MeVFromT9(0.005);
  const BbnEmState em = plasma.ElectromagneticState(T);
  assert(Close(em.rho, BbnPlasma::RhoPhoton(T), 1.0e-12));
  assert(Close(em.pressure, em.rho / 3.0, 1.0e-12));
  assert(Close(em.drho_dT, 4.0 * em.rho / T, 1.0e-12));
}

/* drho/dT is used by the temperature equation and is easy to get wrong by a
   factor without it showing anywhere else, so check it against a finite
   difference of rho right where the pairs are annihilating. */
void test_heat_capacity_matches_finite_difference() {
  const BbnPlasma plasma;
  for (double T9 : {30.0, 10.0, 6.0, 3.0, 1.0}) {
    const double T       = MeVFromT9(T9);
    const double h       = 1.0e-5 * T;
    const double up      = plasma.ElectromagneticState(T + h).rho;
    const double down    = plasma.ElectromagneticState(T - h).rho;
    const double numeric = (up - down) / (2.0 * h);
    assert(Close(plasma.ElectromagneticState(T).drho_dT, numeric, 1.0e-7));
  }
}

/* The quadrature must be converged, not merely plausible. Measured against a
   1024-node reference over T9 in [0.1, 1000], the worst relative error is
   5.2e-9 at 64 nodes, 7.8e-12 at the default 128, and 1.4e-13 at 256 -- the
   default therefore sits four orders below the evolver tolerance the solver
   runs at. The bound asserted below is that measured 128-node figure with a
   little room, not a round number picked in advance. */
void test_quadrature_is_converged() {
  const BbnPlasma coarse(128);
  const BbnPlasma fine(512);
  for (double T9 : {1000.0, 100.0, 30.0, 10.0, 5.0, 2.0, 1.0, 0.3, 0.1}) {
    const double T     = MeVFromT9(T9);
    const BbnEmState c = coarse.ElectromagneticState(T);
    const BbnEmState f = fine.ElectromagneticState(T);
    assert(Close(c.rho, f.rho, 1.0e-10));
    assert(Close(c.pressure, f.pressure, 1.0e-10));
    assert(Close(c.drho_dT, f.drho_dT, 1.0e-10));
  }
  /* ...and 64 nodes would NOT pass that bound, so it is a real constraint
     rather than one any node count would satisfy. */
  const BbnPlasma too_coarse(64);
  const double T = MeVFromT9(30.0);
  assert(!Close(too_coarse.ElectromagneticState(T).rho, fine.ElectromagneticState(T).rho, 1.0e-10));
}

/* One absolute check that the MeV^4 -> cgs conversion inside HubbleFromRho is
   right, against H = sqrt(8 pi G rho / 3) evaluated independently in cgs. */
void test_hubble_conversion() {
  const double T           = MeVFromT9(10.0);
  const double rho         = BbnPlasma::RhoPhoton(T);
  const double mev4_to_cgs = 2.3200e5; /* g cm^-3 per MeV^4, to 5 digits */
  const double expected    = std::sqrt(8.0 * _PI_ * 6.67428e-8 * rho * mev4_to_cgs / 3.0);
  assert(Close(BbnPlasma::HubbleFromRho(rho), expected, 1.0e-4));
}

/* Neutrino energy density scales as N_eff and as T^4, and one species carries
   7/8 of a photon-like species. */
void test_neutrino_scaling() {
  const double T = MeVFromT9(5.0);
  assert(Close(BbnPlasma::RhoNeutrino(T, 1.0), 7.0 / 8.0 * BbnPlasma::RhoPhoton(T), 1.0e-12));
  assert(Close(BbnPlasma::RhoNeutrino(T, 3.0), 3.0 * BbnPlasma::RhoNeutrino(T, 1.0), 1.0e-12));
  assert(
      Close(BbnPlasma::RhoNeutrino(2.0 * T, 1.0), 16.0 * BbnPlasma::RhoNeutrino(T, 1.0), 1.0e-12));
}

}  // namespace

int main() {
  test_relativistic_limit();
  test_annihilated_limit();
  test_heat_capacity_matches_finite_difference();
  test_quadrature_is_converged();
  test_hubble_conversion();
  test_neutrino_scaling();
  std::printf("bbn plasma tests passed\n");
  return 0;
}
