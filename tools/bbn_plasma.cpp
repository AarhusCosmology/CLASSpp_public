#include "bbn_plasma.h"

#include <cmath>

#include "bbn_nuclides.h"
#include "constants.h"
#include "quadrature.h"

namespace {

/** hbar*c in MeV cm, and the resulting (1 MeV)^3 in cm^-3. */
constexpr double kHbarCMeVcm = 1.973269804e-11;
const double kMeV3ToInvCm3   = 1.0 / (kHbarCMeVcm * kHbarCMeVcm * kHbarCMeVcm);

/** 1 MeV/c^2 in grams: 1e6 eV in joules, divided by c^2, times 1000 g/kg. */
const double kMeVToGram = 1.0e6 * _eV_ / (_c_ * _c_) * 1.0e3;

/** 1 MeV^4 in g/cm^3 -- the single place the natural-unit plasma meets cgs. */
const double kMeV4ToGramPerCm3 = kMeV3ToInvCm3 * kMeVToGram;

/** Newton's constant in cgs (cm^3 g^-1 s^-2), from the SI value CLASS uses. */
const double kNewtonCgs = _G_ * 1.0e3;

/** Riemann zeta(3), for the photon number density. */
constexpr double kZeta3 = 1.2020569031595942854;

/** Fermi-Dirac occupation 1/(e^x + 1), written so that large x underflows to
 *  zero rather than overflowing. The argument is always non-negative here (it is
 *  E/T for a massive particle), so the exp(-x) form is unconditionally safe --
 *  which matters because these translation units are built with -ffast-math,
 *  where a materialised infinity is undefined rather than merely large. */
inline double FermiDirac(double x) {
  const double e = std::exp(-x);
  return e / (1.0 + e);
}

/** -df/dx = e^x/(e^x+1)^2, in the same overflow-free form. */
inline double FermiDiracSlope(double x) {
  const double e = std::exp(-x);
  const double d = 1.0 + e;
  return e / (d * d);
}

}  // namespace

BbnPlasma::BbnPlasma(int nodes, double u_max) : u_(nodes), w_(nodes) {
  std::vector<double> mu(nodes), weight(nodes);
  quadrature_gauss_legendre(mu.data(), weight.data(), nodes, 1.0e-14);

  /* Map the [-1,1] rule onto [0, u_max]. */
  for (int i = 0; i < nodes; ++i) {
    u_[i] = 0.5 * u_max * (mu[i] + 1.0);
    w_[i] = 0.5 * u_max * weight[i];
  }
}

double BbnPlasma::RhoPhoton(double T_mev) {
  return (_PI_ * _PI_ / 15.0) * T_mev * T_mev * T_mev * T_mev;
}

double BbnPlasma::NumberDensityPhoton(double T_mev) {
  const double n_natural = 2.0 * kZeta3 / (_PI_ * _PI_) * T_mev * T_mev * T_mev;
  return n_natural * kMeV3ToInvCm3;
}

double BbnPlasma::RhoNeutrino(double T_nu_mev, double n_eff) {
  /* Per effective species (neutrino plus antineutrino, one helicity each):
     rho = 7/8 * 2 * pi^2/30 * T^4 = 7/8 * pi^2/15 * T^4. */
  const double t2 = T_nu_mev * T_nu_mev;
  return n_eff * (7.0 / 8.0) * (_PI_ * _PI_ / 15.0) * t2 * t2;
}

double BbnPlasma::HubbleFromRho(double rho_mev4) {
  const double rho_cgs = rho_mev4 * kMeV4ToGramPerCm3;
  return std::sqrt(8.0 * _PI_ * kNewtonCgs * rho_cgs / 3.0);
}

double BbnPlasma::InvCm3ToMeV3() {
  return 1.0 / kMeV3ToInvCm3;
}

BbnEmState BbnPlasma::ElectromagneticState(double T_mev) const {
  const double t2 = T_mev * T_mev;
  const double t3 = t2 * T_mev;
  const double t4 = t3 * T_mev;

  /* Photons: analytic. */
  BbnEmState state;
  state.rho      = (_PI_ * _PI_ / 15.0) * t4;
  state.pressure = state.rho / 3.0;
  state.drho_dT  = 4.0 * state.rho / T_mev;

  /* Electron-positron pairs, four internal degrees of freedom (two spins, two
     charges), zero chemical potential:

        rho     = (2/pi^2) T^4 Int du u^2 eps       f(eps)
        p       = (2/pi^2) T^4 Int du u^4/(3 eps)   f(eps)
        drho/dT = (2/pi^2) T^3 Int du u^2 eps^2   (-f'(eps))

     with u = p/T and eps = sqrt(u^2 + z^2), z = m_e/T. All three share eps and
     the exponential, so they are accumulated in one pass. */
  const double z  = kElectronMassMeV / T_mev;
  const double z2 = z * z;

  double integral_rho = 0.0;
  double integral_p   = 0.0;
  double integral_c   = 0.0;

  for (std::size_t i = 0; i < u_.size(); ++i) {
    const double u   = u_[i];
    const double u2  = u * u;
    const double eps = std::sqrt(u2 + z2);
    const double f   = FermiDirac(eps);
    const double fp  = FermiDiracSlope(eps);

    integral_rho += w_[i] * u2 * eps * f;
    integral_p   += w_[i] * u2 * u2 / (3.0 * eps) * f;
    integral_c   += w_[i] * u2 * eps * eps * fp;
  }

  const double prefactor  = 2.0 / (_PI_ * _PI_);
  state.rho              += prefactor * t4 * integral_rho;
  state.pressure         += prefactor * t4 * integral_p;
  state.drho_dT          += prefactor * t3 * integral_c;

  return state;
}
