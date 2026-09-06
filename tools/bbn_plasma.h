#pragma once

/** @file bbn_plasma.h Thermodynamics of the BBN plasma and the expansion rate.
 *
 *  Design: docs/superpowers/specs/2026-09-06-bbn-nuclear-network-design.md section 3
 *
 *  CLASS's own background cannot be used here: it contains no electrons or
 *  positrons and pins T_gamma = T_cmb/a with a fixed N_eff at all times, so its H
 *  at T ~ MeV is wrong by a large factor. This is a self-contained plasma, with a
 *  rho_extra seam (see BbnExpansion) for anything the caller wants to add.
 *
 *  Units in this file: temperatures and energy densities in MeV and MeV^4,
 *  number densities in cm^-3, H in s^-1. The MeV^4 -> cgs conversion happens in
 *  exactly one place, HubbleFromRho().
 */

#include <vector>

/** Thermodynamic state of the sector that is thermally COUPLED to the photons:
 *  photons plus electron-positron pairs. Baryons are added by the solver, which
 *  is the only thing that knows the composition. Neutrinos are decoupled and are
 *  deliberately not here -- they contribute to H but not to this energy balance. */
struct BbnEmState {
  double rho;      /**< MeV^4 */
  double pressure; /**< MeV^4 */
  double drho_dT;  /**< MeV^3 */
};

/** Photon + electron-positron thermodynamics.
 *
 *  The electron-positron integrals are done by Gauss-Legendre quadrature on a
 *  fixed grid built once. The integrands are analytic in the integration
 *  variable, so the rule converges geometrically and a modest node count reaches
 *  machine precision; bbn_plasma_test measures that rather than assuming it. */
class BbnPlasma {
 public:
  /** @param nodes    quadrature nodes for the e+- integrals
   *  @param u_max    upper limit in units of p/T; the integrand falls as e^{-u} */
  explicit BbnPlasma(int nodes = 128, double u_max = 60.0);

  /** Photons plus electron-positron pairs at photon temperature T (MeV).
   *
   *  The pairs are treated with zero chemical potential. The true value is set by
   *  charge neutrality and is of order the baryon-to-photon ratio, ~1e-9, so it
   *  is negligible in every quantity computed here. */
  BbnEmState ElectromagneticState(double T_mev) const;

  /** Energy density of Nnu effective neutrino species at temperature T_nu (MeV).
   *
   *  Neutrinos are treated as INSTANTANEOUSLY decoupled: T_nu is propagated as
   *  1/a by the caller, never re-thermalized. Incomplete decoupling -- the
   *  ~0.044 of N_eff that comes from e+- entropy leaking into neutrinos during
   *  annihilation -- is therefore not modelled. That is an energy TRANSFER rather
   *  than an addition, so it barely changes H; what it does change is the final
   *  photon-to-neutrino temperature ratio, and hence eta. */
  static double RhoNeutrino(double T_nu_mev, double n_eff);

  /** Photon energy density alone, MeV^4. */
  static double RhoPhoton(double T_mev);

  /** Photon number density, cm^-3. */
  static double NumberDensityPhoton(double T_mev);

  /** H in s^-1 from a total energy density in MeV^4. */
  static double HubbleFromRho(double rho_mev4);

  /** Conversion used by the solver to express baryon energy densities. */
  static double InvCm3ToMeV3();

 private:
  std::vector<double> u_; /**< quadrature abscissae in p/T */
  std::vector<double> w_; /**< quadrature weights */
};
