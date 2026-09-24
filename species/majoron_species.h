#pragma once

#include <string_view>

#include "ultra_relativistic.h"

/**
 * Neutrinos coupled to a majoron (Escudero & Witte, arXiv:1909.04044 and 2103.03249; the
 * Majoron model of the H0 Olympics, arXiv:2107.10291 sec. 2.3.6). As in Escudero & Witte,
 * the three neutrinos and the eV-scale majoron phi are one massless fluid.
 *
 * Background: their fit to the solved thermodynamics (1909.04044, supplementary material),
 *
 *   rho / T_gamma^4 = alpha_1 - 10^alpha_2 / [1 + alpha_3 (1 + x) e^{-(x - 1)/alpha_4}],
 *
 * x = T_gamma/m_phi, with alpha_i(Gamma_eff). Majorons thermalize with the neutrinos near
 * T ~ m_phi and decay back into them, raising N_eff by up to ~0.11. The pressure follows
 * from continuity.
 *
 * Perturbations: UR's hierarchy (w = c_s^2 = 1/3, as they take it), with every multipole
 * l >= 2 relaxed by inverse decays and decays, dF_l/dtau -= a Gamma F_l, at the H0
 * Olympics' energy-transfer rate for zero neutrino chemical potential:
 *
 *   Gamma = Gamma_phi (m_phi/T_nu)^3 K_2(m_phi/T_nu),   Gamma_phi = lambda^2 m_phi / (16 pi),
 *   Gamma_eff = (lambda / 4e-14)^2 (0.1 eV / m_phi).
 *
 * The species follows UR's approximations. The radiation streaming approximation assumes
 * free streaming; it starts after recombination, when the interaction has switched off for
 * m_phi above ~0.3 eV. Extra free-streaming radiation (the paper's Delta N_ur) is N_ur.
 */
class MajoronSpecies : public UltraRelativisticSpecies {
 public:
  static constexpr std::string_view kTypeName = "majoron";

  MajoronSpecies(const background& pba, double m_phi, double Gamma_eff);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  /** Stiff while the interaction is on. */
  bool SupportsExplicitPerturbationEvolver() const override {
    return false;
  }

  double GetOmega0() const override;
  double GetRadiationOmega0() const override;
  /** The fit's density, for NEDE's trigger: w is not constant while the majorons decay. */
  double BackgroundDensityOverH0Sq(double a, double H0) const override;

  void RegisterBackgroundIndices(int& index_bg) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  double P(const double* pvecback) const override;
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;

  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;

 private:
  /** rho/T_gamma^4 of the fit and its first two derivatives in x = T_gamma/m_phi. */
  void Fit(double x, double* F, double* dF, double* d2F) const;
  /** a Gamma in 1/Mpc. */
  double CollisionRate(double a) const;

  double Gamma_phi_;       ///< majoron decay width, eV
  double x0_;              ///< T_gamma/m_phi today
  double H0_;              ///< 1/Mpc
  double Omega_nu_per_F_;  ///< density today per unit of F: Omega0_g 15/pi^2
  double alpha1_, ten_alpha2_, alpha3_, alpha4_;
};
