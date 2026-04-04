#pragma once

#include "base_species.h"

struct background;

/**
 * Ultra-relativistic neutrinos / relics (UR) species.
 *
 * Background:  rho_ur = Omega0_ur * H0^2 / a_rel^4,
 *              p_ur   = rho_ur / 3,
 *              dp/dloga = -4/3 * rho_ur.
 * Perturbations: full Boltzmann hierarchy in l (delta, theta, shear, l3..l_max).
 *   Supports UFA (fluid approximation) and RSA (radiation streaming).
 *   Non-standard ceff2_ur and cvis2_ur terms from ppt->three_ceff2_ur etc.
 */
class UltraRelativisticSpecies : public BaseSpecies {
public:
  explicit UltraRelativisticSpecies(const background& pba);

  // ── Background ─────────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void ComputeBackground(double a_rel, const double* pvecback_B,
                         double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  /** Registers delta_ur, theta_ur, shear_ur, and (when UFA is off) l3..l_max_ur. */
  void RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  double Delta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;

private:
  double Omega0_ur_;
  double H0_;
};
