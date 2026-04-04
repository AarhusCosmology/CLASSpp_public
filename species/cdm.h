#pragma once

#include "base_species.h"

struct background;

/**
 * Cold Dark Matter (CDM) species.
 *
 * Background:  rho_cdm = Omega0_cdm * H0^2 / a_rel^3,  p = 0
 * Perturbations: pressureless fluid, no anisotropic stress.
 *   Newtonian gauge: evolves delta_cdm and theta_cdm.
 *   Synchronous gauge: evolves delta_cdm only (theta = 0 by gauge choice).
 */
class CDMSpecies : public BaseSpecies {
public:
  explicit CDMSpecies(const background& pba);

  // ── Background ─────────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void ComputeBackground(double a_rel, const double* pvecback_B,
                         double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
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
  double Omega0_cdm_;
  double H0_;
  int index_bg_rho_cdm_ = -1;
};

