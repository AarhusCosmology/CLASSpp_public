#pragma once

#include "base_species.h"

struct background;

/**
 * Decaying Cold Dark Matter (DCDM) species.
 *
 * Background:
 *   rho_dcdm is integrated via ODE; BackgroundDerivs is called from background_derivs_member().
 *   ComputeBackground reads rho_dcdm from pvecback_B.
 *
 * Perturbations (scalar):
 *   dy[delta_dcdm] = -(theta_dcdm + metric_continuity)
 *   dy[theta_dcdm] = -a_prime_over_a * theta_dcdm + metric_euler
 */
class DCDMSpecies : public BaseSpecies {
public:
  explicit DCDMSpecies(const background& pba);

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override { bgm_ = bgm; }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B,
                         double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy,
                        const double* pvecback) override;
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

  int bi_rho_index() const { return index_bi_rho_dcdm_; }

private:
  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;

  int index_bg_rho_dcdm_ = -1;
  int index_bi_rho_dcdm_ = -1;
};
