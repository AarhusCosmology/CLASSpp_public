#pragma once
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"
#include "../tools/non_cold_dark_matter.h"
#include <memory>

class BackgroundModule;

/**
 * Non-Cold Dark Matter: massive neutrinos and warm/hot dark matter.
 * Wraps the existing NonColdDarkMatter class which handles all N_ncdm species.
 */
class NCDMSpecies : public BaseSpecies {
public:
  NCDMSpecies(std::shared_ptr<NonColdDarkMatter> ncdm,
              const background* pba,
              const BackgroundModule* bgm)
    : BaseSpecies("NCDM", EnergyType::Other),
      ncdm_(std::move(ncdm)), pba_(pba), bgm_(bgm) {}

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override { bgm_ = bgm; }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    double rho = 0.;
    for (int n = 0; n < pba_->N_ncdm; ++n) rho += pvecback[index_bg_rho_ + n];
    return rho;
  }
  double P(const double* pvecback) const override {
    double p = 0.;
    for (int n = 0; n < pba_->N_ncdm; ++n) p += pvecback[index_bg_p_ + n];
    return p;
  }
  double DpDloga(const double* pvecback) const override {
    // dp/dloga for NCDM: sum_n (pseudo_p_n - 5*p_n)  [see CLASS IV paper eq. A6]
    double dp = 0.;
    for (int n = 0; n < pba_->N_ncdm; ++n)
      dp += pvecback[index_bg_pseudo_p_ + n] - 5. * pvecback[index_bg_p_ + n];
    return dp;
  }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                   const perturb_workspace* ppw, int gauge) override;
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  double Delta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;

  int bg_number_index()   const { return index_bg_number_; }
  int bg_pseudo_p_index() const { return index_bg_pseudo_p_; }

private:
  std::shared_ptr<NonColdDarkMatter> ncdm_;
  const background*     pba_;
  const BackgroundModule* bgm_;

  // Background indices (N_ncdm slots each)
  int index_bg_number_    = -1;
  int index_bg_pseudo_p_  = -1;
  // index_bg_rho_ and index_bg_p_ are the base-class protected members

  // Integration indices for decaying NCDM distribution function
  int index_bi_lnf_decay_dr1_           = -1;
  int index_bi_dlnfdlnq_separate_decay_ = -1;

  // Background indices for decay-dr lnf/dlnfdlnq slots
  int index_bg_lnf_decay_dr1_    = -1;
  int index_bg_dlnfdlnq_decay_   = -1;
  int index_bg_dlnfdlnq_sep_     = -1;

  // Perturbation indices
  int index_pt_psi0_ = -1;
};
