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
  NCDMSpecies(int ncdm_id, std::shared_ptr<NonColdDarkMatter> ncdm,
              const background* pba,
              const BackgroundModule* bgm)
    : BaseSpecies("NCDM_" + std::to_string(ncdm_id), EnergyType::Other),
      ncdm_id_(ncdm_id), ncdm_(std::move(ncdm)), pba_(pba), bgm_(bgm) {}

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override { bgm_ = bgm; }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double DpDloga(const double* pvecback) const override {
    // dp/dloga for NCDM component: (pseudo_p - 5*p)  [see CLASS IV paper eq. A6]
    return pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_];
  }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr, int& index_pt,
                                   const perturb_workspace* ppw, int gauge) override;
  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  double Delta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv, const double* y, const double* pvecback, const perturb_workspace* ppw) const override;

  int ncdm_id() const { return ncdm_id_; }
  int bg_number_index()   const { return index_bg_number_; }
  int bg_pseudo_p_index() const { return index_bg_pseudo_p_; }
  int bg_lnf_index()      const { return index_bg_lnf_decay_dr1_; }
  int bg_dlnfdlnq_index() const { return index_bg_dlnfdlnq_decay_; }
  int bg_dlnfdlnq_sep_index() const { return index_bg_dlnfdlnq_sep_; }

  int bi_lnf_index()      const { return index_bi_lnf_decay_dr1_; }
  int bi_dlnfdlnq_sep_index() const { return index_bi_dlnfdlnq_separate_decay_; }

private:
  int ncdm_id_;
  std::shared_ptr<NonColdDarkMatter> ncdm_;
  const background*     pba_;
  const BackgroundModule* bgm_;

  // Background indices (single slot each)
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
