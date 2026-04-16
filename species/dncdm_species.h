#pragma once
#include <memory>

#include "../species/base_species.h"
#include "../tools/non_cold_dark_matter.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;
class BackgroundColumnWriter;

/**
 * Decaying Non-Cold Dark Matter (DNCDM).
 * Handles the background and perturbation evolution for a decaying NCDM species.
 */
class DNCDMSpecies : public BaseSpecies {
 public:
  DNCDMSpecies(int ncdm_id,
               std::shared_ptr<NonColdDarkMatter> ncdm,
               const background* pba,
               const BackgroundModule* bgm)
      : BaseSpecies("DNCDM_" + std::to_string(ncdm_id), EnergyType::Other), ncdm_id_(ncdm_id),
        ncdm_(std::move(ncdm)), pba_(pba), bgm_(bgm) {}

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double DpDloga(const double* pvecback) const override {
    return pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_];
  }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  double Delta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  bool IsFreestreaming() const override {
    return true;
  }
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  int ncdm_id() const {
    return ncdm_id_;
  }
  const NonColdDarkMatter& ncdm() const {
    return *ncdm_;
  }
  int bg_number_index() const {
    return index_bg_number_;
  }
  int bg_pseudo_p_index() const {
    return index_bg_pseudo_p_;
  }
  int bg_lnf_index() const {
    return index_bg_lnf_decay_dr1_;
  }
  int bg_dlnfdlnq_index() const {
    return index_bg_dlnfdlnq_decay_;
  }
  int bg_dlnfdlnq_sep_index() const {
    return index_bg_dlnfdlnq_sep_;
  }

  int bi_lnf_index() const {
    return index_bi_lnf_decay_dr1_;
  }
  int bi_dlnfdlnq_sep_index() const {
    return index_bi_dlnfdlnq_separate_decay_;
  }

 private:
  int ncdm_id_;
  std::shared_ptr<NonColdDarkMatter> ncdm_;
  const background* pba_;
  const BackgroundModule* bgm_;

  // Background indices (single slot each)
  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;

  // Integration indices for decaying NCDM distribution function
  int index_bi_lnf_decay_dr1_           = -1;
  int index_bi_dlnfdlnq_separate_decay_ = -1;

  // Background indices for decay-dr lnf/dlnfdlnq slots
  int index_bg_lnf_decay_dr1_  = -1;
  int index_bg_dlnfdlnq_decay_ = -1;
  int index_bg_dlnfdlnq_sep_   = -1;

  // Perturbation indices
  int index_pt_psi0_ = -1;
};
