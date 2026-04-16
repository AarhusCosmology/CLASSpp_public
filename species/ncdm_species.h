#pragma once
#include <memory>

#include "../species/base_species.h"
#include "../tools/non_cold_dark_matter.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Non-Cold Dark Matter: massive neutrinos and warm/hot dark matter.
 * Wraps the existing NonColdDarkMatter class which handles all N_ncdm species.
 */
class NCDMSpecies : public BaseSpecies {
 public:
  NCDMSpecies(int ncdm_id,
              std::shared_ptr<NonColdDarkMatter> ncdm,
              const background* pba,
              const BackgroundModule* bgm)
      : BaseSpecies("NCDM_" + std::to_string(ncdm_id), EnergyType::Other), ncdm_id_(ncdm_id),
        ncdm_(std::move(ncdm)), pba_(pba), bgm_(bgm) {}

  bool IsFreestreaming() const override {
    return true;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
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
  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
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

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  int ncdm_id() const {
    return ncdm_id_;
  }
  int bg_number_index() const {
    return index_bg_number_;
  }
  int bg_pseudo_p_index() const {
    return index_bg_pseudo_p_;
  }

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      enum file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

 private:
  int ncdm_id_;
  std::shared_ptr<NonColdDarkMatter> ncdm_;
  const background* pba_;
  const BackgroundModule* bgm_;

  // Background indices (single slot each)
  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;
  // index_bg_rho_ and index_bg_p_ are the base-class protected members

  // Perturbation indices
  int index_pt_psi0_ = -1;
};
