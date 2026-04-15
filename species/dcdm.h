#pragma once
#include <memory>

#include "../species/base_species.h"
#include "background.h"

class BackgroundModule;

/**
 * Decaying Cold Dark Matter (DCDM).
 * Background density stored in the ODE integration vector.
 */
class DCDMSpecies : public BaseSpecies {
 public:
  explicit DCDMSpecies(const background& pba);

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;

  int bi_rho_index() const {
    return index_bi_rho_dcdm_;
  }
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

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
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

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
  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;

  // Integration indices
  int index_bi_rho_dcdm_ = -1;

  // Background indices
  int index_bg_rho_dcdm_ = -1;
};
