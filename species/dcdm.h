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
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;
    int idx_theta = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  DCDMSpecies(const background& pba,
              double omega0_dcdmdr,
              double Gamma_dcdm,
              double Omega_ini_dcdm);

  double GetOmega0() const override {
    return Omega0_dcdmdr_;
  }

  double Gamma_dcdm() const {
    return Gamma_dcdm_;
  }
  double Omega_ini_dcdm() const {
    return Omega_ini_dcdm_;
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;

  int bi_rho_index() const {
    return index_bi_rho_dcdm_;
  }
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double RhoDotOverRho(const double* pvecback, double a_prime_over_a) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;
  int transfer_delta_index() const {
    return index_tp_delta_;
  }
  int transfer_theta_index() const {
    return index_tp_theta_;
  }

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

 private:
  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
  double Omega0_dcdmdr_;
  double Gamma_dcdm_     = 0.;
  double Omega_ini_dcdm_ = 0.;

  // Integration indices
  int index_bi_rho_dcdm_ = -1;

  // Background indices
  int index_bg_rho_dcdm_ = -1;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
