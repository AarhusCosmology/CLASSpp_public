#pragma once
#include <string_view>

#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"
#include "species_build_context.h"

class BackgroundModule;

/**
 * Scalar field: phi and phi' integrated via ODE.
 * Potential V(phi) = exp(-lambda*phi) * ((phi-B)^alpha + A).
 */
class ScalarFieldSpecies : public BaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "scalar_field";

  ScalarFieldSpecies(const background& pba,
                     double omega0_scf,
                     std::vector<double> scf_parameters,
                     int scf_tuning_index,
                     bool attractor_ic_scf,
                     double phi_ini_scf,
                     double phi_prime_ini_scf);

  double GetOmega0() const override {
    return Omega0_scf_;
  }

  const std::vector<double>& scf_parameters() const {
    return scf_parameters_;
  }
  int scf_tuning_index() const {
    return scf_tuning_index_;
  }
  bool attractor_ic_scf() const {
    return attractor_ic_scf_;
  }
  double phi_ini_scf() const {
    return phi_ini_scf_;
  }
  double phi_prime_ini_scf() const {
    return phi_prime_ini_scf_;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_phi_scf_       = index_bg++;
    index_bg_phi_prime_scf_ = index_bg++;
    index_bg_V_scf_         = index_bg++;
    index_bg_dV_scf_        = index_bg++;
    index_bg_ddV_scf_       = index_bg++;
    index_bg_rho_           = index_bg++;
    index_bg_p_             = index_bg++;
    index_bg_p_prime_scf_   = index_bg++;
  }

  void RegisterIntegrationIndices(int& index_bi) override {
    index_bi_phi_scf_       = index_bi++;
    index_bi_phi_prime_scf_ = index_bi++;
  }

  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;

  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;

  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Perturbations ────────────────────────────────────────────────────────
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_phi = -1;  // delta_phi
    // Synchronous: delta_phi'. Newtonian: q = delta_phi' - phi'_bg (psi + 3 phi).
    int idx_phi_prime = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;
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
                     const perturb_parameters_and_workspace& ppaw) const override;

  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;

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

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Shooting hooks ───────────────────────────────────────────────────────
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

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

  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;
  void FinalizeBackground(double a, double H, const double* pvecback_B, double* pvecback) override;

  int bi_phi_index() const {
    return index_bi_phi_scf_;
  }
  int bi_phi_prime_index() const {
    return index_bi_phi_prime_scf_;
  }

 private:
  double V_scf(double phi) const;
  double dV_scf(double phi) const;
  double ddV_scf(double phi) const;

  ShootingTarget shooting_target_{};  // unknown_param empty => no shooting target
  bool needs_shooting_ = false;       // true iff the direct unknown was absent (we guessed)

  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
  int index_tp_delta_          = -1;  // #309 transfer-source slot
  int index_tp_theta_          = -1;

  int index_bg_phi_scf_       = -1;
  int index_bg_phi_prime_scf_ = -1;
  int index_bg_V_scf_         = -1;
  int index_bg_dV_scf_        = -1;
  int index_bg_ddV_scf_       = -1;
  int index_bg_p_prime_scf_   = -1;
  int index_bi_phi_scf_       = -1;
  int index_bi_phi_prime_scf_ = -1;
  double Omega0_scf_          = 0.;
  std::vector<double> scf_parameters_;
  int scf_tuning_index_     = 0;
  bool attractor_ic_scf_    = true;
  double phi_ini_scf_       = 1.;
  double phi_prime_ini_scf_ = 1.;
};
