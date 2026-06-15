#pragma once
#include "background.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dcdm.h"
#include "species_build_context.h"

class BackgroundModule;

/**
 * DCDM_DR_Species: composite for Decaying Cold Dark Matter + its decay radiation.
 *
 * EnergyType::Other so the background loop splits rho correctly:
 *   rho_m += rho - 3p  (= rho_dcdm)
 *   rho_r += 3p        (= rho_dr)
 *
 * BackgroundDerivs override adds the DCDM->DR decay source on top of children.
 * AddCouplingDerivs adds the DR Boltzmann l=0,1 coupling source terms from DCDM.
 */
class DCDM_DR_Species : public CompositeSpecies {
 public:
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    DCDMSpecies::PerturbLayout dcdm;
    DarkRadiationSpecies::PerturbLayout dr;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  DCDM_DR_Species(const background* pba,
                  const BackgroundModule* bgm,
                  double omega0_dcdmdr,
                  double Gamma_dcdm,
                  double Omega_ini_dcdm);

  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;

  // Override to add DCDM->DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

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

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  double Delta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  double Theta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  /** Matter (cold-DCDM only) contribution. DR is radiation and contributes 0. */
  double MatterRhoDelta(const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const override;
  double MatterRhoPlusPTheta(const perturb_vector* pv,
                             const double* y,
                             const double* pvecback,
                             const perturb_workspace* ppw) const override;

  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      enum file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  // Typed accessors so callers can capture child indices
  DCDMSpecies& dcdm() {
    return *dcdm_;
  }
  DarkRadiationSpecies& dr() {
    return *dr_sp_;
  }
  const DCDMSpecies& dcdm() const {
    return *dcdm_;
  }
  const DarkRadiationSpecies& dr() const {
    return *dr_sp_;
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

 private:
  DCDMSpecies* dcdm_           = nullptr;  // non-owning pointer into children_
  DarkRadiationSpecies* dr_sp_ = nullptr;  // non-owning pointer into children_
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;

  ShootingTarget shooting_target_{};  // unknown_param empty => no shooting target
  bool needs_shooting_ = false;       // true iff the direct unknown was absent (we guessed)
};
