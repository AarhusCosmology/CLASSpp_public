#pragma once
#include <memory>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dncdm_species.h"
#include "parser.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"

class BackgroundModule;

/**
 * DNCDM_DR_Species: composite for one flavor of Decaying Non-Cold Dark Matter + its decay radiation.
 */
class DNCDM_DR_Species : public CompositeSpecies {
 public:
  // ── PerturbLayout ──────────────────────────────────────────────────────────
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    NCDMBaseSpecies::PerturbLayout dncdm;
    DarkRadiationSpecies::PerturbLayout dr;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  // Takes ownership of a pre-built DNCDMSpecies (from DNCDMSpecies::CreateAll)
  DNCDM_DR_Species(std::unique_ptr<DNCDMSpecies> dncdm,
                   const background* pba,
                   const BackgroundModule* bgm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Shooter hooks ──────────────────────────────────────────────────────────
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;

  // Override to add DNCDM->DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    dncdm_->WriteBackgroundColumnTitles(w);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    dncdm_->WriteBackgroundData(pvecback, w);
  }

  /** Total density fraction of the whole decaying sector today (matter + decay
   *  radiation) = the closure-reserved Omega_dncdmdr. Mirrors DCDMSpecies::GetOmega0().
   *  Surfaces the pinned/shot combined rather than summing children, because at
   *  Pass-1 closure time the emergent DR is not yet integrated. */
  double GetOmega0() const override;

  DNCDMSpecies& dncdm() {
    return *dncdm_;
  }
  const DNCDMSpecies& dncdm() const {
    return *dncdm_;
  }
  DarkRadiationSpecies& dr() {
    return *dr_sp_;
  }
  const DarkRadiationSpecies& dr() const {
    return *dr_sp_;
  }

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

  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& layout,
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

  double DeltaRho(const BaseSpecies::PerturbLayout& layout,
                  const perturb_vector* pv,
                  const double* y,
                  const double* pvecback,
                  const perturb_workspace* ppw) const override;

  double RhoPlusPTheta(const BaseSpecies::PerturbLayout& layout,
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

  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

 private:
  void AddCouplingDerivs(const PerturbLayout& my,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw);

  DNCDMSpecies* dncdm_         = nullptr;
  DarkRadiationSpecies* dr_sp_ = nullptr;
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;
};
