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
  /** The parent is a momentum-resolved NCDM hierarchy. */
  bool HasNcdm() const override {
    return true;
  }

  // ── PerturbLayout ──────────────────────────────────────────────────────────
  // Uses the generic CompositeSpecies::PerturbLayout (one owning sub-layout per
  // child in children_ order): the base-class child loops — including the
  // TallyStressEnergy/DelegateTally hot path — index child_layouts, so every
  // composite layout must carry it (#358). Typed views below.
  enum ChildIndex { kDncdm = 0, kDr = 1 };  // children_ order, set in the ctor
  static const NCDMBaseSpecies::PerturbLayout& dncdm_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kDncdm]);
  }
  static const DarkRadiationSpecies::PerturbLayout& dr_layout(
      const BaseSpecies::PerturbLayout& my) {
    return static_cast<const DarkRadiationSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kDr]);
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
  // PerturbDerivs (children + AddCouplingDerivs), ICs, stress-energy and
  // approximation-switch copies all use the generic CompositeSpecies child
  // loops. Registration stays overridden for its DR-first slot ordering.
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void PerturbTensorDerivs(const BaseSpecies::PerturbLayout& layout,
                           double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) const override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  DNCDMSpecies* dncdm_         = nullptr;
  DarkRadiationSpecies* dr_sp_ = nullptr;
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;
};
