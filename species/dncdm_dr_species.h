#pragma once
#include <memory>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "dncdm_decay_radiation_species.h"
#include "dncdm_species.h"
#include "parser.h"
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
    DNCDM_DecayRadiationSpecies::PerturbLayout dr;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  // Takes ownership of a pre-built DNCDMSpecies (from DNCDMSpecies::CreateAll)
  DNCDM_DR_Species(std::unique_ptr<DNCDMSpecies> dncdm,
                   const background* pba,
                   const BackgroundModule* bgm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;

  // Override to add DNCDM->DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    dncdm_->WriteBackgroundColumnTitles(w);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    dncdm_->WriteBackgroundData(pvecback, w);
  }

  DNCDMSpecies& dncdm() {
    return *dncdm_;
  }
  const DNCDMSpecies& dncdm() const {
    return *dncdm_;
  }
  DNCDM_DecayRadiationSpecies& dr() {
    return *dr_sp_;
  }
  const DNCDM_DecayRadiationSpecies& dr() const {
    return *dr_sp_;
  }

  // ── Perturbations ──────────────────────────────────────────────────────────
  // NOTE: RegisterPerturbationIndices is intentionally NOT overridden here.
  // The module calls each child's Register directly (DR child first, then DNCDM child via
  // the NCDM block) to preserve the pv->y slot ordering established before Task 22.
  // See perturbations_module.cpp around the DNCDM_DR block for the direct-call site.

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

 private:
  void AddCouplingDerivs(const PerturbLayout& my,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw);

  DNCDMSpecies* dncdm_                = nullptr;
  DNCDM_DecayRadiationSpecies* dr_sp_ = nullptr;
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;
};
