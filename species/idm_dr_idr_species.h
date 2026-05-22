#pragma once
#include "background.h"
#include "composite_species.h"
#include "idm_dr.h"
#include "idr.h"
#include "species_build_context.h"

/**
 * IDM_DR_IDR_Species: composite for interacting dark matter + interacting dark radiation.
 *
 * Children handle free-streaming terms; this composite's AddCouplingDerivs
 * adds the momentum-exchange and TCA terms that couple IDM_DR to IDR.
 */
class IDM_DR_IDR_Species : public CompositeSpecies {
 public:
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    IDM_DRSpecies::PerturbLayout idm_dr;
    IDRSpecies::PerturbLayout idr;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  explicit IDM_DR_IDR_Species(const background& pba);

  IDM_DRSpecies& idm_dr() {
    return *idm_dr_;
  }
  IDRSpecies& idr() {
    return *idr_;
  }
  const IDM_DRSpecies& idm_dr() const {
    return *idm_dr_;
  }
  const IDRSpecies& idr() const {
    return *idr_;
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

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

  /** Matter (IDM_DR only) contribution. IDR is radiation and contributes 0. */
  double MatterRhoDelta(const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const override;
  double MatterRhoPlusPTheta(const perturb_vector* pv,
                             const double* y,
                             const double* pvecback,
                             const perturb_workspace* ppw) const override;

  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
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

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

 private:
  IDM_DRSpecies* idm_dr_ = nullptr;
  IDRSpecies* idr_       = nullptr;
  const background& pba_;
};
