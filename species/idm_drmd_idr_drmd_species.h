#pragma once
#include "background.h"
#include "composite_species.h"
#include "idm_drmd.h"
#include "idr_drmd.h"
#include "species_build_context.h"

class IDM_DRMD_IDR_DRMD_Species : public CompositeSpecies {
 public:
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    IDM_DRMDSpecies::PerturbLayout idm_drmd;
    IDR_DRMDSpecies::PerturbLayout idr_drmd;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  explicit IDM_DRMD_IDR_DRMD_Species(const background& pba);

  IDM_DRMDSpecies& idm_drmd() {
    return *idm_drmd_;
  }
  IDR_DRMDSpecies& idr_drmd() {
    return *idr_drmd_;
  }
  const IDM_DRMDSpecies& idm_drmd() const {
    return *idm_drmd_;
  }
  const IDR_DRMDSpecies& idr_drmd() const {
    return *idr_drmd_;
  }

  /** Per-sub-species presence, captured at construction (== pba->has_* at build
   *  time). The composite is created on has_idm_drmd || has_idr_drmd, so callers
   *  needing the individual flags must use these rather than count(...). */
  bool has_idm_drmd() const {
    return has_idm_drmd_;
  }
  bool has_idr_drmd() const {
    return has_idr_drmd_;
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

  /** Matter (IDM_DRMD only) contribution. IDR_DRMD is radiation and contributes 0. */
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

  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
    CompositeSpecies::SetBackgroundModule(bgm);
  }

 private:
  IDM_DRMDSpecies* idm_drmd_ = nullptr;
  IDR_DRMDSpecies* idr_drmd_ = nullptr;
  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
  bool has_idm_drmd_           = false;
  bool has_idr_drmd_           = false;
};
