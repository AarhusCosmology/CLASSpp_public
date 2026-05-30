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

  IDM_DR_IDR_Species(
      const background& pba, double omega0_idm_dr, double omega0_idr, double T_idr, int l_max_idr);

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

  /** Per-sub-species presence, captured at construction (== pba->has_* at build
   *  time). The composite is created on has_idm_dr || has_idr, so callers needing
   *  the individual flags must use these rather than count("IDM_DR_IDR"). */
  bool has_idm_dr() const {
    return has_idm_dr_;
  }
  bool has_idr() const {
    return has_idr_;
  }

  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
    CompositeSpecies::SetPerturbs(ppt);
  }

  /** IDR l>=3 multipoles not needed in sources when rsa_idr is off, idr is
   *  free-streaming, and tca_idm_dr is off. */
  void MarkUsedInSources(const BaseSpecies::PerturbLayout& layout,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;

  std::optional<double> GetParam(const std::string& name) const override;

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
  bool has_idm_dr_     = false;
  bool has_idr_        = false;
  const perturbs* ppt_ = nullptr;
};
