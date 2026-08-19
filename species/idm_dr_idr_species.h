#pragma once
#include <string_view>

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
  static constexpr std::string_view kTypeName = "idm_dr_idr";

  // Uses the generic CompositeSpecies::PerturbLayout (one owning sub-layout per
  // child in children_ order): the base-class child loops — including the
  // TallyStressEnergy/DelegateTally hot path — index child_layouts, so every
  // composite layout must carry it (#358). Typed views below.
  enum ChildIndex { kIdmDr = 0, kIdr = 1 };  // children_ order, set in the ctor
  static const IDM_DRSpecies::PerturbLayout& idm_dr_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const IDM_DRSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kIdmDr]);
  }
  static const IDRSpecies::PerturbLayout& idr_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const IDRSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kIdr]);
  }

  IDM_DR_IDR_Species(const background& pba,
                     double omega0_idm_dr,
                     double omega0_idr,
                     double T_idr,
                     int l_max_idr,
                     double a_idm_dr,
                     double nindex_idm_dr,
                     double m_idm,
                     double b_idr,
                     int idr_nature,
                     std::vector<double> alpha_idm_dr,
                     std::vector<double> beta_idr);

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

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

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
  // Registration, derivs (children + AddCouplingDerivs), stress-energy,
  // sync->Newtonian and approximation-switch copies all use the generic
  // CompositeSpecies child loops. Only composite-specific logic is overridden.
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
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

  void PrintVariables(PerturbColumnWriter& writer,
                      const BaseSpecies::PerturbLayout* base,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  IDM_DRSpecies* idm_dr_ = nullptr;
  IDRSpecies* idr_       = nullptr;
  const background& pba_;
  bool has_idm_dr_     = false;
  bool has_idr_        = false;
  const perturbs* ppt_ = nullptr;

  int index_tp_delta_idm_dr_ = -1;  // #309 transfer-source slots (composite-owned)
  int index_tp_theta_idm_dr_ = -1;
  int index_tp_delta_idr_    = -1;
  int index_tp_theta_idr_    = -1;
};
