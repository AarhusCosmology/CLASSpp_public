#pragma once
#include <memory>
#include <string>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "../species/species_build_context.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

class NCDMSpecies : public NCDMBaseSpecies {
 public:
  // New input path: parameters are read from PFC under <instance_name>.<field>.
  NCDMSpecies(FileContent* pfc,
              const std::string& instance_name,
              const NcdmSettings& settings,
              const background* pba,
              const BackgroundModule* bgm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  bool IsFreestreaming() const override {
    return true;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double DpDloga(const double* pvecback) const override {
    return pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_];
  }

  // ── Perturbations ────────────────────────────────────────────────────────

  // Layout-based scalar register (implementation writes both layout and legacy pv arrays).
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  /** Legacy scalar register: no-op — dual-written by the layout-based path above. */
  void RegisterPerturbationIndices(perturb_vector* /*pv*/,
                                   const precision* /*ppr*/,
                                   int& /*index_pt*/,
                                   const perturb_workspace* /*ppw*/,
                                   int /*gauge*/) override {}

  // Layout-based PerturbDerivs.
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  /** Legacy PerturbDerivs: no-op — superseded by layout-based path above. */
  void PerturbDerivs(double /*tau*/,
                     const double* /*y*/,
                     double* /*dy*/,
                     const perturb_parameters_and_workspace& /*ppaw*/) override {}

  // Layout-based FillSources.
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) override;

  /** Legacy FillSources: no-op — superseded by layout-based path above. */
  void FillSources(const double* /*y*/,
                   const double* /*dy*/,
                   PerturbSourceContext& /*ctx*/) override {}

  // Layout-based ApplyInitialConditions.
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  /** Legacy ApplyInitialConditions: no-op — superseded by layout-based path above. */
  void ApplyInitialConditions(double* /*y*/, const PerturbIcContext& /*ctx*/) override {}

  // Layout-based stress-energy observables. Legacy versions throw — all
  // callers go through the layout-based path.
  double Delta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  double Delta(const perturb_vector*,
               const double*,
               const double*,
               const perturb_workspace*) const override;

  double Theta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  double Theta(const perturb_vector*,
               const double*,
               const double*,
               const perturb_workspace*) const override;

  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;

  double DeltaP(const perturb_vector*,
                const double*,
                const double*,
                const perturb_workspace*) const override;

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  double RhoPlusPShear(const perturb_vector*,
                       const double*,
                       const double*,
                       const perturb_workspace*) const override;

  // FA-collapse switch hook.
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void WriteOutputColumns(PerturbColumnWriter&,
                          const PerturbationsModule&,
                          enum file_format,
                          TransferColumnSection) const override;
  void PrintVariables(PerturbColumnWriter&,
                      double,
                      const double*,
                      const PerturbationsModule&,
                      const perturb_workspace*) const override;
  void WriteTensorOutputColumnTitles(char* tensor_titles) const override;

  /** Per-species offset into the NCDM source-slot block, assigned in lex order
   *  by PerturbationsModule when source indices are allocated. Used to compute
   *  this species' (delta, theta) slot as `index_tp_delta_ncdm1_ + source_slot_`
   *  and to label output columns. -1 until the module assigns it. */
  void SetSourceSlot(int s) {
    source_slot_ = s;
  }
  int source_slot() const {
    return source_slot_;
  }

  int bg_number_index() const {
    return index_bg_number_;
  }
  int bg_pseudo_p_index() const {
    return index_bg_pseudo_p_;
  }

 protected:
  double GetDlnf0Dlnq(int iq, const double* /*pvecback*/) const override {
    return dlnf0_dlnq_[iq];
  }

  int source_slot_ = -1;  // assigned in lex order by PerturbationsModule
  const background* pba_;

  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;
  int index_pt_psi0_     = -1;
};
