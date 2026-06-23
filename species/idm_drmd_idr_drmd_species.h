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

  IDM_DRMD_IDR_DRMD_Species(const background& pba,
                            double omega0_idm_drmd,
                            double omega0_idr_drmd,
                            double f_idm_drmd,
                            double G_over_aH_drmd,
                            double delta_Neff_drmd,
                            double z_stop);

  double f_idm_drmd() const {
    return f_idm_drmd_;
  }
  double G_over_aH_drmd() const {
    return G_over_aH_drmd_;
  }
  double delta_Neff_drmd() const {
    return delta_Neff_drmd_;
  }
  double z_stop() const {
    return z_stop_;
  }

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

  int bg_G_over_aH_index() const {
    return index_bg_G_over_aH_drmd_;
  }
  double f_idr_drmd() const {
    return f_idr_drmd_;
  }
  double Gamma0_drmd_ic() const {
    return Gamma0_drmd_ic_;
  }
  double z_dec_drmd() const {
    return z_dec_drmd_;
  }

  void RegisterBackgroundIndices(int& index_bg) override;
  void ComputeIdmDrmd(
      double a, double rho_idm_over_rho_idr, double* Rint, double* csp2, double* Gint) const;
  void InitializeDrmdBackground(double rho_tot, double H, double a, const double* pvecback);
  void FinalizeBackground(double a, double H, const double* pvecback_B, double* pvecback) override;
  void ProcessBackgroundTable(const double* background_table,
                              int n_rows,
                              int row_stride,
                              const double* z_table) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

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
                     const perturb_parameters_and_workspace& ppaw) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

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
  double f_idm_drmd_           = 0.;
  double G_over_aH_drmd_       = 0.;
  double delta_Neff_drmd_      = 0.;
  double z_stop_               = 0.;

  int index_tp_delta_idm_drmd_ = -1;  // #309 transfer-source slots (composite-owned)
  int index_tp_theta_idm_drmd_ = -1;
  int index_tp_delta_idr_drmd_ = -1;
  int index_tp_theta_idr_drmd_ = -1;

  int index_bg_G_over_aH_drmd_ = -1;    // bg-table slot (this species owns it)
  double Gamma0_drmd_ic_       = 0.;    // interaction rate today, computed at IC time
  double f_idr_drmd_           = 0.;    // idr/rho_tot fraction at IC (verbose diagnostic)
  double z_dec_drmd_           = -1.0;  // decoupling redshift (G_over_aH closest to 1); -1 = none
  double G_over_aH_tmp_        = 1e20;  // running |G_over_aH - 1| scan tracker
};
