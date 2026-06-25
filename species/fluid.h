#pragma once

#include "background.h"
#include "base_species.h"
#include "species_build_context.h"

class BackgroundModule;

/**
 * Dark energy fluid species with evolving equation of state w(a).
 *
 * Background:
 *   rho_fld is ODE-integrated via BackgroundDerivs called from background_derivs_member().
 *   ComputeBackground reads rho_fld from pvecback_B.
 *
 * Perturbations (true fluid):
 *   dy[delta_fld] = -(1+w)(theta_fld + metric_continuity)
 *   dy[theta_fld] = -(1-3*cs2)*H'*theta_fld + cs2*k^2/(1+w)*delta_fld + metric_euler
 */
class FluidSpecies : public BaseSpecies {
 public:
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;
    int idx_theta = -1;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  FluidSpecies(const background& pba,
               double omega0_fld,
               equation_of_state fluid_eos,
               double w0_fld,
               double wa_fld,
               double cs2_fld,
               double Omega_EDE);

  double GetOmega0() const override {
    return Omega0_fld_;
  }

  equation_of_state fluid_eos() const {
    return fluid_eos_;
  }
  double w0_fld() const {
    return w0_fld_;
  }
  double wa_fld() const {
    return wa_fld_;
  }
  double cs2_fld() const {
    return cs2_fld_;
  }
  double Omega_EDE() const {
    return Omega_EDE_;
  }
  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  /** Returns w_fld from pvecback. Used by perturbations to read the cached value. */
  double W(const double* pvecback) const {
    return pvecback[index_bg_w_fld_];
  }

  /**
   * Fluid equation-of-state evaluation. Owns w_fld(a), dw_fld/da(a), and the
   * analytic integral used by background initial conditions. Implementation
   * moved here from BackgroundModule::background_w_fld because this is pure
   * fluid physics.
   */
  int ComputeWFld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) const;

  // ── Perturbations ──────────────────────────────────────────────────────────

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;

  // Layout-based signatures: do the real work.
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
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

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

  int bi_rho_index() const {
    return index_bi_rho_fld_;
  }
  int bg_w_index() const {
    return index_bg_w_fld_;
  }
  int bg_dw_over_da_index() const {
    return index_bg_dw_over_da_fld_;
  }

 protected:
  const BackgroundModule* bgm_ = nullptr;
  int index_bg_rho_fld_        = -1;
  int index_bg_w_fld_          = -1;
  int index_bg_dw_over_da_fld_ = -1;
  double cs2_fld_              = 1.;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;

 private:
  const background& pba_;
  double Omega0_fld_;

  equation_of_state fluid_eos_ = CLP;
  double w0_fld_               = -1.;
  double wa_fld_               = 0.;
  double Omega_EDE_            = 0.;

  int index_bi_rho_fld_ = -1;
};
