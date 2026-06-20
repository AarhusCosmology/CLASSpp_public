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
    int idx_Gamma = -1;  // PPF dynamical variable (only populated when use_ppf == _TRUE_)
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
               double Omega_EDE,
               bool use_ppf,
               double c_gamma_over_c_fld);

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
  bool use_ppf() const {
    return use_ppf_;
  }
  double c_gamma_over_c_fld() const {
    return c_gamma_over_c_fld_;
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
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

  /**
   * Compute PPF fluid contribution. Writes ppw->delta_rho_fld,
   * rho_plus_p_theta_fld, delta_p_fld, Gamma_prime_fld (module then
   * accumulates into rho/theta/p totals). Only called when use_ppf() == _TRUE_.
   * Moved from PerturbationsModule — PPF is fluid-specific physics that depends
   * on the rest of the universe (rho_plus_p_theta, rho_plus_p_shear, delta_rho
   * from non-fluid species must be fully populated before calling).
   */
  void ComputePpf(double k,
                  double a,
                  double a_prime_over_a,
                  const precision* ppr,
                  const double* y,
                  perturb_workspace* ppw) const;

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
                     const perturb_parameters_and_workspace& ppaw) override;
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) override;
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

  // Layout-based read virtuals: do the real work.
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

 private:
  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
  double Omega0_fld_;

  equation_of_state fluid_eos_ = CLP;
  double w0_fld_               = -1.;
  double wa_fld_               = 0.;
  double cs2_fld_              = 1.;
  double Omega_EDE_            = 0.;
  bool use_ppf_                = true;
  double c_gamma_over_c_fld_   = 0.4;

  int index_bg_rho_fld_        = -1;
  int index_bg_w_fld_          = -1;
  int index_bg_dw_over_da_fld_ = -1;
  int index_bi_rho_fld_        = -1;

  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
