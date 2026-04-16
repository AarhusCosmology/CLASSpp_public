#pragma once

#include "base_species.h"

struct background;
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
  explicit FluidSpecies(const background& pba);

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;
  bool RequiresDeferredPerturbDerivs() const override {
    return true;
  }
  bool RequiresDeferredBackground() const override {
    return true;
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  /** Called by BackgroundModule::background_functions() before ComputeBackground().
   *  Writes w_fld and dw/da (already computed by BackgroundModule) into pvecback
   *  using FluidSpecies's private indices. */
  void WriteWFld(double w_fld, double dw_over_da_fld, double* pvecback) const;

  /** Returns w_fld from pvecback. Used by perturbations to read the cached value. */
  double W(const double* pvecback) const {
    return pvecback[index_bg_w_fld_];
  }

  // ── Perturbations ──────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;
  double Delta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double DeltaP(const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

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

  int index_bg_rho_fld_        = -1;
  int index_bg_w_fld_          = -1;
  int index_bg_dw_over_da_fld_ = -1;
  int index_bi_rho_fld_        = -1;
};
