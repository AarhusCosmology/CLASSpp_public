#pragma once
#include "background.h"
#include "base_species.h"
#include "perturbations.h"

/**
 * IDM_DRMD: Interacting Dark Matter in the DRMD model.
 */
class IDM_DRMDSpecies : public BaseSpecies {
 public:
  IDM_DRMDSpecies(const background& pba, double omega0_idm_drmd)
      : BaseSpecies("IDM_DRMD", EnergyType::Matter), pba_(pba), Omega0_idm_drmd_(omega0_idm_drmd) {}

  double GetOmega0() const override {
    return Omega0_idm_drmd_;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = Omega0_idm_drmd_ * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* /*pvecback*/) const override {
    return 0.;
  }
  double DpDloga(const double* /*pvecback*/) const override {
    return 0.;
  }

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;  ///< index_pt_delta_idm_drmd
    int idx_theta = -1;  ///< index_pt_theta_idm_drmd
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  // ── Perturbation index registration ────────────────────────────────────────

  /** Layout-based override (primary path). */
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  // ── PerturbDerivs ──────────────────────────────────────────────────────────

  /** Layout-based PerturbDerivs (primary path). */
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  /** Legacy PerturbDerivs. */
  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  // ── Stress-energy observables ───────────────────────────────────────────────

  /** Layout-based Delta (primary). */
  double Delta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  /** Legacy Delta: unreachable — composite IDM_DRMD_IDR_DRMD always dispatches
      the layout-aware Delta() on its IDM_DRMD child. */
  double Delta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  /** Layout-based Theta (primary). */
  double Theta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  /** Legacy Theta: unreachable — see Delta. */
  double Theta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double DeltaP(const perturb_vector* /*pv*/,
                const double* /*y*/,
                const double* /*pvecback*/,
                const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }
  double RhoPlusPShear(const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  // ── Switch-copy hook ────────────────────────────────────────────────────────

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

 private:
  const background& pba_;
  double Omega0_idm_drmd_;
};
