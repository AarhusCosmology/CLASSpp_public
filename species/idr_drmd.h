#pragma once
#include "background.h"
#include "base_species.h"
#include "perturbations.h"

/**
 * IDR_DRMD: Dark Radiation in the DRMD model.
 */
class IDR_DRMDSpecies : public BaseSpecies {
 public:
  IDR_DRMDSpecies(const background& pba, double omega0_idr_drmd)
      : BaseSpecies("IDR_DRMD", EnergyType::Radiation), pba_(pba),
        Omega0_idr_drmd_(omega0_idr_drmd) {}

  double GetOmega0() const override {
    return Omega0_idr_drmd_;
  }

  /** IDR_DRMD is dark radiation: contributes its Omega0 to the early-time radiation density. */
  double GetRadiationOmega0() const override {
    return Omega0_idr_drmd_;
  }

  bool IsFreestreaming() const override {
    return true;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = Omega0_idr_drmd_ * pba_.H0 * pba_.H0 /
                              (a_rel * a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_rho_] / 3.;
  }
  double PPrime(double a,
                double H,
                const double* /*pvecback_B*/,
                const double* pvecback) const override {
    return a * H * (-4. / 3. * pvecback[index_bg_rho_]);
  }

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;  ///< index_pt_delta_idr_drmd
    int idx_theta = -1;  ///< index_pt_theta_idr_drmd
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  // ── Perturbation index registration ────────────────────────────────────────

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  // ── PerturbDerivs ──────────────────────────────────────────────────────────

  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  // ── Stress-energy observables ───────────────────────────────────────────────

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

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& /*layout*/,
                       const perturb_vector* /*pv*/,
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
  double Omega0_idr_drmd_;
};
