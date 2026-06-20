#pragma once
#include "background.h"
#include "base_species.h"
#include "perturbations.h"

/**
 * IDM_DR: Dark matter interacting with dark radiation.
 * Similar to CDM but with interaction terms in perturbations.
 */
class IDM_DRSpecies : public BaseSpecies {
 public:
  IDM_DRSpecies(const background& pba,
                double omega0_idm_dr,
                double a_idm_dr      = 0.,
                double nindex_idm_dr = 4.,
                double m_idm         = 1.e11)
      : BaseSpecies("IDM_DR", EnergyType::Matter), pba_(pba), Omega0_idm_dr_(omega0_idm_dr),
        a_idm_dr_(a_idm_dr), nindex_idm_dr_(nindex_idm_dr), m_idm_(m_idm) {}

  double a_idm_dr() const {
    return a_idm_dr_;
  }
  double nindex_idm_dr() const {
    return nindex_idm_dr_;
  }
  double m_idm() const {
    return m_idm_;
  }

  double GetOmega0() const override {
    return Omega0_idm_dr_;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = Omega0_idm_dr_ * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* /*pvecback*/) const override {
    return 0.;
  }

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;  ///< index_pt_delta_idm_dr
    int idx_theta = -1;  ///< index_pt_theta_idm_dr
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

  double DeltaP(const BaseSpecies::PerturbLayout& /*layout*/,
                const perturb_vector* /*pv*/,
                const double* /*y*/,
                const double* /*pvecback*/,
                const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }
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
  double Omega0_idm_dr_;
  double a_idm_dr_      = 0.;
  double nindex_idm_dr_ = 4.;
  double m_idm_         = 1.e11;
};
