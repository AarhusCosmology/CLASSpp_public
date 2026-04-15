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
  explicit IDM_DRSpecies(const background& pba)
      : BaseSpecies("IDM_DR", EnergyType::Matter), pba_(pba) {}

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_idm_dr * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel);
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

  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  double Delta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return y[pv->index_pt_delta_idm_dr];
  }
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return y[pv->index_pt_theta_idm_dr];
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

 private:
  const background& pba_;
};
