#pragma once
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

/** Baryons: rho ~ a^{-3}. Two perturbation variables: delta_b, theta_b. */
class BaryonsSpecies : public BaseSpecies {
public:
  explicit BaryonsSpecies(const background& pba)
    : BaseSpecies("Baryons", EnergyType::Matter), pba_(pba) {}

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/,
                          double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_b * pba_.H0 * pba_.H0
                              / (a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override { return pvecback[index_bg_rho_]; }
  double P(const double* /*pvecback*/) const override { return 0.; }
  double DpDloga(const double* /*pvecback*/) const override { return 0.; }

  void RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  double Delta(const perturb_vector* pv, const double* y, const double* /*pvecback*/) const override {
    return y[pv->index_pt_delta_b];
  }
  double Theta(const perturb_vector* pv, const double* y, const double* /*pvecback*/) const override {
    return y[pv->index_pt_theta_b];
  }
  double DeltaP(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const override { return 0.; }
  double RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const override { return 0.; }

private:
  const background& pba_;
};
