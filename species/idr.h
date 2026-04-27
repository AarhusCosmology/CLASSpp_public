#pragma once
#include "background.h"
#include "base_species.h"
#include "perturbations.h"

class ThermodynamicsModule;  // forward declaration

/**
 * IDR: Interacting Dark Radiation.
 * Coupled to IDM_DR.
 */
class IDRSpecies : public BaseSpecies {
 public:
  explicit IDRSpecies(const background& pba)
      : BaseSpecies("IDR", EnergyType::Radiation), pba_(pba) {}

  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override {
    thm_ = thm;
  }
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
  }

  /**
   * IDR shear under TCA with IDM_DR. Returns the TCA shear prediction if all
   * TCA guards are satisfied; otherwise returns 0. Called by RhoPlusPShear
   * (when shear_idr is not in the y-vector) and by external code that needs
   * the TCA shear (perturb_vector_init, IDM_DR_IDR_Species::PrintVariables).
   */
  double TcaShearIdr(const perturb_vector* pv, const double* y, const perturb_workspace* ppw) const;

  double GetOmega0() const override {
    return pba_.Omega0_idr;
  }

  bool IsFreestreaming() const override {
    return true;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_idr * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_rho_] / 3.;
  }
  double DpDloga(const double* pvecback) const override {
    return -4. / 3. * pvecback[index_bg_rho_];
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
    return (pv->index_pt_delta_idr >= 0) ? y[pv->index_pt_delta_idr] : 0.;
  }
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return (pv->index_pt_theta_idr >= 0) ? y[pv->index_pt_theta_idr] : 0.;
  }
  double DeltaP(const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* /*ppw*/) const override {
    return (pv->index_pt_delta_idr >= 0) ? pvecback[index_bg_rho_] * y[pv->index_pt_delta_idr] / 3.
                                         : 0.;
  }
  double RhoPlusPShear(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

 private:
  const background& pba_;
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_             = nullptr;
};
