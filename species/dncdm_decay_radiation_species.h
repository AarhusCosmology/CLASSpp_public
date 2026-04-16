#pragma once
#include <memory>

#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Dark Radiation from a single decaying NCDM flavor.
 */
class DNCDM_DecayRadiationSpecies : public BaseSpecies {
 public:
  DNCDM_DecayRadiationSpecies(int ncdm_id, const background* pba, const BackgroundModule* bgm)
      : BaseSpecies("DNCDM_DecayRadiation_" + std::to_string(ncdm_id), EnergyType::Radiation),
        ncdm_id_(ncdm_id), pba_(pba), bgm_(bgm) {}

  bool IsFreestreaming() const override {
    return true;
  }

  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

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

  int bi_rho_index() const {
    return index_bi_rho_;
  }
  int bg_rho_index() const {
    return index_bg_rho_;
  }
  int pt_F0_index() const {
    return index_pt_F0_;
  }

 private:
  int ncdm_id_;
  const background* pba_;
  const BackgroundModule* bgm_;

  int index_bi_rho_ = -1;
  int index_pt_F0_  = -1;
};