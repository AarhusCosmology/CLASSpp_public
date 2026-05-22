#pragma once
#include <memory>
#include <string>

#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Dark Radiation from a single decaying NCDM flavor.
 */
class DNCDM_DecayRadiationSpecies : public BaseSpecies {
 public:
  // ── PerturbLayout ──────────────────────────────────────────────────────────
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_F0 = -1;  ///< base offset for the DR multipole hierarchy in pv->y
    int l_max  = -1;  ///< max multipole (pv->l_max_dr)
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  DNCDM_DecayRadiationSpecies(const std::string& parent_name,
                              const background* pba,
                              const BackgroundModule* bgm)
      : BaseSpecies("DNCDM_DecayRadiation_" + parent_name, EnergyType::Radiation), pba_(pba),
        bgm_(bgm) {}

  /** Decay product: no direct Omega0 input, starts at zero. */
  double GetOmega0() const override {
    return 0.0;
  }

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

  // ── Layout-based perturbation interface ────────────────────────────────────
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  /** Legacy scalar register: no-op — dual-written by layout-based path above. */
  void RegisterPerturbationIndices(perturb_vector* /*pv*/,
                                   const precision* /*ppr*/,
                                   int& /*index_pt*/,
                                   const perturb_workspace* /*ppw*/,
                                   int /*gauge*/) override {}

  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  /** Legacy PerturbDerivs: no-op — composite routes through layout-based path. */
  void PerturbDerivs(double /*tau*/,
                     const double* /*y*/,
                     double* /*dy*/,
                     const perturb_parameters_and_workspace& /*ppaw*/) override {}

  double Delta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  /** Legacy Delta: no-op — composite routes through layout-based path. */
  double Delta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double Theta(const BaseSpecies::PerturbLayout& layout,
               const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;

  /** Legacy Theta: no-op — composite routes through layout-based path. */
  double Theta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;

  /** Legacy DeltaP: no-op — composite routes through layout-based path. */
  double DeltaP(const perturb_vector* /*pv*/,
                const double* /*y*/,
                const double* /*pvecback*/,
                const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  /** Legacy RhoPlusPShear: no-op — composite routes through layout-based path. */
  double RhoPlusPShear(const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  int bi_rho_index() const {
    return index_bi_rho_;
  }
  int bg_rho_index() const {
    return index_bg_rho_;
  }

 private:
  const background* pba_;
  const BackgroundModule* bgm_;

  int index_bi_rho_ = -1;
};