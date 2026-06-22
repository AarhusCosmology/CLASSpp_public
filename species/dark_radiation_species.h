#pragma once
#include <memory>
#include <string>

#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Dark Radiation: a single, self-contained free-streaming decay-radiation channel.
 *
 * One instance == one decay channel. A composite that emits N decay-radiation
 * products simply owns N instances (DR is not shared between channels). The decay
 * *source* is injected by the parent composite (in BackgroundDerivs / AddCouplingDerivs);
 * this species only carries dilution + the free-streaming Boltzmann hierarchy.
 */
class DarkRadiationSpecies : public BaseSpecies {
 public:
  // ── PerturbLayout ──────────────────────────────────────────────────────────
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_F0 = -1;  ///< base offset for the DR multipole hierarchy in pv->y
    int l_max  = -1;  ///< max multipole (ppr->l_max_dr)
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  DarkRadiationSpecies(const std::string& name, const background* pba, const BackgroundModule* bgm)
      : BaseSpecies(name, EnergyType::Radiation), pba_(pba), bgm_(bgm) {}

  /** Decay product: no direct Omega0 input, starts at zero. */
  double GetOmega0() const override {
    return 0.0;
  }

  bool IsFreestreaming() const override {
    return true;
  }

  // ── Background ──────────────────────────────────────────────────────────────
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
  double PPrime(double a,
                double H,
                const double* /*pvecback_B*/,
                const double* pvecback) const override {
    return a * H * (-4. / 3. * pvecback[index_bg_rho_]);
  }

  // ── Perturbations ────────────────────────────────────────────────────────────
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

  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;

  /** Re-seed the DR multipole hierarchy from the gauge-shifted IC. decay_corr is
      supplied by the owning composite (= aΓ·ρ_parent/ρ_dr; 0 when ρ_dr == 0). */
  void PerturbNewtonianReseed(const PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx,
                              double decay_corr) const;

  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

  int bg_rho_index() const {
    return index_bg_rho_;
  }
  int bi_rho_index() const {
    return index_bi_rho_;
  }
  double DarkRadiationRhoToday(const double* pvecback_integration) const override {
    return pvecback_integration[index_bi_rho_];
  }
  int pt_F0_index() const {
    return index_pt_F0_;
  }

  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;
  int transfer_delta_index() const {
    return index_tp_delta_;
  }
  int transfer_theta_index() const {
    return index_tp_theta_;
  }

 private:
  const background* pba_;
  const BackgroundModule* bgm_;

  int index_bi_rho_   = -1;
  int index_pt_F0_    = -1;
  int index_tp_delta_ = -1;  // #309 transfer-source slot
  int index_tp_theta_ = -1;
};
