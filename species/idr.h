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

  // ── PerturbLayout ──────────────────────────────────────────────────────────

  /**
   * Per-species perturbation slot layout for IDR.
   *
   * When idr_nature == idr_free_streaming AND TCA is off:
   *   idx_delta, idx_theta, idx_shear, idx_l3 are set; l_max == pba.l_max_idr.
   * When idr_nature == idr_fluid OR TCA is on (shear not registered):
   *   idx_delta, idx_theta set; idx_shear == idx_l3 == -1; l_max == -1.
   * When RSA is on (nothing registered):
   *   all fields are -1.
   */
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;  ///< index_pt_delta_idr
    int idx_theta = -1;  ///< index_pt_theta_idr
    int idx_shear = -1;  ///< index_pt_shear_idr (free-streaming + tca_off only)
    int idx_l3    = -1;  ///< index_pt_l3_idr    (free-streaming + tca_off + l_max>=3)
    int l_max     = -1;  ///< pba.l_max_idr when full hierarchy, else -1
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  /**
   * IDR shear under TCA with IDM_DR. Returns the TCA shear prediction if all
   * TCA guards are satisfied; otherwise returns 0. Called by the layout-aware
   * RhoPlusPShear (when shear_idr is not in the y-vector) and by external code
   * that needs the TCA shear (perturb_vector_init, IDM_DR_IDR_Species::PrintVariables).
   *
   * `layout` must be the IDR child's PerturbLayout (nested inside the
   * IDM_DR_IDR composite layout). We read layout.idx_shear / layout.idx_theta
   * to drive the "shear not in y" detection and the theta_idr value.
   */
  double TcaShearIdr(const PerturbLayout& layout,
                     const double* y,
                     const perturb_workspace* ppw) const;

  // ── Perturbation index registration ────────────────────────────────────────

  /** Layout-based override (primary path, called by module when IDR is migrated). */
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  /** Legacy override: stubbed — composite IDM_DR_IDR always dispatches the
      layout-aware overload on its IDR child. */
  void RegisterPerturbationIndices(perturb_vector* /*pv*/,
                                   const precision* /*ppr*/,
                                   int& /*index_pt*/,
                                   const perturb_workspace* /*ppw*/,
                                   int /*gauge*/) override {}

  // ── PerturbDerivs ──────────────────────────────────────────────────────────

  /** Layout-based PerturbDerivs (primary path). */
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  /** Legacy PerturbDerivs: delegates to layout-based path via pv indices. */
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
  /** Legacy Delta: unreachable — IDR is only accessed via IDM_DR_IDR composite's
      layout-aware Delta(). Stub kept for vtable completeness. */
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

  /** Layout-based DeltaP (primary). */
  double DeltaP(const BaseSpecies::PerturbLayout& layout,
                const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;
  /** Legacy DeltaP: unreachable — see Delta. */
  double DeltaP(const perturb_vector* /*pv*/,
                const double* /*y*/,
                const double* /*pvecback*/,
                const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  /** Layout-based RhoPlusPShear (primary). */
  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;
  /** Legacy RhoPlusPShear: unreachable — see Delta. */
  double RhoPlusPShear(const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

  // ── Switch-copy hook (for Phase C / Task 25) ───────────────────────────────

  /**
   * Copy IDR perturbations from old layout to new layout across an approximation switch.
   * Handles RSA (rsa_idr off->on) and TCA (tca_idm_dr on->off) transitions.
   * The TCA on->off switch seeds shear_idr from TcaShearIdr; that special case
   * is handled by the module (Task 23/25) since it needs IDM_DR_IDR_Species context.
   * This method covers the simpler slot-by-slot copies (RSA and UFA transitions).
   */
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;

 private:
  const background& pba_;
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_             = nullptr;
};
