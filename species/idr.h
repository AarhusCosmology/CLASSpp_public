#pragma once
#include <vector>

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
  IDRSpecies(const background& pba,
             double omega0_idr,
             bool has_sibling_idm_dr          = false,
             double T_idr                     = 0.,
             int l_max_idr                    = 0,
             double b_idr                     = 0.,
             int idr_nature                   = 0,
             std::vector<double> alpha_idm_dr = {},
             std::vector<double> beta_idr     = {})
      : BaseSpecies("IDR", EnergyType::Radiation), pba_(pba), Omega0_idr_(omega0_idr),
        has_sibling_idm_dr_(has_sibling_idm_dr), T_idr_(T_idr), l_max_idr_(l_max_idr),
        b_idr_(b_idr), idr_nature_(idr_nature), alpha_idm_dr_(std::move(alpha_idm_dr)),
        beta_idr_(std::move(beta_idr)) {}

  bool has_sibling_idm_dr() const {
    return has_sibling_idm_dr_;
  }

  double T_idr() const {
    return T_idr_;
  }
  int l_max_idr() const {
    return l_max_idr_;
  }

  double b_idr() const {
    return b_idr_;
  }
  int idr_nature() const {
    return idr_nature_;
  }
  const std::vector<double>& alpha_idm_dr() const {
    return alpha_idm_dr_;
  }
  const std::vector<double>& beta_idr() const {
    return beta_idr_;
  }

  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override {
    thm_ = thm;
  }
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
  }

  double GetOmega0() const override {
    return Omega0_idr_;
  }

  /** IDR is dark radiation: contributes its Omega0 to the early-time radiation density. */
  double GetRadiationOmega0() const override {
    return Omega0_idr_;
  }

  bool IsFreestreaming() const override {
    return true;
  }

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = Omega0_idr_ * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel * a_rel);
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

  /**
   * Per-species perturbation slot layout for IDR.
   *
   * When idr_nature == idr_free_streaming AND TCA is off:
   *   idx_delta, idx_theta, idx_shear, idx_l3 are set; l_max == IDRSpecies::l_max_idr().
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
    int l_max     = -1;  ///< IDRSpecies::l_max_idr() when full hierarchy, else -1
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

  double RhoPlusPShear(const BaseSpecies::PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

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
  double Omega0_idr_;
  bool has_sibling_idm_dr_;
  double T_idr_                    = 0.;
  int l_max_idr_                   = 0;
  const ThermodynamicsModule* thm_ = nullptr;
  const perturbs* ppt_             = nullptr;
  double b_idr_                    = 0.;
  int idr_nature_                  = 0;
  std::vector<double> alpha_idm_dr_;
  std::vector<double> beta_idr_;
};
