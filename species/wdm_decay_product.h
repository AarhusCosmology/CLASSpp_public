#pragma once
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "background.h"

class BackgroundModule;
class DCDMSpecies;

/**
 * Warm decay products of cold dark matter (arXiv:2606.14849).
 *
 * A DCDM parent decays χ → χ₁ + χ₂ into two identical daughters, each with
 * mass m_d = ε·m_χ/2 and kick velocity v = √(1−ε²)·c. All daughters are born
 * with the same physical momentum, so the comoving distribution f(q,τ) builds
 * up at the moving cutoff q_cut(τ) = a·kQKick. f and g ≡ ∂f/∂ln q are
 * integrated per momentum bin in the background ODE (injection via a
 * normalized Gaussian kernel in ln q — the energy-injection sum rule
 * Σᵢ dqᵢ qᵢ² εᵢ Jᵢ · factor/a⁴ = aΓρ_dcdm holds exactly by construction).
 *
 * Internal conventions: T_wdm = T_cmb, deg = 1, dimensionless kick momentum
 * kQKick = p/T₀ = 10 (only v_kick is physical; the parent mass drops out).
 * Dimensionless daughter mass M = kQKick·ε/v.
 *
 * Perturbations (Task 4) evolve UNNORMALIZED multipoles ψ_ℓ = (Δf)_ℓ, which
 * stay regular through f = 0. Synchronous gauge only; no tensor slots; no
 * fluid approximation.
 */
class WdmDecayProductSpecies : public NCDMBaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "dcdm_wdm";
  static constexpr double kQKick              = 10.0;  // p/T0, internal convention

  WdmDecayProductSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm);

  double Gamma() const {
    return Gamma_;
  }
  double vkick() const {
    return vkick_;
  }
  double epsilon_retention() const {
    return eps_;
  }
  const std::vector<double>& dq() const {
    return dq_;
  }
  const std::vector<double>& u() const {
    return u_;
  }

  /** Wired by DCDM_WDM_Species; used by RhoDotOverRho / FillSources. */
  void SetParent(const DCDMSpecies* dcdm) {
    parent_ = dcdm;
  }

  const std::optional<double>& Omega_ini_pending() const {
    return Omega_ini_pending_;
  }
  const std::optional<double>& Omega_combined_pending() const {
    return Omega_combined_pending_;
  }

  /**
   * Injection source Jᵢ = dfᵢ/dτ (and optionally dJᵢ/dln q) for the current
   * (a, ρ_dcdm). Normalized Gaussian kernel in u = ln q at u_cut = ln(a·kQKick):
   *   Jᵢ = A·Gᵢ/D,  Gᵢ = exp(−(uᵢ−u_cut)²/2σ²),  σ = kernel_width·Δu,
   *   D  = Σⱼ dqⱼ qⱼ² εⱼ(a) Gⱼ,  A = aΓρ_dcdm·a⁴/factor_.
   * A smooth Gaussian gate ramps the injection on over the first few kernel
   * widths of the grid (pre-grid decays, ~1e-6 of the total for defaults, are
   * suppressed/dropped): with the gate at 1 the energy sum rule is exact.
   * The kernel argument is clamped (G = 0 beyond ~11σ): under -ffast-math,
   * exp of extreme arguments and denormal intermediates are not IEEE-safe.
   * dJdlnq may be nullptr. Background-thread only (uses mutable scratch).
   */
  void FillInjection(double a, double rho_dcdm, double* J, double* dJdlnq) const;

  // ── Background (Task 2) ────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;
  double RhoDotOverRho(const double* pvecback, double a_prime_over_a) const override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Decay-product neutrality overrides ─────────────────────────────────────
  // The daughter starts empty: it contributes nothing to the early-time
  // radiation budget, fnu, N_eff, or the flatness budget (the combined sector
  // density is owned by the DCDM_WDM composite).
  double GetOmega0() const override {
    return 0.;
  }
  double NeutrinoOmega0() const override {
    return 0.;
  }
  double NeffContribution(double /*z*/) const override {
    return 0.;
  }
  bool IsFreestreaming() const override {
    return false;
  }
  double BackgroundAIni(double a_proposed, double /*tol*/) const override {
    return a_proposed;  // empty at early times; no relativistic-start constraint
  }
  void CheckUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const override {}
  bool IsUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const override {
    return true;
  }
  void PrintNeffInfo() const override;
  void PrintMassInfo() const override;

  // ── Perturbations (Task 4) ─────────────────────────────────────────────────
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
                     const perturb_parameters_and_workspace& ppaw) const override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  /** Tensor modes: no slots registered — the daughter's tensor anisotropic
   *  stress is neglected (documented scope limit; it is empty until late times). */
  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& /*layout*/,
                                         perturb_vector* /*pv*/,
                                         const precision* /*ppr*/,
                                         int& /*index_pt*/,
                                         const perturb_workspace* /*ppw*/,
                                         int /*gauge*/) override {}
  /** Synchronous gauge only (guarded at CreateAll); base NCDM transform would
   *  zero the injected ICs via GetDlnf0Dlnq() == 0. */
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& /*layout*/,
                                     double* /*y*/,
                                     const PerturbIcContext& /*ctx*/) override {}

  // ── Index accessors for the composite ──────────────────────────────────────
  int bg_f_index() const {
    return index_bg_f_;
  }
  int bg_dfdlnq_index() const {
    return index_bg_dfdlnq_;
  }
  int bg_inj_index() const {
    return index_bg_inj_;
  }
  int bg_number_index() const {
    return index_bg_number_;
  }
  int bi_f_index() const {
    return index_bi_f_;
  }
  int bi_dfdlnq_index() const {
    return index_bi_dfdlnq_;
  }

 protected:
  /** Only consumed by tensor derivs / base gauge transform, both disabled. */
  double GetDlnf0Dlnq(int /*iq*/, const double* /*pvecback*/) const override {
    return 0.;
  }

 private:
  /** ndf15 error weights are max(|y|, threshold=1e-15); components starting at
   *  EXACTLY zero force absolute step errors <= rtol*1e-15 while f climbs ~15
   *  decades at each bin's injection onset — micro-steps at ~n_bins onsets make
   *  the background solve crawl. So f bins start at this seed instead
   *  (subtracted again in ComputeBackground; its would-be density is ~1e-9 of
   *  rho_gamma and never enters any physics). */
  static constexpr double kFSeed = 1e-10;

  const background* pba_;
  const DCDMSpecies* parent_ = nullptr;

  double Gamma_ = 0.;  // decay rate, CLASS units (Mpc^-1 after both conversions)
  double vkick_ = 0.;  // kick velocity in units of c
  double eps_   = 0.;  // total mass retention (m_d = eps*m_parent/2)

  int n_bins_          = 96;
  double q_min_ratio_  = 1e-4;
  double kernel_width_ = 1.0;  // sigma in units of the log-grid spacing du_
  int l_max_input_     = -1;   // -1 → ppr->l_max_ncdm

  double du_ = 0.;          // uniform spacing in u = ln q
  std::vector<double> u_;   // ln q_i
  std::vector<double> dq_;  // q_i * du_ (midpoint rule in ln q)

  std::optional<double> Omega_ini_pending_;
  std::optional<double> Omega_combined_pending_;

  // Background indices (per-bin arrays are n_bins_ long)
  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;
  int index_bg_f_        = -1;
  int index_bg_dfdlnq_   = -1;
  int index_bg_inj_      = -1;
  int index_bi_f_        = -1;
  int index_bi_dfdlnq_   = -1;

  mutable std::vector<double> scratch_G_;  // kernel weights (background thread only)
};
