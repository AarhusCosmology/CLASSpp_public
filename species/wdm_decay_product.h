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
 * integrated per momentum bin in the background ODE. Injection is a
 * conservative cell-integrated Gaussian (erf) energy deposit in ln q: the
 * energy-injection sum rule Σᵢ dqᵢ qᵢ² εᵢ Jᵢ · factor/a⁴ = aΓρ_dcdm holds
 * exactly for every cutoff position, and the source is smooth in time (both
 * evolvers integrate it).
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

  /** Static-sparsity seed for the injection source. ndf15's numjac derives the
   *  Jacobian sparsity pattern from exact zeros and locks it permanently once it
   *  repeats (tools/evolver_ndf15.cpp: pattern from fabs(dFdy)!=0, locked after
   *  trust_sparse repeats, never re-derived). A moving compact source would lock
   *  a too-narrow pattern early and corrupt the grouped Jacobian when the
   *  footprint sweeps on. Flooring every weight (renormalized, Σw stays exactly
   *  1) keeps every injection row structurally coupled to rho_dcdm at all times
   *  at a physically invisible level (~1e-12 of the instantaneous deposit).
   *  This trick only works because J is LINEAR in rho_dcdm (the floor is a
   *  fixed fraction of the same A that scales with rho_dcdm) — a nonlinear
   *  coupling would decouple the seeded pattern from the finite-difference
   *  Jacobian numjac actually probes, defeating the pattern-lock. */
  static constexpr double kSparsityFloor = 1e-12;

  // ── Fiducial-cosmology placement helpers ───────────────────────────────────
  // Placement-only fixed LCDM background — NOT physics. Maps a decay time t (Mpc)
  // to the scale factor a at which those decays inject, so the momentum grid is a
  // function of Gamma alone (invariant across the run's cosmology). See spec
  // docs/superpowers/specs/2026-07-18-dcdm-wdm-adaptive-momentum-grid-design.md.
  static constexpr double kFidOmegaM = 0.31;
  static constexpr double kFidOmegaR = 9.2e-5;            // photons+nu, T_cmb=2.7255, Neff=3.044
  static constexpr double kFidH      = 0.674;             // h
  static double FiducialCosmicTime(double a);             // t(a) [Mpc], radiation+matter
  static double FiducialScaleFactorAtTime(double t_mpc);  // inverse (monotone)

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
  const std::vector<double>& u_edge() const {
    return u_edge_;
  }
  /** Injection kernel width at u_cut: the analytic local quantile-bin width
   *  Δu_loc = [(F_hi−F_lo)/N] · H_fid(a) e^{+Γ t_fid(a)} / Γ  (a = e^{u_cut}/kQKick),
   *  clamped to [kSigmaMin, kSigmaMax]. Smooth in u_cut; evaluated in log form so
   *  no overflow/inf is ever materialized (-ffast-math discipline). Floor = the
   *  background-table resolvability bound; cap = the placement-bias budget. */
  double SigmaAt(double u_cut) const;

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
   * (a, ρ_dcdm). Conservative cell-integrated Gaussian (erf) energy deposit in
   * u = ln q at u_cut = ln(a·kQKick): the kick's energy A = aΓρ_dcdm·a⁴/factor_
   * is split by GaussWeights over the cells (σ = SigmaAt(u_cut)), off-grid tails
   * clamped into the edge bins, so Σ Wᵢ = 1 identically and the exact sum rule
   * Σ dqᵢ qᵢ² εᵢ Jᵢ = A holds for every u_cut. dJ/dln q is the full derivative of
   * J(u) = A·φ_σ(u−u_cut)/(q³ε): the cell-averaged raw-kernel derivative (normal-pdf
   * differences at the cell edges, smooth, no edge special cases) MINUS the
   * measure's own log-derivative (3 + q²/ε²)·J, which supplies the g-channel's
   * ρ-weighted moment (Σ dq q² ε dJdlnq ≈ −3A non-relativistically, by integration
   * by parts) that drives the daughter's fluid-limit metric-driven growth.
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

  /** Build the injection-adapted quantile grid (equal injected decays per bin);
   *  fills q_/q_bg_/u_/dq_/u_edge_. Depends on Gamma + the fiducial H only. */
  void BuildInjectionAdaptedGrid();
  /** Cell-integrated Gaussian (erf) energy-deposit weights for a kick at u_cut:
   *  w[i] = Φ((u_edge[i+1]−u_cut)/σ) − Φ((u_edge[i]−u_cut)/σ), plus the below-/
   *  above-grid tail mass clamped into the corresponding edge bin, so Σ w = 1
   *  identically (energy conserved for every u_cut, no D-renormalization) and
   *  the source is C-infinity in time. */
  void GaussWeights(double u_cut, double sigma, double* w) const;

  const background* pba_;
  const DCDMSpecies* parent_ = nullptr;

  double Gamma_ = 0.;  // decay rate, CLASS units (Mpc^-1 after both conversions)
  double vkick_ = 0.;  // kick velocity in units of c
  double eps_   = 0.;  // total mass retention (m_d = eps*m_parent/2)

  int n_bins_        = 96;
  double q_edge_tol_ = 1e-3;  // decayed-fraction trimmed off each end (grid F-span);
                              // the trimmed tail is clamped into the edge bin, not dropped
  int l_max_input_ = -1;      // -1 → ppr->l_max_ncdm
  std::optional<double> q_min_ratio_floor_;  // optional hard floor on q_lo/q_kick

  static constexpr double kSigmaMin = 0.03;  // ≈ 4 background-table samples (7e-3 in ln a)
  static constexpr double kSigmaMax = 0.25;  // bounds the wide-bin placement bias
  double F_lo_ = 0.;  // realized decayed-fraction span of the grid (0,0 = fallback grid)
  double F_hi_ = 0.;

  std::vector<double> u_;       // ln q_i (bin centres = q(F_mid[i]))
  std::vector<double> dq_;      // cell widths q_edge[i+1] - q_edge[i]
  std::vector<double> u_edge_;  // N+1 cell edges in u (F-quantile edges)

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

  mutable std::vector<double> scratch_w_;  // erf energy-deposit weights (background thread only)
};
