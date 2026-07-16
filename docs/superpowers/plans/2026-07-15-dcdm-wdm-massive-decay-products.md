# DCDM → Massive Decay Products (dcdm_wdm) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement dark matter decaying to two identical massive daughters (arXiv:2606.14849) as a `dcdm_wdm` composite species, validated by unit tests + a physics notebook.

**Architecture:** Composite `DCDM_WDM_Species` (mirrors `DNCDM_DR_Species`) couples the existing `DCDMSpecies` parent to a new `WdmDecayProductSpecies : NCDMBaseSpecies` daughter. The daughter's phase-space distribution f(qᵢ,τ) and its log-slope g(qᵢ,τ)=∂f/∂ln q are integrated per momentum bin in the background ODE, sourced by a normalized Gaussian kernel at the moving injection cutoff q_cut(τ)=a·q_kick (energy conservation exact by construction). Perturbations evolve **unnormalized** multipoles ψ_ℓ(qᵢ)≡(Δf)_ℓ (regular at f=0): the standard CLASS ncdm hierarchy with dln f₀/dln q replaced by g, plus an injection source Jᵢ·δ_dcdm on ℓ=0.

**Tech Stack:** C++17, CMake (+Makefile shim), classy Python wrapper (auto-generated cclassy.pxd — never hand-edit), jupyter/matplotlib.

**Spec:** `docs/superpowers/specs/2026-07-15-dcdm-wdm-massive-decay-products-design.md`

## Global Constraints

- Branch: `ddm-massive-decay-products` (already created). NEVER `git add -A` (in-source CMake/Xcode artifacts get swept in); stage explicit paths only.
- Subagents: never `git checkout` / `git reset` (shared working tree).
- C++ only, no `extern "C"` / `#ifdef __cplusplus` guards. Match surrounding code style (2-space indent, `// ──` section banners, doc comments on class methods).
- New test executables go in **BOTH** `CMakeLists.txt` (add_executable + the foreach list) and `Makefile` `TEST_TARGETS` — CI runs `make test` and only builds targets named there.
- Errors from species constructors/factories: `throw std::invalid_argument(...)` (existing pattern), NOT class_test.
- Synchronous gauge only (guard at CreateAll); scalar modes only (tensor slots not registered); no fluid approximation for the daughter.
- Internal conventions: T_wdm = T_cmb (T_=1), deg=1, q_kick ≡ p/T₀ = 10.0, M = q_kick·ε/v_kick.
- Never require bit-identical outputs in verification; use relative tolerances.
- Build: `make class` (wraps `cmake --build build/cmake --target class --parallel`). Tests: `make test-dcdm-wdm` then run `./build/cmake/test-dcdm-wdm` from repo root (tests assume CWD = repo root).

---

### Task 1: Daughter species — construction, momentum grid, injection kernel

**Files:**
- Create: `species/wdm_decay_product.h`
- Create: `species/wdm_decay_product.cpp`
- Create: `species/dcdm_wdm_test.cpp`
- Modify: `CMakeLists.txt` (CLASS_SPECIES_FILES list ~line 60-101; test section lines 212-234)
- Modify: `Makefile` (TEST_TARGETS list, lines 4-18)

**Interfaces:**
- Produces class `WdmDecayProductSpecies : NCDMBaseSpecies` with:
  - ctor `WdmDecayProductSpecies(FileContent* pfc, const std::string& instance_name, const NcdmSettings& settings, const background* pba, const BackgroundModule* bgm)`
  - `static constexpr std::string_view kTypeName = "dcdm_wdm";` and `static constexpr double kQKick = 10.0;`
  - `double Gamma() const`, `double vkick() const`, `double epsilon_retention() const`, `const std::vector<double>& dq() const`, `const std::vector<double>& u() const`
  - `void FillInjection(double a, double rho_dcdm, double* J, double* dJdlnq) const` (dJdlnq may be nullptr)
  - `void SetParent(const DCDMSpecies* dcdm)`
  - `const std::optional<double>& Omega_ini_pending() const`, `const std::optional<double>& Omega_combined_pending() const`
  - index accessors (used by Task 2/3/4): `bg_f_index()`, `bg_dfdlnq_index()`, `bg_inj_index()`, `bg_number_index()`, `bi_f_index()`, `bi_dfdlnq_index()`
- Consumes: `NCDMBaseSpecies` DeferInit ctor, `SpeciesInput`, `SetDegAndFactor`.

- [ ] **Step 1: Write the header** `species/wdm_decay_product.h`:

```cpp
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
   * Zero while q_cut is below the grid (pre-grid decays dropped; negligible).
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
  double PPrime(double a, double H, const double* pvecback_B, const double* pvecback) const override;
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
  const background* pba_;
  const DCDMSpecies* parent_ = nullptr;

  double Gamma_ = 0.;  // decay rate, CLASS units (Mpc^-1 after both conversions)
  double vkick_ = 0.;  // kick velocity in units of c
  double eps_   = 0.;  // total mass retention (m_d = eps*m_parent/2)

  int n_bins_          = 96;
  double q_min_ratio_  = 1e-4;
  double kernel_width_ = 1.0;  // sigma in units of the log-grid spacing du_
  int l_max_input_     = -1;   // -1 → ppr->l_max_ncdm

  double du_ = 0.;             // uniform spacing in u = ln q
  std::vector<double> u_;      // ln q_i
  std::vector<double> dq_;     // q_i * du_ (midpoint rule in ln q)

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
```

- [ ] **Step 2: Write the constructor + FillInjection + prints** in `species/wdm_decay_product.cpp`. Task 2 appends the background methods; for now add placeholder empty bodies ONLY for the pure-virtual/override methods the header declares that are implemented in later tasks (`RegisterBackgroundIndices` … `WriteOutputColumns`) so the file links — each placeholder body is `{}` (or `return {};` / `return 0.;`) with a `// Task N` comment:

```cpp
#include "wdm_decay_product.h"

#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "dcdm.h"
#include "perturbations_module.h"
#include "species/species_input.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

WdmDecayProductSpecies::WdmDecayProductSpecies(FileContent* pfc,
                                               const std::string& instance_name,
                                               const NcdmSettings& settings,
                                               const background* pba,
                                               const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings, DeferInit{}),
      pba_(pba) {
  bgm_ = bgm;
  SpeciesInput input(pfc, instance_name);

  // NCDM-family keys with no meaning here: the distribution is built from the
  // decay, not from a thermal PSD. Reject loudly instead of silently ignoring.
  for (const char* key : {"m", "T", "deg", "ksi", "Omega", "omega", "max_q",
                          "quadrature_strategy", "momenta_bins_bg", "use_psd_file",
                          "psd_filename", "psd_parameters"}) {
    if (input.get<std::string>(key).has_value()) {
      throw std::invalid_argument("species '" + instance_name + "': parameter '" +
                                  std::string(key) + "' is not supported by dcdm_wdm");
    }
  }

  // Internal thermodynamic convention: T_wdm = T_cmb, deg = 1, no chemical potential.
  T_       = 1.0;
  ksi_     = 0.;
  m_in_eV_ = 0.;          // the daughter mass in eV is not physical (kQKick is a convention)
  SetDegAndFactor(1.0);   // recomputes factor_ with T_ = 1

  // Decay rate: exactly one of Gamma [km/s/Mpc], log10Gamma, lifetime [yr], log10lifetime.
  auto Gamma_value         = input.get<double>("Gamma");
  auto log10Gamma_value    = input.get<double>("log10Gamma");
  auto lifetime_value      = input.get<double>("lifetime");
  auto log10lifetime_value = input.get<double>("log10lifetime");
  int n_provided = (int) Gamma_value.has_value() + (int) log10Gamma_value.has_value() +
                   (int) lifetime_value.has_value() + (int) log10lifetime_value.has_value();
  if (n_provided != 1) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': specify exactly one of Gamma, log10Gamma, lifetime, log10lifetime");
  }
  double Gamma_raw = 0.;
  if (Gamma_value) {
    Gamma_raw = *Gamma_value;
  }
  else if (log10Gamma_value) {
    Gamma_raw = std::pow(10., *log10Gamma_value);
  }
  else if (lifetime_value) {
    Gamma_raw = 1. / *lifetime_value / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  else {  // log10lifetime
    Gamma_raw = 1. / std::pow(10., *log10lifetime_value) / (365 * 24 * 60 * 60) *
                _Mpc_over_m_ * 1e-3;
  }
  Gamma_ = Gamma_raw * (1.e3 / _c_);

  // Kinematics: exactly one of epsilon (mass retention) or vkick (v/c).
  auto eps_opt   = input.get<double>("epsilon");
  auto vkick_opt = input.get<double>("vkick");
  if (eps_opt.has_value() == vkick_opt.has_value()) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of epsilon, vkick");
  }
  if (eps_opt) {
    eps_ = *eps_opt;
    if (eps_ < 0. || eps_ > 1. - 1e-9) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': epsilon must be in [0, 1 - 1e-9]; for smaller kicks "
                                  "input vkick directly");
    }
    vkick_ = std::sqrt((1. - eps_) * (1. + eps_));
  }
  else {
    vkick_ = *vkick_opt;
    if (vkick_ < 1e-8 || vkick_ > 1.) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': vkick must be in [1e-8, 1] (in units of c); for "
                                  "vkick -> 0 the model is plain CDM");
    }
    eps_ = std::sqrt((1. - vkick_) * (1. + vkick_));
  }
  // p/m_d = v/sqrt(1-v^2) = v/eps  =>  M = m_d/T0 = kQKick * eps / v.
  M_ = kQKick * eps_ / vkick_;

  // Momentum grid: uniform in u = ln q on [ln q_lo, ln q_kick + 3 sigma], so the
  // kernel keeps a 3 sigma margin above q_cut(a=1) = kQKick for any bin count:
  //   du (n_bins - 3 kernel_width) = ln(kQKick / q_lo).
  n_bins_       = input.get_or("momenta_bins", 96);
  q_min_ratio_  = input.get_or("q_min_ratio", 1e-4);
  kernel_width_ = input.get_or("kernel_width", 1.0);
  l_max_input_  = input.get_or("l_max", -1);
  if (n_bins_ < 8 || n_bins_ > 4096) {
    throw std::invalid_argument("species '" + instance_name + "': momenta_bins out of range");
  }
  if (!(q_min_ratio_ > 0. && q_min_ratio_ < 0.5)) {
    throw std::invalid_argument("species '" + instance_name + "': q_min_ratio must be in (0, 0.5)");
  }
  if (!(kernel_width_ >= 0.2 && kernel_width_ <= 5.)) {
    throw std::invalid_argument("species '" + instance_name +
                                "': kernel_width must be in [0.2, 5] bins");
  }

  const double q_lo = kQKick * q_min_ratio_;
  du_               = std::log(kQKick / q_lo) / (n_bins_ - 3. * kernel_width_);
  q_.resize(n_bins_);
  u_.resize(n_bins_);
  dq_.resize(n_bins_);
  for (int i = 0; i < n_bins_; ++i) {
    u_[i]  = std::log(q_lo) + (i + 0.5) * du_;
    q_[i]  = std::exp(u_[i]);
    dq_[i] = q_[i] * du_;
  }
  q_bg_ = q_;
  w_.assign(n_bins_, 0.);
  w_bg_.assign(n_bins_, 0.);        // refreshed from f each ComputeBackground
  dlnf0_dlnq_.assign(n_bins_, 0.);  // unused (tensor slots disabled)
  scratch_G_.assign(n_bins_, 0.);

  // Normalization: Omega_ini/omega_ini (parent initial abundance, dcdm convention)
  // XOR Omega_dcdmwdm/omega_dcdmwdm (combined sector today). During a shooting
  // build DoShooting legitimately writes <instance>.Omega_dcdmwdm alongside the
  // user's initial key.
  auto Omega_ini_opt  = input.get<double>("Omega_ini");
  auto omega_ini_opt  = input.get<double>("omega_ini");
  auto Omega_comb_opt = input.get<double>("Omega_dcdmwdm");
  auto omega_comb_opt = input.get<double>("omega_dcdmwdm");
  if (Omega_ini_opt && omega_ini_opt) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify only one of Omega_ini, omega_ini");
  }
  if (Omega_comb_opt && omega_comb_opt) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify only one of Omega_dcdmwdm, omega_dcdmwdm");
  }
  const bool has_initial  = Omega_ini_opt.has_value() || omega_ini_opt.has_value();
  const bool has_combined = Omega_comb_opt.has_value() || omega_comb_opt.has_value();
  const bool in_shooting  = (pfc != nullptr) && pfc->is_shooting;
  if (has_initial && has_combined && !in_shooting) {
    throw std::invalid_argument("species '" + instance_name +
                                "': Omega_ini conflicts with Omega_dcdmwdm — choose one");
  }
  if (!has_initial && !has_combined) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify Omega_ini/omega_ini or Omega_dcdmwdm/omega_dcdmwdm");
  }
  if (Omega_ini_opt)
    Omega_ini_pending_ = *Omega_ini_opt;
  else if (omega_ini_opt)
    Omega_ini_pending_ = *omega_ini_opt / settings.h / settings.h;
  if (Omega_comb_opt)
    Omega_combined_pending_ = *Omega_comb_opt;
  else if (omega_comb_opt)
    Omega_combined_pending_ = *omega_comb_opt / settings.h / settings.h;
}

// ─────────────────────────────────────────────────────────────────────────────
// FillInjection
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::FillInjection(double a,
                                           double rho_dcdm,
                                           double* J,
                                           double* dJdlnq) const {
  const int N        = n_bins_;
  const double q_cut = a * kQKick;
  const double u_cut = std::log(q_cut);
  const double sigma = kernel_width_ * du_;

  auto zero_all = [&]() {
    for (int i = 0; i < N; ++i) {
      J[i] = 0.;
      if (dJdlnq != nullptr)
        dJdlnq[i] = 0.;
    }
  };

  // Injection is dropped while the cutoff is below the grid: the decayed
  // fraction before a = q_min_ratio is negligible by construction (guarded in
  // DCDM_WDM_Species::CreateAll for very short lifetimes).
  if (u_cut < u_[0] || rho_dcdm <= 0. || Gamma_ == 0.) {
    zero_all();
    return;
  }

  double D = 0.;
  for (int i = 0; i < N; ++i) {
    const double x       = (u_[i] - u_cut) / sigma;
    const double G       = std::exp(-0.5 * x * x);
    scratch_G_[i]        = G;
    const double epsilon = std::sqrt(q_[i] * q_[i] + a * a * M_ * M_);
    D                   += dq_[i] * q_[i] * q_[i] * epsilon * G;
  }
  if (D <= 0.) {
    zero_all();
    return;
  }
  const double a2 = a * a;
  const double A  = a * Gamma_ * rho_dcdm * a2 * a2 / factor_;
  for (int i = 0; i < N; ++i) {
    J[i] = A * scratch_G_[i] / D;
    if (dJdlnq != nullptr)
      dJdlnq[i] = J[i] * (u_cut - u_[i]) / (sigma * sigma);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Info prints
// ─────────────────────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::PrintNeffInfo() const {
  printf(" -> dcdm_wdm daughter '%s': %d momentum bins, vkick = %g c (epsilon = %g), "
         "Gamma = %g Mpc^-1\n",
         name().c_str(),
         n_bins_,
         vkick_,
         eps_,
         Gamma_);
}

void WdmDecayProductSpecies::PrintMassInfo() const {
  // The daughter mass in eV is not physical here (internal kQKick convention);
  // report the physical kinematics instead.
  printf(" -> dcdm_wdm daughter '%s': m_d/m_parent = %g, vkick = %g c\n",
         name().c_str(),
         eps_ / 2.,
         vkick_);
}
```

- [ ] **Step 3: Write the failing test** `species/dcdm_wdm_test.cpp`:

```cpp
// Unit test for WdmDecayProductSpecies (dcdm_wdm daughter): construction guards,
// momentum grid, injection-kernel normalization (the exact energy sum rule), and
// (from Task 2 onward) the background interface.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "background.h"
#include "parser.h"
#include "species/wdm_decay_product.h"

static NcdmSettings TestSettings() {
  NcdmSettings s{};
  s.h          = 0.67;
  s.T_cmb      = 2.7255;
  s.tol_ncdm   = 1e-3;
  s.tol_ncdm_bg = 1e-5;
  s.tol_M_ncdm = 1e-8;
  return s;
}

static FileContent BaseFc() {
  FileContent fc;
  fc.set("ddm.type", "dcdm_wdm");
  fc.set("ddm.Gamma", "100");
  fc.set("ddm.vkick", "0.1");
  fc.set("ddm.Omega_ini", "0.05");
  return fc;
}

static bool Throws(FileContent& fc) {
  background pba{};
  pba.H0 = 2.2e-4;
  try {
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    return false;
  }
  catch (const std::invalid_argument&) {
    return true;
  }
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  // ── Construction + grid ─────────────────────────────────────────────────────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    assert(sp.q_size() == 96);
    // M = kQKick * eps / v with eps = sqrt(1 - v^2)
    const double v = 0.1, eps = std::sqrt(1. - v * v);
    assert(std::fabs(sp.M() - 10.0 * eps / v) < 1e-12 * sp.M());
    // Grid top edge: u.back() + du/2 = ln(kQKick) + 3*sigma (3-sigma margin at a=1)
    const double du = sp.u()[1] - sp.u()[0];
    assert(std::fabs((sp.u().back() + 0.5 * du) - (std::log(10.0) + 3.0 * du)) < 1e-10);
    // Grid bottom edge at q_min_ratio: u.front() - du/2 = ln(kQKick * 1e-4)
    assert(std::fabs((sp.u().front() - 0.5 * du) - std::log(10.0 * 1e-4)) < 1e-10);
  }
  {
    FileContent fc;
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.Gamma", "100");
    fc.set("ddm.epsilon", "0.6");
    fc.set("ddm.Omega_ini", "0.05");
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    assert(std::fabs(sp.vkick() - 0.8) < 1e-14);
    assert(std::fabs(sp.M() - 10.0 * 0.6 / 0.8) < 1e-12);
    assert(sp.GetOmega0() == 0.);  // decay product reserves nothing itself
  }

  // ── Guards ──────────────────────────────────────────────────────────────────
  {
    FileContent fc = BaseFc();
    fc.set("ddm.epsilon", "0.5");  // both epsilon and vkick -> reject
    assert(Throws(fc));
  }
  {
    FileContent fc;
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.Gamma", "100");
    fc.set("ddm.epsilon", "1.0");  // too close to 1
    fc.set("ddm.Omega_ini", "0.05");
    assert(Throws(fc));
  }
  {
    FileContent fc = BaseFc();
    fc.set("ddm.m", "0.1");  // thermal-NCDM key -> reject
    assert(Throws(fc));
  }
  {
    FileContent fc;  // no Gamma at all
    fc.set("ddm.type", "dcdm_wdm");
    fc.set("ddm.vkick", "0.1");
    fc.set("ddm.Omega_ini", "0.05");
    assert(Throws(fc));
  }
  {
    FileContent fc = BaseFc();
    fc.set("ddm.Omega_dcdmwdm", "0.1");  // conflicts with Omega_ini outside shooting
    assert(Throws(fc));
  }

  // ── Injection kernel: exact energy sum rule ────────────────────────────────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();
    std::vector<double> J(N), dJ(N);
    for (double a : {2e-4, 1e-3, 0.05, 0.3, 1.0}) {
      const double rho_dcdm = 0.7;  // arbitrary
      sp.FillInjection(a, rho_dcdm, J.data(), dJ.data());
      double sum = 0.;
      for (int i = 0; i < N; ++i) {
        const double q       = sp.q()[i];
        const double epsilon = std::sqrt(q * q + a * a * sp.M() * sp.M());
        sum                 += sp.dq()[i] * q * q * epsilon * J[i];
      }
      const double lhs = sum * sp.factor() / (a * a * a * a);
      const double rhs = a * sp.Gamma() * rho_dcdm;
      assert(std::fabs(lhs - rhs) < 1e-12 * rhs);
    }
    // Below the grid: no injection.
    sp.FillInjection(5e-5, 0.7, J.data(), dJ.data());
    for (int i = 0; i < N; ++i)
      assert(J[i] == 0.);
  }

  std::printf("dcdm_wdm daughter test passed\n");
  return 0;
}
```

- [ ] **Step 4: Register the build targets.**
  In `CMakeLists.txt`: add `species/wdm_decay_product.cpp` to the `CLASS_SPECIES_FILES` list (alphabetical position, near `species/ultra_relativistic.cpp`); add `add_executable(test-dcdm-wdm species/dcdm_wdm_test.cpp)` next to the other tests AND add `test-dcdm-wdm` to the `foreach(_t IN ITEMS ...)` list.
  In `Makefile`: add `test-dcdm-wdm \` to `TEST_TARGETS`.

- [ ] **Step 5: Build + run the test.**
  Run: `make test-dcdm-wdm && ./build/cmake/test-dcdm-wdm`
  Expected: `dcdm_wdm daughter test passed`. Iterate on compile errors (FileContent API details, include paths) until green. Also run `make class` to confirm the library still builds.

- [ ] **Step 6: Commit**

```bash
git add species/wdm_decay_product.h species/wdm_decay_product.cpp species/dcdm_wdm_test.cpp CMakeLists.txt Makefile
git commit -m "dcdm_wdm daughter species: construction, momentum grid, injection kernel"
```

---

### Task 2: Daughter background interface

**Files:**
- Modify: `species/wdm_decay_product.cpp` (replace Task-1 placeholder bodies)
- Modify: `species/dcdm_wdm_test.cpp` (extend)

**Interfaces:**
- Consumes: Task 1 members. Produces working `RegisterBackgroundIndices`, `RegisterIntegrationIndices`, `SetBackgroundInitialConditions`, `ComputeBackground`, `PPrime`, `RhoDotOverRho`, `WriteBackgroundColumnTitles/Data`.
- pvecback slot order (contiguous, claimed in this order): number, rho (=index_bg_rho_), p (=index_bg_p_), pseudo_p, f[N], dfdlnq[N], inj[N]. bi slots: f[N], dfdlnq[N]. The composite (Task 3) fills the inj[N] pvecback columns; ComputeBackground here does NOT.

- [ ] **Step 1: Implement the background methods** (replace placeholders):

```cpp
// ── Background ───────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_   = index_bg++;
  index_bg_rho_      = index_bg++;
  index_bg_p_        = index_bg++;
  index_bg_pseudo_p_ = index_bg++;
  index_bg_f_        = index_bg;
  index_bg          += n_bins_;
  index_bg_dfdlnq_   = index_bg;
  index_bg          += n_bins_;
  index_bg_inj_      = index_bg;  // filled by DCDM_WDM_Species (needs rho_dcdm)
  index_bg          += n_bins_;
}

void WdmDecayProductSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_f_       = index_bi;
  index_bi         += n_bins_;
  index_bi_dfdlnq_  = index_bi;
  index_bi         += n_bins_;
}

void WdmDecayProductSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  for (int i = 0; i < n_bins_; ++i) {
    ctx.pvecback_integration[index_bi_f_ + i]      = 0.;
    ctx.pvecback_integration[index_bi_dfdlnq_ + i] = 0.;
  }
}

void WdmDecayProductSpecies::ComputeBackground(double a,
                                               const double* pvecback_B,
                                               double* pvecback) {
  const double z = 1. / a - 1.;
  for (int i = 0; i < n_bins_; ++i) {
    // Clamp tiny negative excursions from the ODE (f is a density).
    const double f                     = std::max(pvecback_B[index_bi_f_ + i], 0.);
    w_bg_[i]                           = f * dq_[i];
    pvecback[index_bg_f_ + i]          = f;
    pvecback[index_bg_dfdlnq_ + i]     = pvecback_B[index_bi_dfdlnq_ + i];
  }
  double n_wdm, rho_wdm, p_wdm, pseudo_p_wdm;
  ComputeMomenta(z, &n_wdm, &rho_wdm, &p_wdm, nullptr, &pseudo_p_wdm);
  pvecback[index_bg_number_]   = n_wdm;
  pvecback[index_bg_rho_]      = rho_wdm;
  pvecback[index_bg_p_]        = p_wdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_wdm;
  // index_bg_inj_ columns are written by DCDM_WDM_Species::ComputeBackground
  // after the child loop (they need rho_dcdm).
}

double WdmDecayProductSpecies::PPrime(double a,
                                      double H,
                                      const double* /*pvecback_B*/,
                                      const double* pvecback) const {
  // Standard NCDM redshift part + injection contribution to dP/dtau:
  //   dP_inj = factor/(3 a^4) * sum dq q^4/eps * J
  double dp_inj = 0.;
  for (int i = 0; i < n_bins_; ++i) {
    const double q       = q_[i];
    const double q2      = q * q;
    const double epsilon = std::sqrt(q2 + a * a * M_ * M_);
    dp_inj              += dq_[i] * q2 * q2 / epsilon * pvecback[index_bg_inj_ + i];
  }
  const double a2 = a * a;
  return a * H * (pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_]) +
         factor_ / (3. * a2 * a2) * dp_inj;
}

double WdmDecayProductSpecies::RhoDotOverRho(const double* pvecback,
                                             double a_prime_over_a) const {
  const double rho = Rho(pvecback);
  if (rho <= 0.)
    return 0.;
  const double a   = pvecback[bgm_->index_bg_a_];
  const double inj = (parent_ != nullptr) ? a * Gamma_ * parent_->Rho(pvecback) : 0.;
  return (-3. * a_prime_over_a * (rho + P(pvecback)) + inj) / rho;
}

// ── Background output columns ────────────────────────────────────────────────

void WdmDecayProductSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  const std::string& nm = name();
  w.Add("(.)rho_wdm_" + nm, 0.);
  w.Add("(.)p_wdm_" + nm, 0.);
  for (int i = 0; i < n_bins_; ++i)
    w.Add("f_" + nm + "[" + std::to_string(i) + "]", 0.);
}

void WdmDecayProductSpecies::WriteBackgroundData(const double* pvecback,
                                                 BackgroundColumnWriter& w) const {
  const std::string& nm = name();
  w.Add("(.)rho_wdm_" + nm, Rho(pvecback));
  w.Add("(.)p_wdm_" + nm, P(pvecback));
  for (int i = 0; i < n_bins_; ++i)
    w.Add("f_" + nm + "[" + std::to_string(i) + "]", pvecback[index_bg_f_ + i]);
}
```

Note: `ComputeMomenta` signature is `(z, n, rho, p, drho_dM, pseudo_p)` — pass `nullptr` for drho_dM as shown.

- [ ] **Step 2: Extend the test** — append to `main()` in `species/dcdm_wdm_test.cpp` (before the final printf):

```cpp
  // ── Background plumbing: indices, ICs, ComputeMomenta via w_bg_ ────────────
  {
    FileContent fc = BaseFc();
    WdmDecayProductSpecies sp(&fc, "ddm", TestSettings(), &pba, nullptr);
    const int N = sp.q_size();

    int index_bg = 0, index_bi = 0;
    sp.RegisterBackgroundIndices(index_bg);
    sp.RegisterIntegrationIndices(index_bi);
    assert(index_bg == 4 + 3 * N);
    assert(index_bi == 2 * N);

    std::vector<double> bi(index_bi, 0.), bg(index_bg, 0.);
    BackgroundICContext icc{};
    icc.a_ini                 = 1e-14;
    icc.pvecback_integration  = bi.data();
    sp.SetBackgroundInitialConditions(icc);

    // Empty distribution: rho = p = 0 exactly.
    sp.ComputeBackground(0.5, bi.data(), bg.data());
    assert(sp.Rho(bg.data()) == 0.);
    assert(sp.P(bg.data()) == 0.);

    // Populate one bin by hand and check rho against the closed-form quadrature.
    const int iq    = N / 2;
    const double a  = 0.5;
    bi[sp.bi_f_index() + iq] = 3.14;
    sp.ComputeBackground(a, bi.data(), bg.data());
    const double q       = sp.q()[iq];
    const double epsilon = std::sqrt(q * q + a * a * sp.M() * sp.M());
    const double rho_exp = sp.factor() / std::pow(a, 4) * sp.dq()[iq] * q * q * epsilon * 3.14;
    assert(std::fabs(sp.Rho(bg.data()) - rho_exp) < 1e-12 * rho_exp);
  }
```

- [ ] **Step 3: Build + run.**
  Run: `make test-dcdm-wdm && ./build/cmake/test-dcdm-wdm`
  Expected: `dcdm_wdm daughter test passed`.

- [ ] **Step 4: Commit**

```bash
git add species/wdm_decay_product.cpp species/dcdm_wdm_test.cpp
git commit -m "dcdm_wdm daughter: background distribution integration + momenta"
```

---

### Task 3: DCDM_WDM composite, factory registration, shooting, background golden test

**Files:**
- Create: `species/dcdm_wdm_species.h`
- Create: `species/dcdm_wdm_species.cpp`
- Modify: `species/all_species.h` (include + factory row)
- Modify: `CMakeLists.txt` (`species/dcdm_wdm_species.cpp` into CLASS_SPECIES_FILES)
- Modify: `species/dcdm_wdm_test.cpp` (gauge guard + full-background golden test)

**Interfaces:**
- Produces `DCDM_WDM_Species : CompositeSpecies` with `kTypeName = "dcdm_wdm"`, `static std::vector<Named> CreateAll(const SpeciesBuildContext&)`, children order `{kDcdm=0, kWdm=1}`, typed layout accessors `dcdm_layout(...)`/`wdm_layout(...)`, accessors `dcdm()`/`wdm()`.
- The composite's `GetOmega0()` returns the pinned/shot combined `Omega_dcdmwdm` (fallback: child sum). Shooting: mode (a) `Omega_ini` given → fixed-point on `<inst>.Omega_dcdmwdm`; mode (b) `Omega_dcdmwdm` given → shoot `<inst>.Omega_ini`.
- Consumes Task 1/2 daughter API + `DCDMSpecies` (ctor `(const background&, omega0_dcdmdr, Gamma, Omega_ini)`, `Rho`, `Gamma_dcdm()`, `bi_rho_index()`, `transfer_delta_index()`, `transfer_theta_index()`).

- [ ] **Step 1: Write `species/dcdm_wdm_species.h`:**

```cpp
#pragma once
#include <memory>
#include <string_view>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "dcdm.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"
#include "wdm_decay_product.h"

class BackgroundModule;

/**
 * DCDM_WDM_Species: composite for cold dark matter decaying to two identical
 * massive daughters (arXiv:2606.14849). Children: DCDMSpecies parent (index
 * kDcdm) + WdmDecayProductSpecies daughter (kWdm).
 *
 * EnergyType::Other so the background loop splits rho correctly:
 *   rho_m += rho - 3p, rho_r += 3p.
 *
 * Background wiring: ComputeBackground writes the daughter's injection columns
 * after the child loop (needs rho_dcdm); BackgroundDerivs adds df_i/dtau = J_i
 * and dg_i/dtau = dJ_i/dlnq. AddCouplingDerivs adds the perturbation injection
 * source J_i * delta_dcdm to the daughter's psi_0 bins.
 *
 * Synchronous gauge only (guard in CreateAll). The l=1 parent-velocity source
 * is omitted: theta_dcdm == 0 identically in synchronous gauge.
 */
class DCDM_WDM_Species : public CompositeSpecies {
 public:
  static constexpr std::string_view kTypeName = "dcdm_wdm";

  // Uses the generic CompositeSpecies::PerturbLayout (one owning sub-layout per
  // child in children_ order) — required by the TallyStressEnergy/DelegateTally
  // hot path (#358). Typed views below.
  enum ChildIndex { kDcdm = 0, kWdm = 1 };
  static const DCDMSpecies::PerturbLayout& dcdm_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const DCDMSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kDcdm]);
  }
  static const NCDMBaseSpecies::PerturbLayout& wdm_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kWdm]);
  }

  /** Takes ownership of a pre-built daughter (from CreateAll). */
  DCDM_WDM_Species(std::unique_ptr<WdmDecayProductSpecies> wdm,
                   const background* pba,
                   const BackgroundModule* bgm,
                   double omega0_combined,
                   double Omega_ini_dcdm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  void SetBackgroundModule(const BackgroundModule* bgm) override;

  // Injection wiring (see class doc)
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  /** Combined sector density today (pinned/shot Omega_dcdmwdm), mirroring
   *  DNCDM_DR_Species::GetOmega0: at closure time the daughter's density is
   *  not yet integrated, so the child sum alone would under-reserve. */
  double GetOmega0() const override;

  // ── Shooter hooks ──────────────────────────────────────────────────────────
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  // Registration, derivs (children + AddCouplingDerivs), ICs, stress-energy and
  // switch copies use the generic CompositeSpecies child loops. FillSources is
  // overridden to also write the parent's transfer sources (DCDMSpecies itself
  // has none — mirrors DCDM_DR).
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  DCDMSpecies& dcdm() {
    return *dcdm_;
  }
  const DCDMSpecies& dcdm() const {
    return *dcdm_;
  }
  WdmDecayProductSpecies& wdm() {
    return *wdm_;
  }
  const WdmDecayProductSpecies& wdm() const {
    return *wdm_;
  }

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  DCDMSpecies* dcdm_           = nullptr;  // non-owning pointers into children_
  WdmDecayProductSpecies* wdm_ = nullptr;
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;

  ShootingTarget shooting_target_{};
  bool needs_shooting_ = false;

  // Scratch for BackgroundDerivs (background integration is single-threaded).
  mutable std::vector<double> scratch_J_, scratch_dJ_;
};
```

- [ ] **Step 2: Write `species/dcdm_wdm_species.cpp`:**

```cpp
#include "dcdm_wdm_species.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "perturbations_module.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

DCDM_WDM_Species::DCDM_WDM_Species(std::unique_ptr<WdmDecayProductSpecies> wdm_arg,
                                   const background* pba,
                                   const BackgroundModule* bgm,
                                   double omega0_combined,
                                   double Omega_ini_dcdm)
    : CompositeSpecies(wdm_arg->name(), BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm) {
  auto dcdm = std::make_unique<DCDMSpecies>(*pba, omega0_combined, wdm_arg->Gamma(),
                                            Omega_ini_dcdm);
  dcdm_ = dcdm.get();
  wdm_  = wdm_arg.get();
  wdm_->SetParent(dcdm_);
  // Child order matters: DCDM first so its rho is in pvecback before the
  // daughter's injection columns are written (ComputeBackground below).
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(wdm_arg));
  scratch_J_.assign(wdm_->q_size(), 0.);
  scratch_dJ_.assign(wdm_->q_size(), 0.);
}

void DCDM_WDM_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

// ─────────────────────────────────────────────────────────────────────────────
// Background wiring
// ─────────────────────────────────────────────────────────────────────────────

void DCDM_WDM_Species::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  CompositeSpecies::ComputeBackground(a, pvecback_B, pvecback);
  // Injection columns need rho_dcdm, so they are written after the child loop.
  wdm_->FillInjection(a, dcdm_->Rho(pvecback), pvecback + wdm_->bg_inj_index(), nullptr);
}

void DCDM_WDM_Species::BackgroundDerivs(double tau,
                                        const double* y,
                                        double* dy,
                                        const double* pvecback) {
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);  // dcdm: rho' (dilution+decay)

  const double a = pvecback[bgm_->index_bg_a_];
  wdm_->FillInjection(a, dcdm_->Rho(pvecback), scratch_J_.data(), scratch_dJ_.data());
  const int N     = wdm_->q_size();
  const int bi_f  = wdm_->bi_f_index();
  const int bi_g  = wdm_->bi_dfdlnq_index();
  for (int i = 0; i < N; ++i) {
    dy[bi_f + i] = scratch_J_[i];
    dy[bi_g + i] = scratch_dJ_[i];
  }
}

void DCDM_WDM_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm_" + name(), 0.);
  wdm_->WriteBackgroundColumnTitles(w);
}

void DCDM_WDM_Species::WriteBackgroundData(const double* pvecback,
                                           BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm_" + name(), dcdm_->Rho(pvecback));
  wdm_->WriteBackgroundData(pvecback, w);
}

// ─────────────────────────────────────────────────────────────────────────────
// Closure + shooter hooks
// ─────────────────────────────────────────────────────────────────────────────

double DCDM_WDM_Species::GetOmega0() const {
  if (wdm_->Omega_combined_pending().has_value())
    return *wdm_->Omega_combined_pending();
  return CompositeSpecies::GetOmega0();  // pre-shooting discovery fallback
}

std::vector<ShootingTarget> DCDM_WDM_Species::GetShootingTargets() const {
  if (needs_shooting_)
    return {shooting_target_};
  return {};
}

void DCDM_WDM_Species::ComputeShootingGuess(const SpeciesBuildContext& /*ctx*/,
                                            std::vector<double>& guess,
                                            std::vector<double>& dxdy) const {
  const std::string& t = shooting_target_.target_name;
  const double tv      = shooting_target_.target_value;
  if (t == name() + ".Omega_dcdmwdm_fixedpoint") {
    // Initial mode: seed the combined-today reserve with Omega_ini (matter
    // scaling would make them equal; kicks only redshift the kinetic part).
    guess.push_back(tv);
    dxdy.push_back(1.0);
  }
  else {
    // Combined mode: shoot Omega_ini; for slow decays combined-today ~ Omega_ini.
    // For fast decays (Gamma >> H0) the kinetic fraction has redshifted away and
    // combined-today -> eps * Omega_ini, so seed accordingly.
    const double gamma_over_H0 = wdm_->Gamma() / pba_->H0;
    const double eps           = std::max(wdm_->epsilon_retention(), 1e-3);
    const double seed          = (gamma_over_H0 > 1.) ? tv / eps : tv;
    guess.push_back(seed);
    dxdy.push_back(seed / std::max(tv, 1e-30));
  }
}

double DCDM_WDM_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                 const ShootingTarget& target) const {
  const double* bg      = ctx.bg_today;
  const double H0       = ctx.pba->H0;
  const double combined = (dcdm_->Rho(bg) + wdm_->Rho(bg)) / (H0 * H0);
  if (target.target_name == name() + ".Omega_dcdmwdm_fixedpoint") {
    // Initial mode: drive the reserved Omega_dcdmwdm (= GetOmega0()) to the
    // integrated combined density.
    return -combined + GetOmega0();
  }
  // Combined mode: drive the integrated combined to the requested target.
  return combined - target.target_value;
}

// ─────────────────────────────────────────────────────────────────────────────
// Factory
// ─────────────────────────────────────────────────────────────────────────────

std::vector<Named> DCDM_WDM_Species::CreateAll(const SpeciesBuildContext& ctx) {
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  std::vector<Named> result;
  if (instances.empty())
    return result;

  // Synchronous gauge only (the l=1 injection source and the gauge transform of
  // the injected distribution are not implemented; the reference paper works in
  // synchronous gauge). Mirror the Type3 guard: any gauge value starting with
  // "new" (case-insensitive) is Newtonian.
  if (auto gauge = ctx.pfc->get<std::string>("gauge")) {
    std::string g = *gauge;
    std::transform(g.begin(), g.end(), g.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (g.rfind("new", 0) == 0) {
      throw std::invalid_argument(
          "dcdm_wdm supports the synchronous gauge only (arXiv:2606.14849 implementation)");
    }
  }

  result.reserve(instances.size());
  for (const auto& nm : instances) {
    (void) ctx.pfc->get<std::string>(nm + ".type");  // mark consumed

    auto wdm = std::make_unique<WdmDecayProductSpecies>(ctx.pfc, nm, *ctx.ncdm_settings,
                                                        ctx.pba, ctx.bgm);

    const bool has_initial  = wdm->Omega_ini_pending().has_value();
    const bool has_combined = wdm->Omega_combined_pending().has_value();
    const double Omega_ini  = wdm->Omega_ini_pending().value_or(0.);
    double omega0_combined  = wdm->Omega_combined_pending().value_or(0.);

    ShootingTarget target{};
    bool needs_shooting = false;
    const bool in_shooting = ctx.pfc->is_shooting;
    if (has_initial) {
      // Mode (a): fixed-point shoot of the combined-today reserve. During a
      // shooting iteration the shot value is already present (has_combined).
      target         = {nm + ".Omega_dcdmwdm_fixedpoint", nm + ".Omega_dcdmwdm",
                        has_combined ? omega0_combined : Omega_ini};
      needs_shooting = true;
      if (!has_combined)
        omega0_combined = Omega_ini;  // discovery-build seed for the reserve
    }
    else if (has_combined && !in_shooting) {
      // Mode (b): shoot the parent's initial abundance to hit combined-today.
      target         = {nm + ".Omega_dcdmwdm", nm + ".Omega_ini", omega0_combined};
      needs_shooting = true;
    }

    // Mode (b) discovery build needs an Omega_ini seed for the DCDM child.
    double Omega_ini_dcdm = Omega_ini;
    if (!has_initial) {
      std::vector<double> g, d;
      auto tmp = std::make_unique<DCDM_WDM_Species>(
          std::make_unique<WdmDecayProductSpecies>(ctx.pfc, nm, *ctx.ncdm_settings, ctx.pba,
                                                   ctx.bgm),
          ctx.pba, ctx.bgm, omega0_combined, /*Omega_ini_dcdm=*/omega0_combined);
      tmp->shooting_target_ = target;
      tmp->needs_shooting_  = true;
      tmp->ComputeShootingGuess(ctx, g, d);
      Omega_ini_dcdm = g[0];
    }

    auto composite = std::make_unique<DCDM_WDM_Species>(std::move(wdm), ctx.pba, ctx.bgm,
                                                        omega0_combined, Omega_ini_dcdm);
    // Seed the reserve for the discovery build in mode (a): GetOmega0() falls
    // back to the child sum (= omega0_combined via the DCDM child) until the
    // shooter pins <inst>.Omega_dcdmwdm.
    composite->shooting_target_ = target;
    composite->needs_shooting_  = needs_shooting;
    result.push_back({nm, std::move(composite)});
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbation coupling: injection source on the daughter's psi_0 bins
// ─────────────────────────────────────────────────────────────────────────────

void DCDM_WDM_Species::AddCouplingDerivs(double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) const {
  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector* pv     = ppw->pv.get();

  const auto& dcdm_lay = dcdm_layout(*pv->species_layouts[collection_index_]);
  const auto& wdm_lay  = wdm_layout(*pv->species_layouts[collection_index_]);
  if (dcdm_lay.idx_delta < 0 || wdm_lay.q_size <= 0 || wdm_lay.index_per_q.empty())
    return;

  const double* pvecback  = ppw->pvecback.data();
  const double delta_dcdm = y[dcdm_lay.idx_delta];
  const int inj0          = wdm_->bg_inj_index();
  for (int iq = 0; iq < wdm_lay.q_size; ++iq)
    dy[wdm_lay.index_per_q[iq]] += pvecback[inj0 + iq] * delta_dcdm;
}

// ─────────────────────────────────────────────────────────────────────────────
// Sources / output (parent slots; the daughter writes its own via the generic
// child loop)
// ─────────────────────────────────────────────────────────────────────────────

void DCDM_WDM_Species::FillSources(const BaseSpecies::PerturbLayout& base,
                                   const double* y,
                                   const double* dy,
                                   PerturbSourceContext& ctx) const {
  CompositeSpecies::FillSources(base, y, dy, ctx);  // daughter's d_wdm/t_wdm

  PerturbationsModule* p_mod = ctx.p_mod;
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;
  const auto& dcdm_lay = dcdm_layout(base);

  if (dcdm_->transfer_delta_index() >= 0) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, dcdm_->transfer_delta_index(),
                          ctx.index_tau, ctx.index_k,
                          y[dcdm_lay.idx_delta] +
                              (3. * ctx.a_prime_over_a + ctx.a * dcdm_->Gamma_dcdm()) *
                                  ctx.theta_over_k2);  // N-body gauge correction
  }
  if (dcdm_->transfer_theta_index() >= 0) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, dcdm_->transfer_theta_index(),
                          ctx.index_tau, ctx.index_k,
                          y[dcdm_lay.idx_theta] + ctx.theta_shift);
  }
}

void DCDM_WDM_Species::WriteOutputColumns(PerturbColumnWriter& w,
                                          const PerturbationsModule& mod,
                                          file_format fmt,
                                          BaseSpecies::TransferColumnSection section) const {
  dcdm_->WriteOutputColumns(w, mod, fmt, section);
  wdm_->WriteOutputColumns(w, mod, fmt, section);
}
```

- [ ] **Step 3: Register the factory.** In `species/all_species.h`: add `#include "dcdm_wdm_species.h"` (alphabetical) and the row `SpeciesFactoryEntry{DCDM_WDM_Species::kTypeName, &DCDM_WDM_Species::CreateAll},` right after the `DNCDMSpecies` row. In `CMakeLists.txt` add `species/dcdm_wdm_species.cpp` to CLASS_SPECIES_FILES.

- [ ] **Step 4: Extend the test** — append to `species/dcdm_wdm_test.cpp`. Add includes `#include "source/cosmology.h"` (check the actual include path used by other consumers: `#include "cosmology.h"` works — source/ is on the include path), `#include "species/dcdm_wdm_species.h"`, `#include <cstring>`, `#include <sstream>`, `#include <vector>`. Then:

```cpp
  // ── Gauge guard ─────────────────────────────────────────────────────────────
  {
    for (const char* g : {"newtonian", "Newtonian", "new", "newt"}) {
      FileContent fc = BaseFc();
      fc.set("gauge", g);
      SpeciesBuildContext sctx{};
      sctx.pfc = &fc;
      bool threw = false;
      try {
        DCDM_WDM_Species::CreateAll(sctx);
      }
      catch (const std::exception& e) {
        threw = std::string(e.what()).find("synchronous gauge only") != std::string::npos;
      }
      assert(threw);
    }
  }

  // ── Full background run + exact energy-conservation golden ────────────────
  {
    FileContent fc = BaseFc();  // Gamma=100 km/s/Mpc, vkick=0.1, Omega_ini=0.05
    fc.set("h", "0.67");
    fc.set("Omega_cdm", "0.20");
    fc.set("output", "");  // background-only
    Cosmology cosmo(fc);
    auto& bgm = cosmo.GetBackgroundModule();

    // Shooting must have pinned the combined density: reserve == integrated.
    const double reserve = bgm->GetOmega0Species("ddm");
    assert(reserve > 0.01 && reserve < 0.06);

    // Pull the full background table via titles/data.
    std::string titles;
    bgm->background_output_titles(titles);
    std::vector<std::string> cols;
    {
      std::stringstream ss(titles);
      std::string t;
      while (std::getline(ss, t, '\t'))
        if (!t.empty())
          cols.push_back(t);
    }
    auto col = [&](const std::string& nm) {
      for (size_t i = 0; i < cols.size(); ++i)
        if (cols[i] == nm)
          return (int) i;
      std::printf("column '%s' not found\n", nm.c_str());
      assert(false);
      return -1;
    };
    const int n_titles = (int) cols.size();
    std::vector<double> data((size_t) n_titles * bgm->bt_size_);
    bgm->background_output_data(n_titles, data.data());

    const int ic_z    = col("z");
    const int ic_tau  = col("conf. time [Mpc]");
    const int ic_dcdm = col("(.)rho_dcdm_ddm");
    const int ic_wdm  = col("(.)rho_wdm_ddm");

    // rho_wdm(a) must equal the injection-time integral
    //   Int dtau' a' Gamma rho_dcdm(a') (a'/a)^3 sqrt(eps^2 + (1-eps^2)(a'/a)^2)
    // (each daughter is born with energy m_parent/2 and redshifts kinetically).
    const double v = 0.1, eps2 = 1. - v * v;
    const double Gamma = 100. * 1.e3 / _c_;  // same double conversion as the species
    const int n_rows = bgm->bt_size_;
    auto row = [&](int r, int c) { return data[(size_t) r * n_titles + c]; };
    for (int target_row : {n_rows / 2, (3 * n_rows) / 4, n_rows - 1}) {
      const double a_t = 1. / (1. + row(target_row, ic_z));
      double integral  = 0.;
      for (int r = 1; r <= target_row; ++r) {
        const double dtau = row(r, ic_tau) - row(r - 1, ic_tau);
        auto integrand    = [&](int rr) {
          const double ap = 1. / (1. + row(rr, ic_z));
          const double x  = ap / a_t;
          // energy ratio to birth energy m_parent/2: sqrt(eps^2 + (1-eps^2) x^2),
          // mass part eps2 constant, kinetic part (1-eps2) redshifts as x^2
          return ap * Gamma * row(rr, ic_dcdm) * x * x * x *
                 std::sqrt(eps2 + (1. - eps2) * x * x);
        };
        integral += 0.5 * (integrand(r - 1) + integrand(r)) * dtau;
      }
      const double rho_wdm = row(target_row, ic_wdm);
      if (rho_wdm > 0.) {
        const double rel = std::fabs(integral - rho_wdm) / rho_wdm;
        std::printf("energy conservation at z=%.3g: rel. dev. %.2e\n",
                    row(target_row, ic_z), rel);
        assert(rel < 5e-3);
      }
    }
  }
```

Also verify the exact background title strings at runtime — if `col()` fails on "z" or "conf. time [Mpc]", print `titles` and adapt to the actual names used by BackgroundModule::background_output_titles.

- [ ] **Step 5: Build + run.**
  Run: `make test-dcdm-wdm && ./build/cmake/test-dcdm-wdm` (from repo root)
  Expected: prints relative deviations ≲ few×10⁻³ and `dcdm_wdm daughter test passed`. Also `make class && ./class explanatory.ini` still runs (species absent → zero impact).
  If shooting fails to converge, debug `GetShootingTargets`/residual signs against the DNCDM_DR pattern (`species/dncdm_dr_species.cpp:308-364`).

- [ ] **Step 6: Run the full test suite:** `make test`
  Expected: all existing tests still pass.

- [ ] **Step 7: Commit**

```bash
git add species/dcdm_wdm_species.h species/dcdm_wdm_species.cpp species/all_species.h species/wdm_decay_product.h species/wdm_decay_product.cpp CMakeLists.txt species/dcdm_wdm_test.cpp
git commit -m "dcdm_wdm composite: factory, shooting, background injection wiring + energy-conservation golden"
```

---

### Task 4: Daughter perturbations + coupling

**Files:**
- Modify: `species/wdm_decay_product.cpp` (replace remaining placeholders)
- Modify: `species/dcdm_wdm_test.cpp` (extend: perturbed run smoke test)

**Interfaces:**
- Consumes: layout type `NCDMBaseSpecies::PerturbLayout` (l_max, q_size, index_per_q), `PerturbScalarContext` fields (`k, a2, metric_continuity, metric_euler, metric_shear, cotKgen`), `ppw->s_l`, pvecback columns from Task 2/3.
- Equations (synchronous gauge; ψ_ℓ = unnormalized (Δf)_ℓ per bin; g = ∂f₀/∂ln q from pvecback):
  - dψ₀ = −(qk/ε)ψ₁ + (metric_continuity/3)·g  [+ J·δ_dcdm from composite]
  - dψ₁ = (qk/3ε)(ψ₀ − 2s₂ψ₂) − (ε·metric_euler/(3qk))·g
  - dψ₂ = (qk/5ε)(2s₂ψ₁ − 3s₃ψ₃) − s₂·(2/15)·metric_shear·g
  - 3≤ℓ<l_max standard; truncation standard. (Identical structure to `DNCDMSpecies::PerturbDerivs`, with `dlnf0_dlnq` → `g`.)

- [ ] **Step 1: Implement** (replace placeholders in `species/wdm_decay_product.cpp`):

```cpp
// ── Perturbations ────────────────────────────────────────────────────────────

void WdmDecayProductSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                         perturb_vector* /*pv*/,
                                                         const precision* ppr,
                                                         int& index_pt,
                                                         const perturb_workspace* /*ppw*/,
                                                         int /*gauge*/) {
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);

  // No fluid approximation: always the full Boltzmann hierarchy (DNCDM precedent).
  layout.l_max  = (l_max_input_ > 0) ? l_max_input_ : ppr->l_max_ncdm;
  layout.q_size = q_size();

  layout.index_per_q.clear();
  layout.index_per_q.reserve(layout.q_size);
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  index_pt += layout.total_size();
}

void WdmDecayProductSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                           double /*tau*/,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback.data();
  const double M_wdm     = M_;
  const int lmax         = layout.l_max;

  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q = q_[iq];
    // Unnormalized slope g = df0/dlnq: multiplies exactly the terms that carry
    // dlnf0/dlnq in the normalized (Psi) hierarchy.
    const double g = pvecback[index_bg_dfdlnq_ + iq];

    const double epsilon        = std::sqrt(q * q + a2 * M_wdm * M_wdm);
    const double qk_div_epsilon = k * q / epsilon;
    const int idx               = layout.index_per_q[iq];

    dy[idx]     = -qk_div_epsilon * y[idx + 1] + metric_continuity * g / 3.;
    dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) -
                  epsilon * metric_euler / (3. * q * k) * g;
    dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                  s_l[2] * metric_shear * 2. / 15. * g;
    for (int l = 3; l < lmax; ++l)
      dy[idx + l] = qk_div_epsilon / (2. * l + 1.) *
                    (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
    dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] -
                     (1. + lmax) * k * cotKgen * y[idx + lmax];
  }
}

void WdmDecayProductSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;  // isocurvature: daughter starts empty and unperturbed

  // Bins already populated at tau_ini (very short lifetimes) carry the parent's
  // adiabatic density contrast; superhorizon injection shares delta_dcdm.
  const double* pvecback     = ctx.ppw->pvecback.data();
  const double delta_dcdm_ic = 0.75 * ctx.delta_g_ic;
  for (int iq = 0; iq < layout.q_size; ++iq)
    y[layout.index_per_q[iq]] = pvecback[index_bg_f_ + iq] * delta_dcdm_ic;
  // psi_{l>=1} start at zero.
}

BaseSpecies::StressEnergyContribution WdmDecayProductSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_p_];
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return se;

  const double a  = ppw->scalar_ctx.a;
  const double k  = ppw->scalar_ctx.k;
  const double a2 = a * a;

  double s_drho = 0., s_theta = 0., s_dp = 0., s_shear = 0.;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q       = q_[iq];
    const double q2      = q * q;
    const double epsilon = std::sqrt(q2 + a2 * M_ * M_);
    const int idx        = layout.index_per_q[iq];
    const double dqv     = dq_[iq];
    s_drho  += dqv * q2 * epsilon * y[idx];
    s_theta += dqv * q2 * q * y[idx + 1];
    s_dp    += dqv * q2 * q2 / epsilon * y[idx];
    s_shear += dqv * q2 * q2 / epsilon * y[idx + 2];
  }
  const double fac    = factor_ / (a2 * a2);
  se.delta_rho        = fac * s_drho;
  se.rho_plus_p_theta = fac * k * s_theta;
  se.delta_p          = fac * s_dp / 3.;
  se.rho_plus_p_shear = fac * (2. / 3.) * s_shear;
  return se;
}

void WdmDecayProductSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y,
    double* new_y,
    const PerturbSwitchContext& /*ctx*/) const {
  // The daughter has no approximation of its own: the layout shape never
  // changes, so this is always a slot-by-slot copy. (NOTE: NCDMBaseSpecies has
  // no default copy — without this override the hierarchy would be zeroed at
  // every TCA/RSA switch.)
  const auto& old_l = static_cast<const NCDMBaseSpecies::PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const NCDMBaseSpecies::PerturbLayout&>(new_base);
  if (old_l.q_size <= 0 || new_l.q_size <= 0)
    return;
  for (int iq = 0; iq < new_l.q_size; ++iq) {
    const int o = old_l.index_per_q[iq];
    const int n = new_l.index_per_q[iq];
    for (int l = 0; l <= new_l.l_max; ++l)
      new_y[n + l] = old_y[o + l];
  }
}

void WdmDecayProductSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                         const double* y,
                                         const double* /*dy*/,
                                         PerturbSourceContext& ctx) const {
  if (ctx.index_md != ctx.p_mod->index_md_scalars_)
    return;
  if (index_tp_delta_ < 0 && index_tp_theta_ < 0)
    return;

  const auto& layout     = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  const double* pvecback = ctx.ppw->pvecback.data();
  const double rho       = Rho(pvecback);
  const double p         = P(pvecback);
  const double a         = ctx.a;
  const double a2        = ctx.a2;
  const double k         = ctx.k;

  double s_drho = 0., s_theta = 0.;
  if (layout.q_size > 0 && !layout.index_per_q.empty()) {
    for (int iq = 0; iq < layout.q_size; ++iq) {
      const double q       = q_[iq];
      const double q2      = q * q;
      const double epsilon = std::sqrt(q2 + a2 * M_ * M_);
      const int idx        = layout.index_per_q[iq];
      s_drho  += dq_[iq] * q2 * epsilon * y[idx];
      s_theta += dq_[iq] * q2 * q * y[idx + 1];
    }
  }
  const double fac = factor_ / (a2 * a2);

  if (index_tp_delta_ >= 0) {
    // delta with N-body gauge correction: delta_Nb = delta - (rhodot/rho) theta_tot/k^2
    const double delta = (rho > 0.) ? fac * s_drho / rho : 0.;
    const double src   = (rho > 0.)
                             ? delta - RhoDotOverRho(pvecback, ctx.a_prime_over_a) *
                                           ctx.theta_over_k2
                             : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp_delta_, ctx.index_tau,
                              ctx.index_k, src);
  }
  if (index_tp_theta_ >= 0) {
    const double rho_plus_p = rho + p;
    const double theta      = (rho_plus_p > 0.) ? fac * k * s_theta / rho_plus_p : 0.;
    const double src        = (rho_plus_p > 0.) ? theta + ctx.theta_shift : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp_theta_, ctx.index_tau,
                              ctx.index_k, src);
  }
}

void WdmDecayProductSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                                const PerturbationsModule& mod,
                                                file_format fmt,
                                                BaseSpecies::TransferColumnSection section) const {
  if (fmt != file_format::class_format)
    return;
  const perturbs* ppt = mod.GetPerturbs();
  if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
    w.Add("d_wdm_" + name(), index_tp_delta_, index_tp_delta_ >= 0);
  if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
    w.Add("t_wdm_" + name(), index_tp_theta_, index_tp_theta_ >= 0);
}
```

Note: `RegisterTransferSourceIndices` and `MarkUsedInSources` are inherited from NCDMBaseSpecies unchanged (delta/theta slots; ℓ>2 masked from sources).

- [ ] **Step 2: Extend the test** — a full C_ℓ smoke run. Append to `main()`:

```cpp
  // ── Perturbed run smoke test: Cls finite, spectra suppressed vs LCDM-ish ───
  {
    FileContent fc = BaseFc();
    fc.set("h", "0.67");
    fc.set("Omega_cdm", "0.20");
    fc.set("output", "tCl,pCl,lCl,mPk");
    fc.set("lensing", "yes");
    fc.set("l_max_scalars", "600");   // keep the smoke test fast
    fc.set("P_k_max_h/Mpc", "1.");
    fc.set("ddm.momenta_bins", "32"); // fast
    Cosmology cosmo(fc);
    auto& spm = cosmo.GetSpectraModule();
    (void) spm;  // reaching here without a throw = the full pipeline ran
    std::printf("full perturbed dcdm_wdm pipeline ran\n");
  }
```

Check `cosmology.h` for the exact getter name (`GetSpectraModule`). If the run throws, capture the message — typical first-run issues: uninitialized layout in tensor modes (ensure `output` above has no tensors), source slot mismatches in FillSources, NaNs from division by rho=0 (guards above).

- [ ] **Step 3: Build + run.**
  Run: `make test-dcdm-wdm && ./build/cmake/test-dcdm-wdm`
  Expected: all asserts pass including the smoke run. Watch runtime — the smoke test with 32 bins should stay under ~1-2 min.

- [ ] **Step 4: Run full suite + ΛCDM non-regression.**
  Run: `make test`
  Run: `./class explanatory.ini` — must be unaffected (species absent).

- [ ] **Step 5: Commit**

```bash
git add species/wdm_decay_product.cpp species/dcdm_wdm_test.cpp
git commit -m "dcdm_wdm daughter perturbations: unnormalized psi hierarchy, injection source, sources/output"
```

---

### Task 5: Physics limit validation script

**Files:**
- Create: `scripts/validate_dcdm_wdm.py`
- Create: `test/dotsyntax_dcdm_wdm.ini` (needs `git add -f` — test .ini files are gitignored)

**Interfaces:**
- Consumes the classy Python wrapper (build via `make classy-pip-dev`; import `from classy import Class`; parameters passed as the same strings as .ini keys, e.g. `{"ddm.type": "dcdm_wdm", ...}`).
- Produces a script that exits non-zero on tolerance failure, printing a table of checks. This is the pre-notebook quantitative gate.

- [ ] **Step 1: Write `test/dotsyntax_dcdm_wdm.ini`:**

```ini
# DCDM -> two massive decay products (arXiv:2606.14849), dot-syntax example.
# f = Omega_ini_ddm / (Omega_ini_ddm + Omega_cdm) is the unstable DM fraction.
h = 0.6736
omega_b = 0.02237
Omega_cdm = 0.234
ddm.type = dcdm_wdm
ddm.Gamma = 100         # km/s/Mpc  (Gamma^-1 ~ 9.78/h Gyr / 100 ~ 0.1 Gyr... see notebook)
ddm.vkick = 0.02        # kick velocity in units of c; alternative: ddm.epsilon
ddm.Omega_ini = 0.026   # parent initial abundance (dcdm convention)
ddm.momenta_bins = 96
output = tCl,pCl,lCl,mPk
lensing = yes
gauge = synchronous
```

- [ ] **Step 2: Write `scripts/validate_dcdm_wdm.py`:**

```python
#!/usr/bin/env python3
"""Physics-limit validation for the dcdm_wdm species (arXiv:2606.14849).

Checks (all ratios vs a reference run):
  1. vkick -> 0 (1e-6): C_ell^TT and P(k) match LCDM with the same total
     omega_m budget to < 0.1% (daughters are effectively CDM).
  2. vkick = 1 (epsilon = 0, massless daughters): C_ell^TT matches the existing
     dcdm_dr implementation with identical Gamma/Omega_ini to < 0.3%.
  3. f -> 0 (Omega_ini -> 1e-6 of the DM): identical to LCDM to < 0.05%.
Run AFTER `make classy-pip-dev`. Exits non-zero on failure.
"""
import sys

import numpy as np
from classy import Class

BASE = {
    "h": 0.6736,
    "omega_b": 0.02237,
    "output": "tCl,pCl,lCl,mPk",
    "lensing": "yes",
    "l_max_scalars": 2500,
    "P_k_max_h/Mpc": 1.0,
    "gauge": "synchronous",
}
OMEGA_DM = 0.26  # total (stable + unstable) early-time DM density fraction


def run(extra):
    cosmo = Class()
    cosmo.set({**BASE, **extra})
    cosmo.compute()
    cls = cosmo.lensed_cl(2000)
    kk = np.logspace(-3, 0, 200)
    pk = np.array([cosmo.pk(k * cosmo.h(), 0.0) for k in kk])
    cosmo.struct_cleanup()
    cosmo.empty()
    return cls, kk, pk


def ddm_pars(f, Gamma_kms_Mpc, vkick, bins=96):
    return {
        "Omega_cdm": OMEGA_DM * (1.0 - f),
        "ddm.type": "dcdm_wdm",
        "ddm.Gamma": Gamma_kms_Mpc,
        "ddm.vkick": vkick,
        "ddm.Omega_ini": OMEGA_DM * f,
        "ddm.momenta_bins": bins,
    }


def maxrel(a, b, lmin=2):
    a, b = np.asarray(a)[lmin:], np.asarray(b)[lmin:]
    return np.max(np.abs(a / b - 1.0))


failures = []


def check(name, value, tol):
    status = "PASS" if value < tol else "FAIL"
    print(f"  {name:55s} {value:10.3e}  (tol {tol:g})  {status}")
    if value >= tol:
        failures.append(name)


print("== reference LCDM ==")
lcdm_cls, kk, lcdm_pk = run({"Omega_cdm": OMEGA_DM})

print("== check 1: vkick -> 0 is LCDM ==")
cls, _, pk = run(ddm_pars(f=0.3, Gamma_kms_Mpc=100.0, vkick=1e-6))
check("vkick=1e-6: max |dTT/TT| (l<=2000)", maxrel(cls["tt"], lcdm_cls["tt"]), 1e-3)
check("vkick=1e-6: max |dP/P| (k<=1 h/Mpc)", maxrel(pk, lcdm_pk, 0), 2e-3)

print("== check 2: vkick=1 (massless daughters) matches dcdm_dr ==")
cls_wdm, _, pk_wdm = run(ddm_pars(f=0.3, Gamma_kms_Mpc=100.0, vkick=1.0))
cls_dr, _, pk_dr = run({
    "Omega_cdm": OMEGA_DM * 0.7,
    "Omega_ini_dcdm": OMEGA_DM * 0.3,
    "Gamma_dcdm": 100.0,
})
check("vkick=1 vs dcdm_dr: max |dTT/TT|", maxrel(cls_wdm["tt"], cls_dr["tt"]), 3e-3)
check("vkick=1 vs dcdm_dr: max |dP/P|", maxrel(pk_wdm, pk_dr, 0), 5e-3)

print("== check 3: f -> 0 is LCDM ==")
cls, _, pk = run(ddm_pars(f=1e-5, Gamma_kms_Mpc=100.0, vkick=0.1))
check("f=1e-5: max |dTT/TT|", maxrel(cls["tt"], lcdm_cls["tt"]), 5e-4)

print("== check 4: physical suppression present at finite kick ==")
cls, _, pk = run(ddm_pars(f=0.1, Gamma_kms_Mpc=200.0, vkick=0.02))
supp = pk[-1] / lcdm_pk[-1] - 1.0
print(f"  P(k=1)/P_LCDM - 1 = {supp:+.3f} (expect clearly negative)")
if not (supp < -0.005):
    failures.append("suppression sign")

if failures:
    print("FAILED:", failures)
    sys.exit(1)
print("all dcdm_wdm physics-limit checks passed")
```

- [ ] **Step 3: Build classy and run.**
  Run: `make classy-pip-dev` then `python3 scripts/validate_dcdm_wdm.py`
  Expected: `all dcdm_wdm physics-limit checks passed`. Notes for debugging:
  - Check 2 compares two *different discretizations* of the same physics (dark-radiation integrated hierarchy vs massless daughter on 96 bins); if it fails at 0.3%, examine whether deviations shrink with `ddm.momenta_bins = 192` (discretization) or not (physics bug).
  - The Cl^TE zero-crossing trap does not apply (TT/P(k) only) — do NOT add blind max-rel-diff checks on TE.
  - Tolerances may need mild retuning (documented in commit message) but relative deviations must SHRINK when momenta_bins doubles — that is the acceptance criterion for discretization-limited checks.

- [ ] **Step 4: Smoke-run the ini:** `./class test/dotsyntax_dcdm_wdm.ini` — completes without error.

- [ ] **Step 5: Commit**

```bash
git add scripts/validate_dcdm_wdm.py
git add -f test/dotsyntax_dcdm_wdm.ini
git commit -m "dcdm_wdm: physics-limit validation script + example ini"
```

---

### Task 6: Validation notebook (paper Fig. 3 recreation)

**Files:**
- Create: `notebooks/dcdm-wdm-massive-decay-products.ipynb`

**Interfaces:** consumes the built classy module and the validated feature.

- [ ] **Step 1: Write the notebook** with this structure (follow the style of `notebooks/370-axion-species.ipynb`: markdown intro with physics + equations, house plot style, dataviz-skill-conformant colors):
  1. **Intro (markdown):** the model (f, Γ, ε / v_kick), the paper reference, what CLASSpp does differently (per-bin ODE injection vs integral equation), parameter mapping: `f = Omega_ini_ddm/(Omega_ini_ddm + Omega_cdm)`, `Γ[km/s/Mpc] ↔ Γ⁻¹[Gyr]` conversion, internal conventions (q_kick=10, T_wdm=T_cmb).
  2. **Background sanity:** run one model (f=0.3, Γ⁻¹=1 Gyr, v=0.05); plot ρ_dcdm, ρ_wdm, w_wdm(a) from `cosmo.get_background()`; overlay the exact injection-integral for ρ_wdm (same formula as the C++ golden test) — sub-percent agreement.
  3. **Limits (plots):** vkick→0 vs ΛCDM; vkick=1 vs dcdm_dr; residual panels.
  4. **Fig. 3 recreation:** for f=0.1 and Γ⁻¹ ∈ {0.1, 1, 10} Gyr with v_kick at the paper's 1σ values read off their Fig. 2 left panel (approx: v/c ≈ 0.03 for 0.1 Gyr, 0.011 for 1 Gyr, 0.02 for 10 Gyr — state clearly these are read off the figure), three stacked panels: ΔP(k)/P, ΔC_ℓ^TT/C_ℓ, ΔC_L^φφ/C_L vs ΛCDM. Overlay k_fs = 2π/λ_fs and ℓ_fs = πχ*/λ_fs vertical dashed lines using paper Eq. (1) computed from the run's background. State the expected outcome from the paper: suppression onset at k_fs, C_L^φφ deviations reaching ~5–10%, low-ℓ TT changes small.
  5. **Convergence:** ΔP(k)/P for momenta_bins ∈ {48, 96, 192} at fixed model — curves must visually coincide (quote max spread).
  6. **Notes/limitations (markdown):** synchronous only, no tensor contribution, no fluid approximation (runtime scales with bins × l_max), q_min_ratio window (very short lifetimes unsupported), unequal-mass daughters not implemented (single-daughter symmetric case per the paper's main text).
- Conversion helper to include: `Gamma_km_s_Mpc = (1/Gyr_lifetime) * 977.79/ (h*100) ...` — derive properly in the notebook: Γ[km/s/Mpc] = (1 Gyr/Γ⁻¹) × 977.792 (since 1 km/s/Mpc = 1/977.792 Gyr⁻¹).
- Keep total notebook runtime manageable: `l_max_scalars=2500`, `momenta_bins=96` (48 for the scan panels if needed), note timings.

- [ ] **Step 2: Execute the notebook end-to-end** (`jupyter nbconvert --to notebook --execute --inplace notebooks/dcdm-wdm-massive-decay-products.ipynb` with a generous `--ExecutePreprocessor.timeout`), confirm all cells ran and figures render (spot-check the physics claims of §4).

- [ ] **Step 3: Commit**

```bash
git add notebooks/dcdm-wdm-massive-decay-products.ipynb
git commit -m "dcdm_wdm: validation notebook recreating arXiv:2606.14849 Fig. 3 physics panels"
```

---

## Self-review checklist (done at plan time)

- Spec coverage: §2.1 inputs → Task 1 ctor; §2.2 normalization → Task 1; §2.3 background/injection → Tasks 1-3; §2.4 hierarchy/ICs/stress-energy/scope guards → Task 4 (+gauge guard Task 3); §3 files/factory → Tasks 1,3; §5 tests → Tasks 1-5 + notebook Task 6. Fluid approximation deferred (documented).
- Known risks called out inline: background title-string names (Task 3 step 4), FillSources conventions, shooting sign conventions, tolerance retuning rule (Task 5).
- Type consistency: `FillInjection(a, rho_dcdm, J, dJdlnq)` used identically in Tasks 1-3; index accessors named `bg_f_index/bg_dfdlnq_index/bg_inj_index/bi_f_index/bi_dfdlnq_index` throughout; layout accessors `dcdm_layout/wdm_layout`.
