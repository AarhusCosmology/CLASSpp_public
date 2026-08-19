#pragma once
#include <memory>
#include <string>
#include <vector>

#include "background.h"
#include "common.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "decay_transition_kernel.h"
#include "dncdm_species.h"
#include "species/species_build_context.h"

class BackgroundModule;

/**
 * DNCDMProxySpecies: a cheap stand-in for the ν_H ↔ ν_l + φ composite
 * (`dr_representation = proxy`).
 *
 * WHY. The exact scheme (DNCDMInvSpecies) carries a resolved PSD for BOTH daughters
 * in the perturbations — the overwhelming majority of a run's variables — and applies
 * the linearised collision operator to every one of them at every k. At high Γ with a
 * light parent that is orders of magnitude too slow for parameter estimation. The
 * proxy keeps the physics the CMB actually sees and drops the rest:
 *
 *   BACKGROUND — kept essentially exact. The parent keeps its own q-grid; the two
 *   daughters keep PSDs on a COARSE grid, and the shipped DecayTransitionKernel drives
 *   all three. Nothing is re-derived, so number conservation, energy conservation and
 *   detailed balance still hold on the grid to machine precision — those identities
 *   are the kernel's, and are pinned in decay_kernel_test.cpp.
 *
 *   PERTURBATIONS — replaced by a relaxation-time (RTA) closure. The daughters become
 *   ONE integrated massless hierarchy each (DarkRadiationSpecies) instead of one
 *   hierarchy per momentum bin, and the linearised collision operator becomes two
 *   relaxation terms:
 *
 *     ℓ = 0, 1  exchange:  each species' contrast relaxes toward the sector's common
 *               contrast at the rate the collision is actually reprocessing that
 *               species. Exactly conservative in δρ and (ρ+P)θ by construction (the
 *               target is the λ-weighted mean), so the Einstein equations never see a
 *               source they should not.
 *
 *     ℓ ≥ 2     damping:   anisotropic stress and higher moments are damped at
 *               Γ_T,ℓ = a Γ⁰ (ρ_H/ρ_sec) [β_ℓ r₃ + α_ℓ r₅], with β_ℓ = ℓ(ℓ+1)/6,
 *               α_ℓ = (3ℓ⁴+2ℓ³−11ℓ²+6ℓ)/32 (arXiv:2203.09075 eq. 16) and r₃, r₅ set
 *               by `dr_rta_form`.
 *
 * Under the DEFAULT form the two exponents are not rivals: the γ⁻⁵ piece is what
 * survives AT detailed balance (arXiv:2011.01502 §5, re-derived in arXiv:2203.09075),
 * and the γ⁻³ piece is the Hannestad–Raffelt kinematic rate, which survives only in so
 * far as the sector is OUT of balance — which is why its amplitude carries ε_ne, the
 * background's own departure-from-balance order parameter (net/gross parent loss).
 * arXiv:2203.09075 takes the stronger position that the γ⁻³ rate is not an ab initio
 * result at all, so `dr_rta_form = copw` drops r₃ and keeps only that paper's γ⁻⁵
 * term. See RtaForm and TransportRate.
 *
 * Children order (the child_layouts contract, #358): kParent=0, kFermion=1, kBoson=2
 * — the same order as DNCDMInvSpecies, so tooling that walks the sector does not have
 * to branch on which representation built it.
 *
 * Synchronous gauge only, no tensors, no fluid approximation — the same guards as the
 * exact composite, for the same reasons. Calibration and the accuracy-vs-cost ladder:
 * docs/superpowers/specs/2026-08-14-dncdm-proxy-representation.md.
 */
class DNCDMProxySpecies : public CompositeSpecies {
 public:
  /** The parent keeps its momentum-resolved hierarchy; only the daughters collapse. */
  bool HasNcdm() const override {
    return true;
  }

  enum ChildIndex { kParent = 0, kFermion = 1, kBoson = 2 };

  /** Layout: parent gets the NCDM per-q layout, both daughters the DR layout. */
  static const NCDMBaseSpecies::PerturbLayout& parent_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kParent]);
  }
  static const DarkRadiationSpecies::PerturbLayout& dr_layout(const BaseSpecies::PerturbLayout& my,
                                                              int child) {
    return static_cast<const DarkRadiationSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[child]);
  }

  /** `q_d` / `dq_d` are the shared daughter grid (log trapezoid, half-weight
   *  endpoints), built by Create so the kernel can take a view of a buffer this
   *  object owns for its whole life. */
  DNCDMProxySpecies(std::unique_ptr<DNCDMSpecies> parent,
                    std::unique_ptr<DarkRadiationSpecies> fermion,
                    std::unique_ptr<DarkRadiationSpecies> boson,
                    std::vector<double> q_d,
                    std::vector<double> dq_d,
                    DecayTransitionKernel::Config cfg,
                    const background* pba,
                    const BackgroundModule* bgm);

  /** Build from a pre-built parent. Called by DNCDM_DR_Species::CreateAll when
   *  dr_representation = proxy. */
  static Named Create(std::unique_ptr<DNCDMSpecies> parent, const SpeciesBuildContext& ctx);

  DNCDMSpecies& parent() {
    return *parent_;
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  /** The collision's per-bin self-coupling — the background system's whole
   *  stiffness (K/eps -> a*Gamma once the parent is non-relativistic). Reported so
   *  the exponential evolver integrates it exactly. */
  void BackgroundDerivsDiagonal(double tau,
                                const double* y,
                                double* diag,
                                const double* pvecback) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

  /** The RTA's Jacobian diagonal, -(nu_i + Gamma_T,l) per state variable. Reported
   *  so the exponential evolver can integrate the relaxation exactly rather than
   *  resolving it with the step size. */
  void PerturbDerivsDiagonal(const BaseSpecies::PerturbLayout& layout,
                             double tau,
                             const double* y,
                             double* diag,
                             const perturb_parameters_and_workspace& ppaw) const override;

  /** k_output_values time series: the sector's total anisotropic stress
   *  a⁴Π_νφ ≡ a⁴·Σ_{i∈{H,l,φ}} (ρ̄_i+p̄_i)σ_i, plus its per-child decomposition.
   *
   *  Same column names and same construction as DNCDMInvSpecies::PrintVariables,
   *  so a proxy run and an exact run are readable by one analysis script — which
   *  is the only way the closure can be checked on the quantity it approximates.
   *  Composite-owned for the same reason: the children carry no valid
   *  collection_index_, so their layouts are reached through THIS composite's
   *  nested layout. */
  void PrintVariables(PerturbColumnWriter& writer,
                      const BaseSpecies::PerturbLayout* base,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  /** `dr_rta_form`: which transport rate the ℓ ≥ 2 damping uses.
   *
   *  kPowers (default) is the fitted two-power form Γ_T/Γ = C₃ ε_ne^n₃ γ⁻³ + C₅ γ⁻⁵
   *  with a non-relativistic shut-off. It is a CALIBRATION of this code's own runs,
   *  accurate to a few per cent across the measured window, and it is what the
   *  shipped defaults are tuned for.
   *
   *  kCOPW is the analytic prescription of Chen, Oldengott, Pierobon & Wong
   *  (arXiv:2203.09075) — eq. (13) with 𝓕 from eq. (14). A SINGLE γ⁻⁵ term carrying
   *  the paper's own (1/12) amplitude: no γ⁻³ piece (that paper argues against one)
   *  and no fitted normalisation, so `dr_rta_C3` defaults to 0 and `dr_rta_C5` to 1
   *  under this form, and `dr_rta_n3` / `dr_rta_vshut` are ignored outright.
   *
   *  It is the one form derived from first principles rather than fitted, which is
   *  why it is worth being able to run, but it is derived AT detailed balance and
   *  does not reproduce the measured γ-dependence inside the converged window; see
   *  TransportRate. */
  enum class RtaForm { kPowers, kCOPW };

  RtaForm rta_form() const {
    return rta_form_;
  }

  /** The transport rate Γ_T,ℓ (conformal, 1/Mpc) at the current background row. */
  double TransportRate(int l, const double* pvecback) const;

  /** Where the two derived quantities TransportRate depends on sit in pvecback:
   *  the sector's total energy density and the departure-from-balance parameter. */
  int bg_rho_sec_index() const {
    return index_bg_rho_sec_;
  }
  int bg_eps_ne_index() const {
    return index_bg_eps_ne_;
  }

  /** Same, with the scale factor passed in rather than read out of pvecback through
   *  the BackgroundModule. That indirection is the only thing in this function that
   *  needs a live module, so splitting it here is what lets the unit test check the
   *  rate against arXiv:2203.09075 eq. (13) on a hand-filled background row. */
  double TransportRate(int l, double a, const double* pvecback) const;

  /** 𝓕(x) of arXiv:2203.09075 eq. (14) — ½e⁻ˣ[−1 + x − eˣ(x²−2)Γ(0,x)] — via that
   *  paper's own two-branch approximation, which avoids the incomplete gamma
   *  function. Past x ~ 40 the parent is so non-relativistic that inverse decays are
   *  shut and the whole term is irrelevant.
   *
   *  𝓕 ~ −ln x as x → 0, which is physical and must not be floored: the rate carries
   *  x⁵𝓕(x), so the product still vanishes. Public for the unit test, which holds
   *  both branches to the paper's error bars against an independent Γ(0,x). */
  static double CurlyF(double x);

  /** α_ℓ of arXiv:2203.09075 eq. (16), the ℓ-dependence of the γ⁻⁵ term.
   *  α_0 = α_1 = 0 (the ℓ ≤ 1 collision integrals vanish identically by
   *  energy/momentum conservation), α_2 = 1, and it grows as ℓ⁴ — small angular
   *  scales are wiped out fastest at a fixed decay opening angle. */
  static constexpr double AlphaL(int l) {
    const double x = static_cast<double>(l);
    return (3. * x * x * x * x + 2. * x * x * x - 11. * x * x + 6. * x) / 32.;
  }

  /** Stored -> bare occupation for the parent leg (#385): the kernel's band factor
   *  and (1±f) coefficients need bare per-dof occupation, the parent's slot holds
   *  deg_H-suppressed stored occupation. g_H = 2 is the parent's spin dof.
   *  Reads GetDeg() LIVE — shooting updates it between iterations. */
  double KappaStoredToBare() const {
    return parent_->GetDeg() * 8. * _PI_ * _PI_ * _PI_ / 2.;
  }

  /** Density normalisation for a BARE-occupation grid moment: rho = BareFactor()
   *  * a^-4 * sum dq q^2 eps f. Equals DrPsdSpecies' g = 2 fermion factor by
   *  construction; the boson takes half of it (the kernel folds the
   *  g_H·g_l/g_phi = 2 ratio into its phi deposit, so the same f_phi with half the
   *  factor is the physical energy density). */
  double BareFactor() const {
    return parent_->factor() / KappaStoredToBare();
  }

  /** Comoving number of a daughter, in the kernel's bare grid units (no factor):
   *  sum dq q^2 f, halved for the boson. */
  double DaughterNumber(const double* f, bool boson) const;
  /** Physical energy density of a daughter from its PSD. */
  double DaughterRho(const double* f, double a, bool boson) const;
  /** n = BareFactorNumber(a) * sum dq q^2 f: divides a published number column back
   *  into the bare grid moment the gross rate is expressed in. */
  double BareFactorNumber(double a) const;

 private:
  /** Gather the three background PSDs into bare-occupation buffers, run the
   *  kernel, and hand back the per-species number/energy sources. */
  void ApplyKernelBackgroundDerivs(double a, const double* y, double* dy, double a_prime_over_a);

  /** Fill the derived background quantities (gross decay rate, ε_ne, sector energy)
   *  into pvecback. Called from ComputeBackground, so the perturbation side reads a
   *  row consistent with the state it was computed from.
   *
   *  NOT optional output: ε_ne is the departure-from-balance order parameter that
   *  multiplies the γ⁻³ term, and ρ_sec sets the parent's energy fraction, so
   *  TransportRate reads both on every call. */
  void ComputeDerivedBackground(double a, const double* pvecback_B, double* pvecback);

  /** Ceiling on every relaxation rate, as a multiple of max(k, a'/a).
   *
   *  nu_H -> a*Gamma once the parent is non-relativistic, which at Gamma = 1e11 is
   *  ~10^6 times the expansion rate for the ~90% of the run in which the parent is
   *  already extinct and contributes nothing to any observable. Left uncapped that
   *  is a stiff mode the explicit evolver has to resolve step by step. Capping at a
   *  large multiple of the fastest rate the mode actually resolves leaves the
   *  slaved solution unchanged — the same trade tight coupling makes for the
   *  photon-baryon fluid — and keeps the system non-stiff. Raise it (and watch the
   *  step count) if a cell is ever suspected of being cap-limited. */
  static constexpr double kRateCapFactor = 20.;

  /** Default `dr_rate_cap`: the background collision's per-bin rate K/eps is
   *  limited to this multiple of a'/a.
   *
   *  Scanned at Γ = 10⁷, m = 0.3, dr_N_q = 96 (etd, background only). The cap is a
   *  no-op until it is well below the rate it is capping:
   *
   *      cap      0       1e5      1e4      1e3      1e2
   *      steps  67296    62337    26659    13710     5170
   *      E_com  1.254042 1.254042 1.254042 1.254042 1.254118
   *      |ΔH/H|    —     3.6e-15  3.6e-15  8.9e-13  9.9e-09
   *
   *  1e3 buys 5× fewer steps for 1e-12 in H, and 1e2 is still only 1e-8 — so this
   *  default has a decade of margin. It is a RELATIVE cap, so re-scan it if a cell
   *  at much higher Γ matters: the uncapped ratio K/(eps·a'/a) grows with Γ. */
  static constexpr double kDefaultRateCap = 1e3;

  /** Below this fraction of the sector's energy the parent is dropped from the
   *  coupling entirely. Its exchange rate is a*Gamma once it is non-relativistic —
   *  the stiffest eigenvalue in the system — while its weight in every observable
   *  has underflowed to ~1e-60. */
  static constexpr double kParentFloor = 1e-12;

  /** Stack-buffer bound for the per-multipole scratch in AddCouplingDerivs. The
   *  buffers are per-call because the k-loop is threaded; a member would be a data
   *  race. l_max_ncdm/l_max_dr never approach this. */
  static constexpr int kLMaxScratch = 128;

  // Non-owning: children_ owns them (same pattern as DNCDMInvSpecies).
  DNCDMSpecies* parent_          = nullptr;
  DarkRadiationSpecies* fermion_ = nullptr;
  DarkRadiationSpecies* boson_   = nullptr;

  const background* pba_       = nullptr;
  const BackgroundModule* bgm_ = nullptr;

  // Daughter momentum grid (shared by both daughters — they are produced in pairs
  // at q2 + q3 = eps1, so a common grid is what makes the two-bin deposit exact in
  // energy as well as number). DECLARED BEFORE kernel_: the kernel holds a
  // non-owning GridView into these buffers and is built from them.
  std::vector<double> q_d_, dq_d_;
  std::unique_ptr<DecayTransitionKernel> kernel_;

  double f_ini_l_ = 1.0, f_ini_phi_ = 0.0;

  /** Fitted transport-rate coefficients (`dr_rta_C3` / `dr_rta_C5` / `dr_rta_n3` /
   *  `dr_rta_vshut`). These are a CALIBRATION, not a derivation -- see TransportRate
   *  -- so they stay inputs: re-fit them if the measurement improves. n3_ and vshut_
   *  are exponents of the fitted form only; kCOPW ignores them. */
  double C3_ = 0.2802, C5_ = 2.5228, n3_ = 0.8714, vshut_ = 1.7995;

  RtaForm rta_form_ = RtaForm::kPowers;

  // Composite-owned background integration state: the two daughter PSDs.
  int index_bi_f_l_   = -1;
  int index_bi_f_phi_ = -1;
  // Published background quantities.
  int index_bg_rate_    = -1;  // gross comoving decay rate (conformal, 1/Mpc)
  int index_bg_eps_ne_  = -1;  // net/gross: 1 = free decay, 0 = detailed balance
  int index_bg_rho_sec_ = -1;  // total sector energy density
  int index_bg_nu_H_    = -1;  // exchange rates, uncapped (the cap needs k)
  int index_bg_nu_l_    = -1;
  int index_bg_nu_phi_  = -1;
  int index_bg_n_l_     = -1;  // daughter number densities, for the output columns
  int index_bg_n_phi_   = -1;

  // Hot-path scratch (BackgroundDerivs is called by the evolver, so no allocation).
  std::vector<double> fH_bare_, df_H_, df_l_, df_phi_;
  std::vector<double> diag_H_, diag_l_, diag_phi_;
  /** The background row of the call in flight, for the extinction cut. Set and
   *  cleared by BackgroundDerivs; the background integration is single-threaded. */
  const double* pvecback_row_ = nullptr;
  /** Last eps_ne computed on a STORED row. Trial rows publish this rather than pay
   *  a kernel pass for a number only the stored table is ever read for. */
  double eps_ne_cached_ = 1.0;
  /** State the kernel's cached transitions were prepared at, so the diagonal can
   *  skip a redundant pass when it provably matches. NaN = nothing cached. */
  double prep_a_ = -1., prep_sum_ = 0.;
};
