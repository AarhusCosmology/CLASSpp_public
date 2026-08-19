#pragma once
#include <algorithm>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

#include "../species/ncdm_base_species.h"
#include "../species/species_build_context.h"
#include "background.h"
#include "decay_transition_kernel.h"  // Statistics

class BackgroundModule;

/**
 * Massless momentum-resolved dark-radiation daughter of the ν_H ↔ ν_l + φ system
 * (arXiv:2011.01502; design §5 "dr_psd_species"). ONE class plays both daughters —
 * statistics = Fermion (Pauli) or Boson (Bose) is a construction parameter, so there
 * is exactly one daughter equation to get right.
 *
 * The comoving PSD f(q,τ) is integrated directly per momentum bin (positivity floor
 * kFFloor, not ln f — gains dominate at low q). ∂f̄/∂lnq is splined live from f̄ vs
 * ln q every ComputeBackground, from f̄ ITSELF and not ln f̄: a decay-filled PSD is a
 * bump on the floor plateau, and the ln of that cliff rings.
 *
 * Perturbations are the UNNORMALIZED F_ℓ(q) ≡ f̄·Ψ_ℓ = δf hierarchy — deliberately NOT
 * the normalized Ψ = δf/f̄ that DNCDMSpecies and the rest of the NCDM family use. A
 * daughter is sourced from EXACTLY zero abundance, where f̄̇/f̄ is smooth but enormous
 * at early times; both the −(f̄̇/f̄)Ψ dilution and the collision operator's 1/f̄ divide
 * by it, so the normalized form diverges. In F-space those terms cancel structurally
 * (derivation in dr_psd_species.cpp's perturbation section). The free-streaming
 * ℓ,ℓ±1 coefficients are unchanged (same form as DNCDMSpecies with M_ = 0, ε = q);
 * the metric driver uses ∂f̄/∂lnq instead of ∂lnf̄/∂lnq and the moment weights drop
 * their f̄. Synchronous gauge only; no fluid approximation.
 *
 * Standalone it is a free-streaming massless PSD species ("all species are equal"):
 * with no source (BackgroundDerivs is a no-op) comoving f is constant and ρ redshifts
 * as radiation ∝ a⁻⁴. The composite DNCDMInvSpecies injects the kernel source into the
 * f-slots; the sector's density is reserved by the composite, so a daughter's own
 * GetOmega0() is 0 (decay-product neutrality, mirrors WdmDecayProductSpecies).
 *
 * TWO momentum grids, coupled by EXACT INTEGER SUBSAMPLING (`dr_bg_refine` = m): the
 * background grid must be fine to resolve the decay-injection band, but the Boltzmann
 * hierarchy costs N_pt·(l_max+1) variables per k-mode, so the perturbation grid must
 * be allowed to be coarser. q_[j] is read straight out of q_bg_[j·m] — never
 * interpolated, never re-derived from a fresh log-spaced formula. Both reasons are
 * load-bearing: (1) splining background quantities across the daughter's injection
 * cliff rings, and exact subsampling removes that by construction, because a PT-grid
 * point IS a BG-grid array entry; (2) every background column is readable at PT index
 * j as BG column j·m with no reconstruction. Points span N−1 intervals, so
 * N_bg−1 = m·(N_pt−1), NOT N_bg = m·N_pt. m = 1 (the default) reproduces the
 * single-grid case bit-for-bit.
 */
class DrPsdSpecies : public NCDMBaseSpecies {
 public:
  static constexpr std::string_view kTypeName = "dr_psd";
  static constexpr double kFFloor = 1e-100;  // positivity floor on f (gains dominate at low q)

  /** kinematic_q_max: the highest momentum this daughter will ever be FED, from the
   *  parent's decay kinematics (DecayTransitionKernel::MaxDaughterMomentum at a = 1).
   *  Used ONLY to choose the default dr_q_max, so a composite-built daughter is sized
   *  to cover its own emission band instead of inheriting a fixed 1e2 that silently
   *  clamps -- and thereby leaks energy -- as soon as a*m/2 outgrows it. 0 (the
   *  default) means "no parent, no bound known", which is the standalone dr_psd case:
   *  a daughter with no decay channel feeding it has no kinematic requirement, and
   *  the historical default stands. An explicit dr_q_max ALWAYS wins over this; the
   *  composite warns once if the user's choice sits below the band (see
   *  DNCDMInvSpecies::Create). */
  DrPsdSpecies(FileContent* pfc,
               const std::string& instance_name,
               const NcdmSettings& settings,
               const background* pba,
               const BackgroundModule* bgm,
               Statistics stat,
               double kinematic_q_max = 0.);

  /** Margin on the kinematic bound when defaulting dr_q_max. The Gauss-Legendre
   *  nodes approach the band edge without reaching it (measured: 0.966-1.000 of the
   *  bound across the relativistic and non-relativistic regimes), so a few percent
   *  is enough to guarantee the top node stays bracketed. */
  static constexpr double kQMaxMargin = 1.05;

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);  // standalone dr_psd

  // ── Index / grid accessors for the composite (Phase 3/4) ───────────────────
  Statistics statistics() const {
    return stat_;
  }
  /** PT-grid index -> BG-grid index. The perturbation grid is always an exact
   *  integer SUBSAMPLE of the background grid (never a spline: background
   *  quantities ring across the daughter's ~87-decade injection cliff, and never an
   *  independently spaced grid: two exp() formulas agree only to rounding). Which
   *  subsample differs by sampling, and the difference is exactly one index:
   *
   *    log    q_bg[i] = q_min*exp(i*h), so index 0 IS the endpoint q_min and
   *           q_bg[j*m] is uniform in ln q from that endpoint  ->  j*m
   *    linear q_bg[i] = (i+1)*h with q=0 an IMPLICIT zero-weight point, so index 0
   *           sits at h, not at 0. q_bg[j*m] would leave a first cell of width h
   *           among cells of width m*h -- not a trapezoid rule at all. Subsampling
   *           from the other end, q_bg[(j+1)*m - 1] = (j+1)*(m*h), is the canonical
   *           qm_trapz grid of step m*h  ->  (j+1)*m - 1
   *
   *  Hence the divisibility requirement also differs: N_bg-1 divisible by refine for
   *  log (N points span N-1 intervals from the endpoint), N_bg divisible by refine
   *  for linear (N cells of width h from the implicit zero). Both keep the last PT
   *  point ON the last BG point, so the exact q_max endpoint is inherited. */
  int bg_index(int iq) const {
    return (q_sampling_ == QSampling::Linear) ? (iq + 1) * refine_ - 1 : iq * refine_;
  }

  int bi_f_index() const {
    return index_bi_f_;
  }
  int bg_f_index() const {
    return index_bg_f_;
  }
  int bg_dfdlnq_index() const {
    return index_bg_dfdlnq_;
  }
  const std::vector<double>& q_bg() const {
    return q_bg_;
  }
  const std::vector<double>& dq_bg() const {
    return dq_bg_;
  }
  const std::vector<double>& dq() const {
    return dq_;
  }
  /** Background-to-perturbation grid refinement m: PT bin j is BG bin j·m. Any read
   *  of a per-bin BACKGROUND column at a PT loop index must be scaled by this. */
  int refine() const {
    return refine_;
  }

  /** Pin this daughter's T to the parent's (#385 kernel-unit consistency: the
   *  kernel boundary conversion assumes every species reports occupation at the
   *  SAME dimensionless q = p/T). Re-runs SetDegAndFactor(GetDeg()) so factor_
   *  picks up the new T_. Called by the DNCDMInvSpecies constructor for both
   *  daughters, before anything reads factor_. */
  void PinTemperature(double T);

  // ── Background ─────────────────────────────────────────────────────────────
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
  /** Comoving-number density n (already accumulated by ComputeMomenta in
   *  ComputeBackground). Published so the exact transition-network identity
   *  n_H + ½(n_φ + n_l) = const (arXiv:2011.01502 fig. 11) can be read straight
   *  off background.dat for all three species, as the reference implementation
   *  does — instead of being reconstructed from the per-bin f columns. */
  double Number(const double* pvecback) const {
    return pvecback[index_bg_number_];
  }
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override {
    return a * H * (pvecback[index_bg_pseudo_p_] - 5. * pvecback[index_bg_p_]);
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Decay-product neutrality overrides (mirror wdm_decay_product.h:148-167) ──
  // The daughter's sector density is owned by the composite; standalone it is a
  // massless, always-relativistic free-streaming bath.
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
    return true;
  }
  double BackgroundAIni(double a_proposed, double /*tol*/) const override {
    return a_proposed;  // massless: always relativistic, no early-start constraint
  }
  void CheckUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const override {}
  bool IsUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const override {
    return true;
  }

  // ── Reduced-moment perturbation representation (dr_reduced_moments) ─────────
  /** Carry n_moments q-MOMENTS of delta-f per multipole instead of the full q-grid.
   *
   *  WHY IT IS A DROP-IN. For a MASSLESS daughter eps = q, so the free-streaming
   *  coefficient qk/eps is exactly k -- it carries no q at all. Every q-moment of the
   *  hierarchy therefore obeys the SAME hierarchy, with no closure error (design note
   *  M3, a theorem rather than a fit). And all four stress-energy sums carry the same
   *  weight q^3 once eps = q, so the metric only ever reads the ENERGY moment at
   *  l = 0, 1, 2. Replacing the q index by a moment index thus changes three things and
   *  nothing else: how many slots there are, what the metric driver is (the projected
   *  d f-bar / d ln q rather than its pointwise value), and how the moments are read
   *  back out.
   *
   *  THE BASIS IS FROZEN, and that is not an optimisation. psi_j is built from
   *  w = f-bar(1 -/+ f-bar), which evolves; a basis that tracked `a` would make the
   *  state variables mean something different at every step, and dm_j/dtau would
   *  acquire an integral of (dpsi_j/dtau) delta-f that this hierarchy does not have.
   *  Measured cost of freezing at Gamma=1e7 (l=2 collision rate, scale-relative p90):
   *  at n_moments = 2 it is a factor 2.5, at 4 it is 1.5, at 6-8 it is ~1.2 -- which is
   *  why 6 is the recommended working point rather than the 2 the design proposed.
   *
   *  `psi` is n_moments * q_size, row-major in the moment index; `gnorm` is the
   *  (diagonal) Gram. Both come from ReducedCollisionOperator so there is exactly one
   *  Stieltjes recurrence in the build. */
  void SetReducedBasis(int n_moments,
                       const std::vector<double>& psi,
                       const std::vector<double>& gnorm,
                       double q_ref,
                       double alpha0);
  int reduced_moments() const {
    return reduced_moments_;
  }
  bool reduced() const {
    return reduced_moments_ > 0;
  }
  /** The energy moment at one multipole, read straight out of the state vector.
   *  Defined here because StressEnergy and FillSources both need it and a second
   *  spelling of the (m_1 + alpha0 m_0) combination is exactly the kind of thing that
   *  drifts. */
  double ReducedEnergyAtL(const NCDMBaseSpecies::PerturbLayout& layout,
                          const double* y,
                          int l) const {
    if (reduced_moments_ < 2 || l > layout.l_max)
      return 0.;
    const double m0 = y[layout.index_per_q[0] + l];
    const double m1 = y[layout.index_per_q[1] + l];
    return red_q_ref_ * (m1 + red_alpha0_ * m0);
  }
  /** Projection of a q-resolved grid function onto the frozen moments. */
  void ReducedProject(const double* fq, double* m) const;
  /** The metric driver in reduced coordinates: D_j = Sum dq q^2 psi_j (df-bar/dlnq). */
  void ReducedDriver(const double* pvecback, double* D) const;
  /** Cap on n_moments, so the hot paths can use a stack array. The measured accuracy
   *  ceiling of the reduction is 9 (reduced_operator_test), so 12 is headroom, not a
   *  limit anyone will meet. */
  static constexpr int kMaxReducedMoments = 12;

  // ── Perturbations ──────────────────────────────────────────────────────────
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

  /** Tensor modes: no slots (daughter tensor anisotropic stress neglected;
   *  documented scope limit, empty until late times — mirrors dcdm_wdm). */
  void RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& /*layout*/,
                                         perturb_vector* /*pv*/,
                                         const precision* /*ppr*/,
                                         int& /*index_pt*/,
                                         const perturb_workspace* /*ppw*/,
                                         int /*gauge*/) override {}
  /** Synchronous gauge only (guarded by the composite / standalone factory). */
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& /*layout*/,
                                     double* /*y*/,
                                     const PerturbIcContext& /*ctx*/) override {}

  /** name() suffixed by statistics ("_l"/"_phi", the paper's ν_l/φ notation).
   *  The composite constructs both daughters under the SAME instance name (they
   *  are not separately dot-addressable), so background/output COLUMN NAMES
   *  must disambiguate here or the fermion and boson columns collide.
   *  Public because DNCDMInvSpecies names its children's per-q diagnostic columns
   *  and must use this one definition rather than re-spelling the suffix. */
  std::string ColumnLabel() const {
    return name() + (stat_ == Statistics::Fermion ? "_l" : "_phi");
  }

 protected:
  /** The DNCDM seam. dr_psd publishes ∂f̄/∂lnq (its F-hierarchy driver), so the
   *  base class's NORMALIZED contract is recovered by dividing by f̄. Reached only
   *  from NCDMBaseSpecies' tensor hierarchy, for which dr_psd registers no slots —
   *  the scalar path reads the ∂f̄/∂lnq column directly and never normalizes. */
  double GetDlnf0Dlnq(int iq, const double* pvecback) const override {
    const int ibg = bg_index(iq);  // PT loop index -> BG column (exact subsample)
    return pvecback[index_bg_dfdlnq_ + ibg] / std::max(pvecback[index_bg_f_ + ibg], kFFloor);
  }

 private:
  /** Parse the statistics dot-key into the enum; class_test_severe on anything
   *  other than "fermion"/"boson". Used by the standalone CreateAll. */
  static Statistics ParseStatistics(FileContent* pfc, const std::string& instance_name);

  /** Fill q_bg_/dq_bg_ (log-trapezoid, half-weight endpoints, thesis §3), then the
   *  perturbation grid q_/dq_ as the exact refine_-subsample of it. */
  void BuildLogTrapezoidGrids();

  /** Fill q_bg_/dq_bg_ with the UNIFORM trapezoid grid of quadrature.c's `qm_trapz`
   *  (dr_q_sampling = linear):
   *
   *      h = q_max/N,  q_i = (i+1)h  (i = 0 … N−1),  dq_i = h,  last point half-weight
   *
   *  i.e. a trapezoid over [0, q_max] carrying q=0 as an implicit zero-weight point
   *  (q²f vanishes there), so the FIRST grid point takes full weight and only the
   *  last is halved. This is byte-identical to `get_qsampling_manual(..., qm_trapz)`
   *  with qmin==0 in both reference implementations, and to the published
   *  `np.linspace(0, q_max, N+1)[1:]` grid used for the arXiv:2011.01502 background
   *  figures — which is the point: it is what makes a like-for-like comparison
   *  against those codes possible, and it is the grid the neutrino-lasing
   *  calculation is historically done on (thousands of uniform bins).
   *
   *  `dr_q_min` is IGNORED in this mode: the grid is anchored at zero by the
   *  convention above, and its first point is q_max/N by construction. */
  void BuildLinearTrapezoidGrids();

  /** Analytic thermal seed f0(q) = 1/(exp(q − ksi) ∓ 1) (− fermion, + boson),
   *  used only when an initial abundance is requested (else the daughter starts
   *  empty at kFFloor). Regular for ksi < q_min in the boson case. */
  double ThermalF0(double q) const;

  // Frozen reduced-moment basis; empty unless SetReducedBasis was called. See the
  // public block above for why it is frozen and what the coordinates mean.
  int reduced_moments_ = 0;
  std::vector<double> red_psi_;    // reduced_moments_ * q_size, row-major in the moment
  std::vector<double> red_gnorm_;  // diagonal Gram
  double red_q_ref_  = 1.;
  double red_alpha0_ = 0.;

  const background* pba_ = nullptr;
  Statistics stat_       = Statistics::Fermion;
  /** Initial occupation amplitude A: f(q, a_ini) = A · f_eq(q). 0 = starts empty.
   *  Set PER DAUGHTER from `dr_f_ini_l` (fermion ν_l) / `dr_f_ini_phi` (boson φ);
   *  falls back to 1 for both when a legacy thermal key (Omega_ini / dr_T / dr_ksi)
   *  is present, else 0. The published background figures use A = [0, 1, 1]
   *  (φ empty, ν_l full Fermi-Dirac), i.e. dr_f_ini_phi=0, dr_f_ini_l=1. */
  double f_ini_ = 0.;

  /** dr_q_sampling: momentum-grid family. Log is the repo default (wide dynamic
   *  range, needed once daughters are injected far from the thermal peak); Linear
   *  reproduces the reference codes' `qm_trapz` grid for cross-code comparison. */
  enum class QSampling { Log, Linear };
  QSampling q_sampling_ = QSampling::Log;

  // Log-trapezoid grid parameters. Only the BACKGROUND grid has free endpoints;
  // the perturbation grid inherits them exactly (q_[0]=q_bg_[0], q_.back()=q_bg_
  // .back()) and differs only in step: h_pt = refine_·h_bg.
  // (Linear sampling ignores q_min_bg_ — see BuildLinearTrapezoidGrids.)
  double q_min_bg_ = 1e-3, q_max_bg_ = 1e2;
  int N_q_bg_ = 100;
  int refine_ = 1;    // dr_bg_refine: N_bg−1 = refine_·(N_pt−1); 1 = one grid
  int N_q_    = 100;  // derived in BuildLogTrapezoidGrids, never read from input

  std::vector<double> dq_;     // pt-grid cell widths (half-weight endpoints, step h_pt)
  std::vector<double> dq_bg_;  // bg-grid cell widths

  int index_bg_number_   = -1;
  int index_bg_pseudo_p_ = -1;
  int index_bg_f_        = -1;
  int index_bg_dfdlnq_   = -1;
  int index_bi_f_        = -1;
  // Running severity of the negative occupations ComputeBackground sees on STORED
  // rows, i.e. before kFFloor hides them; they drive the one-per-decade warning
  // there and nothing else. Mutable because ComputeBackground is logically const on
  // the species; same pattern and same reasoning as DNCDMSpecies' negative_f_rows_.
  mutable double worst_negative_rel_     = 0.;  // |f|/peak of the worst excursion
  mutable double announced_negative_rel_ = 0.;

 public:
};
