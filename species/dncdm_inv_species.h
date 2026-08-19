#pragma once
#include <memory>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "decay_transition_kernel.h"
#include "dncdm_species.h"
#include "dr_psd_species.h"
#include "reduced_collision_operator.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"

class BackgroundModule;

/**
 * DNCDMInvSpecies: the ν_H ↔ ν_l + φ composite (arXiv:2011.01502; design §5).
 * Owns a collision-owned parent DNCDMSpecies (ν_H), a fermion DrPsdSpecies (ν_l),
 * a boson DrPsdSpecies (φ), and DecayTransitionKernels that derive all three
 * background derivatives (and the per-ℓ collision operator) from ONE discrete
 * transition network each, so number/energy conservation and detailed balance hold
 * on the grid to machine precision. TWO kernel instances, one per daughter-grid
 * resolution: the background network on the daughters' fine grid, the perturbation
 * network on their coarse dr_bg_refine-subsample. Only the BACKGROUND kernel is a
 * member; the perturbation ones live per-k in Scratch (thread safety).
 * Both read the daughters' background occupations at FULL resolution; only the
 * deposit/hierarchy grid is coarsened (DecayTransitionKernel's ctor doc says why
 * that distinction is not optional).
 *
 * All three species carry the unnormalized F = δf as their perturbation variable
 * (#386). The parent used to carry the normalized Ψ_H, which divides by an f̄_H
 * that decays without bound and made high-Γ runs unintegrable.
 *
 * The parent keeps ownership of its PSD (grid, ComputeMomenta, perturbation
 * hierarchy); the composite only supplies its background RHS by writing the kernel
 * source into the parent's f-slots in BackgroundDerivs — the same seam
 * DCDM_WDM_Species uses. With inverse decays OFF the factory builds DNCDM_DR_Species
 * instead, so the decay-only path is untouched.
 *
 * Children order (the child_layouts contract, #358): kParent=0, kFermion=1,
 * kBoson=2. Synchronous gauge only; no tensors, no fluid approximation for the
 * coupled trio (all guarded at Create). The parent must be a single momentum grid
 * (q_size == q_size_bg), also guarded at Create.
 */
class DNCDMInvSpecies : public CompositeSpecies {
 public:
  /** Parent and both daughters are momentum-resolved NCDM hierarchies. */
  bool HasNcdm() const override {
    return true;
  }

  enum ChildIndex { kParent = 0, kFermion = 1, kBoson = 2 };  // children_ order

  static const NCDMBaseSpecies::PerturbLayout& parent_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kParent]);
  }
  static const NCDMBaseSpecies::PerturbLayout& fermion_layout(
      const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kFermion]);
  }
  static const NCDMBaseSpecies::PerturbLayout& boson_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kBoson]);
  }

  DNCDMInvSpecies(std::unique_ptr<DNCDMSpecies> parent,
                  std::unique_ptr<DrPsdSpecies> fermion,
                  std::unique_ptr<DrPsdSpecies> boson,
                  DecayTransitionKernel::Config cfg,
                  const background* pba,
                  const BackgroundModule* bgm);

  /** Build the composite from a pre-built parent (from DNCDMSpecies::CreateAll).
   *  Reads inverse_decays / quantum_statistics / dr_* keys off the parent's
   *  instance; applies the input-time guards; constructs the two daughters and
   *  the kernel. Called by DNCDM_DR_Species::CreateAll when inverse_decays=yes. */
  static Named Create(std::unique_ptr<DNCDMSpecies> parent, const SpeciesBuildContext& ctx);

  DNCDMSpecies& parent() {
    return *parent_;
  }
  DrPsdSpecies& fermion() {
    return *fermion_;
  }
  DrPsdSpecies& boson() {
    return *boson_;
  }

  // (2*pi)^3, as a compile-time constant: KappaStoredToBare is called every
  // background/perturbation RHS evaluation (hot path), so the folding must not
  // depend on whether the compiler treats pow(2*_PI_, 3) as constant-foldable.
  static constexpr double kTwoPiCubed = 8. * _PI_ * _PI_ * _PI_;

  /** Kernel-boundary conversion (#385, occupation-suppression convention):
   *  bare-occupation = kappa * stored-occupation (parent leg only; the daughters
   *  already store bare occupation). g_H = 2 (design's parent spin dof). Reads
   *  GetDeg() LIVE at every call — shooting updates deg_H between iterations, so
   *  this must never be cached across a background/perturbation RHS. Public: the
   *  tests use it (dncdm_inv_test.cpp's perturbation-conservation reconstruction). */
  double KappaStoredToBare() const {
    return parent_->GetDeg() * kTwoPiCubed / 2.;
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  /** The collision's per-bin self-coupling, for the exponential evolver. This is the
   *  stiff part of the whole background system: diag_H -> -a*Gamma, which is what
   *  sets rkdp45's step (h*lambda = 3.0 against a stability radius of 3.3) for the
   *  ~90% of the run that happens after the parent is extinct. */
  void BackgroundDerivsDiagonal(double tau,
                                const double* y,
                                double* diag,
                                const double* pvecback) override;

  /** Writes the kernel background source into the parent+daughter f-slots of dy.
   *  Split out of BackgroundDerivs so it can be driven with an explicit a in the
   *  conservation unit test (which has no BackgroundModule to read a from). */
  void ApplyKernelBackgroundDerivs(double a,
                                   const double* y,
                                   double* dy,
                                   double a_prime_over_a = 0.);

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Closure + shooter hooks (combined 3-child reserve; mirrors DNCDM_DR) ────
  double GetOmega0() const override;
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

  // ── Perturbations: the per-ℓ collision coupling (design §4.3) ───────────────
  /** Runs after each child's free-streaming PerturbDerivs (CompositeSpecies two-
   *  phase contract). Adds the kernel per-ℓ collision operator (ℓ ≤ l_max_dncdm_col).
   *  All THREE species evolve the unnormalized F = δf (#386), so the operator output
   *  is used verbatim and no −(f̄̇/f̄)Ψ dilution term appears anywhere: in F-space it
   *  cancels structurally. Synchronous, scalar only. */
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

  /** Per-k working memory (BaseSpecies::PerturbScratch). Holds this composite's
   *  own perturbation kernel and the contiguous F/dF gather buffers.
   *
   *  The kernel is here, not shared, because DecayTransitionKernel is stateful by
   *  design: PrepareTransitions caches the whole transition geometry for one
   *  (a, f-bar) and ApplyPerturbationOperator advances a rolling per-node Legendre
   *  pair, which is what makes the caller's l-loop O(L) instead of O(L^2). All of
   *  that is per (k, tau). Two threads sharing one kernel corrupt each other's
   *  recurrence — observed as the strict-ascending-l contract firing with
   *  "expected l=3, got l=0", which is the assertion doing its job.
   *
   *  A kernel costs a few small vectors and one Gauss-Legendre setup to build, and
   *  it holds only non-owning GridViews into the species' grids, so one per k-mode
   *  is nothing next to a k-mode integration. */
  struct Scratch : BaseSpecies::PerturbScratch {
    std::unique_ptr<DecayTransitionKernel> kernel_pt;
    std::vector<double> fH_gather;
    std::vector<double> F_H, F_l, F_phi;
    std::vector<double> dF_H, dF_l, dF_phi;
    /** Collision-only dy, for the transport-rate diagnostic in PrintVariables.
     *  Sized lazily to pv->pt_size (the layout is not known at CreatePerturbScratch
     *  time) and used ONLY at k_output_values sampling points, never on the RHS. */
    std::vector<double> dy_coll;
    /** Per-bin Jacobian diagonals for the perturbation leg. Separate from the
     *  members BackgroundDerivsDiagonal uses: the background is single-threaded,
     *  the perturbations are not, and one k-mode must not see another's. Sized on
     *  the PERTURBATION grids, which is what kernel_pt's CollisionDiagonal writes. */
    std::vector<double> diag_H_pt, diag_l_pt, diag_phi_pt;
    /** [bin][l] gather/scatter buffers for the all-l kernel sweep, l contiguous within
     *  a bin. Sized lazily on first use like dy_coll, because the multipole count is a
     *  layout property not known at CreatePerturbScratch time. */
    std::vector<double> FH_all, Fl_all, Fphi_all, dFH_all, dFl_all, dFphi_all;
    /** Reduced-representation scratch: one moment vector and one grid-sized buffer for
     *  the reconstruct/project round trip. Per-workspace like everything else here —
     *  the k-loop is threaded and these are written on the RHS. */
    std::vector<double> red_m, red_grid;
    /** Tabulated path: the S-vector the mat-vec reads and the one it writes. */
    std::vector<double> red_z, red_dz;
  };
  std::unique_ptr<BaseSpecies::PerturbScratch> CreatePerturbScratch() const override;

  /** Pure worker behind AddCouplingDerivs: gathers the background PSDs from the
   *  interpolated columns in pvecback, runs the kernel once, and scatters the
   *  collision + background-evolution terms into dy. Split out (like
   *  ApplyKernelBackgroundDerivs) so the identity test can drive it with explicit
   *  layouts and a hand-built pvecback (no perturb_workspace). l_max_col already
   *  clamped to the children's l_max by the caller. */
  void ApplyKernelPerturbDerivs(double a,
                                int l_max_col,
                                const double* pvecback,
                                const NCDMBaseSpecies::PerturbLayout& p_lay,
                                const NCDMBaseSpecies::PerturbLayout& f_lay,
                                const NCDMBaseSpecies::PerturbLayout& b_lay,
                                const double* y,
                                double* dy,
                                Scratch& scratch,
                                double k = 0.) const;

  /** Same split as ApplyKernelPerturbDerivs, for the Jacobian diagonal: gathers the
   *  background PSDs from pvecback, runs ONE kernel pass, and writes each bin's
   *  d(dy_i)/d(y_i) into every l slot of that bin.
   *
   *  The diagonal is l-INDEPENDENT (pinned in decay_kernel_test): the parent's loss
   *  is a same-bin local term with no angular factor, and the daughters' loss was
   *  lumped so it couples a bin only to itself. One CollisionDiagonal call therefore
   *  serves the whole hierarchy.
   *
   *  Public for the same reason as the derivs worker: the test drives it with
   *  explicit layouts and a hand-built pvecback, with no perturb_workspace. */
  void ApplyKernelPerturbDiagonal(double a,
                                  int l_max_col,
                                  const double* pvecback,
                                  const NCDMBaseSpecies::PerturbLayout& p_lay,
                                  const NCDMBaseSpecies::PerturbLayout& f_lay,
                                  const NCDMBaseSpecies::PerturbLayout& b_lay,
                                  double* diag,
                                  Scratch& scratch) const;

  void PerturbDerivsDiagonal(const BaseSpecies::PerturbLayout& layout,
                             double tau,
                             const double* y,
                             double* diag,
                             const perturb_parameters_and_workspace& ppaw) const override;

  /** Builds the frozen reduced basis, and the tabulated operator when both are asked
   *  for, once the background table exists (BaseSpecies hook). No-op in the exact
   *  scheme (`dr_reduced_moments = 0`). */
  void ProcessBackgroundTable(const double* background_table,
                              int n_rows,
                              int row_stride,
                              const double* z_table) override;

  /** `dr_reduced_moments`: carry N q-moments per daughter instead of the full grid.
   *  0 (default) leaves the exact scheme untouched in every respect. */
  void set_reduced_moments(int n) {
    reduced_moments_ = n;
  }
  int reduced_moments() const {
    return reduced_moments_;
  }
  /** True once the frozen basis exists. Both conditions matter: the key can be set
   *  before ProcessBackgroundTable has run, and the perturbation RHS must not consult a
   *  basis that is not there yet. */
  bool reduced() const {
    return reduced_moments_ > 0 && reduced_op_ != nullptr;
  }
  void set_reduced_table_l_max(int l) {
    reduced_table_l_max_ = l;
  }
  /** True once the tabulated operator exists and may be used on the hot path.
   *
   *  The table is built whenever the reduction is on. The matrix-free path survives
   *  only as the fallback for the cases where the table cannot serve a request — a
   *  hierarchy asking for more multipoles than were tabulated, and a background on
   *  which the collision is never active enough to bound a window.
   *
   *  ⚠ MEMORY. One row is (l_max+1)·S² doubles — ~110 KB at S = 28 and 18 multipoles
   *  — and the window holds ~2000 background rows at Γ = 10⁷, so ~220 MB. That is
   *  one run's worth; it is NOT free to run 32 MCMC chains on a node. Thinning it by
   *  striding the rows was measured to save no wall time at all while costing
   *  4.1e-3 → 5.4e-1 in C_ℓ^TT across Γ = 10⁸→10¹⁰, so the row spacing is not the
   *  lever to pull — dr_N_q is. */
  bool tabulated() const {
    return reduced() && !reduced_table_.empty();
  }

  /** k_output_values time series: the sector's total anisotropic stress
   *  a⁴Π_νφ ≡ a⁴·Σ_{i∈{H,l,φ}} (ρ̄_i+p̄_i)σ_i (paper/thesis eq. 6.19), plus its
   *  per-child decomposition. Composite-owned like Type3Species::PrintVariables:
   *  the children carry no valid collection_index_, so their layouts are reached
   *  through THIS composite's nested layout. */
  void PrintVariables(PerturbColumnWriter& writer,
                      const BaseSpecies::PerturbLayout* base,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  // ── Test introspection (not on the hot path) ───────────────────────────────
  /** The parent's cached PERTURBATION-kernel background PSD rate f̄̇_H from the last
   *  ApplyKernelPerturbDerivs. Every species' Δdy is now already its F-space rate, so
   *  this is not needed to reconstruct one; what the test uses it for is the
   *  representation-independent decay-only statement, that pure decay leaves
   *  Ψ = F/f̄ invariant because dF/F equals f̄̇/f̄ bin by bin. Reads the PT kernel,
   *  not the background one: the comparison only holds against the same discrete
   *  network the ℓ-operator linearizes. */
  const std::vector<double>& KernelDfBgH(const Scratch& scratch) const {
    return scratch.kernel_pt->df_bg_H();
  }

  // ── Test introspection (not on the hot path) ───────────────────────────────
  /** Bare kernel moments of the derivs currently held in dy (read from the
   *  parent+daughter f-slots). The Fable-review conservation identities are
   *  N_H+N_l=0, 2N_H+N_phi=0, 2E_H+2E_l+E_phi=0 (×2 boson dof folded into df_phi).
   *  Precondition: dy's PARENT leg must be in STORED units (i.e. exactly what
   *  ApplyKernelBackgroundDerivs writes back, /kappa'd) — this method internally
   *  re-scales a scratch copy by kappa (#385) before handing it to the kernel,
   *  which itself only ever works in bare units; the daughter legs need no such
   *  conversion, they are already bare. */
  DecayTransitionKernel::Moments ConservationMoments(double a, const double* dy) const;
  /** Off-grid energy clamped by the BACKGROUND kernel (fine daughter grids). */
  double ClampedEnergyResidual() const {
    return kernel_->clamped_energy_residual();
  }
  /** Energy the BACKGROUND kernel's lumped daughter loss deliberately misplaces, in
   *  the same bare units ConservationMoments returns, so that
   *  2 E_H + 2 E_l + E_phi == SplitEnergyResidual() to machine precision. See
   *  DecayTransitionKernel::split_energy_residual() for why it is not zero. */
  double SplitEnergyResidual() const {
    return kernel_->split_energy_residual();
  }
  /** Energy / momentum the PERTURBATION kernel's lumped daughter loss misplaces, so
   *  that 2eH+2el+ephi and 2momH+2moml+momphi close. Reset by PrepareTransitions, so
   *  both survive a full l = 0,1,2,... sweep. */
  double SplitEnergyResidualPert(const Scratch& scratch) const {
    return scratch.kernel_pt->split_energy_residual_pert();
  }
  double SplitMomentumResidualPert(const Scratch& scratch) const {
    return scratch.kernel_pt->split_momentum_residual_pert();
  }
  /** Off-grid energy clamped by the PERTURBATION kernel (coarse daughter grids).
   *  Separate from the background one: the two kernels hold independent scratch, so
   *  after ApplyKernelPerturbDerivs only THIS one is fresh. */
  double ClampedEnergyResidualPert(const Scratch& scratch) const {
    return scratch.kernel_pt->clamped_energy_residual();
  }

 private:
  /** Builds the frozen reduced-moment basis and hands it to both daughters. Called from
   *  ProcessBackgroundTable, which is the first hook with a solved background — the basis
   *  is a functional of the daughter PSDs, so it cannot be built any earlier.
   *
   *  The reference row is the first at which the parent's comoving number has fallen to
   *  1/e of its initial value, i.e. the centre of the window the reduction has to be
   *  accurate in. A threshold crossing rather than an extremum of any rate: n_H falls
   *  monotonically through 1/e once and only once, where every rate-based criterion
   *  picks late times, when the parent's occupation sits on its floor and every
   *  log-derivative of it is floor noise. Anything hand-picked would be the
   *  (Gamma, m)-dependent tuning this whole scheme exists to avoid. */
  void BuildReducedBasis(const double* background_table, int n_rows, int row_stride);

  /** Assembles M~_l on every background row of the collision window and stores it.
   *  Called from ProcessBackgroundTable after BuildReducedBasis, which it needs: the
   *  table is the reduction of the operator ONTO the frozen basis, so the basis has to
   *  exist first.
   *
   *  The window is the CONTIGUOUS row range over which the collision diagonal rises
   *  above `kTableActiveRate * aH` and the parent still exists (kTableMinAbundance).
   *  Contiguous and not merely "the active rows",
   *  because the lookup interpolates between neighbours and a hole in the middle of the
   *  range would be interpolated ACROSS rather than skipped. Rows inside the window that
   *  are individually quiet are tabulated anyway; their blocks are near zero and cost
   *  only storage. */
  void BuildReducedTable(const double* background_table, int n_rows, int row_stride);

  /** Activity threshold, in units of the expansion rate, that bounds the window. The
   *  collision falls ten or more orders outside its ~1.5 decades (design note M1), so
   *  anything in this region is a truncation far below every other error in the scheme
   *  -- but it IS a truncation, and Apply drops to exactly zero outside, so it must be
   *  loose enough that the discontinuity is unobservable. */
  static constexpr double kTableActiveRate = 1e-4;

  /** Late cut: the parent's comoving number as a fraction of its initial value. Below
   *  this the collision has essentially nothing left to act on and every row past it is
   *  pure storage. This, not the rate, is what makes the window ~1.5 decades wide
   *  instead of the whole history -- the rate alone does NOT bound it on the late side,
   *  because the parent's diagonal tends to a*Gamma while aH falls, so rate/aH is still
   *  rising at a = 1. */
  static constexpr double kTableMinAbundance = 1e-6;

  DNCDMSpecies* parent_  = nullptr;
  DrPsdSpecies* fermion_ = nullptr;
  DrPsdSpecies* boson_   = nullptr;
  int reduced_moments_   = 0;
  /** `background_verbose`. The reduced path builds its basis and tabulates its operator
   *  during the background solve, and reports what it built at > 1. Warnings do not
   *  consult this. */
  int background_verbose_ = 0;
  /** Kept alive for the whole run: it owns the FROZEN basis (including the entropy
   *  weight, which the daughters do not carry), and the perturbation RHS calls its
   *  Reconstruct / Project on the hot path. Only const methods are used there, so it is
   *  shared across k-threads safely — Assemble is NOT called after construction. */
  std::unique_ptr<ReducedCollisionOperator> reduced_op_;
  /** ReducedCollisionOperator holds a non-owning reference to a kernel. The hot path
   *  only calls Reconstruct/Project, which never touch it, but a dangling reference
   *  would still be latent UB — so the kernel the basis was built against is kept
   *  alive here rather than being a local in BuildReducedBasis. It costs a few small
   *  vectors and is built once per run. */
  std::unique_ptr<DecayTransitionKernel> reduced_basis_kernel_;
  /** Multipole range the table must cover, resolved at Create time from the precision
   *  structure. ProcessBackgroundTable runs BEFORE RegisterPerturbationIndices, so the
   *  layouts do not exist yet and l_max cannot be read off them; it has to come from
   *  ppr. The RHS re-checks its own L against the table and falls back to the
   *  matrix-free path if the table is short, so a mismatch costs speed, never
   *  correctness. */
  int reduced_table_l_max_ = -1;
  /** M~_l over the window. Written once in ProcessBackgroundTable and read-only on the
   *  threaded k-loop thereafter, so it is shared rather than per-workspace. */
  ReducedOperatorTable reduced_table_;
  // TWO kernel instances over the SAME parent grid and the same Config, differing
  // only in which daughter grids they see. The daughters' background grid must stay
  // fine (it resolves the decay-injection band and carries the exact number/energy
  // conservation of the background RHS), while the perturbation hierarchy has to be
  // able to run on a coarser grid — N_pt·(l_max+1) variables per k-mode per daughter
  // is the whole cost of a run. Two instances rather than one parameterized kernel
  // because DecayTransitionKernel binds its grids (and all its per-node scratch) at
  // construction; it is deliberately left untouched.
  std::unique_ptr<DecayTransitionKernel> kernel_;  // daughters at q_bg()/dq_bg()
  // Kernel Config, kept so CreatePerturbScratch can build the per-k perturbation
  // kernels. There is deliberately NO shared perturbation kernel member: the only
  // instances live in per-workspace Scratch, so a racy one cannot be reached.
  DecayTransitionKernel::Config cfg_;
  std::vector<double> df_H_, df_l_, df_phi_;  // background hot-path scratch (sized in ctor)
  // Same, for BackgroundDerivsDiagonal. Separate from df_* because the evolver may
  // want both at one state and reusing the buffers would make the order load-bearing.
  std::vector<double> diag_H_, diag_l_, diag_phi_;
  // Parent-only kernel-boundary scratch (#385): fH_bare_[i] = kappa * y[parent slot
  // i], the bare occupation ApplyKernelBackgroundDerivs hands to kernel_ in place
  // of the raw (stored) f-slot pointer. Ctor-sized (parent q_size), never resized.
  std::vector<double> fH_bare_;

  const background* pba_       = nullptr;
  const BackgroundModule* bgm_ = nullptr;
};
