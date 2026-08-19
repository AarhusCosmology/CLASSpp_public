#pragma once
#include <cmath>
#include <cstdint>
#include <memory>
#include <vector>

// Non-owning view of a STATIC momentum grid. The kernel never owns these; the
// species owns the buffers (parent reuses its DNCDM grid, daughters own their
// log-trapezoid grids). q ascending, size n; dq = trapezoid cell widths.
struct GridView {
  const double* q  = nullptr;
  const double* dq = nullptr;
  int n            = 0;
};

// Statistics sign for the (1±f) factors: Fermion -> Pauli blocking (1-f);
// Boson -> Bose enhancement (1+f).
enum class Statistics : std::uint8_t { Fermion, Boson };

/** Occupation-number bounds shared by the kernel and the reduced operator.
 *
 *  kFMax keeps 1-f representable for a fermion, so the Pauli factor never underflows
 *  to exactly zero on a bin that is merely close to full. kFMin is the floor an
 *  unfilled daughter bin is read at; it is ~50 decades below any occupation the model
 *  produces, so it never binds on physics. */
inline constexpr double kFMin = 1e-100;
inline constexpr double kFMax = 1.0 - 1e-15;  // 1-f must stay representable

/**
 * Conservative discrete transition network for nu_H <-> nu_l + phi (massless
 * daughters). Discretizes the joint band integral ONCE per RHS evaluation and
 * derives all three background derivatives AND the per-l perturbation collision
 * operator from the SAME discrete transitions, so number+energy conservation and
 * detailed balance hold on the grid to machine precision for any resolution.
 * Physics: arXiv:2011.01502 eqs. (6.2)-(6.9); design spec section 4.
 *
 * Continuum band factor (design 4.1), with the cubic f_H f_l f_phi cancelling:
 *   Lambda(q1,q2,q3) = f_l(q2) f_phi(q3) - f_H(q1) + f_H(q1)[f_l(q2) - f_phi(q3)]
 *                      |--- inv ---|      |dec|     |--------- qs (linear) ---------|
 * Kinematics (massless daughters): eps1 = sqrt(q1^2+a^2 m^2); daughter band
 * q2 in [(eps1-q1)/2, (eps1+q1)/2]; partner q3 = eps1 - q2 (so q2*+q3* = eps1
 * exactly, making the two-bin linear deposit conserve number AND energy per
 * transition).
 */
class DecayTransitionKernel {
 public:
  struct Config {
    bool inverse_decays     = true;   // include the l+phi->H repopulation (f_l f_phi) term
    bool quantum_statistics = false;  // include the (1+-f) linear term f_H(f_l-f_phi)

    /** Stiffness cap, in units of the expansion rate: the fastest per-bin collision
     *  rate K/eps is limited to max_rate * (a'/a). 0 (the default) disables it.
     *
     *  Same device, and same argument, as the interacting-neutrino sectors
     *  (ncdm_interacting_species.cpp, ultra_relativistic.cpp): a process far faster
     *  than the expansion has already reached its attractor within a Hubble time. The
     *  margin is much thinner here, because this rate does NOT run far above H at the
     *  epoch that matters. Only DNCDMProxySpecies uses it (see `dr_rate_cap`); the
     *  exact method leaves it off.
     *
     *  MUST stay ONE global scale on K, never per bin: every deposit is linear in K,
     *  so a uniform factor leaves the conservation identities and detailed balance
     *  untouched, where a per-bin cap would rescale the legs of one transition
     *  differently and break them. */
    double max_rate = 0.;
  };

  struct Moments {
    double N_H, N_l, N_phi, E_H, E_l, E_phi;
  };

  /** fermion/boson are the daughters' STATE grids: where the perturbation F lives,
   *  where deposits land, and (by default) where the background f is read.
   *
   *  fermion_bg/boson_bg override ONLY where the background occupations are sampled.
   *  Pass them when the caller's state grid is coarser than the grid the background
   *  was integrated on. This is a CORRECTNESS requirement, not an optimisation: the
   *  daughters' PSD has an injection front, above which f falls off a cliff, and the
   *  gather below is a two-bin LINEAR interpolation (it has to be, so gather and
   *  scatter are transposes and the conservation identities hold). Interpolating
   *  across a cell that straddles the front returns a value set by the grid rather
   *  than by the physics -- orders of magnitude high -- and because the inverse-decay
   *  rate goes as f_l f_phi / f_H against a parent that has decayed to ~1e-16, that
   *  error turns into a fake repopulation. Sampling f on the background grid removes
   *  it for one extra binary search per node and nothing in the ODE dimension.
   *
   *  Empty views (the default) mean "same as the state grid", which reproduces the
   *  single-grid behaviour bit for bit. */
  DecayTransitionKernel(GridView parent,
                        GridView fermion,
                        GridView boson,
                        Statistics fermion_stat,
                        Statistics boson_stat,
                        Config cfg,
                        GridView fermion_bg = {},
                        GridView boson_bg   = {});

  /** Highest daughter momentum the transition network can ever emit, given the
   *  parent's top grid momentum q1_max and a*m at the LATEST scale factor of interest
   *  (a = 1 for a run that ends today).
   *
   *  Both daughters sit on the band q2 in [(eps1-q1)/2, (eps1+q1)/2], q3 = eps1 - q2,
   *  so both are bounded by (eps1+q1)/2, which increases with q1. CLOSED FORM and
   *  SATURATING -- q* ~ a*m/2 grows linearly in a and stops at a = 1 -- so sizing a
   *  daughter grid is a calculation, not an open problem.
   *
   *  Undersizing it is an ENERGY LEAK rather than a resolution trade-off: TwoBinSplit
   *  clamps the off-grid target into the edge bin, conserving NUMBER but filing the
   *  energy at q_edge instead of q*, and it is silent because the identities the
   *  kernel asserts are number identities. clamped_energy_residual() is the detector.
   *  Callers add a small margin (the composite uses 5%) because the Gauss-Legendre
   *  nodes approach the band edge without reaching it. */
  static double MaxDaughterMomentum(double q1_max, double a_m) {
    return 0.5 * (std::sqrt(q1_max * q1_max + a_m * a_m) + q1_max);
  }

  /** Where one transition's daughter momentum q* is placed on the daughter grid:
   *  two consecutive weights starting at bin `j0`. `clamped` says q* fell off the
   *  grid and the energy was filed at the edge bin rather than at q*.
   *
   *  EVERY consumer (background scatter, lumped loss, positivity decomposition and
   *  the perturbation gather/scatter) uses this one object, and the gather must use
   *  the SAME weights as the scatter -- that transpose is what makes the discrete
   *  number and energy identities exact. Three invariants, all asserted in
   *  decay_kernel_test:
   *      sum_e w_e       == 1     (number)
   *      sum_e w_e q_je  == q*    (energy, per transition)
   *      w_e             >= 0     (an empty bin may only fill; see the lumped loss)
   */
  struct DepositStencil {
    static constexpr int kWidth = 2;
    int j0                      = 0;
    double w[kWidth]            = {0., 0.};
    bool clamped                = false;
  };

  /** Build the two-bin deposit stencil for q* on grid g.
   *
   *  `hint`, when non-null, is read as the bracket index this search should start its
   *  outward hunt from and written back with the one it found, so a caller stepping
   *  q_star monotonically pays O(1) per call instead of a fresh bisection. Purely a
   *  search accelerator: the stencil is unique, so it does not affect the result. */
  static void BuildDeposit(const GridView& g,
                           double q_star,
                           DepositStencil& out,
                           int* hint = nullptr);

  // Background RHS. Zeroes df_* then accumulates. K = a^2 m Gamma (design 4.1).
  // f_*/df_* sized to the matching grid's n. Pure; hot-path scratch preallocated.
  // a_prime_over_a drives the Config::max_rate stiffness cap; 0 (the default)
  // leaves the rate uncapped, which is what the kernel unit tests want.
  void ComputeBackgroundDerivs(double a,
                               double m,
                               double Gamma,
                               const double* f_H,
                               const double* f_l,
                               const double* f_phi,
                               double* df_H,
                               double* df_l,
                               double* df_phi,
                               double a_prime_over_a = 0.) const;

  // Diagnostics: N_i = factor-free Sum dq q^2 df_i, E_i = Sum dq q^2 eps_i df_i.
  // The kernel returns the BARE per-dof grid moments; the g-factor spin weights
  // are applied by the caller (the g_H*g_l/g_phi = 2 boson ratio is already
  // folded into df_phi at deposit time, paper eq. 4.14).
  Moments ComputeMoments(
      double a, double m, const double* df_H, const double* df_l, const double* df_phi) const;

  // Precompute the transition geometry (node momenta q2*/q3*, deposit bins+split
  // weights on both daughter grids, gathered background f, Lambda_s, and the angle
  // cosines cos_a*/cos_b*/cos_g*) for the current (a,m,Gamma,f) into internal
  // scratch. Background derivs call this internally; the perturbation module calls
  // it once per RHS so ApplyPerturbationOperator can run for every l without
  // recomputing kinematics. Resets clamped_energy_residual().
  void PrepareTransitions(double a,
                          double m,
                          double Gamma,
                          const double* f_H,
                          const double* f_l,
                          const double* f_phi,
                          double a_prime_over_a = 0.) const;

  // Per-l perturbation collision operator on the SAME network (design 4.3). Needs a
  // prior PrepareTransitions for the current (a,m,f), and MUST then be called for
  // l = 0,1,2,... in strict ascending sequence (enforced): the angular factors are
  // advanced by a rolling per-node recurrence rather than rebuilt from P_0, which is
  // what makes the caller's l-loop cost O(L) instead of O(L^2).
  //
  // PURE UNNORMALIZED F-SPACE: in F_i = f_i Psi_i = delta-f_i, out
  // dF_i = (dF_i/dtau)^(1)_{C,l}. It applies NO 1/f_i normalization -- a species whose
  // own state variable is normalized (the parent's Psi_H) converts on BOTH sides in
  // the caller. That is load-bearing: a daughter filling from zero abundance has an
  // enormous f-bar-dot/f-bar at early times, so a 1/f_i inside the operator would blow
  // the hierarchy up, whereas F = 0 is the exact statement "no perturbation to an
  // empty distribution". It also omits the -(f-bar-dot_i/f_i)Psi background-evolution
  // term; the composite adds that for the PARENT only, from df_bg_*(), so the
  // pure-decay piece cancels exactly in decay-only mode (paper 4.3).
  //
  // Zeroes dF_* then accumulates, gathering and scattering with transpose weights so
  // the l=0 number/energy and l=1 momentum identities hold to machine precision.
  void ApplyPerturbationOperator(int l,
                                 const double* F_H,
                                 const double* F_l,
                                 const double* F_phi,
                                 double* dF_H,
                                 double* dF_l,
                                 double* dF_phi) const;

  /** Rewind the rolling multipole recurrence WITHOUT rebuilding the transition
   *  geometry, so the next ApplyPerturbationOperator call must again be l = 0.
   *
   *  The only state the recurrence carries between sweeps is next_l_, and the l = 0
   *  call re-seeds P_0 for every node, so this is a legitimate one-line rewind rather
   *  than a partial reset with a hidden precondition. It exists so that assembling the
   *  operator as an explicit matrix can hoist PrepareTransitions out of the per-column
   *  loop. Nothing else may be assumed reset. */
  void ResetMultipoleRecurrence() const {
    next_l_ = 0;
  }

  /** The same operator for ALL multipoles l = 0..l_max in one pass, with the loops
   *  inverted: node outer, l inner. Arrays are [bin][stride], l contiguous within a bin.
   *
   *  The per-l entry point above re-walks the whole node structure once per multipole,
   *  which measured as the dominant cost of the module and behaved as a memory-bound
   *  loop (removing arithmetic from the inner loops changed nothing). Here each node's
   *  data is loaded once and reused across all multipoles, the background-only
   *  coefficients c_H/c_l/c_phi are formed once per node, and the rolling Legendre pair
   *  lives in registers.
   *
   *  IDENTICAL RESULTS, not merely equivalent: for a given (bin, l) the contributions
   *  still arrive in ascending i then ascending n, so the summation order is unchanged.
   *  It replaces the whole l-loop, starting from the post-PrepareTransitions state and
   *  leaving next_l_ where an l_max-long sequence of per-l calls would have.
   *
   *  WITHIN a node, l is the INNERMOST dimension and every loop body is stride-1 over
   *  it. That is why the Legendre factors are evaluated into small arrays up front
   *  rather than rolled in registers: the recurrence is sequential in l and would
   *  otherwise be a vectorisation barrier in the hottest loop in the module.
   *
   *  The dF_* buffers MUST NOT alias the F_* ones (they are declared __restrict, and
   *  the function zeroes dF_* before accumulating, so an alias would destroy its own
   *  input). Callers pass distinct per-workspace scratch. */
  void ApplyPerturbationOperatorAllL(int l_max,
                                     int stride,
                                     const double* __restrict F_H,
                                     const double* __restrict F_l,
                                     const double* __restrict F_phi,
                                     double* __restrict dF_H,
                                     double* __restrict dF_l,
                                     double* __restrict dF_phi) const;

  /** Per-bin decomposition of a daughter's background RHS with respect to that bin's
   *  OWN occupation. Given the parent's f, the daughters' RHS is LINEAR in each
   *  daughter's f-vector, so for daughter bin k
   *
   *      df_k = gain_k + diag_k * f_k,
   *
   *  with gain_k the transitions' f-independent part. There is no off-diagonal term:
   *  the lumped loss couples a bin only to itself, which is exactly what makes the
   *  operator positivity-preserving.
   *
   *  WHY IT EXISTS. gain_k is the RHS the bin would see at f_k = 0. If it is negative,
   *  the operator drives an empty bin to a NEGATIVE occupation -- a property of the
   *  right-hand side, not of the integrator, so no step-size control, state clamp or
   *  change of variable can prevent it. That is the distinction between "the solver
   *  overshot" and "the discretisation is not positivity-preserving".
   *
   *  Requires a prior PrepareTransitions for the state being diagnosed. */
  struct PositivityTerms {
    std::vector<double> gain, diag;
    void Resize(int n) {
      gain.assign(n, 0.);
      diag.assign(n, 0.);
    }
  };
  void PositivityDecomposition(PositivityTerms& fermion_out, PositivityTerms& boson_out) const;

  /** Per-bin d(df_i)/d(f_i) of the collision RHS, for all three species at once.
   *
   *  This is the stiff part of the network, and it is EXACTLY the diagonal: Lambda is
   *  linear in each species' own occupation, and the daughter loss is lumped precisely
   *  so that no off-diagonal survives. Writing the collision as
   *
   *      df_i = gain_i + diag_i * f_i,
   *
   *  diag_i is what an exponential / IMEX step integrates exactly instead of resolving.
   *  For the parent d(Lambda)/d(f_H) = -1 + (f_l - f_phi) under quantum statistics and
   *  -1 without, so diag_H -> -K/eps1 -> -a*Gamma, which measures as the dominant
   *  eigenvalue of the whole network over the epochs that carry the steps.
   *
   *  No allocation and no gain term: this is a hot path, called once per attempted
   *  step, unlike PositivityDecomposition which is a throttled diagnostic. The two
   *  necessarily share structure, so decay_kernel_test pins them against a
   *  finite-difference Jacobian rather than against each other.
   *
   *  Units: invariant under the caller's stored/bare occupation rescaling, since f and
   *  df are scaled by the same constant. Requires a prior PrepareTransitions. */
  void CollisionDiagonal(double* diag_H, double* diag_l, double* diag_phi) const;

  // Cached background PSD rates df_i/dtau from the last PrepareTransitions (SAME
  // discrete rates the operator linearizes). The composite forms -(f-bar-dot/f)Psi
  // from these — never a separately-interpolated background column — which is what
  // makes the parent's decay term cancel exactly.
  const std::vector<double>& df_bg_H() const {
    return df_bg_H_;
  }

  // Positivity floor for the callers' 1/f conversions (mirrors DrPsdSpecies::kFFloor).
  // The operator itself never divides by f; the composite uses this to guard the
  // parent's Psi_H <-> F_H round trip.
  static constexpr double kFFloor = 1e-100;

  // Energy mis-deposited by off-grid clamping since the last PrepareTransitions.
  // Zero iff every node's q2*/q3* land inside the daughter grids (no clamping).
  double clamped_energy_residual() const {
    return clamped_energy_residual_;
  }

  /** Energy misplaced by the LUMPED daughter loss since the last PrepareTransitions,
   *  in the same bare-moment units as ComputeMoments, signed so that
   *
   *      2 E_H + 2 E_l + E_phi == split_energy_residual()
   *
   *  holds to machine precision. Sibling of clamped_energy_residual(): both name an
   *  energy the network does not place at q*, and neither is allowed to be silent.
   *
   *  WHY THERE IS ONE AT ALL. A two-bin deposit that conserves both number and energy
   *  is two equations in two unknowns, so its weights are FORCED. Applied to a loss,
   *  which is proportional to the interpolated occupation, those weights bill a bin
   *  for particles its neighbour holds -- a negative off-diagonal, and an empty bin
   *  driven below zero. Exact per-transition energy placement and positivity are
   *  therefore mutually exclusive here, and this branch buys positivity. The total
   *  number removed is unchanged, so every number identity and the whole parent leg
   *  are untouched; only WHERE in q the energy leaves differs, by O(dq), and only
   *  where f varies sharply across a cell. */
  double split_energy_residual() const {
    return split_energy_residual_;
  }

  /** The same misplacement for the PERTURBATION operator, so that its two q-weighted
   *  identities close:
   *
   *      l = 0:  2 eH + 2 el + ephi        == split_energy_residual_pert()
   *      l = 1:  2 momH + 2 moml + momphi  == split_momentum_residual_pert()
   *
   *  Kept as two accumulators rather than one per-call value because callers run a
   *  whole l = 0,1,2,... sweep between PrepareTransitions calls, which would overwrite
   *  a single one. Both are reset by PrepareTransitions and written only by the l = 0
   *  and l = 1 calls; l >= 2 has no conservation identity to book against.
   *
   *  The l = 0 NUMBER identities (numH+numl, 2numH+numphi) stay EXACTLY zero: lumping
   *  moves where a deposit lands, never how much. */
  double split_energy_residual_pert() const {
    return split_energy_residual_pert_;
  }
  double split_momentum_residual_pert() const {
    return split_momentum_residual_pert_;
  }

  // -- Test introspection (do not use on the hot path) -----------------------
  // Legendre P_l(x) by upward recurrence from P_0, the reference form of the operator's
  // angular factor (exposed so the test can pin that P_l is built by recurrence, never
  // read from a table). The operator itself runs the identical recurrence incrementally
  // across its l = 0,1,2,... calls, one step per call, returning the same values to
  // within the last bit (see the ulp note at that loop).
  static double LegendreP(int l, double x);
  /** Nodes laid out for parent bin i by the last PrepareTransitions. The count is not
   *  fixed: the union partition cuts the emission band wherever either daughter grid
   *  has a cell edge inside it, and how much of the grids the band covers depends on
   *  a*m. */
  int node_count(int i) const {
    return node_cnt_[i];
  }
  /** Quadrature weight of node s of parent bin i, i.e. the dq2 the transition carries.
   *  Sum over s is q1, the band width. */
  double node_weight(int i, int s) const {
    return wn_[node_off_[i] + s];
  }
  /** Kinematics of node s of parent bin i, from the last PrepareTransitions. */
  void NodeKinematics(
      int i, int s, double& q2, double& q3, double& cos_a, double& cos_b, double& cos_g) const;

 private:
  GridView parent_, fermion_, boson_;
  // Where f_l / f_phi are READ (see the ctor doc). Default-initialized to the state
  // grids, so every path that does not opt in behaves exactly as before.
  GridView fermion_bg_, boson_bg_;
  // False when both bg views alias the state grids: then the deposit split already
  // IS the gather split and the extra binary searches are skipped entirely, which is
  // what keeps the single-grid path bit-identical.
  bool separate_bg_grids_ = false;
  Statistics fermion_stat_, boson_stat_;
  Config cfg_;

  // Per-node scratch, sized parent_.n * node_stride_ in the ctor (allocation-free RHS).
  mutable std::vector<double> q2_star_, q3_star_;        // node momenta, q2*+q3* = eps1
  mutable std::vector<double> wn_;                       // node weight 0.5(qhi-qlo)*gl_w
  mutable std::vector<double> lambda_s_;                 // Lambda_s (band factor)
  mutable std::vector<double> fl_gather_, fphi_gather_;  // gathered f_l(q2*), f_phi(q3*)
  // Per-node lumping fraction, one value per node and daughter (see LumpTheta).
  mutable std::vector<double> theta_l_, theta_phi_;

  /** The lumping fraction theta for one node/daughter.
   *
   *  Full lumping (theta = 1) is what buys positivity -- the loss charged to a bin must
   *  vanish as THAT bin empties -- but it is also a one-signed energy error, and the
   *  dominant term in the daughter-grid resolution requirement. Its size per node is a
   *  covariance over the stencil,
   *
   *      sum_e w_e q_e d_e = -lambda(1-lambda) (q_{j+1}-q_j) (f_j - f_{j+1}),
   *
   *  so it vanishes on a FLAT stencil and is worst across a steep cell, while
   *  positivity only cares about the stencil's MINIMUM. Those demands do not conflict,
   *  and theta is the blend: it replaces the deviation d_e = f_e - f_gathered by
   *  theta*d_e, theta = 1 being the fully lumped loss and theta = 0 the energy-exact
   *  two-bin split.
   *
   *  theta MUST be one value for the whole stencil, because sum_e w_e (theta d_e) = 0
   *  needs theta to come out of the sum -- and that identity is what keeps dN, the
   *  parent's leg and every number identity exactly unchanged for any theta. A
   *  per-EDGE theta would break all of them; it is the one real constraint here.
   *
   *  g <= 0 exempts the node outright: the "loss" is then a gain, so its off-diagonal
   *  is positive and the exact split is already Metzler. */
  double LumpTheta(double g, const double* w, int n, int j0, const double* f, double f_gath) const {
    // Nodes with g <= 0 are exempt: there the "loss" is a gain, its off-diagonal is
    // POSITIVE, and Metzler is satisfied by the exact split already.
    if (!(g > 0.))
      return 0.;
    // theta = 1 - H/A, the harmonic mean of the stencil over its arithmetic mean.
    // H <= A always (AM-HM), with equality iff the stencil is flat, so theta is in
    // [0,1] by construction and needs no clamp on physical input.
    //
    // Harmonic rather than 1 - min/A, which was the obvious first choice and is
    // measurably worse: across a stencil f(1 +- delta) the min form gives
    // theta = delta while this gives theta = delta^2, so it lumps far less on a smooth
    // stretch, and at a front it still saturates because H is dominated by the
    // smallest entry. It is also SMOOTH, which min is not -- min has a kink where the
    // argmin switches, and a kink in the RHS is what stalled earlier attempts. That
    // smoothness is what makes d(theta)/df tractable for the analytic Jacobian.
    double S = 0.;
    for (int e = 0; e < n; ++e) {
      const double fe = f[j0 + e];
      if (!(fe > 0.))
        return 1.;  // an exactly empty bin: lump fully, which is where positivity needs it
      S += w[e] / fe;
    }
    if (!(S > 0.) || !(f_gath > 0.))
      return 1.;
    const double r = 1. / (S * f_gath);  // H / A
    return r >= 1. ? 0. : (r <= 0. ? 1. : 1. - r);
  }

  /** d(theta)/d(f_k) for the harmonic limiter above. CollisionDiagonal needs it
   *  because theta is a function of the STATE, so the analytic diagonal must carry
   *
   *      d(Lambda_k)/df_k = g [ theta + (1-theta) dG/df_k + (f_k - G) d(theta)/df_k ].
   *
   *  With H = 1/S and A = f_gath:  d(theta)/df_k = (w_k/A) (H/A - H^2/f_k^2).
   *
   *  It vanishes identically on a flat stencil (H = A = f_k), as it must, and stays
   *  finite as f_k -> 0 (H -> f_k/w_k, so the whole expression tends to -1/(w_k A)).
   *  Returns 0 when the limiter is not the active scheme, where theta is constant. */
  double LumpDTheta(
      double g, const double* w, int n, int j0, const double* f, double f_gath, int e) const {
    if (!(g > 0.) || !(f_gath > 0.))
      return 0.;
    double S = 0.;
    for (int t = 0; t < n; ++t) {
      const double ft = f[j0 + t];
      if (!(ft > 0.))
        return 0.;  // saturated at theta = 1, locally constant
      S += w[t] / ft;
    }
    if (!(S > 0.))
      return 0.;
    const double H  = 1. / S;
    const double fk = f[j0 + e];
    return (w[e] / f_gath) * (H / f_gath - H * H / (fk * fk));
  }
  mutable std::vector<double> cos_a_, cos_b_, cos_g_;  // emission-angle cosines (Phase 4)
  // Rolling Legendre state per node, one (P_{l-2}, P_{l-1}) pair per cosine. The
  // operator advances these by ONE recurrence step per call instead of rebuilding
  // P_l from P_0 (which cost O(l) per node per call, i.e. O(L^2) over the caller's
  // l = 0..L loop and ~97% of the perturbation RHS). The arithmetic is the same
  // relation in the same operand order as LegendreP, agreeing with it to the last bit;
  // what it requires in exchange is that l advance by exactly 1 per call, tracked by
  // next_l_ and enforced in ApplyPerturbationOperator.
  mutable std::vector<double> Pa_prev_, Pa_curr_, Pb_prev_, Pb_curr_, Pg_prev_, Pg_curr_;
  mutable int next_l_ = 0;  // next l the operator will accept
  // Node layout. Nodes for parent i occupy [node_off_[i], node_off_[i] + node_cnt_[i]).
  // The stride is the worst case, so the offsets never move; only the COUNT varies,
  // because the number of sub-intervals in the union partition depends on how much of
  // the daughter grids the band covers and therefore on a*m.
  int node_stride_ = 0;
  mutable std::vector<int> node_off_, node_cnt_;
  mutable std::vector<double> breakpoints_;  // scratch for the union partition

  mutable std::vector<int> j_fermion_, k_boson_;  // low deposit bin per daughter
  // Deposit weights, kWidth per node. The SAME array is the gather weight:
  // transpose consistency is what makes the discrete identities exact.
  mutable std::vector<double> w_fermion_, w_boson_;
  mutable std::vector<double> eps1_;  // per parent bin, size parent_.n
  // Background PSD rates cached by PrepareTransitions (sized to each grid in ctor);
  // ComputeBackgroundDerivs copies these out, the operator reuses them for -f-dot/f.
  mutable std::vector<double> df_bg_H_, df_bg_l_, df_bg_phi_;
  // Daughter background f resampled onto the DEPOSIT (state) grid. Identical to the
  // caller's f_l/f_phi unless separate bg views were given; in that case the state
  // grid is an exact integer subsample of the bg grid (DrPsdSpecies::bg_index), so
  // the two-bin gather at a state momentum lands exactly on a bg node and this is a
  // lookup, not an interpolation. Needed because the lumped loss reads the occupation
  // of the bin it deposits INTO, which the bg-grid array cannot index.
  mutable std::vector<double> fl_state_, fphi_state_;
  // 1/(q^2 dq) on the two daughter STATE grids, formed once at construction because
  // the grids never move. The scatter divides by this measure once per (node,
  // deposit edge, MULTIPOLE), and the compiler cannot hoist that out of the l loop:
  // the measure is read through fermion_.q / boson_.q, which it cannot prove are not
  // aliased by the dF_l / dF_phi it is storing into on the same iteration. So the
  // innermost body of the hottest loop in the module carried three loads, two
  // multiplies and a divide that depend on nothing but the bin index.
  std::vector<double> inv_meas_f_, inv_meas_b_;
  // Per-node, per-deposit-edge occupation DEVIATION f[k] - (state-grid gather), the
  // one background quantity the lumped scheme needs beyond the gathers the operator
  // already loads. Every lumped coefficient is a gathered one plus a multiple of
  // this, so storing the deviation rather than the four derived coefficients halves
  // the extra traffic on the l-loop, which is the hottest path in the module.
  // sum_e w_e (deviation) == 0 by construction -- that identity is what keeps number
  // conservation exact.
  mutable std::vector<double> d_l_edge_, d_phi_edge_;  // size total*kWidth
  // Companions to d_*_edge_ for the limiter, same layout. dev_* is the RAW deviation
  // f_k - G (d_*_edge_ is theta times this) and dth_* is d(theta)/df_e. The
  // perturbation operator needs both to carry the (f_k - G) delta-theta term that
  // makes it the exact Jacobian of the limited loss.
  mutable std::vector<double> dev_l_edge_, dev_phi_edge_;
  mutable std::vector<double> dth_l_edge_, dth_phi_edge_;
  // Parent background f cached by PrepareTransitions (the caller owns the buffer;
  // valid for the duration of one RHS, i.e. until the next PrepareTransitions). Only
  // f_H is kept: the operator needs it for the quantum-statistics coefficients, while
  // the node-gathered daughter f's already live in fl_gather_/fphi_gather_.
  mutable const double* fH_bg_                 = nullptr;
  mutable double k_coeff_                      = 0.;  // a^2 m Gamma
  mutable double am2_                          = 0.;  // a^2 m^2
  mutable double clamped_energy_residual_      = 0.;
  mutable double split_energy_residual_        = 0.;
  mutable double split_energy_residual_pert_   = 0.;
  mutable double split_momentum_residual_pert_ = 0.;

  // Two-bin linear split of q_star on g: (j, weight lambda into j, 1-lambda into
  // j+1) with lambda q[j] + (1-lambda) q[j+1] = q_star (momentum-exact). Off-grid
  // q_star clamps fully into the nearest edge bin (number-preserving) and sets
  // clamped = true so the caller can book the energy mismatch.
  //
  // `hint` is a previously returned j to start an outward hunt from, or -1 for a plain
  // bisection. It changes the SEARCH ONLY -- the bracket is unique, so the returned
  // (j, lambda, clamped) are identical either way. It pays because the node momenta
  // arrive MONOTONICALLY (q2* ascending, q3* = eps1 - q2* descending), making a hunt
  // from the previous bracket O(1) amortised where the bisection was a measurable
  // share of a high-Gamma run.
  static void TwoBinSplit(
      const GridView& g, double q_star, int& j, double& lambda, bool& clamped, int hint = -1);
};
