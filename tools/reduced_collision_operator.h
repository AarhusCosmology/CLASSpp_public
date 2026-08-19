#pragma once
#include <cstdint>
#include <vector>

#include "decay_transition_kernel.h"

/**
 * Galerkin reduction of the per-multipole collision operator for
 * nu_H <-> nu_l + phi onto a low-dimensional daughter basis.
 *
 * WHAT IT IS FOR. The daughters carry N_q x (l_max+1) perturbation variables each --
 * the overwhelming majority of a shipped run's state -- and the grid they need scales
 * with Gamma. The parent does not, and a fluid closure on it is measured to fail. So
 * the reduction to attempt is daughters only.
 *
 * The reduced operator is built AS A CONGRUENCE OF THE EXACT ONE rather than
 * re-derived:
 *
 *      M~_l  =  P . M_l . R
 *
 * with M_l the operator DecayTransitionKernel::ApplyPerturbationOperator already
 * applies, R a reconstruction from daughter moments to a daughter grid function, and
 * P the moment map back. Every coefficient is read out of the tested kernel, so a
 * disagreement cannot come from a second implementation of the physics. Three
 * consequences follow, and they are the whole design -- break any of them and the
 * guarantees go with it:
 *
 *  1. CONSERVATION IS INHERITED. The sector's conserved functionals are moments
 *     (number Sum dq q^2 df; energy and, at l = 1, momentum Sum dq q^3 df), and those
 *     ARE the retained coordinates, so w^T P = w^T exactly. Whatever conservation the
 *     exact operator has, the reduced one has, for ANY daughter DOF count >= 2.
 *
 *  2. THE QUADRATURE GRID DECOUPLES FROM THE STATE DIMENSION. M_l is assembled on a
 *     daughter grid that exists only inside this class; the ODE state is the moment
 *     vector. The lumped daughter loss misplaces momentum by O(dq), and that residual
 *     drives the spurious positive l = 1 eigenvalue -- so it can be driven down by
 *     refining a grid the integrator never sees. That separation of accuracy from
 *     cost, not an analytic identity, is what breaks the N ~ Gamma^0.25 rule.
 *
 *  3. DISSIPATIVITY IS INHERITED IF P AND R ARE ADJOINT. The linearised operator is
 *     symmetric negative semi-definite in the entropy inner product
 *     <df, df'> = Sum dq q^2 df df' / w(q), w = f(1-f) for a fermion and f(1+f) for a
 *     boson. Reconstructing in that weight makes M~ a congruence of M, so
 *     max Re lambda(M~) <= max Re lambda(M) and the reduction cannot manufacture a
 *     growing mode. DissipativityBound() measures it.
 *
 * THE BASIS. R reconstructs df_i(q) = w_i(q) Sum_j c_j psi_j(q). Under the entropy
 * weight, n_moments = 2 spans {w, w q}, which is EXACTLY the space of equilibrium
 * perturbations (df = (df/dmu) dmu + (df/dT) dT), so the basis annihilates them
 * identically. That matters because the occupations here are extreme, making Lambda a
 * small residual of large terms: a basis that breaks the cancellation injects error
 * proportional to the large terms. EquilibriumNullVectors against a detailed-balance
 * background is the machine-precision test of it.
 *
 * UNITS. Everything here is BARE, because that is what the kernel speaks. A caller
 * holding the parent leg in the species' STORED units must apply
 * kappa = deg (2pi)^3/2 on both sides first.
 *
 * Derivation and the measurements behind the basis and DOF choices:
 * docs/superpowers/specs/2026-08-12-dncdm-reduced-perturbation-model-design.md.
 */
class ReducedCollisionOperator {
 public:
  struct Config {
    /** Daughter degrees of freedom per species. 2 is the design's proposal and the
     *  smallest that can span the equilibrium manifold; 3 and 4 are the honest
     *  convergence dimension (the exact scheme is the limit). */
    int n_moments = 2;
  };

  /** `kernel` must be constructed on `parent`/`fermion`/`boson` and is driven, not
   *  owned; the caller keeps it alive and calls PrepareTransitions on it before
   *  Assemble. The grids are the ones the kernel's perturbation state lives on --
   *  here they are pure quadrature, so they may be far finer than any grid the
   *  integrator could afford. */
  ReducedCollisionOperator(const DecayTransitionKernel& kernel,
                           GridView parent,
                           GridView fermion,
                           GridView boson,
                           Statistics fermion_stat,
                           Statistics boson_stat,
                           Config cfg);

  /** Builds the basis (weights, q_ref, Gram inverse) for one background state. Must
   *  be called with the same (a, m, f) the kernel was prepared with; nothing here
   *  re-derives the kernel's own coefficients, but the BASIS is a functional of the
   *  background and would otherwise be stale. */
  void SetBackground(double a, double m, const double* f_H, const double* f_l, const double* f_phi);

  /** Point this operator at THIS row's background but keep `src`'s FROZEN basis.
   *
   *  Assembling a table row by row must not rebuild the basis at each row: the moments
   *  are coordinates, and a basis that tracked `a` would make row n's coordinates mean
   *  something different from row n+1's, so interpolating between them -- or evolving a
   *  state through them -- would be meaningless. Only the OPERATOR may vary with a.
   *
   *  This is the reason SetBackground and this are separate calls rather than one with a
   *  flag: the failure it prevents is silent. Everything still runs, every identity in
   *  this class still holds row by row, and only the table is wrong. */
  void AdoptFrozenBasis(double a, double m, const double* f_H, const ReducedCollisionOperator& src);

  int n_moments() const {
    return cfg_.n_moments;
  }
  /** Reduced state dimension: q-resolved parent + n_moments per daughter. */
  int size() const {
    return parent_.n + 2 * cfg_.n_moments;
  }
  /** Full (exact-operator) state dimension, for the side-by-side controls. */
  int full_size() const {
    return parent_.n + fermion_.n + boson_.n;
  }

  /** Assemble M~_l for l = 0..l_max. `out` is resized to (l_max+1)*size()*size();
   *  block l starts at l*size()*size() and is ROW-MAJOR, so element (i, j) of M~_l
   *  is out[l*S*S + i*S + j] with S = size().
   *
   *  Cost is size() applications of the exact operator per background state -- about
   *  20 columns at the shipped parent grid, against the 118-602 a full-space
   *  precompute pays. It is k-independent, so it amortises over the whole k-loop. */
  void Assemble(int l_max, std::vector<double>& out) const;

  /** The same assembly for the EXACT operator, so every test below can be run on
   *  both without a second harness. `out` is (l_max+1)*full_size()^2, row-major. */
  void AssembleExact(int l_max, std::vector<double>& out) const;

  /** Residual of the sector's conservation identities, as the relative size of the
   *  left-null-vector product ||w^T M|| / (||w|| ||M||). Zero means the identity
   *  holds on this operator at this resolution.
   *
   *  The identities are the kernel's own (decay_kernel_test): N_H + N_l = 0 and
   *  2 N_H + N_phi = 0 (number, l = 0), 2 E_H + 2 E_l + E_phi = 0 (energy, l = 0),
   *  and the same q-weighting at l = 1 is momentum -- the identity whose discrete
   *  violation owns the spurious eigenvalue. `momentum` is only meaningful for the
   *  l = 1 block, the other three only for l = 0; all four are computed regardless
   *  so a caller can watch them across l. */
  struct ConservationResidual {
    double number_fermion = 0.;
    double number_boson   = 0.;
    double energy         = 0.;
    double momentum       = 0.;
  };
  ConservationResidual Conservation(const double* M_l) const;
  ConservationResidual ConservationExact(const double* M_l) const;

  /** max Re lambda(M) <= this, rigorously: the largest eigenvalue of the symmetric
   *  part of M in the entropy inner product. A value <= 0 PROVES the operator cannot
   *  grow any mode -- which is the property the shipped discrete operator lacks at
   *  l = 0 and l = 1, and the property a reduced operator has to keep at every
   *  number of degrees of freedom for its grid requirement to stop scaling with
   *  Gamma.
   *
   *  Returned in the same units as the operator (1/Mpc), so compare it against a*H. */
  double DissipativityBound(const double* M_l) const;
  /** The same bound for the EXACT operator, in the same physical norm. This is the
   *  control the reduced number is meaningless without.
   *
   *  The congruence P = R-adjoint makes the reduced Rayleigh quotient at z equal the
   *  exact one at R z, so the reduced numerical range is a SUBSET of the exact one and
   *  therefore
   *                DissipativityBound() <= DissipativityBoundExact()
   *  identically, at every background and every resolution. That inequality -- not
   *  "the reduced bound is negative" -- is the milestone: it says the reduction cannot
   *  manufacture a growing mode, so whatever the exact operator's defect is, refining
   *  the QUADRATURE grid (which the reduced state does not pay for) removes it. */
  double DissipativityBoundExact(const double* M_l) const;

  /** The three equilibrium perturbations, in reduced coordinates: two chemical-
   *  potential shifts obeying mu_H = mu_l + mu_phi and one temperature shift. On a
   *  background in detailed balance the exact operator annihilates all three, so
   *  ||M~ v|| measures whether the reduction preserves detailed balance.
   *
   *  Each vector is size(); `out` is resized to 3*size(). */
  void EquilibriumNullVectors(std::vector<double>& out) const;
  /** Same three vectors on the full grid, size 3*full_size(). */
  void EquilibriumNullVectorsExact(std::vector<double>& out) const;

  /** Reconstruction and projection, exposed so an operator-against-operator probe can
   *  drive the exact kernel with a df that lies in the reduced span (where the two
   *  must agree to machine precision) or with a real df (where the difference is the
   *  closure error being measured). `daughter` is 0 for the fermion, 1 for the boson.
   *  `coeff` has n_moments entries, `df` has that daughter's grid size. */
  void Reconstruct(int daughter, const double* coeff, double* df) const;
  void Project(int daughter, const double* df, double* coeff) const;

  /** The entropy weight f(1 -/+ f) on that daughter's grid -- the basis weight, the
   *  metric of the H-theorem inner product, and (times {1, q}) the equilibrium
   *  perturbations themselves. Zero wherever the daughter is unpopulated, which is
   *  most of the grid early on. */
  const std::vector<double>& weight(int daughter) const {
    return basis_[daughter].w;
  }
  double q_ref(int daughter) const {
    return basis_[daughter].q_ref;
  }

  /** What SetBackground had to do to the supplied occupations to make a weight out of
   *  them. `n_clamped` counts bins where the value used differs from the one passed --
   *  a fermion past saturation, or a negative occupation from either statistic --
   *  and `f_max` is the peak occupation seen.
   *
   *  This exists because the clamp is otherwise silent, and what it clamps is a defect
   *  in the BACKGROUND solve, not in this class: a fermion bin at f > 1 means the
   *  discrete gather has broken Pauli blocking. Clamping keeps the reduction
   *  well-posed; reporting is what stops that standing in for a fix. A healthy
   *  background must give n_clamped = 0, or the report carries no information. */
  struct WeightClamp {
    int n_clamped = 0;
    double f_max  = 0.;
  };
  WeightClamp weight_clamp(int daughter) const {
    return basis_[daughter].clamp;
  }
  /** The frozen basis, for the species that has to CARRY these moments as its state.
   *
   *  The reduced daughter's perturbation variables are m_j = Sum dq q^2 psi_j df, so the
   *  species needs the same psi_j this operator projects with -- not a second
   *  construction of them. One Stieltjes recurrence, one set of polynomials, both
   *  consumers reading it: a duplicated basis that drifted by a rounding would break the
   *  conservation identities in a way no test here would catch.
   *
   *  psi is n_moments * n_grid, row-major in the moment index. */
  const std::vector<double>& psi(int daughter) const {
    return basis_[daughter].psi;
  }
  const std::vector<double>& gnorm(int daughter) const {
    return basis_[daughter].gnorm;
  }
  /** psi_1 = t - alpha0, so the ENERGY functional is q_ref*(m_1 + alpha0*m_0). Every
   *  consumer that reads energy out of the moments needs this, and getting it wrong is
   *  silent: the hierarchy still runs, it just conserves the wrong thing. */
  double alpha0(int daughter) const {
    return basis_[daughter].alpha0;
  }

 private:
  /** One daughter's basis.
   *
   *  The test functions psi_j are ORTHOGONAL polynomials in the measure dq q^2 w(q),
   *  built by the Stieltjes three-term recurrence rather than as monomials
   *  (q/q_ref)^j. Monomials are what the first version used, and they run out of
   *  conditioning around n_moments ~ 8 -- the Gram matrix is a Hankel moment matrix
   *  and the Cholesky fails outright. The recurrence never forms t^j and stays well
   *  conditioned, so the DOF axis extends as far as the daughter grid supports. It
   *  also makes the Gram DIAGONAL, so reconstruction is a scaling, not a solve.
   *
   *  psi_0 = 1 and psi_1 = t - alpha0 by construction, so span{psi_0, psi_1} =
   *  span{1, t}: number and energy stay exactly representable, which is what makes
   *  conservation survive the reduction. Their COORDINATES change -- energy is
   *  q_ref*(m_1 + alpha0*m_0), not q_ref*m_1 -- and Conservation() carries that. */
  struct Basis {
    std::vector<double> w;       ///< the entropy weight f(1 -/+ f) on that daughter's grid
    std::vector<double> psi;     ///< n_moments * n_grid, the orthogonal test functions
    std::vector<double> gnorm;   ///< n_moments, <phi_j, phi_j> -- the Gram is diagonal
    std::vector<double> metric;  ///< n_moments^2, the induced entropy metric
    double q_ref  = 1.;
    double alpha0 = 0.;  ///< psi_1 = t - alpha0, so t = psi_1 + alpha0 psi_0
    WeightClamp clamp;   ///< what SetBackground had to clamp to build w
  };

  /** Metric of the ENTROPY inner product on the reduced state, as a dense SPD matrix:
   *  diag(g_H dq q^2 / w_H) on the parent, and g_i * G^-1 on each daughter block.
   *  G is the (diagonal) Gram of the orthogonalised basis, so the daughter blocks are
   *  a scaling. This is the form the entropy weight induces, and it is what makes P
   *  and R adjoint -- i.e. what makes the reduction a congruence.
   *
   *  The degeneracy weights (2, 2, 1) are the kernel's own, folded so the conservation
   *  identities read as they do in decay_kernel_test. */
  void BuildMetric(std::vector<double>& S) const;
  void BuildMetricExact(std::vector<double>& S) const;

  const DecayTransitionKernel* kernel_;
  GridView parent_, fermion_, boson_;
  Statistics fermion_stat_, boson_stat_;
  Config cfg_;

  double a_ = 0., m_ = 0.;
  std::vector<double> eps_;  ///< parent sqrt(q^2 + a^2 m^2)
  std::vector<double> w_H_;  ///< parent entropy weight, for the metric only
  Basis basis_[2];           ///< 0 = fermion, 1 = boson
};

/**
 * M~_l tabulated over background rows, applied as a dense mat-vec.
 *
 * The reduced operator is k-INDEPENDENT: a functional of the background alone.
 * Applying it matrix-free -- reconstruct onto the quadrature grid, sweep the kernel,
 * project back -- therefore repeats k-independent work at every k-mode and every
 * step, and leaves the per-RHS cost scaling with dr_N_q just as the exact scheme's
 * did. Assembling it once per background row costs size() kernel applications per row
 * and turns the k-loop's collision into an S x S mat-vec, S = n_parent + 2 n_moments.
 *
 * That is what takes dr_N_q off the hot path, and dr_N_q is where the error actually
 * is: the internal quadrature grid moves P(k) more than the moment count does, by
 * well over an order of magnitude. Tabulated, refining it costs only the precompute,
 * which amortises over the whole k-loop.
 *
 * WHAT IT APPROXIMATES. The assembly is exact (reduced_operator_test pins M~ z against
 * P M R z at machine precision), so the only new approximation is interpolation in
 * ln a between rows -- and the rows are the background table's own, which are dense,
 * on an operator that is smooth in ln a.
 *
 * OUTSIDE THE TABULATED RANGE THE OPERATOR IS TAKEN AS ZERO. The table covers the
 * collision window, and a caller that picks the window too narrow gets a
 * DISCONTINUITY, not a slow degradation. Choose it on a norm threshold, never by eye.
 */
class ReducedOperatorTable {
 public:
  /** Discards any existing table and sets the block geometry. */
  void Reset(int size, int l_max);

  /** Append one background row. Rows must arrive in ASCENDING a -- the lookup is a
   *  binary search and does not sort. `M` is the (l_max+1)*S*S row-major block that
   *  ReducedCollisionOperator::Assemble produces, and is copied. */
  void AddRow(double a, const std::vector<double>& M);

  int size() const {
    return size_;
  }
  int l_max() const {
    return l_max_;
  }
  int n_rows() const {
    return static_cast<int>(a_.size());
  }
  bool empty() const {
    return a_.empty();
  }
  /** Bytes held, for the precompute's own report -- this is the resource the design
   *  trades the k-loop's time for, so it should never be silent. */
  size_t bytes() const {
    return M_.size() * sizeof(double);
  }

  /** dz += M~_l(a) z, with M~ linear in ln a between bracketing rows.
   *
   *  ACCUMULATES rather than assigns, because the collision is one term of a
   *  right-hand side that free streaming has already written. Outside the tabulated
   *  range it adds nothing at all. */
  void Apply(double a, int l, const double* z, double* dz) const;

  /** The interpolated operator's DIAGONAL, written (not accumulated) into `out[size]`.
   *
   *  The exponential evolver integrates diag(M) exactly and treats the rest explicitly,
   *  so it needs this every RHS.  Reading it off the table matters for more than tidiness:
   *  taking it from the kernel instead would force a PrepareTransitions per step, which is
   *  most of the work the tabulation exists to remove.
   *
   *  It is also a better diagonal than the matrix-free path can offer.  There the
   *  daughters have no per-moment diagonal at all -- the collision does not act diagonally
   *  on them -- so each moment is handed the Rayleigh quotient of the per-bin diagonal, a
   *  defensible integrating-factor choice but a guess.  Here it is the actual diagonal of
   *  the actual operator the ODE sees.  Zero outside the tabulated range. */
  void Diagonal(double a, int l, double* out) const;

 private:
  /** Index of the last row with a_[i] <= a, or -1 if a is below the first row.
   *
   *  O(1) when the rows are uniform in ln a, which is the case that actually occurs:
   *  CLASS's background table is uniform in ln a (measured, constant ratio 4.0566 per
   *  200 rows to five digits) and striding preserves that, so the index is one
   *  subtraction and one multiply.  Uniformity is DETECTED in AddRow rather than
   *  assumed -- a non-uniform table still works, just via the binary search. */
  int LowerRow(double a) const;

  int size_  = 0;
  int l_max_ = -1;
  std::vector<double> a_;     ///< ascending
  std::vector<double> ln_a_;  ///< cached: the interpolation variable
  /** Set while the ln-a spacing is uniform to within kUniformTol; cleared by the first
   *  row that is not.  inv_dlna_ is only meaningful while it holds. */
  bool uniform_ln_a_                  = true;
  double inv_dlna_                    = 0.;
  static constexpr double kUniformTol = 1e-9;
  std::vector<double> M_;  ///< n_rows * (l_max+1) * size^2, row-major throughout
};
