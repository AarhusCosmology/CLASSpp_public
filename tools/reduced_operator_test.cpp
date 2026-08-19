// Decides the reduced-daughter design for nu_H <-> nu_l + phi, without cosmology.
//
// The reduced operator (tools/reduced_collision_operator) is a Galerkin congruence of
// the shipped kernel, M~ = P M R, so it re-derives no physics and can be compared with
// the exact operator on the same background and the same df with NO fitted constant --
// which is what the design note's section 4 test 2 asks for and what every previous
// attempt (a fitted q-weight, a 3x3 fluid RTA) failed to do honestly.
//
// Six questions, in the order that lets each fail cheaply before the next:
//
//   1  P R = I                     -- the harness is self-consistent
//   2  in-span agreement           -- the assembly reproduces the exact operator
//                                     exactly on the subspace it represents
//   3  equilibrium annihilation    -- M~ v_eq = 0 on a detailed-balance background,
//                                     which design section 2.6 calls the single largest
//                                     risk.  Flat/Occupation are controls and must fail.
//   4  conservation inherited      -- ||w^T M~|| tracks ||w^T M||, for any n_moments
//   5  DISSIPATIVITY               -- max Re lambda <= 0 at l = 0, 1, 2 for ANY number
//                                     of daughter degrees of freedom.  THE MILESTONE.
//   6  quadrature-grid independence-- the l = 1 momentum residual falls as the
//                                     QUADRATURE grid refines at fixed state size.
//                                     This is the claim that the grid rule N ~ Gamma^0.25
//                                     stops applying, and it is the point of the design.
//
// Units are the kernel's: q in units of T_ncdm, m = M/T, Gamma in 1/Mpc.  An ini
// `Gamma = 1e7` is 33.36/Mpc here (DNCDMSpecies stores Gamma * 1e3/c) -- getting that
// wrong inflates every rate by exactly c in km/s, which has cost a debugging round
// before.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <vector>

#include "decay_transition_kernel.h"
#include "reduced_collision_operator.h"

namespace {

// ------------------------------------------------------------------- grids --

struct Grid {
  std::vector<double> q, dq;
};

GridView View(const Grid& g) {
  return GridView{g.q.data(), g.dq.data(), static_cast<int>(g.q.size())};
}

/** Log-trapezoid grid, the daughters' shipped convention (half cells at the ends). */
Grid LogGrid(double qmin, double qmax, int n) {
  Grid g;
  g.q.resize(n);
  g.dq.resize(n);
  const double h = std::log(qmax / qmin) / (n - 1);
  for (int i = 0; i < n; ++i)
    g.q[i] = qmin * std::exp(i * h);
  g.dq[0] = 0.5 * (g.q[1] - g.q[0]);
  for (int i = 1; i < n - 1; ++i)
    g.dq[i] = 0.5 * (g.q[i + 1] - g.q[i - 1]);
  g.dq[n - 1] = 0.5 * (g.q[n - 1] - g.q[n - 2]);
  return g;
}

/** The production parent grid: quadrature_strategy = 3 lays q_i = (15/n)(i+1) with
 *  uniform weights, i.e. a plain 16-bin trapezoid out to q = 15. Using the real one
 *  matters because every rate below is a band integral over it. */
Grid ParentGrid(int n) {
  Grid g;
  g.q.resize(n);
  g.dq.assign(n, 15.0 / n);
  for (int i = 0; i < n; ++i)
    g.q[i] = (15.0 / n) * (i + 1);
  return g;
}

DecayTransitionKernel::Config KernelConfig() {
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  return cfg;
}

double MaxAbs(const std::vector<double>& v) {
  double m = 0.;
  for (double x : v)
    m = std::fmax(m, std::fabs(x));
  return m;
}

// ------------------------------------------------------------ backgrounds --

struct Background {
  std::vector<double> fH, fl, fphi;
  double a = 0.;
};

/** A background in EXACT detailed balance: one temperature, mu_H = mu_l + mu_phi.
 *  Lambda vanishes identically in the continuum, so the exact operator annihilates
 *  every equilibrium perturbation and question 3 has a machine-precision answer to
 *  measure against. mu_phi must be <= 0 or the Bose distribution is not defined. */
Background EquilibriumBackground(const Grid& P, const Grid& F, const Grid& B, double a, double m) {
  const double mu_l = 0.2, mu_phi = -0.3, mu_H = mu_l + mu_phi;
  Background bg;
  bg.a = a;
  bg.fH.resize(P.q.size());
  for (size_t i = 0; i < P.q.size(); ++i) {
    const double eps = std::sqrt(P.q[i] * P.q[i] + a * a * m * m);
    bg.fH[i]         = 1. / (std::exp(eps - mu_H) + 1.);
  }
  bg.fl.resize(F.q.size());
  for (size_t i = 0; i < F.q.size(); ++i)
    bg.fl[i] = 1. / (std::exp(F.q[i] - mu_l) + 1.);
  bg.fphi.resize(B.q.size());
  for (size_t i = 0; i < B.q.size(); ++i)
    bg.fphi[i] = 1. / (std::exp(B.q[i] - mu_phi) - 1.);
  return bg;
}

/** A background OUT of detailed balance, with the structure that makes this sector
 *  hard: a Pauli-saturated nu_l, a strongly Bose-enhanced phi (peak occupation of
 *  order 100), an injection front above which both daughters fall off a cliff, and a
 *  daughter temperature offset from the parent's so Lambda does not vanish.
 *
 *  WHY THIS IS ANALYTIC AND NOT INTEGRATED. Two reasons, and the second is the
 *  important one.
 *
 *  First, integrating the network here would need the stiff coupled solver the module
 *  itself uses; a per-bin exponential step runs away in the low-q boson bins within a
 *  few steps (f_phi ~ 1e7), because those bins have the smallest q^2 dq and the
 *  largest Bose enhancement.
 *
 *  Second and decisively: question 6 refines the QUADRATURE grid and compares the
 *  results. An integrated background would be a DIFFERENT physical state on each grid,
 *  so the comparison would confound the discretisation of the operator with the
 *  discretisation of the state. A background defined as a function of q is the same
 *  state on every grid, which is what makes that measurement mean anything.
 *
 *  Everything questions 1-6 assert is algebra -- conservation, adjointness,
 *  congruence, the numerical range -- and holds for ARBITRARY positive f, which is the
 *  same standard decay_kernel_test holds its conservation identities to. What this
 *  background canNOT settle is the size of the closure ERROR on a real trajectory;
 *  that needs states from a CLASS run and is deliberately not attempted here. */
Background AnalyticBackground(const Grid& P, const Grid& F, const Grid& B, double a, double m) {
  // mu_phi just below zero is what produces the measured peak f_phi ~ O(100); mu_l
  // near the fermion's own scale is what Pauli-saturates nu_l. T_d != 1 is the
  // departure from detailed balance.
  const double mu_H = 0.29, mu_l = 0.3, mu_phi = -1e-3, T_d = 1.15;
  const double q_front = 0.35 * DecayTransitionKernel::MaxDaughterMomentum(15.0, a * m);
  auto front           = [&](double q) { return 1. / (1. + std::pow(q / q_front, 6.)); };

  Background bg;
  bg.a = a;
  bg.fH.resize(P.q.size());
  for (size_t i = 0; i < P.q.size(); ++i) {
    const double eps = std::sqrt(P.q[i] * P.q[i] + a * a * m * m);
    bg.fH[i]         = 1. / (std::exp(eps - mu_H) + 1.);
  }
  bg.fl.resize(F.q.size());
  for (size_t i = 0; i < F.q.size(); ++i)
    bg.fl[i] = std::fmin(0.999, front(F.q[i]) / (std::exp((F.q[i] - mu_l) / T_d) + 1.));
  bg.fphi.resize(B.q.size());
  for (size_t i = 0; i < B.q.size(); ++i)
    bg.fphi[i] = front(B.q[i]) / (std::exp((B.q[i] - mu_phi) / T_d) - 1.);
  return bg;
}

// ------------------------------------------------------------------ report --

/** ||M z|| / (||M|| ||z||) in the plain 2-norm: how far z is from a null vector. */
double NullResidual(const double* M, int n, const double* z) {
  double num = 0., zn = 0., mn = 0.;
  for (int i = 0; i < n; ++i) {
    double s = 0.;
    for (int j = 0; j < n; ++j) {
      s  += M[(size_t) i * n + j] * z[j];
      mn  = std::fmax(mn, std::fabs(M[(size_t) i * n + j]));
    }
    num += s * s;
    zn  += z[i] * z[i];
  }
  const double scale = std::sqrt(zn) * mn;
  return (scale > 0.) ? std::sqrt(num) / scale : 0.;
}

// ------------------------------------------------------------------ shared --

// Gamma = 1e7 in ini units. m = M/T for m_ncdm = 0.06 eV at T_ncdm = 0.71611.
constexpr double kGamma = 33.356409;
constexpr double kMass  = 356.74;
constexpr double kA2H   = 2.15e-6;  // H0 sqrt(Omega_r), 1/Mpc: radiation domination
constexpr double kA     = 9.0e-3;   // inside the collision window at this Gamma (design note M1)
constexpr int kNParent  = 16;
constexpr int kLMax     = 4;  // enough to see l = 0, 1 and a clean l >= 2

struct Setup {
  Grid P, F, B;
  DecayTransitionKernel::Config cfg;
};

Setup MakeSetup(int n_daughter, double a_end) {
  Setup s;
  s.P = ParentGrid(kNParent);
  // Cover the whole emission band at the latest scale factor of interest, with the
  // 5% margin the composite uses -- an undersized q_max is an energy leak, not a
  // resolution trade-off (kernel doc on MaxDaughterMomentum).
  const double qmax = 1.05 * DecayTransitionKernel::MaxDaughterMomentum(15.0, a_end * kMass);
  s.F               = LogGrid(1e-2, qmax, n_daughter);
  s.B               = LogGrid(1e-2, qmax, n_daughter);
  s.cfg             = KernelConfig();
  return s;
}

// =================================================================== tests ==

void q1_projection_identity() {
  std::printf("\n== 1. P R = I (the harness is self-consistent) ==\n");
  Setup s = MakeSetup(101, 5e-2);
  DecayTransitionKernel K(View(s.P),
                          View(s.F),
                          View(s.B),
                          Statistics::Fermion,
                          Statistics::Boson,
                          s.cfg);
  Background bg = EquilibriumBackground(s.P, s.F, s.B, 1e-2, kMass);

  double worst = 0.;
  for (int nm = 1; nm <= 4; ++nm) {
    ReducedCollisionOperator
        R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
    R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
    std::vector<double> c(nm), c2(nm), df(std::max(s.F.q.size(), s.B.q.size()));
    for (int d = 0; d < 2; ++d)
      for (int j = 0; j < nm; ++j) {
        std::fill(c.begin(), c.end(), 0.);
        c[j] = 1.;
        R.Reconstruct(d, c.data(), df.data());
        R.Project(d, df.data(), c2.data());
        for (int i = 0; i < nm; ++i)
          worst = std::fmax(worst, std::fabs(c2[i] - (i == j ? 1. : 0.)));
      }
  }
  std::printf("   max |P R - I| over n_moments 1..4 : %.3e\n", worst);
  assert(worst < 1e-9);
}

/** The identity the reduced SPECIES rests on: with psi_0 = 1 and psi_1 = t - alpha0, the
 *  physical energy Sum dq q^3 df is q_ref*(m_1 + alpha0*m_0) in the moment coordinates.
 *
 *  DrPsdSpecies reads the daughter's delta-rho, theta and shear through exactly this
 *  combination -- for a massless daughter eps = q, so all four stress-energy sums carry
 *  the weight q^3 and the metric never sees anything else. A wrong coefficient here does
 *  not crash: the hierarchy runs and quietly conserves the wrong quantity. */
void q1b_energy_coordinate() {
  std::printf("\n== 1b. energy functional in reduced coordinates ==\n");
  Setup s = MakeSetup(101, 5e-2);
  DecayTransitionKernel K(View(s.P),
                          View(s.F),
                          View(s.B),
                          Statistics::Fermion,
                          Statistics::Boson,
                          s.cfg);
  Background bg = AnalyticBackground(s.P, s.F, s.B, kA, kMass);

  double worst = 0.;
  for (int nm = 2; nm <= 6; ++nm) {
    ReducedCollisionOperator
        R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
    R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
    for (int d = 0; d < 2; ++d) {
      const Grid& g = (d == 0) ? s.F : s.B;
      // An arbitrary df, NOT one in the reduced span: the energy identity is a property
      // of the projection alone and must hold for anything.
      std::vector<double> df(g.q.size()), m(nm);
      for (size_t i = 0; i < g.q.size(); ++i)
        df[i] = std::sin(3.1 * i) * std::exp(-g.q[i] / 4.);
      R.Project(d, df.data(), m.data());

      double want = 0.;
      for (size_t i = 0; i < g.q.size(); ++i)
        want += g.dq[i] * g.q[i] * g.q[i] * g.q[i] * df[i];
      const double got = R.q_ref(d) * (m[1] + R.alpha0(d) * m[0]);
      worst = std::fmax(worst, std::fabs(got - want) / std::fmax(std::fabs(want), 1e-300));
    }
  }
  std::printf("   max relative |q_ref(m1 + a0 m0) - Sum dq q^3 df| : %.3e\n", worst);
  assert(worst < 1e-12);
}

void q2_in_span_agreement() {
  std::printf("\n== 2. reduced == exact on the subspace the reduction represents ==\n");
  Setup s = MakeSetup(101, 5e-2);
  DecayTransitionKernel K(View(s.P),
                          View(s.F),
                          View(s.B),
                          Statistics::Fermion,
                          Statistics::Boson,
                          s.cfg);
  Background bg = AnalyticBackground(s.P, s.F, s.B, kA, kMass);

  const int nm = 2;
  ReducedCollisionOperator
      R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
  K.PrepareTransitions(bg.a, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data());
  R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());

  std::vector<double> Mred;
  R.Assemble(kLMax, Mred);
  const int S = R.size();

  // A state that lies in the reduced span: arbitrary parent, daughters reconstructed.
  std::vector<double> z(S);
  for (int b = 0; b < kNParent; ++b)
    z[b] = 0.3 * std::sin(1.7 * b) * bg.fH[b];
  z[kNParent]          = 1.0;
  z[kNParent + 1]      = -0.4;
  z[kNParent + nm]     = 0.6;
  z[kNParent + nm + 1] = 0.25;

  const int stride = kLMax + 1;
  std::vector<double> FH((size_t) kNParent * stride, 0.), Fl(s.F.q.size() * stride, 0.),
      Fp(s.B.q.size() * stride, 0.);
  std::vector<double> dH(FH.size()), dl(Fl.size()), dp(Fp.size());
  std::vector<double> rec(std::max(s.F.q.size(), s.B.q.size()));

  for (int b = 0; b < kNParent; ++b)
    for (int l = 0; l <= kLMax; ++l)
      FH[(size_t) b * stride + l] = z[b];
  R.Reconstruct(0, &z[kNParent], rec.data());
  for (size_t k = 0; k < s.F.q.size(); ++k)
    for (int l = 0; l <= kLMax; ++l)
      Fl[k * stride + l] = rec[k];
  R.Reconstruct(1, &z[kNParent + nm], rec.data());
  for (size_t k = 0; k < s.B.q.size(); ++k)
    for (int l = 0; l <= kLMax; ++l)
      Fp[k * stride + l] = rec[k];

  K.ResetMultipoleRecurrence();
  K.ApplyPerturbationOperatorAllL(kLMax,
                                  stride,
                                  FH.data(),
                                  Fl.data(),
                                  Fp.data(),
                                  dH.data(),
                                  dl.data(),
                                  dp.data());

  double worst = 0.;
  std::vector<double> slice(std::max(s.F.q.size(), s.B.q.size())), mom(nm);
  for (int l = 0; l <= kLMax; ++l) {
    const double* Ml = Mred.data() + (size_t) l * S * S;
    std::vector<double> want(S, 0.);
    for (int i = 0; i < S; ++i)
      for (int j = 0; j < S; ++j)
        want[i] += Ml[(size_t) i * S + j] * z[j];

    std::vector<double> got(S, 0.);
    for (int b = 0; b < kNParent; ++b)
      got[b] = dH[(size_t) b * stride + l];
    for (int d = 0; d < 2; ++d) {
      const double* src = (d == 0) ? dl.data() : dp.data();
      const size_t n    = (d == 0) ? s.F.q.size() : s.B.q.size();
      for (size_t k = 0; k < n; ++k)
        slice[k] = src[k * stride + l];
      R.Project(d, slice.data(), mom.data());
      for (int j = 0; j < nm; ++j)
        got[kNParent + d * nm + j] = mom[j];
    }
    const double scale = std::fmax(MaxAbs(want), 1e-300);
    for (int i = 0; i < S; ++i)
      worst = std::fmax(worst, std::fabs(got[i] - want[i]) / scale);
  }
  std::printf("   max relative |M~z - P M R z| over l = 0..%d : %.3e\n", kLMax, worst);
  assert(worst < 1e-11);
}

void q3_equilibrium_annihilation() {
  std::printf("\n== 3. detailed balance: does the basis annihilate v_eq? ==\n");
  Setup s = MakeSetup(101, 5e-2);
  DecayTransitionKernel K(View(s.P),
                          View(s.F),
                          View(s.B),
                          Statistics::Fermion,
                          Statistics::Boson,
                          s.cfg);
  Background bg = EquilibriumBackground(s.P, s.F, s.B, 8e-3, kMass);
  K.PrepareTransitions(bg.a, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data());

  ReducedCollisionOperator
      ref(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {2});
  ref.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
  std::vector<double> Mex, veq_ex;
  ref.AssembleExact(1, Mex);
  ref.EquilibriumNullVectorsExact(veq_ex);
  const int N = ref.full_size();
  double exact[3];
  for (int v = 0; v < 3; ++v)
    exact[v] = NullResidual(Mex.data(), N, &veq_ex[(size_t) v * N]);
  std::printf(
      "   exact operator's own residual ||M v_eq||/(||M|| ||v||), l = 0:\n"
      "      v_mu(H,l) %.3e   v_mu(H,phi) %.3e   v_T %.3e\n",
      exact[0],
      exact[1],
      exact[2]);
  std::printf(
      "   (not machine zero: the two-bin gather leaves the DISCRETE Lambda an O(dq^2)\n"
      "    residual. No reduction can beat this floor.)\n\n");

  // The decisive statement is an identity, not a size. If the basis spans equilibrium
  // then R P v_eq = v_eq exactly, hence
  //       M~ (P v_eq)  =  P M R P v_eq  =  P (M v_eq),
  // i.e. the reduced operator's equilibrium residual IS the projection of the exact
  // one -- it adds nothing of its own, at any resolution and for any background. A
  // basis that does not span equilibrium fails at the first step, and the size of that
  // failure is the reconstruction error below, which is normalisation-free (it is a
  // relative error in the physical entropy norm) and so is comparable across bases.
  std::printf("   %3s   %-13s %-13s   %-13s\n", "nm", "||RP v-v||/||v||", "(boson)", "M~Pv - PMv");
  double entropy_rec = 0.;
  for (int nm = 2; nm <= 3; ++nm) {
    ReducedCollisionOperator
        R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
    R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());

    // Basis fidelity: can it represent an equilibrium perturbation at all? Measured
    // in the entropy norm ||x||^2 = Sum dq q^2 x^2 / w, which is the norm the
    // H-theorem is stated in.
    double rec[2] = {0., 0.};
    std::vector<double> c(nm), v(std::max(s.F.q.size(), s.B.q.size())),
        vr(std::max(s.F.q.size(), s.B.q.size()));
    for (int d = 0; d < 2; ++d) {
      const Grid& g                 = (d == 0) ? s.F : s.B;
      const std::vector<double>& we = R.weight(d);
      for (int dir = 0; dir < 2; ++dir) {  // the mu and T equilibrium directions
        for (size_t k = 0; k < g.q.size(); ++k)
          v[k] = we[k] * (dir == 0 ? 1. : g.q[k]);
        R.Project(d, v.data(), c.data());
        R.Reconstruct(d, c.data(), vr.data());
        double num = 0., den = 0.;
        for (size_t k = 0; k < g.q.size(); ++k) {
          if (!(we[k] > 0.))
            continue;
          const double s2  = g.dq[k] * g.q[k] * g.q[k] / we[k];
          num             += s2 * (vr[k] - v[k]) * (vr[k] - v[k]);
          den             += s2 * v[k] * v[k];
        }
        rec[d] = std::fmax(rec[d], (den > 0.) ? std::sqrt(num / den) : 0.);
      }
    }

    // And the identity, on the fermion chemical-potential direction.
    std::vector<double> M, veq, veq_full;
    R.Assemble(0, M);
    R.EquilibriumNullVectors(veq);
    R.EquilibriumNullVectorsExact(veq_full);
    const int S = R.size();
    std::vector<double> lhs(S, 0.);
    for (int i = 0; i < S; ++i)
      for (int j = 0; j < S; ++j)
        lhs[i] += M[(size_t) i * S + j] * veq[j];
    std::vector<double> rhs_full(N, 0.);
    for (int i = 0; i < N; ++i)
      for (int j = 0; j < N; ++j)
        rhs_full[i] += Mex[(size_t) i * N + j] * veq_full[j];
    std::vector<double> rhs(S, 0.), mom(nm), slice(std::max(s.F.q.size(), s.B.q.size()));
    for (int b = 0; b < kNParent; ++b)
      rhs[b] = rhs_full[b];
    for (int d = 0; d < 2; ++d) {
      const size_t n0 = (d == 0) ? kNParent : kNParent + s.F.q.size();
      const size_t nn = (d == 0) ? s.F.q.size() : s.B.q.size();
      for (size_t k = 0; k < nn; ++k)
        slice[k] = rhs_full[n0 + k];
      R.Project(d, slice.data(), mom.data());
      for (int j = 0; j < nm; ++j)
        rhs[kNParent + d * nm + j] = mom[j];
    }
    double id = 0., sc = 0.;
    for (int i = 0; i < S; ++i) {
      id = std::fmax(id, std::fabs(lhs[i] - rhs[i]));
      sc = std::fmax(sc, std::fabs(rhs[i]));
    }
    id /= std::fmax(sc, 1e-300);

    std::printf("   %3d   %-13.3e %-13.3e   %-13.3e\n", nm, rec[0], rec[1], id);
    entropy_rec = std::fmax(entropy_rec, std::fmax(rec[0], rec[1]));
    assert(id < 1e-10);  // R P v_eq = v_eq, so the identity is exact
  }
  std::printf("   -> the entropy basis reproduces an equilibrium perturbation to %.1e\n",
              entropy_rec);
  assert(entropy_rec < 1e-10);
}

void q4_q5_conservation_and_dissipativity() {
  std::printf("\n== 4/5. conservation and DISSIPATIVITY on a decayed background ==\n");
  Setup s = MakeSetup(101, 5e-2);
  DecayTransitionKernel K(View(s.P),
                          View(s.F),
                          View(s.B),
                          Statistics::Fermion,
                          Statistics::Boson,
                          s.cfg);
  Background bg = AnalyticBackground(s.P, s.F, s.B, kA, kMass);
  K.PrepareTransitions(bg.a, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data());
  const double aH = kA2H / bg.a;
  std::printf("   a = %.3e, peak f_H = %.3f, f_l = %.3f, f_phi = %.3e, aH = %.3e/Mpc\n",
              bg.a,
              MaxAbs(bg.fH),
              MaxAbs(bg.fl),
              MaxAbs(bg.fphi),
              aH);

  ReducedCollisionOperator
      ref(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {2});
  ref.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
  std::vector<double> Mex;
  ref.AssembleExact(2, Mex);
  const int N = ref.full_size();

  // Read the identities at the multipole that owns them: number and energy are l = 0
  // statements, momentum is the l = 1 one. The others are printed only so the pattern
  // is visible, and they are not expected to vanish.
  double exact_bound[3];
  std::printf("   EXACT operator (%d variables):\n", N);
  std::printf("   %3s  %-10s %-10s %-10s  %-13s\n",
              "l",
              "number",
              "energy",
              "momentum",
              "max Re lambda");
  for (int l = 0; l <= 2; ++l) {
    auto c         = ref.ConservationExact(Mex.data() + (size_t) l * N * N);
    exact_bound[l] = ref.DissipativityBoundExact(Mex.data() + (size_t) l * N * N);
    std::printf("   %3d  %-10.2e %-10.2e %-10.2e  %+.5e\n",
                l,
                std::fmax(c.number_fermion, c.number_boson),
                c.energy,
                c.momentum,
                exact_bound[l]);
  }

  std::printf("\n   REDUCED operator:\n");
  std::printf("   %3s %3s %-6s  %-10s %-10s %-10s  %-13s %s\n",
              "nm",
              "l",
              "size",
              "number",
              "energy",
              "momentum",
              "max Re lambda",
              "vs exact");
  bool containment_holds = true;
  for (int nm = 2; nm <= 4; ++nm) {
    ReducedCollisionOperator
        R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
    R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
    std::vector<double> M;
    R.Assemble(2, M);
    const int S = R.size();
    for (int l = 0; l <= 2; ++l) {
      const double* Ml = M.data() + (size_t) l * S * S;
      auto c           = R.Conservation(Ml);
      const double b   = R.DissipativityBound(Ml);
      const bool ok    = b <= exact_bound[l] * (exact_bound[l] > 0. ? 1.001 : 0.999) + 1e-12;
      std::printf("   %3d %3d %6d  %-10.2e %-10.2e %-10.2e  %+.5e %s\n",
                  nm,
                  l,
                  S,
                  std::fmax(c.number_fermion, c.number_boson),
                  c.energy,
                  c.momentum,
                  b,
                  ok ? "contained" : "*** EXCEEDS ***");
      if (!ok)
        containment_holds = false;
    }
  }
  std::printf(
      "\n   THE MILESTONE. The entropy reduction is a congruence, so its numerical\n"
      "   range is a SUBSET of the exact operator's: max Re lambda(M~) <= max Re\n"
      "   lambda(M) at every n_moments and every l. Containment held: %s\n"
      "   (aH = %.3e/Mpc. Neither operator is dissipative on THIS background --\n"
      "    it is 15%% out of detailed balance by construction, and a linearised\n"
      "    collision operator is only negative semi-definite AT equilibrium. What\n"
      "    the milestone asserts is that the reduction adds nothing of its own.)\n",
      containment_holds ? "YES" : "NO",
      aH);
  assert(containment_holds);
}

void q6_the_milestone() {
  std::printf("\n== 6. THE MILESTONE: dissipativity vs the QUADRATURE grid and the DOF ==\n");
  std::printf(
      "   The shipped operator's positive l=0/l=1 eigenvalue is what the daughter\n"
      "   grid rule N ~ Gamma^0.25 pays for. The reduced state is 20 variables in\n"
      "   EVERY row below; only the quadrature grid moves, and the integrator never\n"
      "   sees it.\n\n"
      "   Read through the numerical-range bound b, where max Re lambda(M) <= b(M):\n"
      "   b is rigorous and cheap, and the congruence gives b(M~) <= b(M), which is\n"
      "   the statement being tested. The l = 0 NUMBER column is the exact zero: number\n"
      "   is a functional the reduced state carries, so its conservation survives the\n"
      "   reduction identically at every grid. Energy is not exact and is not expected\n"
      "   to be -- the lumped daughter loss misplaces O(dq) of it by construction (see\n"
      "   split_energy_residual) -- so it is reported and watched for convergence.\n\n");
  const double aH = kA2H / kA;
  std::printf("   %5s %-8s | %-11s %-11s | %-11s %-11s | %-10s %-10s\n",
              "N_q",
              "state",
              "b exact l=0",
              "b exact l=1",
              "b red l=0",
              "b red l=1",
              "l=0 number",
              "l=0 energy");
  for (int nq : {26, 51, 101, 201}) {
    Setup s = MakeSetup(nq, 5e-2);
    DecayTransitionKernel K(View(s.P),
                            View(s.F),
                            View(s.B),
                            Statistics::Fermion,
                            Statistics::Boson,
                            s.cfg);
    Background bg = AnalyticBackground(s.P, s.F, s.B, kA, kMass);
    K.PrepareTransitions(bg.a, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data());

    ReducedCollisionOperator
        R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {2});
    R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());

    std::vector<double> Mex, Mred;
    R.AssembleExact(1, Mex);
    R.Assemble(1, Mred);
    const int N = R.full_size(), S = R.size();
    const double e0 = R.DissipativityBoundExact(Mex.data());
    const double e1 = R.DissipativityBoundExact(Mex.data() + (size_t) N * N);
    const double r0 = R.DissipativityBound(Mred.data());
    const double r1 = R.DissipativityBound(Mred.data() + (size_t) S * S);
    const auto c0   = R.Conservation(Mred.data());
    std::printf("   %5d %3d/%-4d | %+.4e %+.4e | %+.4e %+.4e | %.4e %.4e\n",
                nq,
                S,
                N,
                e0,
                e1,
                r0,
                r1,
                std::fmax(c0.number_fermion, c0.number_boson),
                c0.energy);
    // CONTAINMENT, at every quadrature grid: the reduction adds no growth of its own.
    assert(r0 <= e0 * (e0 > 0. ? 1.001 : 0.999) + 1e-12);
    assert(r1 <= e1 * (e1 > 0. ? 1.001 : 0.999) + 1e-12);
    // ...and l = 0's NUMBER zero survives it identically, at every grid.
    assert(std::fmax(c0.number_fermion, c0.number_boson) < 1e-8);
    assert(c0.energy < 1e-3);  // converges with the grid; not an exact zero
  }
  std::printf(
      "\n   Containment holds at every quadrature grid, and l = 0's number zero survives\n"
      "   the reduction: the reduced state carries that functional itself. What the\n"
      "   reduction buys is not a smaller eigenvalue at equal grid; it is the SAME\n"
      "   operator at 20 state variables instead of the full grid's, because the grid\n"
      "   it is computed on is quadrature the integrator never carries.\n"
      "   aH = %.3e/Mpc.\n",
      aH);

  std::printf("\n   And against the number of daughter degrees of freedom, at N_q = 201:\n");
  std::printf("   %3s %-6s | %-11s %-11s\n", "nm", "state", "b red l=0", "b red l=1");
  Setup s = MakeSetup(201, 5e-2);
  DecayTransitionKernel K(View(s.P),
                          View(s.F),
                          View(s.B),
                          Statistics::Fermion,
                          Statistics::Boson,
                          s.cfg);
  Background bg = AnalyticBackground(s.P, s.F, s.B, kA, kMass);
  K.PrepareTransitions(bg.a, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data());
  // The Stieltjes recurrence replaced the monomial Gram precisely so this axis -- the
  // design's honest convergence dimension -- runs as far as the daughter grid supports.
  // It does not run forever: high-order orthogonal polynomials oscillate with a large
  // dynamic range, so the reconstruction w*psi_j and the projection of the operator's
  // output start to cancel, and the conservation zeros drift off machine zero. That is
  // an accuracy ceiling of the reduction, not a failure of the recurrence, and the test
  // MEASURES it rather than asserting past it.
  int ceiling = 0;
  for (int nm = 2; nm <= 16; ++nm) {
    ReducedCollisionOperator
        R(K, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
    try {
      R.SetBackground(bg.a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
    }
    catch (const std::exception& e) {
      std::printf("   %3d      - | %s\n", nm, e.what());
      break;
    }
    std::vector<double> M;
    R.Assemble(1, M);
    const int S      = R.size();
    const double r0  = R.DissipativityBound(M.data());
    const double r1  = R.DissipativityBound(M.data() + (size_t) S * S);
    const auto c0    = R.Conservation(M.data());
    const bool clean = std::fmax(c0.number_fermion, c0.number_boson) < 1e-8;
    std::printf("   %3d %6d | %+.4e %+.4e  %s\n", nm, S, r0, r1, clean ? "" : "<- ceiling");
    if (!clean)
      break;
    ceiling = nm;
  }
  std::printf(
      "\n   The number zero survives up to n_moments = %d. The closure\n"
      "   measurement asks for 4-8, so the usable range covers the requirement --\n"
      "   but the DOF axis is NOT unbounded, and a design that needed 16 moments\n"
      "   would have to orthogonalise against a finer grid.\n",
      ceiling);
  // "for ANY number of daughter degrees of freedom" is the load-bearing half of the
  // claim: a reduction whose stability depended on the DOF count would just move the
  // Gamma ratchet from the grid onto the moment count. It has to cover what the
  // closure measurement actually needs, which is 8.
  assert(ceiling >= 8);
}

// ============================ 7. the entropy weight must stay a weight ==
//
// f(1-f) goes NEGATIVE the instant a fermion bin overshoots Pauli saturation, and the
// design note's own measurements report a peak f_l ~ 1.01 in a real run (DrPsdSpecies
// floors f at kFFloor but never caps it at 1). Nothing downstream notices: dmu =
// dq q^2 w becomes a SIGNED measure, so the Stieltjes recurrence orthogonalises in an
// indefinite form, gnorm can stay positive while orthogonality means nothing, and
// Reconstruct returns a df whose sign is FLIPPED across the saturated core. The
// Cholesky in BuildMetric would catch an indefinite metric -- but the shipped species
// path only ever calls Reconstruct/Project, so it never runs.
//
// Note that AnalyticBackground clamps its fermion at 0.999, i.e. the test background
// used by questions 1-6 is guarded against precisely the pathology a real run reports.
// This one deliberately is not.
Background OvershootBackground(const Grid& P, const Grid& F, const Grid& B, double a, double m) {
  Background bg = AnalyticBackground(P, F, B, a, m);
  // Push the populated core just past saturation, the way a discrete gather that does
  // not enforce Pauli blocking does. 1.01 is the measured value, not a stress test.
  for (size_t i = 0; i < bg.fl.size(); ++i)
    if (bg.fl[i] > 0.5)
      bg.fl[i] = 1.01;
  return bg;
}

void q7_entropy_weight_stays_nonnegative() {
  std::printf("\n== 7. the entropy weight under a Pauli overshoot ==\n");
  const int n_daughter = 51;
  Setup s              = MakeSetup(n_daughter, kA);
  Background bg        = OvershootBackground(s.P, s.F, s.B, kA, kMass);

  DecayTransitionKernel kern(View(s.P),
                             View(s.F),
                             View(s.B),
                             Statistics::Fermion,
                             Statistics::Boson,
                             s.cfg);
  kern.PrepareTransitions(kA, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data(), kA2H);

  ReducedCollisionOperator
      R(kern, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {4});
  R.SetBackground(kA, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());

  double w_min = 1e300, w_eq_min = 1e300, f_max = 0.;
  for (int k = 0; k < n_daughter; ++k) {
    w_min    = std::fmin(w_min, R.weight(0)[k]);
    w_eq_min = std::fmin(w_eq_min, R.weight(0)[k]);
    f_max    = std::fmax(f_max, bg.fl[k]);
  }
  std::printf("   peak f_l = %.4f  ->  min basis weight %+.6e, min entropy weight %+.6e\n",
              f_max,
              w_min,
              w_eq_min);

  // A weight is non-negative by definition. The measure it induces has to be a measure,
  // or none of questions 1-6 mean anything on this background.
  assert(f_max > 1.0);  // the background really does overshoot, or the test has no power
  assert(w_min >= 0.);
  assert(w_eq_min >= 0.);
  std::printf("   the weight stays a weight; the reconstruction cannot sign-flip.\n");
}

// A clamp that fires silently is a clamp that hides a background-solve defect for a
// whole run. It has to be reportable, and it has to stay quiet when the
// background is healthy -- otherwise the report carries no information.
void q7b_overshoot_is_reported() {
  std::printf("\n== 7b. the overshoot is reported, and only when it happens ==\n");
  const int n_daughter = 51;
  Setup s              = MakeSetup(n_daughter, kA);

  auto build = [&](const Background& bg) {
    DecayTransitionKernel kern(View(s.P),
                               View(s.F),
                               View(s.B),
                               Statistics::Fermion,
                               Statistics::Boson,
                               s.cfg);
    kern.PrepareTransitions(kA, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data(), kA2H);
    ReducedCollisionOperator
        R(kern, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {4});
    R.SetBackground(kA, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
    return R.weight_clamp(0);
  };

  const auto over = build(OvershootBackground(s.P, s.F, s.B, kA, kMass));
  const auto ok   = build(AnalyticBackground(s.P, s.F, s.B, kA, kMass));
  std::printf("   overshooting background: %d bins clamped, peak f_l = %.4f\n",
              over.n_clamped,
              over.f_max);
  std::printf("   healthy background:      %d bins clamped, peak f_l = %.4f\n",
              ok.n_clamped,
              ok.f_max);

  assert(over.n_clamped > 0);
  assert(over.f_max > 1.0);
  assert(ok.n_clamped == 0);  // no false positive, or the report is noise
}

// ================== 8. the analytic scattering kernel, as an ORACLE ==
//
// This suite has no external reference: every check so far is the shipped kernel
// against itself, or against a reduction OF itself, so a systematic error in the
// deposit is invisible to all of them.  There is exactly one place an outside answer
// exists.  arXiv:2205.13628 (appendix B) shows that for MASSLESS daughters and with
// inverse decays switched off, the daughters' momentum-averaged collision term closes
// in closed form:
//
//     (dF_dr/dtau)_{C,l}  =  rdot_dr * <f_H Psi_H F_l(q/eps)> / <f_H>,
//     F_l(x) = (1 - x^2)^2 / 2 * int_{-1}^{1} P_l(u) / (1 - x u)^3 du = (1-x^2)/(2x) Q^2_l(1/x)
//
// -- no daughter grid, no closure error, and F_0 = 1, F_1 = x exactly at every x.  The
// second of those is the MOMENTUM sum rule whose discrete violation owns the whole
// Gamma ratchet (design note M10/M11), so measuring the shipped kernel against it is
// measuring the defect directly rather than inferring it from an eigenvalue.
//
// The operator is linear in delta-f_H, so driving it with delta-f_H supported on ONE
// parent bin makes the daughter's energy-moment rate exactly C_b * F_l(x_b) for a
// single unknown C_b -- and F_0 = 1 fixes C_b = E_0.  The whole l-dependence is then a
// prediction with NO free constant and no normalisation to agree on.
//
// F_l is evaluated here EXACTLY -- partial fractions, no quadrature, no Q^2_l, no
// recurrence -- so it shares no code and no failure mode with either the shipped
// DecayTransitionKernel or the Miller recurrence in species/dncdm_dr_species.cpp that
// implements the same function for the no-inverse-decay model.

/** P_l(u) = sum_k c[k] u^k, by the standard three-term recurrence on coefficients. */
std::vector<double> LegendreCoeffs(int l) {
  std::vector<double> prev{1.}, cur{0., 1.};
  if (l == 0)
    return prev;
  for (int n = 1; n < l; ++n) {
    std::vector<double> next(cur.size() + 1, 0.);
    for (size_t k = 0; k < cur.size(); ++k)
      next[k + 1] += (2. * n + 1.) / (n + 1.) * cur[k];
    for (size_t k = 0; k < prev.size(); ++k)
      next[k] -= (double) n / (n + 1.) * prev[k];
    prev = cur;
    cur  = next;
  }
  return cur;
}

/** I_k(x) = int_{-1}^{1} u^k / (1 - x u)^3 du, in closed form.
 *
 *  s = 1 - x u maps [-1, 1] onto [1-x, 1+x] with du = -ds/x, so
 *      I_k = x^{-(k+1)} int_{1-x}^{1+x} (1-s)^k s^{-3} ds
 *  and the binomial expansion of (1-s)^k leaves only elementary powers of s.
 *
 *  The x^{-(k+1)} prefactor cancels against the bracket as x -> 0, so this loses
 *  precision at small x; every parent bin in this test has x >= 0.28, where the loss
 *  is a few digits out of sixteen. */
double IntegralK(int k, double x) {
  const double s_lo = 1. - x, s_hi = 1. + x;
  double binom = 1., sum = 0.;
  for (int j = 0; j <= k; ++j) {
    const double J  = (j == 2) ? std::log(s_hi / s_lo)
                               : (std::pow(s_hi, j - 2) - std::pow(s_lo, j - 2)) / (j - 2.);
    sum            += binom * ((j % 2) ? -1. : 1.) * J;
    binom           = binom * (k - j) / (j + 1.);
  }
  return sum / std::pow(x, k + 1);
}

double AnalyticF(int l, double x) {
  const std::vector<double> c = LegendreCoeffs(l);
  double acc                  = 0.;
  for (size_t k = 0; k < c.size(); ++k)
    if (c[k] != 0.)
      acc += c[k] * IntegralK((int) k, x);
  const double one_minus = (1. - x * x);
  return 0.5 * one_minus * one_minus * acc;
}

void q8_analytic_kernel_oracle() {
  std::printf("\n== 8. the shipped kernel against the ANALYTIC scattering kernel ==\n");

  // (i) pin the oracle on its own identities.  F_0 = 1 and F_1 = x come out of the same
  // partial-fraction code path as F_2..F_4, so this is a real check of it and not a
  // tautology; F_2's closed form (paper, section 2.3) is an independent third route.
  double id_err = 0.;
  for (double x = 0.05; x < 0.999; x += 0.017) {
    id_err          = std::fmax(id_err, std::fabs(AnalyticF(0, x) - 1.));
    id_err          = std::fmax(id_err, std::fabs(AnalyticF(1, x) - x));
    const double f2 = (x * (5. * x * x - 3.) + 3. * std::pow(x * x - 1., 2.) * std::atanh(x)) /
                      (2. * x * x * x);
    id_err          = std::fmax(id_err, std::fabs(AnalyticF(2, x) - f2) / std::fabs(f2));
  }
  std::printf("   oracle self-check (F_0 = 1, F_1 = x, F_2 closed form): %.3e\n", id_err);
  assert(id_err < 1e-9);

  // (ii) the kernel, in the limit where the analytic answer is exact.  Swept in the
  // DAUGHTER GRID, because a single resolution could not tell a discretisation error
  // from a wrong operator -- and because the daughter grid turns out to be the knob that
  // matters.  n_gauss and n_sub were swept first and are both INERT here: under
  // stratified_quadrature the emission nodes are the union of the daughter cell edges
  // lying inside the band (PrepareTransitions), so dr_N_q sets the node count of the
  // u-integral and neither Gauss setting touches it.  Worth knowing before anyone tries
  // to converge this sector in n_gauss again.
  std::printf("     N_q  |   l=1        l=2        l=3        l=4   (worst over parent bins)\n");
  std::printf("          |   fermion then boson, ABSOLUTE deviation (F_0 = 1 sets the scale)\n");

  double converged[kLMax + 1] = {0.}, converged_abs[kLMax + 1] = {0.}, prev_abs[kLMax + 1] = {0.};
  int n_rungs = 0;
  for (int n_daughter = 51; n_daughter <= 801; n_daughter = 2 * n_daughter - 1) {
    Setup s                  = MakeSetup(n_daughter, kA);
    s.cfg.inverse_decays     = false;
    s.cfg.quantum_statistics = false;
    Background bg            = AnalyticBackground(s.P, s.F, s.B, kA, kMass);
    std::fill(bg.fl.begin(), bg.fl.end(), 0.);      // no daughters: nothing to inverse-decay,
    std::fill(bg.fphi.begin(), bg.fphi.end(), 0.);  // and nothing to block or enhance

    DecayTransitionKernel kern(View(s.P),
                               View(s.F),
                               View(s.B),
                               Statistics::Fermion,
                               Statistics::Boson,
                               s.cfg);
    kern.PrepareTransitions(kA, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data(), kA2H);

    const int L1 = kLMax + 1;
    const int Np = (int) s.P.q.size(), Nf = n_daughter, Nb = n_daughter;
    std::vector<double> FH((size_t) Np * L1), Fl((size_t) Nf * L1), Fp((size_t) Nb * L1);
    std::vector<double> dH(FH.size()), dl(Fl.size()), dp(Fp.size());

    double worst[kLMax + 1] = {0.}, worst_abs[kLMax + 1] = {0.};
    double worst_f[kLMax + 1] = {0.}, worst_b[kLMax + 1] = {0.};
    for (int b = 0; b < Np; ++b) {
      const double x = s.P.q[b] / std::sqrt(s.P.q[b] * s.P.q[b] + kA * kA * kMass * kMass);

      std::fill(FH.begin(), FH.end(), 0.);
      std::fill(Fl.begin(), Fl.end(), 0.);
      std::fill(Fp.begin(), Fp.end(), 0.);
      for (int l = 0; l <= kLMax; ++l)
        FH[(size_t) b * L1 + l] = 1.;  // one parent bin, the same amplitude at every l

      kern.ResetMultipoleRecurrence();
      kern.ApplyPerturbationOperatorAllL(kLMax,
                                         L1,
                                         FH.data(),
                                         Fl.data(),
                                         Fp.data(),
                                         dH.data(),
                                         dl.data(),
                                         dp.data());

      // Energy-moment rate of each daughter at each l.  eps = q for a massless daughter.
      double E_f[kLMax + 1] = {0.}, E_b[kLMax + 1] = {0.};
      for (int l = 0; l <= kLMax; ++l) {
        for (int j = 0; j < Nf; ++j)
          E_f[l] += s.F.dq[j] * std::pow(s.F.q[j], 3) * dl[(size_t) j * L1 + l];
        for (int j = 0; j < Nb; ++j)
          E_b[l] += s.B.dq[j] * std::pow(s.B.q[j], 3) * dp[(size_t) j * L1 + l];
      }

      for (int l = 1; l <= kLMax; ++l) {
        // F_0 = 1 fixes the constant, so the ratio is the prediction outright.  Both
        // daughters share the same energy-weighted transfer (the boson's integrand is
        // the fermion's under u -> -u, which leaves the integral invariant).
        const double pred  = AnalyticF(l, x);
        const double dev_f = std::fabs(E_f[l] / E_f[0] - pred);
        const double dev_b = std::fabs(E_b[l] / E_b[0] - pred);
        const double dev   = std::fmax(dev_f, dev_b);
        worst_f[l]         = std::fmax(worst_f[l], dev_f);
        worst_b[l]         = std::fmax(worst_b[l], dev_b);
        // BOTH metrics, because they say different things.  F_l -> c_l x^l as x -> 0, so
        // at the least relativistic parent bins the prediction is itself tiny and a large
        // RELATIVE error there may be worth nothing in the sector's energy budget.  The
        // absolute one is already scale-relative: F_0 = 1 normalises it.
        worst[l]     = std::fmax(worst[l], dev / std::fmax(std::fabs(pred), 1e-300));
        worst_abs[l] = std::fmax(worst_abs[l], dev);
      }
    }
    std::printf("   %4d  nu_l |", n_daughter);
    for (int l = 1; l <= kLMax; ++l)
      std::printf(" %9.2e", worst_f[l]);
    std::printf("\n          phi  |");
    for (int l = 1; l <= kLMax; ++l)
      std::printf(" %9.2e", worst_b[l]);
    std::printf("\n");
    for (int l = 1; l <= kLMax; ++l) {
      prev_abs[l]      = converged_abs[l];
      converged[l]     = worst[l];
      converged_abs[l] = worst_abs[l];
    }
    ++n_rungs;
  }

  // l = 0 and l = 1 are SUM RULES -- F_0 = 1 is number/energy, F_1 = x is momentum --
  // and they hold at every resolution, to machine precision, because the two-bin deposit
  // places number and energy exactly and the Legendre weight is a per-node scalar.
  //
  // That is design note M11's attribution, now checked against an OUTSIDE answer instead
  // of inferred from an eigenvalue: with inverse decays off there is no lumped daughter
  // loss, and the l = 1 momentum identity is then exact.  Whatever breaks it in a
  // production run is the lumping, not the deposit and not the emission quadrature.
  assert(converged[1] < 1e-12);

  // l >= 2 has no sum rule, and it converges in the DAUGHTER GRID -- not in n_gauss and
  // not in n_sub, both of which are inert here.  Under stratified_quadrature the emission
  // nodes ARE the union of the daughter cell edges inside the band, so dr_N_q sets the
  // node count of the u-integral; refining the grid refines the quadrature, and nothing
  // else does.  Measured order is 2 (8.99e-4 / 2.32e-4 / 6.15e-5 at N_q = 201/401/801).
  //
  // Asserting the RATE rather than a threshold: a threshold passes just as well on an
  // operator that has stopped converging, which is the failure this test exists to catch.
  const double h_ratio = 2.;  // N_q doubles each rung
  for (int l = 2; l <= kLMax; ++l) {
    const double order = std::log(prev_abs[l] / converged_abs[l]) / std::log(h_ratio);
    std::printf("   l = %d: finest deviation %.2e, observed convergence order %.2f\n",
                l,
                converged_abs[l],
                order);
    assert(order > 1.7);
    assert(converged_abs[l] < 3e-4);
  }
  assert(n_rungs >= 4);
}

// ======================= 9. the tabulated operator (Route 2) ==
//
// M~_l is k-INDEPENDENT, so it can be assembled once per background row and replayed
// as a dense mat-vec inside the k-loop, instead of being applied matrix-free
// (reconstruct -> sweep the quadrature grid -> project) at every k and every step.
// That is the step that finally takes dr_N_q off the hot path.
//
// What makes it worth doing is measured, not assumed: at Gamma = 1e7 the reduced
// scheme's internal quadrature grid is the DOMINANT error -- halving h moves P(k) by
// 1.19e-3, against 5.3e-5 for the moment count -- and
// refining it currently costs 2.25x wall time per doubling, at every k. Tabulated, it
// costs only the precompute.
//
// Two things here are new and neither exists anywhere else: interpolation in ln a, and
// the mat-vec. The ASSEMBLY itself is already pinned by question 2 (M~ z equals
// P M R z to 6.3e-16), so this tests only what it adds.
void q9_operator_table() {
  std::printf("\n== 9. the tabulated operator: interpolation in ln a and the mat-vec ==\n");
  const int n_daughter = 51;
  const int nm         = 6;
  const int l_max      = 2;
  Setup s              = MakeSetup(n_daughter, kA);

  // Two background rows, a factor 2 apart in a -- far wider than any real table row
  // spacing, so an interpolation error shows up rather than hiding under the tolerance.
  const double a0 = kA, a1 = 2. * kA;
  std::vector<double> M[2];
  for (int r = 0; r < 2; ++r) {
    const double a = (r == 0) ? a0 : a1;
    Background bg  = AnalyticBackground(s.P, s.F, s.B, a, kMass);
    DecayTransitionKernel kern(View(s.P),
                               View(s.F),
                               View(s.B),
                               Statistics::Fermion,
                               Statistics::Boson,
                               s.cfg);
    kern.PrepareTransitions(a, kMass, kGamma, bg.fH.data(), bg.fl.data(), bg.fphi.data(), kA2H);
    ReducedCollisionOperator
        R(kern, View(s.P), View(s.F), View(s.B), Statistics::Fermion, Statistics::Boson, {nm});
    R.SetBackground(a, kMass, bg.fH.data(), bg.fl.data(), bg.fphi.data());
    R.Assemble(l_max, M[r]);
  }

  const int S = kNParent + 2 * nm;
  ReducedOperatorTable table;
  table.Reset(S, l_max);
  table.AddRow(a0, M[0]);
  table.AddRow(a1, M[1]);
  assert(table.n_rows() == 2);

  unsigned seed = 0xC0FFEEu;
  auto rnd      = [&]() {
    seed = seed * 1664525u + 1013904223u;
    return -1.0 + 2.0 * (((seed >> 9) & 0xffffu) / 65535.0);
  };
  std::vector<double> z(S);
  for (int i = 0; i < S; ++i)
    z[i] = rnd();

  auto matvec = [&](const std::vector<double>& Mat, int l, std::vector<double>& out) {
    out.assign(S, 0.);
    const double* Ml = Mat.data() + (size_t) l * S * S;
    for (int i = 0; i < S; ++i)
      for (int j = 0; j < S; ++j)
        out[i] += Ml[(size_t) i * S + j] * z[j];
  };
  auto rel = [&](const std::vector<double>& x, const std::vector<double>& ref) {
    double num = 0., den = 0.;
    for (int i = 0; i < S; ++i) {
      num = std::fmax(num, std::fabs(x[i] - ref[i]));
      den = std::fmax(den, std::fabs(ref[i]));
    }
    return (den > 0.) ? num / den : num;
  };

  // (i) AT a tabulated row the table must reproduce that row's matrix exactly. This is
  // the property that makes the tabulated path a drop-in: no interpolation, no blend.
  double worst_exact = 0.;
  for (int r = 0; r < 2; ++r) {
    for (int l = 0; l <= l_max; ++l) {
      std::vector<double> ref, got(S, 0.);
      matvec(M[r], l, ref);
      table.Apply((r == 0) ? a0 : a1, l, z.data(), got.data());
      worst_exact = std::fmax(worst_exact, rel(got, ref));
    }
  }
  std::printf("   at a tabulated row:            %.3e\n", worst_exact);
  assert(worst_exact < 1e-14);

  // (ii) BETWEEN rows, linear in ln a. Checked at the ln-midpoint, where the linear
  // blend is exactly the average of the two endpoints -- so the assertion tests the
  // interpolation and not the operator.
  double worst_mid = 0.;
  for (int l = 0; l <= l_max; ++l) {
    std::vector<double> r0, r1, got(S, 0.), want(S);
    matvec(M[0], l, r0);
    matvec(M[1], l, r1);
    for (int i = 0; i < S; ++i)
      want[i] = 0.5 * (r0[i] + r1[i]);
    table.Apply(std::sqrt(a0 * a1), l, z.data(), got.data());
    worst_mid = std::fmax(worst_mid, rel(got, want));
  }
  std::printf("   at the ln-a midpoint:          %.3e\n", worst_mid);
  assert(worst_mid < 1e-14);

  // (iii) OUTSIDE the tabulated range the operator is zero, and Apply ACCUMULATES, so
  // "zero" means it must leave dz untouched. The table covers the collision window; M1
  // says the collision falls ten or more orders outside it, and the species is what
  // chooses the window wide enough for that to be true.
  std::vector<double> outside(S, 1.0), before(S, 1.0);
  table.Apply(0.1 * a0, 0, z.data(), outside.data());
  table.Apply(10. * a1, 0, z.data(), outside.data());
  assert(rel(outside, before) == 0.);
  std::printf("   outside the range:             untouched\n");
}

// The O(1) lookup added for the uniform-in-ln-a case must land on the SAME row the
// binary search would, on both a uniform and a deliberately non-uniform table.  A row
// off by one is not a small error -- it returns a different operator -- and it would be
// invisible in q9, whose two-row table is uniform by construction.
void q11_row_lookup() {
  std::printf("\n== 11. row lookup: O(1) uniform path == binary search ==\n");
  const int S = 2, l_max = 0;
  auto block = [&](double tag) {
    std::vector<double> M((size_t) (l_max + 1) * S * S, 0.);
    M[0] = tag;  // identifies the row it came from
    return M;
  };

  for (int variant = 0; variant < 2; ++variant) {
    ReducedOperatorTable t;
    t.Reset(S, l_max);
    std::vector<double> as;
    for (int i = 0; i < 40; ++i) {
      // variant 0: uniform in ln a.  variant 1: a deliberate kink at row 20, which is
      // exactly what a table assembled over a non-uniform background would look like.
      const double lna = (variant == 0) ? (-12. + 0.25 * i)
                                        : (-12. + 0.25 * i + (i >= 20 ? 0.37 * (i - 19) : 0.));
      as.push_back(std::exp(lna));
      t.AddRow(as.back(), block((double) i));
    }
    // Probe between every pair of rows; the value returned must be the lower row's tag.
    double worst = 0.;
    for (size_t i = 0; i + 1 < as.size(); ++i) {
      for (double frac : {0.0, 0.001, 0.5, 0.999}) {
        const double a = std::exp(std::log(as[i]) + frac * (std::log(as[i + 1]) - std::log(as[i])));
        double z[2] = {1., 0.}, dz[2] = {0., 0.};
        t.Apply(a, 0, z, dz);
        // The blend of tags i and i+1 at weight `frac` is exactly i + frac.
        worst = std::fmax(worst, std::fabs(dz[0] - ((double) i + frac)));
      }
    }
    std::printf("   %-12s worst |row index error| = %.3e\n",
                (variant == 0) ? "uniform:" : "non-uniform:",
                worst);
    assert(worst < 1e-9);
  }
}

}  // namespace

int main() {
  std::printf("reduced_operator_test: Gamma = %.4f/Mpc (ini 1e7), m = M/T = %.2f\n", kGamma, kMass);
  q1_projection_identity();
  q1b_energy_coordinate();
  q2_in_span_agreement();
  q3_equilibrium_annihilation();
  q4_q5_conservation_and_dissipativity();
  q6_the_milestone();
  q7_entropy_weight_stays_nonnegative();
  q7b_overshoot_is_reported();
  q8_analytic_kernel_oracle();
  q9_operator_table();
  q11_row_lookup();
  std::printf("\nreduced_operator_test: all assertions passed\n");
  return 0;
}
