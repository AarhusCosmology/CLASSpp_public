// Unit tests for DecayTransitionKernel (tools/decay_transition_kernel): the pure
// conservative transition network for nu_H <-> nu_l + phi (massless daughters).
// The load-bearing property is machine-precision number+energy conservation for
// ARBITRARY positive f arrays; the equilibrium residual is the only soft
// (convergence-ratio) criterion. Physics: arXiv:2011.01502; design spec section 4.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>

#include "decay_transition_kernel.h"
#include "quadrature.h"  // same Gauss-Legendre rule the kernel ctor uses (P_l reference)

namespace {

struct Grid {
  std::vector<double> q, dq;
};

// Log-spaced trapezoid grid with explicit q_min (thesis grid convention): half
// cells at the endpoints so Sum dq[i] f[i] approximates the integral.
Grid MakeGrid(double qmin, double qmax, int N) {
  Grid g;
  g.q.resize(N);
  g.dq.resize(N);
  const double h = std::log(qmax / qmin) / (N - 1);
  for (int i = 0; i < N; ++i)
    g.q[i] = qmin * std::exp(i * h);
  g.dq[0] = 0.5 * (g.q[1] - g.q[0]);
  for (int i = 1; i < N - 1; ++i)
    g.dq[i] = 0.5 * (g.q[i + 1] - g.q[i - 1]);
  g.dq[N - 1] = 0.5 * (g.q[N - 1] - g.q[N - 2]);
  return g;
}

GridView View(const Grid& g) {
  return GridView{g.q.data(), g.dq.data(), static_cast<int>(g.q.size())};
}

double Scale(double a, double b) {
  return std::fmax(std::fmax(std::fabs(a), std::fabs(b)), 1e-300);
}
double Scale3(double a, double b, double c) {
  return std::fmax(std::fmax(std::fabs(a), std::fabs(b)), std::fmax(std::fabs(c), 1e-300));
}

double MaxAbs(const std::vector<double>& a,
              const std::vector<double>& b,
              const std::vector<double>& c) {
  double m = 0.;
  for (double v : a)
    m = std::fmax(m, std::fabs(v));
  for (double v : b)
    m = std::fmax(m, std::fabs(v));
  for (double v : c)
    m = std::fmax(m, std::fabs(v));
  return m;
}

// MOMENT CONVENTION (Fable review -- load-bearing): ComputeMoments returns BARE
// per-dof grid moments N_i = Sum dq q^2 df_i, E_i = Sum dq q^2 eps_i df_i. The x2
// boson dof ratio (g_H*g_l/g_phi) is folded into df_phi (paper eq. 4.14). With
// g_H = g_l = 2, g_phi = 1 one decay H -> l + phi gives bare moments N_l = -N_H,
// N_phi = -2*N_H, and App. C reduces to EXACTLY:
//   N_H + N_l == 0            (fermion-leg number)
//   2*N_H + N_phi == 0        (boson-leg number)
//   2*E_H + 2*E_l + E_phi == 0  (g-weighted energy; per transition
//                                2 eps1 - 2 q2* - 2 q3* == 0 since q2*+q3* = eps1)
// These hold for ARBITRARY f -- a failure is a gather/scatter transpose or a
// 1/(q^2 dq) conversion bug, never a tolerance to loosen.
// `balanced` is a free parameter here on purpose: the discrete number and energy
// identities are properties of the DEPOSIT weights and hold for ANY value of the band
// factor, so a gather change must leave them at round-off. That is the check that
// Config::balanced_gather really did only move a coefficient.
void test_conservation(bool balanced = false) {
  Grid P      = MakeGrid(1e-3, 1e2, 64);
  Grid F      = MakeGrid(1e-4, 2e2, 64);  // wider than the parent bands: no off-grid clamp
  Grid B      = MakeGrid(1e-4, 2e2, 64);
  const int N = 64;

  std::mt19937 rng(12345u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::vector<double> fH(N), fl(N), fphi(N);
  for (auto& x : fH)
    x = U(rng);
  for (auto& x : fl)
    x = U(rng);
  for (auto& x : fphi)
    x = U(rng);

  const double a = 0.5, m = 2.0, Gamma = 0.3;  // A = a^2 m^2 = 1 -> bands in [1e-4,2e2]
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);

  std::vector<double> dfH(N), dfl(N), dfphi(N);
  K.ComputeBackgroundDerivs(a,
                            m,
                            Gamma,
                            fH.data(),
                            fl.data(),
                            fphi.data(),
                            dfH.data(),
                            dfl.data(),
                            dfphi.data());
  assert(K.clamped_energy_residual() == 0.);  // grids cover the bands

  auto M             = K.ComputeMoments(a, m, dfH.data(), dfl.data(), dfphi.data());
  const double num_f = std::fabs(M.N_H + M.N_l) / Scale(M.N_H, M.N_l);
  const double num_b = std::fabs(2. * M.N_H + M.N_phi) / Scale(2. * M.N_H, M.N_phi);
  // ENERGY, post-positivity-fix: the identity is exact up to the energy the lumped
  // daughter loss deliberately misplaces, and THAT is booked exactly. So the test
  // is not weakened to a tolerance -- it still pins a machine-precision statement,
  // just the right one: no energy goes missing unaccounted. The residual itself is
  // reported because this f is uniform random per bin, i.e. the roughest input the
  // scheme can be handed; see test_daughter_positivity for its convergence.
  const double escale = Scale3(2. * M.E_H, 2. * M.E_l, M.E_phi);
  const double ene    = std::fabs(2. * M.E_H + 2. * M.E_l + M.E_phi - K.split_energy_residual()) /
                        escale;
  std::printf(
      "conservation (random f)%s: |N_H+N_l|=%.2e |2N_H+N_phi|=%.2e "
      "|2E_H+2E_l+E_phi - booked|=%.2e (booked %.2e of scale)\n",
      balanced ? " (balanced_gather)" : "",
      num_f,
      num_b,
      ene,
      std::fabs(K.split_energy_residual()) / escale);
  assert(num_f < 1e-12);
  assert(num_b < 1e-12);
  assert(ene < 1e-12);
  assert(std::fabs(M.N_H) > 0.);  // guard against a trivially-zero kernel "passing"
}

// Equilibrium residual (Fable review -- replaces the unachievable hard
// detailed-balance tolerance). For FD parent+fermion, BE boson with common T and
// mu_H = mu_l + mu_phi the continuum Lambda == 0 on eps1 = q2*+q3*; on the grid
// the gathered f_l, f_phi are linear interpolants so Lambda_s = O(dq^2), tangent
// to the conservation manifold (N/E still exact). Normalize by the decay-only
// rate (same f/grids, resolution-independent) and test the convergence ratio.
double parent_only_resid = 0.;
double EquilibResidual(int Nd, bool balanced = false) {
  Grid P         = MakeGrid(1e-3, 1e2, 64);  // parent grid FIXED across the refinement
  Grid F         = MakeGrid(1e-4, 2e2, Nd);
  Grid B         = MakeGrid(1e-4, 2e2, Nd);
  const double a = 0.5, m = 2.0, Gamma = 0.3, A = a * a * m * m;
  const double T = 1.0, mu_l = -0.5, mu_phi = -2.0,
               mu_H = mu_l + mu_phi;  // mu_phi < q_min: BE regular
  auto fH_of   = [&](double q) { return 1. / (std::exp((std::sqrt(q * q + A) - mu_H) / T) + 1.); };
  auto fl_of   = [&](double q) { return 1. / (std::exp((q - mu_l) / T) + 1.); };
  auto fphi_of = [&](double q) { return 1. / (std::exp((q - mu_phi) / T) - 1.); };

  std::vector<double> fH(64), fl(Nd), fphi(Nd);
  for (int i = 0; i < 64; ++i)
    fH[i] = fH_of(P.q[i]);
  for (int i = 0; i < Nd; ++i) {
    fl[i]   = fl_of(F.q[i]);
    fphi[i] = fphi_of(B.q[i]);
  }

  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  DecayTransitionKernel Kfull(View(P),
                              View(F),
                              View(B),
                              Statistics::Fermion,
                              Statistics::Boson,
                              cfg);
  std::vector<double> dfH(64), dfl(Nd), dfphi(Nd);
  Kfull.ComputeBackgroundDerivs(a,
                                m,
                                Gamma,
                                fH.data(),
                                fl.data(),
                                fphi.data(),
                                dfH.data(),
                                dfl.data(),
                                dfphi.data());

  auto M             = Kfull.ComputeMoments(a, m, dfH.data(), dfl.data(), dfphi.data());
  const double resid = MaxAbs(dfH, dfl, dfphi);
  // Split the residual: the PARENT leg is df_H = dN/(q1^2 dq1) with dN = C Lambda,
  // so it isolates the band factor alone. The DAUGHTER legs additionally carry the
  // lumped loss, whose per-bin deviation term survives even at Lambda == 0.
  std::vector<double> empty;
  parent_only_resid = MaxAbs(dfH, empty, empty);

  DecayTransitionKernel::Config cfg0;
  cfg0.inverse_decays     = false;
  cfg0.quantum_statistics = false;
  DecayTransitionKernel Kdec(View(P),
                             View(F),
                             View(B),
                             Statistics::Fermion,
                             Statistics::Boson,
                             cfg0);
  std::vector<double> gH(64), gl(Nd), gphi(Nd);
  Kdec.ComputeBackgroundDerivs(a,
                               m,
                               Gamma,
                               fH.data(),
                               fl.data(),
                               fphi.data(),
                               gH.data(),
                               gl.data(),
                               gphi.data());
  const double scale = MaxAbs(gH, gl, gphi);

  // N/E conservation is EXACT regardless of resolution -- assert here too, but
  // against the DECAY-ONLY rate rather than against the net moments themselves.
  // At equilibrium the net moments are the cancellation this function exists to
  // drive to zero, so self-scaling would ask round-off to cancel against ~0 and fail
  // on success. The decay-only rate is the same order as the individual node
  // contributions that must cancel, which is the honest scale for "did the sum of
  // large terms come out conserved".
  auto Md         = Kdec.ComputeMoments(a, m, gH.data(), gl.data(), gphi.data());
  const double ns = Scale(Md.N_H, Md.N_l);
  const double es = Scale3(2. * Md.E_H, 2. * Md.E_l, Md.E_phi);
  assert(std::fabs(M.N_H + M.N_l) < 1e-12 * ns);
  assert(std::fabs(2. * M.N_H + M.N_phi) < 1e-12 * ns);
  assert(std::fabs(2. * M.E_H + 2. * M.E_l + M.E_phi - Kfull.split_energy_residual()) < 1e-12 * es);
  return resid / scale;
}

void test_equilibrium() {
  // DETAILED BALANCE. In equilibrium the continuum band factor vanishes identically,
  // so what is left is pure discretisation error and it must CONVERGE. A quadrature
  // that smoothed the PSD by distorting the band factor would show up here as a
  // residual that stops shrinking.
  const double r64  = EquilibResidual(64);
  const double p64  = parent_only_resid;
  const double r128 = EquilibResidual(128);
  const double p128 = parent_only_resid;
  std::printf(
      "equilibrium residual: r64=%.3e r128=%.3e ratio=%.2f (need > 2.5); "
      "parent leg alone %.3e / %.3e\n",
      r64,
      r128,
      r64 / r128,
      p64,
      p128);
  assert(r128 < r64 / 2.5);

  // Config::balanced_gather: the band factor is built from the log-odds, which are
  // EXACTLY linear in q for this state, so the parent leg -- which is Lambda alone --
  // must come back at round-off and STAY there under refinement. That flatness is the
  // claim: a scheme that merely converged faster would still fall with Nd.
  const double b64   = EquilibResidual(64, true);
  const double bp64  = parent_only_resid;
  const double b128  = EquilibResidual(128, true);
  const double bp128 = parent_only_resid;
  std::printf(
      "   balanced_gather: r64=%.3e r128=%.3e; parent leg alone %.3e / %.3e "
      "(linear gather: %.3e / %.3e)\n",
      b64,
      b128,
      bp64,
      bp128,
      p64,
      p128);
  // The parent leg is Lambda with nothing else on it; the daughter legs additionally
  // carry the lumped loss deviation, which survives at Lambda == 0 and is NOT what
  // this flag addresses, so the well-balanced claim is asserted where it is made.
  // These are residuals already normalised by the decay-only rate, so 1e-15 is a
  // round-off claim, and asserting it at BOTH resolutions is the flatness claim.
  assert(bp64 < 1e-15);
  assert(bp128 < 1e-15);
}

// Mode reduction: inverse+qs off => Lambda = -f_H(q1) => the parent decay term
// reduces to the exact continuum diagonal df_H(q1_i) = -K/eps1 * f_H(q1_i)
// (matches dncdm_species.cpp diagonal decay; K = a^2 m Gamma).
void test_mode_reduction() {
  Grid P = MakeGrid(1e-3, 1e2, 64);
  Grid F = MakeGrid(1e-4, 2e2, 48);
  Grid B = MakeGrid(1e-4, 2e2, 48);

  std::mt19937 rng(777u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::vector<double> fH(64), fl(48), fphi(48);
  for (auto& x : fH)
    x = U(rng);
  for (auto& x : fl)
    x = U(rng);
  for (auto& x : fphi)
    x = U(rng);

  const double a = 0.5, m = 2.0, Gamma = 0.3, A = a * a * m * m, K = a * a * m * Gamma;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = false;
  cfg.quantum_statistics = false;
  // Sampled on purpose: the parent's decay term is -K/eps1*f_H analytically, so it
  // must NOT depend on how the daughter side lays out its nodes. Pinning the older
  // scheme here keeps that a statement about the parent leg rather than about the
  // default quadrature.
  DecayTransitionKernel Kd(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
  std::vector<double> dfH(64), dfl(48), dfphi(48);
  Kd.ComputeBackgroundDerivs(a,
                             m,
                             Gamma,
                             fH.data(),
                             fl.data(),
                             fphi.data(),
                             dfH.data(),
                             dfl.data(),
                             dfphi.data());

  double maxrel = 0.;
  for (int i = 0; i < 64; ++i) {
    const double eps1   = std::sqrt(P.q[i] * P.q[i] + A);
    const double expect = -K / eps1 * fH[i];
    maxrel              = std::fmax(maxrel, std::fabs(dfH[i] - expect) / Scale(expect, dfH[i]));
  }
  std::printf("mode reduction: max rel dev of parent decay term = %.2e\n", maxrel);
  assert(maxrel < 1e-12);
}

// RHS continuity across band-edge node crossings. Sweep a finely so the band
// edges (eps1 +- q1)/2 and the interior Gauss nodes cross daughter grid points;
// the RHS must be continuous (kinks only C0 as a node crosses a grid point --
// design 4.2). The 2020-2024 branch's #1 bug was a DISCONTINUOUS jump here
// a step would blow the finite-difference ratio up.
//
// CRITERION: the da-REFINEMENT ratio, not a magnitude threshold. m = |df(a+da) -
// df(a)|/(da*scale) tends to |d(df)/da|/scale for a continuous RHS -- FLAT in da --
// but reads J/(da*scale) for a jump J, which DOUBLES every time da halves. So
// m(da/2)/m(da) is ~1 when continuous and ~2 when not, and it needs no calibration.
//
// It replaces an absolute `m < 50`, which was a proxy for the same thing and broke
// the moment stratified quadrature became the default: stratified measures 125 with
// NO discontinuity whatsoever (verified -- 16x refinement of da moves it 125.191 ->
// 125.316, i.e. flat to 0.1%). It is legitimately steeper because it actually
// delivers the band-edge emission into the low-q boson bins that the sampled nodes
// mostly stepped over: at the maximum, df_phi[5] = 1.4e6 against a `scale` of
// 2.2e5, and the kink at a ~ 0.2112 is parent bins 60/61/62 walking their band's
// lower edge (eps1-q1)/2 across boson cell 5. A magnitude bound cannot tell that
// from a defect; the refinement ratio can.
//
// PER-BIN, deliberately: a jump in bin j doubles bin j's OWN ratio no matter how
// large a different bin's honest derivative is. Measured sensitivity, by injecting
// a step of known size and checking this assertion fires -- sampled catches a jump
// of 1e-3*scale, stratified 1e-1*scale (its low-q boson bins have a large genuine
// d/da for the injected jump to hide behind). Globally-maxed instead of per-bin,
// those floors are 10x worse. Running BOTH quadratures is therefore load-bearing
// for sensitivity and not merely for coverage: the sampled sweep is the sharp one.
void test_rhs_continuity() {
  Grid P = MakeGrid(1e-3, 1e2, 64);
  Grid F = MakeGrid(1e-4, 2e2, 40);  // coarse daughter grid: nodes cross bins during the sweep
  Grid B = MakeGrid(1e-4, 2e2, 40);
  const int Np = 64, Nd = 40;

  std::mt19937 rng(2024u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::vector<double> fH(Np), fl(Nd), fphi(Nd);
  for (auto& x : fH)
    x = U(rng);
  for (auto& x : fl)
    x = U(rng);
  for (auto& x : fphi)
    x = U(rng);

  const double m = 2.0, Gamma = 0.3;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  // Both quadratures, and stratified is the one that needs it most: its breakpoint
  // SET changes discretely as the band sweeps (a daughter cell edge enters or leaves
  // the band), which is a fresh opportunity for the jump this test exists to catch.
  // Each such event contributes a piece of zero width, so the RHS stays C0 -- that is
  // the claim under test, not an assumption.
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);

  std::vector<double> dfH(Np), dfl(Nd), dfphi(Nd), pH(Np), pl(Nd), pphi(Nd);
  const double a0 = 0.2, a1 = 0.8;
  K.ComputeBackgroundDerivs(0.5,
                            m,
                            Gamma,
                            fH.data(),
                            fl.data(),
                            fphi.data(),
                            dfH.data(),
                            dfl.data(),
                            dfphi.data());
  const double scale = MaxAbs(dfH, dfl, dfphi);

  // Per-bin max of |df(a+da) - df(a)|/(da*scale) over the sweep; layout is
  // [parent | fermion | boson] so one index names a bin across all three arrays.
  auto sweep = [&](double da) {
    K.ComputeBackgroundDerivs(a0,
                              m,
                              Gamma,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              pH.data(),
                              pl.data(),
                              pphi.data());
    std::vector<double> per_bin(Np + 2 * Nd, 0.);
    for (double a = a0 + da; a < a1; a += da) {
      K.ComputeBackgroundDerivs(a,
                                m,
                                Gamma,
                                fH.data(),
                                fl.data(),
                                fphi.data(),
                                dfH.data(),
                                dfl.data(),
                                dfphi.data());
      for (int i = 0; i < Np; ++i)
        per_bin[i] = std::fmax(per_bin[i], std::fabs(dfH[i] - pH[i]) / (da * scale));
      for (int j = 0; j < Nd; ++j) {
        per_bin[Np + j]      = std::fmax(per_bin[Np + j], std::fabs(dfl[j] - pl[j]) / (da * scale));
        per_bin[Np + Nd + j] = std::fmax(per_bin[Np + Nd + j],
                                         std::fabs(dfphi[j] - pphi[j]) / (da * scale));
      }
      pH   = dfH;
      pl   = dfl;
      pphi = dfphi;
    }
    return per_bin;
  };

  const double da   = 2e-4;
  const auto coarse = sweep(da);
  const auto fine   = sweep(0.5 * da);
  // Bins whose motion is at roundoff carry no information and a noisy ratio.
  const double kFloor = 1e-6;
  double worst_R = 0., worst_m = 0.;
  int at = -1, n_live = 0;
  for (size_t i = 0; i < coarse.size(); ++i) {
    worst_m = std::fmax(worst_m, coarse[i]);
    if (coarse[i] <= kFloor)
      continue;
    ++n_live;
    const double R = fine[i] / coarse[i];
    if (R > worst_R) {
      worst_R = R;
      at      = static_cast<int>(i);
    }
  }
  std::printf(
      "rhs continuity: worst per-bin da-refinement ratio = %.3f at bin %d "
      "(need < 1.5; 2.0 = jump), largest |d(df)/da|/scale = %.2f over %d live bins\n",
      worst_R,
      at,
      worst_m,
      n_live);
  assert(worst_R < 1.5);  // ~2 would mean a step, i.e. a discontinuous RHS
  assert(n_live > Np);    // else the sweep moved almost nothing and the test is vacuous
}

// Perturbation operator identities (design 4.3). ApplyPerturbationOperator is a pure
// UNNORMALIZED F-space operator: in F_i = f_i Psi_i = delta-f_i, out (dF_i/dtau)^(1)_{C,l}.
// No 1/f_i normalization anywhere (the caller applies it only where its own state
// variable is normalized -- the parent). The output obeys the SAME g-weighted moment
// identities as the background network, read DIRECTLY off the returned arrays:
//   l=0: num_i = Sum dq q^2 dF_i,0 conserves number (num_H+num_l==0,
//        2num_H+num_phi==0) and energy (2e_H+2e_l+e_phi==0, e uses eps_i);
//   l=1: mom_i = Sum dq q^3 dF_i,1 conserves momentum (2mom_H+2mom_l+mom_phi==0),
//        because q1 = q2* cos_a* + q3* cos_b* holds exactly at every node.
// These hold for ARBITRARY F (equilibrium or not) — a failure is a gather/scatter
// transpose or angular-projection bug, never a tolerance to loosen.
//
// Run under both quadratures, and for the same reason as test_conservation: the
// operator must remain the exact Jacobian of the background network under EITHER
// node layout, and the l=1 momentum identity in particular leans on
// q1 = q2* cos_a* + q3* cos_b* holding at every node -- which the stratified
// breakpoints generate by a different route than the Gauss-Legendre ones.
void test_perturbation_operator(bool balanced = false) {
  Grid P      = MakeGrid(1e-3, 1e2, 64);
  Grid F      = MakeGrid(1e-4, 2e2, 64);  // wider than the parent bands: no off-grid clamp
  Grid B      = MakeGrid(1e-4, 2e2, 64);
  const int N = 64;

  std::mt19937 rng(98765u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::uniform_real_distribution<double> Upm(-1.0, 1.0);
  std::vector<double> fH(N), fl(N), fphi(N), FH(N), Fl(N), Fphi(N);
  for (auto& x : fH)
    x = U(rng);
  for (auto& x : fl)
    x = U(rng);
  for (auto& x : fphi)
    x = U(rng);
  // F = delta-f, an INDEPENDENT perturbation array (not tied to f): the identities
  // must hold for arbitrary F, so seeding it from f would weaken the test.
  for (auto& x : FH)
    x = Upm(rng);
  for (auto& x : Fl)
    x = Upm(rng);
  for (auto& x : Fphi)
    x = Upm(rng);

  const double a = 0.5, m = 2.0, Gamma = 0.3, A = a * a * m * m;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  // The chain rule changes what the operator gathers, never how it scatters, so every
  // identity below must hold unchanged under it. That is the point of running this
  // whole test a second time rather than adding a separate one.
  cfg.balanced_gather = balanced;
  cfg.balanced_pert   = balanced;
  cfg.lumped_loss     = !balanced;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());

  // P_l built by recurrence (never a table): spot-check P_2 = 1.5 cos^2 - 0.5 at a node.
  {
    double q2, q3, ca, cb, cg;
    K.NodeKinematics(10, 2, q2, q3, ca, cb, cg);
    assert(std::fabs(DecayTransitionKernel::LegendreP(2, ca) - (1.5 * ca * ca - 0.5)) < 1e-13);
    assert(std::fabs(DecayTransitionKernel::LegendreP(0, cb) - 1.0) < 1e-14);
    assert(std::fabs(DecayTransitionKernel::LegendreP(1, cg) - cg) < 1e-14);
  }

  // l=0: number & energy conservation of the perturbed operator.
  {
    std::vector<double> dH(N), dl(N), dphi(N);
    K.ApplyPerturbationOperator(0,
                                FH.data(),
                                Fl.data(),
                                Fphi.data(),
                                dH.data(),
                                dl.data(),
                                dphi.data());
    double numH = 0, numl = 0, numphi = 0, eH = 0, el = 0, ephi = 0;
    for (int i = 0; i < N; ++i) {
      const double q = P.q[i], w = P.dq[i] * q * q, eps = std::sqrt(q * q + A);
      numH += w * dH[i];
      eH   += w * eps * dH[i];
    }
    for (int j = 0; j < N; ++j) {
      const double q = F.q[j], w = F.dq[j] * q * q;
      numl += w * dl[j];
      el   += w * q * dl[j];
    }
    for (int k = 0; k < N; ++k) {
      const double q = B.q[k], w = B.dq[k] * q * q;
      numphi += w * dphi[k];
      ephi   += w * q * dphi[k];
    }
    // NUMBER stays exactly conserved under the lumped operator -- lumping moves where
    // a deposit lands, never how much -- so those two are still asserted at zero. The
    // q-weighted ENERGY identity carries the O(dq) misplacement, booked per call.
    const double nf = std::fabs(numH + numl) / Scale(numH, numl);
    const double nb = std::fabs(2. * numH + numphi) / Scale(2. * numH, numphi);
    const double bk = K.split_energy_residual_pert();
    const double en = std::fabs(2. * eH + 2. * el + ephi - bk) / Scale3(2. * eH, 2. * el, ephi);
    std::printf(
        "pert l=0: |numH+numl|=%.2e |2numH+numphi|=%.2e "
        "|2eH+2el+ephi - booked|=%.2e (booked %.2e)\n",
        nf,
        nb,
        en,
        std::fabs(bk) / Scale3(2. * eH, 2. * el, ephi));
    assert(nf < 1e-11);
    assert(nb < 1e-11);
    assert(en < 1e-11);
    assert(bk != 0.);              // else the energy assertion is vacuous
    assert(std::fabs(numH) > 0.);  // guard against a trivially-zero operator "passing"
  }

  // l=1: momentum conservation of the perturbed operator.
  {
    std::vector<double> dH(N), dl(N), dphi(N);
    K.ApplyPerturbationOperator(1,
                                FH.data(),
                                Fl.data(),
                                Fphi.data(),
                                dH.data(),
                                dl.data(),
                                dphi.data());
    double momH = 0, moml = 0, momphi = 0;
    for (int i = 0; i < N; ++i) {
      const double q  = P.q[i];
      momH           += P.dq[i] * q * q * q * dH[i];
    }
    for (int j = 0; j < N; ++j) {
      const double q  = F.q[j];
      moml           += F.dq[j] * q * q * q * dl[j];
    }
    for (int k = 0; k < N; ++k) {
      const double q  = B.q[k];
      momphi         += B.dq[k] * q * q * q * dphi[k];
    }
    // Same first-moment structure as the l=0 energy identity, same booked residual.
    const double bk = K.split_momentum_residual_pert();
    const double mm = std::fabs(2. * momH + 2. * moml + momphi - bk) /
                      Scale3(2. * momH, 2. * moml, momphi);
    std::printf("pert l=1: |2momH+2moml+momphi - booked|=%.2e (booked %.2e)\n",
                mm,
                std::fabs(bk) / Scale3(2. * momH, 2. * moml, momphi));
    assert(mm < 1e-10);
    assert(bk != 0.);
    assert(std::fabs(momH) > 0.);  // guard against a trivially-zero operator
  }

  // Decay-only reduction: inverse+qs off => the parent operator is pure decay
  // dF_H = -(K/eps1) F_H (design 4.3 cancellation source: after the caller's /f_H it
  // is exactly minus the parent's own f-bar-dot/f-bar). Daughters still receive the
  // parent's perturbation (dN_l ~ -s_H P_l(cos_a)), only the DAUGHTER->anything legs
  // vanish (c_l = c_phi = 0).
  {
    DecayTransitionKernel::Config cfg0;
    cfg0.inverse_decays     = false;
    cfg0.quantum_statistics = false;  // parent leg is analytic: mode-independent
    DecayTransitionKernel Kd(View(P),
                             View(F),
                             View(B),
                             Statistics::Fermion,
                             Statistics::Boson,
                             cfg0);
    Kd.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
    std::vector<double> dH(N), dl(N), dphi(N);
    // The operator advances its per-node Legendre recurrence one l per call, so l=2 is
    // only reachable through l=0 and l=1 (whose outputs are discarded here); calling it
    // cold at l=2 is an API violation and now throws. The parent's decay-only result is
    // l-independent anyway (-K/eps1 F_H carries no angular factor), which is exactly
    // what makes it the cancellation source -- so l=2 is still the meaningful probe:
    // it confirms no P_l leaked into the parent leg.
    for (int lp = 0; lp <= 2; ++lp)
      Kd.ApplyPerturbationOperator(lp,
                                   FH.data(),
                                   Fl.data(),
                                   Fphi.data(),
                                   dH.data(),
                                   dl.data(),
                                   dphi.data());
    double maxrel = 0.;
    for (int i = 0; i < N; ++i) {
      const double eps1   = std::sqrt(P.q[i] * P.q[i] + A);
      const double expect = -(a * a * m * Gamma) / eps1 * FH[i];
      maxrel              = std::fmax(maxrel, std::fabs(dH[i] - expect) / Scale(expect, dH[i]));
    }
    std::printf("pert decay-only: max rel dev of parent operator = %.2e\n", maxrel);
    assert(maxrel < 1e-12);
  }

  // Sequential-call contract. The operator carries a rolling per-node (P_{l-2},P_{l-1})
  // Legendre state, so (a) an out-of-order l must be rejected rather than silently
  // evaluate the wrong multipole, and (b) PrepareTransitions must restart the
  // recurrence, making the l-sequence a pure function of the prepared state: the SAME
  // sweep repeated after a re-prepare has to come back BIT-identical, which is what
  // pins the cache reset and the prev/curr shift.
  {
    const int LMAX = 5;
    std::vector<double> sweep1, sweep2;
    for (int pass = 0; pass < 2; ++pass) {
      K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
      bool threw = false;
      std::vector<double> dH(N), dl(N), dphi(N);
      try {  // cold kernel expects l=0; anything else is API misuse
        K.ApplyPerturbationOperator(1,
                                    FH.data(),
                                    Fl.data(),
                                    Fphi.data(),
                                    dH.data(),
                                    dl.data(),
                                    dphi.data());
      }
      catch (const std::invalid_argument&) {
        threw = true;
      }
      assert(threw);
      std::vector<double>& out = pass == 0 ? sweep1 : sweep2;
      for (int l = 0; l <= LMAX; ++l) {
        K.ApplyPerturbationOperator(l,
                                    FH.data(),
                                    Fl.data(),
                                    Fphi.data(),
                                    dH.data(),
                                    dl.data(),
                                    dphi.data());
        out.insert(out.end(), dH.begin(), dH.end());
        out.insert(out.end(), dl.begin(), dl.end());
        out.insert(out.end(), dphi.begin(), dphi.end());
      }
      // Re-running an already-consumed l is equally rejected.
      bool threw_repeat = false;
      try {
        K.ApplyPerturbationOperator(LMAX,
                                    FH.data(),
                                    Fl.data(),
                                    Fphi.data(),
                                    dH.data(),
                                    dl.data(),
                                    dphi.data());
      }
      catch (const std::invalid_argument&) {
        threw_repeat = true;
      }
      assert(threw_repeat);
    }
    assert(sweep1.size() == sweep2.size());
    bool identical = true;
    for (std::size_t i = 0; i < sweep1.size(); ++i)
      identical &= (sweep1[i] == sweep2[i]);
    std::printf("pert l-sweep: out-of-order rejected, l=0..%d repeat bit-identical = %s\n",
                LMAX,
                identical ? "yes" : "NO");
    assert(identical);
  }

  // Angular factor at HIGH l. The l=0/l=1 identities above cannot see an error in P_l
  // (P_0 = 1, and the momentum identity uses cos_a*/cos_b* only through l=1), so the
  // rolling recurrence needs its own anchor at l >= 2. Decay-only (c_l = c_phi = 0) with
  // ONLY the parent perturbed leaves the fermion deposit per node at
  // dN_l = B_n s_H,n P_l(cos_a*_n), s_H,n = -F_H[i]; the two-bin split is
  // number-preserving, so the fermion NUMBER moment collapses to a closed form,
  //   Sum_j dq q^2 dF_l[j] = Sum_n B_n F_H[i] P_l(cos_a*_n),
  // and likewise the boson moment with cos_b* and the x2 dof factor. B_n = K dq1
  // (q1/eps1) wn is rebuilt here from the kernel's own node census and weights, and
  // P_l from the public from-scratch LegendreP -- so a stale cache, a shift-order slip
  // or an off-by-one in l shows up as an O(1) deviation. What this pins is the ANGULAR
  // factor, not the quadrature: the node list is taken as given and only P_l(cos) and
  // the l-shift are under test.
  {
    DecayTransitionKernel::Config cfg0;
    cfg0.inverse_decays     = false;
    cfg0.quantum_statistics = false;
    DecayTransitionKernel Kd(View(P),
                             View(F),
                             View(B),
                             Statistics::Fermion,
                             Statistics::Boson,
                             cfg0);
    Kd.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());

    std::vector<double> zero(N, 0.), dH(N), dl(N), dphi(N);
    double worst = 0.;
    for (int l = 0; l <= 12; ++l) {
      Kd.ApplyPerturbationOperator(l,
                                   FH.data(),
                                   zero.data(),
                                   zero.data(),
                                   dH.data(),
                                   dl.data(),
                                   dphi.data());
      double got_l = 0., got_phi = 0.;
      for (int j = 0; j < N; ++j)
        got_l += F.dq[j] * F.q[j] * F.q[j] * dl[j];
      for (int k = 0; k < N; ++k)
        got_phi += B.dq[k] * B.q[k] * B.q[k] * dphi[k];
      double ref_l = 0., ref_phi = 0.;
      for (int i = 0; i < N; ++i) {
        const double q1 = P.q[i], eps1 = std::sqrt(q1 * q1 + A);
        for (int sN = 0; sN < Kd.node_count(i); ++sN) {
          double q2, q3, ca, cb, cg;
          Kd.NodeKinematics(i, sN, q2, q3, ca, cb, cg);
          const double Bn  = (a * a * m * Gamma) * P.dq[i] * (q1 / eps1) * Kd.node_weight(i, sN);
          ref_l           += Bn * FH[i] * DecayTransitionKernel::LegendreP(l, ca);
          ref_phi         += 2.0 * Bn * FH[i] * DecayTransitionKernel::LegendreP(l, cb);
        }
      }
      worst = std::fmax(worst, std::fabs(got_l - ref_l) / Scale(got_l, ref_l));
      worst = std::fmax(worst, std::fabs(got_phi - ref_phi) / Scale(got_phi, ref_phi));
    }
    std::printf("pert P_l(l=0..12): max rel dev of daughter number moments = %.2e\n", worst);
    assert(worst < 1e-12);
  }
}

// CollisionDiagonal must equal the true d(df_i)/d(f_i) of ComputeBackgroundDerivs,
// because an exponential / IMEX step integrates exactly that and carries the rest
// explicitly. If the analytic diagonal drifts from the operator, the scheme is
// stable but converges to the WRONG attractor -- a silent physics error, not a
// crash. So it is pinned against central differences on the operator itself.
//
// The diagonal's DOMINANCE is deliberately reported and not asserted. It is a
// property of the trajectory, not of the operator: measured on the real solution at
// a = 0.3 the off-diagonal sits 5-9 orders below the diagonal, but that is because
// the emission band has swept out to q ~ a*m/2 while the daughters still occupy
// q ~ 1-6, so there is nothing at the band momenta for the inverse decay to couple
// to. Populate the daughters uniformly instead -- physically impossible, but a
// legal input -- and the ratio comes back to ~1 at EVERY a*m (swept 1 to 1e4).
// So an assertion here would either encode an unphysical state or silently depend
// on the daughter profile a future change might alter. What the scheme needs to be
// CORRECT is the diagonal matching the operator, which is asserted; whether it is
// FAST is a question about trajectories, which the timings answer.
void test_collision_diagonal(bool balanced = false) {
  Grid P       = MakeGrid(1e-3, 1e2, 24);
  Grid F       = MakeGrid(1e-4, 2e2, 48);
  Grid B       = MakeGrid(1e-4, 2e2, 48);
  const int Np = 24, Nd = 48, N = Np + 2 * Nd;

  std::mt19937 rng(31337u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::vector<double> y(N);
  for (int i = 0; i < Np; ++i)
    y[i] = U(rng);
  for (int j = 0; j < 2 * Nd; ++j)
    y[Np + j] = U(rng);

  const double a = 0.5, m = 2.0, Gamma = 0.3;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);

  auto rhs = [&](const std::vector<double>& yy, std::vector<double>& out) {
    K.ComputeBackgroundDerivs(a,
                              m,
                              Gamma,
                              yy.data(),
                              yy.data() + Np,
                              yy.data() + Np + Nd,
                              out.data(),
                              out.data() + Np,
                              out.data() + Np + Nd);
  };

  // Analytic diagonal (needs the transitions prepared for THIS state).
  std::vector<double> d0(N);
  rhs(y, d0);
  std::vector<double> dH(Np), dl(Nd), dphi(Nd);
  K.CollisionDiagonal(dH.data(), dl.data(), dphi.data());
  std::vector<double> analytic(N);
  for (int i = 0; i < Np; ++i)
    analytic[i] = dH[i];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + j] = dl[j];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + Nd + j] = dphi[j];

  // Central-difference Jacobian: diagonal for the comparison, full for the premise.
  std::vector<double> Jrow(N), dp(N), dm(N), offsum(N, 0.), numeric(N);
  double worst = 0., dscale = 0.;
  int at = -1;
  for (int c = 0; c < N; ++c) {
    const double h         = std::fmax(1e-7 * std::fabs(y[c]), 1e-13);
    std::vector<double> yp = y, ym = y;
    yp[c] += h;
    ym[c] -= h;
    rhs(yp, dp);
    rhs(ym, dm);
    for (int r = 0; r < N; ++r) {
      const double J = (dp[r] - dm[r]) / (2 * h);
      if (r == c)
        numeric[r] = J;
      else
        offsum[r] += std::fabs(J);
    }
  }
  for (int i = 0; i < N; ++i)
    dscale = std::fmax(dscale, std::fabs(numeric[i]));
  for (int i = 0; i < N; ++i) {
    const double e = std::fabs(analytic[i] - numeric[i]) / Scale(numeric[i], dscale);
    if (e > worst) {
      worst = e;
      at    = i;
    }
  }
  std::printf("collision diagonal%s: max rel dev from finite differences = %.2e at index %d\n",
              balanced ? " (balanced_gather)" : "",
              worst,
              at);
  assert(worst < 1e-6);  // central differences on a smooth RHS; not a machine-eps claim
  assert(dscale > 0.);   // else the comparison is vacuous

  // Reported, not asserted (see the header comment): the late-time state the solver
  // actually spends its steps in -- parent extinct, and daughters confined to the
  // momenta they physically occupy rather than smeared over the whole grid, which is
  // what makes the band's deposit region empty.
  for (int i = 0; i < Np; ++i)
    y[i] *= 1e-60;
  for (int j = 0; j < Nd; ++j) {
    const double lq = std::log(F.q[j] / 2.8);
    y[Np + j]       = 1. / (std::exp(F.q[j]) + 1.) + 0.3 * std::exp(-lq * lq);
    const double lb = std::log(B.q[j] / 2.8);
    y[Np + Nd + j]  = 0.5 * std::exp(-lb * lb);
  }
  rhs(y, d0);
  // Grab the analytic diagonal for THIS state before the sweep below re-prepares the
  // transitions. This is the interesting state for balanced_gather: the chain factor
  // is G(1-+G)/(f(1-+f)), so it is O(1) on the smooth random state above and
  // ENORMOUS here, where the Gaussian daughter profile leaves bins many decades
  // below their neighbours. If the linearisation is wrong anywhere it is here.
  K.CollisionDiagonal(dH.data(), dl.data(), dphi.data());
  for (int i = 0; i < Np; ++i)
    analytic[i] = dH[i];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + j] = dl[j];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + Nd + j] = dphi[j];
  std::fill(offsum.begin(), offsum.end(), 0.);
  for (int c = 0; c < N; ++c) {
    const double h         = std::fmax(1e-7 * std::fabs(y[c]), 1e-13);
    std::vector<double> yp = y, ym = y;
    yp[c] += h;
    ym[c] -= h;
    rhs(yp, dp);
    rhs(ym, dm);
    for (int r = 0; r < N; ++r) {
      const double J = (dp[r] - dm[r]) / (2 * h);
      if (r == c)
        numeric[r] = J;
      else
        offsum[r] += std::fabs(J);
    }
  }
  double Dmax = 0., Offmax = 0.;
  for (int i = 0; i < N; ++i) {
    Dmax   = std::fmax(Dmax, std::fabs(numeric[i]));
    Offmax = std::fmax(Offmax, offsum[i]);
  }
  double worst2 = 0.;
  int at2       = -1;
  for (int i = 0; i < N; ++i) {
    const double e = std::fabs(analytic[i] - numeric[i]) / Scale(numeric[i], Dmax);
    if (e > worst2) {
      worst2 = e;
      at2    = i;
    }
  }
  std::printf(
      "   extinct parent, physical daughter profile: max|J_ii| = %.3e, "
      "max row off-diagonal = %.3e (ratio %.2e -- reported, not asserted); "
      "diagonal dev = %.2e at %d\n",
      Dmax,
      Offmax,
      Offmax / Dmax,
      worst2,
      at2);
  assert(Dmax > 0.);  // the report must not be vacuous
  // worst2 is REPORTED, not asserted, and the linear path is the reason: with the
  // parent at 1e-60 the central difference's step floors at 1e-13, which is a huge
  // RELATIVE perturbation, so the numeric column there is the difference quotient of
  // a different state. It reads ~3e-1 on the shipped scheme too. What it is good for
  // is comparing the two gathers against each other on the same unsound estimator.

  // DYNAMIC RANGE, on a state the finite difference can still resolve. The chain
  // factor G(1-+G)/(f(1-+f)) is 1 by construction on a flat stencil, so a smooth
  // random state barely exercises it; here the daughters are log-uniform over eight
  // decades, which makes the factor range over many orders while every bin stays far
  // enough above the 1e-13 step floor for the quotient to mean something.
  std::mt19937 rng2(90210u);
  std::uniform_real_distribution<double> LU(-8., -0.05);
  for (int i = 0; i < Np; ++i)
    y[i] = U(rng2);
  for (int j = 0; j < 2 * Nd; ++j)
    y[Np + j] = std::pow(10., LU(rng2));
  rhs(y, d0);
  K.CollisionDiagonal(dH.data(), dl.data(), dphi.data());
  for (int i = 0; i < Np; ++i)
    analytic[i] = dH[i];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + j] = dl[j];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + Nd + j] = dphi[j];
  // A RELATIVE step, and a coarse one. The block above can floor at 1e-13 because its
  // state is O(1); here a bin may sit at 1e-8, where an absolute floor is a 1e-5
  // relative kick. And the step has to be coarse: df_c is dominated by transitions out
  // of far larger bins, so dp[c] - dm[c] is a cancellation and the estimator's
  // round-off branch takes over below h/y ~ 1e-4. Measured on the linear path,
  // h/y = 1e-4/1e-5/1e-6/1e-7 gives 5.3e-6/8.3e-5/4.2e-4/1.8e-3 -- growing as h
  // SHRINKS, i.e. the5.3e-6 is the estimator, not the diagonal.
  for (int c = 0; c < N; ++c) {
    const double h         = 1e-4 * std::fabs(y[c]);
    std::vector<double> yp = y, ym = y;
    yp[c] += h;
    ym[c] -= h;
    rhs(yp, dp);
    rhs(ym, dm);
    numeric[c] = (dp[c] - dm[c]) / (2 * h);
  }
  double dscale3 = 0., worst3 = 0.;
  int at3 = -1;
  for (int i = 0; i < N; ++i)
    dscale3 = std::fmax(dscale3, std::fabs(numeric[i]));
  for (int i = 0; i < N; ++i) {
    const double e = std::fabs(analytic[i] - numeric[i]) / Scale(numeric[i], dscale3);
    if (e > worst3) {
      worst3 = e;
      at3    = i;
    }
  }
  std::printf("   log-uniform daughters over 8 decades: max rel dev = %.2e at index %d\n",
              worst3,
              at3);
  // 1e-4 leaves ~20x over the estimator's own floor and still catches a dropped chain
  // factor, which is an O(1) relative error, not a small one.
  assert(worst3 < 1e-4);
  assert(dscale3 > 0.);
}

// Per-l diagonal of the PERTURBATION operator, against the background's
// CollisionDiagonal.
//
// PREMISE UNDER TEST. An exponential integrator on the perturbation path
// integrates d(dF_i)/d(F_i) exactly and everything else explicitly. That is the
// right object only if ApplyPerturbationOperator's diagonal really is the
// background network's diagonal -- which it should be, the operator being the
// exact Jacobian of that network (440ec131) -- and it is only CHEAP if the
// diagonal is l-INDEPENDENT, so one CollisionDiagonal call can serve the whole
// hierarchy instead of one per multipole.
//
// The operator is LINEAR, so column j is exactly what it returns for the unit
// vector e_j. No finite differences: the only tolerance here is roundoff between
// two independent code paths computing the same analytic quantity, which is why
// this asserts far tighter than test_collision_diagonal's 1e-6.
//
// The strict-ascending-l contract (the rolling Legendre recurrence) means column
// j at level l costs one PrepareTransitions plus a replay of ll = 0..l. Same
// shape as DumpCollisionSpectrum, which is where this construction comes from.
void test_perturbation_diagonal(bool balanced = false) {
  Grid P       = MakeGrid(1e-3, 1e2, 16);
  Grid F       = MakeGrid(1e-4, 2e2, 24);  // wider than the parent bands: no off-grid clamp
  Grid B       = MakeGrid(1e-4, 2e2, 24);
  const int Np = 16, Nd = 24, N = Np + 2 * Nd;
  const int LMAX = 4;

  std::mt19937 rng(4242u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::vector<double> fH(Np), fl(Nd), fphi(Nd);
  for (auto& x : fH)
    x = U(rng);
  for (auto& x : fl)
    x = U(rng);
  for (auto& x : fphi)
    x = U(rng);

  const double a = 0.5, m = 2.0, Gamma = 0.3;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  cfg.balanced_pert      = balanced;
  cfg.lumped_loss        = !balanced;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);

  // Analytic diagonal, taken ONCE from the background decomposition. The same
  // vector is compared at every l, so an l-dependent diagonal fails here.
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
  std::vector<double> dH(Np), dl(Nd), dphi(Nd);
  K.CollisionDiagonal(dH.data(), dl.data(), dphi.data());
  std::vector<double> analytic(N);
  for (int i = 0; i < Np; ++i)
    analytic[i] = dH[i];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + j] = dl[j];
  for (int j = 0; j < Nd; ++j)
    analytic[Np + Nd + j] = dphi[j];

  double scale = 0.;
  for (int i = 0; i < N; ++i)
    scale = std::fmax(scale, std::fabs(analytic[i]));
  assert(scale > 0.);  // else every comparison below is vacuous

  std::vector<double> FH(Np), Fl(Nd), Fphi(Nd), oH(Np), ol(Nd), ophi(Nd);
  double worst = 0., worst_parent = 0.;
  int at = -1, at_l = -1;

  for (int l = 0; l <= LMAX; ++l) {
    for (int c = 0; c < N; ++c) {
      K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
      std::fill(FH.begin(), FH.end(), 0.);
      std::fill(Fl.begin(), Fl.end(), 0.);
      std::fill(Fphi.begin(), Fphi.end(), 0.);
      if (c < Np)
        FH[c] = 1.;
      else if (c < Np + Nd)
        Fl[c - Np] = 1.;
      else
        Fphi[c - Np - Nd] = 1.;
      for (int ll = 0; ll <= l; ++ll)
        K.ApplyPerturbationOperator(ll,
                                    FH.data(),
                                    Fl.data(),
                                    Fphi.data(),
                                    oH.data(),
                                    ol.data(),
                                    ophi.data());
      const double got = (c < Np) ? oH[c] : (c < Np + Nd) ? ol[c - Np] : ophi[c - Np - Nd];
      const double e   = std::fabs(got - analytic[c]) / Scale(analytic[c], scale);
      if (e > worst) {
        worst = e;
        at    = c;
        at_l  = l;
      }
      if (c < Np)
        worst_parent = std::fmax(worst_parent, e);
    }
  }

  std::printf(
      "perturbation diagonal: max rel dev from CollisionDiagonal = %.2e at index %d "
      "l=%d; PARENT leg alone = %.2e (l=0..%d)\n",
      worst,
      at,
      at_l,
      worst_parent,
      LMAX);

  // The parent leg carries the stiff eigenvalue (diag_H = -K/eps1 -> -a*Gamma) and
  // is the whole point of the exercise, so it is asserted separately and first:
  // if the daughters turn out to be l-dependent, the parent-only diagonal is still
  // a correct and useful L (the evolver forms N = f - Ly, so the fixed point is
  // exact for any L).
  assert(worst_parent < 1e-9);
  assert(worst < 1e-9);
}

// The daughter grids must COVER the emission band, or TwoBinSplit clamps the deposit
// into the edge bin: number survives (the split still sums to dN) but the energy is
// filed at q_edge instead of q*, so the sector silently loses (1 - q_max/q*) of every
// decay. That is a real defect this branch shipped -- with q_max = 300 against a band
// reaching a*M/2 = 2973, the daughters received 40% of the parent's energy loss at
// z=3 and 20% at z=1. MaxDaughterMomentum is the a-priori bound that makes it
// impossible to size a grid wrongly; this test pins the bound to the band the kernel
// actually builds, so the two cannot drift apart.
void test_band_bound() {
  // Bound vs. the real node momenta, over parent grids and a*m values spanning the
  // relativistic (a*m << q1) and non-relativistic (a*m >> q1) regimes.
  for (double q1_max : {1e1, 1e2}) {
    Grid P = MakeGrid(1e-3, q1_max, 32);
    Grid F = MakeGrid(1e-6, 1e6, 32);  // deliberately huge: no clamping anywhere
    Grid B = MakeGrid(1e-6, 1e6, 32);
    DecayTransitionKernel::Config cfg;
    DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
    const int NP = static_cast<int>(P.q.size());
    std::vector<double> fH(NP, 0.5), fl(32, 0.), fphi(32, 0.);
    for (double am : {1e-3, 1e-1, 1.0, 1e1, 1e3, 1e5}) {
      const double a = 1.0, m = am;
      const double bound = DecayTransitionKernel::MaxDaughterMomentum(P.q.back(), a * m);
      double worst       = 0.;
      // The node list is laid out per call, so it has to exist before it can be walked.
      K.PrepareTransitions(a, m, 1.0, fH.data(), fl.data(), fphi.data());
      for (int i = 0; i < NP; ++i) {
        for (int sN = 0; sN < K.node_count(i); ++sN) {
          double q2, q3, ca, cb, cg;
          K.NodeKinematics(i, sN, q2, q3, ca, cb, cg);
          assert(q2 <= bound && q3 <= bound);  // the bound is a bound
          worst = std::fmax(worst, std::fmax(q2, q3));
        }
      }
      // ... and a USEFUL one: never more than the top parent cell above the band.
      assert(worst > 0.5 * bound);
      std::printf("band bound: q1_max=%.1e a*m=%.1e bound=%.6e worst node=%.6e (%.3f)\n",
                  q1_max,
                  am,
                  bound,
                  worst,
                  worst / bound);
    }
  }

  // Regression: a grid sized to the bound clamps nothing; one sized just under it does.
  // This is the defect in its smallest reproducible form.
  {
    Grid P         = MakeGrid(1e-3, 1e1, 32);
    const int N    = 32;
    const double a = 1.0, m = 6e3, Gamma = 0.3;  // non-relativistic parent, band at ~a*m/2
    const double bound = DecayTransitionKernel::MaxDaughterMomentum(P.q.back(), a * m);
    std::vector<double> fH(N, 0.3), fl(N, 0.), fphi(N, 0.);
    std::vector<double> dfH(N), dfl(N), dfphi(N);
    DecayTransitionKernel::Config cfg;

    Grid Fok = MakeGrid(1e-6, 1.05 * bound, N);
    DecayTransitionKernel Kok(View(P),
                              View(Fok),
                              View(Fok),
                              Statistics::Fermion,
                              Statistics::Boson,
                              cfg);
    Kok.ComputeBackgroundDerivs(a,
                                m,
                                Gamma,
                                fH.data(),
                                fl.data(),
                                fphi.data(),
                                dfH.data(),
                                dfl.data(),
                                dfphi.data());
    assert(Kok.clamped_energy_residual() == 0.);

    Grid Fbad = MakeGrid(1e-6, 0.1 * bound, N);
    DecayTransitionKernel Kbad(View(P),
                               View(Fbad),
                               View(Fbad),
                               Statistics::Fermion,
                               Statistics::Boson,
                               cfg);
    Kbad.ComputeBackgroundDerivs(a,
                                 m,
                                 Gamma,
                                 fH.data(),
                                 fl.data(),
                                 fphi.data(),
                                 dfH.data(),
                                 dfl.data(),
                                 dfphi.data());
    assert(Kbad.clamped_energy_residual() != 0.);
    std::printf("band bound: undersized grid leaks energy residual %.6e (as it must)\n",
                Kbad.clamped_energy_residual());
  }
}

// POSITIVITY OF THE DAUGHTER LEGS. An occupation that is zero can only go up: no
// process removes a particle from a bin that holds none. That is a property the
// RIGHT-HAND SIDE must have, and if it does not, no integrator can rescue it --
// step-size control, state clamps and changes of variable were all tried against
// this and all three stalled, because each holds the state against an RHS that
// points out of the physical region.
//
// The state below is the defect's smallest form: an injection FRONT, one populated
// boson bin surrounded by empty ones, with f_l > f_H so that inverse decays are net
// ABSORBING. The exact two-bin deposit makes the loss proportional to the
// INTERPOLATED occupation lambda f_k + (1-lambda) f_{k+1} and then charges a share
// of it back to bin k -- so an empty bin is billed for its neighbour's particles.
// Measured in production at Gamma=1e8: boson bin 6 saw gain +1.52, off-diagonal
// -1.59, i.e. RHS = -6.7e-2 at f = 0, and the state ran to -1.13.
//
// Conservation is asserted on the SAME state, because the fix must not buy
// positivity with number: the lumped loss removes exactly the same TOTAL number
// (sum_e w_e f_{k+e} is the gathered value by definition), only from the bins that
// actually hold the particles.
// `lumped`/`balanced` span the 2x2 because positivity is exactly what distinguishes
// the four schemes, and three of them CLAIM it. (linear, unlumped) is the known-broken
// corner and is run as a control: it must actually go negative, or the front states
// below are not exercising the defect and the other three pass for nothing.
void test_daughter_positivity(bool balanced = false, bool lumped = true) {
  double worst_neg = 0.;
  // Only the corner without a positivity argument is allowed to violate. For the other
  // three, `assert` fires on the spot; here we also record the worst excursion so the
  // control can be shown to be a real one.
  const bool must_hold = lumped || balanced;
  auto check           = [&](double df) {
    if (df < worst_neg)
      worst_neg = df;
    assert(!must_hold || df >= 0.);
  };
  Grid P       = MakeGrid(1e-3, 1e2, 24);
  Grid F       = MakeGrid(1e-4, 2e2, 48);
  Grid B       = MakeGrid(1e-4, 2e2, 48);
  const int Np = 24, Nd = 48;

  const double a = 0.5, m = 2.0, Gamma = 0.3;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  cfg.lumped_loss        = lumped;
  // Both quadratures: positivity is a property of the LUMPED loss, which is what
  // makes the operator Metzler, and lumping is independent of where the nodes sit.
  // It is checked under stratified as well because that is now the shipped path --
  // and because stratified spreads a band over many more bins, so far more bins are
  // near-empty-but-touched, which is exactly the state that used to go negative.
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
  std::vector<double> dfH(Np), dfl(Nd), dfphi(Nd);

  // (a) boson front: one populated bin, everything else empty. f_l > f_H makes
  //     dLambda/df_phi = f_l - f_H > 0, i.e. a net drain -- the production case.
  for (int spike = 12; spike < Nd - 12; spike += 7) {
    std::vector<double> fH(Np, 0.05), fl(Nd, 0.5), fphi(Nd, 0.);
    fphi[spike] = 2.0;
    K.ComputeBackgroundDerivs(a,
                              m,
                              Gamma,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dfH.data(),
                              dfl.data(),
                              dfphi.data());
    for (int k = 0; k < Nd; ++k) {
      if (fphi[k] != 0.)
        continue;
      check(dfphi[k]);  // an empty bin may only fill
    }
    auto M = K.ComputeMoments(a, m, dfH.data(), dfl.data(), dfphi.data());
    assert(std::fabs(M.N_H + M.N_l) < 1e-12 * Scale(M.N_H, M.N_l));
    assert(std::fabs(2. * M.N_H + M.N_phi) < 1e-12 * Scale(2. * M.N_H, M.N_phi));
  }

  // (b) fermion front. dLambda/df_l = f_phi + f_H >= 0 always, so the fermion leg
  //     is a drain in EVERY configuration -- it is only spared in production
  //     because a thermally-seeded nu_l has no front to fall off.
  for (int spike = 12; spike < Nd - 12; spike += 7) {
    std::vector<double> fH(Np, 0.05), fl(Nd, 0.), fphi(Nd, 0.5);
    fl[spike] = 0.9;
    K.ComputeBackgroundDerivs(a,
                              m,
                              Gamma,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dfH.data(),
                              dfl.data(),
                              dfphi.data());
    for (int k = 0; k < Nd; ++k) {
      if (fl[k] != 0.)
        continue;
      check(dfl[k]);
    }
    auto M = K.ComputeMoments(a, m, dfH.data(), dfl.data(), dfphi.data());
    assert(std::fabs(M.N_H + M.N_l) < 1e-12 * Scale(M.N_H, M.N_l));
    assert(std::fabs(2. * M.N_H + M.N_phi) < 1e-12 * Scale(2. * M.N_H, M.N_phi));
  }

  // (c) The energy the lumped split misplaces is BOOKED, never silent: the residual
  //     accessor must reproduce the energy identity's defect to machine precision.
  //     (Sibling of clamped_energy_residual(); both are zero for smooth f.)
  std::vector<double> fH(Np, 0.05), fl(Nd), fphi(Nd, 0.);
  for (int k = 0; k < Nd; ++k)
    fl[k] = 1. / (std::exp(F.q[k]) + 1.);
  fphi[26] = 2.0;
  K.ComputeBackgroundDerivs(a,
                            m,
                            Gamma,
                            fH.data(),
                            fl.data(),
                            fphi.data(),
                            dfH.data(),
                            dfl.data(),
                            dfphi.data());
  auto M              = K.ComputeMoments(a, m, dfH.data(), dfl.data(), dfphi.data());
  const double ident  = 2. * M.E_H + 2. * M.E_l + M.E_phi;
  const double booked = K.split_energy_residual();
  // With the exact split nothing is misplaced, so the residual is identically zero and
  // the energy identity closes on its own -- assert THAT rather than skipping, since
  // "no residual" is the substantive claim for that scheme. With lumping on, a nonzero
  // booking is required or the identity below is vacuous (a constant f lumps exactly).
  assert(lumped ? booked != 0. : booked == 0.);
  assert(std::fabs(ident - booked) < 1e-12 * Scale3(2. * M.E_H, 2. * M.E_l, M.E_phi));
  std::printf("positivity: energy identity %.6e == booked split residual %.6e\n", ident, booked);

  // (d) The residual is a DISCRETISATION error, not a bias: for smooth f it shrinks
  //     with the daughter grid. First order, because the misplacement per node is
  //     lambda(1-lambda)(q_k - q_k+1)(f_k - f_k+1) -- one power of dq. That is the
  //     price of positivity, and it is the whole price; number stays exact.
  double prev = 0.;
  for (int Nref : {48, 96, 192}) {
    Grid Fr = MakeGrid(1e-4, 2e2, Nref);
    DecayTransitionKernel Kr(View(P),
                             View(Fr),
                             View(Fr),
                             Statistics::Fermion,
                             Statistics::Boson,
                             cfg);
    std::vector<double> flr(Nref), fphir(Nref), dflr(Nref), dfphir(Nref);
    for (int k = 0; k < Nref; ++k) {
      flr[k]   = 1. / (std::exp(Fr.q[k]) + 1.);
      fphir[k] = 1. / (std::exp(Fr.q[k] + 2.) - 1.);
    }
    Kr.ComputeBackgroundDerivs(a,
                               m,
                               Gamma,
                               fH.data(),
                               flr.data(),
                               fphir.data(),
                               dfH.data(),
                               dflr.data(),
                               dfphir.data());
    auto Mr            = Kr.ComputeMoments(a, m, dfH.data(), dflr.data(), dfphir.data());
    const double scale = Scale3(2. * Mr.E_H, 2. * Mr.E_l, Mr.E_phi);
    const double rel   = std::fabs(Kr.split_energy_residual()) / scale;
    // Booked exactly, at every resolution.
    assert(std::fabs(2. * Mr.E_H + 2. * Mr.E_l + Mr.E_phi - Kr.split_energy_residual()) <
           1e-12 * scale);
    std::printf("positivity: N=%3d smooth-f split residual %.3e (of energy scale)%s\n",
                Nref,
                rel,
                prev > 0. ? "" : "");
    if (prev > 0.)
      assert(rel < 0.75 * prev);  // shrinking with resolution
    prev = rel;
  }

  std::printf("positivity [%s gather, %s loss]: worst empty-bin RHS = %.3e (%s)\n",
              balanced ? "balanced" : "linear",
              lumped ? "lumped" : "exact",
              worst_neg,
              must_hold ? "asserted >= 0" : "CONTROL: expected < 0");
  // The control must actually break, or the front states are not reaching the defect
  // and the three positive results above are vacuous.
  assert(must_hold || worst_neg < 0.);
}

// The deposit stencil's three invariants, for EVERY order. These are the
// properties the whole conservative-network design rests on, so a new deposit rule
// is only admissible if it keeps all three exactly -- they are not tolerances.
void test_deposit_moments() {
  Grid D           = MakeGrid(1e-4, 20.0, 200);  // a production-sized daughter grid
  double worst_num = 0., worst_ene = 0., most_negative = 0.;
  // Sweep q* densely across the interior, including exactly-on-node values.
  for (int i = 3; i < 196; ++i) {
    for (double frac : {0.0, 0.013, 0.25, 0.5, 0.749, 0.99}) {
      const double q_star = D.q[i] + frac * (D.q[i + 1] - D.q[i]);
      DecayTransitionKernel::DepositStencil st;
      DecayTransitionKernel::BuildDeposit(View(D), q_star, st);
      double sw = 0., sq = 0.;
      for (int e = 0; e < DecayTransitionKernel::DepositStencil::kWidth; ++e) {
        sw            += st.w[e];
        sq            += st.w[e] * D.q[st.j0 + e];
        most_negative  = std::fmin(most_negative, st.w[e]);
      }
      worst_num = std::fmax(worst_num, std::fabs(sw - 1.0));
      worst_ene = std::fmax(worst_ene, std::fabs(sq - q_star) / q_star);
    }
  }
  std::printf("deposit: |sum w - 1|=%.2e  |sum w q - q*|/q*=%.2e  min w=%.2e\n",
              worst_num,
              worst_ene,
              most_negative);
  assert(worst_num < 1e-14);    // number, exactly
  assert(worst_ene < 1e-13);    // energy, exactly
  assert(most_negative >= 0.);  // positivity: an empty bin may only fill
}

// THE BAND PARTITION. The emission band of parent bin i is cut at every fermion cell
// edge inside it AND at every image eps1 - q_k of a boson cell edge, and one midpoint
// node is placed per resulting sub-interval. Two things have to hold, and neither is
// visible to the conservation identities:
//
//   * the node weights over a band must sum to EXACTLY q1, its width. A dropped
//     sub-interval under-counts the decay rate while leaving every PER-TRANSITION
//     identity satisfied, so nothing else in this file would see it;
//   * the node count must equal the number of sub-intervals the union partition
//     implies, so a merge that skipped one stream's edges is caught structurally
//     rather than as a small numerical drift.
//
// Then the quadrature itself: one midpoint node per cell is only enough if the
// integrand is smooth across a cell, so refining BOTH daughter grids (which refines
// the partition) must leave the total decay rate essentially unchanged.
void test_band_partition() {
  const double a = 0.5, m = 2.0, Gamma = 0.3;

  auto total_rate = [&](int ND) {
    Grid P       = MakeGrid(1e-3, 1e2, 16);
    Grid F       = MakeGrid(1e-4, 2e2, ND);
    Grid B       = MakeGrid(1e-4, 2e2, ND);
    const int NP = 16;
    std::vector<double> fH(NP), fl(ND), fphi(ND);
    for (int i = 0; i < NP; ++i)
      fH[i] = 1. / (std::exp(P.q[i]) + 1.);
    for (int i = 0; i < ND; ++i) {
      fl[i]   = 0.10 / (std::exp(F.q[i]) + 1.);
      fphi[i] = 0.10 / (std::exp(B.q[i]) + 1.);
    }
    DecayTransitionKernel::Config cfg;
    cfg.inverse_decays     = true;
    cfg.quantum_statistics = true;
    DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
    std::vector<double> dH(NP), dl(ND), dp(ND);
    K.ComputeBackgroundDerivs(a,
                              m,
                              Gamma,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dH.data(),
                              dl.data(),
                              dp.data());

    // Band weight and node census, per parent bin, against the partition rebuilt here
    // from the two grids by a completely separate route.
    double worst_w = 0.;
    int worst_n    = 0;
    const double A = a * a * m * m;
    for (int i = 0; i < NP; ++i) {
      const double q1   = P.q[i];
      const double eps1 = std::sqrt(q1 * q1 + A);
      const double qlo = 0.5 * (eps1 - q1), qhi = qlo + q1;
      int cuts = 0;
      for (int k = 0; k < ND; ++k) {
        if (F.q[k] > qlo && F.q[k] < qhi)
          ++cuts;
        const double img = eps1 - B.q[k];
        if (img > qlo && img < qhi)
          ++cuts;
      }
      const int n = K.node_count(i);
      worst_n     = std::max(worst_n, std::abs(n - (cuts + 1)));
      double sw   = 0.;
      for (int sN = 0; sN < n; ++sN)
        sw += K.node_weight(i, sN);
      worst_w = std::fmax(worst_w, std::fabs(sw / q1 - 1.));
    }
    std::printf(
        "band partition (N_daughter=%3d): |sum wn / q1 - 1| = %.2e, node-count "
        "mismatch = %d\n",
        ND,
        worst_w,
        worst_n);
    assert(worst_w < 1e-13);  // the band is covered exactly, with nothing dropped
    assert(worst_n == 0);     // and by exactly the sub-intervals the union implies

    auto M = K.ComputeMoments(a, m, dH.data(), dl.data(), dp.data());
    return M.N_H;
  };

  const double r128  = total_rate(128);
  const double r512  = total_rate(512);
  const double drift = std::fabs(r512 / r128 - 1.);
  std::printf("total decay rate N_H, daughter grid 128 -> 512: %.3e relative\n", drift);
  assert(std::fabs(r128) > 0.);
  assert(drift < 1e-3);  // one midpoint node per cell is already converged
}

// Separate background grids (dr_bg_refine > 1), the layout the per-k perturbation
// kernels are built on. The band is partitioned against the STATE grids, because that
// is where dF is written, while f is gathered from the finer background grids -- so
// the gathered f has kinks inside a state cell that one node per cell does not
// resolve exactly. That costs O(dq^2) accuracy and NOTHING in conservation, which is
// per transition; this test pins that distinction.
void test_separate_bg_grids(bool balanced = false) {
  const int NP = 16, NST = 64, NBG = 64 * 3;  // refine = 3
  Grid P         = MakeGrid(1e-3, 1e2, NP);
  Grid FS        = MakeGrid(1e-4, 2e2, NST);
  Grid BS        = MakeGrid(1e-4, 2e2, NST);
  Grid FB        = MakeGrid(1e-4, 2e2, NBG);
  Grid BB        = MakeGrid(1e-4, 2e2, NBG);
  const double a = 0.5, m = 2.0, Gamma = 0.3;

  std::vector<double> fH(NP), fl(NBG), fphi(NBG);
  for (int i = 0; i < NP; ++i)
    fH[i] = 1. / (std::exp(P.q[i]) + 1.);
  for (int i = 0; i < NBG; ++i) {
    fl[i]   = 0.10 / (std::exp(FB.q[i]) + 1.);
    fphi[i] = 0.10 / (std::exp(BB.q[i]) + 1.);
  }

  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  // The production geometry for the chain rule, and the one place its two halves are
  // read off DIFFERENT grids: G(1-+G) comes from the band factor's BACKGROUND-grid
  // gather (the accurate estimate of f at q*), while 1/(f(1-+f)) is on the STATE grid,
  // because that is where F lives and so where Psi = F/f can be formed at all. The
  // diagonal check below is what pins that pairing -- at refine = 1 the two grids
  // coincide and it cannot see a mix-up.
  cfg.balanced_gather = balanced;
  cfg.balanced_pert   = balanced;
  cfg.lumped_loss     = !balanced;
  DecayTransitionKernel K(View(P),
                          View(FS),
                          View(BS),
                          Statistics::Fermion,
                          Statistics::Boson,
                          cfg,
                          View(FB),
                          View(BB));

  std::vector<double> dH(NP), dl(NST), dp(NST);
  K.ComputeBackgroundDerivs(a,
                            m,
                            Gamma,
                            fH.data(),
                            fl.data(),
                            fphi.data(),
                            dH.data(),
                            dl.data(),
                            dp.data());
  auto M          = K.ComputeMoments(a, m, dH.data(), dl.data(), dp.data());
  const double nf = std::fabs(M.N_H + M.N_l) / Scale(M.N_H, M.N_l);
  const double nb = std::fabs(2. * M.N_H + M.N_phi) / Scale(2. * M.N_H, M.N_phi);
  const double en = std::fabs(2. * M.E_H + 2. * M.E_l + M.E_phi - K.split_energy_residual()) /
                    Scale3(2. * M.E_H, 2. * M.E_l, M.E_phi);
  std::printf(
      "separate bg grids (refine=3): |N_H+N_l|=%.2e "
      "|2N_H+N_phi|=%.2e |energy - booked|=%.2e\n",
      nf,
      nb,
      en);
  assert(nf < 1e-12);
  assert(nb < 1e-12);

  // The PERTURBATION operator's implied diagonal against CollisionDiagonal, ON THE
  // SEPARATE-GRID LAYOUT. This is the invariant that matters for ETD: the exponential
  // step subtracts CollisionDiagonal from the operator, so if the two disagree the
  // split is inconsistent. It is checked here and not only on the single-grid layout
  // because every lumping quantity (theta, the deviation, d(theta)/df) is defined
  // against the STATE grid, while the tempting fl_gather_[n] is the BACKGROUND-grid
  // one -- a distinction that is invisible at refine = 1 and wrong above it. That is
  // the trap the lumped-loss comment in PrepareTransitions records, and one I walked
  // into while writing the limiter's Jacobian.
  {
    std::vector<double> dgH(NP), dgl(NST), dgp(NST);
    K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
    K.CollisionDiagonal(dgH.data(), dgl.data(), dgp.data());

    // ⚠ SCALE-RELATIVE, NOT POINTWISE. This used to divide by |dgl[k]| floored at
    // 1e-300, i.e. not floored at all, and the daughter grid has bins that are empty to
    // within double precision: the worst offender here sits at |dgl| = 2.6e-15 against a
    // max of 1.14, and the bin the LINEAR gather reports at 4.3e-38. Dividing round-off
    // by round-off is not a test -- it passed on clang and failed on gcc at 1.9e-10 for
    // no reason either compiler was wrong about, and the balanced gather makes it worse
    // by construction, since its chain factor is DESIGNED to be maximally sensitive
    // exactly where a bin empties.
    //
    // Flooring at 1e-12 of the operator's own scale keeps every discrepancy that could
    // matter to ETD -- which subtracts this diagonal from the operator, so what counts is
    // the deviation against the operator's magnitude, not against a bin that holds
    // nothing. The real disagreements this test exists to catch (a lumping quantity read
    // off the wrong grid) are O(1) relative and clear the floor by twelve orders.
    double dgmax = 0.;
    for (int k = 0; k < NST; ++k)
      dgmax = std::fmax(dgmax, std::fabs(dgl[k]));
    assert(dgmax > 0.);  // else every ratio below is 0/0 and the test is vacuous

    // Column k of the operator at l = 0: F = e_k, read back dF_l[k].
    double worst = 0.;
    int worst_k  = -1;
    for (int k = 0; k < NST; ++k) {
      std::vector<double> FH(NP, 0.), Fl(NST, 0.), Fp(NST, 0.);
      std::vector<double> dH2(NP), dl2(NST), dp2(NST);
      Fl[k] = 1.;
      K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
      K.ApplyPerturbationOperator(0,
                                  FH.data(),
                                  Fl.data(),
                                  Fp.data(),
                                  dH2.data(),
                                  dl2.data(),
                                  dp2.data());
      const double sc  = std::fmax(std::fabs(dgl[k]), 1e-12 * dgmax);
      const double rel = std::fabs(dl2[k] - dgl[k]) / sc;
      if (rel > worst) {
        worst   = rel;
        worst_k = k;
      }
    }
    std::printf(
        "   pert diagonal vs CollisionDiagonal on refine=3: worst rel dev = %.2e at bin %d\n",
        worst,
        worst_k);
    assert(worst < 1e-10);
  }
  assert(en < 1e-12);
  assert(std::fabs(M.N_H) > 0.);

  // And the perturbation operator on the same layout must stay conservative too.
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
  std::vector<double> FH(NP), Fl(NST, 0.), Fp(NST, 0.);
  for (int i = 0; i < NP; ++i)
    FH[i] = 0.05 * fH[i];
  std::vector<double> gH(NP), gl(NST), gp(NST);
  K.ApplyPerturbationOperator(0, FH.data(), Fl.data(), Fp.data(), gH.data(), gl.data(), gp.data());
  auto Mp         = K.ComputeMoments(a, m, gH.data(), gl.data(), gp.data());
  const double pf = std::fabs(Mp.N_H + Mp.N_l) / Scale(Mp.N_H, Mp.N_l);
  std::printf("   perturbation l=0 on the same layout: |numH+numl|=%.2e\n", pf);
  assert(pf < 1e-12);
}

// The all-multipole sweep against the per-l entry point, term by term.
//
// WHY THIS EXISTS. ApplyPerturbationOperatorAllL is ~70% of a Gamma=1e8 run and is
// the only form the production path calls (dncdm_inv_species.cpp), yet every other
// test in this file drives the per-l form -- so until now the hot function had no
// direct test at all and an optimisation could silently change the physics. The
// per-l form is the reference: it is the one the conservation, diagonal and
// equilibrium tests pin, so agreement with it inherits all of those.
//
// The tolerance is relative-to-scale rather than bit-identical ON PURPOSE. The two
// forms accumulate each (bin, l) in the same order, so today they agree exactly, but
// requiring that would forbid any reassociation in the sweep -- including hoisting a
// loop-invariant reciprocal out of the l loop, which is a legitimate optimisation
// with no physical content. 1e-12 of the array's own scale is far tighter than any
// real defect could hide under and still leaves that room.
void test_all_l_matches_per_l(bool balanced = false) {
  Grid P       = MakeGrid(1e-3, 1e2, 48);
  Grid F       = MakeGrid(1e-4, 2e2, 64);
  Grid B       = MakeGrid(1e-4, 2e2, 64);
  const int NP = 48, ND = 64;
  const int L = 17;  // l_max_dncdm_col's default, i.e. the production sweep length

  std::mt19937 rng(31337u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::uniform_real_distribution<double> Upm(-1.0, 1.0);
  std::vector<double> fH(NP), fl(ND), fphi(ND);
  for (auto& x : fH)
    x = U(rng);
  for (auto& x : fl)
    x = U(rng);
  for (auto& x : fphi)
    x = U(rng);

  const int stride = L + 1;
  std::vector<double> FH(NP * stride), Fl(ND * stride), Fphi(ND * stride);
  for (auto& x : FH)
    x = Upm(rng);
  for (auto& x : Fl)
    x = Upm(rng);
  for (auto& x : Fphi)
    x = Upm(rng);

  const double a = 0.5, m = 2.0, Gamma = 0.3;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  cfg.balanced_pert      = balanced;
  cfg.lumped_loss        = !balanced;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);

  // Reference: the per-l form, one multipole at a time into contiguous slices, then
  // transposed into the [bin][stride] layout the sweep writes.
  std::vector<double> refH(NP * stride, 0.), refl(ND * stride, 0.), refp(ND * stride, 0.);
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
  {
    std::vector<double> sFH(NP), sFl(ND), sFp(ND), sdH(NP), sdl(ND), sdp(ND);
    for (int l = 0; l <= L; ++l) {
      for (int i = 0; i < NP; ++i)
        sFH[i] = FH[(size_t) i * stride + l];
      for (int j = 0; j < ND; ++j)
        sFl[j] = Fl[(size_t) j * stride + l];
      for (int k = 0; k < ND; ++k)
        sFp[k] = Fphi[(size_t) k * stride + l];
      K.ApplyPerturbationOperator(l,
                                  sFH.data(),
                                  sFl.data(),
                                  sFp.data(),
                                  sdH.data(),
                                  sdl.data(),
                                  sdp.data());
      for (int i = 0; i < NP; ++i)
        refH[(size_t) i * stride + l] = sdH[i];
      for (int j = 0; j < ND; ++j)
        refl[(size_t) j * stride + l] = sdl[j];
      for (int k = 0; k < ND; ++k)
        refp[(size_t) k * stride + l] = sdp[k];
    }
  }
  const double ref_resid_e = K.split_energy_residual_pert();
  const double ref_resid_m = K.split_momentum_residual_pert();

  // The sweep. PrepareTransitions again: it resets next_l_, which the sweep requires.
  std::vector<double> gH(NP * stride, 0.), gl(ND * stride, 0.), gp(ND * stride, 0.);
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
  K.ApplyPerturbationOperatorAllL(L,
                                  stride,
                                  FH.data(),
                                  Fl.data(),
                                  Fphi.data(),
                                  gH.data(),
                                  gl.data(),
                                  gp.data());

  auto worst_rel = [](const std::vector<double>& x, const std::vector<double>& y) {
    double scale = 0.;
    for (double v : y)
      scale = std::fmax(scale, std::fabs(v));
    scale    = std::fmax(scale, 1e-300);
    double w = 0.;
    for (size_t i = 0; i < x.size(); ++i)
      w = std::fmax(w, std::fabs(x[i] - y[i]) / scale);
    return w;
  };
  const double wH = worst_rel(gH, refH);
  const double wl = worst_rel(gl, refl);
  const double wp = worst_rel(gp, refp);
  // The booked split residuals are part of the contract too: they feed the energy
  // identity in test_perturbation_operator, so a sweep that got them wrong would
  // look conservative here and be wrong there.
  //
  // Scaled against the operator's own l = 0 ENERGY, not against the residual itself:
  // with the exact split (lumped_loss off, which is what balanced_gather brings) the
  // misplacement is identically zero, so both implementations book round-off and a
  // residual-relative comparison would be dividing 1e-17 by 1e-17. The identity the
  // residual serves is 2eH + 2el + ephi == booked, so the energy is its natural scale.
  double e_scale = 0.;
  for (int j = 0; j < ND; ++j)
    e_scale += std::fabs(F.dq[j] * F.q[j] * F.q[j] * F.q[j] * refl[(size_t) j * stride]);
  const double we = std::fabs(K.split_energy_residual_pert() - ref_resid_e) /
                    Scale3(K.split_energy_residual_pert(), ref_resid_e, e_scale);
  const double wm = std::fabs(K.split_momentum_residual_pert() - ref_resid_m) /
                    Scale3(K.split_momentum_residual_pert(), ref_resid_m, e_scale);
  std::printf(
      "all-l vs per-l: dF_H=%.2e dF_l=%.2e dF_phi=%.2e "
      "resid_e=%.2e resid_m=%.2e\n",
      wH,
      wl,
      wp,
      we,
      wm);
  assert(wH < 1e-12);
  assert(wl < 1e-12);
  assert(wp < 1e-12);
  assert(we < 1e-12);
  assert(wm < 1e-12);
  // Guard against a trivially-zero operator "agreeing" with itself.
  double mag = 0.;
  for (double v : refl)
    mag = std::fmax(mag, std::fabs(v));
  assert(mag > 0.);
}

// THE OPERATOR IS A JACOBIAN, AND OF WHICH BAND FACTOR.
//
// ApplyPerturbationOperator at l = 0 (every Legendre factor is 1) is a DIRECTIONAL
// DERIVATIVE of ComputeBackgroundDerivs: the same network, linearised. So it can be
// pinned against a central difference of the background RHS along the perturbation --
// two entirely independent code paths, one of them not differentiated by hand at all.
//
// This is the test the balanced gather did not have. test_perturbation_diagonal pins the
// operator against CollisionDiagonal, but both were written from the same derivation, so
// a wrong gather in BOTH agrees with itself. Against finite differences it does not:
// with balanced_gather on and the chain rule off, the operator gathers F linearly and is
// the exact Jacobian of a DIFFERENT discrete model than the one the background
// integrates -- measured 1.1e-1 here, against 8.8e-12 on the linear path.
//
// The perturbation is seeded as F = f * Psi with Psi a smooth O(1) function, which is
// what the hierarchy delivers and what makes f + h F stay positive for the difference.
// The step is scanned rather than fixed: the quotient is a cancellation, so it has a
// truncation branch and a round-off branch and the minimum over the scan is the signal.
void test_perturbation_jacobian(bool balanced, bool chain) {
  const int N = 48;
  Grid P      = MakeGrid(1e-3, 1e2, N);
  Grid F      = MakeGrid(1e-4, 2e2, N);
  Grid B      = MakeGrid(1e-4, 2e2, N);

  std::mt19937 rng(98765u);
  std::uniform_real_distribution<double> U(0.01, 0.9);
  std::vector<double> fH(N), fl(N), fphi(N), FH(N), Fl(N), Fphi(N);
  for (int i = 0; i < N; ++i) {
    fH[i]   = U(rng);
    fl[i]   = U(rng);
    fphi[i] = U(rng);
  }
  for (int i = 0; i < N; ++i) {
    FH[i]   = fH[i] * 0.37 * std::sin(1.3 * std::log(P.q[i]));
    Fl[i]   = fl[i] * 0.29 * std::cos(0.9 * std::log(F.q[i]));
    Fphi[i] = fphi[i] * 0.41 * std::sin(0.7 * std::log(B.q[i]) + 0.3);
  }

  const double a = 0.5, m = 2.0, Gamma = 0.3;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = balanced;
  cfg.balanced_pert      = chain;
  cfg.lumped_loss        = !balanced;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);

  std::vector<double> oH(N), ol(N), ophi(N);
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
  K.ApplyPerturbationOperator(0,
                              FH.data(),
                              Fl.data(),
                              Fphi.data(),
                              oH.data(),
                              ol.data(),
                              ophi.data());

  auto rel = [](const std::vector<double>& x, const std::vector<double>& y) {
    double num = 0., den = 0.;
    for (size_t i = 0; i < x.size(); ++i) {
      num = std::fmax(num, std::fabs(x[i] - y[i]));
      den = std::fmax(den, std::fabs(y[i]));
    }
    return den > 0. ? num / den : num;
  };

  std::vector<double> pH(N), pl(N), pphi(N);
  std::vector<double> aH(N), al(N), aphi(N), bH(N), bl(N), bphi(N), gH(N), gl(N), gphi(N);
  double best = 1e300, best_h = 0.;
  for (double h : {1e-3, 3e-4, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6}) {
    for (int i = 0; i < N; ++i) {
      pH[i]   = fH[i] + h * FH[i];
      pl[i]   = fl[i] + h * Fl[i];
      pphi[i] = fphi[i] + h * Fphi[i];
    }
    K.ComputeBackgroundDerivs(a,
                              m,
                              Gamma,
                              pH.data(),
                              pl.data(),
                              pphi.data(),
                              aH.data(),
                              al.data(),
                              aphi.data());
    for (int i = 0; i < N; ++i) {
      pH[i]   = fH[i] - h * FH[i];
      pl[i]   = fl[i] - h * Fl[i];
      pphi[i] = fphi[i] - h * Fphi[i];
    }
    K.ComputeBackgroundDerivs(a,
                              m,
                              Gamma,
                              pH.data(),
                              pl.data(),
                              pphi.data(),
                              bH.data(),
                              bl.data(),
                              bphi.data());
    for (int i = 0; i < N; ++i) {
      gH[i]   = (aH[i] - bH[i]) / (2 * h);
      gl[i]   = (al[i] - bl[i]) / (2 * h);
      gphi[i] = (aphi[i] - bphi[i]) / (2 * h);
    }
    best = std::fmin(best, std::fmax(rel(oH, gH), std::fmax(rel(ol, gl), rel(ophi, gphi))));
    if (best == std::fmax(rel(oH, gH), std::fmax(rel(ol, gl), rel(ophi, gphi))))
      best_h = h;
  }
  std::printf("pert Jacobian vs finite differences (%s gather%s): %.2e at h=%.0e\n",
              balanced ? "balanced" : "linear",
              chain ? " + chain rule" : "",
              best,
              best_h);
  // Guard against a vacuous pass: the operator must not be identically zero.
  double mag = 0.;
  for (int i = 0; i < N; ++i)
    mag = std::fmax(mag, std::fabs(ol[i]));
  assert(mag > 0.);

  if (balanced && !chain) {
    // CONTROL. Without the chain rule the operator differentiates the linear gather
    // while the background integrates the balanced one. If this ever stops failing, the
    // two assertions below have stopped measuring anything.
    assert(best > 1e-3);
  }
  else {
    assert(best < 1e-9);
  }
}

// The chain factor diverges as a bin empties -- dG/df_e = w_e (f_other/f_e)^{w_other} in
// the dilute limit -- which is a real property of a geometric mean and the reason this
// was left unshipped. ChainFactor caps it. This pins that the cap holds on the state that
// provokes it (an injection front: f falling off a cliff mid-grid) and, just as
// importantly, that the operator stays FINITE and CONSERVATIVE there.
void test_chain_factor_cap() {
  const int N = 40;
  Grid P      = MakeGrid(1e-3, 1e2, N);
  Grid F      = MakeGrid(1e-4, 2e2, N);
  Grid B      = MakeGrid(1e-4, 2e2, N);

  // Injection front: the daughters are filled up to a cut and empty above it, spanning
  // ~80 decades across two bins -- the geometry that makes G/f_e ~ 1e40.
  std::vector<double> fH(N), fl(N), fphi(N), FH(N), Fl(N), Fphi(N);
  for (int i = 0; i < N; ++i) {
    fH[i]   = 0.4 * std::exp(-P.q[i] / 3.0);
    fl[i]   = (i < N / 2) ? 0.3 * std::exp(-F.q[i]) : 1e-90 * std::exp(-(double) (i - N / 2) * 4.);
    fphi[i] = (i < N / 2) ? 0.3 * std::exp(-B.q[i]) : 1e-90 * std::exp(-(double) (i - N / 2) * 4.);
  }
  // Psi ~ O(1) everywhere, so F inherits the same cliff: the physical statement that an
  // empty bin has an empty perturbation.
  for (int i = 0; i < N; ++i) {
    FH[i]   = 0.3 * fH[i];
    Fl[i]   = -0.2 * fl[i];
    Fphi[i] = 0.5 * fphi[i];
  }

  const double a = 0.5, m = 2.0, Gamma = 0.3, A = a * a * m * m;
  DecayTransitionKernel::Config cfg;
  cfg.inverse_decays     = true;
  cfg.quantum_statistics = true;
  cfg.balanced_gather    = true;
  cfg.balanced_pert      = true;
  cfg.lumped_loss        = false;
  DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
  K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());

  std::vector<double> dH(N), dl(N), dphi(N);
  K.ApplyPerturbationOperator(0,
                              FH.data(),
                              Fl.data(),
                              Fphi.data(),
                              dH.data(),
                              dl.data(),
                              dphi.data());

  double worst = 0.;
  for (int i = 0; i < N; ++i) {
    assert(std::isfinite(dH[i]) && std::isfinite(dl[i]) && std::isfinite(dphi[i]));
    worst = std::fmax(worst, std::fabs(dl[i]));
  }
  // The l = 0 number identities do not care what the gather did -- all three legs are
  // built from the same node source terms -- so they must still be exact ON THIS STATE.
  // That is what says the cap did not quietly break conservation to buy finiteness.
  double numH = 0, numl = 0, numphi = 0;
  for (int i = 0; i < N; ++i) {
    const double q = P.q[i], w = P.dq[i] * q * q;
    (void) A;
    numH += w * dH[i];
  }
  for (int j = 0; j < N; ++j)
    numl += F.dq[j] * F.q[j] * F.q[j] * dl[j];
  for (int k = 0; k < N; ++k)
    numphi += B.dq[k] * B.q[k] * B.q[k] * dphi[k];
  const double nf = std::fabs(numH + numl) / Scale(numH, numl);
  const double nb = std::fabs(2. * numH + numphi) / Scale(2. * numH, numphi);
  std::printf(
      "chain cap on a front state: max|dF_l| = %.2e, |numH+numl| = %.2e, "
      "|2numH+numphi| = %.2e\n",
      worst,
      nf,
      nb);
  assert(nf < 1e-11);
  assert(nb < 1e-11);
  assert(worst > 0.);  // else the finiteness assertions above are vacuous
}

// REGRESSION: the cap must reach the BACKGROUND diagonal, not just the operator.
//
// The cap used to be gated on chain_pert_ (= balanced_gather && balanced_pert) while
// ChainFactor is CONSUMED under chain_diag_, which is a superset: a background kernel
// has separate_bg_grids_ == false at the dr_bg_refine = 1 production setting, so it
// carries the chain rule on balanced_gather ALONE. Plain `balanced_gather` therefore ran
// CollisionDiagonal with inv_cap = 0 and reached diag = 6.8e84 on an emptying bin, which
// etd cannot integrate -- hpc_ratchet's cmb_G05.0_m0.3_x_bal_q2x died in 3 s, and a
// dr_N_q scan at that cell was clean at 83..103 EXCEPT exactly 93, a knife-edge on where
// the deposit lands. test_chain_factor_cap above never caught it because it sets
// balanced_pert = true, i.e. the one configuration in which the cap was already applied.
//
// The invariant asserted here is exact rather than a magnitude bound: balanced_pert
// changes the PERTURBATION gather and nothing the background diagonal reads, so the two
// kernels must return the SAME diagonal. Before the fix they differed by ~1e36.
void test_chain_cap_reaches_background_diagonal() {
  const int N = 40;
  Grid P      = MakeGrid(1e-3, 1e2, N);
  Grid F      = MakeGrid(1e-4, 2e2, N);
  Grid B      = MakeGrid(1e-4, 2e2, N);

  // Same injection front as test_chain_factor_cap: daughters filled to a cut and empty
  // above it, which is the geometry that drives X = G(1-+G)/(f_e(1-+f_e)) to ~1e40.
  std::vector<double> fH(N), fl(N), fphi(N);
  for (int i = 0; i < N; ++i) {
    fH[i]   = 0.4 * std::exp(-P.q[i] / 3.0);
    fl[i]   = (i < N / 2) ? 0.3 * std::exp(-F.q[i]) : 1e-90 * std::exp(-(double) (i - N / 2) * 4.);
    fphi[i] = (i < N / 2) ? 0.3 * std::exp(-B.q[i]) : 1e-90 * std::exp(-(double) (i - N / 2) * 4.);
  }
  const double a = 0.5, m = 2.0, Gamma = 0.3;

  auto diagonal_with = [&](bool balanced_pert,
                           std::vector<double>& dH,
                           std::vector<double>& dl,
                           std::vector<double>& dphi) {
    DecayTransitionKernel::Config cfg;
    cfg.inverse_decays     = true;
    cfg.quantum_statistics = true;
    cfg.balanced_gather    = true;
    cfg.balanced_pert      = balanced_pert;
    cfg.lumped_loss        = false;
    DecayTransitionKernel K(View(P), View(F), View(B), Statistics::Fermion, Statistics::Boson, cfg);
    K.PrepareTransitions(a, m, Gamma, fH.data(), fl.data(), fphi.data());
    K.CollisionDiagonal(dH.data(), dl.data(), dphi.data());
    return K.chain_factor_max();
  };

  std::vector<double> dH0(N), dl0(N), dphi0(N), dH1(N), dl1(N), dphi1(N);
  const double raw = diagonal_with(false, dH0, dl0, dphi0);  // `bal`: the shipped config
  diagonal_with(true, dH1, dl1, dphi1);                      // `chn`: always capped

  // The raw worst X really does reach the runaway regime on this state, else every
  // assertion below is vacuous and the test would keep passing through a reintroduction.
  assert(raw > 1e20);

  double worst = 0.;
  for (int i = 0; i < N; ++i) {
    assert(std::isfinite(dH0[i]) && std::isfinite(dl0[i]) && std::isfinite(dphi0[i]));
    worst =
        std::fmax(worst,
                  std::fmax(std::fabs(dH0[i] - dH1[i]),
                            std::fmax(std::fabs(dl0[i] - dl1[i]), std::fabs(dphi0[i] - dphi1[i]))));
  }
  double scale = 0.;
  for (int i = 0; i < N; ++i)
    scale = std::fmax(scale,
                      std::fmax(std::fabs(dH1[i]),
                                std::fmax(std::fabs(dl1[i]), std::fabs(dphi1[i]))));
  std::printf(
      "   background diagonal, balanced_pert off vs on: worst abs dev = %.2e "
      "(scale %.2e, raw chain factor %.2e)\n",
      worst,
      scale,
      raw);
  assert(scale > 0.);  // a trivially-zero diagonal would agree for the wrong reason
  assert(worst <= 1e-12 * scale);
}

}  // namespace

int main() {
  test_deposit_moments();
  test_band_partition();
  test_separate_bg_grids();
  test_separate_bg_grids(true);
  test_conservation();
  test_conservation(true);
  test_equilibrium();
  test_rhs_continuity();
  test_perturbation_operator();
  test_perturbation_operator(true);  // balanced gather + chain rule: same identities
  test_all_l_matches_per_l();
  test_all_l_matches_per_l(true);
  test_perturbation_jacobian(false, false);  // linear gather: already a Jacobian
  test_perturbation_jacobian(true, false);   // control: balanced band factor, linear gather
  test_perturbation_jacobian(true, true);    // the claim
  test_chain_factor_cap();
  test_chain_cap_reaches_background_diagonal();
  test_daughter_positivity();              // shipped: linear gather, lumped loss
  test_daughter_positivity(false, false);  // control: the defect lumping exists for
  test_daughter_positivity(true, true);    // balanced gather, lumping still on
  test_daughter_positivity(true, false);   // balanced gather, EXACT split -- the claim
  test_collision_diagonal();
  test_collision_diagonal(true);
  test_perturbation_diagonal();
  test_perturbation_diagonal(true);
  test_mode_reduction();
  test_band_bound();
  std::printf("decay kernel test passed\n");
  return 0;
}
