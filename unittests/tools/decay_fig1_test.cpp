// Hannestad & Madsen 1997 (astro-ph/9702125) fig 1 through DecayHarness, and the
// balanced-gather stall it exposed (issue #389, step 2).
//
// Two pieces:
//   1. The kernel's floored log-odds gather (Config::gather_floor) is Lipschitz in an
//      empty daughter bin, where the unfloored one goes like f^(1-w). That Hoelder
//      dependence on noise-level bins is what stalled ndf15 at m_H = 2 MeV,
//      tau = 0.2 s; the control runs the same state with the floor off.
//   2. Fig 1 itself with the harness's DEFAULT scheme, against the paper's digitized
//      curves (scripts/data/hannestad_madsen_1997/fig1.csv; the needed columns are
//      embedded below). Every run is under a watchdog, so a returning stall fails
//      instead of hanging CI.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <future>
#include <vector>

#include "decay_harness.h"
#include "decay_transition_kernel.h"

namespace {

// ── 1. the gather's regularity in an empty bin ────────────────────────────────────

struct Grid {
  std::vector<double> q, dq;
  GridView View() const {
    return GridView{q.data(), dq.data(), static_cast<int>(q.size())};
  }
};

Grid LogGrid(double qmin, double qmax, int n) {
  Grid g;
  g.q.resize(n);
  g.dq.resize(n);
  const double h = std::log(qmax / qmin) / (n - 1);
  for (int i = 0; i < n; ++i) {
    g.q[i] = qmin * std::exp(i * h);
  }
  g.dq[0] = 0.5 * (g.q[1] - g.q[0]);
  for (int i = 1; i < n - 1; ++i) {
    g.dq[i] = 0.5 * (g.q[i + 1] - g.q[i - 1]);
  }
  g.dq[n - 1] = 0.5 * (g.q[n - 1] - g.q[n - 2]);
  return g;
}

// The stall's shape, built by hand: a light parent remnant whose emission band reaches a
// boson injection front, i.e. f_phi alive below q = 10 and exactly empty above. The
// parents' right-hand side is then differentiated with respect to the first empty bin,
// with steps of 1e-27, 1e-24 and 1e-21 -- all below the harness's abstol of 1e-15, where
// ndf15 neither controls the occupation nor can keep it from changing sign.
//
// Returns max over parent bins of |slope(eps) / slope(1e-21)| for eps = 1e-27, 1e-24, with
// slope the forward difference. A Lipschitz gather gives 1; a gather going like f^(1-w)
// gives up to (1e6)^w, i.e. orders of magnitude.
double EmptyBinSlopeRatio(double floor) {
  const Grid P                      = LogGrid(0.1, 20., 16);
  const Grid D                      = LogGrid(1e-2, 60., 40);
  DecayTransitionKernel::Config cfg = DecayHarness::Config::DefaultKernel();
  cfg.fermion_mass_ratio            = 0.1;
  cfg.boson_mass_ratio              = 0.0025;
  cfg.gather_floor                  = floor;
  const DecayTransitionKernel K(P.View(),
                                D.View(),
                                D.View(),
                                Statistics::Fermion,
                                Statistics::Boson,
                                cfg);

  const double a = 7.;
  const double m = 2.;  // a m = 14: T_nu = 0.14 MeV at m_H = 2 MeV, the stall's epoch
  std::vector<double> fH(P.q.size());
  std::vector<double> fl(D.q.size());
  std::vector<double> fphi(D.q.size());
  for (size_t i = 0; i < P.q.size(); ++i) {
    fH[i] = 1e-12 * std::exp(-std::sqrt(P.q[i] * P.q[i] + a * a * m * m) / 14.);
  }
  int dead = -1;
  for (size_t j = 0; j < D.q.size(); ++j) {
    fl[j]   = 1. / (std::exp(D.q[j]) + 1.);
    fphi[j] = D.q[j] < 10. ? 1. / std::expm1(D.q[j]) : 0.;
    if (dead < 0 && fphi[j] == 0.) {
      dead = static_cast<int>(j);
    }
  }
  assert(dead > 0);

  std::vector<double> dH0(P.q.size());
  std::vector<double> dH1(P.q.size());
  std::vector<double> dl(D.q.size());
  std::vector<double> dp(D.q.size());
  const auto slopes = [&](double eps) {
    K.ComputeBackgroundDerivs(a,
                              m,
                              1.,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dH0.data(),
                              dl.data(),
                              dp.data());
    fphi[dead] = eps;
    K.ComputeBackgroundDerivs(a,
                              m,
                              1.,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dH1.data(),
                              dl.data(),
                              dp.data());
    fphi[dead] = 0.;
    std::vector<double> s(P.q.size());
    for (size_t i = 0; i < P.q.size(); ++i) {
      s[i] = (dH1[i] - dH0[i]) / eps;
    }
    return s;
  };
  const auto fine   = slopes(1e-27);
  const auto mid    = slopes(1e-24);
  const auto coarse = slopes(1e-21);
  double worst      = 0.;
  int reached       = 0;
  for (size_t i = 0; i < P.q.size(); ++i) {
    if (coarse[i] == 0. && fine[i] == 0.) {
      continue;  // this parent's band does not reach the empty bin
    }
    ++reached;
    worst = std::fmax(worst, std::fabs(fine[i] / coarse[i]));
    worst = std::fmax(worst, std::fabs(mid[i] / coarse[i]));
  }
  // Guard, not skip: the state must actually put the empty bin inside some band.
  assert(reached > 0);
  return worst;
}

void test_floor_makes_gather_lipschitz() {
  const double with    = EmptyBinSlopeRatio(1e-18);
  const double without = EmptyBinSlopeRatio(0.);
  std::printf(
      "empty-bin slope ratio (1e-27 vs 1e-21 step): floor 1e-18 -> %.6f, "
      "no floor (control) -> %.3e\n",
      with,
      without);
  // 1.0018 measured: the forward difference at a 1e-27 step resolves the slope to ~1e-3,
  // against 7e5 for the control.
  assert(std::fabs(with - 1.) < 1e-2);
  assert(without > 1e2);
}

// The floored chain factor is the exact derivative of the floored gather:
// CollisionDiagonal against a central difference of the background right-hand side, on a
// resolved thermal state well above the floor (so the chain cap does not bind) but with a
// floor big enough to matter if EtaSlopeFloored were wrong (d = 1e-3).
void test_floored_diagonal_matches_jacobian() {
  const Grid P                      = LogGrid(0.1, 20., 16);
  const Grid D                      = LogGrid(1e-2, 40., 48);
  DecayTransitionKernel::Config cfg = DecayHarness::Config::DefaultKernel();
  cfg.fermion_mass_ratio            = 0.1;
  cfg.boson_mass_ratio              = 0.0025;
  cfg.gather_floor                  = 1e-3;
  const DecayTransitionKernel K(P.View(),
                                D.View(),
                                D.View(),
                                Statistics::Fermion,
                                Statistics::Boson,
                                cfg);
  const double a = 1.;
  const double m = 2.;
  std::vector<double> fH(P.q.size());
  std::vector<double> fl(D.q.size());
  std::vector<double> fphi(D.q.size());
  for (size_t i = 0; i < P.q.size(); ++i) {
    fH[i] = 0.5 / (std::exp(std::sqrt(P.q[i] * P.q[i] + a * a * m * m) / 1.3) + 1.);
  }
  for (size_t j = 0; j < D.q.size(); ++j) {
    fl[j]   = 0.7 / (std::exp(D.q[j] / 1.1) + 1.);
    fphi[j] = 0.3 / std::expm1(D.q[j] / 0.9 + 0.05);
  }
  std::vector<double> dH(P.q.size());
  std::vector<double> dl(D.q.size());
  std::vector<double> dp(D.q.size());
  std::vector<double> diagH(P.q.size());
  std::vector<double> diagl(D.q.size());
  std::vector<double> diagp(D.q.size());
  K.ComputeBackgroundDerivs(a,
                            m,
                            1.,
                            fH.data(),
                            fl.data(),
                            fphi.data(),
                            dH.data(),
                            dl.data(),
                            dp.data());
  K.CollisionDiagonal(diagH.data(), diagl.data(), diagp.data());
  double worst = 0.;
  double scale = 0.;
  for (size_t j = 0; j < D.q.size(); ++j) {
    scale = std::fmax(scale, std::fmax(std::fabs(diagl[j]), std::fabs(diagp[j])));
  }
  const auto fd = [&](std::vector<double>& f, size_t j, const std::vector<double>& out) {
    const double f0 = f[j];
    const double h  = 1e-6 * f0;
    f[j]            = f0 + h;
    K.ComputeBackgroundDerivs(a,
                              m,
                              1.,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dH.data(),
                              dl.data(),
                              dp.data());
    const double up = out[j];
    f[j]            = f0 - h;
    K.ComputeBackgroundDerivs(a,
                              m,
                              1.,
                              fH.data(),
                              fl.data(),
                              fphi.data(),
                              dH.data(),
                              dl.data(),
                              dp.data());
    const double dn = out[j];
    f[j]            = f0;
    return (up - dn) / (2. * h);
  };
  for (size_t j = 0; j < D.q.size(); ++j) {
    if (fl[j] > 1e-6) {
      worst = std::fmax(worst, std::fabs(fd(fl, j, dl) - diagl[j]) / scale);
    }
    if (fphi[j] > 1e-6) {
      worst = std::fmax(worst, std::fabs(fd(fphi, j, dp) - diagp[j]) / scale);
    }
  }
  std::printf("floored CollisionDiagonal vs central difference: %.2e (of max |diag|)\n", worst);
  assert(worst < 1e-6);
}

// ── 2. fig 1 against the digitized paper ──────────────────────────────────────────

// q = p/T_nu (labelled p/T_gamma in the paper) and q^2 f_F for tau = 0.2, 1, 20 s,
// verbatim from scripts/data/hannestad_madsen_1997/fig1.csv (q in [0.1, 12.75]; the
// q = 13 row is empty there).
constexpr double kQ[]     = {0.1,  0.25,  0.5,  0.75,  1,    1.25,  1.5,  1.75, 2,     2.25, 2.5,
                             2.75, 3,     3.25, 3.5,   3.75, 4,     4.25, 4.5,  4.75,  5,    5.25,
                             5.5,  5.75,  6,    6.25,  6.5,  6.75,  7,    7.25, 7.5,   7.75, 8,
                             8.25, 8.5,   8.75, 9,     9.25, 9.5,   9.75, 10,   10.25, 10.5, 10.75,
                             11,   11.25, 11.5, 11.75, 12,   12.25, 12.5, 12.75};
constexpr double kTau02[] = {0.009, 0.051, 0.179, 0.367, 0.585, 0.800, 0.976, 1.094, 1.147,
                             1.140, 1.087, 1.004, 0.903, 0.796, 0.689, 0.590, 0.499, 0.418,
                             0.348, 0.288, 0.237, 0.194, 0.158, 0.128, 0.104, 0.084, 0.067,
                             0.054, 0.043, 0.035, 0.028, 0.022, 0.018, 0.014, 0.011, 0.009,
                             0.007, 0.006, 0.004, 0.004, 0.003, 0.002, 0.002, 0.001, 0.001,
                             0.001, 0.001, 0.001, 0.001, 0.000, 0.000, 0.000};
constexpr double kTau1[]  = {0.009, 0.048, 0.166, 0.324, 0.496, 0.662, 0.805, 0.913, 0.982,
                             1.011, 1.004, 0.966, 0.905, 0.830, 0.748, 0.663, 0.580, 0.502,
                             0.429, 0.365, 0.309, 0.259, 0.217, 0.180, 0.149, 0.124, 0.102,
                             0.084, 0.069, 0.057, 0.046, 0.038, 0.031, 0.025, 0.021, 0.017,
                             0.014, 0.011, 0.009, 0.007, 0.006, 0.005, 0.004, 0.003, 0.003,
                             0.002, 0.002, 0.001, 0.001, 0.001, 0.001, 0.001};
constexpr double kTau20[] = {0.005, 0.030, 0.105, 0.203, 0.306, 0.402, 0.484, 0.548, 0.593,
                             0.620, 0.632, 0.631, 0.621, 0.603, 0.581, 0.555, 0.528, 0.500,
                             0.471, 0.444, 0.417, 0.391, 0.366, 0.343, 0.320, 0.299, 0.279,
                             0.259, 0.241, 0.224, 0.207, 0.192, 0.177, 0.163, 0.150, 0.137,
                             0.126, 0.115, 0.104, 0.095, 0.086, 0.077, 0.070, 0.062, 0.056,
                             0.050, 0.044, 0.039, 0.035, 0.031, 0.027, 0.024};
constexpr int kNQ         = sizeof(kQ) / sizeof(kQ[0]);
static_assert(sizeof(kTau02) / sizeof(double) == kNQ);
static_assert(sizeof(kTau1) / sizeof(double) == kNQ);
static_assert(sizeof(kTau20) / sizeof(double) == kNQ);

struct Fig1Run {
  double max_dy = 0.;  // max |Delta(q^2 f_F)| at the digitized points
  double peak   = 0.;  // our peak height of q^2 f_F
  double drift  = 0.;  // worst of the two number invariants' relative drift
  double energy = 0.;  // energy budget defect
  double secs   = 0.;
};

// One fig 1 run at the paper's parameters, with the harness's default kernel scheme and
// grids. Run on a worker thread under a watchdog: a stall is the failure this test exists
// to catch, and a hung CI job reports nothing.
Fig1Run RunFig1(double tau_s, const double* paper) {
  DecayHarness::Config c;
  c.m_H_mev         = 2.;
  c.tau_s           = tau_s;
  c.r_l             = 0.1;
  c.r_phi           = 0.0025;
  c.T_start_mev     = 5.;
  c.n_eff_spectator = 1.;
  // The scheme under test is the harness default; only the floor is asserted so a change
  // of default cannot silently move this test off what it guards.
  assert(c.kernel.balanced_gather && c.kernel.inverse_decays && c.kernel.quantum_statistics);
  assert(c.kernel.gather_floor > 0.);

  const auto t0            = std::chrono::steady_clock::now();
  auto job                 = std::async(std::launch::async, [c, paper]() {
    DecayHarness h(c);
    std::vector<double> q;
    std::vector<double> f;
    h.Run({}, [&](const DecayHarness::Snapshot& s) {
      q = h.daughter_q();
      f.assign(s.f_l, s.f_l + q.size());
    });
    Fig1Run r;
    for (size_t j = 0; j < q.size(); ++j) {
      r.peak = std::fmax(r.peak, q[j] * q[j] * f[j]);
    }
    for (int k = 0; k < kNQ; ++k) {
      const double x = kQ[k];
      const auto it  = std::upper_bound(q.begin(), q.end(), x);
      assert(it != q.begin() && it != q.end());
      const size_t hi = static_cast<size_t>(it - q.begin());
      const size_t lo = hi - 1;
      // ln f linear in q between log-spaced nodes, as scripts/hannestad_fig1.py does.
      assert(f[lo] > 0. && f[hi] > 0.);
      const double t  = (x - q[lo]) / (q[hi] - q[lo]);
      const double fx = std::exp((1. - t) * std::log(f[lo]) + t * std::log(f[hi]));
      r.max_dy        = std::fmax(r.max_dy, std::fabs(x * x * fx - paper[k]));
    }
    const auto& d = h.diagnostics();
    r.drift       = std::fmax(std::fabs(d.number_defect_l), std::fabs(d.number_defect_phi));
    r.energy      = d.energy_defect;
    return r;
  });
  constexpr int kWatchdogS = 120;  // the stall never finished in 600 s; a run takes ~3 s
  if (job.wait_for(std::chrono::seconds(kWatchdogS)) != std::future_status::ready) {
    std::fprintf(stderr,
                 "FAIL: fig 1 at tau = %g s did not finish in %d s -- the balanced-gather "
                 "stall is back\n",
                 tau_s,
                 kWatchdogS);
    std::fflush(stderr);
    std::_Exit(1);
  }
  Fig1Run r = job.get();
  r.secs    = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf(
      "fig 1, tau = %4g s: max |d(q^2 f)| = %.4f, peak %.3f, number drift %.1e, "
      "energy defect %.1e, %.1f s\n",
      tau_s,
      r.max_dy,
      r.peak,
      r.drift,
      r.energy,
      r.secs);
  assert(r.drift < 1e-10);
  assert(std::fabs(r.energy) < 1e-5);
  return r;
}

// THE LIFETIMES. The caption says tau_H = 0.2, 1, 20 s. The curves are instead those of
// the paper's own lifetime parameter alpha = (m_H/3 MeV)^2 (tau/1 s) = 0.1, 1, 10, i.e.
// tau = 0.225, 2.25, 22.5 s at m_H = 2 MeV: one rule, three curves, nothing fitted. The
// caption's 0.2 and 20 s are those values rounded; its "1 s" is not, and misses the
// dotted curve by 0.113 at every resolution. Controls that fail: the caption's tau
// (0.0036 / 0.112 / 0.020) and alpha with a 3 MeV mass scale (tau = 0.1, 1, 10 s: 0.020 /
// 0.112 / 0.145).
//
// Bounds, measured at this grid (parent 48, daughter 96): alpha = 0.1, 1, 10 give 0.0024
// / 0.0038 / 0.0026 (<= 0.0019 at 2x), inside the digitization's own error (p95 0.007).
// Asserted at 0.01. An 11% change of lifetime breaks that for alpha = 1 (tau = 2.0 / 2.5
// s: 0.022 / 0.015) and alpha = 10 (20 / 25 s: 0.022 / 0.014); the relativistic alpha =
// 0.1 curve is softer (0.2 / 0.25 / 0.3 s: 0.003 / 0.006 / 0.014). tau = 1 s against the
// dotted curve is asserted ABOVE 0.05, so the test also pins which reading it encodes.
void test_fig1() {
  const Fig1Run a01 = RunFig1(0.225, kTau02);
  assert(a01.max_dy <= 0.01);
  const Fig1Run a1 = RunFig1(2.25, kTau1);
  assert(a1.max_dy <= 0.01);
  const Fig1Run a10 = RunFig1(22.5, kTau20);
  assert(a10.max_dy <= 0.01);
  const Fig1Run caption = RunFig1(1., kTau1);
  assert(caption.max_dy > 0.05);
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);  // keep the measurements when an assert fires
  test_floor_makes_gather_lipschitz();
  test_floored_diagonal_matches_jacobian();
  test_fig1();
  std::printf("decay_fig1_test: all passed\n");
  return 0;
}
