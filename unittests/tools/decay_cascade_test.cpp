// Unit tests for DecayCascade (tools/decay_harness): Hannestad & Madsen's three-level
// cascade nu_tau -> nu_mu + phi, nu_mu -> nu_e + phi, two DecayTransitionKernel channels
// on one nu_mu grid and one phi grid (issue #389, step 3). Coarse grids: the comparison
// with the paper's figures is on issue #389; this pins the bookkeeping.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

#include "common.h"  // _PI_
#include "decay_harness.h"

namespace {

using C = DecayCascade;

C::Config Coarse() {
  C::Config c;
  c.m_tau_mev = 2.;
  c.m_mu_mev  = 0.5;
  c.tau_tau_s = 0.2025;  // the lasing model of figs 5-7
  c.tau_mu_s  = 1.08;
  c.T_end_mev = 0.02;
  c.tau_n     = 16;
  c.mu_n      = 36;
  c.e_n       = 36;
  c.phi_n     = 48;
  return c;
}

// THE THREE INVARIANTS, AT EVERY SAMPLE OF A WHOLE LASING RUN. Each kernel's identities
// are exact per right-hand side and the network only sums them, so
//
//     2N_tau + 2N_mu + 2N_e,    N_phi + 2N_tau - 2N_e,    E_c - W
//
// may drift only by the integrator's error. Checked along the run, not only at its
// end, so a defect that cancels by the end cannot hide. Then the end state: both
// parents are gone, so N_phi = 4 N_tau,0 + 2 N_mu,0 (phi starts empty), and Pauli
// blocking holds both fermion daughters at or below 1 while phi lases far above it.
void test_invariants_along_run() {
  C h(Coarse());
  std::vector<double> samples;
  for (int k = 1; k < 40; ++k) {
    samples.push_back(5. * std::pow(0.02 / 5., k / 40.));
  }
  const double inv2pi2 = 1. / (2. * _PI_ * _PI_);
  std::vector<double> I_f;
  std::vector<double> I_p;
  std::vector<double> E_W;
  double f_mu_max  = 0.;
  double f_e_max   = 0.;
  double f_phi_max = 0.;
  std::vector<double> f_phi_end;
  const auto max_of = [](const double* f, size_t n, double& m) {
    for (size_t j = 0; j < n; ++j) {
      m = std::fmax(m, f[j]);
    }
  };
  h.Run(samples, [&](const C::Snapshot& s) {
    const double Nt = h.BareNumber(C::kTau, s.f_tau);
    const double Nm = h.BareNumber(C::kMu, s.f_mu);
    const double Ne = h.BareNumber(C::kE, s.f_e);
    const double Np = h.BareNumber(C::kPhi, s.f_phi);
    I_f.push_back(2. * Nt + 2. * Nm + 2. * Ne);
    I_p.push_back(Np + 2. * Nt - 2. * Ne);
    // E_c from the four species' own moments, independently of SectorEnergy.
    const double Ec = inv2pi2 * (2. * h.BareEnergy(C::kTau, s.f_tau, s.T_nu) +
                                 2. * h.BareEnergy(C::kMu, s.f_mu, s.T_nu) +
                                 2. * h.BareEnergy(C::kE, s.f_e, s.T_nu) +
                                 h.BareEnergy(C::kPhi, s.f_phi, s.T_nu));
    E_W.push_back(Ec - s.W);
    max_of(s.f_mu, h.q(C::kMu).size(), f_mu_max);
    max_of(s.f_e, h.q(C::kE).size(), f_e_max);
    max_of(s.f_phi, h.q(C::kPhi).size(), f_phi_max);
    f_phi_end.assign(s.f_phi, s.f_phi + h.q(C::kPhi).size());  // the last call is T_end
  });
  const auto& n0    = h.diagnostics().N_ini;
  const double I_f0 = 2. * n0[C::kTau] + 2. * n0[C::kMu] + 2. * n0[C::kE];
  const double I_p0 = n0[C::kPhi] + 2. * n0[C::kTau] - 2. * n0[C::kE];
  const double Ec0  = h.network().diagnostics().Ec_ini;
  double worst_f    = 0.;
  double worst_p    = 0.;
  double worst_E    = 0.;
  for (size_t k = 0; k < I_f.size(); ++k) {
    worst_f = std::fmax(worst_f, std::fabs(I_f[k] - I_f0) / I_f0);
    worst_p = std::fmax(worst_p, std::fabs(I_p[k] - I_p0) / I_f0);
    worst_E = std::fmax(worst_E, std::fabs(E_W[k] - Ec0) / Ec0);
  }
  assert(I_f.size() == samples.size() + 1);
  const auto& d     = h.diagnostics();
  const double Nphi = d.N_end[C::kPhi];
  const double pred = 4. * d.N_ini[C::kTau] + 2. * d.N_ini[C::kMu];
  std::printf(
      "invariants along the run: fermion %.2e, phi leg %.2e, energy %.2e (end: %.2e / %.2e / "
      "%.2e); N_phi,end = %.6f vs 4 N_tau,0 + 2 N_mu,0 = %.6f; N_tau,end %.1e, N_mu,end %.1e; "
      "max f_mu %.4f, f_e %.4f, f_phi %.3g; %lld steps, %lld RHS\n",
      worst_f,
      worst_p,
      worst_E,
      d.fermion_number_defect,
      d.phi_number_defect,
      d.energy_defect,
      Nphi,
      pred,
      d.N_end[C::kTau],
      d.N_end[C::kMu],
      f_mu_max,
      f_e_max,
      f_phi_max,
      d.stats.steps_accepted,
      d.stats.derivs_evaluations);
  assert(worst_f < 1e-10);
  assert(worst_p < 1e-10);
  assert(worst_E < 1e-5);
  assert(std::fabs(d.fermion_number_defect) < 1e-10);
  assert(std::fabs(d.phi_number_defect) < 1e-10);
  assert(std::fabs(d.energy_defect) < 1e-5);
  // Both parents gone: what is left of them is below the 1e-6 the identity is held to.
  assert(d.N_end[C::kTau] < 1e-8 * pred);
  assert(d.N_end[C::kMu] < 1e-6 * pred);
  assert(std::fabs(Nphi / pred - 1.) < 1e-6);
  assert(f_mu_max <= 1. && f_e_max <= 1.);
  assert(f_phi_max > 1e3);  // the cascade lases: phi is far above any thermal occupation

  // THE LASING PEAK ITSELF. The invariants above hold for any band factor, including a
  // broken one, so pin the physics too: the cold peak of q^2 f_phi (figs 5/6, q ~ 0.01).
  // On this grid it sits at q = 0.0113, height 17.6, i.e. 15% below the converged
  // 20.83 of an independent solver of the paper's equations that shares nothing with the
  // kernel (issue #389). The kernel's own refinement approaches it: 20.57 at its default
  // grids, 20.82 +- 0.04 extrapolated. The paper's 31.3 is reproduced by neither.
  // Bounds: a local maximum must exist in [0.004, 0.03] (no lasing, no peak) and its
  // height must stay within 5% of 17.6.
  const auto& qp = h.q(C::kPhi);
  double peak_q  = 0.;
  double peak_y  = 0.;
  for (size_t j = 1; j + 1 < qp.size(); ++j) {
    const auto y = [&](size_t i) { return qp[i] * qp[i] * f_phi_end[i]; };
    if (qp[j] < 0.004 || qp[j] > 0.03 || !(y(j) >= y(j - 1) && y(j) > y(j + 1))) {
      continue;
    }
    // Parabola in ln q through the three nodes.
    const double x0 = std::log(qp[j - 1]);
    const double x1 = std::log(qp[j]);
    const double x2 = std::log(qp[j + 1]);
    const double l0 = std::log(y(j - 1));
    const double l1 = std::log(y(j));
    const double l2 = std::log(y(j + 1));
    const double c2 = ((l2 - l1) / (x2 - x1) - (l1 - l0) / (x1 - x0)) / (x2 - x0);
    const double c1 = (l1 - l0) / (x1 - x0) - c2 * (x1 + x0);
    const double xv = -c1 / (2. * c2);
    peak_q          = std::exp(xv);
    peak_y          = std::exp(l1 + c1 * (xv - x1) + c2 * (xv * xv - x1 * x1));
  }
  std::printf("lasing peak at T_end: q = %.4f, q^2 f_phi = %.3f (pinned 17.6 +- 5%%)\n",
              peak_q,
              peak_y);
  assert(peak_y > 0.);
  assert(std::fabs(peak_y / 17.6 - 1.) < 0.05);
}

// THE CASCADE WITH K2 SWITCHED OFF IS THE SINGLE CHANNEL. With tau_mu = 1e12 s and
// nu_e empty, the cascade is DecayHarness's nu_tau -> nu_mu + phi with no spectator.
// Built on the harness's own grids (its shared daughter grid for both nu_mu and phi),
// the physics is then identical, so the two must agree to the integrator's tolerance:
// a slip in how the network sums channels, or a channel writing to the wrong species,
// could not. Bins below f = 1e-6 are left out: the phi bins just under the emission
// threshold (f ~ 1e-8) are set by absolute, not relative, tolerance.
void test_k2_off_is_single_channel() {
  DecayHarness::Config hc;
  hc.m_H_mev          = 2.;
  hc.tau_s            = 0.2025;
  hc.r_l              = 0.25;
  hc.r_phi            = 0.0025;
  hc.T_start_mev      = 5.;
  hc.T_end_mev        = 0.02;
  hc.parent_n         = 16;
  hc.daughter_n       = 48;
  hc.n_eff_spectator  = 0.;
  hc.daughter_thermal = true;
  DecayHarness h(hc);
  std::vector<double> fl;
  std::vector<double> fp;
  h.Run({}, [&](const DecayHarness::Snapshot& s) {
    fl.assign(s.f_l, s.f_l + h.daughter_q().size());
    fp.assign(s.f_phi, s.f_phi + h.daughter_q().size());
  });

  C::Config c = Coarse();
  c.tau_tau_s = hc.tau_s;
  c.tau_mu_s  = 1e12;
  c.e_thermal = false;
  c.kernel    = hc.kernel;  // the same scheme, gather floor included
  c.tau_n     = hc.parent_n;
  c.mu_n      = hc.daughter_n;
  c.phi_n     = hc.daughter_n;
  c.mu_q_min  = hc.daughter_q_min;
  c.phi_q_min = hc.daughter_q_min;
  c.mu_q_max  = h.daughter_q().back();
  c.phi_q_max = h.daughter_q().back();
  c.T_end_mev = hc.T_end_mev;
  C k(c);
  std::vector<double> km;
  std::vector<double> kp;
  k.Run({}, [&](const C::Snapshot& s) {
    km.assign(s.f_mu, s.f_mu + k.q(C::kMu).size());
    kp.assign(s.f_phi, s.f_phi + k.q(C::kPhi).size());
  });
  assert(k.q(C::kMu) == h.daughter_q());
  double worst = 0.;
  for (size_t j = 0; j < fl.size(); ++j) {
    if (fl[j] > 1e-6) {
      worst = std::fmax(worst, std::fabs(km[j] / fl[j] - 1.));
    }
    if (fp[j] > 1e-6) {
      worst = std::fmax(worst, std::fabs(kp[j] / fp[j] - 1.));
    }
  }
  std::printf("K2 off vs DecayHarness: worst |df/f| = %.2e over nu_mu and phi\n", worst);
  assert(worst < 1e-4);
}

// SAMPLING BY THE REAL T_gamma. Figs 6 and 7 are labelled by T_gamma, which the plasma
// alone sets (its entropy is conserved), so TnuAtTgamma inverts it in closed form. The
// ODE integrates T_gamma independently; the two must agree at the sample to the
// integrator's tolerance.
void test_tgamma_sampling() {
  C::Config c = Coarse();
  c.T_end_mev = 0.3;
  C h(c);
  std::vector<double> tg = {2.3, 1.1, 0.5};
  std::vector<double> tn;
  for (double t : tg) {
    tn.push_back(h.TnuAtTgamma(t));
  }
  std::vector<double> got;
  h.Run(tn, [&](const C::Snapshot& s) { got.push_back(s.T_gamma); });
  // Samples come back in time order (descending T) and the last call is T_end.
  assert(got.size() == 4);
  double worst = 0.;
  for (size_t i = 0; i < tg.size(); ++i) {
    worst = std::fmax(worst, std::fabs(got[i] / tg[i] - 1.));
  }
  std::printf("T_gamma sampling: T_nu(2.3, 1.1, 0.5) = %.5f, %.5f, %.5f MeV, worst %.1e\n",
              tn[0],
              tn[1],
              tn[2],
              worst);
  assert(worst < 1e-5);
  assert(tn[2] < 0.5 && tn[2] > 0.45);  // annihilation has begun heating the photons
}

// THE WEIGHTS THE KERNEL ASSUMES. Its rates are per internal state with the boson leg
// carrying g_H/g_phi = 2, so a network whose species weights differ would pass the
// kernel's own N/E identities and still break E_c - W. The network must refuse it.
void test_network_rejects_inconsistent_weights() {
  const auto make = [](double gH, double gl, double gphi) {
    DecayNetwork::Setup s;
    s.grids.push_back(DecayNetwork::LogGrid(1e-2, 20., 12));
    s.grids.push_back(DecayNetwork::LogGrid(1e-3, 40., 16, true));
    s.species = {{1., gH, 0}, {0., gl, 1}, {0., gphi, 1}};
    s.channels.emplace_back();
    s.T_start_mev = 5.;
    s.T_end_mev   = 1.;
    return s;
  };
  const auto refused = [&](double gH, double gl, double gphi) {
    try {
      DecayNetwork net(make(gH, gl, gphi));
    }
    catch (const std::invalid_argument&) {
      return true;
    }
    return false;
  };
  assert(!refused(2., 2., 1.));  // control: Hannestad's and the cascade's weights
  assert(!refused(4., 4., 2.));  // an overall scale is consistent
  assert(refused(2., 2., 2.));   // the reviewer's g_phi = 2
  assert(refused(2., 1., 1.));   // fermion leg
  assert(refused(1., 2., 1.));   // parent leg
  std::printf("network weights: (2,2,1) and (4,4,2) accepted, (2,2,2), (2,1,1), (1,2,1) refused\n");
}

}  // namespace

int main() {
  std::setvbuf(stdout, nullptr, _IONBF, 0);  // the numbers must survive a failed assert
  test_network_rejects_inconsistent_weights();
  test_tgamma_sampling();
  test_k2_off_is_single_channel();
  test_invariants_along_run();
  return 0;
}
