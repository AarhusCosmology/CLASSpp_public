// Unit tests for DecayHarness (tools/decay_harness): one DecayTransitionKernel channel
// driven through the MeV era outside CLASS (issue #389, step 2). Each test pins one
// piece of the harness against something it cannot share code with.
#include "decay_harness.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "common.h"  // _PI_

namespace {

// THE PLASMA, with the sector switched off in effect (a lifetime of 1e12 s). The plasma's
// entropy is conserved through annihilation, so with a T_nu = 1 the final ratio is fixed
// by the entropy at the START:
//
//     (T_gamma/T_nu)^3 = s(T_start) / (T_start^3 * 4 pi^2/45),
//
// which is (11/4)^{1/3} only in the limit T_start >> m_e; at 10 MeV the pairs carry
// slightly less entropy than that and the ratio is ~6e-5 lower. Pinned against both: the
// exact one tightly, the textbook one at the size of that correction.
void test_standard_model_limit() {
  DecayHarness::Config c;
  c.m_H_mev     = 1e-6;  // relativistic throughout, so it is a third neutrino
  c.tau_s       = 1e12;
  c.T_start_mev = 10.;
  c.T_end_mev   = 1e-3;
  c.parent_n    = 24;
  c.daughter_n  = 24;
  c.rtol        = 1e-8;
  DecayHarness h(c);
  double ratio = 0.;
  double tH    = 0.;
  h.Run({}, [&](const DecayHarness::Snapshot& s) {
    ratio = s.T_gamma / s.T_nu;
    tH    = s.t * s.H;
  });
  const BbnPlasma plasma;
  const BbnEmState em  = plasma.ElectromagneticState(c.T_start_mev);
  const double s_start = (em.rho + em.pressure) / c.T_start_mev;
  const double exact   = std::cbrt(s_start / (std::pow(c.T_start_mev, 3) * 4. * _PI_ * _PI_ / 45.));
  const double ur      = std::cbrt(11. / 4.);
  std::printf(
      "standard-model limit: T_gamma/T_nu = %.8f (entropy at T_start: %.8f, "
      "(11/4)^(1/3) = %.8f), t*H = %.6f\n",
      ratio,
      exact,
      ur,
      tH);
  assert(std::fabs(ratio / exact - 1.) < 1e-6);
  assert(std::fabs(ratio / ur - 1.) < 1e-4);
  // Radiation domination: t H -> 1/2 once the start-up transient has redshifted away.
  assert(std::fabs(tH - 0.5) < 1e-3);
}

// DECAY ONLY, against an independent quadrature. Without inverse decays or quantum
// statistics each parent bin decays at its own time-dilated rate, so
//
//     ln(f_H(q, t) / f_H(q, 0)) = -Gamma int_0^t (a m / eps(q, a)) dt'.
//
// The integral is formed here by the trapezoid rule over the harness's own samples of
// t(a), which is a different discretisation from the ODE's. This is the test that the
// kernel's conformal-time output, the a = 1/T_nu units and the 1/(aH) conversion to
// ln a are mutually consistent: an error in any of them is an O(1) factor in the rate.
void test_decay_only_rate() {
  DecayHarness::Config c;
  c.m_H_mev                   = 5.;
  c.tau_s                     = 1.;
  c.r_l                       = 0.1;
  c.T_start_mev               = 5.;
  c.T_end_mev                 = 0.3;
  c.parent_n                  = 32;
  c.daughter_n                = 64;
  c.rtol                      = 1e-9;
  c.kernel.inverse_decays     = false;
  c.kernel.quantum_statistics = false;
  DecayHarness h(c);

  std::vector<double> samples;
  const int ns = 1200;
  for (int k = 0; k < ns; ++k) {
    samples.push_back(c.T_start_mev *
                      std::pow(c.T_end_mev / c.T_start_mev, static_cast<double>(k) / (ns - 1)));
  }
  const auto& q = h.parent_q();
  const int nP  = static_cast<int>(q.size());
  std::vector<double> integral(nP, 0.);
  std::vector<double> f0(nP);
  std::vector<double> f_end(nP);
  double t_prev = -1.;
  std::vector<double> rate_prev(nP, 0.);
  h.Run(samples, [&](const DecayHarness::Snapshot& s) {
    const double am = c.m_H_mev / s.T_nu;
    for (int i = 0; i < nP; ++i) {
      const double r = am / std::sqrt(q[i] * q[i] + am * am);
      if (t_prev >= 0.) {
        integral[i] += 0.5 * (r + rate_prev[i]) * (s.t - t_prev);
      }
      else {
        f0[i] = s.f_H[i];
      }
      rate_prev[i] = r;
      f_end[i]     = s.f_H[i];
    }
    t_prev = s.t;
  });
  double worst = 0.;
  for (int i = 0; i < nP; ++i) {
    const double expect = -integral[i] / c.tau_s;
    worst = std::fmax(worst, std::fabs(std::log(f_end[i] / f0[i]) - expect) / std::fabs(expect));
  }
  std::printf(
      "decay only: ln f_H vs -Gamma int (m/E) dt, worst %.2e over %d bins (t_end = %.3g s)\n",
      worst,
      nP,
      t_prev);
  assert(worst < 1e-4);
  assert(t_prev > 2. * c.tau_s);  // the parent really did decay
}

// THE FULL NETWORK, inverse decays and quantum statistics on, massive fermion daughter.
// The kernel's identities are exact per right-hand side, so over a whole run the two
// number invariants and the comoving energy budget may drift only by the integrator's
// error, which is what is asserted. A daughter grid sized by hand too small is the
// control: it must show up in the energy budget, because clamping files energy at the
// grid edge.
void test_full_network_budgets() {
  DecayHarness::Config c;
  c.m_H_mev     = 5.;
  c.tau_s       = 1.;
  c.r_l         = 0.1;
  c.T_start_mev = 5.;
  c.T_end_mev   = 0.05;
  c.parent_n    = 32;
  c.daughter_n  = 64;
  c.rtol        = 1e-8;
  DecayHarness h(c);
  double fl_peak = 0.;
  h.Run({}, [&](const DecayHarness::Snapshot& s) {
    for (size_t j = 0; j < h.daughter_q().size(); ++j) {
      fl_peak = std::fmax(fl_peak, s.f_l[j]);
    }
  });
  const auto& d = h.diagnostics();
  std::printf("full network (r_l = 0.1): number drift %.2e / %.2e, energy defect %.2e\n",
              d.number_defect_l,
              d.number_defect_phi,
              d.energy_defect);
  assert(std::fabs(d.number_defect_l) < 1e-6);
  assert(std::fabs(d.number_defect_phi) < 1e-6);
  assert(std::fabs(d.energy_defect) < 1e-6);
  assert(fl_peak <= 1.);  // Pauli blocking holds the fermion daughter at or below 1

  DecayHarness::Config small = c;
  small.daughter_q_max       = 3.;  // far below the a*m/2 the decays emit at
  DecayHarness hs(small);
  hs.Run({}, [](const DecayHarness::Snapshot&) {});
  std::printf("   control, daughter grid cut at q = 3: energy defect %.2e\n",
              hs.diagnostics().energy_defect);
  assert(std::fabs(hs.diagnostics().energy_defect) > 1e-3);
}

// THE BOTTOM DAUGHTER CELL IS THE RESERVOIR FOR q < q_min. A massive daughter is emitted
// at rest by any parent moving at the rest-frame daughter velocity p*/E* (0.98 at
// m_F/m_H = 0.1), and at T ~ 9 MeV those parents sit in the bulk of the thermal
// distribution. Everything below the grid is filed into bin 0, so bin 0 has to own that
// volume. With the trapezoid's half cell it held ~6x too little, was driven to the Pauli
// limit (f = 1 - 6e-7), and the integrator stalled: this run did not finish in ten
// minutes. Asserted on both symptoms: the bin reads like its neighbour, and it is cheap.
void test_bottom_cell_owns_rest() {
  DecayHarness::Config c;
  c.m_H_mev     = 5.;
  c.tau_s       = 0.2;
  c.r_l         = 0.1;
  c.T_start_mev = 10.;
  c.T_end_mev   = 5.;
  DecayHarness h(c);
  double f0 = 0.;
  double f1 = 0.;
  h.Run({}, [&](const DecayHarness::Snapshot& s) {
    f0 = s.f_l[0];
    f1 = s.f_l[1];
  });
  const auto& st = h.diagnostics().stats;
  std::printf("bottom cell: f_l[0] = %.4f, f_l[1] = %.4f, %lld steps\n", f0, f1, st.steps_accepted);
  assert(std::fabs(f0 / f1 - 1.) < 0.02);
  assert(st.steps_accepted < 5000);
}

}  // namespace

int main() {
  test_bottom_cell_owns_rest();
  test_standard_model_limit();
  test_decay_only_rate();
  test_full_network_budgets();
  return 0;
}
