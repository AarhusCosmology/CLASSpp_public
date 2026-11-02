/** @file decay_harness.cpp Command-line driver for DecayHarness (tools/decay_harness.h):
 *  one decay channel H <-> l + phi through the MeV era, outside CLASS (issue #389).
 *
 *  Arguments are key=value pairs, every one optional except m_H and tau:
 *
 *      decay_harness m_H=5 tau=1 r_l=0.1 [r_phi=0] [T_start=10] [T_end=0.01]
 *                    [samples=2.3,1,0.5] [n_eff_spectator=1] [daughter_thermal=1]
 *                    [inverse_decays=1] [quantum_statistics=1] [balanced_gather=1]
 *                    [parent_n=48] [daughter_n=96] [daughter_q_min=1e-3] [rtol=1e-6]
 *
 *  Masses and temperatures in MeV, tau in s. `samples` are T_nu values at which to
 *  print the spectra in addition to T_end.
 *
 *  Output on stdout: one block per sample, blank-line separated. Each block starts
 *  with `#` lines giving T_nu, T_gamma, t, then one row per daughter bin:
 *
 *      p/T_gamma   f_l   f_phi   (p/T_gamma)^2 f_l   (p/T_gamma)^2 f_phi
 *
 *  then one row per parent bin, prefixed by `parent`. The run's conservation
 *  diagnostics follow on stderr, so a figure script can read stdout untouched.
 */

#include "decay_harness.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<double> ParseList(const std::string& s) {
  std::vector<double> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, ',')) {
    if (!item.empty()) {
      out.push_back(std::stod(item));
    }
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; ++i) {
      const char* eq = std::strchr(argv[i], '=');
      if (!eq) {
        std::fprintf(stderr, "decay_harness: argument '%s' is not key=value\n", argv[i]);
        return 2;
      }
      kv[std::string(argv[i], static_cast<size_t>(eq - argv[i]))] = std::string(eq + 1);
    }
    if (kv.count("m_H") == 0 || kv.count("tau") == 0) {
      std::fprintf(stderr, "decay_harness: m_H=<MeV> and tau=<s> are required\n");
      return 2;
    }

    DecayHarness::Config c;
    std::vector<double> samples;
    const auto num = [&](const char* key, double& dst) {
      if (kv.count(key) != 0) {
        dst = std::stod(kv[key]);
        kv.erase(key);
      }
    };
    const auto integer = [&](const char* key, int& dst) {
      if (kv.count(key) != 0) {
        dst = std::stoi(kv[key]);
        kv.erase(key);
      }
    };
    const auto flag = [&](const char* key, bool& dst) {
      if (kv.count(key) != 0) {
        dst = std::stoi(kv[key]) != 0;
        kv.erase(key);
      }
    };
    num("m_H", c.m_H_mev);
    num("tau", c.tau_s);
    num("r_l", c.r_l);
    num("r_phi", c.r_phi);
    num("T_start", c.T_start_mev);
    num("T_end", c.T_end_mev);
    num("n_eff_spectator", c.n_eff_spectator);
    num("parent_q_min", c.parent_q_min);
    num("parent_q_max", c.parent_q_max);
    num("daughter_q_min", c.daughter_q_min);
    num("daughter_q_max", c.daughter_q_max);
    num("rtol", c.rtol);
    num("gather_floor", c.kernel.gather_floor);
    integer("parent_n", c.parent_n);
    integer("daughter_n", c.daughter_n);
    integer("emission_gauss", c.kernel.emission_gauss);
    flag("daughter_thermal", c.daughter_thermal);
    flag("inverse_decays", c.kernel.inverse_decays);
    flag("quantum_statistics", c.kernel.quantum_statistics);
    flag("balanced_gather", c.kernel.balanced_gather);
    flag("lumped_loss", c.kernel.lumped_loss);
    if (kv.count("samples") != 0) {
      samples = ParseList(kv["samples"]);
      kv.erase("samples");
    }
    // An unread key is a typo, and a typo silently runs the default.
    if (!kv.empty()) {
      std::fprintf(stderr, "decay_harness: unknown argument '%s'\n", kv.begin()->first.c_str());
      return 2;
    }

    DecayHarness h(c);
    const auto& qd = h.daughter_q();
    const auto& qp = h.parent_q();
    h.Run(samples, [&](const DecayHarness::Snapshot& s) {
      std::printf("# T_nu = %.10g MeV\n# T_gamma = %.10g MeV\n# t = %.10g s\n",
                  s.T_nu,
                  s.T_gamma,
                  s.t);
      std::printf("# p/T_gamma  f_l  f_phi  (p/T_gamma)^2 f_l  (p/T_gamma)^2 f_phi\n");
      // q = p/T_nu, so p/T_gamma = q T_nu / T_gamma.
      const double conv = s.T_nu / s.T_gamma;
      for (size_t j = 0; j < qd.size(); ++j) {
        const double x = qd[j] * conv;
        std::printf("%.10e %.10e %.10e %.10e %.10e\n",
                    x,
                    s.f_l[j],
                    s.f_phi[j],
                    x * x * s.f_l[j],
                    x * x * s.f_phi[j]);
      }
      for (size_t i = 0; i < qp.size(); ++i) {
        const double x = qp[i] * conv;
        std::printf("parent %.10e %.10e %.10e\n", x, s.f_H[i], x * x * s.f_H[i]);
      }
      std::printf("\n\n");
    });
    const auto& d = h.diagnostics();
    std::fprintf(stderr,
                 "decay_harness: number drift (fermion leg) %.3e, (boson leg) %.3e, "
                 "energy defect %.3e; %lld steps (%lld rejected), %lld right-hand sides, "
                 "%lld Jacobians\n",
                 d.number_defect_l,
                 d.number_defect_phi,
                 d.energy_defect,
                 d.stats.steps_accepted,
                 d.stats.steps_rejected,
                 d.stats.derivs_evaluations,
                 d.stats.jacobians);
  }
  catch (const std::exception& e) {
    std::fprintf(stderr, "decay_harness: %s\n", e.what());
    return 1;
  }
  return 0;
}
