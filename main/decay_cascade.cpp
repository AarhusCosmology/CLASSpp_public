/** @file decay_cascade.cpp Command-line driver for DecayCascade (tools/decay_harness.h):
 *  Hannestad & Madsen's three-level cascade nu_tau -> nu_mu + phi, nu_mu -> nu_e + phi
 *  through the MeV era, outside CLASS (issue #389, step 3).
 *
 *  Arguments are key=value pairs, every one optional:
 *
 *      decay_cascade [m_tau=2] [m_mu=0.5] [m_e=0.005] [m_phi=0.005]
 *                    [tau_tau=0.2] [tau_mu=1] [T_start=5] [T_end=0.01]
 *                    [samples=2.3,1] [samples_Tgamma=2.3,1.1,0.5] [spectra=1]
 *                    [n_eff_spectator=0] [inverse_decays=1] [quantum_statistics=1]
 *                    [balanced_gather=1] [gather_floor=1e-10] [emission_gauss=2]
 *                    [tau_n=48] [mu_n=112] [e_n=128] [phi_n=160]
 *                    [tau_q_min=1e-2] [mu_q_min=1e-3] [e_q_min=1e-4] [phi_q_min=1e-5]
 *                    [tau_q_max=30] [mu_q_max=0] [e_q_max=0] [phi_q_max=0] [rtol=1e-6]
 *
 *  Masses and temperatures in MeV, lifetimes in s. `samples` are T_nu values and
 *  `samples_Tgamma` real photon temperatures (converted through the plasma's entropy,
 *  DecayNetwork::TnuAtTgamma) at which to print the state in addition to T_end.
 *
 *  Output on stdout: one block per sample, blank-line separated. Each block starts
 *  with `#` lines: T_nu, T_gamma, t, the bare numbers N_s = Sum q^2 dq f per dof, and
 *  n_s / n_gamma with the real T_gamma. Unless spectra=0, one row per bin follows,
 *
 *      <species> q f q^2 f
 *
 *  with q = p/T_nu and species one of tau, mu, e, phi. The conservation diagnostics and
 *  the integrator's counters follow on stderr.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "common.h"
#include "decay_harness.h"

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

constexpr double kZeta3 = 1.2020569031595942;

}  // namespace

int main(int argc, char** argv) {
  try {
    std::map<std::string, std::string> kv;
    for (int i = 1; i < argc; ++i) {
      const char* eq = std::strchr(argv[i], '=');
      if (!eq) {
        std::fprintf(stderr, "decay_cascade: argument '%s' is not key=value\n", argv[i]);
        return 2;
      }
      kv[std::string(argv[i], static_cast<size_t>(eq - argv[i]))] = std::string(eq + 1);
    }

    DecayCascade::Config c;
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
    num("m_tau", c.m_tau_mev);
    num("m_mu", c.m_mu_mev);
    num("m_e", c.m_e_mev);
    num("m_phi", c.m_phi_mev);
    num("tau_tau", c.tau_tau_s);
    num("tau_mu", c.tau_mu_s);
    num("T_start", c.T_start_mev);
    num("T_end", c.T_end_mev);
    num("n_eff_spectator", c.n_eff_spectator);
    num("tau_q_min", c.tau_q_min);
    num("tau_q_max", c.tau_q_max);
    num("mu_q_min", c.mu_q_min);
    num("mu_q_max", c.mu_q_max);
    num("e_q_min", c.e_q_min);
    num("e_q_max", c.e_q_max);
    num("phi_q_min", c.phi_q_min);
    num("phi_q_max", c.phi_q_max);
    num("rtol", c.rtol);
    num("gather_floor", c.kernel.gather_floor);
    integer("tau_n", c.tau_n);
    integer("mu_n", c.mu_n);
    integer("e_n", c.e_n);
    integer("phi_n", c.phi_n);
    integer("emission_gauss", c.kernel.emission_gauss);
    flag("mu_thermal", c.mu_thermal);
    flag("e_thermal", c.e_thermal);
    flag("inverse_decays", c.kernel.inverse_decays);
    flag("quantum_statistics", c.kernel.quantum_statistics);
    flag("balanced_gather", c.kernel.balanced_gather);
    flag("lumped_loss", c.kernel.lumped_loss);
    bool spectra = true;
    flag("spectra", spectra);
    std::vector<double> samples;
    std::vector<double> samples_tg;
    if (kv.count("samples") != 0) {
      samples = ParseList(kv["samples"]);
      kv.erase("samples");
    }
    if (kv.count("samples_Tgamma") != 0) {
      samples_tg = ParseList(kv["samples_Tgamma"]);
      kv.erase("samples_Tgamma");
    }
    // An unread key is a typo, and a typo silently runs the default.
    if (!kv.empty()) {
      std::fprintf(stderr, "decay_cascade: unknown argument '%s'\n", kv.begin()->first.c_str());
      return 2;
    }

    DecayCascade h(c);
    for (double Tg : samples_tg) {
      samples.push_back(h.TnuAtTgamma(Tg));
    }
    const char* names[DecayCascade::kNumSpecies] = {"tau", "mu", "e", "phi"};
    const double g[DecayCascade::kNumSpecies]    = {2., 2., 2., 1.};
    h.Run(samples, [&](const DecayCascade::Snapshot& s) {
      const double* f[DecayCascade::kNumSpecies] = {s.f_tau, s.f_mu, s.f_e, s.f_phi};
      std::printf("# T_nu = %.10g MeV\n# T_gamma = %.10g MeV\n# t = %.10g s\n# W = %.10g\n",
                  s.T_nu,
                  s.T_gamma,
                  s.t,
                  s.W);
      // n_s = g_s/(2 pi^2) T_nu^3 N_s and n_gamma = 2 zeta(3)/pi^2 T_gamma^3.
      const double r3 = (s.T_nu / s.T_gamma) * (s.T_nu / s.T_gamma) * (s.T_nu / s.T_gamma);
      for (int k = 0; k < DecayCascade::kNumSpecies; ++k) {
        const double N = h.BareNumber(k, f[k]);
        std::printf("# N_%s = %.12g\n# E_%s = %.12g\n# n_over_ngamma_%s = %.10g\n",
                    names[k],
                    N,
                    names[k],
                    h.BareEnergy(k, f[k], s.T_nu),
                    names[k],
                    g[k] * N * r3 / (4. * kZeta3));
      }
      if (spectra) {
        for (int k = 0; k < DecayCascade::kNumSpecies; ++k) {
          const auto& q = h.q(k);
          for (size_t j = 0; j < q.size(); ++j) {
            std::printf("%s %.10e %.10e %.10e\n", names[k], q[j], f[k][j], q[j] * q[j] * f[k][j]);
          }
        }
      }
      std::printf("\n\n");
    });
    const auto& d = h.diagnostics();
    std::fprintf(stderr,
                 "decay_cascade: fermion number drift %.3e, phi-leg drift %.3e, energy defect "
                 "%.3e; N_phi,end = %.6f (4 N_tau,0 + 2 N_mu,0 = %.6f); %lld steps (%lld "
                 "rejected), %lld right-hand sides, %lld Jacobians, %lld LU\n",
                 d.fermion_number_defect,
                 d.phi_number_defect,
                 d.energy_defect,
                 d.N_end[DecayCascade::kPhi],
                 4. * d.N_ini[DecayCascade::kTau] + 2. * d.N_ini[DecayCascade::kMu],
                 d.stats.steps_accepted,
                 d.stats.steps_rejected,
                 d.stats.derivs_evaluations,
                 d.stats.jacobians,
                 d.stats.lu_decompositions);
  }
  catch (const std::exception& e) {
    std::fprintf(stderr, "decay_cascade: %s\n", e.what());
    return 1;
  }
  return 0;
}
