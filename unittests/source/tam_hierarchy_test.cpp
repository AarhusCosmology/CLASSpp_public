// The total angular momentum (TAM) photon hierarchy against the optimal one, on converged runs
// (long photon hierarchies, tensor tight coupling ended early), so that what is left is the
// difference between the two bases and nothing else.
//
// - Flat space: the two bases solve the same equations, so they agree to the truncation of
//   the long hierarchies (measured 2026-11-27: at most 7e-7 in TT, EE and BB).
// - Curved space: the optimal basis is an approximation for polarisation (Pitrou, Pereira &
//   Lesgourgues 2005.12119). At Omega_k = 0.1 tensor EE and BB at l = 2-10 move by 4.94% and
//   5.12%, scalar EE by 0.0635%; class_public v3.4.3 (Cyril Pitrou's implementation) gives
//   4.94%, 5.11% and 0.0633% at the same settings. TT stays within 2e-4. A mistake in the
//   TAM equations moves these numbers.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "cosmology.h"
#include "spectra_module.h"

namespace {

using Cls = std::map<std::string, std::vector<double>>;

Cls Spectra(const char* modes, const char* hierarchy, const char* Omega_k, int l_max) {
  FileContent fc;
  /* the cosmology of the cross-code check: the low-l numbers depend on it (scalar EE at
     l = 2-10 is 0.096% with the default tau_reio) */
  fc.set("h", "0.67810");
  fc.set("omega_b", "0.02238280");
  fc.set("omega_cdm", "0.1201075");
  fc.set("tau_reio", "0.05430842");
  fc.set("output", "tCl,pCl");
  fc.set("modes", modes);
  fc.set("hierarchy", hierarchy);
  if (Omega_k) {
    fc.set("Omega_k", Omega_k);
  }
  if (std::string(modes) == "s") {
    fc.set("A_s", "2.100549e-09");
    fc.set("n_s", "0.9660499");
    fc.set("l_max_scalars", std::to_string(l_max));
    fc.set("l_max_g", "50");
    fc.set("l_max_pol_g", "50");
  }
  else {
    fc.set("r", "0.1");
    fc.set("l_max_tensors", std::to_string(l_max));
    fc.set("l_max_g_ten", "40");
    fc.set("l_max_pol_g_ten", "40");
    fc.set("tight_coupling_trigger_tau_c_over_tau_k_ten", "0.002");
    fc.set("tight_coupling_trigger_tau_c_over_tau_h_ten", "0.002");
  }
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  Cosmology cosmology{fc};
  return cosmology.GetSpectraModule()->cl_output(l_max);
}

// max over l in [l_lo, l_hi] of |a/b - 1| for one spectrum
double MaxRel(const Cls& a, const Cls& b, const std::string& key, int l_lo, int l_hi) {
  double dev = 0.;
  for (int l = l_lo; l <= l_hi; ++l) {
    dev = std::max(dev, std::abs(a.at(key)[l] / b.at(key)[l] - 1.));
  }
  return dev;
}

void test_flat_space_bases_agree() {
  for (const char* modes : {"s", "t"}) {
    const int l_max = (std::string(modes) == "s") ? 1000 : 600;
    const Cls opt   = Spectra(modes, "optimal", nullptr, l_max);
    const Cls tam   = Spectra(modes, "tam", nullptr, l_max);
    for (const char* key : {"tt", "ee", "bb"}) {
      if (std::string(modes) == "s" && std::string(key) == "bb") {
        continue;
      }
      const double dev = MaxRel(tam, opt, key, 2, l_max);
      std::printf("flat, modes = %s, %s: max |tam/optimal - 1| = %.2e\n", modes, key, dev);
      assert(dev < 1e-5);
    }
  }
}

void test_curved_space_polarisation_moves() {
  const Cls st_opt  = Spectra("s", "optimal", "0.1", 1000);
  const Cls st_tam  = Spectra("s", "tam", "0.1", 1000);
  const double s_ee = MaxRel(st_tam, st_opt, "ee", 2, 10);
  const double s_tt = MaxRel(st_tam, st_opt, "tt", 2, 100);
  std::printf("Omega_k = 0.1, scalars: EE l = 2-10 %.4f%%, TT %.2e\n", 100. * s_ee, s_tt);
  assert(std::abs(s_ee - 6.35e-4) < 0.3e-4);
  assert(s_tt < 2e-4);

  const Cls t_opt   = Spectra("t", "optimal", "0.1", 600);
  const Cls t_tam   = Spectra("t", "tam", "0.1", 600);
  const double t_ee = MaxRel(t_tam, t_opt, "ee", 2, 10);
  const double t_bb = MaxRel(t_tam, t_opt, "bb", 2, 10);
  const double t_tt = MaxRel(t_tam, t_opt, "tt", 2, 100);
  std::printf("Omega_k = 0.1, tensors: EE l = 2-10 %.3f%%, BB %.3f%%, TT %.2e\n",
              100. * t_ee,
              100. * t_bb,
              t_tt);
  assert(std::abs(t_ee - 0.0494) < 0.001);
  assert(std::abs(t_bb - 0.0512) < 0.001);
  assert(t_tt < 2e-4);
}

}  // namespace

int main() {
  test_flat_space_bases_agree();
  test_curved_space_polarisation_moves();
  std::printf("TAM hierarchy test passed\n");
  return 0;
}
