// The non-linear corrections are computed for every time of the perturbation sources, and the
// times are shared among threads. Each time is computed as it was in the serial loop, so the
// result must not depend on the number of threads, bit for bit, for halofit and HMcode alike.
//
// The serial loop stopped at the latest time at which the correction is not computable (high z,
// where k_NL exceeds the trusted k_max) and left R_NL = 1 at all earlier times. Threads compute
// some of those earlier times before the stop is known, and they are reset afterwards. The times
// run to z = 10, beyond the stop (z ~ 7 with the continued spectrum, z ~ 2 without it).
#include <cassert>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "cosmology.h"
#include "nonlinear_module.h"

namespace {

struct Corrections {
  std::vector<std::vector<double>> nl_corr_density;
  std::vector<double> k_nl;
  int k_size = 0;
};

Corrections Run(const std::vector<std::pair<std::string, std::string>>& settings, int threads) {
  FileContent fc;
  fc.set("output", "mPk");
  fc.set("z_max_pk", "10");
  fc.set("N_ncdm", "1");
  fc.set("m_ncdm", "0.06");
  fc.set("N_ur", "2.0308");
  fc.set("threads", std::to_string(threads));
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  for (const auto& [key, value] : settings) {
    fc.set(key, value);
  }
  Cosmology cosmology{fc};
  const auto& nl = *cosmology.GetNonlinearModule();
  Corrections out{nl.nl_corr_density_, {}, nl.k_size_};
  for (double z : {0., 1., 2., 3.}) {
    double k_nl    = 0.;
    double k_nl_cb = 0.;
    nl.nonlinear_k_nl_at_z(z, &k_nl, &k_nl_cb);
    out.k_nl.push_back(k_nl);
    out.k_nl.push_back(k_nl_cb);
  }
  return out;
}

/* is R_NL = 1 for every k at time slice index_tau of P(k) type index_pk? */
bool IsLinear(const Corrections& c, int index_pk, int index_tau) {
  for (int index_k = 0; index_k < c.k_size; ++index_k) {
    if (c.nl_corr_density[index_pk][index_tau * c.k_size + index_k] != 1.) {
      return false;
    }
  }
  return true;
}

}  // namespace

int main() {
  const std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string>>>>
      cases = {{"halofit", {{"non linear", "halofit"}}},
               {"HMcode 2016", {{"non linear", "hmcode"}, {"hmcode_version", "2016"}}},
               {"HMcode 2020", {{"non linear", "hmcode"}, {"hmcode_version", "2020"}}},
               {"halofit, k_max_extrapolated = 0",
                {{"non linear", "halofit"}, {"k_max_extrapolated", "0"}}}};

  for (const auto& [name, settings] : cases) {
    const Corrections serial   = Run(settings, 1);
    const Corrections parallel = Run(settings, 8);
    const int tau_size         = static_cast<int>(serial.nl_corr_density[0].size()) / serial.k_size;
    const bool same            = (serial.nl_corr_density == parallel.nl_corr_density) &&
                                 (serial.k_nl == parallel.k_nl);
    std::printf("%-32s 1 vs 8 threads: %s; R_NL = 1 at the first time: %s, at the last: %s\n",
                name.c_str(),
                same ? "identical" : "DIFFERENT",
                IsLinear(serial, 0, 0) ? "yes" : "no",
                IsLinear(serial, 0, tau_size - 1) ? "yes" : "no");
    assert(same);
    /* the correction is computed today, and is not computable at the earliest time, z = 10 */
    assert(!IsLinear(serial, 0, tau_size - 1));
    assert(IsLinear(serial, 0, 0));
  }
  std::printf("nonlinear threads test passed\n");
  return 0;
}
