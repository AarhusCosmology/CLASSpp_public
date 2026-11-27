// Tensor tight coupling must not leave the photons out of the gravitational-wave source.
//
// While tight coupling is on, the photon anisotropic stress is 8/5 a^2 rho_g Theta_2^(2), with
// Theta_2^(2) = -4/3 h'/kappa' to leading order in 1/kappa' (C. Pitrou, class_public v3.4.0). If the source
// leaves it out, h feels the photons only once tight coupling ends, and h depends on when that happens
// (tight_coupling_trigger_tau_c_over_tau_k_ten).
//
// With the defaults, the spectra that follow h must agree with tight coupling ended five times earlier
// (trigger 0.01 -> 0.002): TT at every l up to 600, EE and BB above l = 100. Without the term they differ by
// 2e-3; with it by 4e-5 at most.
//
// EE and BB below l = 100 are left out on purpose: they depend on the trigger by 3e-3 to 5e-3 with or without
// the term, through the polarisation source of the modes still tightly coupled during recombination.
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

constexpr int kLmax         = 600;
constexpr int kLminPol      = 101;
constexpr double kTolerance = 2e-4;

struct Case {
  const char* name;
  const char* tensor_method;  // nullptr: the default
  const char* Omega_k;        // nullptr: flat
};

std::map<std::string, std::vector<double>> TensorCl(const Case& c, const char* trigger) {
  FileContent fc;
  fc.set("output", "tCl,pCl");
  fc.set("modes", "t");
  fc.set("l_max_tensors", std::to_string(kLmax));
  if (c.tensor_method) {
    fc.set("tensor method", c.tensor_method);
  }
  if (c.Omega_k) {
    fc.set("Omega_k", c.Omega_k);
  }
  if (trigger) {
    fc.set("tight_coupling_trigger_tau_c_over_tau_k_ten", trigger);
  }
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  Cosmology cosmology{fc};
  return cosmology.GetSpectraModule()->cl_output(kLmax);
}

// max of |TT/TT_tight - 1| over l >= 2, and of |EE/EE_tight - 1|, |BB/BB_tight - 1| over l >= kLminPol
double MaxDeviation(const Case& c) {
  const auto def   = TensorCl(c, nullptr);
  const auto tight = TensorCl(c, "0.002");
  double dev       = 0.;
  int l_at         = 0;
  for (int l = 2; l <= kLmax; ++l) {
    double d = std::abs(def.at("tt")[l] / tight.at("tt")[l] - 1.);
    if (l >= kLminPol) {
      d = std::max({d,
                    std::abs(def.at("ee")[l] / tight.at("ee")[l] - 1.),
                    std::abs(def.at("bb")[l] / tight.at("bb")[l] - 1.)});
    }
    if (d > dev) {
      dev  = d;
      l_at = l;
    }
  }
  std::printf("%s: default vs tight trigger: max deviation %.2e at l = %d\n", c.name, dev, l_at);
  return dev;
}

}  // namespace

int main() {
  const Case cases[] = {{"ur", nullptr, nullptr},
                        {"photons only", "photons", nullptr},
                        {"open", nullptr, "0.02"},
                        {"closed", nullptr, "-0.02"}};
  for (const Case& c : cases) {
    assert(MaxDeviation(c) < kTolerance);
  }
  std::printf("tensor tight coupling test passed\n");
  return 0;
}
