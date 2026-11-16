// Tensor initial conditions to order (k tau)^2: h = h_ini [1 - (k^2+2K) tau^2 / (6 + 8/5 rho_fs/rho_r)], with
// rho_fs the free-streaming species whose anisotropic stress enters the equation of h, and their monopoles
// F_0 = sqrt(6) (h - h_ini). A mode that starts with h' = 0 keeps an error of order (k tau_ini)^2, so that the
// tensor C_l's depend on the starting time.
//
// With the defaults, each tensor spectrum must agree with a start four times earlier in k tau
// (start_large_k_at_tau_h_over_tau_k) and at tau_c/tau_H = 4e-4 for the smallest wavenumbers
// (start_small_k_at_tau_c_over_tau_h). The cases cover the three sources of the stress (ur, ncdm as massless
// species, and the ncdm momentum hierarchy), no stress at all, and curvature. In the third case all neutrinos are
// massive, so that it fails without the initial monopoles of the momentum bins (1.6e-3, as ur without its own).
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

constexpr int kLmax = 300;

// With h' = 0 at the start the deviation is 4e-3 to 7e-3. With the (k tau)^2 term it is 2e-4 to 3e-4 with the
// stress of the neutrinos and 4e-4 without it, largest at l < 10: the smallest wavenumbers start close to
// equality, where the radiation-era series has a correction of relative order om*tau.
constexpr double kTolerance = 1e-3;

struct Case {
  const char* name;
  const char* tensor_method;  // nullptr: the default
  int n_massive;              // 0, 1 (0.06 eV) or 3 (degenerate, 0.02 eV each)
  const char* Omega_k;        // nullptr: flat
};

std::map<std::string, std::vector<double>> TensorCl(const Case& c, bool early) {
  FileContent fc;
  fc.set("output", "tCl,pCl");
  fc.set("modes", "t");
  fc.set("l_max_tensors", std::to_string(kLmax));
  if (c.tensor_method) {
    fc.set("tensor method", c.tensor_method);
  }
  if (c.n_massive == 1) {
    fc.set("N_ncdm", "1");
    fc.set("m_ncdm", "0.06");
    fc.set("N_ur", "2.0328");
  }
  if (c.n_massive == 3) {
    fc.set("N_ncdm", "1");
    fc.set("deg_ncdm", "3");
    fc.set("m_ncdm", "0.02");
    fc.set("N_ur", "0.00641");
  }
  if (c.Omega_k) {
    fc.set("Omega_k", c.Omega_k);
  }
  if (early) {
    fc.set("start_small_k_at_tau_c_over_tau_h", "0.0004");
    fc.set("start_large_k_at_tau_h_over_tau_k", "0.0175");
  }
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  Cosmology cosmology{fc};
  return cosmology.GetSpectraModule()->cl_output(kLmax);
}

// max over l of |TT/TT_ref - 1|, |EE/EE_ref - 1| and |BB/BB_ref - 1|
double MaxDeviation(const Case& c) {
  const auto def   = TensorCl(c, false);
  const auto early = TensorCl(c, true);
  double dev       = 0.;
  int l_at         = 0;
  for (int l = 2; l <= kLmax; ++l) {
    const double d = std::max({std::abs(def.at("tt")[l] / early.at("tt")[l] - 1.),
                               std::abs(def.at("ee")[l] / early.at("ee")[l] - 1.),
                               std::abs(def.at("bb")[l] / early.at("bb")[l] - 1.)});
    if (d > dev) {
      dev  = d;
      l_at = l;
    }
  }
  std::printf("%s: default vs early start: max deviation %.2e at l = %d\n", c.name, dev, l_at);
  return dev;
}

}  // namespace

int main() {
  const Case cases[] = {{"ur", nullptr, 0, nullptr},
                        {"ncdm as massless", nullptr, 1, nullptr},
                        {"ncdm hierarchy", "exact", 3, nullptr},
                        {"photons only", "photons", 0, nullptr},
                        {"open", nullptr, 0, "0.02"},
                        {"closed", nullptr, 0, "-0.02"}};
  for (const Case& c : cases) {
    assert(MaxDeviation(c) < kTolerance);
  }
  std::printf("tensor start test passed\n");
  return 0;
}
