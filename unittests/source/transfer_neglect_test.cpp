// The CMB transfer functions are set to zero where l < (q - delta_k(l)) r_a(tau_rec): beyond
// that q the source sits on the oscillating tail of the Bessel function. delta_k grows with l,
// delta_k(l) = transfer_neglect_delta_k_S_X sqrt(l / (l + transfer_neglect_l_half_S_X)), because
// at low l the integrand dies off within a small fraction of the q a constant cut kept
// (docs/superpowers/specs/2026-10-21-transfer-neglect-per-l-design.md).
//
// With the defaults, what the cut neglects must stay small against the same run with nothing
// neglected. The constant cuts used until 2026-10-21 neglected 2.8e-4 of C_l^TT at the median
// over l <= 2500 and 4.7e-4 at worst; the per-l defaults measure 3.5e-5 and 1.2e-4.
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

constexpr int kLmax = 2500;

std::map<std::string, std::vector<double>> UnlensedCl(bool neglect_nothing) {
  FileContent fc;
  fc.set("output", "tCl,pCl");
  fc.set("l_max_scalars", std::to_string(kLmax));
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  if (neglect_nothing) {
    for (const char* type : {"t0", "t1", "t2", "e"}) {
      fc.set(std::string("transfer_neglect_delta_k_S_") + type, "100");
      fc.set(std::string("transfer_neglect_l_half_S_") + type, "0");
    }
  }
  Cosmology cosmology{fc};
  return cosmology.GetSpectraModule()->cl_output(kLmax);
}

double Median(std::vector<double> v) {
  std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
  return v[v.size() / 2];
}

}  // namespace

int main() {
  const auto cut  = UnlensedCl(false);
  const auto full = UnlensedCl(true);

  std::vector<double> tt;
  double tt_max = 0.;
  double ee_max = 0.;
  double te_max = 0.;
  for (int l = 2; l <= kLmax; ++l) {
    const double tt_full = full.at("tt")[l];
    const double ee_full = full.at("ee")[l];
    tt.push_back(std::abs(cut.at("tt")[l] / tt_full - 1.));
    tt_max = std::max(tt_max, tt.back());
    ee_max = std::max(ee_max, std::abs(cut.at("ee")[l] / ee_full - 1.));
    te_max = std::max(te_max,
                      std::abs(cut.at("te")[l] - full.at("te")[l]) / std::sqrt(tt_full * ee_full));
  }
  const double tt_median = Median(tt);
  std::printf("neglected: TT median %.2e, max %.2e; EE max %.2e; TE max %.2e of sqrt(TT EE)\n",
              tt_median,
              tt_max,
              ee_max,
              te_max);
  assert(tt_median < 1e-4);
  assert(tt_max < 2.5e-4);
  assert(ee_max < 1e-4);
  assert(te_max < 2e-4);
  std::printf("transfer neglect test passed\n");
  return 0;
}
