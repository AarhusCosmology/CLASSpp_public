// In the full-Limber C_L^phiphi, multipole L at wavenumber k comes from comoving distance
// (L + 1/2)/k: a q grid that stops at the last perturbation mode leaves out everything closer
// than that, where halofit's non-linear power is largest. The lensing source is therefore
// continued past the last mode (AsymptoticExtrapolation, times halofit's correction there) up
// to k_max_extrapolated, instead of solving modes that far
// (docs/superpowers/specs/2026-10-21-lensing-source-extrapolation-design.md).
//
// With the defaults, C_L^phiphi must match a run that solves modes to 70/Mpc. Without the
// continuation it was 2% low at high L. The same continuation lets halofit trust the linear
// spectrum beyond the last mode, so nonlinear_min_k_max came down from 5 to 2/Mpc: P_NL must
// also still match the solved modes, up to z = 3.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "cosmology.h"
#include "nonlinear_module.h"
#include "spectra_module.h"

namespace {

constexpr int kLmax = 2500;

std::vector<double> ClPhiPhi(bool solve_to_70) {
  FileContent fc;
  fc.set("output", "tCl,pCl,lCl");
  fc.set("lensing", "yes");
  fc.set("non linear", "halofit");
  fc.set("l_max_scalars", std::to_string(kLmax));
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  if (solve_to_70) {
    fc.set("nonlinear_min_k_max", "70");
  }
  Cosmology cosmology{fc};
  return cosmology.GetSpectraModule()->cl_output(kLmax).at("pp");
}

/* max |P_NL / P_NL(modes solved to 70/Mpc) - 1| for k < 0.67/Mpc and z <= 3, P(k) only */
double PkNonlinearError() {
  auto pk = [](bool solve_to_70) {
    FileContent fc;
    fc.set("output", "mPk");
    fc.set("non linear", "halofit");
    fc.set("P_k_max_h/Mpc", "1");
    fc.set("z_max_pk", "3");
    fc.set("input_verbose", "0");
    fc.set("write warnings", "no");
    fc.set("write parameters", "no");
    if (solve_to_70) {
      fc.set("nonlinear_min_k_max", "70");
      fc.set("k_max_extrapolated", "0");
    }
    Cosmology cosmology{fc};
    const auto& nl = *cosmology.GetNonlinearModule();
    std::vector<double> out;
    for (double z : {0., 1., 2., 3.}) {
      for (double k = 1e-3; k < 0.67; k *= 1.2) {
        double p = 0.;
        nl.nonlinear_pk_at_k_and_z(pk_nonlinear, k, z, nl.index_pk_m_, &p, nullptr);
        out.push_back(p);
      }
    }
    return out;
  };
  const auto solved = pk(true);
  const auto p      = pk(false);
  double worst      = 0.;
  for (size_t i = 0; i < p.size(); ++i) {
    worst = std::max(worst, std::abs(p[i] / solved[i] - 1.));
  }
  return worst;
}

}  // namespace

int main() {
  const double pk_error = PkNonlinearError();
  std::printf("halofit P_NL (z <= 3) against modes solved to 70/Mpc: worst %.2e\n", pk_error);
  assert(pk_error < 1e-3);

  const auto solved = ClPhiPhi(true);
  const auto cl     = ClPhiPhi(false);
  double worst      = 0.;
  int l_worst       = 0;
  for (int l = 2; l <= kLmax; ++l) {
    const double d = std::abs(cl[l] / solved[l] - 1.);
    if (d > worst) {
      worst   = d;
      l_worst = l;
    }
  }
  std::printf("C_L^phiphi against modes solved to 70/Mpc: worst %.2e at L = %d\n", worst, l_worst);
  assert(worst < 3e-3);
  std::printf("lensing extrapolation test passed\n");
  return 0;
}
