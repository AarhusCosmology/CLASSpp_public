// The non-adiabatic initial conditions (cdi, nid, niv) are series in om*tau, om = a rho_m / sqrt(rho_r), in which
// om*tau is the leading order of the mode. With only the adiabatic start-time conditions the small-k modes start
// at om*tau ~ 1.5 (start_small_k_at_tau_c_over_tau_h = 0.006), where the series is wrong at the percent level,
// so start_isocurvature_at_omega_tau starts them no later than om*tau = 0.01.
//
// With the defaults, each mode must then agree with a start four times earlier. Measured 2026-11-05 on
// class_public with the same initial conditions: without the start condition cdi/bi EE at l = 2-5 move by 94%
// and nid by 5%; with it but with the niv shear, l=3 and eta values used until 2026-11, niv still moves by 1e-3
// at low l and 2e-2 at l > 30, because the wrong leading-order shear excites the free-streaming hierarchy.
// With ncdm (m = 0.06 eV), niv C_l^TT at l = 2 changed by a factor 240 between om*tau = 0.01 and 0.0005 at the
// default tol_ncdm_synchronous, since the ncdm momentum did not cancel the photon one (see NCDMSpecies).
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

constexpr int kLmax = 1500;

std::map<std::string, std::vector<double>> UnlensedCl(const char* ic,
                                                      const char* start_omega_tau,
                                                      bool ncdm) {
  FileContent fc;
  fc.set("output", "tCl,pCl");
  fc.set("ic", ic);
  if (ncdm) {
    fc.set("N_ncdm", "1");
    fc.set("m_ncdm", "0.06");
    fc.set("N_ur", "2.0328");
  }
  fc.set("l_max_scalars", std::to_string(kLmax));
  if (start_omega_tau) {
    fc.set("start_isocurvature_at_omega_tau", start_omega_tau);
  }
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  Cosmology cosmology{fc};
  return cosmology.GetSpectraModule()->cl_output(kLmax);
}

// max over l of |TT/TT_ref - 1|, |EE/EE_ref - 1| and |TE - TE_ref| / sqrt(TT_ref EE_ref)
double MaxDeviation(const char* ic, bool ncdm = false) {
  const auto def   = UnlensedCl(ic, nullptr, ncdm);
  const auto early = UnlensedCl(ic, "0.0025", ncdm);
  double dev       = 0.;
  int l_at         = 0;
  for (int l = 2; l <= kLmax; ++l) {
    const double tt = early.at("tt")[l];
    const double ee = early.at("ee")[l];
    const double d = std::max({std::abs(def.at("tt")[l] / tt - 1.),
                               std::abs(def.at("ee")[l] / ee - 1.),
                               std::abs(def.at("te")[l] - early.at("te")[l]) / std::sqrt(tt * ee)});
    if (d > dev) {
      dev  = d;
      l_at = l;
    }
  }
  std::printf("%s%s: default vs start at om*tau = 0.0025: max deviation %.2e at l = %d\n",
              ic,
              ncdm ? " + ncdm" : "",
              dev,
              l_at);
  return dev;
}

}  // namespace

int main() {
  assert(MaxDeviation("cdi") < 3e-4);
  assert(MaxDeviation("nid") < 3e-4);
  assert(MaxDeviation("niv") < 1e-3);
  // niv needs the photon and neutrino momenta to cancel at order tau^0: with ncdm, its momentum comes from the
  // perturbation quadrature, so the initial Psi_l must match the background density that sets the photon velocity.
  assert(MaxDeviation("niv", true) < 1e-3);
  std::printf("isocurvature start test passed\n");
  return 0;
}
