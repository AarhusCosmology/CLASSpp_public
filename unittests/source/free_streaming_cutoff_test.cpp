// Halofit and HMcode integrate the linear P(k) continued past the last perturbation mode up to
// k_max_extrapolated, in the high-k form of cold matter. The free-streaming cutoff of warm dark
// matter is not of that form: with 1 keV thermal warm dark matter as a quarter of the dark
// matter, the default modes (to 2/Mpc) put P_NL at z = 3 5e-2 off. So a species whose cutoff
// lies there has its modes solved up to k_max_extrapolated, and P_NL is that of the solved modes
// (docs/superpowers/specs/2026-10-22-nonlinear-followups-design.md).
//
// Pinned: 1 keV triggers it, and P_NL then matches a run that solves to 70/Mpc without the
// continuation exactly; 10 keV warm dark matter (whose cutoff starts near the end of the
// continued range), 0.06 eV neutrinos and a run without a non-linear method do not.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "cosmology.h"
#include "nonlinear_module.h"
#include "perturbations_module.h"

namespace {

using Settings = std::vector<std::pair<std::string, std::string>>;

/* thermal warm dark matter of mass m_kev as a quarter of the dark matter */
Settings WarmDarkMatter(double m_kev) {
  const double omega_ncdm = 0.03;
  const double T_ncdm     = 0.71611 * std::cbrt(omega_ncdm * 93.14 / (m_kev * 1e3));
  return {{"N_ncdm", "1"},
          {"omega_ncdm", std::to_string(omega_ncdm)},
          {"T_ncdm", std::to_string(T_ncdm)},
          {"omega_cdm", "0.09"}};
}

struct Result {
  double k_last = 0.;
  std::vector<double> pk_nl; /* z = 3, k < 1/Mpc; empty without a non-linear method */
};

Result Run(const Settings& model, const Settings& extra) {
  FileContent fc;
  fc.set("output", "mPk");
  fc.set("P_k_max_1/Mpc", "1");
  fc.set("z_max_pk", "3");
  fc.set("non linear", "halofit");
  fc.set("input_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  for (const auto& [key, value] : model) {
    fc.set(key, value);
  }
  for (const auto& [key, value] : extra) {
    fc.set(key, value);
  }
  Cosmology cosmology{fc};
  const auto& pt = *cosmology.GetPerturbationsModule();
  Result out;
  out.k_last = pt.k_[pt.index_md_scalars_].back();
  if (!extra.empty() && (extra[0].first == "non linear") && extra[0].second.empty()) {
    return out;
  }
  const auto& nl = *cosmology.GetNonlinearModule();
  for (double k = 1e-3; k < 1.; k *= 1.2) {
    double p = 0.;
    nl.nonlinear_pk_at_k_and_z(pk_nonlinear, k, 3., nl.index_pk_m_, &p, nullptr);
    out.pk_nl.push_back(p);
  }
  return out;
}

}  // namespace

int main() {
  const Settings neutrinos = {{"N_ncdm", "1"}, {"m_ncdm", "0.06"}, {"N_ur", "2.0308"}};

  const Result wdm1     = Run(WarmDarkMatter(1.), {});
  const Result wdm1_ref = Run(WarmDarkMatter(1.),
                              {{"nonlinear_min_k_max", "70"}, {"k_max_extrapolated", "0"}});
  const Result wdm10    = Run(WarmDarkMatter(10.), {});
  const Result nu       = Run(neutrinos, {});
  const Result wdm1_lin = Run(WarmDarkMatter(1.), {{"non linear", ""}});

  std::printf("last mode: 1 keV %.3g, 10 keV %.3g, 0.06 eV neutrinos %.3g, 1 keV linear %.3g\n",
              wdm1.k_last,
              wdm10.k_last,
              nu.k_last,
              wdm1_lin.k_last);
  assert(wdm1.k_last >= 70.);
  assert(wdm10.k_last < 10.);
  assert(nu.k_last < 10.);
  assert(wdm1_lin.k_last < 10.);

  assert(wdm1.pk_nl.size() == wdm1_ref.pk_nl.size());
  double worst = 0.;
  for (size_t i = 0; i < wdm1.pk_nl.size(); ++i) {
    worst = std::max(worst, std::abs(wdm1.pk_nl[i] / wdm1_ref.pk_nl[i] - 1.));
  }
  std::printf("1 keV P_NL(z = 3) against modes solved to 70/Mpc: worst %.2e\n", worst);
  assert(worst < 1e-12);
  std::printf("free-streaming cutoff test passed\n");
  return 0;
}
