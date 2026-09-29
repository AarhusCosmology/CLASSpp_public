/* Shooting: scaled tolerances, adoption of the converged build, and the decaying-NCDM
   budget closure. docs/superpowers/specs/2026-09-29-shooting-reuse-and-dncdm-budget-design.md

   Every expectation below is a physical identity or an input, never a number computed by
   the code under test: 100*theta_s must equal the requested value, a flat universe must
   have H(z=0) = H0, and a combined-mode sector must carry the density it was asked for. */

#include <cmath>
#include <cstdio>
#include <exception>
#include <stdexcept>
#include <string>

#include "background_module.h"
#include "cosmology.h"
#include "input_module.h"
#include "thermodynamics_module.h"

namespace {

int failures = 0;

void expect(bool ok, const char* what) {
  if (!ok) {
    ++failures;
    std::printf("  FAILED: %s\n", what);
  }
}

template <class F>
void run(const char* name, F&& test) {
  std::printf("%s\n", name);
  try {
    test();
  }
  catch (const std::exception& e) {
    ++failures;
    std::printf("  FAILED: threw: %s\n", e.what());
  }
}

FileContent quiet() {
  FileContent fc;
  fc.set("input_verbose", "0");
  fc.set("background_verbose", "0");
  fc.set("thermodynamics_verbose", "0");
  fc.set("write warnings", "no");
  fc.set("write parameters", "no");
  return fc;
}

// A = [0,1,1] A1/NO-like sector: m = 0.05 eV parent, one thermal nu_l daughter, phi empty,
// next to a stable 0.0087 eV spectator; N_ur completes N_eff = 3.044.
FileContent decay_sector(const char* representation, const char* Gamma, bool one_grid = true) {
  FileContent fc = quiet();
  fc.set("H0", "67.96");
  fc.set("omega_b", "0.0223");
  fc.set("omega_cdm", "0.118");
  fc.set("N_ur", "0.0044");
  fc.set("nu2.type", "ncdm_standard");
  fc.set("nu2.m", "0.0087");
  fc.set("nu2.T", "0.71611");
  fc.set("nu2.deg", "1");
  fc.set("d.type", "ncdm_decay_dr");
  fc.set("d.m", "0.05");
  fc.set("d.T", "0.71611");
  fc.set("d.Gamma", Gamma);
  // One momentum grid for background and perturbations: every representation evolves the
  // parent's distribution on it (the integrated one is refused otherwise, see below).
  if (one_grid) {
    fc.set("d.quadrature_strategy", "3");
    fc.set("d.momenta_bins", "16");
  }
  if (std::string(representation) == "proxy") {
    fc.set("d.dr_representation", "proxy");
    fc.set("d.dr_f_ini_l", "1");
    fc.set("d.dr_f_ini_phi", "0");
    fc.set("evolver_background", "3");
  }
  return fc;
}

double H_today_over_H0(Cosmology& cosmology) {
  const BackgroundModulePtr& bgm = cosmology.GetBackgroundModule();
  const double* today = bgm->background_table_.data() + (bgm->bt_size_ - 1) * bgm->bg_size_;
  return today[bgm->index_bg_H_] / cosmology.GetInputModule()->background_.H0;
}

double sector_Omega_today(Cosmology& cosmology) {
  const BackgroundModulePtr& bgm = cosmology.GetBackgroundModule();
  const double* today      = bgm->background_table_.data() + (bgm->bt_size_ - 1) * bgm->bg_size_;
  const InputModulePtr& im = cosmology.GetInputModule();
  const double H0          = im->background_.H0;
  return im->all_species_.at("d")->Rho(today) / (H0 * H0);
}

// Before: DoShooting passed 1e-3 absolute to Newton, and a guess within 1e-3 of the target
// in 100*theta_s was accepted (3-4e-5 off in practice, ~0.1 sigma of Planck).
void theta_s_shoot_lands_on_target() {
  FileContent fc = quiet();
  fc.set("100*theta_s", "1.04200");
  fc.set("omega_b", "0.02237");
  fc.set("omega_cdm", "0.1200");
  Cosmology cosmology{fc};
  const ThermodynamicsModulePtr& thm = cosmology.GetThermodynamicsModule();
  const double theta                 = 100. * thm->rs_rec_ / thm->ra_rec_;
  std::printf("  100*theta_s = %.8f (asked 1.04200)\n", theta);
  expect(std::fabs(theta - 1.04200) < 2e-6, "std::fabs(theta - 1.04200) < 2e-6");
}

// Before: a deg-normalised sector reported Omega0 = 0 to the closure, so Omega_Lambda took its
// share and the universe over-closed: H(0)/H0 - 1 = 5.9e-4 at Gamma = 0, 1.5e-5 at high Gamma.
void deg_mode_closes_the_budget(const char* representation, const char* Gamma) {
  FileContent fc = decay_sector(representation, Gamma);
  Cosmology cosmology{fc};
  const double excess = H_today_over_H0(cosmology) - 1.;
  std::printf("  %-10s Gamma = %-6s H(0)/H0 - 1 = %+.2e\n", representation, Gamma, excess);
  expect(std::fabs(excess) < 1e-6, "std::fabs(excess) < 1e-6");
}

// Before: DoShooting returned a fresh InputModule and the Cosmology rebuilt the background the
// converged residual had already built. At Gamma = 0 the stable-sector seed is exact, so the
// first evaluation converges and must become the production build; the thermodynamics module
// is then built lazily on top of the adopted background.
void converged_evaluation_is_adopted() {
  FileContent fc = decay_sector("integrated", "0");
  Cosmology cosmology{fc};
  const InputModulePtr& im = cosmology.GetInputModule();
  std::printf("  evaluations = %d, converged build kept = %d\n",
              im->shooting_evaluations_,
              static_cast<int>(im->shooting_build_kept_));
  expect(im->shooting_build_kept_, "im->shooting_build_kept_");
  expect(std::fabs(H_today_over_H0(cosmology) - 1.) < 1e-6,
         "std::fabs(H_today_over_H0(cosmology) - 1.) < 1e-6");
  const ThermodynamicsModulePtr& thm = cosmology.GetThermodynamicsModule();
  expect(std::isfinite(thm->rs_rec_) && thm->rs_rec_ > 0.,
         "std::isfinite(thm->rs_rec_) && thm->rs_rec_ > 0.");
}

// Where the proxy's background is expensive (Gamma12 >= 1e-6 here, 2-57 s), the high-Gamma seed
// must be within tolerance: one background, and that one is the production build.
void high_gamma_proxy_closes_in_one_build() {
  FileContent fc = decay_sector("proxy", "1e6");
  Cosmology cosmology{fc};
  const InputModulePtr& im = cosmology.GetInputModule();
  const double excess      = H_today_over_H0(cosmology) - 1.;
  std::printf("  proxy Gamma = 1e6: evaluations = %d, kept = %d, H(0)/H0 - 1 = %+.2e\n",
              im->shooting_evaluations_,
              static_cast<int>(im->shooting_build_kept_),
              excess);
  expect(im->shooting_evaluations_ == 1, "im->shooting_evaluations_ == 1");
  expect(im->shooting_build_kept_, "im->shooting_build_kept_");
  expect(std::fabs(excess) < 1e-6, "std::fabs(excess) < 1e-6");
}

// Before: the proxy had no shooting hooks (Omega_dncdmdr was never shot: at Gamma = 0 its
// unshot guess gave N_eff = 97.8), and the integrated representation's guess below Gamma = H0
// was ~94x too large, sending Newton to N_eff = -5251.
void combined_mode_hits_its_target(const char* representation) {
  FileContent fc = decay_sector(representation, "0");
  fc.set("d.omega_dncdmdr", "0.000540");
  Cosmology cosmology{fc};
  const double h     = cosmology.GetInputModule()->background_.h;
  const double asked = 0.000540 / (h * h);
  const double got   = sector_Omega_today(cosmology);
  std::printf("  %-10s Omega_dncdmdr: asked %.8e, got %.8e, H(0)/H0 - 1 = %+.2e, %d evaluations\n",
              representation,
              asked,
              got,
              H_today_over_H0(cosmology) - 1.,
              cosmology.GetInputModule()->shooting_evaluations_);
  expect(std::fabs(got - asked) < 1e-6, "std::fabs(got - asked) < 1e-6");
  // At Gamma = 0 the stable-sector guess is exact, so one build. The old guess (~94x too
  // large) took three, probing unphysical deg (N_eff = -5249) on the way.
  expect(cosmology.GetInputModule()->shooting_evaluations_ == 1,
         "cosmology.GetInputModule()->shooting_evaluations_ == 1");
  expect(std::fabs(H_today_over_H0(cosmology) - 1.) < 1e-6,
         "std::fabs(H_today_over_H0(cosmology) - 1.) < 1e-6");
}

// Before: the proxy never shot the initial-abundance fixed point either.
void initial_mode_proxy_closes_the_budget() {
  FileContent fc = decay_sector("proxy", "0");
  fc.set("d.Neff_ini", "1.0132");
  Cosmology cosmology{fc};
  const double excess = H_today_over_H0(cosmology) - 1.;
  std::printf("  proxy Neff_ini: H(0)/H0 - 1 = %+.2e\n", excess);
  expect(std::fabs(excess) < 1e-6, "std::fabs(excess) < 1e-6");
}

// Before: with the default quadrature the integrated representation sampled its parent on 9
// background points but evolved ln f on the 3 perturbation points, so 6 bins never decayed and
// kept sourcing decay radiation: Omega_DR = 0.70 today at Gamma = 1e6 (the parent holds 1.2e-3).
// psd and proxy refuse the mismatch; the integrated representation must too, as a SEVERE
// (configuration) error.
void integrated_representation_refuses_two_grids() {
  FileContent fc = decay_sector("integrated", "1e6", /*one_grid=*/false);
  bool refused   = false;
  try {
    Cosmology cosmology{fc};
    cosmology.GetBackgroundModule();
  }
  catch (const std::invalid_argument& e) {
    refused = std::string(e.what()).find("q_size") != std::string::npos;
    std::printf("  refused: %s\n", e.what());
  }
  expect(refused, "refused");
}

}  // namespace

int main() {
  run("theta_s shoot lands on target", theta_s_shoot_lands_on_target);
  run("deg mode closes the budget: integrated, Gamma = 0",
      [] { deg_mode_closes_the_budget("integrated", "0"); });
  run("deg mode closes the budget: integrated, Gamma = 1e6",
      [] { deg_mode_closes_the_budget("integrated", "1e6"); });
  run("deg mode closes the budget: proxy, Gamma = 0",
      [] { deg_mode_closes_the_budget("proxy", "0"); });
  run("converged evaluation is adopted", converged_evaluation_is_adopted);
  run("high-Gamma proxy closes in one build", high_gamma_proxy_closes_in_one_build);
  run("combined mode hits its target: integrated",
      [] { combined_mode_hits_its_target("integrated"); });
  run("combined mode hits its target: proxy", [] { combined_mode_hits_its_target("proxy"); });
  run("initial mode, proxy", initial_mode_proxy_closes_the_budget);
  run("integrated representation refuses two grids", integrated_representation_refuses_two_grids);
  std::printf(failures ? "shooting_test: %d FAILED\n" : "shooting_test: all passed\n", failures);
  return failures ? 1 : 0;
}
