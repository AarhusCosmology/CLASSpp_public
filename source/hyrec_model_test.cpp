#include "hyrec_model.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <string>

#include "parser.h"
#include "precision.h"

namespace {

/* Planck-2018-like, matching what a default CLASS run hands the model. */
constexpr double kT_cmb    = 2.7255;
constexpr double kObh2     = 0.02238;
constexpr double kOcbh2    = 0.14240;
constexpr double kOkh2     = 0.;
constexpr double kNeff     = 3.046;
constexpr double kYHe      = 0.2454;
constexpr double kFHe      = 0.0813; /* helium by number, relative to hydrogen */
constexpr double kNH0      = 0.1904; /* hydrogen nuclei today, m^-3 */
constexpr double kZInitial = 20000.;

/* The vendored table directory, as ResolveDataPaths would hand it over. */
std::string DefaultPath() {
  return std::string(__CLASSDIR__) + "/hyrec/";
}

HyrecModel MakeModel(const std::string& path = DefaultPath()) {
  return HyrecModel(path, kT_cmb, kObh2, kOcbh2, kOkh2, kNeff, kYHe, kFHe, kNH0, kZInitial);
}

/* A state deep in hydrogen recombination, where SWIFT's fitting function is the
   branch that runs. Densities and rates are the ones CLASS would supply. */
RecombinationState StateAt(double z, double x_H, double x_He) {
  RecombinationState s = {};
  s.z                  = z;
  s.x_H                = x_H;
  s.x_He               = x_He;
  s.x                  = x_H + kFHe * x_He;
  s.n_H                = kNH0 * (1. + z) * (1. + z) * (1. + z);
  s.T_rad              = kT_cmb * (1. + z);
  s.T_mat              = s.T_rad;
  /* H(z) in s^-1, matter+radiation, near enough for a unit test. */
  s.H = 2.2e-18 * std::sqrt(0.31 * std::pow(1. + z, 3) + 9.2e-5 * std::pow(1. + z, 4));
  s.hydrogen_frozen = false;
  s.helium_ode      = true;
  return s;
}

const EnergyDeposition kNoInjection = {0., 0., 0., 0., 0.};

/* The vendored tables load, and the model reports which code and mode it is. */
void test_tables_load() {
  HyrecModel model = MakeModel();
  assert(std::string(model.Name()).find("HYREC-2") != std::string::npos);
  assert(std::string(model.Name()).find("SWIFT") != std::string::npos);
}

/* HYREC-2 builds its file names by strcat-ing onto hyrec_path, so a directory
   given without the trailing separator would look for "<dir>Alpha_inf.dat".
   precision::parse appends it; check on the real directory, so the assertion is
   that the tables load, not merely that a string ends in a slash. */
void test_hyrec_path_without_trailing_separator_still_loads() {
  FileContent fc;
  fc.set("hyrec_path", std::string(__CLASSDIR__) + "/hyrec");

  precision ppr;
  ppr.parse(fc);
  assert(ppr.hyrec_path == DefaultPath());

  HyrecModel model = MakeModel(ppr.hyrec_path);
  assert(std::string(model.Name()).find("HYREC-2") != std::string::npos);
}

/* A missing data directory must surface as a CLASS exception, never as the
   exit(1) that HYREC-2 would otherwise reach on its own. */
void test_missing_tables_throw_rather_than_exit() {
  bool threw = false;
  try {
    HyrecModel bad = MakeModel("/nonexistent-hyrec-path/");
  }
  catch (const std::exception&) {
    threw = true;
  }
  assert(threw);
}

/* The derivative must be finite everywhere SWIFT runs, and must drive x_H back
   towards equilibrium: at z = 1600 an under-ionized gas re-ionizes (dx_H/dz < 0,
   so x_H rises as z falls), while by z = 900 the same x_H is far above
   equilibrium and recombines (dx_H/dz > 0). Pinning the flip pins the sign
   convention -- both the dlna -> dz change of variable and its direction. */
void test_hydrogen_relaxes_towards_equilibrium() {
  HyrecModel model = MakeModel();

  for (double z = 1600.; z >= 700.; z -= 25.) {
    IonisationDerivatives d = model.Derivatives(StateAt(z, 0.5, 1e-8), kNoInjection, 0.);
    assert(std::isfinite(d.dx_H_dz));
  }

  assert(model.Derivatives(StateAt(1600., 0.5, 1e-8), kNoInjection, 0.).dx_H_dz < 0.);
  assert(model.Derivatives(StateAt(900., 0.5, 1e-8), kNoInjection, 0.).dx_H_dz > 0.);

  /* Once recombination is under way, the further x_H sits above equilibrium the
     faster it falls. */
  const double at_09 = model.Derivatives(StateAt(1100., 0.9, 1e-8), kNoInjection, 0.).dx_H_dz;
  const double at_05 = model.Derivatives(StateAt(1100., 0.5, 1e-8), kNoInjection, 0.).dx_H_dz;
  const double at_01 = model.Derivatives(StateAt(1100., 0.1, 1e-8), kNoInjection, 0.).dx_H_dz;
  assert(at_09 > at_05 && at_05 > at_01 && at_01 > 0.);
}

/* Frozen hydrogen is the driver's call, not the model's: CLASS holds x_H on its
   Saha branch in the early phases and the model must report exactly zero. */
void test_frozen_hydrogen_reports_zero() {
  HyrecModel model     = MakeModel();
  RecombinationState s = StateAt(1400., 0.9, 0.5);
  s.hydrogen_frozen    = true;
  assert(model.Derivatives(s, kNoInjection, 0.).dx_H_dz == 0.);
}

/* Helium: evolved only while CLASS is integrating it and while HYREC-2's own
   window applies. Outside either, the derivative is exactly zero rather than a
   division by a vanishing neutral-helium fraction. */
void test_helium_window() {
  HyrecModel model = MakeModel();

  /* Inside the window: HeII recombining, with hydrogen still on its Saha branch
     -- which is exactly what CLASS's helium phase hands the model. */
  RecombinationState inside = StateAt(2200., 0.999999, 0.5);
  inside.hydrogen_frozen    = true;
  IonisationDerivatives d   = model.Derivatives(inside, kNoInjection, 0.);
  assert(std::isfinite(d.dx_He_dz));
  assert(d.dx_He_dz != 0.);

  /* CLASS is holding helium at Saha (its analytic phase, where x_He can reach 2
     and HYREC-2's HeII -> HeI routine does not apply). */
  RecombinationState saha = inside;
  saha.helium_ode         = false;
  assert(model.Derivatives(saha, kNoInjection, 0.).dx_He_dz == 0.);

  /* Fully ionized helium: no neutral helium for the escape probability. */
  RecombinationState ionized = StateAt(2200., 0.999999, 1.);
  ionized.hydrogen_frozen    = true;
  assert(model.Derivatives(ionized, kNoInjection, 0.).dx_He_dz == 0.);

  /* Helium finished: below HYREC-2's XHEII_MIN it stops tracking. */
  RecombinationState done = StateAt(1400., 0.9, 1e-9);
  assert(model.Derivatives(done, kNoInjection, 0.).dx_He_dz == 0.);
}

/* Outside HYREC-2's effective-rate tables the multilevel atom has no data, and
   calling it anyway is a hard error rather than a graceful fallback. Both ends
   must fall back to Peebles instead: kB*T_rad drops under TR_MIN = 0.004 eV at
   low z, and rises over TR_MAX = 0.4 eV above z ~ 1700. */
void test_peebles_fallback_outside_the_tables() {
  HyrecModel model = MakeModel();

  for (double z : {12., 5., 1.}) {
    IonisationDerivatives d = model.Derivatives(StateAt(z, 1e-4, 1e-9), kNoInjection, 0.);
    assert(std::isfinite(d.dx_H_dz));
  }

  for (double z : {1800., 2200., 3000., 6000.}) {
    IonisationDerivatives d = model.Derivatives(StateAt(z, 0.999, 1e-9), kNoInjection, 0.);
    assert(std::isfinite(d.dx_H_dz));
  }
}

/* Energy injection reaches HYREC-2 as a deposition rate and must speed
   ionization, i.e. reduce the net recombination rate. */
void test_injection_slows_recombination() {
  HyrecModel model     = MakeModel();
  RecombinationState s = StateAt(1100., 0.3, 1e-8);

  const double without = model.Derivatives(s, kNoInjection, 0.).dx_H_dz;

  const EnergyDeposition dep = {0.15, 0.35, 0.02, 0.32, 0.14};
  const double with          = model.Derivatives(s, dep, 1e-24).dx_H_dz;

  assert(std::isfinite(with));
  assert(with < without);
}

/* A units slip -- cm^-3 for m^-3, eV for K -- would not change the sign, so check
   the magnitude too: dx_H/dz is O(1e-2) in the thick of recombination. */
void test_derivative_magnitude_is_sane() {
  HyrecModel hyrec     = MakeModel();
  RecombinationState s = StateAt(1100., 0.3, 1e-8);
  const double d_hyrec = hyrec.Derivatives(s, kNoInjection, 0.).dx_H_dz;
  assert(d_hyrec > 0.);
  assert(d_hyrec < 1.); /* dx_H/dz is O(1e-3) here; a unit slip would blow this */
  assert(d_hyrec > 1e-8);
}

}  // namespace

int main() {
  test_tables_load();
  test_hyrec_path_without_trailing_separator_still_loads();
  test_missing_tables_throw_rather_than_exit();
  test_hydrogen_relaxes_towards_equilibrium();
  test_frozen_hydrogen_reports_zero();
  test_helium_window();
  test_peebles_fallback_outside_the_tables();
  test_injection_slows_recombination();
  test_derivative_magnitude_is_sane();

  std::printf("hyrec model tests passed\n");
  return 0;
}
