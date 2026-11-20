// Cosmology::GetDerivedParameter, the one lookup behind classy's get_current_derived_parameters
// (docs/superpowers/specs/2026-11-20-classy-package-design.md, section 3).
//
// Modules answer their own names. Species answer theirs through GetParam, so a name only a
// species knows needs no other code. A module's name keeps its meaning when a species answers
// it too (NEDE's z_star). Two species answering one name is an error that names them. Legacy
// species names keep classy's placeholders in runs without the species.
#include <cassert>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>

#include "background_module.h"
#include "cosmology.h"
#include "nonlinear_module.h"
#include "parser.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "thermodynamics_module.h"

namespace {

FileContent Inputs(const std::map<std::string, std::string>& extra) {
  FileContent fc;
  fc.set("h", "0.6781");
  fc.set("omega_b", "0.02238");
  fc.set("omega_cdm", "0.1201");
  fc.set("output", "");  // no perturbations to integrate
  for (const auto& [key, value] : extra) {
    fc.set(key, value);
  }
  return fc;
}

// Equal to an ULP or two. Where the library and the test both evaluate an expression,
// -ffast-math may order it differently in the two, and GCC and clang differ in that.
bool Close(double a, double b, double rel = 1e-15) {
  return std::fabs(a - b) <= rel * std::fabs(b);
}

// The message of the severe error the lookup throws, or "" if it answered.
std::string SevereMessage(Cosmology& cosmo, const std::string& name) {
  try {
    cosmo.GetDerivedParameter(name);
  }
  catch (const std::invalid_argument& error) {
    return error.what();
  }
  return "";
}

bool Contains(const std::string& text, const std::string& part) {
  return text.find(part) != std::string::npos;
}

void ModuleNames() {
  FileContent fc = Inputs({});
  Cosmology cosmo(fc);
  const double h = cosmo.GetInputModule()->background_.h;
  assert(cosmo.GetDerivedParameter("h") == h);
  // omega_m is Omega_m h^2 (classy returned Omega_m / h^2 until 2026-11-20).
  assert(Close(cosmo.GetDerivedParameter("omega_m"), cosmo.GetDerivedParameter("Omega_m") * h * h));
  assert(cosmo.GetDerivedParameter("rs_d_h") == cosmo.GetDerivedParameter("rs_d") * h);

  const auto& thermo = *cosmo.GetThermodynamicsModule();
  assert(Close(cosmo.GetDerivedParameter("100*theta_s"), 100. * thermo.rs_rec_ / thermo.ra_rec_));
  assert(cosmo.GetDerivedParameter("z_star") == thermo.z_star_);
  assert(Close(cosmo.GetDerivedParameter("exp_m_2_tau_As"),
               std::exp(-2 * thermo.tau_reionization_) * cosmo.GetPrimordialModule()->A_s_));

  // Refusals, with their reasons.
  assert(Contains(SevereMessage(cosmo, "no_such_parameter"),
                  "was not recognized as a derived parameter"));
  assert(Contains(SevereMessage(cosmo, "D_H"), "BBN network"));
  assert(Contains(SevereMessage(cosmo, "100*theta_d"), "damping scale"));
  assert(Contains(SevereMessage(cosmo, "nosuchspecies.alpha"), "no species 'nosuchspecies'"));
  assert(Contains(SevereMessage(cosmo, "CDM.alpha"), "does not report 'alpha'"));
  assert(cosmo.GetDerivedParameter("CDM.f_peak") > 0.);  // from the table, for any species

  // Legacy species names: present species answer, absent ones give classy's placeholders.
  const auto& background = *cosmo.GetBackgroundModule();
  assert(cosmo.GetDerivedParameter("Omega_Lambda") == background.GetOmega0Species("Lambda"));
  assert(cosmo.GetDerivedParameter("Omega0_lambda") == background.GetOmega0Species("Lambda"));
  assert(cosmo.GetDerivedParameter("Omega0_cdm") == background.GetOmega0Species("CDM"));
  assert(cosmo.GetDerivedParameter("Omega0_fld") == 0.);
  assert(cosmo.GetDerivedParameter("N_dg") == 0.);
  assert(cosmo.GetDerivedParameter("z_dec_drmd") == -1.);
  std::printf("module names, refusals and placeholders: ok\n");
}

void DampingScaleWhenComputed() {
  FileContent fc = Inputs({{"compute damping scale", "yes"}});
  Cosmology cosmo(fc);
  const auto& thermo = *cosmo.GetThermodynamicsModule();
  assert(Close(cosmo.GetDerivedParameter("100*theta_d"), 100. * thermo.rd_rec_ / thermo.ra_rec_));
  std::printf("100*theta_d with the damping scale: ok\n");
}

// a_c is in no table: the pheno-axion fluid answers it, and that is all it takes.
void SpeciesNameNeedsNoOtherCode() {
  FileContent fc = Inputs({{"fluid_equation_of_state", "pheno_axion"},
                           {"n_pheno_axion", "3"},
                           {"log10_axion_ac", "-3.5"},
                           {"fraction_fld_ac", "0.1"},
                           {"Theta_initial_fld", "2.8"}});
  Cosmology cosmo(fc);
  // To an ULP: GCC's pow and the species' own evaluation of 10^-3.5 need not agree exactly.
  assert(Close(cosmo.GetDerivedParameter("a_c"), std::pow(10., -3.5)));
  assert(cosmo.GetDerivedParameter("a_c") == cosmo.GetDerivedParameter("Fluid.a_c"));
  assert(cosmo.GetDerivedParameter("Omega0_fld") ==
         cosmo.GetBackgroundModule()->GetOmega0Species("Fluid"));
  std::printf("a species name by itself: ok\n");
}

// NEDE answers z_star (its transition); the plain name stays the recombination's.
void ModuleNameKeepsItsMeaning() {
  FileContent fc = Inputs({{"f_NEDE", "0.1"}, {"m_NEDE", "316.227766"}});
  Cosmology cosmo(fc);
  const double nede_z_star = cosmo.GetDerivedParameter("NEDE.z_star");
  assert(cosmo.GetDerivedParameter("z_star") == cosmo.GetThermodynamicsModule()->z_star_);
  assert(nede_z_star != cosmo.GetDerivedParameter("z_star"));
  std::printf("z_star with NEDE: ok\n");
}

// IDR inside IDM_DR_IDR answers N_dg and xi_idr; the composite answers its coupling.
void CompositeAnswersForItsParts() {
  FileContent fc = Inputs({{"Omega_idm_dr", "0.12"}, {"xi_idr", "0.3"}, {"a_idm_dr", "1e-4"}});
  Cosmology cosmo(fc);
  const auto& input      = cosmo.GetInputModule()->background_;
  const double Omega_idr = cosmo.GetBackgroundModule()->GetOmega0Species("IDR");
  assert(Omega_idr > 0.);
  // To an ULP: -ffast-math may order the arithmetic differently here and in the species.
  assert(Close(cosmo.GetDerivedParameter("N_dg"),
               Omega_idr / input.Omega0_g * 8. / 7. * std::pow(11. / 4., 4. / 3.)));
  assert(Close(cosmo.GetDerivedParameter("xi_idr"), 0.3));
  assert(cosmo.GetDerivedParameter("a_dark") == 1e-4);
  assert(Close(cosmo.GetDerivedParameter("Gamma_0_nadm"),
               1e-4 * (4. / 3.) * (input.h * input.h * Omega_idr)));
  std::printf("IDR inside IDM_DR_IDR: ok\n");
}

// Two greybody instances both answer alpha: the plain name is ambiguous, the dotted one is not.
void TwoInstancesAnsweringOneName() {
  std::map<std::string, std::string> extra;
  for (const auto& [key, alpha] : {std::pair{"gb1", "1.0"}, std::pair{"gb2", "2.5"}}) {
    const std::string k               = key;
    extra[k + ".type"]                = "ncdm_greybody";
    extra[k + ".parameterization"]    = "direct";
    extra[k + ".alpha"]               = alpha;
    extra[k + ".x"]                   = "1.0";
    extra[k + ".q0"]                  = "1.0";
    extra[k + ".m"]                   = "0.06";
    extra[k + ".deg"]                 = "1.0";
    extra[k + ".momenta_bins"]        = "7";
    extra[k + ".quadrature_strategy"] = "2";
    extra[k + ".max_q"]               = "18.0";
  }
  FileContent fc = Inputs(extra);
  Cosmology cosmo(fc);
  const std::string message = SevereMessage(cosmo, "alpha");
  assert(Contains(message, "'gb1'") && Contains(message, "'gb2'") &&
         Contains(message, "'<species>.alpha'"));
  assert(cosmo.GetDerivedParameter("gb1.alpha") == 1.0);
  assert(cosmo.GetDerivedParameter("gb2.alpha") == 2.5);
  std::printf("two instances answering alpha: ok\n");
}

// The checks classy made before calling these, now made by the modules.
void ChecksMovedFromClassy() {
  FileContent fc = Inputs({});
  Cosmology cosmo(fc);
  bool refused = false;
  try {
    cosmo.GetNonlinearModule()->GetSigma(8., 0., cosmo.GetNonlinearModule()->index_pk_m_);
  }
  catch (const std::invalid_argument& error) {
    refused = Contains(error.what(), "add mPk");
  }
  assert(refused);
  refused = false;
  try {
    cosmo.GetPerturbationsModule()->GetScalarTransferK();
  }
  catch (const std::invalid_argument& error) {
    refused = Contains(error.what(), "No transfer functions were computed");
  }
  assert(refused);

  // Flat: (chi(z2) - chi(z1))/(1 + z2); zero when the source is not behind the observer.
  const auto& background = *cosmo.GetBackgroundModule();
  const int index_chi    = background.index_bg_conf_distance_;
  const double chi1      = background.GetBackgroundAtZ(1.)[index_chi];
  const double chi2      = background.GetBackgroundAtZ(2.)[index_chi];
  assert(Close(background.GetAngularDistanceFromTo(1., 2.), (chi2 - chi1) / 3.));
  assert(background.GetAngularDistanceFromTo(2., 1.) == 0.);
  std::printf("checks moved from classy: ok\n");
}

}  // namespace

int main() {
  ChecksMovedFromClassy();
  ModuleNames();
  DampingScaleWhenComputed();
  SpeciesNameNeedsNoOtherCode();
  ModuleNameKeepsItsMeaning();
  CompositeAnswersForItsParts();
  TwoInstancesAnsweringOneName();
  return 0;
}
