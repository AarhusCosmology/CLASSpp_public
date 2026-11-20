#include "cosmology.h"

#include <cmath>
#include <map>
#include <set>
#include <string>

#include "background_module.h"
#include "lensing_module.h"
#include "nonlinear_module.h"
#include "output_module.h"
#include "perturbations_module.h"
#include "primordial_module.h"
#include "spectra_module.h"
#include "thermodynamics_module.h"
#include "transfer_module.h"

InputModulePtr& Cosmology::GetInputModule() {
  if (!shot_) {
    ShootingResult shoot = InputModule::DoShooting(input_module_ptr_);
    if (shoot.converged) {
      // The solver's last evaluation was built from exactly the resolved file content, so it
      // IS this cosmology: adopt it whole, with every module the residual already built.
      // docs/superpowers/specs/2026-09-29-shooting-reuse-and-dncdm-budget-design.md (B)
      // All but the Bessel cache, which belongs to whoever constructed this cosmology.
      std::shared_ptr<BesselTableCache> bessel_cache = std::move(bessel_cache_);
      *this                                          = std::move(*shoot.converged);
      bessel_cache_                                  = std::move(bessel_cache);
    }
    else {
      input_module_ptr_ = std::move(shoot.input);
    }
    shot_ = true;
  }
  return input_module_ptr_;
}

BackgroundModulePtr& Cosmology::GetBackgroundModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!background_module_ptr_) {
    background_module_ptr_ = BackgroundModulePtr(new BackgroundModule(GetInputModule()));
  }
  return background_module_ptr_;
}

ThermodynamicsModulePtr& Cosmology::GetThermodynamicsModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!thermodynamics_module_ptr_) {
    thermodynamics_module_ptr_ = ThermodynamicsModulePtr(
        new ThermodynamicsModule(GetInputModule(), GetBackgroundModule()));
  }
  return thermodynamics_module_ptr_;
}

PerturbationsModulePtr& Cosmology::GetPerturbationsModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!perturbations_module_ptr_) {
    perturbations_module_ptr_ = PerturbationsModulePtr(
        new PerturbationsModule(GetInputModule(),
                                GetBackgroundModule(),
                                GetThermodynamicsModule()));
  }
  return perturbations_module_ptr_;
}

PrimordialModulePtr& Cosmology::GetPrimordialModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!primordial_module_ptr_) {
    /** If sigma8 was input, compute local pm and nl module here, compute sigma8, update As and continue*/
    if (GetInputModule()->primordial_.sigma8 > 0) {
      auto pm = PrimordialModulePtr(
          new PrimordialModule(GetInputModule(), GetPerturbationsModule()));
      auto nl =
          NonlinearModule(GetInputModule(), GetBackgroundModule(), GetPerturbationsModule(), pm);
      const double sigma8 = nl.GetSigma8(nl.index_pk_m_);
      // The sigma8 solve rescales A_s once, before the primordial module below reads it.
      // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
      const_cast<primordial*>(&GetInputModule()->primordial_)->A_s *=
          pow(GetInputModule()->primordial_.sigma8 / sigma8, 2);
    }
    primordial_module_ptr_ = PrimordialModulePtr(
        new PrimordialModule(GetInputModule(), GetPerturbationsModule()));
  }
  return primordial_module_ptr_;
}

NonlinearModulePtr& Cosmology::GetNonlinearModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!nonlinear_module_ptr_) {
    nonlinear_module_ptr_ = NonlinearModulePtr(new NonlinearModule(GetInputModule(),
                                                                   GetBackgroundModule(),
                                                                   GetPerturbationsModule(),
                                                                   GetPrimordialModule()));
  }
  return nonlinear_module_ptr_;
}

TransferModulePtr& Cosmology::GetTransferModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!transfer_module_ptr_) {
    transfer_module_ptr_ = TransferModulePtr(new TransferModule(GetInputModule(),
                                                                GetBackgroundModule(),
                                                                GetThermodynamicsModule(),
                                                                GetPerturbationsModule(),
                                                                GetNonlinearModule(),
                                                                bessel_cache_));
  }
  return transfer_module_ptr_;
}

SpectraModulePtr& Cosmology::GetSpectraModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!spectra_module_ptr_) {
    spectra_module_ptr_ = SpectraModulePtr(new SpectraModule(GetInputModule(),
                                                             GetPerturbationsModule(),
                                                             GetPrimordialModule(),
                                                             GetNonlinearModule(),
                                                             GetTransferModule()));
  }
  return spectra_module_ptr_;
}

LensingModulePtr& Cosmology::GetLensingModule() {
  GetInputModule();  // may adopt a converged shooting build, filling the pointers below
  if (!lensing_module_ptr_) {
    lensing_module_ptr_ = LensingModulePtr(new LensingModule(GetInputModule(), GetSpectraModule()));
  }
  return lensing_module_ptr_;
}

namespace {
// Names that species answer, which callers ask for also in runs without that species.
// classy returned these values in such runs before the species answered for themselves,
// so they still do. Any other name that nothing answers is an error.
const std::map<std::string, double>& PlaceholdersForAbsentSpecies() {
  static const std::map<std::string, double> placeholders = {
      {"Omega_Lambda", 0.},
      {"Omega0_lambda", 0.},
      {"Omega0_fld", 0.},
      {"Omega0_cdm", 0.},
      {"N_dg", 0.},
      {"xi_idr", 0.},
      {"a_dark", 0.},
      {"Gamma_0_nadm", 0.},
      {"f_idr_drmd", 0.},
      {"z_dec_drmd", -1.},  // -1: no decoupling
  };
  return placeholders;
}

// Names read from the input parameters. classy answered them without building any module, and
// so does this: asking for h must not compute the background.
const DerivedParameterTable<background>& InputParameters() {
  static const DerivedParameterTable<background> table = {
      {"h", [](const background& b) { return b.h; }},
      {"H0", [](const background& b) { return b.h * 100; }},
      {"Omega_g", [](const background& b) { return b.Omega0_g; }},
      {"Omega_b", [](const background& b) { return b.Omega0_b; }},
      {"omega_b", [](const background& b) { return b.Omega0_b * b.h * b.h; }},
      {"Omega0_k", [](const background& b) { return b.Omega0_k; }},
      {"T_cmb", [](const background& b) { return b.T_cmb; }},
  };
  return table;
}

// The distances also offered in Mpc/h, under the name with "_h" appended.
bool IsHScaledDistance(const std::string& name) {
  static const std::set<std::string> names =
      {"rs_rec_h", "ds_rec_h", "ra_rec_h", "da_rec_h", "ds_d_h", "rs_d_h"};
  return names.count(name) > 0;
}

template <class Module>
auto FindReader(const std::string& name) -> const
    typename DerivedParameterTable<Module>::mapped_type* {
  const auto& table = Module::DerivedParameters();
  const auto it     = table.find(name);
  return it == table.end() ? nullptr : &it->second;
}
}  // namespace

double Cosmology::GetDerivedParameter(const std::string& name) {
  if (IsHScaledDistance(name)) {
    return GetDerivedParameter(name.substr(0, name.size() - 2)) * GetInputModule()->background_.h;
  }
  if (name == "exp_m_2_tau_As") {
    return std::exp(-2 * GetDerivedParameter("tau_reio")) * GetDerivedParameter("A_s");
  }
  if (const size_t dot = name.rfind('.'); dot != std::string::npos) {
    return GetBackgroundModule()->GetInstanceParameter(name,
                                                       name.substr(0, dot),
                                                       name.substr(dot + 1));
  }
  if (const auto it = InputParameters().find(name); it != InputParameters().end()) {
    return it->second(GetInputModule()->background_);
  }
  // A module's name means the module's quantity, whatever the species say: NEDE also answers
  // z_star (its transition), which stays reachable as "<instance>.z_star".
  if (const auto* reader = FindReader<BackgroundModule>(name)) {
    return (*reader)(*GetBackgroundModule());
  }
  if (const auto* reader = FindReader<ThermodynamicsModule>(name)) {
    return (*reader)(*GetThermodynamicsModule());
  }
  if (const auto* reader = FindReader<PrimordialModule>(name)) {
    return (*reader)(*GetPrimordialModule());
  }
  if (const auto* reader = FindReader<NonlinearModule>(name)) {
    return (*reader)(*GetNonlinearModule());
  }
  const auto answers = GetBackgroundModule()->SpeciesAnswers(name);
  if (answers.size() == 1) {
    return answers.front().second;
  }
  if (answers.size() > 1) {
    std::string keys;
    for (const auto& answer : answers) {
      keys += (keys.empty() ? "'" : ", '") + answer.first + "'";
    }
    class_stop_severe(
        "%s is answered by more than one species (%s); ask for one of them as "
        "'<species>.%s'",
        name.c_str(),
        keys.c_str(),
        name.c_str());
  }
  const auto& placeholders = PlaceholdersForAbsentSpecies();
  if (const auto it = placeholders.find(name); it != placeholders.end()) {
    return it->second;
  }
  class_stop_severe("%s was not recognized as a derived parameter", name.c_str());
}
