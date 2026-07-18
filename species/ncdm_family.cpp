#include "species/ncdm_family.h"

#include <stdexcept>
#include <string_view>

#include "errors.h"
#include "species/axion_ncdm_species.h"
#include "species/dncdm_species.h"
#include "species/greybody_ncdm_species.h"
#include "species/ncdm_interacting_species.h"
#include "species/ncdm_species.h"
#include "species/species_input.h"

namespace {

// Exactly the types whose perturbation code consults index_ap_ncdmfa:
// NCDMSpecies (and its subclasses ncdm_greybody / ncdm_axion, which inherit
// its perturbation code), NCDMInteractingSpecies, and the ncdm_decay_dr
// composite. dcdm_wdm's daughter has no fluid approximation and is absent
// on purpose.
constexpr std::string_view kFluidApproximationConsumers[] = {
    NCDMSpecies::kTypeName,
    GreyBodyNCDMSpecies::kTypeName,
    AxionNCDMSpecies::kTypeName,
    NCDMInteractingSpecies::kTypeName,
    DNCDMSpecies::kTypeName,
};

bool ConsumesFluidApproximation(std::string_view type) {
  for (const auto known : kFluidApproximationConsumers) {
    if (known == type) {
      return true;
    }
  }
  return false;
}

}  // namespace

void SynthesiseNcdmFluidApproximation(FileContent* pfc) {
  if (!pfc) {
    throw std::logic_error("SynthesiseNcdmFluidApproximation: null FileContent*");
  }

  // Reject the dot key on species that ignore the NCDM fluid approximation;
  // silently accepting it is the same bug class as ignoring it on ncdm_axion
  // (#376). Keys with no matching <name>.type are left unread for the generic
  // unused-parameter warning.
  constexpr std::string_view kSuffix = ".fluid_approximation";
  std::vector<std::string> fa_keys;
  pfc->for_each([&](const std::string& name, const std::string& /*value*/, bool /*read*/) {
    if ((name.size() > kSuffix.size()) &&
        (name.compare(name.size() - kSuffix.size(), kSuffix.size(), kSuffix) == 0)) {
      fa_keys.push_back(name);
    }
  });
  for (const auto& key : fa_keys) {
    const std::string instance = key.substr(0, key.size() - kSuffix.size());
    const auto type            = pfc->get<std::string>(instance + ".type");
    class_test_severe(type && !ConsumesFluidApproximation(*type),
                      "species '%s' (type '%s') does not use the NCDM fluid approximation; remove "
                      "'%s'",
                      instance.c_str(),
                      type->c_str(),
                      key.c_str());
  }

  std::vector<std::string> instances;
  for (const auto type : kFluidApproximationConsumers) {
    const auto found = pfc->instances_with("type", std::string(type));
    instances.insert(instances.end(), found.begin(), found.end());
  }
  SynthesiseIdenticalScalarField(pfc,
                                 instances,
                                 "fluid_approximation",
                                 "ncdm_fluid_approximation",
                                 "dot-syntax NCDM-family species");
}
