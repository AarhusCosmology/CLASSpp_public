#include "species/species_input.h"

#include <cstdlib>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "species/all_species.h"

SpeciesInput::SpeciesInput(FileContent* pfc, std::string instance_name)
    : pfc_(pfc), instance_name_(std::move(instance_name)) {
  if (!pfc_) {
    throw std::invalid_argument("SpeciesInput: null FileContent*");
  }
  if (instance_name_.empty()) {
    throw std::invalid_argument("SpeciesInput: empty instance name");
  }
}

std::string SpeciesInput::qualify(const std::string& field) const {
  return instance_name_ + "." + field;
}

namespace {

int ParseIntValue(const std::string& value, const std::string& field_name) {
  char* end   = nullptr;
  long parsed = std::strtol(value.c_str(), &end, 10);
  if ((end == value.c_str()) || (*end != '\0')) {
    throw std::invalid_argument("dot-syntax field '" + field_name + "' must be an integer, got '" +
                                value + "'");
  }
  return static_cast<int>(parsed);
}

struct DotFieldMap {
  std::string_view dot;
  std::string_view legacy;
};

struct SingleInstanceSpec {
  std::string_view type;
  std::vector<DotFieldMap> fields;
};

// Module-owned: which named species are single-instance, and how each clean
// dot-field maps to its legacy key. Keyed by the species' canonical kTypeName.
// ScalarField is intentionally absent (reserved for a multi-instance future).
const std::vector<SingleInstanceSpec>& SingleInstanceTable() {
  static const std::vector<SingleInstanceSpec> table = {
      {PhotonsSpecies::kTypeName, {{"Omega", "Omega_g"}, {"T_cmb", "T_cmb"}}},
      {BaryonsSpecies::kTypeName, {{"Omega", "Omega_b"}}},
      {CDMSpecies::kTypeName, {{"Omega", "Omega_cdm"}}},
      {UltraRelativisticSpecies::kTypeName,
       {{"N", "N_ur"}, {"Omega", "Omega_ur"}, {"omega", "omega_ur"}}},
      {LambdaSpecies::kTypeName, {{"Omega", "Omega_Lambda"}}},
      {FluidSpecies::kTypeName,
       {{"Omega", "Omega_fld"},
        {"w0", "w0_fld"},
        {"wa", "wa_fld"},
        {"cs2", "cs2_fld"},
        {"Omega_EDE", "Omega_EDE"},
        {"c_gamma_over_c", "c_gamma_over_c_fld"},
        {"use_ppf", "use_ppf"},
        {"equation_of_state", "fluid_equation_of_state"}}},
  };
  return table;
}

}  // namespace

std::vector<std::string> CollectInstanceFieldValues(FileContent* pfc,
                                                    const std::vector<std::string>& instances,
                                                    const std::string& field) {
  std::vector<std::string> values(instances.size());
  for (size_t i = 0; i < instances.size(); ++i) {
    SpeciesInput input(pfc, instances[i]);
    if (auto v = input.get<std::string>(field))
      values[i] = *v;
  }
  return values;
}

bool AnyInstanceFieldValue(const std::vector<std::string>& values) {
  for (const auto& value : values) {
    if (!value.empty()) {
      return true;
    }
  }
  return false;
}

bool SynthesiseIdenticalScalarField(FileContent* pfc,
                                    const std::vector<std::string>& instances,
                                    const std::string& dot_field,
                                    const std::string& legacy_key,
                                    const std::string& species_description) {
  const std::vector<std::string> values = CollectInstanceFieldValues(pfc, instances, dot_field);
  if (!AnyInstanceFieldValue(values)) {
    return false;
  }

  for (const auto& value : values) {
    if (value.empty()) {
      throw std::invalid_argument("dot-syntax field '" + dot_field +
                                  "' must be specified for all " + species_description +
                                  " if provided for any of them");
    }
  }

  const std::string& first = values.front();
  for (const auto& value : values) {
    if (value != first) {
      throw std::invalid_argument(dot_field + " must be identical for all " + species_description);
    }
  }

  auto existing_value = pfc->get<std::string>(legacy_key);
  if (existing_value.has_value() && *existing_value != first) {
    throw std::invalid_argument("input sets both '" + legacy_key + "' and dot-syntax field '" +
                                dot_field + "' with different values");
  }

  pfc->set(legacy_key, first);
  return true;
}

void TranslateSingleInstanceDotSyntax(FileContent* pfc) {
  if (!pfc) {
    throw std::invalid_argument("TranslateSingleInstanceDotSyntax: null FileContent*");
  }
  for (const auto& spec : SingleInstanceTable()) {
    const std::string type(spec.type);
    const auto instances = pfc->instances_with("type", type);
    if (instances.empty()) {
      continue;
    }
    if (instances.size() > 1) {
      std::string joined;
      for (size_t i = 0; i < instances.size(); ++i) {
        joined += (i ? ", " : "") + instances[i];
      }
      throw std::invalid_argument("species type '" + type + "' is single-instance but was given " +
                                  std::to_string(instances.size()) + " times (" + joined + ")");
    }
    SpeciesInput input(pfc, instances.front());
    (void) input.get<std::string>("type");  // consume "<name>.type"
    for (const auto& fm : spec.fields) {
      const std::string dot(fm.dot);
      auto value = input.get<std::string>(dot);  // consumes "<name>.<dot>" if present
      if (!value) {
        continue;
      }
      const std::string legacy(fm.legacy);
      auto existing = pfc->get<std::string>(legacy);
      if (existing && *existing != *value) {
        throw std::invalid_argument("input sets both legacy key '" + legacy + "' and dot-syntax '" +
                                    instances.front() + "." + dot + "' with different values");
      }
      pfc->set(
          legacy,
          *value);  // set() resets the read flag so the downstream consumer still sees this key
    }
  }
}
