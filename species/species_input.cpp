#include "species/species_input.h"

#include <cstdlib>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "errors.h"
#include "species/all_species.h"

SpeciesInput::SpeciesInput(FileContent* pfc, std::string instance_name)
    : pfc_(pfc), instance_name_(std::move(instance_name)) {
  if (!pfc_) {
    throw std::logic_error("SpeciesInput: null FileContent*");
  }
  class_test_severe(instance_name_.empty(), "SpeciesInput: empty instance name");
}

std::string SpeciesInput::qualify(const std::string& field) const {
  return instance_name_ + "." + field;
}

namespace {

int ParseIntValue(const std::string& value, const std::string& field_name) {
  char* end   = nullptr;
  long parsed = std::strtol(value.c_str(), &end, 10);
  class_test_severe((end == value.c_str()) || (*end != '\0'),
                    "dot-syntax field '%s' must be an integer, got '%s'",
                    field_name.c_str(),
                    value.c_str());
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
    class_test_severe(value.empty(),
                      "dot-syntax field '%s' must be specified for all %s if provided for any of "
                      "them",
                      dot_field.c_str(),
                      species_description.c_str());
  }

  const std::string& first = values.front();
  for (const auto& value : values) {
    class_test_severe(value != first,
                      "%s must be identical for all %s",
                      dot_field.c_str(),
                      species_description.c_str());
  }

  auto existing_value = pfc->get<std::string>(legacy_key);
  class_test_severe(existing_value.has_value() && *existing_value != first,
                    "input sets both '%s' and dot-syntax field '%s' with different values",
                    legacy_key.c_str(),
                    dot_field.c_str());

  pfc->set(legacy_key, first);
  return true;
}

void TranslateSingleInstanceDotSyntax(FileContent* pfc) {
  if (!pfc) {
    throw std::logic_error("TranslateSingleInstanceDotSyntax: null FileContent*");
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
      class_stop_severe("species type '%s' is single-instance but was given %zu times (%s)",
                        type.c_str(),
                        instances.size(),
                        joined.c_str());
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
      class_test_severe(existing && *existing != *value,
                        "input sets both legacy key '%s' and dot-syntax '%s.%s' with different "
                        "values",
                        legacy.c_str(),
                        instances.front().c_str(),
                        dot.c_str());
      pfc->set(
          legacy,
          *value);  // set() resets the read flag so the downstream consumer still sees this key
    }
  }
}
