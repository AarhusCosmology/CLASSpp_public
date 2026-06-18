#include "species/species_input.h"

#include <cstdlib>
#include <stdexcept>

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
