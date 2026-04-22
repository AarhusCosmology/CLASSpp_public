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

bool SpeciesInput::read_double(const std::string& field, double& out) {
  return pfc_->read_double(qualify(field), out);
}
bool SpeciesInput::read_int(const std::string& field, int& out) {
  return pfc_->read_int(qualify(field), out);
}
bool SpeciesInput::read_string(const std::string& field, std::string& out) {
  return pfc_->read_string(qualify(field), out);
}
bool SpeciesInput::read_list_of_doubles(const std::string& field, std::vector<double>& out) {
  return pfc_->read_list_of_doubles(qualify(field), out);
}

double SpeciesInput::required_double(const std::string& field) {
  double v = 0.;
  if (!read_double(field, v)) {
    throw std::invalid_argument("species '" + instance_name_ + "': missing required field '" +
                                field + "'");
  }
  return v;
}
int SpeciesInput::required_int(const std::string& field) {
  int v = 0;
  if (!read_int(field, v)) {
    throw std::invalid_argument("species '" + instance_name_ + "': missing required field '" +
                                field + "'");
  }
  return v;
}
std::string SpeciesInput::required_string(const std::string& field) {
  std::string v;
  if (!read_string(field, v)) {
    throw std::invalid_argument("species '" + instance_name_ + "': missing required field '" +
                                field + "'");
  }
  return v;
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
    input.read_string(field, values[i]);
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

std::string CsvWithDefaults(const std::vector<std::string>& values,
                            const std::string& default_value) {
  std::string csv;
  for (size_t i = 0; i < values.size(); ++i) {
    if (i) {
      csv += ", ";
    }
    csv += values[i].empty() ? default_value : values[i];
  }
  return csv;
}

std::string CsvForPsdFilenames(const std::vector<std::string>& use_psd_file_values,
                               const std::vector<std::string>& filename_values,
                               const std::string& dot_flag_field,
                               const std::string& dot_filename_field,
                               const std::string& legacy_filename_key) {
  std::string csv;
  for (size_t i = 0; i < use_psd_file_values.size(); ++i) {
    const std::string flag = use_psd_file_values[i].empty() ? "0" : use_psd_file_values[i];
    if (ParseIntValue(flag, dot_flag_field) == 0) {
      continue;
    }
    if (filename_values[i].empty()) {
      throw std::invalid_argument("dot-syntax field '" + dot_filename_field +
                                  "' must be provided whenever '" + dot_flag_field +
                                  "' is nonzero while synthesising legacy key '" +
                                  legacy_filename_key + "'");
    }
    if (!csv.empty()) {
      csv += ", ";
    }
    csv += filename_values[i];
  }
  return csv;
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

  std::string existing_value;
  if (pfc->read_string(legacy_key, existing_value) && existing_value != first) {
    throw std::invalid_argument("input sets both '" + legacy_key + "' and dot-syntax field '" +
                                dot_field + "' with different values");
  }

  pfc->set(legacy_key, first);
  return true;
}
