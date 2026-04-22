#pragma once

#include <string>
#include <vector>

#include "parser.h"

/**
 * Per-instance read wrapper over FileContent for dot-syntax species input.
 * Prefixes every key with "<instance_name>." before delegating to the
 * underlying FileContent. Successful reads mark the fully-qualified key as
 * read so unread_parameters() remains accurate.
 */
class SpeciesInput {
 public:
  SpeciesInput(FileContent* pfc, std::string instance_name);

  const std::string& instance_name() const {
    return instance_name_;
  }

  bool read_double(const std::string& field, double& out);
  bool read_int(const std::string& field, int& out);
  bool read_string(const std::string& field, std::string& out);
  bool read_list_of_doubles(const std::string& field, std::vector<double>& out);

  double required_double(const std::string& field);
  int required_int(const std::string& field);
  std::string required_string(const std::string& field);

 private:
  std::string qualify(const std::string& field) const;

  FileContent* pfc_;
  std::string instance_name_;
};

std::vector<std::string> CollectInstanceFieldValues(FileContent* pfc,
                                                    const std::vector<std::string>& instances,
                                                    const std::string& field);

bool AnyInstanceFieldValue(const std::vector<std::string>& values);

std::string CsvWithDefaults(const std::vector<std::string>& values,
                            const std::string& default_value);

std::string CsvForPsdFilenames(const std::vector<std::string>& use_psd_file_values,
                               const std::vector<std::string>& filename_values,
                               const std::string& dot_flag_field,
                               const std::string& dot_filename_field,
                               const std::string& legacy_filename_key);

bool SynthesiseIdenticalScalarField(FileContent* pfc,
                                    const std::vector<std::string>& instances,
                                    const std::string& dot_field,
                                    const std::string& legacy_key,
                                    const std::string& species_description);
