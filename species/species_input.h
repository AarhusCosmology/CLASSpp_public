#pragma once

#include <optional>
#include <stdexcept>
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

  /** Typed accessor: the value if "<instance>.<field>" is present (marking it
   *  read), or std::nullopt if absent. Delegates to FileContent::get<T>. */
  template <class T>
  std::optional<T> get(const std::string& field) const {
    return pfc_->get<T>(qualify(field));
  }

  /** Read-with-default: get<T>(field) if present, else @p fallback. */
  template <class T>
  T get_or(const std::string& field, T fallback) const {
    if (auto v = get<T>(field))
      return *v;
    return fallback;
  }

  /** Required accessor: throws std::invalid_argument if the field is absent. */
  template <class T>
  T require(const std::string& field) const {
    if (auto v = get<T>(field))
      return *v;
    throw std::invalid_argument("species '" + instance_name_ + "': missing required field '" +
                                field + "'");
  }

 private:
  std::string qualify(const std::string& field) const;

  FileContent* pfc_;
  std::string instance_name_;
};

std::vector<std::string> CollectInstanceFieldValues(FileContent* pfc,
                                                    const std::vector<std::string>& instances,
                                                    const std::string& field);

bool AnyInstanceFieldValue(const std::vector<std::string>& values);

bool SynthesiseIdenticalScalarField(FileContent* pfc,
                                    const std::vector<std::string>& instances,
                                    const std::string& dot_field,
                                    const std::string& legacy_key,
                                    const std::string& species_description);
