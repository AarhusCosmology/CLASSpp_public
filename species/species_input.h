#pragma once

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "errors.h"
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

  /** Yes/no flag on the repo's convention: the value must begin with 'y'/'Y' or
   *  'n'/'N', absent falls back to @p fallback, and anything else is REJECTED --
   *  same rule input_module.cpp applies to its own yes/no keys. Rejecting rather
   *  than reading an unrecognised value as "no" is the point: these flags select
   *  physics, and a typo that silently turns one off is not something a run reports.
   *  Structural (decidable from the text, not from a value a sampler varies), hence
   *  severe. */
  bool get_flag(const std::string& field, bool fallback) const {
    auto v = get<std::string>(field);
    if (!v)
      return fallback;
    const char c = v->empty() ? '\0' : (*v)[0];
    class_test_severe(c != 'y' && c != 'Y' && c != 'n' && c != 'N',
                      "species '%s': field '%s' must begin with 'y' or 'n', not '%s'",
                      instance_name_.c_str(),
                      field.c_str(),
                      v->c_str());
    return c == 'y' || c == 'Y';
  }

  /** Required accessor: throws std::invalid_argument if the field is absent. */
  template <class T>
  T require(const std::string& field) const {
    if (auto v = get<T>(field))
      return *v;
    class_stop_severe("species '%s': missing required field '%s'",
                      instance_name_.c_str(),
                      field.c_str());
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

/**
 * Translate dot-syntax for SINGLE-INSTANCE legacy species into their legacy
 * keys, in place. For each recognised single-instance type T (photons, baryons,
 * cdm, ur, lambda, fluid), finds the at-most-one instance N with "N.type == T",
 * rewrites each known "N.<dot_field>" to its legacy key via FileContent::set,
 * and marks the dot entries read. The instance name N is discarded, so output
 * is identical to the legacy-key form. Throws std::invalid_argument if a type
 * appears more than once, or if a translated legacy key is already present with
 * a different value. Unknown "N.<field>" entries are left untouched (and unread,
 * so the usual unrecognised-parameter warning still fires).
 *
 * Copies at most once per FileContent (FileContent::single_instance_translated).
 * A reparse of the same content or a copy of it consumes "N.type" and the dot
 * fields but leaves the legacy keys alone, so a later override of a legacy key
 * wins: pk_eq's effective w0_fld/wa_fld on a dot-syntax fluid used to abort here
 * as "sets both legacy key and dot-syntax" (see
 * docs/superpowers/specs/2026-09-22-dot-syntax-translate-once-design.md).
 */
void TranslateSingleInstanceDotSyntax(FileContent* pfc);

/**
 * Reject every dot-syntax species instance that no factory built (#430).
 *
 * Call once, after every factory in kAllSpeciesFactories has run. Rejects:
 *   - "N.type" that no factory consumed: a type name that is not a species type,
 *     a species type that is configured by legacy keys only, or an N that is not
 *     a legal instance name;
 *   - "N.<field>" with no "N.type" line at all (the type line misspelled or
 *     forgotten).
 * Either way the species would otherwise be silently absent while the run
 * reports success. Relies on the factory contract stated in all_species.h: a
 * factory reads "N.type" for exactly the instances it builds. Structural, hence
 * severe; every offending instance is named in one message.
 *
 * Unread fields of an instance that WAS built are not checked here.
 *
 * Design: docs/superpowers/specs/2026-09-21-unbuilt-species-instances-design.md
 */
void RejectUnbuiltSpeciesInstances(const FileContent& fc);
