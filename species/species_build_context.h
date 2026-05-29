#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "base_species.h"

class FileContent;
struct background;
struct precision;
class BackgroundModule;
struct NcdmSettings;
class SpeciesCollection;

/**
 * Pre-resolved Omega0 values for the coupled-species cluster (CDM, IDM_DR,
 * IDR, IDM_DRMD, IDR_DRMD, DCDM_DR). These cannot be parsed independently by
 * each species: f_idm_dr/f_idm_drmd subtract from CDM, IDR is derived from
 * T_idr*Omega0_g, and the synchronous-gauge CDM minimum applies after all of
 * those adjustments. InputModule computes the budget once via
 * ReadCoupledOmegaBudget() before ConstructSpecies runs and passes it through
 * SpeciesBuildContext; each affected CreateAll reads its slot.
 *
 * std::nullopt means "absent" — the corresponding CreateAll returns an empty
 * result. A present-but-zero value also means absent (preserved for legacy
 * parser semantics like Omega_cdm=0 with no synchronous-gauge minimum kick).
 */
struct SpeciesOmegaBudget {
  std::optional<double> cdm;
  std::optional<double> dcdmdr;
  std::optional<double> idm_dr;
  std::optional<double> idr;
  std::optional<double> idm_drmd;
  std::optional<double> idr_drmd;
};

/**
 * Inputs every species' static CreateAll factory needs.
 * Bundled into one struct to keep factory signatures uniform.
 */
struct SpeciesBuildContext {
  FileContent* pfc;
  const background* pba;
  const precision* ppr;
  const NcdmSettings* ncdm_settings;  // non-null
  const BackgroundModule* bgm;        // nullptr at species-construction time
  const SpeciesCollection* all_species = nullptr;
  // ^ ConstructSpecies passes a non-null pointer for the whole build, but the
  // collection is filled incrementally — earlier factories in Pass 1 see fewer
  // entries than later ones. nullptr in pre-construction / standalone contexts
  // (e.g. DoShooting's per-species shooting-guess loop).

  // Pre-resolved Omega budget for the coupled-species cluster.  Populated by
  // InputModule before ConstructSpecies runs and threaded through to DoShooting's
  // per-species guess loop; default nullptr for standalone callers (tests etc.).
  const SpeciesOmegaBudget* omega_budget = nullptr;

  // Set by ConstructSpecies when running Pass 2 (closure species). Carries the
  // budget-closure value Omega0 that the closure species must adopt instead of
  // reading its Omega from the input file. std::nullopt for Pass 1 and for
  // non-closure species.
  std::optional<double> omega0_closure_override;
};

/**
 * One species sector produced by a CreateAll factory.
 * `key` is the SpeciesCollection insertion key.
 */
struct Named {
  std::string key;
  std::unique_ptr<BaseSpecies> species;
};
