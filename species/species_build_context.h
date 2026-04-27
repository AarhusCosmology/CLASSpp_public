#pragma once

#include <memory>
#include <string>
#include <vector>

#include "base_species.h"

class FileContent;
struct background;
struct precision;
class BackgroundModule;
struct NcdmSettings;

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
};

/**
 * One species sector produced by a CreateAll factory.
 * `key` is the SpeciesCollection insertion key.
 */
struct Named {
  std::string key;
  std::unique_ptr<BaseSpecies> species;
};
