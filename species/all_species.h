#pragma once
/**
 * Aggregator: include every species class and declare the constexpr factory
 * list iterated by InputModule::ConstructSpecies.
 *
 * Adding a new species: write the class, include its header here, add one
 * row to kAllSpeciesFactories below.
 */
#include <array>
#include <string_view>

#include "baryons.h"
#include "base_species.h"
#include "cdm.h"
#include "dark_radiation_species.h"
#include "dcdm.h"
#include "dcdm_dr_species.h"
#include "dncdm_dr_species.h"
#include "fluid.h"
#include "greybody_ncdm_species.h"
#include "idm_dr.h"
#include "idm_dr_idr_species.h"
#include "idm_drmd.h"
#include "idm_drmd_idr_drmd_species.h"
#include "idr.h"
#include "idr_drmd.h"
#include "lambda.h"
#include "ncdm_interacting_species.h"
#include "ncdm_species.h"
#include "photons.h"
#include "scalar_field.h"
#include "species_build_context.h"
#include "ultra_relativistic.h"

struct SpeciesFactoryEntry {
  std::string_view name;
  std::vector<Named> (*create_all)(const SpeciesBuildContext&);
};

inline constexpr std::array kAllSpeciesFactories = {
    SpeciesFactoryEntry{"Photons", &PhotonsSpecies::CreateAll},
    SpeciesFactoryEntry{"Baryons", &BaryonsSpecies::CreateAll},
    SpeciesFactoryEntry{"CDM", &CDMSpecies::CreateAll},
    SpeciesFactoryEntry{"UR", &UltraRelativisticSpecies::CreateAll},
    SpeciesFactoryEntry{"DCDM_DR", &DCDM_DR_Species::CreateAll},
    SpeciesFactoryEntry{"NCDM", &NCDMSpecies::CreateAll},
    SpeciesFactoryEntry{"NCDMGreyBody", &GreyBodyNCDMSpecies::CreateAll},
    SpeciesFactoryEntry{"DNCDM_DR", &DNCDM_DR_Species::CreateAll},
    SpeciesFactoryEntry{"NCDMInt", &NCDMInteractingSpecies::CreateAll},
    SpeciesFactoryEntry{"IDM_DR_IDR", &IDM_DR_IDR_Species::CreateAll},
    SpeciesFactoryEntry{"IDM_DRMD_IDR_DRMD", &IDM_DRMD_IDR_DRMD_Species::CreateAll},
    SpeciesFactoryEntry{"Lambda", &LambdaSpecies::CreateAll},
    SpeciesFactoryEntry{"Fluid", &FluidSpecies::CreateAll},
    SpeciesFactoryEntry{"ScalarField", &ScalarFieldSpecies::CreateAll},
};

/** Maps the closure-species enum to the matching factory entry name. */
inline constexpr std::string_view ClosureSpeciesName(ClosureSpecies cs) {
  switch (cs) {
    case ClosureSpecies::Lambda:
      return "Lambda";
    case ClosureSpecies::Fluid:
      return "Fluid";
    case ClosureSpecies::ScalarField:
      return "ScalarField";
    case ClosureSpecies::None:
      return "";
  }
  return "";
}
