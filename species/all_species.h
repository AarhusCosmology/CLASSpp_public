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

#include "axion_ncdm_species.h"
#include "baryons.h"
#include "base_species.h"
#include "cdm.h"
#include "dark_radiation_species.h"
#include "dcdm.h"
#include "dcdm_dr_species.h"
#include "dcdm_wdm_species.h"
#include "dncdm_dr_species.h"
#include "dncdm_inv_species.h"
#include "dr_psd_species.h"
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
#include "type3_species.h"
#include "ultra_relativistic.h"

struct SpeciesFactoryEntry {
  std::string_view name;
  std::vector<Named> (*create_all)(const SpeciesBuildContext&);
};

// NOTE: each entry name is the species' kTypeName (the .ini `type =` value).
// This is deliberately distinct from a species' collection registration key
// (CamelCase, e.g. all_species_.at("Lambda")), which lives in each CreateAll.
inline constexpr std::array kAllSpeciesFactories = {
    SpeciesFactoryEntry{PhotonsSpecies::kTypeName, &PhotonsSpecies::CreateAll},
    SpeciesFactoryEntry{BaryonsSpecies::kTypeName, &BaryonsSpecies::CreateAll},
    SpeciesFactoryEntry{CDMSpecies::kTypeName, &CDMSpecies::CreateAll},
    SpeciesFactoryEntry{UltraRelativisticSpecies::kTypeName, &UltraRelativisticSpecies::CreateAll},
    SpeciesFactoryEntry{DCDM_DR_Species::kTypeName, &DCDM_DR_Species::CreateAll},
    SpeciesFactoryEntry{NCDMSpecies::kTypeName, &NCDMSpecies::CreateAll},
    SpeciesFactoryEntry{GreyBodyNCDMSpecies::kTypeName, &GreyBodyNCDMSpecies::CreateAll},
    SpeciesFactoryEntry{AxionNCDMSpecies::kTypeName, &AxionNCDMSpecies::CreateAll},
    SpeciesFactoryEntry{DNCDMSpecies::kTypeName, &DNCDM_DR_Species::CreateAll},
    SpeciesFactoryEntry{DCDM_WDM_Species::kTypeName, &DCDM_WDM_Species::CreateAll},
    SpeciesFactoryEntry{DrPsdSpecies::kTypeName, &DrPsdSpecies::CreateAll},
    SpeciesFactoryEntry{NCDMInteractingSpecies::kTypeName, &NCDMInteractingSpecies::CreateAll},
    SpeciesFactoryEntry{IDM_DR_IDR_Species::kTypeName, &IDM_DR_IDR_Species::CreateAll},
    SpeciesFactoryEntry{IDM_DRMD_IDR_DRMD_Species::kTypeName,
                        &IDM_DRMD_IDR_DRMD_Species::CreateAll},
    SpeciesFactoryEntry{LambdaSpecies::kTypeName, &LambdaSpecies::CreateAll},
    SpeciesFactoryEntry{FluidSpecies::kTypeName, &FluidSpecies::CreateAll},
    SpeciesFactoryEntry{ScalarFieldSpecies::kTypeName, &ScalarFieldSpecies::CreateAll},
    SpeciesFactoryEntry{Type3Species::kTypeName, &Type3Species::CreateAll},
};

/** Maps the closure-species enum to the matching factory entry name. */
inline constexpr std::string_view ClosureSpeciesName(ClosureSpecies cs) {
  switch (cs) {
    case ClosureSpecies::Lambda:
      return LambdaSpecies::kTypeName;
    case ClosureSpecies::Fluid:
      return FluidSpecies::kTypeName;
    case ClosureSpecies::ScalarField:
      return ScalarFieldSpecies::kTypeName;
    case ClosureSpecies::None:
      return "";
  }
  return "";
}
