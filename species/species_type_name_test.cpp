// Pins each species' canonical type string and the kAllSpeciesFactories
// registry to the agreed set, so a rename or a missing override is caught.
#include <cassert>
#include <set>
#include <string>

#include "species/all_species.h"

int main() {
  // Per-species constants (the inherited-shadowing cases matter: GreyBody and
  // NCDMInt derive from NCDMSpecies and must declare their own kTypeName).
  assert(PhotonsSpecies::kTypeName == "photons");
  assert(BaryonsSpecies::kTypeName == "baryons");
  assert(CDMSpecies::kTypeName == "cdm");
  assert(UltraRelativisticSpecies::kTypeName == "ur");
  assert(DCDM_DR_Species::kTypeName == "dcdm_dr");
  assert(NCDMSpecies::kTypeName == "ncdm_standard");
  assert(GreyBodyNCDMSpecies::kTypeName == "ncdm_greybody");
  assert(DNCDMSpecies::kTypeName == "ncdm_decay_dr");
  assert(NCDMInteractingSpecies::kTypeName == "ncdm_self_interacting");
  assert(IDM_DR_IDR_Species::kTypeName == "idm_dr_idr");
  assert(IDM_DRMD_IDR_DRMD_Species::kTypeName == "idm_drmd_idr_drmd");
  assert(LambdaSpecies::kTypeName == "lambda");
  assert(FluidSpecies::kTypeName == "fluid");
  assert(ScalarFieldSpecies::kTypeName == "scalar_field");
  assert(Type3Species::kTypeName == "cdm_scf_momentum");

  // The factory registry carries exactly these type strings.
  std::set<std::string> names;
  for (const auto& e : kAllSpeciesFactories) {
    names.insert(std::string(e.name));
  }
  const std::set<std::string> expected = {
      "photons",
      "baryons",
      "cdm",
      "ur",
      "dcdm_dr",
      "ncdm_standard",
      "ncdm_greybody",
      "ncdm_decay_dr",
      "ncdm_self_interacting",
      "idm_dr_idr",
      "idm_drmd_idr_drmd",
      "lambda",
      "fluid",
      "scalar_field",
      "cdm_scf_momentum",
  };
  assert(names == expected);

  // ClosureSpeciesName routes through the same strings.
  assert(ClosureSpeciesName(ClosureSpecies::Lambda) == "lambda");
  assert(ClosureSpeciesName(ClosureSpecies::Fluid) == "fluid");
  assert(ClosureSpeciesName(ClosureSpecies::ScalarField) == "scalar_field");
  return 0;
}
