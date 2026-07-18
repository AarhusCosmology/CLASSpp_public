#include "base_species.h"

#include "ncdm_base_species.h"

BudgetBucket BudgetBucketOf(const BaseSpecies& s) {
  // NCDM-family species carry EnergyType::Other (they tally into both rho_r and
  // rho_m), so route them by family, not energy_type. Same dynamic_cast idiom as
  // GetNcdmSpecies.
  if (dynamic_cast<const NCDMBaseSpecies*>(&s))
    return BudgetBucket::Ncdm;
  switch (s.energy_type()) {
    case BaseSpecies::EnergyType::Radiation:
      return BudgetBucket::Radiation;
    case BaseSpecies::EnergyType::Matter:
      return BudgetBucket::NonRelativistic;
    case BaseSpecies::EnergyType::DarkEnergy:
    case BaseSpecies::EnergyType::Other:
      return BudgetBucket::Other;
  }
  return BudgetBucket::Other;  // unreachable; silences -Wreturn-type
}
