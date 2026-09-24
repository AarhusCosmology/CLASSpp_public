#include "base_species.h"

#include "ncdm_base_species.h"
#include "perturbations.h"

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

int BaseSpecies::ApproximationRegime(const perturb_workspace* ppw, int which) const {
  return ppw->approx[ppw->species_approx_index[collection_index_] + which];
}

double BaseSpecies::BackgroundDensityOverH0Sq(double a, double /*H0*/) const {
  switch (energy_type()) {
    case EnergyType::Radiation:
      return GetOmega0() / (a * a * a * a);
    case EnergyType::Matter:
      return GetOmega0() / (a * a * a);
    case EnergyType::DarkEnergy:
      return GetOmega0();
    case EnergyType::Other:
      break;
  }
  class_stop_severe(
      "species '%s' cannot give its density before the background is solved, which NEDE "
      "needs to place its transition",
      name().c_str());
}
