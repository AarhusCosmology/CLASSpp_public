#include "base_species.h"

#include "perturbations.h"

double BaseSpecies::MatterRhoDelta(const perturb_vector* pv,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  if (!IsMatterSpecies())
    return 0.;
  return Rho(pvecback) * Delta(*pv->species_layouts.at(collection_index_), pv, y, pvecback, ppw);
}

double BaseSpecies::MatterRhoPlusPTheta(const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  if (!IsMatterSpecies())
    return 0.;
  return (Rho(pvecback) + P(pvecback)) *
         Theta(*pv->species_layouts.at(collection_index_), pv, y, pvecback, ppw);
}
