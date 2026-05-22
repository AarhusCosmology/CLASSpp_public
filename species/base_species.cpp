#include "base_species.h"

#include "perturbations.h"

double BaseSpecies::MatterRhoDelta(const perturb_vector* pv,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  if (!IsMatterSpecies())
    return 0.;
  if (collection_index_ < pv->species_layouts.size())
    return Rho(pvecback) * Delta(*pv->species_layouts[collection_index_], pv, y, pvecback, ppw);
  return Rho(pvecback) * Delta(pv, y, pvecback, ppw);
}

double BaseSpecies::MatterRhoPlusPTheta(const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  if (!IsMatterSpecies())
    return 0.;
  if (collection_index_ < pv->species_layouts.size())
    return (Rho(pvecback) + P(pvecback)) *
           Theta(*pv->species_layouts[collection_index_], pv, y, pvecback, ppw);
  return (Rho(pvecback) + P(pvecback)) * Theta(pv, y, pvecback, ppw);
}
