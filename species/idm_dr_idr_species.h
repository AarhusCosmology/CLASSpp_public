#pragma once
#include "composite_species.h"
#include "idm_dr.h"
#include "idr.h"
#include "background.h"

/**
 * IDM_DR_IDR_Species: composite for interacting dark matter + interacting dark radiation.
 *
 * Children handle free-streaming terms; this composite's AddCouplingDerivs
 * adds the momentum-exchange and TCA terms that couple IDM_DR to IDR.
 */
class IDM_DR_IDR_Species : public CompositeSpecies {
public:
  explicit IDM_DR_IDR_Species(const background& pba);

  IDM_DRSpecies&       idm_dr()       { return *idm_dr_; }
  IDRSpecies&          idr()          { return *idr_; }
  const IDM_DRSpecies& idm_dr() const { return *idm_dr_; }
  const IDRSpecies&    idr()    const { return *idr_; }

protected:
  void AddCouplingDerivs(double tau, const double* y, double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

private:
  IDM_DRSpecies* idm_dr_ = nullptr;
  IDRSpecies*    idr_    = nullptr;
  const background& pba_;
};
