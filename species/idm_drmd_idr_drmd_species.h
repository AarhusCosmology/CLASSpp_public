#pragma once
#include "composite_species.h"
#include "idm_drmd.h"
#include "idr_drmd.h"
#include "background.h"

class IDM_DRMD_IDR_DRMD_Species : public CompositeSpecies {
public:
  explicit IDM_DRMD_IDR_DRMD_Species(const background& pba);

  IDM_DRMDSpecies&       idm_drmd()       { return *idm_drmd_; }
  IDR_DRMDSpecies&       idr_drmd()       { return *idr_drmd_; }
  const IDM_DRMDSpecies& idm_drmd() const { return *idm_drmd_; }
  const IDR_DRMDSpecies& idr_drmd() const { return *idr_drmd_; }

protected:
  void AddCouplingDerivs(double tau, const double* y, double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

private:
  IDM_DRMDSpecies* idm_drmd_ = nullptr;
  IDR_DRMDSpecies* idr_drmd_ = nullptr;
  const background& pba_;
};
