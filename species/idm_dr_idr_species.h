#pragma once
#include "background.h"
#include "composite_species.h"
#include "idm_dr.h"
#include "idr.h"

/**
 * IDM_DR_IDR_Species: composite for interacting dark matter + interacting dark radiation.
 *
 * Children handle free-streaming terms; this composite's AddCouplingDerivs
 * adds the momentum-exchange and TCA terms that couple IDM_DR to IDR.
 */
class IDM_DR_IDR_Species : public CompositeSpecies {
 public:
  explicit IDM_DR_IDR_Species(const background& pba);

  IDM_DRSpecies& idm_dr() {
    return *idm_dr_;
  }
  IDRSpecies& idr() {
    return *idr_;
  }
  const IDM_DRSpecies& idm_dr() const {
    return *idm_dr_;
  }
  const IDRSpecies& idr() const {
    return *idr_;
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      enum file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

 private:
  IDM_DRSpecies* idm_dr_ = nullptr;
  IDRSpecies* idr_       = nullptr;
  const background& pba_;
};
