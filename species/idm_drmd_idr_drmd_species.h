#pragma once
#include "background.h"
#include "composite_species.h"
#include "idm_drmd.h"
#include "idr_drmd.h"

class IDM_DRMD_IDR_DRMD_Species : public CompositeSpecies {
 public:
  explicit IDM_DRMD_IDR_DRMD_Species(const background& pba);

  IDM_DRMDSpecies& idm_drmd() {
    return *idm_drmd_;
  }
  IDR_DRMDSpecies& idr_drmd() {
    return *idr_drmd_;
  }
  const IDM_DRMDSpecies& idm_drmd() const {
    return *idm_drmd_;
  }
  const IDR_DRMDSpecies& idr_drmd() const {
    return *idr_drmd_;
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

  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
    CompositeSpecies::SetBackgroundModule(bgm);
  }

 private:
  IDM_DRMDSpecies* idm_drmd_ = nullptr;
  IDR_DRMDSpecies* idr_drmd_ = nullptr;
  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
};
