#pragma once
#include <memory>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "dncdm_decay_radiation_species.h"
#include "dncdm_species.h"
#include "parser.h"

class BackgroundModule;

/**
 * DNCDM_DR_Species: composite for one flavor of Decaying Non-Cold Dark Matter + its decay radiation.
 */
class DNCDM_DR_Species : public CompositeSpecies {
 public:
  // Takes ownership of a pre-built DNCDMSpecies (from DNCDMSpecies::CreateAll)
  DNCDM_DR_Species(std::unique_ptr<DNCDMSpecies> dncdm,
                   const background* pba,
                   const BackgroundModule* bgm);

  static std::vector<std::unique_ptr<DNCDM_DR_Species>> CreateAll(FileContent* pfc,
                                                                  const NcdmSettings& settings,
                                                                  const background* pba,
                                                                  const BackgroundModule* bgm);

  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;

  // Override to add DNCDM->DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    dncdm_->WriteBackgroundColumnTitles(w);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    dncdm_->WriteBackgroundData(pvecback, w);
  }

  DNCDMSpecies& dncdm() {
    return *dncdm_;
  }
  DNCDM_DecayRadiationSpecies& dr() {
    return *dr_sp_;
  }

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

 private:
  int ncdm_id_;
  DNCDMSpecies* dncdm_                = nullptr;
  DNCDM_DecayRadiationSpecies* dr_sp_ = nullptr;
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;
};
