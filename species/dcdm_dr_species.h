#pragma once
#include "background.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dcdm.h"

class BackgroundModule;

/**
 * DCDM_DR_Species: composite for Decaying Cold Dark Matter + its decay radiation.
 *
 * EnergyType::Other so the background loop splits rho correctly:
 *   rho_m += rho - 3p  (= rho_dcdm)
 *   rho_r += 3p        (= rho_dr)
 *
 * BackgroundDerivs override adds the DCDM->DR decay source on top of children.
 * AddCouplingDerivs adds the DR Boltzmann l=0,1 coupling source terms from DCDM.
 */
class DCDM_DR_Species : public CompositeSpecies {
 public:
  DCDM_DR_Species(std::shared_ptr<DarkRadiation> dr,
                  const background* pba,
                  const BackgroundModule* bgm);

  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;

  // Override to add DCDM->DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;

  // Typed accessors so callers can capture child indices
  DCDMSpecies& dcdm() {
    return *dcdm_;
  }
  DarkRadiationSpecies& dr() {
    return *dr_sp_;
  }
  const DCDMSpecies& dcdm() const {
    return *dcdm_;
  }
  const DarkRadiationSpecies& dr() const {
    return *dr_sp_;
  }

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) override;

 private:
  DCDMSpecies* dcdm_           = nullptr;  // non-owning pointer into children_
  DarkRadiationSpecies* dr_sp_ = nullptr;  // non-owning pointer into children_
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;
};
