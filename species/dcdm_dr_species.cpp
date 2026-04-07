#include "dcdm_dr_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include <cmath>

DCDM_DR_Species::DCDM_DR_Species(std::shared_ptr<DarkRadiation> dr,
                                   const background* pba,
                                   const BackgroundModule* bgm)
  : CompositeSpecies("DCDM_DR", BaseSpecies::EnergyType::Other)
  , pba_(pba), bgm_(bgm)
{
  auto dcdm = std::make_unique<DCDMSpecies>(*pba);
  auto dr_sp = std::make_unique<DarkRadiationSpecies>(dr, pba, bgm);
  dcdm_  = dcdm.get();
  dr_sp_ = dr_sp.get();
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(dr_sp));
}

void DCDM_DR_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

void DCDM_DR_Species::BackgroundDerivs(double tau, const double* y,
                                        double* dy, const double* pvecback) {
  // Children handle their own dilution terms
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);

  // DCDM->DR decay source for first DR channel
  const double a = pvecback[bgm_->index_bg_a_];
  dy[dr_sp_->bi_rho_dr_species_index()] +=
      a * pba_->Gamma_dcdm * pvecback[dcdm_->bg_rho_index()];
}

void DCDM_DR_Species::AddCouplingDerivs(double /*tau*/, const double* y, double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*    ppw = ppaw.ppw;
  const perturb_vector*       pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  if (pv->index_pt_delta_dcdm < 0 || pv->index_pt_F0_dr_species < 0) return;

  const double* pvecback = ppw->pvecback;
  const double a = ctx.a;
  const double k = ctx.k;

  const int base = pv->index_pt_F0_dr_species;  // first DR channel, index_dr=0
  const double rprime_dr = pba_->Gamma_dcdm
                           * pvecback[dcdm_->bg_rho_index()]
                           * std::pow(a, 5) / (pba_->H0 * pba_->H0);

  const double delta_dcdm = y[pv->index_pt_delta_dcdm];
  const double theta_dcdm = y[pv->index_pt_theta_dcdm];

  // Add DCDM source to DR l=0 and l=1
  const double dl0 = rprime_dr * (delta_dcdm + ctx.metric_euler / (k * k));
  const double dl1 = rprime_dr / k * theta_dcdm;
  dy[base + 0] += dl0;
  dy[base + 1] += dl1;

  // Keep sum slots consistent
  dy[pv->index_pt_F0_dr_sum + 0] += dl0;
  dy[pv->index_pt_F0_dr_sum + 1] += dl1;
}
