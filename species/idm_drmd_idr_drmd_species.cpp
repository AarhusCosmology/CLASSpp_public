#include "idm_drmd_idr_drmd_species.h"
#include "background_module.h"
#include "perturbations_module.h"

IDM_DRMD_IDR_DRMD_Species::IDM_DRMD_IDR_DRMD_Species(const background& pba)
  : CompositeSpecies("IDM_DRMD_IDR_DRMD", BaseSpecies::EnergyType::Other)
  , pba_(pba)
{
  auto idm = std::make_unique<IDM_DRMDSpecies>(pba);
  auto idr = std::make_unique<IDR_DRMDSpecies>(pba);
  idm_drmd_ = idm.get();
  idr_drmd_ = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}

void IDM_DRMD_IDR_DRMD_Species::AddCouplingDerivs(
    double /*tau*/, const double* y, double* dy,
    const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace*    ppw = ppaw.ppw;
  const perturb_vector*       pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  auto* bgm = ppaw.perturbations_module->GetBackgroundModule().get();
  const double* pvecback = ppw->pvecback;

  const double rho_idm_drmd = idm_drmd_->Rho(pvecback);
  const double rho_idr_drmd = idr_drmd_->Rho(pvecback);

  // Guard against zero densities to avoid division by zero in background_idm_drmd.
  if (rho_idm_drmd <= 0. || rho_idr_drmd <= 0.) return;

  double Rint, csp2, Gint;
  class_call(bgm->background_idm_drmd(
      ctx.a,
      rho_idm_drmd / rho_idr_drmd,
      &Rint, &csp2, &Gint),
    bgm->error_message_, bgm->error_message_);

  const double theta_idm = y[pv->index_pt_theta_idm_drmd];
  const double theta_idr = y[pv->index_pt_theta_idr_drmd];
  const double delta_idr = y[pv->index_pt_delta_idr_drmd];

  if (ppw->approx[ppw->index_ap_tca_idm_drmd] == (int)tca_idm_drmd_on) {
    // TCA on: ASSIGN theta_idm_drmd (replaces the free-streaming term written by child)
    double GdDelta = 3.*csp2*(ctx.a_prime_over_a*theta_idm + ctx.k2*delta_idr/4.);
    dy[pv->index_pt_theta_idm_drmd] = 0.25*ctx.k2*delta_idr + ctx.metric_euler
                                     - GdDelta * Rint;
    // IDR_DRMD: subtract coupling from free-streaming equation already written
    dy[pv->index_pt_theta_idr_drmd] -= GdDelta * Rint;
  } else {
    // TCA off: add coupling increments
    dy[pv->index_pt_theta_idm_drmd] += Gint*(theta_idr - theta_idm);
    dy[pv->index_pt_theta_idr_drmd] -= Gint * Rint * (theta_idr - theta_idm);
  }
}
