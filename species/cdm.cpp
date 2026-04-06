#include "cdm.h"
#include "background.h"
#include "perturbations.h"
#include "perturbations_module.h"

CDMSpecies::CDMSpecies(const background& pba)
    : BaseSpecies("CDM", EnergyType::Matter),
      Omega0_cdm_(pba.Omega0_cdm), H0_(pba.H0) {}

void CDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_cdm_, _TRUE_, index_bg, 1);
}

void CDMSpecies::ComputeBackground(double a_rel, const double* /*pvecback_B*/,
                                    double* pvecback) {
  pvecback[index_bg_rho_cdm_] = Omega0_cdm_ * H0_ * H0_ / (a_rel * a_rel * a_rel);
}

double CDMSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_cdm_];
}

double CDMSpecies::P(const double* /*pvecback*/) const { return 0.; }

double CDMSpecies::DpDloga(const double* /*pvecback*/) const { return 0.; }

void CDMSpecies::RegisterPerturbationIndices(perturb_vector* pv, const precision* /*ppr*/,
                                              int& index_pt,
                                              const perturb_workspace* /*ppw*/,
                                              int gauge) {
  class_define_index(pv->index_pt_delta_cdm, _TRUE_, index_pt, 1);

  /* theta_cdm is a dynamical variable only in Newtonian gauge (gauge==0);
     in synchronous gauge it is zero by definition. Sentinel -1 signals absent. */
  pv->index_pt_theta_cdm = -1;
  class_define_index(pv->index_pt_theta_cdm, (gauge == 0), index_pt, 1);
}

void CDMSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;
  const int gauge = (int)ppaw.perturbations_module->GetPerturbs()->gauge;

  if (gauge == 0) { /* newtonian */
    dy[pv->index_pt_delta_cdm] =
        -(y[pv->index_pt_theta_cdm] + ctx.metric_continuity);
    dy[pv->index_pt_theta_cdm] =
        -ctx.a_prime_over_a * y[pv->index_pt_theta_cdm] + ctx.metric_euler;
  } else { /* synchronous: theta_cdm = 0 by gauge choice */
    dy[pv->index_pt_delta_cdm] = -ctx.metric_continuity;
  }
}

double CDMSpecies::Delta(const perturb_vector* pv, const double* y, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return y[pv->index_pt_delta_cdm];
}

double CDMSpecies::Theta(const perturb_vector* pv, const double* y, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return (pv->index_pt_theta_cdm >= 0) ? y[pv->index_pt_theta_cdm] : 0.;
}

double CDMSpecies::DeltaP(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return 0.;
}

double CDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return 0.;
}
