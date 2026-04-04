#include "fluid.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "background_module.h"
#include "background.h"
#include <cmath>

FluidSpecies::FluidSpecies(const background& pba)
    : BaseSpecies("Fluid", EnergyType::DarkEnergy), pba_(pba) {}

void FluidSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_fld_,        _TRUE_, index_bg, 1);
  class_define_index(index_bg_w_fld_,          _TRUE_, index_bg, 1);
  class_define_index(index_bg_dw_over_da_fld_, _TRUE_, index_bg, 1);
}

void FluidSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_rho_fld_, _TRUE_, index_bi, 1);
}

void FluidSpecies::ComputeBackground(double /*a_rel*/, const double* pvecback_B,
                                      double* pvecback) {
  /* w_fld and dw_over_da_fld are pre-computed by BackgroundModule::background_functions()
     (with checked error handling) and written to pvecback before this call. */
  pvecback[index_bg_rho_fld_] = pvecback_B[index_bi_rho_fld_];
}

void FluidSpecies::BackgroundDerivs(double /*tau*/, const double* y, double* dy,
                                     const double* pvecback) {
  const double a     = pvecback[bgm_->index_bg_a_];
  const double H     = pvecback[bgm_->index_bg_H_];
  const double w_fld = pvecback[index_bg_w_fld_];
  /** rho' = -3*a*H*(1+w)*rho */
  dy[index_bi_rho_fld_] = -3.*a*H*(1. + w_fld)*y[index_bi_rho_fld_];
}

double FluidSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_fld_];
}

double FluidSpecies::P(const double* pvecback) const {
  return pvecback[index_bg_w_fld_] * pvecback[index_bg_rho_fld_];
}

double FluidSpecies::DpDloga(const double* pvecback) const {
  const double w_fld = pvecback[index_bg_w_fld_];
  const double dw_over_da_fld = pvecback[index_bg_dw_over_da_fld_];
  const double a = pvecback[bgm_->index_bg_a_];
  return (a * dw_over_da_fld - 3. * (1. + w_fld) * w_fld) * pvecback[index_bg_rho_fld_];
}

void FluidSpecies::RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                                const perturb_workspace* /*ppw*/,
                                                int /*gauge*/) {
  if (pba_.use_ppf == _FALSE_) {
    pv->index_pt_delta_fld = -1;
    pv->index_pt_theta_fld = -1;
    class_define_index(pv->index_pt_delta_fld, _TRUE_, index_pt, 1);
    class_define_index(pv->index_pt_theta_fld, _TRUE_, index_pt, 1);
  } else {
    pv->index_pt_delta_fld = -1;
    pv->index_pt_theta_fld = -1;
    class_define_index(pv->index_pt_Gamma_fld, _TRUE_, index_pt, 1);
  }
}

void FluidSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                  const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector*    pv     = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  const double a                = ctx.a;
  const double a_prime_over_a   = ctx.a_prime_over_a;
  const double k2               = ctx.k2;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;

  if (pba_.use_ppf == _FALSE_) {
    const double w_fld = ppw->pvecback[index_bg_w_fld_];
    const double dw_over_da_fld = ppw->pvecback[index_bg_dw_over_da_fld_];
    const double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;
    const double ca2 = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;
    const double cs2 = pba_.cs2_fld;

    dy[pv->index_pt_delta_fld] =
        -(1. + w_fld) * (y[pv->index_pt_theta_fld] + metric_continuity)
        - 3. * (cs2 - w_fld) * a_prime_over_a * y[pv->index_pt_delta_fld]
        - 9. * (1. + w_fld) * (cs2 - ca2) * a_prime_over_a * a_prime_over_a
              * y[pv->index_pt_theta_fld] / k2;

    dy[pv->index_pt_theta_fld] =
        -(1. - 3. * cs2) * a_prime_over_a * y[pv->index_pt_theta_fld]
        + cs2 * k2 / (1. + w_fld) * y[pv->index_pt_delta_fld]
        + metric_euler;
  } else {
    dy[pv->index_pt_Gamma_fld] = ppw->Gamma_prime_fld;
  }
}

double FluidSpecies::Delta(const perturb_vector* pv, const double* y, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return (pv->index_pt_delta_fld >= 0) ? y[pv->index_pt_delta_fld] : 0.;
}
double FluidSpecies::Theta(const perturb_vector* pv, const double* y, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return (pv->index_pt_theta_fld >= 0) ? y[pv->index_pt_theta_fld] : 0.;
}
double FluidSpecies::DeltaP(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const { return 0.; }
double FluidSpecies::RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const { return 0.; }
