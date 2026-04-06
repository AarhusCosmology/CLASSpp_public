#include "dcdm.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "background.h"
#include "background_module.h"

DCDMSpecies::DCDMSpecies(const background& pba)
    : BaseSpecies("DCDM", EnergyType::Matter), pba_(pba) {}

void DCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_dcdm_, _TRUE_, index_bg, 1);
}

void DCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_rho_dcdm_, _TRUE_, index_bi, 1);
}

void DCDMSpecies::ComputeBackground(double /*a_rel*/, const double* pvecback_B,
                                     double* pvecback) {
  pvecback[index_bg_rho_dcdm_] = pvecback_B[index_bi_rho_dcdm_];
}

void DCDMSpecies::BackgroundDerivs(double /*tau*/, const double* y, double* dy,
                                    const double* pvecback) {
  const double a   = pvecback[bgm_->index_bg_a_];
  const double H   = pvecback[bgm_->index_bg_H_];
  const double rho = y[index_bi_rho_dcdm_];
  /** rho' = -a*(3H + Gamma) * rho */
  dy[index_bi_rho_dcdm_] = -a * (3.*H + pba_.Gamma_dcdm) * rho;
}

double DCDMSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_dcdm_];
}

double DCDMSpecies::P(const double* /*pvecback*/) const { return 0.; }
double DCDMSpecies::DpDloga(const double* /*pvecback*/) const { return 0.; }

void DCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv, const precision* /*ppr*/,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  class_define_index(pv->index_pt_delta_dcdm, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_dcdm, _TRUE_, index_pt, 1);
}

void DCDMSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw   = ppaw.ppw;
  const perturb_vector* pv       = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[pv->index_pt_delta_dcdm] = -(y[pv->index_pt_theta_dcdm] + ctx.metric_continuity)
      - ctx.a * pba_.Gamma_dcdm / ctx.k2 * ctx.metric_euler;
  dy[pv->index_pt_theta_dcdm] = -ctx.a_prime_over_a * y[pv->index_pt_theta_dcdm]
      + ctx.metric_euler;
}

double DCDMSpecies::Delta(const perturb_vector* pv, const double* y, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return y[pv->index_pt_delta_dcdm];
}
double DCDMSpecies::Theta(const perturb_vector* pv, const double* y, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const {
  return y[pv->index_pt_theta_dcdm];
}
double DCDMSpecies::DeltaP(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const { return 0.; }
double DCDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const { return 0.; }
