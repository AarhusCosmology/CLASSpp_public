#include "dcdm.h"

#include "background.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

DCDMSpecies::DCDMSpecies(const background& pba)
    : BaseSpecies("DCDM", EnergyType::Matter), pba_(pba) {}

void DCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_dcdm_, _TRUE_, index_bg, 1);
  index_bg_rho_ = index_bg_rho_dcdm_;
}

void DCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_rho_dcdm_, _TRUE_, index_bi, 1);
}

void DCDMSpecies::SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {
  pvecback_integration[index_bi_rho_dcdm_] = pba_.Omega_ini_dcdm * std::pow(pba_.H0, 2) *
                                             std::pow(1.0 / a_rel, 3);
}

void DCDMSpecies::ComputeBackground(double /*a_rel*/, const double* pvecback_B, double* pvecback) {
  pvecback[index_bg_rho_dcdm_] = pvecback_B[index_bi_rho_dcdm_];
}

void DCDMSpecies::BackgroundDerivs(double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const double* pvecback) {
  const double a   = pvecback[bgm_->index_bg_a_];
  const double H   = pvecback[bgm_->index_bg_H_];
  const double rho = y[index_bi_rho_dcdm_];
  /** rho' = -a*(3H + Gamma) * rho */
  dy[index_bi_rho_dcdm_] = -a * (3. * H + pba_.Gamma_dcdm) * rho;
}

double DCDMSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_dcdm_];
}

double DCDMSpecies::P(const double* /*pvecback*/) const {
  return 0.;
}
double DCDMSpecies::DpDloga(const double* /*pvecback*/) const {
  return 0.;
}

void DCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                              const precision* /*ppr*/,
                                              int& index_pt,
                                              const perturb_workspace* /*ppw*/,
                                              int /*gauge*/) {
  class_define_index(pv->index_pt_delta_dcdm, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_dcdm, _TRUE_, index_pt, 1);
}

void DCDMSpecies::PerturbDerivs(double /*tau*/,
                                const double* y,
                                double* dy,
                                const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[pv->index_pt_delta_dcdm] = -(y[pv->index_pt_theta_dcdm] + ctx.metric_continuity) -
                                ctx.a * pba_.Gamma_dcdm / ctx.k2 * ctx.metric_euler;
  dy[pv->index_pt_theta_dcdm] = -ctx.a_prime_over_a * y[pv->index_pt_theta_dcdm] + ctx.metric_euler;
}

double DCDMSpecies::Delta(const perturb_vector* pv,
                          const double* y,
                          const double* /*pvecback*/,
                          const perturb_workspace* /*ppw*/) const {
  return y[pv->index_pt_delta_dcdm];
}
double DCDMSpecies::Theta(const perturb_vector* pv,
                          const double* y,
                          const double* /*pvecback*/,
                          const perturb_workspace* /*ppw*/) const {
  return y[pv->index_pt_theta_dcdm];
}
double DCDMSpecies::DeltaP(const perturb_vector* /*pv*/,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  return 0.;
}
double DCDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/,
                                  const double* /*y*/,
                                  const double* /*pvecback*/,
                                  const perturb_workspace* /*ppw*/) const {
  return 0.;
}

void DCDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  perturb_vector* pv = ctx.ppw->pv;
  if (pv->index_pt_delta_dcdm >= 0)
    y[pv->index_pt_delta_dcdm] = 3. / 4. * ctx.delta_g_ic;
}

void DCDMSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                     const PerturbationsModule& mod,
                                     enum file_format fmt,
                                     BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_dcdm", mod.index_tp_delta_dcdm_, pba->has_dcdm);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_dcdm", mod.index_tp_theta_dcdm_, pba->has_dcdm);
  }
}

void DCDMSpecies::PrintVariables(PerturbColumnWriter& w,
                                 double /*tau*/,
                                 const double* y,
                                 const PerturbationsModule& mod,
                                 const perturb_workspace* ppw) const {
  const background* pba = mod.GetBackground();
  if (pba->has_dcdm != _TRUE_)
    return;

  double delta_dcdm = 0., theta_dcdm = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double k           = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;

    delta_dcdm = y[pv->index_pt_delta_dcdm];
    theta_dcdm = y[pv->index_pt_theta_dcdm];

    // Synchronous → Newtonian gauge conversion for DCDM (matter + decay term)
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_dcdm         += alpha * (-a * pba_.Gamma_dcdm - 3. * a * H);
      theta_dcdm         += k * k * alpha;
    }
  }

  w.Add("delta_dcdm", delta_dcdm, true);
  w.Add("theta_dcdm", theta_dcdm, true);
}
