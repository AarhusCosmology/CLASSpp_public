#include "fluid.h"

#include <cmath>

#include "background.h"
#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

FluidSpecies::FluidSpecies(const background& pba)
    : BaseSpecies("Fluid", EnergyType::DarkEnergy), pba_(pba) {}

void FluidSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_fld_, _TRUE_, index_bg, 1);
  index_bg_rho_ = index_bg_rho_fld_;
  class_define_index(index_bg_w_fld_, _TRUE_, index_bg, 1);
  class_define_index(index_bg_dw_over_da_fld_, _TRUE_, index_bg, 1);
}

void FluidSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_rho_fld_, _TRUE_, index_bi, 1);
}

void FluidSpecies::ComputeBackground(double /*a_rel*/, const double* pvecback_B, double* pvecback) {
  /* w_fld and dw_over_da_fld are pre-computed by BackgroundModule::background_functions()
     (with checked error handling) and written to pvecback before this call. */
  pvecback[index_bg_rho_fld_] = pvecback_B[index_bi_rho_fld_];
}

void FluidSpecies::BackgroundDerivs(double /*tau*/,
                                    const double* y,
                                    double* dy,
                                    const double* pvecback) {
  const double a     = pvecback[bgm_->index_bg_a_];
  const double H     = pvecback[bgm_->index_bg_H_];
  const double w_fld = pvecback[index_bg_w_fld_];
  /** rho' = -3*a*H*(1+w)*rho */
  dy[index_bi_rho_fld_] = -3. * a * H * (1. + w_fld) * y[index_bi_rho_fld_];
}

double FluidSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_fld_];
}

double FluidSpecies::P(const double* pvecback) const {
  return pvecback[index_bg_w_fld_] * pvecback[index_bg_rho_fld_];
}

double FluidSpecies::DpDloga(const double* pvecback) const {
  const double w_fld          = pvecback[index_bg_w_fld_];
  const double dw_over_da_fld = pvecback[index_bg_dw_over_da_fld_];
  const double a              = pvecback[bgm_->index_bg_a_];
  return (a * dw_over_da_fld - 3. * (1. + w_fld) * w_fld) * pvecback[index_bg_rho_fld_];
}

void FluidSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_fld", 0.);
  w.Add("(.)w_fld", 0.);
}

void FluidSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_fld", pvecback[index_bg_rho_fld_]);
  w.Add("(.)w_fld", pvecback[index_bg_w_fld_]);
}

void FluidSpecies::WriteWFld(double w_fld, double dw_over_da_fld, double* pvecback) const {
  pvecback[index_bg_w_fld_]          = w_fld;
  pvecback[index_bg_dw_over_da_fld_] = dw_over_da_fld;
}

void FluidSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                               const precision* /*ppr*/,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  if (pba_.use_ppf == _FALSE_) {
    pv->index_pt_delta_fld = -1;
    pv->index_pt_theta_fld = -1;
    class_define_index(pv->index_pt_delta_fld, _TRUE_, index_pt, 1);
    class_define_index(pv->index_pt_theta_fld, _TRUE_, index_pt, 1);
  }
  else {
    pv->index_pt_delta_fld = -1;
    pv->index_pt_theta_fld = -1;
    class_define_index(pv->index_pt_Gamma_fld, _TRUE_, index_pt, 1);
  }
}

void FluidSpecies::PerturbDerivs(double /*tau*/,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double a                 = ctx.a;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double k2                = ctx.k2;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;

  if (pba_.use_ppf == _FALSE_) {
    const double w_fld          = ppw->pvecback[index_bg_w_fld_];
    const double dw_over_da_fld = ppw->pvecback[index_bg_dw_over_da_fld_];
    const double w_prime_fld    = dw_over_da_fld * a_prime_over_a * a;
    const double ca2            = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;
    const double cs2            = pba_.cs2_fld;

    dy[pv->index_pt_delta_fld] = -(1. + w_fld) * (y[pv->index_pt_theta_fld] + metric_continuity) -
                                 3. * (cs2 - w_fld) * a_prime_over_a * y[pv->index_pt_delta_fld] -
                                 9. * (1. + w_fld) * (cs2 - ca2) * a_prime_over_a * a_prime_over_a *
                                     y[pv->index_pt_theta_fld] / k2;

    dy[pv->index_pt_theta_fld] = -(1. - 3. * cs2) * a_prime_over_a * y[pv->index_pt_theta_fld] +
                                 cs2 * k2 / (1. + w_fld) * y[pv->index_pt_delta_fld] + metric_euler;
  }
  else {
    dy[pv->index_pt_Gamma_fld] = ppw->Gamma_prime_fld;
  }
}

void FluidSpecies::FillSources(const double* /*y*/,
                               const double* /*dy*/,
                               PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod  = ctx.p_mod;
  perturb_workspace* ppw      = ctx.ppw;
  const BackgroundModule* bgm = p_mod->GetBackgroundModule().get();
  const double* pvecback      = ppw->pvecback;

  // Fluid sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_fld ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_fld_ == _TRUE_) {
    const double w_fld = W(pvecback);
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_fld_,
                          ctx.index_tau,
                          ctx.index_k,
                          ppw->delta_rho_fld / Rho(pvecback) +
                              3. * ctx.a_prime_over_a * (1. + w_fld) *
                                  ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_fld ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_fld_ == _TRUE_) {
    const double w_fld = W(pvecback);
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_fld_,
                          ctx.index_tau,
                          ctx.index_k,
                          ppw->rho_plus_p_theta_fld / (1. + w_fld) / Rho(pvecback) +
                              ctx.theta_shift);  // N-body gauge correction
  }
}

void FluidSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  if (ctx.p_mod->GetBackground()->use_ppf == _TRUE_)
    return;
  perturb_vector* pv = ctx.ppw->pv;
  if (pv->index_pt_delta_fld < 0 || pv->index_pt_theta_fld < 0)
    return;

  double w_fld, dw_over_da, integral;
  auto* bgm = ctx.p_mod->GetBackgroundModule().get();
  class_call(bgm->background_w_fld(ctx.a, &w_fld, &dw_over_da, &integral),
             bgm->error_message_,
             bgm->error_message_);

  y[pv->index_pt_delta_fld] = -ctx.ktau_two / 4. * (1. + w_fld) *
                              (4. - 3. * ctx.p_mod->GetBackground()->cs2_fld) /
                              (4. - 6. * w_fld + 3. * ctx.p_mod->GetBackground()->cs2_fld) *
                              ctx.ppr->curvature_ini * ctx.s2_squared;

  y[pv->index_pt_theta_fld] = -ctx.k * ctx.ktau_three / 4. * ctx.p_mod->GetBackground()->cs2_fld /
                              (4. - 6. * w_fld + 3. * ctx.p_mod->GetBackground()->cs2_fld) *
                              ctx.ppr->curvature_ini * ctx.s2_squared;
}

void FluidSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                      const PerturbationsModule& mod,
                                      enum file_format fmt,
                                      BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_fld", mod.index_tp_delta_fld_, pba->has_fld);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_fld", mod.index_tp_theta_fld_, pba->has_fld);
  }
}

void FluidSpecies::PrintVariables(PerturbColumnWriter& w,
                                  double /*tau*/,
                                  const double* /*y*/,
                                  const PerturbationsModule& mod,
                                  const perturb_workspace* ppw) const {
  const background* pba = mod.GetBackground();
  if (pba->has_fld != _TRUE_)
    return;

  double delta_rho_fld = 0., rho_plus_p_theta_fld = 0., delta_p_fld = 0.;

  if (!w.IsTitleMode()) {
    delta_rho_fld        = ppw->delta_rho_fld;
    rho_plus_p_theta_fld = ppw->rho_plus_p_theta_fld;
    delta_p_fld          = ppw->delta_p_fld;
  }

  w.Add("delta_rho_fld", delta_rho_fld, true);
  w.Add("rho_plus_p_theta_fld", rho_plus_p_theta_fld, true);
  w.Add("delta_p_fld", delta_p_fld, true);
}

double FluidSpecies::Delta(const perturb_vector* pv,
                           const double* y,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  return (pv->index_pt_delta_fld >= 0) ? y[pv->index_pt_delta_fld] : 0.;
}
double FluidSpecies::Theta(const perturb_vector* pv,
                           const double* y,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  return (pv->index_pt_theta_fld >= 0) ? y[pv->index_pt_theta_fld] : 0.;
}
double FluidSpecies::DeltaP(const perturb_vector* /*pv*/,
                            const double* /*y*/,
                            const double* /*pvecback*/,
                            const perturb_workspace* /*ppw*/) const {
  return 0.;
}
double FluidSpecies::RhoPlusPShear(const perturb_vector* /*pv*/,
                                   const double* /*y*/,
                                   const double* /*pvecback*/,
                                   const perturb_workspace* /*ppw*/) const {
  return 0.;
}
