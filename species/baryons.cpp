#include "baryons.h"

#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

void BaryonsSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                                 const precision* /*ppr*/,
                                                 int& index_pt,
                                                 const perturb_workspace* /*ppw*/,
                                                 int /*gauge*/) {
  class_define_index(pv->index_pt_delta_b, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_b, _TRUE_, index_pt, 1);
}

void BaryonsSpecies::PerturbDerivs(double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double k2                   = ppaw.k * ppaw.k;
  const double a_prime_over_a       = ctx.a_prime_over_a;
  const double metric_continuity    = ctx.metric_continuity;
  const double metric_euler         = ctx.metric_euler;
  const double R                    = ctx.R;
  const double s2_squared           = ctx.s2_squared;
  const double delta_p_b_over_rho_b = ctx.delta_p_b_over_rho_b;

  const double theta_b = y[pv->index_pt_theta_b];
  /* theta_g from context: RSA-corrected photon velocity */
  const double theta_g = ctx.theta_g;
  const double delta_g = ctx.delta_g;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  /* density equation */
  dy[pv->index_pt_delta_b] = -(theta_b + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == (int) tca_off) {
    /* Full equation */
    dy[pv->index_pt_theta_b] = -a_prime_over_a * theta_b + metric_euler +
                               k2 * delta_p_b_over_rho_b + R * dkappa * (theta_g - theta_b);
  }
  else {
    /* TCA on: tight-coupling approximation for theta_b */
    dy[pv->index_pt_theta_b] =
        (-a_prime_over_a * theta_b +
         k2 * (delta_p_b_over_rho_b + R * (delta_g / 4. - s2_squared * ppw->tca_shear_g)) +
         R * ppw->tca_slip) /
            (1. + R) +
        metric_euler;
  }
}

double BaryonsSpecies::DeltaP(const perturb_vector* /*pv*/,
                              const double* /*y*/,
                              const double* pvecback,
                              const perturb_workspace* ppw) const {
  return pvecback[index_bg_rho_] * ppw->scalar_ctx.delta_p_b_over_rho_b;
}

void BaryonsSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv             = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  if (pv->index_pt_delta_b < 0 || pv->index_pt_theta_b < 0)
    return;

  if (ctx.index_ic == mod->index_ic_ad_) {
    y[pv->index_pt_delta_b] = 3. / 4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[pv->index_pt_delta_b] = 3. / 4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_bi_) {
    y[pv->index_pt_delta_b] = ctx.ppr->entropy_ini + 3. / 4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_nid_) {
    y[pv->index_pt_delta_b] = ctx.ppr->entropy_ini * ctx.fracnu / ctx.fracg / 8. * ctx.ktau_two;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_niv_) {
    y[pv->index_pt_delta_b] = 3. / 4. * ctx.delta_g_ic;
    y[pv->index_pt_theta_b] = ctx.theta_g_ic;
  }
}

void BaryonsSpecies::FillSources(const double* y, const double* /*dy*/, PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const perturb_vector* pv   = ppw->pv;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const double a_prime_over_a = ctx.a_prime_over_a;

  // ── delta_b: baryon density transfer ─────────────────────────────────────
  if (p_mod->has_source_delta_b_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_b_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[pv->index_pt_delta_b] +
                              3. * a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_b: baryon velocity transfer ────────────────────────────────────
  if (p_mod->has_source_theta_b_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_b_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[pv->index_pt_theta_b] + ctx.theta_shift);  // N-body gauge correction
  }
}

void BaryonsSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                        const PerturbationsModule& mod,
                                        enum file_format fmt,
                                        BaseSpecies::TransferColumnSection section) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_b", mod.index_tp_delta_b_, true);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_b", mod.index_tp_theta_b_, true);
  }
}

void BaryonsSpecies::PrintVariables(PerturbColumnWriter& w,
                                    double /*tau*/,
                                    const double* y,
                                    const PerturbationsModule& mod,
                                    const perturb_workspace* ppw) const {
  double delta_b = 0., theta_b = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double k           = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;

    delta_b = y[pv->index_pt_delta_b];
    theta_b = y[pv->index_pt_theta_b];

    // Converting synchronous variables to Newtonian
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_b            -= 3. * H * a * alpha;
      theta_b            += k * k * alpha;
    }
  }

  w.Add("delta_b", delta_b, true);
  w.Add("theta_b", theta_b, true);
}
