#include "cdm.h"

#include "background.h"
#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

CDMSpecies::CDMSpecies(const background& pba)
    : BaseSpecies("CDM", EnergyType::Matter), Omega0_cdm_(pba.Omega0_cdm), H0_(pba.H0) {}

void CDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_cdm_, _TRUE_, index_bg, 1);
  index_bg_rho_ = index_bg_rho_cdm_;
}

void CDMSpecies::ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) {
  pvecback[index_bg_rho_cdm_] = Omega0_cdm_ * H0_ * H0_ / (a_rel * a_rel * a_rel);
}

double CDMSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_cdm_];
}

double CDMSpecies::P(const double* /*pvecback*/) const {
  return 0.;
}

double CDMSpecies::DpDloga(const double* /*pvecback*/) const {
  return 0.;
}

void CDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_cdm", 0.);
}

void CDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_cdm", pvecback[index_bg_rho_cdm_]);
}

void CDMSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                             const precision* /*ppr*/,
                                             int& index_pt,
                                             const perturb_workspace* /*ppw*/,
                                             int gauge) {
  class_define_index(pv->index_pt_delta_cdm, _TRUE_, index_pt, 1);

  /* theta_cdm is a dynamical variable only in Newtonian gauge (gauge==0);
     in synchronous gauge it is zero by definition. Sentinel -1 signals absent. */
  pv->index_pt_theta_cdm = -1;
  class_define_index(pv->index_pt_theta_cdm, (gauge == 0), index_pt, 1);
}

void CDMSpecies::PerturbDerivs(double /*tau*/,
                               const double* y,
                               double* dy,
                               const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const int gauge                 = (int) ppaw.perturbations_module->GetPerturbs()->gauge;

  if (gauge == 0) { /* newtonian */
    dy[pv->index_pt_delta_cdm] = -(y[pv->index_pt_theta_cdm] + ctx.metric_continuity);
    dy[pv->index_pt_theta_cdm] = -ctx.a_prime_over_a * y[pv->index_pt_theta_cdm] + ctx.metric_euler;
  }
  else { /* synchronous: theta_cdm = 0 by gauge choice */
    dy[pv->index_pt_delta_cdm] = -ctx.metric_continuity;
  }
}

double CDMSpecies::Delta(const perturb_vector* pv,
                         const double* y,
                         const double* /*pvecback*/,
                         const perturb_workspace* /*ppw*/) const {
  return y[pv->index_pt_delta_cdm];
}

double CDMSpecies::Theta(const perturb_vector* pv,
                         const double* y,
                         const double* /*pvecback*/,
                         const perturb_workspace* /*ppw*/) const {
  return (pv->index_pt_theta_cdm >= 0) ? y[pv->index_pt_theta_cdm] : 0.;
}

double CDMSpecies::DeltaP(const perturb_vector* /*pv*/,
                          const double* /*y*/,
                          const double* /*pvecback*/,
                          const perturb_workspace* /*ppw*/) const {
  return 0.;
}

double CDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/,
                                 const double* /*y*/,
                                 const double* /*pvecback*/,
                                 const perturb_workspace* /*ppw*/) const {
  return 0.;
}

void CDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv             = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  if (pv->index_pt_delta_cdm < 0)
    return;

  if (ctx.index_ic == mod->index_ic_ad_ || ctx.index_ic == mod->index_ic_bi_) {
    y[pv->index_pt_delta_cdm] = 3. / 4. * ctx.delta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[pv->index_pt_delta_cdm] = ctx.ppr->entropy_ini + 3. / 4. * ctx.delta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_nid_) {
    y[pv->index_pt_delta_cdm] = -ctx.ppr->entropy_ini * ctx.fracnu * ctx.fracb / ctx.fracg / 80. *
                                ctx.ktau_two * ctx.om * ctx.tau;
  }
  else if (ctx.index_ic == mod->index_ic_niv_) {
    y[pv->index_pt_delta_cdm] = -ctx.ppr->entropy_ini * 9. / 64. * ctx.fracnu * ctx.fracb /
                                ctx.fracg * ctx.k * ctx.tau * ctx.om * ctx.tau;
  }
}

void CDMSpecies::FillSources(const double* y, const double* /*dy*/, PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const perturb_vector* pv   = ppw->pv;

  const double a_prime_over_a = ctx.a_prime_over_a;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_cdm ──────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_cdm_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_cdm_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[pv->index_pt_delta_cdm] +
                              3. * a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_cdm ──────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_cdm_ == _TRUE_) {
    const double theta_cdm = (pv->index_pt_theta_cdm >= 0) ? y[pv->index_pt_theta_cdm] : 0.;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_cdm_,
                          ctx.index_tau,
                          ctx.index_k,
                          theta_cdm + ctx.theta_shift);  // N-body gauge correction
  }
}

void CDMSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                    const PerturbationsModule& mod,
                                    enum file_format fmt,
                                    BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_cdm", mod.index_tp_delta_cdm_, pba->has_cdm);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_cdm",
            mod.index_tp_theta_cdm_,
            (pba->has_cdm == _TRUE_) && (ppt->gauge != synchronous));
  }
  else if (fmt == camb_format) {
    if (section != TransferColumnSection::velocity)
      w.Add("-T_cdm/k2", mod.index_tp_delta_cdm_, pba->has_cdm);
  }
}

void CDMSpecies::PrintVariables(PerturbColumnWriter& w,
                                double /*tau*/,
                                const double* y,
                                const PerturbationsModule& mod,
                                const perturb_workspace* ppw) const {
  const background* pba = mod.GetBackground();
  if (pba->has_cdm != _TRUE_)
    return;

  double delta_cdm = 0., theta_cdm = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double k           = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;

    delta_cdm           = y[pv->index_pt_delta_cdm];
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      theta_cdm = 0.;
    }
    else {
      theta_cdm = y[pv->index_pt_theta_cdm];
    }

    // Converting synchronous variables to Newtonian
    if (ppt->gauge == synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_cdm          -= 3. * H * a * alpha;
      theta_cdm          += k * k * alpha;
    }
  }

  w.Add("delta_cdm", delta_cdm, true);
  w.Add("theta_cdm", theta_cdm, true);
}
