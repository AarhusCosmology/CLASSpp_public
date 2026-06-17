#include "cdm.h"

#include "background.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

CDMSpecies::CDMSpecies(const background& pba, double omega0_cdm)
    : BaseSpecies("CDM", EnergyType::Matter), Omega0_cdm_(omega0_cdm), H0_(pba.H0) {}

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

void CDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_cdm", 0.);
}

void CDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_cdm", pvecback[index_bg_rho_cdm_]);
}

void CDMSpecies::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity && ctx.gauge != synchronous, index_tp, 1);
}

void CDMSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                             perturb_vector* pv,
                                             const precision* /*ppr*/,
                                             int& index_pt,
                                             const perturb_workspace* /*ppw*/,
                                             int gauge) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = index_pt;
  ++index_pt;

  /* theta_cdm is a dynamical variable only in Newtonian gauge (gauge==0);
     in synchronous gauge it is zero by definition. Sentinel -1 signals absent. */
  layout.idx_theta = -1;
  if (gauge == 0) {
    layout.idx_theta = index_pt;
    ++index_pt;
  }
}

void CDMSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                               double /*tau*/,
                               const double* y,
                               double* dy,
                               const perturb_parameters_and_workspace& ppaw) {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const int gauge                 = (int) ppaw.perturbations_module->GetPerturbs()->gauge;

  if (gauge == 0) { /* newtonian */
    dy[layout.idx_delta] = -(y[layout.idx_theta] + ctx.metric_continuity);
    dy[layout.idx_theta] = -ctx.a_prime_over_a * y[layout.idx_theta] + ctx.metric_euler;
  }
  else { /* synchronous: theta_cdm = 0 by gauge choice */
    dy[layout.idx_delta] = -ctx.metric_continuity;
  }
}

double CDMSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                         const perturb_vector* /*pv*/,
                         const double* y,
                         const double* /*pvecback*/,
                         const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  return y[layout.idx_delta];
}

double CDMSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                         const perturb_vector* /*pv*/,
                         const double* y,
                         const double* /*pvecback*/,
                         const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  return (layout.idx_theta >= 0) ? y[layout.idx_theta] : 0.;
}

double CDMSpecies::DeltaP(const BaseSpecies::PerturbLayout& /*base*/,
                          const perturb_vector* /*pv*/,
                          const double* /*y*/,
                          const double* /*pvecback*/,
                          const perturb_workspace* /*ppw*/) const {
  return 0.;
}

double CDMSpecies::RhoPlusPShear(const BaseSpecies::PerturbLayout& /*base*/,
                                 const perturb_vector* /*pv*/,
                                 const double* /*y*/,
                                 const double* /*pvecback*/,
                                 const perturb_workspace* /*ppw*/) const {
  return 0.;
}

void CDMSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                        double* y,
                                        const PerturbIcContext& ctx) {
  const auto& layout             = static_cast<const PerturbLayout&>(base);
  const PerturbationsModule* mod = ctx.p_mod;
  if (layout.idx_delta < 0)
    return;

  if (ctx.index_ic == mod->index_ic_ad_ || ctx.index_ic == mod->index_ic_bi_) {
    y[layout.idx_delta] = 3. / 4. * ctx.delta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[layout.idx_delta] = ctx.ppr->entropy_ini + 3. / 4. * ctx.delta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_nid_) {
    y[layout.idx_delta] = -ctx.ppr->entropy_ini * ctx.fracnu * ctx.fracb / ctx.fracg / 80. *
                          ctx.ktau_two * ctx.om * ctx.tau;
  }
  else if (ctx.index_ic == mod->index_ic_niv_) {
    y[layout.idx_delta] = -ctx.ppr->entropy_ini * 9. / 64. * ctx.fracnu * ctx.fracb / ctx.fracg *
                          ctx.k * ctx.tau * ctx.om * ctx.tau;
  }
}

void CDMSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                             const double* y,
                             const double* /*dy*/,
                             PerturbSourceContext& ctx) {
  const auto& layout         = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod = ctx.p_mod;

  const double a_prime_over_a = ctx.a_prime_over_a;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_cdm ──────────────────────────────────────────────────────────────
  if (index_tp_delta_ >= 0) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_delta_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[layout.idx_delta] +
                              3. * a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_cdm ──────────────────────────────────────────────────────────────
  if (index_tp_theta_ >= 0) {
    const double theta_cdm = (layout.idx_theta >= 0) ? y[layout.idx_theta] : 0.;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_theta_,
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
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
      w.Add("d_cdm", index_tp_delta_, index_tp_delta_ >= 0);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
      w.Add("t_cdm", index_tp_theta_, index_tp_theta_ >= 0);
  }
  else if (fmt == camb_format) {
    if (section != TransferColumnSection::velocity)
      w.Add("-T_cdm/k2", index_tp_delta_, index_tp_delta_ >= 0);
  }
}

void CDMSpecies::PrintVariables(PerturbColumnWriter& w,
                                double /*tau*/,
                                const double* y,
                                const PerturbationsModule& mod,
                                const perturb_workspace* ppw) const {
  double delta_cdm = 0., theta_cdm = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    const auto& layout = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
    const double k     = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;

    delta_cdm           = y[layout.idx_delta];
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      theta_cdm = 0.;
    }
    else {
      theta_cdm = (layout.idx_theta >= 0) ? y[layout.idx_theta] : 0.;
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

// ── Newtonian-gauge transform ─────────────────────────────────────────────────

void CDMSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                               double* y,
                                               const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}

void CDMSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                               const BaseSpecies::PerturbLayout& new_base,
                                               const double* old_y,
                                               double* new_y,
                                               const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  // theta is registered only in newtonian gauge (idx_theta >= 0); synchronous
  // gauge derives cdm velocity from the metric, so idx_theta == -1 there.
  if (old_l.idx_theta >= 0 && new_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> CDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  // CDM's Omega0 is the resolved cdm slot of the coupled-species budget
  // (parser block in InputModule::ReadCoupledOmegaBudget applied the user's
  // Omega_cdm / omega_cdm, then any f_idm_dr / f_idm_drmd subtractions, then
  // the synchronous-gauge minimum). Absent slot → CDM absent.
  if (!ctx.omega_budget || !ctx.omega_budget->cdm.has_value())
    return result;
  const double omega0 = *ctx.omega_budget->cdm;
  if (omega0 == 0.0)
    return result;
  result.push_back({"CDM", std::make_unique<CDMSpecies>(*ctx.pba, omega0)});
  return result;
}
