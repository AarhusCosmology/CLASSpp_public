#include "idm_drmd_idr_drmd_species.h"

#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

void IDM_DRMD_IDR_DRMD_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr_drmd", 0.);
  w.Add("(.)rho_idm_drmd", 0.);
  w.Add("G_over_aH_drmd", 0.);
}

void IDM_DRMD_IDR_DRMD_Species::WriteBackgroundData(const double* pvecback,
                                                    BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr_drmd", idr_drmd().Rho(pvecback));
  w.Add("(.)rho_idm_drmd", idm_drmd().Rho(pvecback));
  w.Add("G_over_aH_drmd", pvecback[bgm_->index_bg_G_over_aH_drmd_]);
}

void IDM_DRMD_IDR_DRMD_Species::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv             = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  if (ctx.index_ic != mod->index_ic_ad_)
    return;

  const auto& my_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  const auto& idm_drm_lay = my_lay.idm_drmd;
  const auto& idr_drm_lay = my_lay.idr_drmd;

  if (pba_.has_idr_drmd == _TRUE_) {
    if (idr_drm_lay.idx_delta >= 0)
      y[idr_drm_lay.idx_delta] = ctx.delta_g_ic;
    if (idr_drm_lay.idx_theta >= 0)
      y[idr_drm_lay.idx_theta] = ctx.theta_g_ic;
  }

  if (pba_.has_idm_drmd == _TRUE_) {
    if (idm_drm_lay.idx_delta >= 0)
      y[idm_drm_lay.idx_delta] = 3. / 4. * ctx.delta_g_ic;

    if (idm_drm_lay.idx_theta >= 0) {
      if (pba_.has_idr_drmd == _TRUE_) {
        if (ctx.ppw->approx[ctx.ppw->index_ap_tca_idm_drmd] == (int) tca_idm_drmd_on) {
          y[idm_drm_lay.idx_theta] = (idr_drm_lay.idx_theta >= 0) ? y[idr_drm_lay.idx_theta] : 0.;
        }
        else {
          double Rint, csp2, Gint;
          auto* bgm = ctx.p_mod->GetBackgroundModule().get();
          class_call(bgm->background_idm_drmd(ctx.ppw->pvecback[bgm->index_bg_a_],
                                              idm_drmd_->Rho(ctx.ppw->pvecback) /
                                                  idr_drmd_->Rho(ctx.ppw->pvecback),
                                              &Rint,
                                              &csp2,
                                              &Gint),
                     bgm->error_message_,
                     bgm->error_message_);
          y[idm_drm_lay.idx_theta] = Gint / (4. + Gint) *
                                     ((idr_drm_lay.idx_theta >= 0) ? y[idr_drm_lay.idx_theta] : 0.);
        }
      }
      else {
        y[idm_drm_lay.idx_theta] = 0.;
      }
    }
  }
}

void IDM_DRMD_IDR_DRMD_Species::FillSources(const double* y,
                                            const double* /*dy*/,
                                            PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const perturb_vector* pv   = ppw->pv;

  // These sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  auto set_source = [&](int index_tp, double value) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp, ctx.index_tau, ctx.index_k, value);
  };

  // ── IDM_DRMD ──────────────────────────────────────────────────────────────
  const auto& my_src_lay = static_cast<const PerturbLayout&>(
      *pv->species_layouts[collection_index_]);
  const auto& idm_drm_src_lay = my_src_lay.idm_drmd;
  const auto& idr_drm_src_lay = my_src_lay.idr_drmd;
  if (p_mod->has_source_delta_idm_drmd_ == _TRUE_) {
    set_source(p_mod->index_tp_delta_idm_drmd_,
               y[idm_drm_src_lay.idx_delta] +
                   3. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  if (p_mod->has_source_theta_idm_drmd_ == _TRUE_) {
    set_source(p_mod->index_tp_theta_idm_drmd_,
               y[idm_drm_src_lay.idx_theta] + ctx.theta_shift);  // N-body gauge correction
  }

  // ── IDR_DRMD ──────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_idr_drmd_ == _TRUE_) {
    set_source(p_mod->index_tp_delta_idr_drmd_,
               y[idr_drm_src_lay.idx_delta] +
                   4. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // The original perturb_sources_member allocates this slot but never writes it.
  // Write zero explicitly to avoid relying on zero-initialization of the source table.
  if (p_mod->has_source_theta_idr_drmd_ == _TRUE_) {
    set_source(p_mod->index_tp_theta_idr_drmd_, 0.);
  }
}

void IDM_DRMD_IDR_DRMD_Species::WriteOutputColumns(
    PerturbColumnWriter& w,
    const PerturbationsModule& mod,
    enum file_format fmt,
    BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_) {
      w.Add("d_idm_drmd", mod.index_tp_delta_idm_drmd_, pba->has_idm_drmd);
      w.Add("d_idr_drmd", mod.index_tp_delta_idr_drmd_, pba->has_idr_drmd);
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_) {
      w.Add("t_idm_drmd", mod.index_tp_theta_idm_drmd_, pba->has_idm_drmd);
      w.Add("t_idr_drmd", mod.index_tp_theta_idr_drmd_, pba->has_idr_drmd);
    }
  }
}

void IDM_DRMD_IDR_DRMD_Species::PrintVariables(PerturbColumnWriter& w,
                                               double /*tau*/,
                                               const double* y,
                                               const PerturbationsModule& mod,
                                               const perturb_workspace* ppw) const {
  const background* pba = mod.GetBackground();

  double delta_idm_drmd = 0., theta_idm_drmd = 0.;
  double delta_idr_drmd = 0., theta_idr_drmd = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;
    const double k           = ppw->scalar_ctx.k;
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const perturbs* ppt      = mod.GetPerturbs();

    const auto& my_pr_lay = static_cast<const PerturbLayout&>(
        *pv->species_layouts[collection_index_]);
    const auto& idm_drm_pr_lay = my_pr_lay.idm_drmd;
    const auto& idr_drm_pr_lay = my_pr_lay.idr_drmd;

    if (pba->has_idm_drmd == _TRUE_) {
      delta_idm_drmd = y[idm_drm_pr_lay.idx_delta];
      theta_idm_drmd = y[idm_drm_pr_lay.idx_theta];
    }
    if (pba->has_idr_drmd == _TRUE_) {
      delta_idr_drmd = y[idr_drm_pr_lay.idx_delta];
      theta_idr_drmd = y[idr_drm_pr_lay.idx_theta];
    }

    if (ppt->gauge == synchronous) {
      const double alpha = pvecmetric[ppw->index_mt_alpha];
      if (pba->has_idm_drmd == _TRUE_) {
        delta_idm_drmd -= 3. * H * a * alpha;
        theta_idm_drmd += k * k * alpha;
      }
      if (pba->has_idr_drmd == _TRUE_) {
        delta_idr_drmd -= 4. * H * a * alpha;
        theta_idr_drmd += k * k * alpha;
      }
    }
  }

  w.Add("delta_idr_drmd", delta_idr_drmd, pba->has_idr_drmd == _TRUE_);
  w.Add("theta_idr_drmd", theta_idr_drmd, pba->has_idr_drmd == _TRUE_);
  w.Add("delta_idm_drmd", delta_idm_drmd, pba->has_idm_drmd == _TRUE_);
  w.Add("theta_idm_drmd", theta_idm_drmd, pba->has_idm_drmd == _TRUE_);
}

IDM_DRMD_IDR_DRMD_Species::IDM_DRMD_IDR_DRMD_Species(const background& pba)
    : CompositeSpecies("IDM_DRMD_IDR_DRMD", BaseSpecies::EnergyType::Other), pba_(pba) {
  auto idm  = std::make_unique<IDM_DRMDSpecies>(pba);
  auto idr  = std::make_unique<IDR_DRMDSpecies>(pba);
  idm_drmd_ = idm.get();
  idr_drmd_ = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}

// ── Perturbation layout-based overrides ───────────────────────────────────────

void IDM_DRMD_IDR_DRMD_Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                            perturb_vector* pv,
                                                            const precision* ppr,
                                                            int& index_pt,
                                                            const perturb_workspace* ppw,
                                                            int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  idm_drmd_->RegisterPerturbationIndices(my.idm_drmd, pv, ppr, index_pt, ppw, gauge);
  idr_drmd_->RegisterPerturbationIndices(my.idr_drmd, pv, ppr, index_pt, ppw, gauge);
}

void IDM_DRMD_IDR_DRMD_Species::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                              double tau,
                                              const double* y,
                                              double* dy,
                                              const perturb_parameters_and_workspace& ppaw) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  idm_drmd_->PerturbDerivs(my.idm_drmd, tau, y, dy, ppaw);
  idr_drmd_->PerturbDerivs(my.idr_drmd, tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);
}

void IDM_DRMD_IDR_DRMD_Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& /*base*/,
                                                       double* y,
                                                       const PerturbIcContext& ctx) {
  // Children have no IC logic; all ICs live on this composite.
  ApplyInitialConditions(y, ctx);
}

double IDM_DRMD_IDR_DRMD_Species::Delta(const BaseSpecies::PerturbLayout& base,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rho_idm = idm_drmd_->Rho(pvecback);
  const double rho_idr = idr_drmd_->Rho(pvecback);
  const double rho_tot = rho_idm + rho_idr;
  if (rho_tot <= 0.)
    return 0.;
  return (rho_idm * idm_drmd_->Delta(my.idm_drmd, pv, y, pvecback, ppw) +
          rho_idr * idr_drmd_->Delta(my.idr_drmd, pv, y, pvecback, ppw)) /
         rho_tot;
}

double IDM_DRMD_IDR_DRMD_Species::Theta(const BaseSpecies::PerturbLayout& base,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rpp_idm = idm_drmd_->Rho(pvecback) + idm_drmd_->P(pvecback);
  const double rpp_idr = idr_drmd_->Rho(pvecback) + idr_drmd_->P(pvecback);
  const double rpp_tot = rpp_idm + rpp_idr;
  if (rpp_tot <= 0.)
    return 0.;
  return (rpp_idm * idm_drmd_->Theta(my.idm_drmd, pv, y, pvecback, ppw) +
          rpp_idr * idr_drmd_->Theta(my.idr_drmd, pv, y, pvecback, ppw)) /
         rpp_tot;
}

double IDM_DRMD_IDR_DRMD_Species::MatterRhoDelta(const perturb_vector* pv,
                                                 const double* y,
                                                 const double* pvecback,
                                                 const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return idm_drmd_->IsMatterSpecies()
             ? idm_drmd_->Rho(pvecback) * idm_drmd_->Delta(my.idm_drmd, pv, y, pvecback, ppw)
             : 0.;
}

double IDM_DRMD_IDR_DRMD_Species::MatterRhoPlusPTheta(const perturb_vector* pv,
                                                      const double* y,
                                                      const double* pvecback,
                                                      const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return idm_drmd_->IsMatterSpecies() ? (idm_drmd_->Rho(pvecback) + idm_drmd_->P(pvecback)) *
                                            idm_drmd_->Theta(my.idm_drmd, pv, y, pvecback, ppw)
                                      : 0.;
}

double IDM_DRMD_IDR_DRMD_Species::DeltaP(const BaseSpecies::PerturbLayout& base,
                                         const perturb_vector* pv,
                                         const double* y,
                                         const double* pvecback,
                                         const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  // IDM_DRMD has DeltaP == 0; only IDR_DRMD contributes.
  return idr_drmd_->DeltaP(my.idr_drmd, pv, y, pvecback, ppw);
}

double IDM_DRMD_IDR_DRMD_Species::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                                const perturb_vector* pv,
                                                const double* y,
                                                const double* pvecback,
                                                const perturb_workspace* ppw) const {
  // IDM_DRMD and IDR_DRMD both have RhoPlusPShear == 0.
  (void) base;
  (void) pv;
  (void) y;
  (void) pvecback;
  (void) ppw;
  return 0.;
}

void IDM_DRMD_IDR_DRMD_Species::AddCouplingDerivs(double /*tau*/,
                                                  const double* y,
                                                  double* dy,
                                                  const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  auto* bgm              = ppaw.perturbations_module->GetBackgroundModule().get();
  const double* pvecback = ppw->pvecback;

  const double rho_idm_drmd = idm_drmd_->Rho(pvecback);
  const double rho_idr_drmd = idr_drmd_->Rho(pvecback);

  // Guard against zero densities to avoid division by zero in background_idm_drmd.
  if (rho_idm_drmd <= 0. || rho_idr_drmd <= 0.)
    return;

  double Rint, csp2, Gint;
  class_call(bgm->background_idm_drmd(ctx.a, rho_idm_drmd / rho_idr_drmd, &Rint, &csp2, &Gint),
             bgm->error_message_,
             bgm->error_message_);

  const auto& my_acc_lay = static_cast<const PerturbLayout&>(
      *pv->species_layouts[collection_index_]);
  const auto& idm_drm_acc_lay = my_acc_lay.idm_drmd;
  const auto& idr_drm_acc_lay = my_acc_lay.idr_drmd;

  const double theta_idm = y[idm_drm_acc_lay.idx_theta];
  const double theta_idr = y[idr_drm_acc_lay.idx_theta];
  const double delta_idr = y[idr_drm_acc_lay.idx_delta];

  if (ppw->approx[ppw->index_ap_tca_idm_drmd] == (int) tca_idm_drmd_on) {
    // TCA on: ASSIGN theta_idm_drmd (replaces the free-streaming term written by child)
    double GdDelta = 3. * csp2 * (ctx.a_prime_over_a * theta_idm + ctx.k2 * delta_idr / 4.);
    dy[idm_drm_acc_lay.idx_theta] = 0.25 * ctx.k2 * delta_idr + ctx.metric_euler - GdDelta * Rint;
    // IDR_DRMD: subtract coupling from free-streaming equation already written
    dy[idr_drm_acc_lay.idx_theta] -= GdDelta * Rint;
  }
  else {
    // TCA off: add coupling increments
    dy[idm_drm_acc_lay.idx_theta] += Gint * (theta_idr - theta_idm);
    dy[idr_drm_acc_lay.idx_theta] -= Gint * Rint * (theta_idr - theta_idm);
  }
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> IDM_DRMD_IDR_DRMD_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_idm_drmd == _TRUE_ || ctx.pba->has_idr_drmd == _TRUE_) {
    result.push_back({"IDM_DRMD_IDR_DRMD", std::make_unique<IDM_DRMD_IDR_DRMD_Species>(*ctx.pba)});
  }
  return result;
}
