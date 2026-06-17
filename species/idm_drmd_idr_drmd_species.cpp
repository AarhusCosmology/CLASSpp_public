#include "idm_drmd_idr_drmd_species.h"

#include <cmath>

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
  w.Add("G_over_aH_drmd", pvecback[index_bg_G_over_aH_drmd_]);
}

void IDM_DRMD_IDR_DRMD_Species::RegisterBackgroundIndices(int& index_bg) {
  CompositeSpecies::RegisterBackgroundIndices(index_bg);  // children: idm_drmd, idr_drmd
  class_define_index(index_bg_G_over_aH_drmd_, _TRUE_, index_bg, 1);
}

void IDM_DRMD_IDR_DRMD_Species::ComputeIdmDrmd(
    double a, double rho_idm_over_rho_idr, double* Rint, double* csp2, double* Gint) const {
  const double z         = 1.0 / a - 1.0;
  const double R_int_tmp = 3.0 / 4.0 * rho_idm_over_rho_idr;
  *Rint                  = R_int_tmp;
  *csp2                  = 1.0 / 3.0 / (1.0 + R_int_tmp);
  if ((1.0 + z_stop_) / (1.0 + z) > 100)  // avoid exp() overflow
    *Gint = 0.;
  else
    *Gint = Gamma0_drmd_ic_ / R_int_tmp * exp(-(1.0 + z_stop_) / (1.0 + z));
}

void IDM_DRMD_IDR_DRMD_Species::InitializeDrmdBackground(
    double rho_tot, double H, double a, double a_today, const double* pvecback) {
  const double rho_idr = idr_drmd_->Rho(pvecback);
  const double rho_idm = idm_drmd_->Rho(pvecback);
  f_idr_drmd_          = rho_idr / rho_tot;
  Gamma0_drmd_ic_      = 0.;
  if (rho_idm > 0. && rho_idr > 0.)
    Gamma0_drmd_ic_ = 3. / 4. * G_over_aH_drmd_ * rho_idm / rho_idr * a / a_today * H;
}

void IDM_DRMD_IDR_DRMD_Species::FinalizeBackground(double a,
                                                   double H,
                                                   const double* /*pvecback_B*/,
                                                   double* pvecback) {
  double Rint, csp2, Gint;
  const double a_rel = a / pba_.a_today;
  ComputeIdmDrmd(a, idm_drmd_->Rho(pvecback) / idr_drmd_->Rho(pvecback), &Rint, &csp2, &Gint);
  pvecback[index_bg_G_over_aH_drmd_] = Gint / (H * a_rel);
}

void IDM_DRMD_IDR_DRMD_Species::ProcessBackgroundTable(const double* background_table,
                                                       int n_rows,
                                                       int row_stride,
                                                       const double* z_table) {
  // Decoupling redshift: the row where G_over_aH is closest to 1.
  for (int i = 0; i < n_rows; i++) {
    const double g = background_table[i * row_stride + index_bg_G_over_aH_drmd_];
    if (pow(g - 1.0, 2.0) < pow(G_over_aH_tmp_ - 1.0, 2.0)) {
      G_over_aH_tmp_ = g;
      z_dec_drmd_    = z_table[i];
    }
  }
}

void IDM_DRMD_IDR_DRMD_Species::RegisterTransferSourceIndices(int& index_tp,
                                                              const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_idm_drmd_, ctx.wants_density && has_idm_drmd(), index_tp, 1);
  class_define_index(index_tp_delta_idr_drmd_, ctx.wants_density && has_idr_drmd(), index_tp, 1);
  class_define_index(index_tp_theta_idm_drmd_, ctx.wants_velocity && has_idm_drmd(), index_tp, 1);
  class_define_index(index_tp_theta_idr_drmd_, ctx.wants_velocity && has_idr_drmd(), index_tp, 1);
}

void IDM_DRMD_IDR_DRMD_Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                       double* y,
                                                       const PerturbIcContext& ctx) {
  // Children have no IC logic; all ICs live on this composite.
  const PerturbationsModule* mod = ctx.p_mod;
  if (ctx.index_ic != mod->index_ic_ad_)
    return;

  const auto& my_lay      = static_cast<const PerturbLayout&>(base);
  const auto& idm_drm_lay = my_lay.idm_drmd;
  const auto& idr_drm_lay = my_lay.idr_drmd;

  if (has_idr_drmd()) {
    if (idr_drm_lay.idx_delta >= 0)
      y[idr_drm_lay.idx_delta] = ctx.delta_g_ic;
    if (idr_drm_lay.idx_theta >= 0)
      y[idr_drm_lay.idx_theta] = ctx.theta_g_ic;
  }

  if (has_idm_drmd()) {
    if (idm_drm_lay.idx_delta >= 0)
      y[idm_drm_lay.idx_delta] = 3. / 4. * ctx.delta_g_ic;

    if (idm_drm_lay.idx_theta >= 0) {
      if (has_idr_drmd()) {
        if (ctx.ppw->approx[ctx.ppw->index_ap_tca_idm_drmd] == (int) tca_idm_drmd_on) {
          y[idm_drm_lay.idx_theta] = (idr_drm_lay.idx_theta >= 0) ? y[idr_drm_lay.idx_theta] : 0.;
        }
        else {
          double Rint, csp2, Gint;
          auto* bgm = ctx.p_mod->GetBackgroundModule().get();
          ComputeIdmDrmd(ctx.ppw->pvecback[bgm->index_bg_a_],
                         idm_drmd_->Rho(ctx.ppw->pvecback) / idr_drmd_->Rho(ctx.ppw->pvecback),
                         &Rint,
                         &csp2,
                         &Gint);
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

void IDM_DRMD_IDR_DRMD_Species::FillSources(const BaseSpecies::PerturbLayout& base,
                                            const double* y,
                                            const double* /*dy*/,
                                            PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;

  // These sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  auto set_source = [&](int index_tp, double value) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp, ctx.index_tau, ctx.index_k, value);
  };

  // ── IDM_DRMD ──────────────────────────────────────────────────────────────
  const auto& my_src_lay      = static_cast<const PerturbLayout&>(base);
  const auto& idm_drm_src_lay = my_src_lay.idm_drmd;
  const auto& idr_drm_src_lay = my_src_lay.idr_drmd;
  if (index_tp_delta_idm_drmd_ >= 0) {
    set_source(index_tp_delta_idm_drmd_,
               y[idm_drm_src_lay.idx_delta] +
                   3. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  if (index_tp_theta_idm_drmd_ >= 0) {
    set_source(index_tp_theta_idm_drmd_,
               y[idm_drm_src_lay.idx_theta] + ctx.theta_shift);  // N-body gauge correction
  }

  // ── IDR_DRMD ──────────────────────────────────────────────────────────────
  if (index_tp_delta_idr_drmd_ >= 0) {
    set_source(index_tp_delta_idr_drmd_,
               y[idr_drm_src_lay.idx_delta] +
                   4. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // The original perturb_sources_member allocates this slot but never writes it.
  // Write zero explicitly to avoid relying on zero-initialization of the source table.
  if (index_tp_theta_idr_drmd_ >= 0) {
    set_source(index_tp_theta_idr_drmd_, 0.);
  }
}

void IDM_DRMD_IDR_DRMD_Species::WriteOutputColumns(
    PerturbColumnWriter& w,
    const PerturbationsModule& mod,
    enum file_format fmt,
    BaseSpecies::TransferColumnSection section) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers) {
      w.Add("d_idm_drmd", index_tp_delta_idm_drmd_, index_tp_delta_idm_drmd_ >= 0);
      w.Add("d_idr_drmd", index_tp_delta_idr_drmd_, index_tp_delta_idr_drmd_ >= 0);
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers) {
      w.Add("t_idm_drmd", index_tp_theta_idm_drmd_, index_tp_theta_idm_drmd_ >= 0);
      w.Add("t_idr_drmd", index_tp_theta_idr_drmd_, index_tp_theta_idr_drmd_ >= 0);
    }
  }
}

void IDM_DRMD_IDR_DRMD_Species::PrintVariables(PerturbColumnWriter& w,
                                               double /*tau*/,
                                               const double* y,
                                               const PerturbationsModule& mod,
                                               const perturb_workspace* ppw) const {
  double delta_idm_drmd = 0., theta_idm_drmd = 0.;
  double delta_idr_drmd = 0., theta_idr_drmd = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
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

    if (has_idm_drmd()) {
      delta_idm_drmd = y[idm_drm_pr_lay.idx_delta];
      theta_idm_drmd = y[idm_drm_pr_lay.idx_theta];
    }
    if (has_idr_drmd()) {
      delta_idr_drmd = y[idr_drm_pr_lay.idx_delta];
      theta_idr_drmd = y[idr_drm_pr_lay.idx_theta];
    }

    if (ppt->gauge == synchronous) {
      const double alpha = pvecmetric[ppw->index_mt_alpha];
      if (has_idm_drmd()) {
        delta_idm_drmd -= 3. * H * a * alpha;
        theta_idm_drmd += k * k * alpha;
      }
      if (has_idr_drmd()) {
        delta_idr_drmd -= 4. * H * a * alpha;
        theta_idr_drmd += k * k * alpha;
      }
    }
  }

  w.Add("delta_idr_drmd", delta_idr_drmd, has_idr_drmd());
  w.Add("theta_idr_drmd", theta_idr_drmd, has_idr_drmd());
  w.Add("delta_idm_drmd", delta_idm_drmd, has_idm_drmd());
  w.Add("theta_idm_drmd", theta_idm_drmd, has_idm_drmd());
}

IDM_DRMD_IDR_DRMD_Species::IDM_DRMD_IDR_DRMD_Species(const background& pba,
                                                     double omega0_idm_drmd,
                                                     double omega0_idr_drmd,
                                                     double f_idm_drmd,
                                                     double G_over_aH_drmd,
                                                     double delta_Neff_drmd,
                                                     double z_stop)
    : CompositeSpecies("IDM_DRMD_IDR_DRMD", BaseSpecies::EnergyType::Other), pba_(pba),
      f_idm_drmd_(f_idm_drmd), G_over_aH_drmd_(G_over_aH_drmd), delta_Neff_drmd_(delta_Neff_drmd),
      z_stop_(z_stop) {
  has_idm_drmd_ = (omega0_idm_drmd != 0.);
  has_idr_drmd_ = (omega0_idr_drmd != 0.);
  auto idm      = std::make_unique<IDM_DRMDSpecies>(pba, omega0_idm_drmd);
  auto idr      = std::make_unique<IDR_DRMDSpecies>(pba, omega0_idr_drmd);
  idm_drmd_     = idm.get();
  idr_drmd_     = idr.get();
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

void IDM_DRMD_IDR_DRMD_Species::PerturbSynchronousToNewtonian(
    const BaseSpecies::PerturbLayout& base, double* y, const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  idm_drmd_->PerturbSynchronousToNewtonian(my.idm_drmd, y, ctx);
  idr_drmd_->PerturbSynchronousToNewtonian(my.idr_drmd, y, ctx);
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
  const perturb_vector* pv        = ppw->pv.get();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double* pvecback = ppw->pvecback;

  const double rho_idm_drmd = idm_drmd_->Rho(pvecback);
  const double rho_idr_drmd = idr_drmd_->Rho(pvecback);

  // Guard against zero densities to avoid division by zero in ComputeIdmDrmd.
  if (rho_idm_drmd <= 0. || rho_idr_drmd <= 0.)
    return;

  double Rint, csp2, Gint;
  ComputeIdmDrmd(ctx.a, rho_idm_drmd / rho_idr_drmd, &Rint, &csp2, &Gint);

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

void IDM_DRMD_IDR_DRMD_Species::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y,
    double* new_y,
    const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  idm_drmd_->CopyPerturbationsAcrossSwitch(old_l.idm_drmd, new_l.idm_drmd, old_y, new_y, ctx);
  idr_drmd_->CopyPerturbationsAcrossSwitch(old_l.idr_drmd, new_l.idr_drmd, old_y, new_y, ctx);
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> IDM_DRMD_IDR_DRMD_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  // Read both slots from the resolved coupled-species budget.  Missing budget
  // (shooting-guess fallback) or both slots absent → composite is absent.
  if (!ctx.omega_budget)
    return result;
  const double omega0_idm_drmd = ctx.omega_budget->idm_drmd.value_or(0.);
  const double omega0_idr_drmd = ctx.omega_budget->idr_drmd.value_or(0.);
  if (omega0_idm_drmd != 0. || omega0_idr_drmd != 0.) {
    const auto* ci               = ctx.coupled_inputs;
    const double f_idm_drmd      = ci ? ci->f_idm_drmd : 0.;
    const double G_over_aH_drmd  = ci ? ci->G_over_aH_drmd_ini : 0.;
    const double delta_Neff_drmd = ci ? ci->delta_Neff_drmd : 0.;
    const double z_stop          = ci ? ci->z_stop : 0.;
    result.push_back({"IDM_DRMD_IDR_DRMD",
                      std::make_unique<IDM_DRMD_IDR_DRMD_Species>(*ctx.pba,
                                                                  omega0_idm_drmd,
                                                                  omega0_idr_drmd,
                                                                  f_idm_drmd,
                                                                  G_over_aH_drmd,
                                                                  delta_Neff_drmd,
                                                                  z_stop)});
  }
  return result;
}
