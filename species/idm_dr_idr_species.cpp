#include "idm_dr_idr_species.h"

#include <cmath>
#include <optional>
#include <string>

#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics.h"
#include "thermodynamics_module.h"

std::optional<double> IDM_DR_IDR_Species::GetParam(const std::string& name) const {
  if (name == "T_idr")
    return idr().T_idr();
  return std::nullopt;
}

void IDM_DR_IDR_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr", 0.);
  w.Add("(.)rho_idm_dr", 0.);
}

void IDM_DR_IDR_Species::WriteBackgroundData(const double* pvecback,
                                             BackgroundColumnWriter& w) const {
  w.Add("(.)rho_idr", idr().Rho(pvecback));
  w.Add("(.)rho_idm_dr", idm_dr().Rho(pvecback));
}

void IDM_DR_IDR_Species::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv             = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  const perturbs* ppt            = mod->GetPerturbs();
  if (ctx.index_ic != mod->index_ic_ad_)
    return;

  const auto& my_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  const auto& idm_dr_lay = my_lay.idm_dr;
  const auto& idr_lay    = my_lay.idr;

  if (has_idm_dr()) {
    if (idm_dr_lay.idx_delta >= 0)
      y[idm_dr_lay.idx_delta] = 3. / 4. * ctx.delta_g_ic;
    if (idm_dr_lay.idx_theta >= 0)
      y[idm_dr_lay.idx_theta] = ctx.theta_ur;
  }
  if (has_idr()) {
    if (idr_lay.idx_delta >= 0)
      y[idr_lay.idx_delta] = ctx.delta_ur;
    if (idr_lay.idx_theta >= 0)
      y[idr_lay.idx_theta] = ctx.theta_ur;
    if (ppt->idr_nature == idr_free_streaming &&
        (!has_idm_dr() ||
         (ctx.ppw->approx[ctx.ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off))) {
      if (idr_lay.idx_shear >= 0)
        y[idr_lay.idx_shear] = ctx.shear_ur;
      if (idr_lay.idx_l3 >= 0)
        y[idr_lay.idx_l3] = ctx.l3_ur;
    }
  }
}

void IDM_DR_IDR_Species::FillSources(const double* y,
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

  // ── IDM_DR ────────────────────────────────────────────────────────────────
  const auto& my_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  const auto& idm_dr_lay = my_lay.idm_dr;
  const auto& idr_lay    = my_lay.idr;
  if (p_mod->has_source_delta_idm_dr_ == _TRUE_) {
    set_source(p_mod->index_tp_delta_idm_dr_,
               y[idm_dr_lay.idx_delta] +
                   3. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  if (p_mod->has_source_theta_idm_dr_ == _TRUE_) {
    set_source(p_mod->index_tp_theta_idm_dr_,
               y[idm_dr_lay.idx_theta] + ctx.theta_shift);  // N-body gauge correction
  }

  // ── IDR ───────────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_idr_ == _TRUE_) {
    if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_off)
      set_source(p_mod->index_tp_delta_idr_,
                 y[idr_lay.idx_delta] +
                     4. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
    else
      set_source(p_mod->index_tp_delta_idr_,
                 ppw->rsa_delta_idr +
                     4. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  if (p_mod->has_source_theta_idr_ == _TRUE_) {
    if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_off)
      set_source(p_mod->index_tp_theta_idr_,
                 y[idr_lay.idx_theta] + ctx.theta_shift);  // N-body gauge correction
    else
      set_source(p_mod->index_tp_theta_idr_,
                 ppw->rsa_theta_idr + ctx.theta_shift);  // N-body gauge correction
  }
}

void IDM_DR_IDR_Species::WriteOutputColumns(PerturbColumnWriter& w,
                                            const PerturbationsModule& mod,
                                            enum file_format fmt,
                                            BaseSpecies::TransferColumnSection section) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_) {
      w.Add("d_idm_dr", mod.index_tp_delta_idm_dr_, has_idm_dr() ? _TRUE_ : _FALSE_);
      w.Add("d_idr", mod.index_tp_delta_idr_, has_idr() ? _TRUE_ : _FALSE_);
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_) {
      w.Add("t_idm_dr", mod.index_tp_theta_idm_dr_, has_idm_dr() ? _TRUE_ : _FALSE_);
      w.Add("t_idr", mod.index_tp_theta_idr_, has_idr() ? _TRUE_ : _FALSE_);
    }
  }
  else if (fmt == camb_format) {
    if (section != TransferColumnSection::velocity)
      w.Add("-T_idm_dr/k2", mod.index_tp_delta_idm_dr_, _TRUE_);
  }
}

void IDM_DR_IDR_Species::PrintVariables(PerturbColumnWriter& w,
                                        double /*tau*/,
                                        const double* y,
                                        const PerturbationsModule& mod,
                                        const perturb_workspace* ppw) const {
  double delta_idm_dr = 0., theta_idm_dr = 0.;
  double delta_idr = 0., theta_idr = 0., shear_idr = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;
    const double k           = ppw->scalar_ctx.k;
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const perturbs* ppt      = mod.GetPerturbs();

    const auto& my_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
    const auto& idm_dr_lay = my_lay.idm_dr;
    const auto& idr_lay    = my_lay.idr;

    if (has_idm_dr()) {
      delta_idm_dr = y[idm_dr_lay.idx_delta];
      theta_idm_dr = y[idm_dr_lay.idx_theta];
    }

    if (has_idr()) {
      if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_off) {
        delta_idr = y[idr_lay.idx_delta];
        theta_idr = y[idr_lay.idx_theta];
        if (ppt->idr_nature == idr_free_streaming) {
          if (has_idm_dr() && (ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_on)) {
            shear_idr = idr_->TcaShearIdr(idr_lay, y, ppw);
          }
          else {
            shear_idr = y[idr_lay.idx_shear];
          }
        }
      }
      else {
        delta_idr = ppw->rsa_delta_idr;
        theta_idr = ppw->rsa_theta_idr;
        shear_idr = 0.;
      }
    }

    if (ppt->gauge == synchronous) {
      const double alpha = pvecmetric[ppw->index_mt_alpha];
      if (has_idm_dr()) {
        delta_idm_dr -= 3. * H * a * alpha;
        theta_idm_dr += k * k * alpha;
      }
      if (has_idr()) {
        delta_idr -= 4. * H * a * alpha;
        theta_idr += k * k * alpha;
      }
    }
  }

  w.Add("delta_idr", delta_idr, has_idr());
  w.Add("theta_idr", theta_idr, has_idr());
  if (has_idr() && mod.GetPerturbs()->idr_nature == idr_free_streaming)
    w.Add("shear_idr", shear_idr, true);
  w.Add("delta_idm_dr", delta_idm_dr, has_idm_dr());
  w.Add("theta_idm_dr", theta_idm_dr, has_idm_dr());
}

IDM_DR_IDR_Species::IDM_DR_IDR_Species(
    const background& pba, double omega0_idm_dr, double omega0_idr, double T_idr, int l_max_idr)
    : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other), pba_(pba) {
  has_idm_dr_ = (omega0_idm_dr != 0.);
  has_idr_    = (omega0_idr != 0.);
  auto idm    = std::make_unique<IDM_DRSpecies>(pba, omega0_idm_dr);
  auto idr    = std::make_unique<IDRSpecies>(pba, omega0_idr, has_idm_dr_, T_idr, l_max_idr);
  idm_dr_     = idm.get();
  idr_        = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}

// ── Perturbation layout-based overrides ───────────────────────────────────────

void IDM_DR_IDR_Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                     perturb_vector* pv,
                                                     const precision* ppr,
                                                     int& index_pt,
                                                     const perturb_workspace* ppw,
                                                     int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  idm_dr_->RegisterPerturbationIndices(my.idm_dr, pv, ppr, index_pt, ppw, gauge);
  idr_->RegisterPerturbationIndices(my.idr, pv, ppr, index_pt, ppw, gauge);
}

void IDM_DR_IDR_Species::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                       double tau,
                                       const double* y,
                                       double* dy,
                                       const perturb_parameters_and_workspace& ppaw) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  idm_dr_->PerturbDerivs(my.idm_dr, tau, y, dy, ppaw);
  idr_->PerturbDerivs(my.idr, tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);
}

void IDM_DR_IDR_Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& /*base*/,
                                                double* y,
                                                const PerturbIcContext& ctx) {
  // Children have no IC logic; all ICs live on this composite.
  ApplyInitialConditions(y, ctx);
}

void IDM_DR_IDR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                       double* y,
                                                       const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  idm_dr_->PerturbSynchronousToNewtonian(my.idm_dr, y, ctx);
  idr_->PerturbSynchronousToNewtonian(my.idr, y, ctx);
}

double IDM_DR_IDR_Species::Delta(const BaseSpecies::PerturbLayout& base,
                                 const perturb_vector* pv,
                                 const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rho_idm = idm_dr_->Rho(pvecback);
  const double rho_idr = idr_->Rho(pvecback);
  const double rho_tot = rho_idm + rho_idr;
  if (rho_tot <= 0.)
    return 0.;
  return (rho_idm * idm_dr_->Delta(my.idm_dr, pv, y, pvecback, ppw) +
          rho_idr * idr_->Delta(my.idr, pv, y, pvecback, ppw)) /
         rho_tot;
}

double IDM_DR_IDR_Species::Theta(const BaseSpecies::PerturbLayout& base,
                                 const perturb_vector* pv,
                                 const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rpp_idm = idm_dr_->Rho(pvecback) + idm_dr_->P(pvecback);
  const double rpp_idr = idr_->Rho(pvecback) + idr_->P(pvecback);
  const double rpp_tot = rpp_idm + rpp_idr;
  if (rpp_tot <= 0.)
    return 0.;
  return (rpp_idm * idm_dr_->Theta(my.idm_dr, pv, y, pvecback, ppw) +
          rpp_idr * idr_->Theta(my.idr, pv, y, pvecback, ppw)) /
         rpp_tot;
}

double IDM_DR_IDR_Species::MatterRhoDelta(const perturb_vector* pv,
                                          const double* y,
                                          const double* pvecback,
                                          const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return idm_dr_->IsMatterSpecies()
             ? idm_dr_->Rho(pvecback) * idm_dr_->Delta(my.idm_dr, pv, y, pvecback, ppw)
             : 0.;
}

double IDM_DR_IDR_Species::MatterRhoPlusPTheta(const perturb_vector* pv,
                                               const double* y,
                                               const double* pvecback,
                                               const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return idm_dr_->IsMatterSpecies() ? (idm_dr_->Rho(pvecback) + idm_dr_->P(pvecback)) *
                                          idm_dr_->Theta(my.idm_dr, pv, y, pvecback, ppw)
                                    : 0.;
}

double IDM_DR_IDR_Species::DeltaP(const BaseSpecies::PerturbLayout& base,
                                  const perturb_vector* pv,
                                  const double* y,
                                  const double* pvecback,
                                  const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  // IDM_DR has DeltaP == 0; only IDR contributes.
  return idr_->DeltaP(my.idr, pv, y, pvecback, ppw);
}

double IDM_DR_IDR_Species::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                         const perturb_vector* pv,
                                         const double* y,
                                         const double* pvecback,
                                         const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  // IDM_DR has RhoPlusPShear == 0; only IDR contributes.
  return idr_->RhoPlusPShear(my.idr, pv, y, pvecback, ppw);
}

void IDM_DR_IDR_Species::AddCouplingDerivs(double /*tau*/,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  auto* pth_mod            = ppaw.perturbations_module->GetThermodynamicsModule().get();
  const double* pvecback   = ppw->pvecback;
  const double* pvecthermo = ppw->pvecthermo;
  auto* ppt                = ppaw.perturbations_module->GetPerturbs();

  const double dmu_idm_dr = pvecthermo[pth_mod->index_th_dmu_idm_dr_];
  const double rho_idm_dr = idm_dr_->Rho(pvecback);
  const double rho_idr    = idr_->Rho(pvecback);

  // No coupling terms if either sector is absent or coupling rate vanishes.
  // This also guards the Sinv = rho_idr/rho_idm_dr and 1/dmu_idm_dr divisions.
  if (rho_idm_dr <= 0. || rho_idr <= 0. || dmu_idm_dr <= 0.)
    return;

  const auto& my_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  const auto& idm_dr_lay    = my_lay.idm_dr;
  const auto& idr_lay       = my_lay.idr;
  const double Sinv         = 4. / 3. * rho_idr / rho_idm_dr;
  const double theta_idm_dr = y[idm_dr_lay.idx_theta];

  // Under RSA for IDR: IDR Boltzmann hierarchy is not evolved; use RSA-approximated
  // theta_idr to keep the IDM_DR drag correct. Do not write to IDR equations.
  if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on) {
    const double theta_idr    = ppw->rsa_theta_idr;
    dy[idm_dr_lay.idx_theta] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr) -
                                ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                    y[idm_dr_lay.idx_delta];
    return;
  }

  const double theta_idr = (idr_lay.idx_theta >= 0) ? y[idr_lay.idx_theta] : 0.;

  if (ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off) {
    const thermo* pth    = pth_mod->GetThermodynamics();
    const double dmu_idr = pth->b_idr / pth->a_idm_dr * idr().GetOmega0() / idm_dr().GetOmega0() *
                           dmu_idm_dr;

    // IDM_DR velocity coupling
    dy[idm_dr_lay.idx_theta] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr) -
                                ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                    y[idm_dr_lay.idx_delta];

    // IDR velocity coupling
    if (ctx.idr_nature == idr_free_streaming) {
      dy[idr_lay.idx_theta] += dmu_idm_dr * (theta_idm_dr - theta_idr);

      // IDR Compton collision terms in hierarchy l>=2
      const int l_max = idr_lay.l_max;
      for (int l = 2; l <= l_max; l++) {
        dy[idr_lay.idx_delta + l] -= (ppt->alpha_idm_dr[l - 2] * dmu_idm_dr +
                                      ppt->beta_idr[l - 2] * dmu_idr) *
                                     y[idr_lay.idx_delta + l];
      }
    }
  }
  else {
    // TCA on: compute tca_shear and tca_slip locally
    const double delta_idr = (idr_lay.idx_delta >= 0) ? y[idr_lay.idx_delta] : 0.;

    const double tca_shear_idm_dr = 0.5 * 8. / 15. / dmu_idm_dr / ppt->alpha_idm_dr[0] *
                                    (theta_idm_dr + ctx.metric_shear);

    const double tca_slip_idm_dr =
        (pth_mod->GetThermodynamics()->nindex_idm_dr - 2. / (1. + Sinv)) * ctx.a_prime_over_a *
            (theta_idm_dr - theta_idr) +
        1. / (1. + Sinv) / dmu_idm_dr *
            (-ctx.a_prime_over_a * theta_idm_dr +
             ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] * y[idm_dr_lay.idx_delta] +
             ctx.k2 * Sinv * (delta_idr / 4. - tca_shear_idm_dr));

    // ASSIGN (=): TCA replaces the free-streaming velocity written by the children
    dy[idm_dr_lay.idx_theta] = 1. / (1. + Sinv) *
                                   (-ctx.a_prime_over_a * theta_idm_dr +
                                    ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                        y[idm_dr_lay.idx_delta] +
                                    ctx.k2 * Sinv * (delta_idr / 4. - tca_shear_idm_dr)) +
                               ctx.metric_euler + Sinv / (1. + Sinv) * tca_slip_idm_dr;

    dy[idr_lay.idx_theta] = 1. / (1. + Sinv) *
                                (-ctx.a_prime_over_a * theta_idm_dr +
                                 ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                     y[idm_dr_lay.idx_delta] +
                                 ctx.k2 * Sinv * (delta_idr / 4. - tca_shear_idm_dr)) +
                            ctx.metric_euler - 1. / (1. + Sinv) * tca_slip_idm_dr;
  }
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> IDM_DR_IDR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  // Read both slots from the resolved coupled-species budget.  Missing budget
  // (shooting-guess fallback) or both slots absent → composite is absent.
  if (!ctx.omega_budget)
    return result;
  const double omega0_idm_dr = ctx.omega_budget->idm_dr.value_or(0.);
  const double omega0_idr    = ctx.omega_budget->idr.value_or(0.);
  if (omega0_idm_dr != 0. || omega0_idr != 0.) {
    // ── Parse T_idr (same logic as ReadCoupledOmegaBudget, kept local) ─────
    double T_idr_local = 0.;
    double stat_f_idr  = 7. / 8.;
    ctx.pfc->read_double("stat_f_idr", stat_f_idr);

    double N_idr = 0., N_dg = 0., xi_idr = 0.;
    bool flag_N_idr = ctx.pfc->read_double("N_idr", N_idr);
    bool flag_N_dg  = ctx.pfc->read_double("N_dg", N_dg);
    bool flag_xi    = ctx.pfc->read_double("xi_idr", xi_idr);

    if (flag_N_idr)
      T_idr_local = std::pow(N_idr / stat_f_idr * (7. / 8.) / std::pow(11. / 4., 4. / 3.),
                             1. / 4.) *
                    ctx.pba->T_cmb;
    else if (flag_N_dg)
      T_idr_local = std::pow(N_dg / stat_f_idr * (7. / 8.) / std::pow(11. / 4., 4. / 3.), 1. / 4.) *
                    ctx.pba->T_cmb;
    else if (flag_xi)
      T_idr_local = xi_idr * ctx.pba->T_cmb;

    const int l_max_idr_local = ctx.ppr->l_max_idr;

    result.push_back({"IDM_DR_IDR",
                      std::make_unique<IDM_DR_IDR_Species>(*ctx.pba,
                                                           omega0_idm_dr,
                                                           omega0_idr,
                                                           T_idr_local,
                                                           l_max_idr_local)});
  }
  return result;
}

void IDM_DR_IDR_Species::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                       const BaseSpecies::PerturbLayout& new_base,
                                                       const double* old_y,
                                                       double* new_y,
                                                       const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  idm_dr_->CopyPerturbationsAcrossSwitch(old_l.idm_dr, new_l.idm_dr, old_y, new_y, ctx);
  idr_->CopyPerturbationsAcrossSwitch(old_l.idr, new_l.idr, old_y, new_y, ctx);
}
