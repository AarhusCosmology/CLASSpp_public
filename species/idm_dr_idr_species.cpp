#include "idm_dr_idr_species.h"

#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics.h"
#include "thermodynamics_module.h"

std::optional<double> IDM_DR_IDR_Species::GetParam(const std::string& name) const {
  if (name == "T_idr")
    return idr().T_idr();
  if (name == "a_idm_dr")
    return idm_dr().a_idm_dr();
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

void IDM_DR_IDR_Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                double* y,
                                                const PerturbIcContext& ctx) {
  // Children have no IC logic; all ICs live on this composite.
  const PerturbationsModule* mod = ctx.p_mod;
  if (ctx.index_ic != mod->index_ic_ad_)
    return;

  const auto& my_lay     = static_cast<const PerturbLayout&>(base);
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
    if (idr_->idr_nature() == idr_free_streaming &&
        (!has_idm_dr() ||
         (ctx.ppw->approx[ctx.ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off))) {
      if (idr_lay.idx_shear >= 0)
        y[idr_lay.idx_shear] = ctx.shear_ur;
      if (idr_lay.idx_l3 >= 0)
        y[idr_lay.idx_l3] = ctx.l3_ur;
    }
  }
}

void IDM_DR_IDR_Species::FillSources(const BaseSpecies::PerturbLayout& base,
                                     const double* y,
                                     const double* /*dy*/,
                                     PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;

  // These sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  auto set_source = [&](int index_tp, double value) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp, ctx.index_tau, ctx.index_k, value);
  };

  // ── IDM_DR ────────────────────────────────────────────────────────────────
  const auto& my_lay     = static_cast<const PerturbLayout&>(base);
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
    const perturb_vector* pv = ppw->pv.get();
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
        if (idr_->idr_nature() == idr_free_streaming) {
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
  if (has_idr() && idr_->idr_nature() == idr_free_streaming)
    w.Add("shear_idr", shear_idr, true);
  w.Add("delta_idm_dr", delta_idm_dr, has_idm_dr());
  w.Add("theta_idm_dr", theta_idm_dr, has_idm_dr());
}

IDM_DR_IDR_Species::IDM_DR_IDR_Species(const background& pba,
                                       double omega0_idm_dr,
                                       double omega0_idr,
                                       double T_idr,
                                       int l_max_idr,
                                       double a_idm_dr,
                                       double nindex_idm_dr,
                                       double m_idm,
                                       double b_idr,
                                       int idr_nature,
                                       std::vector<double> alpha_idm_dr,
                                       std::vector<double> beta_idr)
    : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other), pba_(pba) {
  has_idm_dr_ = (omega0_idm_dr != 0.);
  has_idr_    = (omega0_idr != 0.);
  auto idm    = std::make_unique<IDM_DRSpecies>(pba, omega0_idm_dr, a_idm_dr, nindex_idm_dr, m_idm);
  auto idr    = std::make_unique<IDRSpecies>(pba,
                                             omega0_idr,
                                             has_idm_dr_,
                                             T_idr,
                                             l_max_idr,
                                             b_idr,
                                             idr_nature,
                                             std::move(alpha_idm_dr),
                                             std::move(beta_idr));
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
  const perturb_vector* pv        = ppw->pv.get();
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
    const double dmu_idr = idr().b_idr() / idm_dr().a_idm_dr() * idr().GetOmega0() /
                           idm_dr().GetOmega0() * dmu_idm_dr;

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
        dy[idr_lay.idx_delta + l] -= (idr_->alpha_idm_dr()[l - 2] * dmu_idm_dr +
                                      idr_->beta_idr()[l - 2] * dmu_idr) *
                                     y[idr_lay.idx_delta + l];
      }
    }
  }
  else {
    // TCA on: compute tca_shear and tca_slip locally
    const double delta_idr = (idr_lay.idx_delta >= 0) ? y[idr_lay.idx_delta] : 0.;

    const double tca_shear_idm_dr = 0.5 * 8. / 15. / dmu_idm_dr / idr_->alpha_idm_dr()[0] *
                                    (theta_idm_dr + ctx.metric_shear);

    const double tca_slip_idm_dr = (idm_dr().nindex_idm_dr() - 2. / (1. + Sinv)) *
                                       ctx.a_prime_over_a * (theta_idm_dr - theta_idr) +
                                   1. / (1. + Sinv) / dmu_idm_dr *
                                       (-ctx.a_prime_over_a * theta_idm_dr +
                                        ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                            y[idm_dr_lay.idx_delta] +
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
    // Re-homed from the monolith idm/idr parse block: an IDM_DR component
    // without any IDR density is unphysical (the interaction needs a partner).
    if (omega0_idm_dr > 0. && omega0_idr == 0.)
      throw std::invalid_argument(
          "You have requested interacting DM with DR, this requires a non-zero density of "
          "interacting DR. Please set either N_idr or xi_idr");

    const double T_idr_local = (ctx.coupled_inputs && ctx.coupled_inputs->T_idr)
                                   ? *ctx.coupled_inputs->T_idr
                                   : 0.;

    const int l_max_idr_local = ctx.ppr->l_max_idr;

    // ── Authoritative parse for IDM/IDR interaction parameters.  This is the
    //    sole site that reads a_idm_dr/nindex_idm_dr/m_idm (onto the IDM_DR
    //    child) and b_idr/idr_nature/alpha_idm_dr/beta_idr (onto the IDR
    //    child).  Thermodynamics and perturbation modules consume these values
    //    via the species' typed accessors. ──
    int input_verbose = 0;
    ctx.pfc->read_int("input_verbose", input_verbose);
    const double h2 = ctx.pba->h * ctx.pba->h;

    // a_idm_dr / a_dark / Gamma_0_nadm — at most one may be set.
    double a_idm_dr      = 0.;
    double nindex_idm_dr = 4.;  // thermo struct default
    double m_idm         = 1.e11;
    double b_idr         = 0.;
    int idr_nature       = idr_free_streaming;  // perturbs struct default (== 0)

    double a_param = 0., a_dark = 0., gamma0 = 0.;
    const bool f_a  = ctx.pfc->read_double("a_idm_dr", a_param);
    const bool f_ad = ctx.pfc->read_double("a_dark", a_dark);
    const bool f_g  = ctx.pfc->read_double("Gamma_0_nadm", gamma0);
    if (((int) f_a + (int) f_ad + (int) f_g) >= 2)
      throw std::invalid_argument(
          "In input file, you can only enter one of a_idm_dr, a_dark or Gamma_0_nadm, choose one");

    if (f_a) {
      a_idm_dr = a_param;
      if (input_verbose > 1)
        printf(
            "You passed a_idm_dr = a_dark = %e, this is equivalent to Gamma_0_nadm = %e in the "
            "NADM notation. \n",
            a_param,
            a_param * (4. / 3.) * (h2 * omega0_idr));
    }
    else if (f_ad) {
      a_idm_dr = a_dark;
      if (input_verbose > 1)
        printf(
            "You passed a_dark = a_idm_dr = %e, this is equivalent to Gamma_0_nadm = %e in the "
            "NADM notation. \n",
            a_dark,
            a_dark * (4. / 3.) * (h2 * omega0_idr));
    }
    else if (f_g) {
      a_idm_dr = gamma0 * (3. / 4.) / (h2 * omega0_idr);
      if (input_verbose > 1)
        printf(
            "You passed Gamma_0_nadm = %e, this is equivalent to a_idm_dr = a_dark = %e in the "
            "ETHOS notation. \n",
            gamma0,
            a_idm_dr);
    }

    if (f_g) {
      // If the user passed Gamma_0_nadm, assume they want nadm parameterisation.
      nindex_idm_dr = 0.;
      idr_nature    = idr_fluid;
      if (input_verbose > 1)
        printf("NADM requested. Defaulting on nindex_idm_dr = %e and idr_nature = fluid \n",
               nindex_idm_dr);
    }
    else {
      // nindex_dark / nindex_idm_dr — at most one may be set.
      double nd_dark = 0., nd_idm = 0.;
      const bool h_nd_dark = ctx.pfc->read_double("nindex_dark", nd_dark);
      const bool h_nd_idm  = ctx.pfc->read_double("nindex_idm_dr", nd_idm);
      if (h_nd_dark && h_nd_idm)
        throw std::invalid_argument(
            "In input file, you can only enter one of nindex_dark, nindex_idm_dr, choose one");
      if (h_nd_dark)
        nindex_idm_dr = nd_dark;
      if (h_nd_idm)
        nindex_idm_dr = nd_idm;

      std::string nature;
      if (ctx.pfc->read_string("idr_nature", nature)) {
        if (nature.find("free_streaming") != std::string::npos ||
            nature.find("Free_Streaming") != std::string::npos ||
            nature.find("Free_streaming") != std::string::npos ||
            nature.find("FREE_STREAMING") != std::string::npos) {
          idr_nature = idr_free_streaming;
        }
        if (nature.find("fluid") != std::string::npos ||
            nature.find("Fluid") != std::string::npos ||
            nature.find("FLUID") != std::string::npos) {
          idr_nature = idr_fluid;
        }
      }
    }

    // m_idm / m_dm — at most one may be set.
    {
      double m1 = 0., m2 = 0.;
      const bool found_m_idm = ctx.pfc->read_double("m_idm", m1);
      const bool found_m_dm  = ctx.pfc->read_double("m_dm", m2);
      if (found_m_idm && found_m_dm)
        throw std::invalid_argument(
            "In input file, you can only enter one of m_idm, m_dm, choose one");
      if (found_m_idm)
        m_idm = m1;
      if (found_m_dm)
        m_idm = m2;
    }

    // b_dark / b_idr — at most one may be set.
    {
      double b1 = 0., b2 = 0.;
      const bool found_b_dark = ctx.pfc->read_double("b_dark", b1);
      const bool found_b_idr  = ctx.pfc->read_double("b_idr", b2);
      if (found_b_dark && found_b_idr)
        throw std::invalid_argument(
            "In input file, you can only enter one of b_dark, b_idr, choose one");
      if (found_b_dark)
        b_idr = b1;
      if (found_b_idr)
        b_idr = b2;
    }

    // Number of l>=2 angular coefficients (one per multipole l = 2..l_max_idr).
    // Clamp at 0: l_max_idr <= 1 means no such multipoles, so the coefficient
    // arrays are empty (guards the (l_max_idr - 1) int->size_t underflow).
    const int n_idr_coeff = (ctx.ppr->l_max_idr > 1) ? ctx.ppr->l_max_idr - 1 : 0;

    // alpha_idm_dr / alpha_dark — list, default 1.5, resized to n_idr_coeff.
    std::vector<double> alpha_idm_dr;
    {
      bool found = ctx.pfc->read_list_of_doubles("alpha_idm_dr", alpha_idm_dr);
      if (!found)
        found = ctx.pfc->read_list_of_doubles("alpha_dark", alpha_idm_dr);
      if (found) {
        const int entries_read = static_cast<int>(alpha_idm_dr.size());
        if (entries_read != n_idr_coeff)
          alpha_idm_dr.resize(n_idr_coeff, alpha_idm_dr[entries_read - 1]);
      }
      else {
        alpha_idm_dr.assign(n_idr_coeff, 1.5);
      }
    }

    // beta_idr / beta_dark — list, default 1.5, resized to n_idr_coeff.
    std::vector<double> beta_idr;
    {
      bool found = ctx.pfc->read_list_of_doubles("beta_idr", beta_idr);
      if (!found)
        found = ctx.pfc->read_list_of_doubles("beta_dark", beta_idr);
      if (found) {
        const int entries_read = static_cast<int>(beta_idr.size());
        if (entries_read != n_idr_coeff)
          beta_idr.resize(n_idr_coeff, beta_idr[entries_read - 1]);
      }
      else {
        beta_idr.assign(n_idr_coeff, 1.5);
      }
    }

    result.push_back({"IDM_DR_IDR",
                      std::make_unique<IDM_DR_IDR_Species>(*ctx.pba,
                                                           omega0_idm_dr,
                                                           omega0_idr,
                                                           T_idr_local,
                                                           l_max_idr_local,
                                                           a_idm_dr,
                                                           nindex_idm_dr,
                                                           m_idm,
                                                           b_idr,
                                                           idr_nature,
                                                           std::move(alpha_idm_dr),
                                                           std::move(beta_idr))});
  }
  return result;
}

// ── MarkUsedInSources ─────────────────────────────────────────────────────────

void IDM_DR_IDR_Species::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                           const perturb_workspace* ppw,
                                           int* used_in_sources) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  /* IDM_DR_IDR has no tensor-mode slots; idx_l3 < 0 also when TCA is on.
     index_ap_rsa_idr and index_ap_tca_idm_dr are only defined in scalar mode,
     so we must not access them in tensor/vector mode. */
  if (my.idr.idx_l3 < 0)
    return;
  /* IDR l>=3 multipoles not needed in sources when rsa_idr is off, idr is
     free-streaming, and tca_idm_dr is off. */
  if (ppw->approx[ppw->index_ap_rsa_idr] != (int) rsa_idr_off)
    return;
  if (idr_->idr_nature() != idr_free_streaming)
    return;
  if (ppw->approx[ppw->index_ap_tca_idm_dr] != (int) tca_idm_dr_off)
    return;

  for (int idx = my.idr.idx_l3; idx <= my.idr.idx_delta + my.idr.l_max; ++idx)
    used_in_sources[idx] = _FALSE_;
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
