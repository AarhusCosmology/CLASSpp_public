#include "idm_dr_idr_species.h"

#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics.h"
#include "thermodynamics_module.h"

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

  if (mod->GetBackground()->has_idm_dr == _TRUE_) {
    if (pv->index_pt_delta_idm_dr >= 0)
      y[pv->index_pt_delta_idm_dr] = 3. / 4. * ctx.delta_g_ic;
    if (pv->index_pt_theta_idm_dr >= 0)
      y[pv->index_pt_theta_idm_dr] = ctx.theta_ur;
  }
  if (mod->GetBackground()->has_idr == _TRUE_) {
    if (pv->index_pt_delta_idr >= 0)
      y[pv->index_pt_delta_idr] = ctx.delta_ur;
    if (pv->index_pt_theta_idr >= 0)
      y[pv->index_pt_theta_idr] = ctx.theta_ur;
    if (ppt->idr_nature == idr_free_streaming &&
        ((mod->GetBackground()->has_idm_dr == _FALSE_) ||
         (ctx.ppw->approx[ctx.ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off))) {
      if (pv->index_pt_shear_idr >= 0)
        y[pv->index_pt_shear_idr] = ctx.shear_ur;
      if (pv->index_pt_l3_idr >= 0)
        y[pv->index_pt_l3_idr] = ctx.l3_ur;
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
  if (p_mod->has_source_delta_idm_dr_ == _TRUE_) {
    set_source(p_mod->index_tp_delta_idm_dr_,
               y[pv->index_pt_delta_idm_dr] +
                   3. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  if (p_mod->has_source_theta_idm_dr_ == _TRUE_) {
    set_source(p_mod->index_tp_theta_idm_dr_,
               y[pv->index_pt_theta_idm_dr] + ctx.theta_shift);  // N-body gauge correction
  }

  // ── IDR ───────────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_idr_ == _TRUE_) {
    if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_off)
      set_source(p_mod->index_tp_delta_idr_,
                 y[pv->index_pt_delta_idr] +
                     4. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
    else
      set_source(p_mod->index_tp_delta_idr_,
                 ppw->rsa_delta_idr +
                     4. * ctx.a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  if (p_mod->has_source_theta_idr_ == _TRUE_) {
    if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_off)
      set_source(p_mod->index_tp_theta_idr_,
                 y[pv->index_pt_theta_idr] + ctx.theta_shift);  // N-body gauge correction
    else
      set_source(p_mod->index_tp_theta_idr_,
                 ppw->rsa_theta_idr + ctx.theta_shift);  // N-body gauge correction
  }
}

void IDM_DR_IDR_Species::WriteOutputColumns(PerturbColumnWriter& w,
                                            const PerturbationsModule& mod,
                                            enum file_format fmt,
                                            BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_) {
      w.Add("d_idm_dr", mod.index_tp_delta_idm_dr_, pba->has_idm_dr);
      w.Add("d_idr", mod.index_tp_delta_idr_, pba->has_idr);
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_) {
      w.Add("t_idm_dr", mod.index_tp_theta_idm_dr_, pba->has_idm_dr);
      w.Add("t_idr", mod.index_tp_theta_idr_, pba->has_idr);
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
  const background* pba = mod.GetBackground();

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

    if (pba->has_idm_dr == _TRUE_) {
      delta_idm_dr = y[pv->index_pt_delta_idm_dr];
      theta_idm_dr = y[pv->index_pt_theta_idm_dr];
    }

    if (pba->has_idr == _TRUE_) {
      if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_off) {
        delta_idr = y[pv->index_pt_delta_idr];
        theta_idr = y[pv->index_pt_theta_idr];
        if (ppt->idr_nature == idr_free_streaming) {
          if ((pba->has_idm_dr == _TRUE_) &&
              (ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_on)) {
            shear_idr = ppw->tca_shear_idm_dr;
          }
          else {
            shear_idr = y[pv->index_pt_shear_idr];
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
      if (pba->has_idm_dr == _TRUE_) {
        delta_idm_dr -= 3. * H * a * alpha;
        theta_idm_dr += k * k * alpha;
      }
      if (pba->has_idr == _TRUE_) {
        delta_idr -= 4. * H * a * alpha;
        theta_idr += k * k * alpha;
      }
    }
  }

  w.Add("delta_idr", delta_idr, pba->has_idr == _TRUE_);
  w.Add("theta_idr", theta_idr, pba->has_idr == _TRUE_);
  if (pba->has_idr == _TRUE_ && mod.GetPerturbs()->idr_nature == idr_free_streaming)
    w.Add("shear_idr", shear_idr, true);
  w.Add("delta_idm_dr", delta_idm_dr, pba->has_idm_dr == _TRUE_);
  w.Add("theta_idm_dr", theta_idm_dr, pba->has_idm_dr == _TRUE_);
}

IDM_DR_IDR_Species::IDM_DR_IDR_Species(const background& pba)
    : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other), pba_(pba) {
  auto idm = std::make_unique<IDM_DRSpecies>(pba);
  auto idr = std::make_unique<IDRSpecies>(pba);
  idm_dr_  = idm.get();
  idr_     = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
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

  const double Sinv         = 4. / 3. * rho_idr / rho_idm_dr;
  const double theta_idm_dr = y[pv->index_pt_theta_idm_dr];

  // Under RSA for IDR: IDR Boltzmann hierarchy is not evolved; use RSA-approximated
  // theta_idr to keep the IDM_DR drag correct. Do not write to IDR equations.
  if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on) {
    const double theta_idr         = ppw->rsa_theta_idr;
    dy[pv->index_pt_theta_idm_dr] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr) -
                                     ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                         y[pv->index_pt_delta_idm_dr];
    return;
  }

  const double theta_idr = (pv->index_pt_theta_idr >= 0) ? y[pv->index_pt_theta_idr] : 0.;

  if (ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off) {
    const thermo* pth    = pth_mod->GetThermodynamics();
    const double dmu_idr = pth->b_idr / pth->a_idm_dr * pba_.Omega0_idr / pba_.Omega0_idm_dr *
                           dmu_idm_dr;

    // IDM_DR velocity coupling
    dy[pv->index_pt_theta_idm_dr] -= Sinv * dmu_idm_dr * (theta_idm_dr - theta_idr) -
                                     ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                         y[pv->index_pt_delta_idm_dr];

    // IDR velocity coupling
    if (ctx.idr_nature == idr_free_streaming) {
      dy[pv->index_pt_theta_idr] += dmu_idm_dr * (theta_idm_dr - theta_idr);

      // IDR Compton collision terms in hierarchy l>=2
      const int l_max = pv->l_max_idr;
      for (int l = 2; l <= l_max; l++) {
        dy[pv->index_pt_delta_idr + l] -= (ppt->alpha_idm_dr[l - 2] * dmu_idm_dr +
                                           ppt->beta_idr[l - 2] * dmu_idr) *
                                          y[pv->index_pt_delta_idr + l];
      }
    }
  }
  else {
    // TCA on: compute tca_shear and tca_slip locally
    const double delta_idr = (pv->index_pt_delta_idr >= 0) ? y[pv->index_pt_delta_idr] : 0.;

    const double tca_shear_idm_dr = 0.5 * 8. / 15. / dmu_idm_dr / ppt->alpha_idm_dr[0] *
                                    (theta_idm_dr + ctx.metric_shear);

    const double tca_slip_idm_dr =
        (pth_mod->GetThermodynamics()->nindex_idm_dr - 2. / (1. + Sinv)) * ctx.a_prime_over_a *
            (theta_idm_dr - theta_idr) +
        1. / (1. + Sinv) / dmu_idm_dr *
            (-ctx.a_prime_over_a * theta_idm_dr +
             ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] * y[pv->index_pt_delta_idm_dr] +
             ctx.k2 * Sinv * (delta_idr / 4. - tca_shear_idm_dr));

    // ASSIGN (=): TCA replaces the free-streaming velocity written by the children
    dy[pv->index_pt_theta_idm_dr] = 1. / (1. + Sinv) *
                                        (-ctx.a_prime_over_a * theta_idm_dr +
                                         ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                             y[pv->index_pt_delta_idm_dr] +
                                         ctx.k2 * Sinv * (delta_idr / 4. - tca_shear_idm_dr)) +
                                    ctx.metric_euler + Sinv / (1. + Sinv) * tca_slip_idm_dr;

    dy[pv->index_pt_theta_idr] = 1. / (1. + Sinv) *
                                     (-ctx.a_prime_over_a * theta_idm_dr +
                                      ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] *
                                          y[pv->index_pt_delta_idm_dr] +
                                      ctx.k2 * Sinv * (delta_idr / 4. - tca_shear_idm_dr)) +
                                 ctx.metric_euler - 1. / (1. + Sinv) * tca_slip_idm_dr;
  }
}
