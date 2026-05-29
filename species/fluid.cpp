#include "fluid.h"

#include <cmath>

#include "background.h"
#include "background_column_writer.h"
#include "background_module.h"
#include "idm_dr_idr_species.h"
#include "idm_drmd_idr_drmd_species.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"

FluidSpecies::FluidSpecies(const background& pba, double omega0_fld)
    : BaseSpecies("Fluid", EnergyType::DarkEnergy), pba_(pba), Omega0_fld_(omega0_fld) {}

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

void FluidSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                               perturb_vector* /*pv*/,
                                               const precision* /*ppr*/,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  if (pba_.use_ppf == _FALSE_) {
    class_define_index(layout.idx_delta, _TRUE_, index_pt, 1);
    class_define_index(layout.idx_theta, _TRUE_, index_pt, 1);
  }
  else {
    class_define_index(layout.idx_Gamma, _TRUE_, index_pt, 1);
  }
}

void FluidSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                 double /*tau*/,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
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

    dy[layout.idx_delta] = -(1. + w_fld) * (y[layout.idx_theta] + metric_continuity) -
                           3. * (cs2 - w_fld) * a_prime_over_a * y[layout.idx_delta] -
                           9. * (1. + w_fld) * (cs2 - ca2) * a_prime_over_a * a_prime_over_a *
                               y[layout.idx_theta] / k2;

    dy[layout.idx_theta] = -(1. - 3. * cs2) * a_prime_over_a * y[layout.idx_theta] +
                           cs2 * k2 / (1. + w_fld) * y[layout.idx_delta] + metric_euler;
  }
  else {
    dy[layout.idx_Gamma] = ppw->Gamma_prime_fld;
  }
}

void FluidSpecies::FillSources(const BaseSpecies::PerturbLayout& /*layout*/,
                               const double* /*y*/,
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

void FluidSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                          double* y,
                                          const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  if (ctx.p_mod->GetBackground()->use_ppf == _TRUE_)
    return;
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return;

  double w_fld, dw_over_da, integral;
  auto* bgm = ctx.p_mod->GetBackgroundModule().get();
  class_call(bgm->background_w_fld(ctx.a, &w_fld, &dw_over_da, &integral),
             bgm->error_message_,
             bgm->error_message_);

  y[layout.idx_delta] = -ctx.ktau_two / 4. * (1. + w_fld) *
                        (4. - 3. * ctx.p_mod->GetBackground()->cs2_fld) /
                        (4. - 6. * w_fld + 3. * ctx.p_mod->GetBackground()->cs2_fld) *
                        ctx.ppr->curvature_ini * ctx.s2_squared;

  y[layout.idx_theta] = -ctx.k * ctx.ktau_three / 4. * ctx.p_mod->GetBackground()->cs2_fld /
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
      w.Add("d_fld", mod.index_tp_delta_fld_, _TRUE_);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_fld", mod.index_tp_theta_fld_, _TRUE_);
  }
}

void FluidSpecies::PrintVariables(PerturbColumnWriter& w,
                                  double /*tau*/,
                                  const double* /*y*/,
                                  const PerturbationsModule& mod,
                                  const perturb_workspace* ppw) const {
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

double FluidSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                           const perturb_vector* /*pv*/,
                           const double* y,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  return (layout.idx_delta >= 0) ? y[layout.idx_delta] : 0.;
}
double FluidSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                           const perturb_vector* /*pv*/,
                           const double* y,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  return (layout.idx_theta >= 0) ? y[layout.idx_theta] : 0.;
}
double FluidSpecies::DeltaP(const BaseSpecies::PerturbLayout& base,
                            const perturb_vector* /*pv*/,
                            const double* y,
                            const double* pvecback,
                            const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  // PPF uses a dedicated path (ComputePpf); this method is not called for
  // PPF — the module skips FluidSpecies in the main loop when use_ppf.
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return 0.;

  const double k2             = ppw->scalar_ctx.k2;
  const double a              = ppw->scalar_ctx.a;
  const double a_prime_over_a = pvecback[bgm_->index_bg_H_] * a;

  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld);
  const double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

  const double rho                  = Rho(pvecback);
  const double delta_rho_fld        = rho * y[layout.idx_delta];
  const double rho_plus_p_theta_fld = (1. + w_fld) * rho * y[layout.idx_theta];
  const double ca2_fld              = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;

  return pba_.cs2_fld * delta_rho_fld +
         (pba_.cs2_fld - ca2_fld) * (3. * a_prime_over_a * rho_plus_p_theta_fld / k2);
}
double FluidSpecies::RhoPlusPShear(const BaseSpecies::PerturbLayout& /*base*/,
                                   const perturb_vector* /*pv*/,
                                   const double* /*y*/,
                                   const double* /*pvecback*/,
                                   const perturb_workspace* /*ppw*/) const {
  return 0.;
}

int FluidSpecies::ComputeWFld(double a,
                              double* w_fld,
                              double* dw_over_da_fld,
                              double* integral_fld) const {
  double Omega_ede          = 0.;
  double dOmega_ede_over_da = 0.;
  double a_eq               = 0.0;

  /** - first, define the function w(a) */
  switch (pba_.fluid_equation_of_state) {
    case CLP:
      *w_fld = pba_.w0_fld + pba_.wa_fld * (1. - a / pba_.a_today);
      break;
    case EDE: {
      // Omega_ede(a) taken from eq. (10) in 1706.00730
      Omega_ede = (Omega0_fld_ - pba_.Omega_EDE * (1. - pow(a, -3. * pba_.w0_fld))) /
                      (Omega0_fld_ + (1. - Omega0_fld_) * pow(a, 3. * pba_.w0_fld)) +
                  pba_.Omega_EDE * (1. - pow(a, -3. * pba_.w0_fld));

      // d Omega_ede / d a taken analytically from the above
      dOmega_ede_over_da = -pba_.Omega_EDE * 3. * pba_.w0_fld * pow(a, -3. * pba_.w0_fld - 1.) /
                               (Omega0_fld_ + (1. - Omega0_fld_) * pow(a, 3. * pba_.w0_fld)) -
                           (Omega0_fld_ - pba_.Omega_EDE * (1. - pow(a, -3. * pba_.w0_fld))) *
                               (1. - Omega0_fld_) * 3. * pba_.w0_fld *
                               pow(a, 3. * pba_.w0_fld - 1.) /
                               pow(Omega0_fld_ + (1. - Omega0_fld_) * pow(a, 3. * pba_.w0_fld), 2) +
                           pba_.Omega_EDE * 3. * pba_.w0_fld * pow(a, -3. * pba_.w0_fld - 1.);

      // find a_equality (needed because EDE tracks first radiation, then matter)
      double Omega_r =
          pba_.Omega0_g *
          (1. +
           3.046 * 7. / 8. *
               pow(4. / 11.,
                   4. /
                       3.));  // assumes LambdaCDM + eventually massive neutrinos so light that they are relativistic at equality; needs to be generalised later on.
      double Omega_m = pba_.Omega0_b;
      if (auto* p = bgm_->all_species_.find("CDM"))
        Omega_m += (*p)->GetOmega0();
      if (auto* p = bgm_->all_species_.find("IDM_DR_IDR"))
        Omega_m += static_cast<const IDM_DR_IDR_Species&>(**p).idm_dr().GetOmega0();
      if (auto* p = bgm_->all_species_.find("IDM_DRMD_IDR_DRMD"))
        Omega_m += static_cast<const IDM_DRMD_IDR_DRMD_Species&>(**p).idm_drmd().GetOmega0();
      if (bgm_->all_species_.count("DCDM_DR"))
        class_stop(bgm_->error_message_,
                   "Early Dark Energy not compatible with decaying Dark Matter because we omitted "
                   "to code the calculation of a_eq in that case, but it would not be difficult to "
                   "add it if necessary, should be a matter of 5 minutes");
      a_eq = Omega_r / Omega_m;  // assumes a flat universe with a=1 today

      // w_ede(a) taken from eq. (11) in 1706.00730
      *w_fld = -dOmega_ede_over_da * a / Omega_ede / 3. / (1. - Omega_ede) + a_eq / 3. / (a + a_eq);
      break;
    }
  }

  /** - then, give the corresponding analytic derivative dw/da (used
      by perturbation equations; we could compute it numerically,
      but with a loss of precision; as long as there is a simple
      analytic expression of the derivative of the previous
      function, let's use it! */
  switch (pba_.fluid_equation_of_state) {
    case CLP:
      *dw_over_da_fld = -pba_.wa_fld / pba_.a_today;
      break;
    case EDE: {
      double d2Omega_ede_over_da2 = 0.;
      *dw_over_da_fld = -d2Omega_ede_over_da2 * a / 3. / (1. - Omega_ede) / Omega_ede -
                        dOmega_ede_over_da / 3. / (1. - Omega_ede) / Omega_ede +
                        dOmega_ede_over_da * dOmega_ede_over_da * a / 3. / (1. - Omega_ede) /
                            (1. - Omega_ede) / Omega_ede +
                        a_eq / 3. / (a + a_eq) / (a + a_eq);
      break;
    }
  }

  /** - finally, give the analytic solution of the following integral:
        \f$ \int_{a}^{a0} da 3(1+w_{fld})/a \f$. This is used in only
        one place, in the initial conditions for the background, and
        with a=a_ini. If your w(a) does not lead to a simple analytic
        solution of this integral, no worry: instead of writing
        something here, the best would then be to leave it equal to
        zero, and then in background_initial_conditions() you should
        implement a numerical calculation of this integral only for
        a=a_ini, using for instance Romberg integration. It should be
        fast, simple, and accurate enough. */
  switch (pba_.fluid_equation_of_state) {
    case CLP:
      *integral_fld = 3. * ((1. + pba_.w0_fld + pba_.wa_fld) * log(pba_.a_today / a) +
                            pba_.wa_fld * (a / pba_.a_today - 1.));
      break;
    case EDE:
      class_stop(bgm_->error_message_,
                 "EDE implementation not finished: to finish it, read the comments in background.c "
                 "just before this line\n");
      break;
  }

  /** note: of course you can generalise these formulas to anything,
      defining new parameters pba->w..._fld. Just remember that so
      far, HyRec explicitely assumes that w(a)= w0 + wa (1-a/a0); but
      Recfast does not assume anything */

  return _SUCCESS_;
}

void FluidSpecies::ComputePpf(double k,
                              double a,
                              double a_prime_over_a,
                              const precision* ppr,
                              const double* y,
                              perturb_workspace* ppw) const {
  const double a2 = a * a;
  const double k2 = k * k;

  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld);
  const double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

  double s2sq               = ppw->s_l[2] * ppw->s_l[2];
  double c_gamma_k_H_square = pow(pba_.c_gamma_over_c_fld * k / a_prime_over_a, 2) * pba_.cs2_fld;
  /** The equation is too stiff for Runge-Kutta when c_gamma_k_H_square is large.
      Use the asymptotic solution Gamma=Gamma'=0 in that case.
  */
  double Gamma_fld;
  if (c_gamma_k_H_square > ppr->c_gamma_k_H_square_max)
    Gamma_fld = 0.;
  else {
    const auto& layout = static_cast<const PerturbLayout&>(
        *ppw->pv->species_layouts[collection_index_]);
    Gamma_fld = y[layout.idx_Gamma];
  }

  double alpha, alpha_prime, metric_euler;
  if (ppw->scalar_ctx.gauge == synchronous) {
    alpha        = (y[ppw->pv->index_pt_eta] +
                    1.5 * a2 / k2 / s2sq *
                        (ppw->delta_rho + 3 * a_prime_over_a / k2 * ppw->rho_plus_p_theta) -
                    Gamma_fld) /
                   a_prime_over_a;
    alpha_prime  = -2. * a_prime_over_a * alpha + y[ppw->pv->index_pt_eta] -
                   4.5 * (a2 / k2) * ppw->rho_plus_p_shear;
    metric_euler = 0.;
  }
  else {
    alpha        = 0.;
    alpha_prime  = 0.;
    metric_euler = k2 * y[ppw->pv->index_pt_phi] - 4.5 * a2 * ppw->rho_plus_p_shear;
  }
  ppw->S_fld = Rho(ppw->pvecback) * (1. + w_fld) * 1.5 * a2 / k2 / a_prime_over_a *
               (ppw->rho_plus_p_theta / ppw->rho_plus_p_tot + k2 * alpha);
  // note that the last terms in the ratio do not include fld, that's correct, it's the whole point of the PPF scheme
  /** We must now check the stiffenss criterion again and set Gamma_prime_fld accordingly. */
  if (c_gamma_k_H_square > ppr->c_gamma_k_H_square_max) {
    ppw->Gamma_prime_fld = 0.;
  }
  else {
    ppw->Gamma_prime_fld = a_prime_over_a * (ppw->S_fld / (1. + c_gamma_k_H_square) -
                                             (1. + c_gamma_k_H_square) * Gamma_fld);
  }
  double Gamma_prime_plus_a_prime_over_a_Gamma = ppw->Gamma_prime_fld + a_prime_over_a * Gamma_fld;
  // delta and theta in both gauges gauge:
  ppw->rho_plus_p_theta_fld =
      Rho(ppw->pvecback) * (1. + w_fld) * ppw->rho_plus_p_theta / ppw->rho_plus_p_tot -
      k2 * 2. / 3. * a_prime_over_a / a2 / (1 + 4.5 * a2 / k2 / s2sq * ppw->rho_plus_p_tot) *
          (ppw->S_fld - Gamma_prime_plus_a_prime_over_a_Gamma / a_prime_over_a);
  ppw->delta_rho_fld = -2. / 3. * k2 * s2sq / a2 * Gamma_fld -
                       3 * a_prime_over_a / k2 * ppw->rho_plus_p_theta_fld;

  /** Now construct the pressure perturbation, see 1903.xxxxx. */
  /** Construct energy density and pressure for DE (_fld) and the rest (_t).
      Also compute derivatives. */
  double rho_fld       = Rho(ppw->pvecback);
  double p_fld         = w_fld * rho_fld;
  double rho_fld_prime = -3 * a_prime_over_a * (rho_fld + p_fld);
  double p_fld_prime   = w_prime_fld * rho_fld - 3 * a_prime_over_a * (1 + w_fld) * p_fld;
  double rho_t         = ppw->pvecback[bgm_->index_bg_rho_tot_] - rho_fld;
  double p_t           = ppw->pvecback[bgm_->index_bg_p_tot_] - p_fld;
  double rho_t_prime   = -3 * a_prime_over_a * (rho_t + p_t);
  double p_t_prime     = ppw->pvecback[bgm_->index_bg_p_tot_prime_] - p_fld_prime;
  /** Compute background quantities X,Y,Z and their derivatives. */
  double X       = c_gamma_k_H_square;
  double X_prime = -2 * X *
                   (a_prime_over_a +
                    ppw->pvecback[bgm_->index_bg_H_prime_] / ppw->pvecback[bgm_->index_bg_H_]);
  double Y       = 4.5 * a2 / k2 / s2sq * (rho_t + p_t);
  double Y_prime = Y * (2. * a_prime_over_a + (rho_t_prime + p_t_prime) / (rho_t + p_t));
  double Z       = 2. / 3. * k2 * ppw->pvecback[bgm_->index_bg_H_] / a;
  double Z_prime = Z * (ppw->pvecback[bgm_->index_bg_H_prime_] / ppw->pvecback[bgm_->index_bg_H_] -
                        a_prime_over_a);
  /** Construct theta_t and its derivative from the Euler equation */
  double theta_t       = ppw->rho_plus_p_theta / ppw->rho_plus_p_tot;
  double theta_t_prime = -a_prime_over_a * theta_t -
                         (p_t_prime * theta_t - k2 * ppw->delta_p + k2 * ppw->rho_plus_p_shear) /
                             ppw->rho_plus_p_tot +
                         metric_euler;
  double S             = ppw->S_fld;
  double S_prime       = -Z_prime / Z * S +
                         1. / Z * (rho_fld_prime + p_fld_prime) * (theta_t + k2 * alpha) +
                         1. / Z * (rho_fld + p_fld) * (theta_t_prime + k2 * alpha_prime);
  /** Analytic derivative of the equation for ppw->rho_plus_p_theta_fld above. */
  double rho_plus_p_theta_fld_prime =
      Z_prime * (S - 1. / (1. + Y) * (S / (1. + 1. / X) + Gamma_fld * X)) +
      Z * (S_prime + Y_prime / (1. + Y * Y + 2 * Y) * (S / (1. + 1. / X) + Gamma_fld * X) -
           1. / (1. + Y) *
               (S_prime / (1. + 1. / X) + S * X_prime / (1. + X * X + 2 * X) +
                ppw->Gamma_prime_fld * X + Gamma_fld * X_prime)) -
      k2 * alpha_prime * (rho_fld + p_fld) - k2 * alpha * (rho_fld_prime + p_fld_prime);

  /** We can finally compute the pressure perturbation using the Euler equation for theta_fld */
  ppw->delta_p_fld = (rho_plus_p_theta_fld_prime + 4 * a_prime_over_a * ppw->rho_plus_p_theta_fld -
                      (rho_fld + p_fld) * metric_euler) /
                     k2;
}

// ── Newtonian-gauge transform ─────────────────────────────────────────────────

void FluidSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                 double* y,
                                                 const PerturbIcContext& ctx) {
  // Non-PPF only: in PPF mode idx_delta/idx_theta are unregistered (== -1) and
  // the helper is a no-op. delta uses the universal ρ̇/ρ shift = -3(1+w)ℋα, which
  // corrects the historical opposite-sign bug.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}

void FluidSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                 const BaseSpecies::PerturbLayout& new_base,
                                                 const double* old_y,
                                                 double* new_y,
                                                 const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  // Standard fluid registers idx_delta/idx_theta; PPF registers idx_Gamma. The
  // unused slots are -1. Vector/tensor pv: all -1 -> no-op.
  if (old_l.idx_delta >= 0 && new_l.idx_delta >= 0)
    new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  if (old_l.idx_theta >= 0 && new_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
  if (old_l.idx_Gamma >= 0 && new_l.idx_Gamma >= 0)
    new_y[new_l.idx_Gamma] = old_y[old_l.idx_Gamma];
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> FluidSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // Closure-Fluid path: ConstructSpecies's Pass 2 hands us the budget-closure value.
  if (ctx.omega0_closure_override.has_value()) {
    const double omega0_fld = *ctx.omega0_closure_override;
    if (omega0_fld != 0.)
      result.push_back({"Fluid", std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld)});
    return result;
  }

  // Non-closure Fluid: read Omega_fld directly from the input file.
  double omega0_fld = 0.;
  if (!ctx.pfc->read_double("Omega_fld", omega0_fld) || omega0_fld == 0.)
    return result;

  result.push_back({"Fluid", std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld)});
  return result;
}
