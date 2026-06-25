#include "ppf_fluid.h"

#include "background_module.h"
#include "perturbations_module.h"

PpfFluid::PpfFluid(const background& pba,
                   double omega0_fld,
                   equation_of_state fluid_eos,
                   double w0_fld,
                   double wa_fld,
                   double cs2_fld,
                   double Omega_EDE,
                   double c_gamma_over_c_fld)
    : FluidSpecies(pba, omega0_fld, fluid_eos, w0_fld, wa_fld, cs2_fld, Omega_EDE),
      c_gamma_over_c_fld_(c_gamma_over_c_fld) {}

void PpfFluid::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                           perturb_vector* /*pv*/,
                                           const precision* /*ppr*/,
                                           int& index_pt,
                                           const perturb_workspace* /*ppw*/,
                                           int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  class_define_index(layout.idx_Gamma, _TRUE_, index_pt, 1);
}

void PpfFluid::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                             double /*tau*/,
                             const double* /*y*/,
                             double* dy,
                             const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout   = static_cast<const PerturbLayout&>(base);
  dy[layout.idx_Gamma] = ppaw.ppw->Gamma_prime_fld;
}

void PpfFluid::FillSources(const BaseSpecies::PerturbLayout& /*base*/,
                           const double* /*y*/,
                           const double* /*dy*/,
                           PerturbSourceContext& ctx) const {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const double* pvecback     = ppw->pvecback.data();
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  if (index_tp_delta_ >= 0) {
    const double w_fld     = W(pvecback);
    const double delta_fld = ppw->delta_rho_fld / Rho(pvecback);
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_delta_,
                          ctx.index_tau,
                          ctx.index_k,
                          delta_fld + 3. * ctx.a_prime_over_a * (1. + w_fld) * ctx.theta_over_k2);
  }
  if (index_tp_theta_ >= 0) {
    const double w_fld     = W(pvecback);
    const double theta_fld = ppw->rho_plus_p_theta_fld / (1. + w_fld) / Rho(pvecback);
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_theta_,
                          ctx.index_tau,
                          ctx.index_k,
                          theta_fld + ctx.theta_shift);
  }
}

void PpfFluid::PrintVariables(PerturbColumnWriter& w,
                              double /*tau*/,
                              const double* /*y*/,
                              const PerturbationsModule& /*mod*/,
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

void PpfFluid::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                             const BaseSpecies::PerturbLayout& new_base,
                                             const double* old_y,
                                             double* new_y,
                                             const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_Gamma >= 0 && new_l.idx_Gamma >= 0)
    new_y[new_l.idx_Gamma] = old_y[old_l.idx_Gamma];
}

void PpfFluid::ComputePpf(double k,
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

  // The PPF closure is defined relative to the rest of the universe. The module's
  // runtime rho_plus_p_tot now INCLUDES this fluid (it flows through the generic
  // loop), so subtract our own background rho+p here.
  const double rho_plus_p_tot_rest = ppw->rho_plus_p_tot - Rho(ppw->pvecback.data()) * (1. + w_fld);

  double s2sq               = ppw->s_l[2] * ppw->s_l[2];
  double c_gamma_k_H_square = pow(c_gamma_over_c_fld_ * k / a_prime_over_a, 2) * cs2_fld_;
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
  if (ppw->scalar_ctx.gauge == static_cast<int>(possible_gauges::synchronous)) {
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
  // S quantity sourcing Gamma_prime evolution in the PPF scheme (eq. 15 of
  // 0808.3125); a ComputePpf-internal temporary, not stored on the workspace.
  const double S_fld = Rho(ppw->pvecback.data()) * (1. + w_fld) * 1.5 * a2 / k2 / a_prime_over_a *
                       (ppw->rho_plus_p_theta / rho_plus_p_tot_rest + k2 * alpha);
  // note that the last terms in the ratio do not include fld, that's correct, it's the whole point of the PPF scheme
  /** We must now check the stiffenss criterion again and set Gamma_prime_fld accordingly. */
  if (c_gamma_k_H_square > ppr->c_gamma_k_H_square_max) {
    ppw->Gamma_prime_fld = 0.;
  }
  else {
    ppw->Gamma_prime_fld = a_prime_over_a * (S_fld / (1. + c_gamma_k_H_square) -
                                             (1. + c_gamma_k_H_square) * Gamma_fld);
  }
  double Gamma_prime_plus_a_prime_over_a_Gamma = ppw->Gamma_prime_fld + a_prime_over_a * Gamma_fld;
  // delta and theta in both gauges gauge:
  ppw->rho_plus_p_theta_fld = Rho(ppw->pvecback.data()) * (1. + w_fld) * ppw->rho_plus_p_theta /
                                  rho_plus_p_tot_rest -
                              k2 * 2. / 3. * a_prime_over_a / a2 /
                                  (1 + 4.5 * a2 / k2 / s2sq * rho_plus_p_tot_rest) *
                                  (S_fld - Gamma_prime_plus_a_prime_over_a_Gamma / a_prime_over_a);
  ppw->delta_rho_fld        = -2. / 3. * k2 * s2sq / a2 * Gamma_fld -
                              3 * a_prime_over_a / k2 * ppw->rho_plus_p_theta_fld;

  /** Now construct the pressure perturbation, see 1903.xxxxx. */
  /** Construct energy density and pressure for DE (_fld) and the rest (_t).
      Also compute derivatives. */
  double rho_fld       = Rho(ppw->pvecback.data());
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
  double theta_t       = ppw->rho_plus_p_theta / rho_plus_p_tot_rest;
  double theta_t_prime = -a_prime_over_a * theta_t -
                         (p_t_prime * theta_t - k2 * ppw->delta_p + k2 * ppw->rho_plus_p_shear) /
                             rho_plus_p_tot_rest +
                         metric_euler;
  double S             = S_fld;
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
