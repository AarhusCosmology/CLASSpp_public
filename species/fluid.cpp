#include "fluid.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "axion_ede_fluid.h"
#include "background.h"
#include "background_module.h"
#include "errors.h"
#include "idm_dr_idr_species.h"
#include "idm_drmd_idr_drmd_species.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "ppf_fluid.h"
#include "species_collection.h"

FluidSpecies::FluidSpecies(const background& pba,
                           double omega0_fld,
                           equation_of_state fluid_eos,
                           double w0_fld,
                           double wa_fld,
                           double cs2_fld,
                           double Omega_EDE)
    : BaseSpecies("Fluid", EnergyType::DarkEnergy), pba_(pba), Omega0_fld_(omega0_fld),
      fluid_eos_(fluid_eos), w0_fld_(w0_fld), wa_fld_(wa_fld), cs2_fld_(cs2_fld),
      Omega_EDE_(Omega_EDE) {}

void FluidSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_fld_, true, index_bg, 1);
  index_bg_rho_ = index_bg_rho_fld_;
  class_define_index(index_bg_w_fld_, true, index_bg, 1);
  class_define_index(index_bg_dw_over_da_fld_, true, index_bg, 1);
}

void FluidSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_rho_fld_, true, index_bi, 1);
}

void FluidSpecies::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void FluidSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  /* rho_fld today */
  const double rho_fld_today = GetOmega0() * pow(pba_.H0, 2);

  /* integrate rho_fld(a) from a_ini to a_0, to get rho_fld(a_ini) given rho_fld(a0). */
  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(ctx.a_ini, &w_fld, &dw_over_da_fld, &integral_fld);

  /* Note: for complicated w_fld(a) functions with no simple
     analytic integral, this is the place were you should compute
     numerically the simple 1d integral [int_{a_ini}^{a_0} 3
     [(1+w_fld)/a] da] (e.g. with the Romberg method?) instead of
     calling ComputeWFld */

  /* rho_fld at initial time */
  ctx.pvecback_integration[bi_rho_index()] = rho_fld_today * exp(integral_fld);
}

void FluidSpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  // Compute w(a) ourselves rather than having the module pre-fill the slots.
  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld);
  pvecback[index_bg_w_fld_]          = w_fld;
  pvecback[index_bg_dw_over_da_fld_] = dw_over_da_fld;
  pvecback[index_bg_rho_fld_]        = pvecback_B[index_bi_rho_fld_];
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

double FluidSpecies::PPrime(double a,
                            double H,
                            const double* /*pvecback_B*/,
                            const double* pvecback) const {
  const double w_fld          = pvecback[index_bg_w_fld_];
  const double dw_over_da_fld = pvecback[index_bg_dw_over_da_fld_];
  return a * H * (a * dw_over_da_fld - 3. * (1. + w_fld) * w_fld) * pvecback[index_bg_rho_fld_];
}

void FluidSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_fld", 0.);
  w.Add("(.)w_fld", 0.);
}

void FluidSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_fld", pvecback[index_bg_rho_fld_]);
  w.Add("(.)w_fld", pvecback[index_bg_w_fld_]);
}

void FluidSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                               perturb_vector* /*pv*/,
                                               const precision* /*ppr*/,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  class_define_index(layout.idx_delta, true, index_pt, 1);
  class_define_index(layout.idx_theta, true, index_pt, 1);
}

void FluidSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                 double /*tau*/,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double a                 = ctx.a;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double k2                = ctx.k2;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;

  const double w_fld          = ppw->pvecback[index_bg_w_fld_];
  const double dw_over_da_fld = ppw->pvecback[index_bg_dw_over_da_fld_];
  const double w_prime_fld    = dw_over_da_fld * a_prime_over_a * a;
  const double ca2            = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;
  const double cs2            = Cs2(k2, a);

  dy[layout.idx_delta] = -(1. + w_fld) * (y[layout.idx_theta] + metric_continuity) -
                         3. * (cs2 - w_fld) * a_prime_over_a * y[layout.idx_delta] -
                         9. * (1. + w_fld) * (cs2 - ca2) * a_prime_over_a * a_prime_over_a *
                             y[layout.idx_theta] / k2;

  dy[layout.idx_theta] = -(1. - 3. * cs2) * a_prime_over_a * y[layout.idx_theta] +
                         cs2 * k2 / (1. + w_fld) * y[layout.idx_delta] + metric_euler;
}

void FluidSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                               const double* y,
                               const double* /*dy*/,
                               PerturbSourceContext& ctx) const {
  const auto& layout         = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const double* pvecback     = ppw->pvecback.data();

  // Fluid sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_fld ─────────────────────────────────────────────────────────────
  if (index_tp_delta_ >= 0) {
    const double w_fld     = W(pvecback);
    const double delta_fld = y[layout.idx_delta];
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_delta_,
                          ctx.index_tau,
                          ctx.index_k,
                          delta_fld + 3. * ctx.a_prime_over_a * (1. + w_fld) *
                                          ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_fld ─────────────────────────────────────────────────────────────
  if (index_tp_theta_ >= 0) {
    const double theta_fld = y[layout.idx_theta];
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_theta_,
                          ctx.index_tau,
                          ctx.index_k,
                          theta_fld + ctx.theta_shift);  // N-body gauge correction
  }
}

void FluidSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                          double* y,
                                          const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return;

  double w_fld, dw_over_da, integral;
  auto* bgm = ctx.p_mod->GetBackgroundModule().get();
  bgm->background_w_fld(ctx.a, &w_fld, &dw_over_da, &integral);
  const double cs2 = Cs2(ctx.k * ctx.k, ctx.a);

  y[layout.idx_delta] = -ctx.ktau_two / 4. * (1. + w_fld) * (4. - 3. * cs2) /
                        (4. - 6. * w_fld + 3. * cs2) * ctx.ppr->curvature_ini * ctx.s2_squared;

  y[layout.idx_theta] = -ctx.k * ctx.ktau_three / 4. * cs2 / (4. - 6. * w_fld + 3. * cs2) *
                        ctx.ppr->curvature_ini * ctx.s2_squared;
}

void FluidSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                      const PerturbationsModule& mod,
                                      file_format fmt,
                                      BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
      w.Add("d_fld", index_tp_delta_, index_tp_delta_ >= 0);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
      w.Add("t_fld", index_tp_theta_, index_tp_theta_ >= 0);
  }
}

void FluidSpecies::PrintVariables(PerturbColumnWriter& w,
                                  const BaseSpecies::PerturbLayout* base,
                                  double /*tau*/,
                                  const double* y,
                                  const PerturbationsModule& mod,
                                  const perturb_workspace* ppw) const {
  double delta_rho_fld = 0., rho_plus_p_theta_fld = 0., delta_p_fld = 0.;

  if (!w.IsTitleMode()) {
    const auto& layout   = static_cast<const PerturbLayout&>(*base);
    const auto se        = StressEnergy(layout, ppw->pv.get(), y, ppw->pvecback.data(), ppw);
    delta_rho_fld        = se.delta_rho;
    rho_plus_p_theta_fld = se.rho_plus_p_theta;
    delta_p_fld          = se.delta_p;
  }

  w.Add("delta_rho_fld", delta_rho_fld, true);
  w.Add("rho_plus_p_theta_fld", rho_plus_p_theta_fld, true);
  w.Add("delta_p_fld", delta_p_fld, true);
}

BaseSpecies::StressEnergyContribution FluidSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;

  // rho, p
  se.rho = pvecback[index_bg_rho_fld_];
  se.p   = pvecback[index_bg_w_fld_] * pvecback[index_bg_rho_fld_];

  // delta_rho: body of DeltaRho
  se.delta_rho = (layout.idx_delta >= 0) ? pvecback[index_bg_rho_fld_] * y[layout.idx_delta] : 0.;

  // rho_plus_p_theta: body of RhoPlusPTheta
  se.rho_plus_p_theta = (layout.idx_theta >= 0)
                            ? (pvecback[index_bg_rho_fld_] +
                               pvecback[index_bg_w_fld_] * pvecback[index_bg_rho_fld_]) *
                                  y[layout.idx_theta]
                            : 0.;

  // delta_p: body of DeltaP (PPF uses a dedicated path; this branch is only reached
  // for non-PPF fluid where idx_delta and idx_theta are registered).
  if (layout.idx_delta < 0 || layout.idx_theta < 0) {
    se.delta_p = 0.;
  }
  else {
    const double k2             = ppw->scalar_ctx.k2;
    const double a              = ppw->scalar_ctx.a;
    const double a_prime_over_a = pvecback[bgm_->index_bg_H_] * a;
    const double cs2            = Cs2(k2, a);

    double w_fld, dw_over_da_fld, integral_fld;
    ComputeWFld(a, &w_fld, &dw_over_da_fld, &integral_fld);
    const double w_prime_fld = dw_over_da_fld * a_prime_over_a * a;

    const double rho                  = pvecback[index_bg_rho_fld_];
    const double delta_rho_fld        = rho * y[layout.idx_delta];
    const double rho_plus_p_theta_fld = (1. + w_fld) * rho * y[layout.idx_theta];
    const double ca2_fld              = w_fld - w_prime_fld / 3. / (1. + w_fld) / a_prime_over_a;

    se.delta_p = cs2 * delta_rho_fld +
                 (cs2 - ca2_fld) * (3. * a_prime_over_a * rho_plus_p_theta_fld / k2);
  }

  // rho_plus_p_shear: always 0 for fluid
  se.rho_plus_p_shear = 0.;

  return se;
}

void FluidSpecies::ComputeWFld(double a,
                               double* w_fld,
                               double* dw_over_da_fld,
                               double* integral_fld) const {
  double Omega_ede          = 0.;
  double dOmega_ede_over_da = 0.;
  double a_eq               = 0.0;

  /** - first, define the function w(a) */
  switch (fluid_eos_) {
    case CLP:
      *w_fld = w0_fld_ + wa_fld_ * (1. - a);
      break;
    case EDE: {
      // Omega_ede(a) taken from eq. (10) in 1706.00730
      Omega_ede = (Omega0_fld_ - Omega_EDE_ * (1. - pow(a, -3. * w0_fld_))) /
                      (Omega0_fld_ + (1. - Omega0_fld_) * pow(a, 3. * w0_fld_)) +
                  Omega_EDE_ * (1. - pow(a, -3. * w0_fld_));

      // d Omega_ede / d a taken analytically from the above
      dOmega_ede_over_da = -Omega_EDE_ * 3. * w0_fld_ * pow(a, -3. * w0_fld_ - 1.) /
                               (Omega0_fld_ + (1. - Omega0_fld_) * pow(a, 3. * w0_fld_)) -
                           (Omega0_fld_ - Omega_EDE_ * (1. - pow(a, -3. * w0_fld_))) *
                               (1. - Omega0_fld_) * 3. * w0_fld_ * pow(a, 3. * w0_fld_ - 1.) /
                               pow(Omega0_fld_ + (1. - Omega0_fld_) * pow(a, 3. * w0_fld_), 2) +
                           Omega_EDE_ * 3. * w0_fld_ * pow(a, -3. * w0_fld_ - 1.);

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
        class_stop(
            "Early Dark Energy not compatible with decaying Dark Matter because we omitted "
            "to code the calculation of a_eq in that case, but it would not be difficult to "
            "add it if necessary, should be a matter of 5 minutes");
      a_eq = Omega_r / Omega_m;  // assumes a flat universe with a=1 today

      // w_ede(a) taken from eq. (11) in 1706.00730
      *w_fld = -dOmega_ede_over_da * a / Omega_ede / 3. / (1. - Omega_ede) + a_eq / 3. / (a + a_eq);
      break;
    }
    case PhenoAxion:
      class_stop(
          "internal error: base FluidSpecies::ComputeWFld reached with PhenoAxion eos "
          "(AxionEDEFluid must override)");
      break;
  }

  /** - then, give the corresponding analytic derivative dw/da (used
      by perturbation equations; we could compute it numerically,
      but with a loss of precision; as long as there is a simple
      analytic expression of the derivative of the previous
      function, let's use it! */
  switch (fluid_eos_) {
    case CLP:
      *dw_over_da_fld = -wa_fld_;
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
    case PhenoAxion:
      class_stop(
          "internal error: base FluidSpecies::ComputeWFld reached with PhenoAxion eos "
          "(AxionEDEFluid must override)");
      break;
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
  switch (fluid_eos_) {
    case CLP:
      *integral_fld = 3. * ((1. + w0_fld_ + wa_fld_) * log(1. / a) + wa_fld_ * (a - 1.));
      break;
    case EDE:
      class_stop(
          "EDE implementation not finished: to finish it, read the comments in background.c "
          "just before this line\n");
      break;
    case PhenoAxion:
      class_stop(
          "internal error: base FluidSpecies::ComputeWFld reached with PhenoAxion eos "
          "(AxionEDEFluid must override)");
      break;
  }

  /** note: of course you can generalise these formulas to anything,
      defining new parameters pba->w..._fld. Just remember that so
      far, HyRec explicitely assumes that w(a)= w0 + wa (1-a/a0); but
      Recfast does not assume anything */
}

bool FluidSpecies::ReachesPhantomDivide() const {
  double w_ini, w_0, dw_over_da, integral;
  ComputeWFld(0., &w_ini, &dw_over_da, &integral);
  ComputeWFld(1., &w_0, &dw_over_da, &integral);
  return (w_ini + 1.) * (w_0 + 1.) <= 0.;
}

bool FluidSpecies::HyrecCplApproximation(double* w0, double* wa) const {
  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(1., &w_fld, &dw_over_da_fld, &integral_fld);
  *w0 = w_fld;
  *wa = -dw_over_da_fld;
  return true;
}

// ── Newtonian-gauge transform ─────────────────────────────────────────────────

void FluidSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                 double* y,
                                                 const PerturbIcContext& ctx) {
  // Non-PPF only: in PPF mode idx_delta/idx_theta are unregistered (== -1) and
  // the helper is a no-op. delta uses the universal ρ̇/ρ shift = -3(1+w)ℋα, which
  // corrects the historical opposite-sign bug.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

void FluidSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                 const BaseSpecies::PerturbLayout& new_base,
                                                 const double* old_y,
                                                 double* new_y,
                                                 const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  // Vector/tensor pv: all -1 -> no-op.
  if (old_l.idx_delta >= 0 && new_l.idx_delta >= 0)
    new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  if (old_l.idx_theta >= 0 && new_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> FluidSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // ── fluid_equation_of_state (string-keyed enum) ────────────────────────────
  equation_of_state fluid_eos = CLP;
  if (auto eos_opt = ctx.pfc->get<std::string>("fluid_equation_of_state")) {
    const std::string& eos_str = *eos_opt;
    if (eos_str.find("pheno_axion") != std::string::npos ||
        eos_str.find("PhenoAxion") != std::string::npos) {
      return CreatePhenoAxion(ctx);
    }
    else if (eos_str.find("CLP") != std::string::npos || eos_str.find("clp") != std::string::npos) {
      fluid_eos = CLP;
    }
    else if (eos_str.find("EDE") != std::string::npos || eos_str.find("ede") != std::string::npos) {
      fluid_eos = EDE;
    }
    else {
      class_stop_severe("incomprehensible input '%s' for the field 'fluid_equation_of_state'",
                        eos_str.c_str());
    }
  }

  // Resolve Omega0_fld: either the closure-budget override or pfc.
  double omega0_fld = 0.;
  if (ctx.omega0_closure_override.has_value()) {
    omega0_fld = *ctx.omega0_closure_override;
  }
  else {
    omega0_fld = ctx.pfc->get_or("Omega_fld", omega0_fld);
  }
  if (omega0_fld == 0.)
    return result;

  // ── numeric params ────────────────────────────────────────────────────────
  double w0_fld    = -1.;
  double wa_fld    = 0.;
  double cs2_fld   = 1.;
  double Omega_EDE = 0.;
  if (fluid_eos == CLP) {
    w0_fld  = ctx.pfc->get_or("w0_fld", w0_fld);
    wa_fld  = ctx.pfc->get_or("wa_fld", wa_fld);
    cs2_fld = ctx.pfc->get_or("cs2_fld", cs2_fld);
  }
  else {  // EDE
    w0_fld    = ctx.pfc->get_or("w0_fld", w0_fld);
    Omega_EDE = ctx.pfc->get_or("Omega_EDE", Omega_EDE);
    cs2_fld   = ctx.pfc->get_or("cs2_fld", cs2_fld);
  }

  // ── PPF flag + sound-speed param ──────────────────────────────────────────
  bool use_ppf              = true;
  double c_gamma_over_c_fld = 0.4;
  if (auto ppf_opt = ctx.pfc->get<std::string>("use_ppf")) {
    const std::string& ppf_str = *ppf_opt;
    use_ppf = (ppf_str.find("y") != std::string::npos || ppf_str.find("Y") != std::string::npos)
                  ? true
                  : false;
  }
  if (use_ppf)
    c_gamma_over_c_fld = ctx.pfc->get_or("c_gamma_over_c_fld", c_gamma_over_c_fld);

  std::unique_ptr<BaseSpecies> sp;
  if (use_ppf) {
    sp = std::make_unique<PpfFluid>(*ctx.pba,
                                    omega0_fld,
                                    fluid_eos,
                                    w0_fld,
                                    wa_fld,
                                    cs2_fld,
                                    Omega_EDE,
                                    c_gamma_over_c_fld);
  }
  else {
    sp = std::make_unique<FluidSpecies>(*ctx.pba,
                                        omega0_fld,
                                        fluid_eos,
                                        w0_fld,
                                        wa_fld,
                                        cs2_fld,
                                        Omega_EDE);
  }
  result.push_back({"Fluid", std::move(sp)});
  return result;
}

std::vector<Named> FluidSpecies::CreatePhenoAxion(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  class_test_severe(ctx.omega0_closure_override.has_value(),
                    "the pheno_axion fluid cannot be the budget-closure species (its density is "
                    "set by a_c and the EDE fraction); let Lambda close the budget");

  // ── n XOR w_fld_f ──────────────────────────────────────────────────────────
  const auto n_opt  = ctx.pfc->get<double>("n_pheno_axion");
  const auto wf_opt = ctx.pfc->get<double>("w_fld_f");
  class_test_severe(n_opt && wf_opt, "give only one of 'n_pheno_axion' and 'w_fld_f'");
  class_test_severe(!n_opt && !wf_opt, "pheno_axion needs 'n_pheno_axion' (or 'w_fld_f')");
  double n, w_f;
  if (n_opt) {
    n = *n_opt;
    class_test(n < 1., "n_pheno_axion must be >= 1");
    w_f = AxionEDEFluid::WFinal(n);
  }
  else {
    w_f = *wf_opt;
    // w_f < 0 would give a derived n = (1+w_f)/(1-w_f) < 1, violating the model
    // requirement n >= 1 (and making the k->0 sound speed w_f negative).
    class_test(w_f < 0. || w_f >= 1.,
               "w_fld_f must lie in [0, 1) so that the derived n = (1+w_f)/(1-w_f) is >= 1");
    n = (1. + w_f) / (1. - w_f);
  }
  const double w_i = ctx.pfc->get_or("w_fld_i", -1.);
  const double nu  = ctx.pfc->get_or("nu_fld", 1.);
  class_test(nu <= 0., "nu_fld must be > 0");
  // w_i < -1 would make w(a) cross -1 on its way to w_f, which the true-fluid
  // equations cannot handle (AxionEDEFluid::ReachesPhantomDivide answers false
  // on the premise w(a) > -1 for all a > 0; w_i = -1 exactly is the designed
  // asymptote and stays safe because w(a) > w_i at every a > 0).
  class_test(w_i < -1., "w_fld_i must be >= -1 (w may not cross the phantom divide)");
  class_test(w_f <= w_i, "pheno_axion needs w_fld_f > w_fld_i");

  // ── a_c XOR log10_axion_ac ─────────────────────────────────────────────────
  const auto ac_opt     = ctx.pfc->get<double>("a_c");
  const auto log_ac_opt = ctx.pfc->get<double>("log10_axion_ac");
  class_test_severe(ac_opt && log_ac_opt, "give only one of 'a_c' and 'log10_axion_ac'");
  class_test_severe(!ac_opt && !log_ac_opt, "pheno_axion needs 'a_c' or 'log10_axion_ac'");
  const double a_c = ac_opt ? *ac_opt : std::pow(10., *log_ac_opt);
  class_test(a_c <= 0. || a_c >= 1., "a_c must lie in (0, 1)");

  // ── Theta_initial_fld (required; enters only the omega_axion calibration) ──
  const auto theta_opt = ctx.pfc->get<double>("Theta_initial_fld");
  class_test_severe(!theta_opt, "pheno_axion requires 'Theta_initial_fld'");
  class_test(*theta_opt <= 0. || *theta_opt >= _PI_,
             "pheno_axion requires Theta_initial_fld in (0, pi)");

  // ── density: exactly one of Omega_fld / Omega_fld_ac / fraction_fld_ac ────
  const auto om0_opt  = ctx.pfc->get<double>("Omega_fld");
  const auto omac_opt = ctx.pfc->get<double>("Omega_fld_ac");
  const auto frac_opt = ctx.pfc->get<double>("fraction_fld_ac");
  const int n_density = int(om0_opt.has_value()) + int(omac_opt.has_value()) +
                        int(frac_opt.has_value());
  class_test_severe(n_density != 1,
                    "pheno_axion needs exactly one of 'Omega_fld', 'Omega_fld_ac', "
                    "'fraction_fld_ac'");

  double omega0_fld = 0.;
  if (om0_opt) {
    omega0_fld = *om0_opt;
  }
  else if (omac_opt) {
    omega0_fld = AxionEDEFluid::OmegaZeroFromOmegaAc(*omac_opt, a_c, nu, w_i, w_f);
  }
  else {
    // fraction f = rho_fld(a_c)/rho_tot(a_c): Omega_fld_ac = Omega_tot_others(a_c) f/(1-f),
    // with the other species scaled to a_c by EnergyType (AxiCLASS input.c:4157-4174;
    // Lambda is excluded there too — negligible at EDE-era a_c).
    const double f_ac = *frac_opt;
    class_test(f_ac <= 0. || f_ac >= 1., "fraction_fld_ac must lie in (0, 1)");
    // Programmer-error invariant: the caller wiring, not user input, controls whether
    // all_species is populated. This can only fire if ConstructSpecies fails to pass it.
    if (ctx.all_species == nullptr)
      throw std::logic_error(
          "fraction_fld_ac needs the species collection (internal: all_species missing)");
    double omega_tot_ac = 0.;
    for (const auto& sp : *ctx.all_species) {
      switch (sp->energy_type()) {
        case BaseSpecies::EnergyType::Radiation:
          omega_tot_ac += sp->GetOmega0() * std::pow(a_c, -4.);
          break;
        case BaseSpecies::EnergyType::Matter:
          omega_tot_ac += sp->GetOmega0() * std::pow(a_c, -3.);
          break;
        default:
          break;
      }
    }
    const double omega_ac = omega_tot_ac * f_ac / (1. - f_ac);
    omega0_fld            = AxionEDEFluid::OmegaZeroFromOmegaAc(omega_ac, a_c, nu, w_i, w_f);
  }
  class_test(omega0_fld <= 0.,
             "pheno_axion needs a positive fluid density: Omega_fld / Omega_fld_ac must be > 0");

  // ── PPF is meaningless here; cs2 is derived, not an input ──────────────────
  if (auto ppf_opt = ctx.pfc->get<std::string>("use_ppf")) {
    class_test_severe(ppf_opt->find("y") != std::string::npos ||
                          ppf_opt->find("Y") != std::string::npos,
                      "use_ppf = yes is incompatible with pheno_axion (w > -1 always; the GDM "
                      "sound speed requires the true fluid equations)");
  }
  class_test_severe(ctx.pfc->get<double>("cs2_fld").has_value(),
                    "cs2_fld is not an input for pheno_axion (the sound speed is the GDM "
                    "formula derived from n, a_c and Theta_initial_fld)");

  result.push_back(
      {"Fluid",
       std::make_unique<AxionEDEFluid>(*ctx.pba, omega0_fld, a_c, n, nu, w_i, w_f, *theta_opt)});
  return result;
}
