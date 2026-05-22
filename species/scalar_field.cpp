#include "scalar_field.h"

#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

double ScalarFieldSpecies::V_scf(double phi) const {
  double lambda = pba_.scf_parameters[0];
  double alpha  = pba_.scf_parameters[1];
  double A      = pba_.scf_parameters[2];
  double B      = pba_.scf_parameters[3];
  return exp(-lambda * phi) * (pow(phi - B, alpha) + A);
}

double ScalarFieldSpecies::dV_scf(double phi) const {
  double lambda = pba_.scf_parameters[0];
  double alpha  = pba_.scf_parameters[1];
  double A      = pba_.scf_parameters[2];
  double B      = pba_.scf_parameters[3];
  double Ve     = exp(-lambda * phi);
  double Vp     = pow(phi - B, alpha) + A;
  double dVe    = -lambda * Ve;
  double dVp    = alpha * pow(phi - B, alpha - 1);
  return dVe * Vp + Ve * dVp;
}

double ScalarFieldSpecies::ddV_scf(double phi) const {
  double lambda = pba_.scf_parameters[0];
  double alpha  = pba_.scf_parameters[1];
  double A      = pba_.scf_parameters[2];
  double B      = pba_.scf_parameters[3];
  double Ve     = exp(-lambda * phi);
  double Vp     = pow(phi - B, alpha) + A;
  double dVe    = -lambda * Ve;
  double dVp    = alpha * pow(phi - B, alpha - 1);
  double ddVe   = lambda * lambda * Ve;
  double ddVp   = alpha * (alpha - 1.) * pow(phi - B, alpha - 2);
  return ddVe * Vp + 2 * dVe * dVp + Ve * ddVp;
}

void ScalarFieldSpecies::ComputeBackground(double a_rel,
                                           const double* pvecback_B,
                                           double* pvecback) {
  const double a                    = a_rel * pba_.a_today;
  const double phi                  = pvecback_B[index_bi_phi_scf_];
  const double phi_prime            = pvecback_B[index_bi_phi_prime_scf_];
  pvecback[index_bg_phi_scf_]       = phi;
  pvecback[index_bg_phi_prime_scf_] = phi_prime;
  pvecback[index_bg_V_scf_]         = V_scf(phi);
  pvecback[index_bg_dV_scf_]        = dV_scf(phi);
  pvecback[index_bg_ddV_scf_]       = ddV_scf(phi);
  pvecback[index_bg_rho_]           = (phi_prime * phi_prime / (2. * a * a) + V_scf(phi)) / 3.;
  pvecback[index_bg_p_]             = (phi_prime * phi_prime / (2. * a * a) - V_scf(phi)) / 3.;
  pvecback[index_bg_p_prime_scf_]   = 0.;
}

void ScalarFieldSpecies::BackgroundDerivs(double /*tau*/,
                                          const double* y,
                                          double* dy,
                                          const double* pvecback) {
  const double a         = pvecback[bgm_->index_bg_a_];
  const double H         = pvecback[bgm_->index_bg_H_];
  const double phi       = y[index_bi_phi_scf_];
  const double phi_prime = y[index_bi_phi_prime_scf_];
  /** phi'' + 2*a*H*phi' + a^2*dV = 0 */
  dy[index_bi_phi_scf_]       = phi_prime;
  dy[index_bi_phi_prime_scf_] = -a * (2. * H * phi_prime + a * dV_scf(phi));
}

double ScalarFieldSpecies::ComputePPrimeAndWrite(double a, double* pvecback) const {
  const double phi_prime          = pvecback[index_bg_phi_prime_scf_];
  const double H                  = pvecback[bgm_->index_bg_H_];
  const double dV                 = pvecback[index_bg_dV_scf_];
  const double p_prime            = phi_prime * (-phi_prime * H / a - 2. / 3. * dV);
  pvecback[index_bg_p_prime_scf_] = p_prime;
  return p_prime;
}

void ScalarFieldSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_scf", 0.);
  w.Add("(.)p_scf", 0.);
  w.Add("(.)p_prime_scf", 0.);
  w.Add("phi_scf", 0.);
  w.Add("phi'_scf", 0.);
  w.Add("V_scf", 0.);
  w.Add("V'_scf", 0.);
  w.Add("V''_scf", 0.);
}

void ScalarFieldSpecies::WriteBackgroundData(const double* pvecback,
                                             BackgroundColumnWriter& w) const {
  w.Add("(.)rho_scf", Rho(pvecback));
  w.Add("(.)p_scf", P(pvecback));
  w.Add("(.)p_prime_scf", pvecback[index_bg_p_prime_scf_]);
  w.Add("phi_scf", pvecback[index_bg_phi_scf_]);
  w.Add("phi'_scf", pvecback[index_bg_phi_prime_scf_]);
  w.Add("V_scf", pvecback[index_bg_V_scf_]);
  w.Add("V'_scf", pvecback[index_bg_dV_scf_]);
  w.Add("V''_scf", pvecback[index_bg_ddV_scf_]);
}

void ScalarFieldSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                     perturb_vector* /*pv*/,
                                                     const precision* /*ppr*/,
                                                     int& index_pt,
                                                     const perturb_workspace* /*ppw*/,
                                                     int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_phi = index_pt;
  ++index_pt;

  layout.idx_phi_prime = index_pt;
  ++index_pt;
}

void ScalarFieldSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                       double /*tau*/,
                                       const double* y,
                                       double* dy,
                                       const perturb_parameters_and_workspace& ppaw) {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  /* Use BackgroundModule's canonical pvecback slots for background quantities */
  auto bgm                  = ppaw.perturbations_module->GetBackgroundModule();
  const double phi_prime_bg = ppw->pvecback[bgm->index_bg_phi_prime_scf_];
  const double ddV_bg       = ppw->pvecback[bgm->index_bg_ddV_scf_];

  const double k2                = ctx.k2;
  const double a2                = ctx.a2;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double metric_continuity = ctx.metric_continuity;

  dy[layout.idx_phi]       = y[layout.idx_phi_prime];
  dy[layout.idx_phi_prime] = -2. * a_prime_over_a * y[layout.idx_phi_prime] -
                             metric_continuity * phi_prime_bg -
                             (k2 + a2 * ddV_bg) * y[layout.idx_phi];
}

void ScalarFieldSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                     const double* y,
                                     const double* /*dy*/,
                                     PerturbSourceContext& ctx) {
  const auto& layout          = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod  = ctx.p_mod;
  perturb_workspace* ppw      = ctx.ppw;
  const BackgroundModule* bgm = p_mod->GetBackgroundModule().get();
  const double* pvecback      = ppw->pvecback;

  // Scalar field sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const double phi_prime_bg = pvecback[bgm->index_bg_phi_prime_scf_];
  const double dV_bg        = pvecback[bgm->index_bg_dV_scf_];
  const double rho_scf      = Rho(pvecback);
  const double p_scf        = P(pvecback);
  const double a2_rel       = ctx.a2_rel;
  const double k2 = ctx.k *
                    ctx.k;  // PerturbSourceContext has no k2 field (unlike PerturbScalarContext)

  // ── delta_scf ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_scf_ == _TRUE_) {
    double delta_rho_scf;
    const perturbs* ppt = p_mod->GetPerturbs();
    if (ppt->gauge == synchronous) {
      delta_rho_scf = 1. / 3. *
                          (1. / a2_rel * phi_prime_bg * y[layout.idx_phi_prime] +
                           dV_bg * y[layout.idx_phi]) +
                      3. * ctx.a_prime_over_a * (1. + p_scf / rho_scf) *
                          ctx.theta_over_k2;  // N-body gauge correction
    }
    else {
      delta_rho_scf =
          1. / 3. *
              (1. / a2_rel * phi_prime_bg * y[layout.idx_phi_prime] + dV_bg * y[layout.idx_phi] -
               1. / a2_rel * phi_prime_bg * phi_prime_bg * ppw->pvecmetric[ppw->index_mt_psi]) +
          3. * ctx.a_prime_over_a * (1. + p_scf / rho_scf) *
              ctx.theta_over_k2;  // N-body gauge correction
    }
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_scf_,
                          ctx.index_tau,
                          ctx.index_k,
                          delta_rho_scf / rho_scf);
  }

  // ── theta_scf ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_scf_ == _TRUE_) {
    const double rho_plus_p_theta_scf = 1. / 3. * k2 / a2_rel * phi_prime_bg * y[layout.idx_phi];

    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_scf_,
                          ctx.index_tau,
                          ctx.index_k,
                          rho_plus_p_theta_scf / (rho_scf + p_scf) +
                              ctx.theta_shift);  // N-body gauge correction
  }
}

void ScalarFieldSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                double* y,
                                                const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_phi >= 0)
    y[layout.idx_phi] = 0.;
  if (layout.idx_phi_prime >= 0)
    y[layout.idx_phi_prime] = 0.;
}

void ScalarFieldSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                            const PerturbationsModule& mod,
                                            enum file_format fmt,
                                            BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_scf", mod.index_tp_delta_scf_, pba->has_scf);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t__scf", mod.index_tp_theta_scf_, pba->has_scf);
  }
}

void ScalarFieldSpecies::PrintVariables(PerturbColumnWriter& w,
                                        double /*tau*/,
                                        const double* y,
                                        const PerturbationsModule& mod,
                                        const perturb_workspace* ppw) const {
  const background* pba = mod.GetBackground();
  if (pba->has_scf != _TRUE_)
    return;

  double delta_scf = 0., theta_scf = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;
    const double k           = ppw->scalar_ctx.k;
    const double H           = pvecback[mod.GetBackgroundModule()->index_bg_H_];
    const double a           = pvecback[mod.GetBackgroundModule()->index_bg_a_];
    const double a2          = ppw->scalar_ctx.a2;
    const perturbs* ppt      = mod.GetPerturbs();

    const double phi_prime_bg = pvecback[mod.GetBackgroundModule()->index_bg_phi_prime_scf_];
    const double dV_bg        = pvecback[mod.GetBackgroundModule()->index_bg_dV_scf_];
    const double rho_scf      = Rho(pvecback);
    const double p_scf        = P(pvecback);

    const auto& layout = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);

    double delta_rho_scf = 0.;
    if (ppt->gauge == synchronous) {
      delta_rho_scf =
          1. / 3. * (1. / a2 * phi_prime_bg * y[layout.idx_phi_prime] + dV_bg * y[layout.idx_phi]);
    }
    else {
      delta_rho_scf = 1. / 3. *
                      (1. / a2 * phi_prime_bg * y[layout.idx_phi_prime] +
                       dV_bg * y[layout.idx_phi] -
                       1. / a2 * phi_prime_bg * phi_prime_bg * pvecmetric[ppw->index_mt_psi]);
    }

    const double rho_plus_p_theta_scf = 1. / 3. * ppw->scalar_ctx.k2 / a2 * phi_prime_bg *
                                        y[layout.idx_phi];

    delta_scf = delta_rho_scf / rho_scf;
    theta_scf = rho_plus_p_theta_scf / (rho_scf + p_scf);

    if (ppt->gauge == synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double w_scf  = p_scf / rho_scf;
      delta_scf          += alpha * (-3. * H * (1. + w_scf));
      theta_scf          += k * k * alpha;
    }
  }

  w.Add("delta_scf", delta_scf, true);
  w.Add("theta_scf", theta_scf, true);
}

double ScalarFieldSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                                 const perturb_vector* pv,
                                 const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  const double rho   = pvecback[index_bg_rho_];
  if (rho == 0.)
    return 0.;
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double dV        = pvecback[index_bg_dV_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;

  double delta_rho = (1. / 3.) *
                     (1. / a2 * phi_prime * y[layout.idx_phi_prime] + dV * y[layout.idx_phi]);

  if (ppw->scalar_ctx.gauge == newtonian) {
    double psi  = y[pv->index_pt_phi] - 4.5 * (a2 / k2) * ppw->rho_plus_p_shear;
    delta_rho  -= (1. / 3.) * (1. / a2) * phi_prime * phi_prime * psi;
  }

  return delta_rho / rho;
}

double ScalarFieldSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                                 const perturb_vector* pv,
                                 const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  const auto& layout      = static_cast<const PerturbLayout&>(base);
  const double rho        = pvecback[index_bg_rho_];
  const double p          = pvecback[index_bg_p_];
  const double rho_plus_p = rho + p;
  if (rho_plus_p == 0.)
    return 0.;
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;
  return (1. / 3.) * k2 / a2 * phi_prime * y[layout.idx_phi] / rho_plus_p;
}

double ScalarFieldSpecies::DeltaP(const BaseSpecies::PerturbLayout& base,
                                  const perturb_vector* /*pv*/,
                                  const double* y,
                                  const double* pvecback,
                                  const perturb_workspace* ppw) const {
  const auto& layout     = static_cast<const PerturbLayout&>(base);
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double dV        = pvecback[index_bg_dV_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;

  double delta_p = (1. / 3.) *
                   (1. / a2 * phi_prime * y[layout.idx_phi_prime] - dV * y[layout.idx_phi]);

  if (ppw->scalar_ctx.gauge == newtonian) {
    double psi  = y[ppw->pv->index_pt_phi] - 4.5 * (a2 / k2) * ppw->rho_plus_p_shear;
    delta_p    -= (1. / 3.) * (1. / a2) * phi_prime * phi_prime * psi;
  }

  return delta_p;
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> ScalarFieldSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_scf == _TRUE_) {
    result.push_back({"ScalarField", std::make_unique<ScalarFieldSpecies>(*ctx.pba)});
  }
  return result;
}
