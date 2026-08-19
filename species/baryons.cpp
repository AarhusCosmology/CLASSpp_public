#include "baryons.h"

#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

// ── RegisterPerturbationIndices ──────────────

void BaryonsSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                 perturb_vector* pv,
                                                 const precision* /*ppr*/,
                                                 int& index_pt,
                                                 const perturb_workspace* /*ppw*/,
                                                 int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = index_pt;
  ++index_pt;

  layout.idx_theta = index_pt;
  ++index_pt;
}

void BaryonsSpecies::RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                       perturb_vector* /*pv*/,
                                                       const precision* /*ppr*/,
                                                       int& index_pt,
                                                       const perturb_workspace* /*ppw*/,
                                                       int /*gauge*/) {
  auto& layout     = static_cast<PerturbLayout&>(base);
  layout.idx_theta = index_pt++; /* v_b^{(1)} */
}

// ── PerturbDerivs ────────────────────────────

void BaryonsSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                   double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double k2                = ppaw.k * ppaw.k;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;
  const double R                 = ctx.R;
  const double s2_squared        = ctx.s2_squared;

  const double theta_b = y[layout.idx_theta];
  /* theta_g from context: RSA-corrected photon velocity */
  const double theta_g = ctx.theta_g;
  const double delta_g = ctx.delta_g;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  /* baryon pressure perturbation / rho (same logic as StressEnergy delta_p branch) */
  double delta_p_b_over_rho_b;
  if (ppt_->has_perturbed_recombination && ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    delta_p_b_over_rho_b = ppw->pvecthermo[thm_->index_th_wb_] *
                           (y[layout.idx_delta] +
                            y[ppw->pv->index_pt_perturbed_recombination_delta_temp]);
  else
    delta_p_b_over_rho_b = ppw->pvecthermo[thm_->index_th_cb2_] * y[layout.idx_delta];

  /* density equation */
  dy[layout.idx_delta] = -(theta_b + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == (int) tca_off) {
    /* Full equation */
    dy[layout.idx_theta] = -a_prime_over_a * theta_b + metric_euler + k2 * delta_p_b_over_rho_b +
                           R * dkappa * (theta_g - theta_b);
  }
  else {
    /* TCA on: tight-coupling approximation for theta_b */
    dy[layout.idx_theta] =
        (-a_prime_over_a * theta_b +
         k2 * (delta_p_b_over_rho_b + R * (delta_g / 4. - s2_squared * ppw->tca_shear_g)) +
         R * ppw->tca_slip) /
            (1. + R) +
        metric_euler;
  }
}

// Direct fused override: one dispatch, no inner virtual calls. Reproduces each
// individual function bit-for-bit (baryons: P ≡ 0, ρ+P = ρ, no shear).
BaseSpecies::StressEnergyContribution BaryonsSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* pv,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const BaryonsSpecies::PerturbLayout&>(base);
  const double rho_b = pvecback[index_bg_rho_];
  StressEnergyContribution se;
  se.rho              = rho_b;
  se.delta_rho        = rho_b * y[layout.idx_delta];
  se.rho_plus_p_theta = rho_b * y[layout.idx_theta];
  if (ppt_->has_perturbed_recombination && ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    se.delta_p = rho_b * ppw->pvecthermo[thm_->index_th_wb_] *
                 (y[layout.idx_delta] + y[pv->index_pt_perturbed_recombination_delta_temp]);
  else
    se.delta_p = rho_b * ppw->pvecthermo[thm_->index_th_cb2_] * y[layout.idx_delta];
  return se;
}

// ── ApplyInitialConditions ───────────────────

void BaryonsSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                            double* y,
                                            const PerturbIcContext& ctx) {
  const auto& layout             = static_cast<const PerturbLayout&>(base);
  const PerturbationsModule* mod = ctx.p_mod;
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return;

  if (ctx.index_ic == mod->index_ic_ad_) {
    y[layout.idx_delta] = 3. / 4. * ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[layout.idx_delta] = 3. / 4. * ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_bi_) {
    y[layout.idx_delta] = ctx.ppr->entropy_ini + 3. / 4. * ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_nid_) {
    y[layout.idx_delta] = ctx.ppr->entropy_ini * ctx.fracnu / ctx.fracg / 8. * ctx.ktau_two;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_niv_) {
    y[layout.idx_delta] = 3. / 4. * ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
}

// ── FillSources ──────────────────────────────

void BaryonsSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                 const double* y,
                                 const double* /*dy*/,
                                 PerturbSourceContext& ctx) const {
  const auto& layout         = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod = ctx.p_mod;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const double a_prime_over_a = ctx.a_prime_over_a;

  // ── delta_b: baryon density transfer ─────────────────────────────────────
  if (index_tp_delta_ >= 0) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_delta_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[layout.idx_delta] +
                              3. * a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_b: baryon velocity transfer ────────────────────────────────────
  if (index_tp_theta_ >= 0) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_theta_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[layout.idx_theta] + ctx.theta_shift);  // N-body gauge correction
  }
}

void BaryonsSpecies::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void BaryonsSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                        const PerturbationsModule& mod,
                                        file_format fmt,
                                        BaseSpecies::TransferColumnSection section) const {
  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
      w.Add("d_b", index_tp_delta_, index_tp_delta_ >= 0);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
      w.Add("t_b", index_tp_theta_, index_tp_theta_ >= 0);
  }
}

void BaryonsSpecies::PrintVariables(PerturbColumnWriter& w,
                                    const BaseSpecies::PerturbLayout* base,
                                    double /*tau*/,
                                    const double* y,
                                    const PerturbationsModule& mod,
                                    const perturb_workspace* ppw) const {
  double delta_b = 0., theta_b = 0.;

  if (!w.IsTitleMode()) {
    const auto& layout       = static_cast<const PerturbLayout&>(*base);
    const double k           = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback.data();
    const double* pvecmetric = ppw->pvecmetric.data();

    delta_b = y[layout.idx_delta];
    theta_b = y[layout.idx_theta];

    // Converting synchronous variables to Newtonian
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == possible_gauges::synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_b            -= 3. * H * a * alpha;
      theta_b            += k * k * alpha;
    }
  }

  w.Add("delta_b", delta_b, true);
  w.Add("theta_b", theta_b, true);
}

// ── Newtonian-gauge transform ─────────────────────────────────────────────────

void BaryonsSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                   double* y,
                                                   const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

void BaryonsSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                   const BaseSpecies::PerturbLayout& new_base,
                                                   const double* old_y,
                                                   double* new_y,
                                                   const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  // Scalar-mode only: idx_delta is registered (>=0) in scalar pv and -1 in
  // vector/tensor pv, so this no-ops outside scalar mode (the vector-mode
  // baryon theta copy stays inline -- out of scope for this PR).
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> BaryonsSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  result.push_back({"Baryons", std::make_unique<BaryonsSpecies>(*ctx.pba)});
  return result;
}
