#include "dark_radiation_species.h"

#include <cmath>

#include "background_module.h"
#include "perturbations_module.h"

// ── Background ─────────────────────────────────────────────────────────────

void DarkRadiationSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_rho_ = index_bg++;
}

void DarkRadiationSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_rho_ = index_bi++;
}

void DarkRadiationSpecies::ComputeBackground(double /*a*/,
                                             const double* pvecback_B,
                                             double* pvecback) {
  pvecback[index_bg_rho_] = pvecback_B[index_bi_rho_];
}

void DarkRadiationSpecies::BackgroundDerivs(double /*tau*/,
                                            const double* y,
                                            double* dy,
                                            const double* pvecback) {
  const double a = pvecback[bgm_->index_bg_a_];
  const double H = pvecback[bgm_->index_bg_H_];

  // Dilution only; decay source added by the parent composite's BackgroundDerivs.
  dy[index_bi_rho_] = -4. * a * H * y[index_bi_rho_];
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DarkRadiationSpecies::RegisterTransferSourceIndices(int& index_tp,
                                                         const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void DarkRadiationSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                       perturb_vector* /*pv*/,
                                                       const precision* ppr,
                                                       int& index_pt,
                                                       const perturb_workspace* /*ppw*/,
                                                       int /*gauge*/) {
  auto& layout   = static_cast<PerturbLayout&>(base);
  layout.idx_F0  = index_pt;
  layout.l_max   = ppr->l_max_dr;
  index_pt_F0_   = index_pt;
  index_pt      += ppr->l_max_dr + 1;
}

void DarkRadiationSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                         double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ctx.k;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;
  const double s2_squared         = ctx.s2_squared;
  const double a                  = ctx.a;
  const double* pvecback          = ppw->pvecback.data();
  const int base_idx              = layout.idx_F0;
  const int lmax                  = layout.l_max;
  double r_dr = pvecback[index_bg_rho_] * (a * a) * (a * a) / (pba_->H0 * pba_->H0);

  // l=0: free-streaming only (coupling source added by parent composite's AddCouplingDerivs)
  dy[base_idx + 0] = -k * y[base_idx + 1] - 4. / 3. * metric_continuity * r_dr;
  // l=1: free-streaming only
  dy[base_idx + 1] = k / 3. * y[base_idx + 0] - 2. / 3. * k * y[base_idx + 2] * s2_squared +
                     4. * metric_euler / (3. * k) * r_dr;
  // l=2
  dy[base_idx + 2] = 8. / 15. * (3. / 4. * k * y[base_idx + 1] + metric_shear * r_dr) -
                     3. / 5. * k * s_l[3] / s_l[2] * y[base_idx + 3];
  // l=3
  {
    int l            = 3;
    dy[base_idx + l] = k / (2. * l + 1.) *
                       (l * s_l[l] * s_l[2] * y[base_idx + l - 1] -
                        (l + 1.) * s_l[l + 1] * y[base_idx + l + 1]);
  }
  // l=4..lmax-1
  for (int l = 4; l < lmax; ++l)
    dy[base_idx + l] = k / (2. * l + 1.) *
                       (l * s_l[l] * y[base_idx + l - 1] -
                        (l + 1.) * s_l[l + 1] * y[base_idx + l + 1]);
  // l=lmax (truncation)
  {
    int l            = lmax;
    dy[base_idx + l] = k * (s_l[l] * y[base_idx + l - 1] - (1. + l) * cotKgen * y[base_idx + l]);
  }
}

/**
 * DR perturbation variables F_l follow the convention of astro-ph/9907388:
 *   delta_rho_dr = rho_dr * F0 / f,  where f = rho_dr * a^4 / H0^2
 * so delta_rho_dr = (H0/a^2)^2 * F0 = rho_dr_over_f * F0.
 */
BaseSpecies::StressEnergyContribution DarkRadiationSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_rho_] / 3.;

  if (layout.idx_F0 < 0)
    return se;
  const double a             = pvecback[bgm_->index_bg_a_];
  const double a2            = a * a;
  const double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  // δρ and (ρ+P)θ also vanish for ρ≤0 (legacy of their old /ρ form); δp and
  // (ρ+P)σ never divided by ρ, so they guard on idx_F0 only.
  const bool rho_pos = pvecback[index_bg_rho_] > 0.;

  se.delta_rho        = rho_pos ? rho_dr_over_f * y[layout.idx_F0] : 0.;
  se.rho_plus_p_theta = rho_pos ? ppw->scalar_ctx.k * rho_dr_over_f * y[layout.idx_F0 + 1] : 0.;
  se.delta_p          = rho_dr_over_f * y[layout.idx_F0] / 3.;
  se.rho_plus_p_shear = 2. / 3. * rho_dr_over_f * y[layout.idx_F0 + 2];

  return se;
}

void DarkRadiationSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                  double* y,
                                                  const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;

  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0)
    return;

  const double* pvecback   = ctx.ppw->pvecback.data();
  const double r_prefactor = std::pow(std::pow(ctx.a, 2) / pba_->H0, 2);
  const double r_dr        = r_prefactor * pvecback[index_bg_rho_];

  y[layout.idx_F0 + 0] = ctx.delta_dr * r_dr;
  if (layout.l_max >= 1)
    y[layout.idx_F0 + 1] = 4. / (3. * ctx.k) * ctx.theta_ur * r_dr;
  if (layout.l_max >= 2)
    y[layout.idx_F0 + 2] = 2. * ctx.shear_ur * r_dr;
  if (layout.l_max >= 3)
    y[layout.idx_F0 + 3] = ctx.l3_ur * r_dr;
}

void DarkRadiationSpecies::PerturbNewtonianReseed(const PerturbLayout& layout,
                                                  double* y,
                                                  const PerturbIcContext& ctx,
                                                  double decay_corr) const {
  if (layout.idx_F0 < 0)
    return;
  const double* pvecback = ctx.ppw->pvecback.data();
  const double r_dr      = std::pow(std::pow(ctx.a, 2) / pba_->H0, 2) * pvecback[index_bg_rho_];
  const double delta_dr  = ctx.delta_dr + (-4. * ctx.a_prime_over_a + decay_corr) * ctx.alpha;
  const double theta_ur  = ctx.theta_ur + ctx.k * ctx.k * ctx.alpha;
  y[layout.idx_F0 + 0]   = delta_dr * r_dr;
  if (layout.l_max >= 1)
    y[layout.idx_F0 + 1] = 4. / (3. * ctx.k) * theta_ur * r_dr;
  if (layout.l_max >= 2)
    y[layout.idx_F0 + 2] = 2. * ctx.shear_ur * r_dr;  // shear/l3 gauge-invariant
  if (layout.l_max >= 3)
    y[layout.idx_F0 + 3] = ctx.l3_ur * r_dr;
}

void DarkRadiationSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y,
    double* new_y,
    const PerturbSwitchContext& /*ctx*/) const {
  const auto& o = static_cast<const PerturbLayout&>(old_base);
  const auto& n = static_cast<const PerturbLayout&>(new_base);
  if (o.idx_F0 < 0 || n.idx_F0 < 0)
    return;
  for (int l = 0; l <= n.l_max; ++l)
    new_y[n.idx_F0 + l] = old_y[o.idx_F0 + l];
}
