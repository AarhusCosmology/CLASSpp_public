#include "dark_radiation_species.h"

#include <cmath>

#include "background_module.h"
#include "dcdm.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

// ── Background ─────────────────────────────────────────────────────────────

void DarkRadiationSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_rho_dr_species_  = index_bg;
  index_bg                 += pba_->N_decay_dr;
  index_bg_rho_             = index_bg++;  // total DR density (from base class protected member)
}

void DarkRadiationSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_rho_dr_species_  = index_bi;
  index_bi                 += pba_->N_decay_dr;
}

void DarkRadiationSpecies::ComputeBackground(double /*a_rel*/,
                                             const double* pvecback_B,
                                             double* pvecback) {
  double rho_dr_total = 0.;
  for (int n = 0; n < pba_->N_decay_dr; ++n) {
    pvecback[index_bg_rho_dr_species_ + n]  = pvecback_B[index_bi_rho_dr_species_ + n];
    rho_dr_total                           += pvecback_B[index_bi_rho_dr_species_ + n];
  }
  pvecback[index_bg_rho_] = rho_dr_total;
}

void DarkRadiationSpecies::BackgroundDerivs(double /*tau*/,
                                            const double* y,
                                            double* dy,
                                            const double* pvecback) {
  const double a = pvecback[bgm_->index_bg_a_];
  const double H = pvecback[bgm_->index_bg_H_];

  for (int n = 0; n < pba_->N_decay_dr; ++n) {
    // Dilution only; DCDM decay source added by DCDM_DR_Species::BackgroundDerivs
    dy[index_bi_rho_dr_species_ + n] = -4. * a * H * y[index_bi_rho_dr_species_ + n];
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DarkRadiationSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                       perturb_vector* pv,
                                                       const precision* /*ppr*/,
                                                       int& index_pt,
                                                       const perturb_workspace* /*ppw*/,
                                                       int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  // DR sum multipoles: F0_dr_sum[l=0..l_max_dr]
  layout.idx_F0_sum       = index_pt;
  pv->index_pt_F0_dr_sum  = index_pt;
  index_pt               += pv->l_max_dr + 1;

  // Per-species multipoles: N_decay_dr * (l_max_dr+1) slots
  layout.idx_F0_species       = index_pt;
  pv->index_pt_F0_dr_species  = index_pt;
  index_pt_F0_dr_species_     = index_pt;
  index_pt                   += pba_->N_decay_dr * (pv->l_max_dr + 1);

  layout.l_max = pv->l_max_dr;
}

void DarkRadiationSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                         double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;
  const double k                  = ctx.k;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;
  const double s2_squared         = ctx.s2_squared;
  const double a                  = ctx.a;

  const double* pvecback = ppw->pvecback;

  const int idx_sum = layout.idx_F0_sum;
  const int idx_sp  = layout.idx_F0_species;
  const int l_max   = layout.l_max;

  // Zero out the running sum
  for (int l = 0; l <= l_max; ++l)
    dy[idx_sum + l] = 0.;

  if (pba_->N_decay_dr > 0) {
    double r_dr = pvecback[index_bg_rho_dr_species_] * (a * a) * (a * a) / (pba_->H0 * pba_->H0);

    // l=0: free-streaming only (coupling source from DCDM added by DCDM_DR_Species::AddCouplingDerivs)
    dy[idx_sp + 0] = -k * y[idx_sp + 1] - 4. / 3. * metric_continuity * r_dr;
    // l=1: free-streaming only
    dy[idx_sp + 1] = k / 3. * y[idx_sp + 0] - 2. / 3. * k * y[idx_sp + 2] * s2_squared +
                     4. * metric_euler / (3. * k) * r_dr;
    // l=2
    dy[idx_sp + 2] = 8. / 15. * (3. / 4. * k * y[idx_sp + 1] + metric_shear * r_dr) -
                     3. / 5. * k * s_l[3] / s_l[2] * y[idx_sp + 3];
    // l=3
    {
      int l          = 3;
      dy[idx_sp + l] = k / (2. * l + 1.) *
                       (l * s_l[l] * s_l[2] * y[idx_sp + l - 1] -
                        (l + 1.) * s_l[l + 1] * y[idx_sp + l + 1]);
    }
    // l=4..l_max-1
    for (int l = 4; l < l_max; ++l)
      dy[idx_sp + l] = k / (2. * l + 1.) *
                       (l * s_l[l] * y[idx_sp + l - 1] - (l + 1.) * s_l[l + 1] * y[idx_sp + l + 1]);
    // l=l_max (truncation)
    {
      int l          = l_max;
      dy[idx_sp + l] = k * (s_l[l] * y[idx_sp + l - 1] - (1. + l) * cotKgen * y[idx_sp + l]);
    }
    // Accumulate into sum
    for (int l = 0; l <= l_max; ++l)
      dy[idx_sum + l] += dy[idx_sp + l];
  }
}

/**
 * DR perturbation variables F_l follow the convention of astro-ph/9907388:
 *   delta_rho_dr = rho_dr * F0 / f,  where f = rho_dr * a^4 / H0^2
 * so delta_rho_dr = (H0/a^2)^2 * F0 = rho_dr_over_f * F0.
 * The factor rho_dr_over_f = H0^2/a^4 converts F_l to physical quantities.
 */
double DarkRadiationSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                                   const perturb_vector* /*pv*/,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0_sum < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[layout.idx_F0_sum] / pvecback[index_bg_rho_];
}

double DarkRadiationSpecies::Delta(const perturb_vector* pv,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* /*ppw*/) const {
  // Legacy path: pv->index_pt_F0_dr_sum is dual-written by RegisterPerturbationIndices
  if (pv->index_pt_F0_dr_sum < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[pv->index_pt_F0_dr_sum] / pvecback[index_bg_rho_];
}

double DarkRadiationSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                                   const perturb_vector* /*pv*/,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0_sum < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  double k             = ppw->scalar_ctx.k;
  // (rho+p)*theta = k * rho_dr_over_f * F1, and (rho+p) = (4/3)*rho_dr
  return 0.75 * k * rho_dr_over_f * y[layout.idx_F0_sum + 1] / pvecback[index_bg_rho_];
}

double DarkRadiationSpecies::Theta(const perturb_vector* pv,
                                   const double* y,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  // Legacy path: pv->index_pt_F0_dr_sum is dual-written by RegisterPerturbationIndices
  if (pv->index_pt_F0_dr_sum < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  double k             = ppw->scalar_ctx.k;
  return 0.75 * k * rho_dr_over_f * y[pv->index_pt_F0_dr_sum + 1] / pvecback[index_bg_rho_];
}

double DarkRadiationSpecies::DeltaP(const BaseSpecies::PerturbLayout& base,
                                    const perturb_vector* /*pv*/,
                                    const double* y,
                                    const double* pvecback,
                                    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0_sum < 0)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[layout.idx_F0_sum] / 3.;
}

double DarkRadiationSpecies::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                           const perturb_vector* /*pv*/,
                                           const double* y,
                                           const double* pvecback,
                                           const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0_sum < 0)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return 2. / 3. * rho_dr_over_f * y[layout.idx_F0_sum + 2];
}

void DarkRadiationSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                  double* y,
                                                  const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;

  const auto& layout       = static_cast<const PerturbLayout&>(base);
  const double* pvecback   = ctx.ppw->pvecback;
  const double r_prefactor = std::pow(std::pow(ctx.a / pba_->a_today, 2) / pba_->H0, 2);

  if (layout.idx_F0_species >= 0) {
    for (int n_dr = 0; n_dr < pba_->N_decay_dr; ++n_dr) {
      const int base_idx        = layout.idx_F0_species + n_dr * (layout.l_max + 1);
      const double r_dr_species = r_prefactor * pvecback[bgm_->index_bg_rho_dr_species_ + n_dr];
      y[base_idx + 0]           = ctx.delta_dr * r_dr_species;
      if (layout.l_max >= 1)
        y[base_idx + 1] = 4. / (3. * ctx.k) * ctx.theta_ur * r_dr_species;
      if (layout.l_max >= 2)
        y[base_idx + 2] = 2. * ctx.shear_ur * r_dr_species;
      if (layout.l_max >= 3)
        y[base_idx + 3] = ctx.l3_ur * r_dr_species;
    }
  }

  if (layout.idx_F0_sum >= 0) {
    const double r_dr_sum = r_prefactor * pvecback[bgm_->index_bg_rho_dr_];
    y[layout.idx_F0_sum]  = ctx.delta_dr * r_dr_sum;
    if (layout.l_max >= 1)
      y[layout.idx_F0_sum + 1] = 4. / (3. * ctx.k) * ctx.theta_ur * r_dr_sum;
    if (layout.l_max >= 2)
      y[layout.idx_F0_sum + 2] = 2. * ctx.shear_ur * r_dr_sum;
    if (layout.l_max >= 3)
      y[layout.idx_F0_sum + 3] = ctx.l3_ur * r_dr_sum;
  }
}

void DarkRadiationSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                              const PerturbationsModule& mod,
                                              enum file_format fmt,
                                              BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_dr", mod.index_tp_delta_dr_, pba->has_dr);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_dr", mod.index_tp_theta_dr_, pba->has_dr);
  }
}

void DarkRadiationSpecies::PrintVariables(PerturbColumnWriter& w,
                                          double /*tau*/,
                                          const double* y,
                                          const PerturbationsModule& mod,
                                          const perturb_workspace* ppw) const {
  if (pba_->has_dr != _TRUE_)
    return;

  // l_max_dr is the same in both modes (ppr->l_max_dr == ppw->pv->l_max_dr at runtime)
  const int l_max_dr = mod.GetPrecision()->l_max_dr;

  // Pre-compute gauge correction scalars (zero in title mode — ppw is nullptr)
  double k = 0., a = 0., H = 0., alpha_corr = 0., decay_corr = 0.;
  const perturb_vector* pv = nullptr;
  if (!w.IsTitleMode()) {
    pv = ppw->pv;
    k  = ppw->scalar_ctx.k;
    a  = ppw->pvecback[bgm_->index_bg_a_];
    H  = ppw->pvecback[bgm_->index_bg_H_];
    if (mod.GetPerturbs()->gauge == synchronous) {
      alpha_corr              = ppw->pvecmetric[ppw->index_mt_alpha];
      const double rho_dcdm   = dcdm_ ? dcdm_->Rho(ppw->pvecback) : 0.;
      const double rho_dr_tot = ppw->pvecback[bgm_->index_bg_rho_dr_];
      decay_corr              = a * pba_->Gamma_dcdm * rho_dcdm / rho_dr_tot;
    }
  }

  // Per-decay-DR-species: delta, theta, shear
  char tmp[64];
  for (int n_dr = 0; n_dr < pba_->N_decay_dr; n_dr++) {
    double delta_dr = 0., theta_dr = 0., shear_dr = 0.;
    if (!w.IsTitleMode()) {
      const double rho_dr_n = ppw->pvecback[bgm_->index_bg_rho_dr_species_ + n_dr];
      const double r_dr     = (a * a / pba_->H0) * (a * a / pba_->H0) * rho_dr_n;
      const int base        = pv->index_pt_F0_dr_species + n_dr * (pv->l_max_dr + 1);
      delta_dr              = y[base + 0] / r_dr + (-4. * a * H + decay_corr) * alpha_corr;
      theta_dr              = y[base + 1] * 3. / 4. * k / r_dr + k * k * alpha_corr;
      shear_dr              = y[base + 2] * 0.5 / r_dr;
    }
    snprintf(tmp, sizeof(tmp), "delta_ncdm[%d]", n_dr);
    w.Add(tmp, delta_dr, true);
    snprintf(tmp, sizeof(tmp), "theta_ncdm[%d]", n_dr);
    w.Add(tmp, theta_dr, true);
    snprintf(tmp, sizeof(tmp), "shear_ncdm[%d]", n_dr);
    w.Add(tmp, shear_dr, true);
  }

  // Raw F multipoles per DR species
  for (int n_dr = 0; n_dr < pba_->N_decay_dr; n_dr++) {
    for (int l = 0; l <= l_max_dr; l++) {
      double F_val = 0.;
      if (!w.IsTitleMode()) {
        const int base = pv->index_pt_F0_dr_species + n_dr * (pv->l_max_dr + 1);
        F_val          = y[base + l];
      }
      snprintf(tmp, sizeof(tmp), "F_dr[%d][%d]", n_dr, l);
      w.Add(tmp, F_val, true);
    }
  }
}
