#include "dncdm_decay_radiation_species.h"

#include <cmath>

#include "background_module.h"
#include "perturbations_module.h"

// ── Background ─────────────────────────────────────────────────────────────

void DNCDM_DecayRadiationSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_rho_ = index_bg++;
}

void DNCDM_DecayRadiationSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_rho_ = index_bi++;
}

void DNCDM_DecayRadiationSpecies::ComputeBackground(double /*a_rel*/,
                                                    const double* pvecback_B,
                                                    double* pvecback) {
  pvecback[index_bg_rho_] = pvecback_B[index_bi_rho_];
}

void DNCDM_DecayRadiationSpecies::BackgroundDerivs(double /*tau*/,
                                                   const double* y,
                                                   double* dy,
                                                   const double* pvecback) {
  const double a = pvecback[bgm_->index_bg_a_];
  const double H = pvecback[bgm_->index_bg_H_];

  // Dilution only; DNCDM decay source added by DNCDM_DR_Species::BackgroundDerivs
  dy[index_bi_rho_] = -4. * a * H * y[index_bi_rho_];
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DNCDM_DecayRadiationSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                              perturb_vector* pv,
                                                              const precision* /*ppr*/,
                                                              int& index_pt,
                                                              const perturb_workspace* /*ppw*/,
                                                              int /*gauge*/) {
  auto& layout   = static_cast<PerturbLayout&>(base);
  layout.idx_F0  = index_pt;
  layout.l_max   = pv->l_max_dr;
  index_pt      += pv->l_max_dr + 1;
}

void DNCDM_DecayRadiationSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
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
  const double* pvecback          = ppw->pvecback;
  const int base_idx              = layout.idx_F0;
  const int lmax                  = layout.l_max;
  double r_dr = pvecback[index_bg_rho_] * (a * a) * (a * a) / (pba_->H0 * pba_->H0);

  // l=0: free-streaming only (coupling source added by DNCDM_DR_Species::AddCouplingDerivs)
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

double DNCDM_DecayRadiationSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                                          const perturb_vector* /*pv*/,
                                          const double* y,
                                          const double* pvecback,
                                          const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[layout.idx_F0] / pvecback[index_bg_rho_];
}

double DNCDM_DecayRadiationSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                                          const perturb_vector* /*pv*/,
                                          const double* y,
                                          const double* pvecback,
                                          const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  double k             = ppw->scalar_ctx.k;
  return 0.75 * k * rho_dr_over_f * y[layout.idx_F0 + 1] / pvecback[index_bg_rho_];
}

double DNCDM_DecayRadiationSpecies::DeltaP(const BaseSpecies::PerturbLayout& base,
                                           const perturb_vector* /*pv*/,
                                           const double* y,
                                           const double* pvecback,
                                           const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[layout.idx_F0] / 3.;
}

double DNCDM_DecayRadiationSpecies::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                                  const perturb_vector* /*pv*/,
                                                  const double* y,
                                                  const double* pvecback,
                                                  const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return 2. / 3. * rho_dr_over_f * y[layout.idx_F0 + 2];
}