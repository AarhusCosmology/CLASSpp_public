#include "dark_radiation_species.h"
#include "background_module.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"
#include <cmath>

// ── Background ─────────────────────────────────────────────────────────────

void DarkRadiationSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_rho_dr_species_ = index_bg;
  index_bg += pba_->N_decay_dr;
  index_bg_rho_ = index_bg++;  // total DR density (from base class protected member)
}

void DarkRadiationSpecies::RegisterIntegrationIndices(int& index_bi) {
  index_bi_rho_dr_species_ = index_bi;
  index_bi += pba_->N_decay_dr;
}

void DarkRadiationSpecies::ComputeBackground(double /*a_rel*/, const double* pvecback_B,
                                              double* pvecback) {
  double rho_dr_total = 0.;
  for (int n = 0; n < pba_->N_decay_dr; ++n) {
    if (dr_) dr_->rho_species_[n] = pvecback_B[index_bi_rho_dr_species_ + n];
    pvecback[index_bg_rho_dr_species_ + n] = pvecback_B[index_bi_rho_dr_species_ + n];
    rho_dr_total += pvecback_B[index_bi_rho_dr_species_ + n];
  }
  pvecback[index_bg_rho_] = rho_dr_total;
}

void DarkRadiationSpecies::BackgroundDerivs(double /*tau*/, const double* y, double* dy,
                                             const double* pvecback) {
  const double a = pvecback[bgm_->index_bg_a_];
  const double H = pvecback[bgm_->index_bg_H_];

  for (int n = 0; n < pba_->N_decay_dr; ++n) {
    // Dilution only; DCDM decay source added by DCDM_DR_Species::BackgroundDerivs
    dy[index_bi_rho_dr_species_ + n] =
        -4. * a * H * y[index_bi_rho_dr_species_ + n];
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DarkRadiationSpecies::RegisterPerturbationIndices(perturb_vector* pv, const precision* /*ppr*/,
                                                        int& index_pt,
                                                        const perturb_workspace* /*ppw*/,
                                                        int /*gauge*/) {
  // DR sum multipoles: F0_dr_sum[l=0..l_max_dr]
  pv->index_pt_F0_dr_sum = index_pt;
  index_pt += pv->l_max_dr + 1;

  // Per-species multipoles: N_decay_dr * (l_max_dr+1) slots
  pv->index_pt_F0_dr_species = index_pt;
  index_pt_F0_dr_species_ = index_pt;
  index_pt += pba_->N_decay_dr * (pv->l_max_dr + 1);
}

void DarkRadiationSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                          const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector*    pv  = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l  = ppw->s_l;
  const double k     = ctx.k;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;
  const double metric_shear      = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;
  const double s2_squared         = ctx.s2_squared;
  const double a = ctx.a;

  const double* pvecback = ppw->pvecback;

  // Zero out the running sum
  for (int l = 0; l <= pv->l_max_dr; ++l)
    dy[pv->index_pt_F0_dr_sum + l] = 0.;

  if (pba_->N_decay_dr > 0) {
    const int base = pv->index_pt_F0_dr_species;
    double r_dr = pvecback[index_bg_rho_dr_species_]
                  * (a * a) * (a * a) / (pba_->H0 * pba_->H0);

    // l=0: free-streaming only (coupling source from DCDM added by DCDM_DR_Species::AddCouplingDerivs)
    dy[base + 0] = -k * y[base + 1] - 4./3. * metric_continuity * r_dr;
    // l=1: free-streaming only
    dy[base + 1] = k/3. * y[base + 0]
                   - 2./3. * k * y[base + 2] * s2_squared
                   + 4. * metric_euler / (3. * k) * r_dr;
    // l=2
    dy[base + 2] = 8./15. * (3./4. * k * y[base + 1] + metric_shear * r_dr)
                   - 3./5. * k * s_l[3] / s_l[2] * y[base + 3];
    // l=3
    {
      int l = 3;
      dy[base + l] = k / (2.*l + 1.) * (l * s_l[l] * s_l[2] * y[base + l - 1]
                                          - (l + 1.) * s_l[l + 1] * y[base + l + 1]);
    }
    // l=4..l_max_dr-1
    for (int l = 4; l < pv->l_max_dr; ++l)
      dy[base + l] = k / (2.*l + 1.) * (l * s_l[l] * y[base + l - 1]
                                          - (l + 1.) * s_l[l + 1] * y[base + l + 1]);
    // l=l_max_dr (truncation)
    {
      int l = pv->l_max_dr;
      dy[base + l] = k * (s_l[l] * y[base + l - 1] - (1. + l) * cotKgen * y[base + l]);
    }
    // Accumulate into sum
    for (int l = 0; l <= pv->l_max_dr; ++l)
      dy[pv->index_pt_F0_dr_sum + l] += dy[base + l];
  }
}

/**
 * DR perturbation variables F_l follow the convention of astro-ph/9907388:
 *   delta_rho_dr = rho_dr * F0 / f,  where f = rho_dr * a^4 / H0^2
 * so delta_rho_dr = (H0/a^2)^2 * F0 = rho_dr_over_f * F0.
 * The factor rho_dr_over_f = H0^2/a^4 converts F_l to physical quantities.
 */
double DarkRadiationSpecies::Delta(const perturb_vector* pv, const double* y,
                                    const double* pvecback, const perturb_workspace* /*ppw*/) const {
  if (pv->index_pt_F0_dr_sum < 0 || pvecback[index_bg_rho_] <= 0.) return 0.;
  double a = pvecback[bgm_->index_bg_a_];
  double a2 = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[pv->index_pt_F0_dr_sum] / pvecback[index_bg_rho_];
}

double DarkRadiationSpecies::Theta(const perturb_vector* pv, const double* y,
                                    const double* pvecback, const perturb_workspace* ppw) const {
  if (pv->index_pt_F0_dr_sum < 0 || pvecback[index_bg_rho_] <= 0.) return 0.;
  double a = pvecback[bgm_->index_bg_a_];
  double a2 = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  double k = ppw->scalar_ctx.k;
  // (rho+p)*theta = k * rho_dr_over_f * F1, and (rho+p) = (4/3)*rho_dr
  return 0.75 * k * rho_dr_over_f * y[pv->index_pt_F0_dr_sum + 1] / pvecback[index_bg_rho_];
}

double DarkRadiationSpecies::DeltaP(const perturb_vector* pv, const double* y,
                                     const double* pvecback, const perturb_workspace* /*ppw*/) const {
  if (pv->index_pt_F0_dr_sum < 0) return 0.;
  double a = pvecback[bgm_->index_bg_a_];
  double a2 = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[pv->index_pt_F0_dr_sum] / 3.;
}

double DarkRadiationSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y,
                                            const double* pvecback, const perturb_workspace* /*ppw*/) const {
  if (pv->index_pt_F0_dr_sum < 0) return 0.;
  double a = pvecback[bgm_->index_bg_a_];
  double a2 = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return 2./3. * rho_dr_over_f * y[pv->index_pt_F0_dr_sum + 2];
}
