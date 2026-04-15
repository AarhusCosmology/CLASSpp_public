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

void DNCDM_DecayRadiationSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                                              const precision* /*ppr*/,
                                                              int& index_pt,
                                                              const perturb_workspace* /*ppw*/,
                                                              int /*gauge*/) {
  index_pt_F0_  = index_pt;
  index_pt     += pv->l_max_dr + 1;  // Reuse l_max_dr parameter for simplicity
}

void DNCDM_DecayRadiationSpecies::PerturbDerivs(double /*tau*/,
                                                const double* y,
                                                double* dy,
                                                const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
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

  const int base = index_pt_F0_;
  double r_dr    = pvecback[index_bg_rho_] * (a * a) * (a * a) / (pba_->H0 * pba_->H0);

  // l=0: free-streaming only (coupling source from DNCDM added by DNCDM_DR_Species::AddCouplingDerivs)
  dy[base + 0] = -k * y[base + 1] - 4. / 3. * metric_continuity * r_dr;
  // l=1: free-streaming only
  dy[base + 1] = k / 3. * y[base + 0] - 2. / 3. * k * y[base + 2] * s2_squared +
                 4. * metric_euler / (3. * k) * r_dr;
  // l=2
  dy[base + 2] = 8. / 15. * (3. / 4. * k * y[base + 1] + metric_shear * r_dr) -
                 3. / 5. * k * s_l[3] / s_l[2] * y[base + 3];
  // l=3
  {
    int l        = 3;
    dy[base + l] = k / (2. * l + 1.) *
                   (l * s_l[l] * s_l[2] * y[base + l - 1] -
                    (l + 1.) * s_l[l + 1] * y[base + l + 1]);
  }
  // l=4..l_max_dr-1
  for (int l = 4; l < pv->l_max_dr; ++l)
    dy[base + l] = k / (2. * l + 1.) *
                   (l * s_l[l] * y[base + l - 1] - (l + 1.) * s_l[l + 1] * y[base + l + 1]);
  // l=l_max_dr (truncation)
  {
    int l        = pv->l_max_dr;
    dy[base + l] = k * (s_l[l] * y[base + l - 1] - (1. + l) * cotKgen * y[base + l]);
  }
}

double DNCDM_DecayRadiationSpecies::Delta(const perturb_vector* pv,
                                          const double* y,
                                          const double* pvecback,
                                          const perturb_workspace* /*ppw*/) const {
  if (index_pt_F0_ < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[index_pt_F0_] / pvecback[index_bg_rho_];
}

double DNCDM_DecayRadiationSpecies::Theta(const perturb_vector* pv,
                                          const double* y,
                                          const double* pvecback,
                                          const perturb_workspace* ppw) const {
  if (index_pt_F0_ < 0 || pvecback[index_bg_rho_] <= 0.)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  double k             = ppw->scalar_ctx.k;
  return 0.75 * k * rho_dr_over_f * y[index_pt_F0_ + 1] / pvecback[index_bg_rho_];
}

double DNCDM_DecayRadiationSpecies::DeltaP(const perturb_vector* pv,
                                           const double* y,
                                           const double* pvecback,
                                           const perturb_workspace* /*ppw*/) const {
  if (index_pt_F0_ < 0)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return rho_dr_over_f * y[index_pt_F0_] / 3.;
}

double DNCDM_DecayRadiationSpecies::RhoPlusPShear(const perturb_vector* pv,
                                                  const double* y,
                                                  const double* pvecback,
                                                  const perturb_workspace* /*ppw*/) const {
  if (index_pt_F0_ < 0)
    return 0.;
  double a             = pvecback[bgm_->index_bg_a_];
  double a2            = a * a;
  double rho_dr_over_f = (pba_->H0 / a2) * (pba_->H0 / a2);
  return 2. / 3. * rho_dr_over_f * y[index_pt_F0_ + 2];
}