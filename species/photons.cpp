#include "photons.h"

#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

// ── Perturbations ─────────────────────────────────────────────────────────────

void PhotonsSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                                 const precision* ppr,
                                                 int& index_pt,
                                                 const perturb_workspace* ppw,
                                                 int /*gauge*/) {
  /* Initialize all photon indices to sentinel -1 so that any early return
     (RSA, TCA) leaves them in a well-defined state that guards can check. */
  pv->index_pt_delta_g = -1;
  pv->index_pt_theta_g = -1;
  pv->index_pt_shear_g = -1;
  pv->index_pt_l3_g    = -1;
  pv->index_pt_pol0_g  = -1;
  pv->index_pt_pol1_g  = -1;
  pv->index_pt_pol2_g  = -1;
  pv->index_pt_pol3_g  = -1;

  /* RSA active: photons handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;

  class_define_index(pv->index_pt_delta_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_g, _TRUE_, index_pt, 1);

  /* TCA active: only delta and theta; shear/polarization are analytic */
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  /* Full hierarchy: shear (l=2), l3..l_max_g, polarization */
  class_define_index(pv->index_pt_shear_g, _TRUE_, index_pt, 1);
  /* l=3 .. l_max_g  → (l_max_g - 2) slots */
  class_define_index(pv->index_pt_l3_g, (pv->l_max_g >= 3), index_pt, pv->l_max_g - 2);

  class_define_index(pv->index_pt_pol0_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol1_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol2_g, _TRUE_, index_pt, 1);
  /* pol l=3 .. l_max_pol_g → (l_max_pol_g - 2) slots */
  class_define_index(pv->index_pt_pol3_g, (pv->l_max_pol_g >= 3), index_pt, pv->l_max_pol_g - 2);
}

void PhotonsSpecies::RegisterVectorPerturbationIndices(perturb_vector* pv,
                                                       int& index_pt,
                                                       const perturb_workspace* ppw,
                                                       int /*gauge*/) {
  pv->index_pt_delta_g = -1;
  pv->index_pt_theta_g = -1;
  pv->index_pt_shear_g = -1;
  pv->index_pt_l3_g    = -1;
  pv->index_pt_pol0_g  = -1;
  pv->index_pt_pol1_g  = -1;
  pv->index_pt_pol2_g  = -1;
  pv->index_pt_pol3_g  = -1;

  /* In vector mode there is no TCA shortcut: either full hierarchy or nothing. */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  class_define_index(pv->index_pt_delta_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_shear_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_l3_g, _TRUE_, index_pt, pv->l_max_g - 2);
  class_define_index(pv->index_pt_pol0_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol1_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol2_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol3_g, _TRUE_, index_pt, pv->l_max_pol_g - 2);
}

void PhotonsSpecies::RegisterTensorPerturbationIndices(perturb_vector* pv,
                                                       int& index_pt,
                                                       const perturb_workspace* ppw,
                                                       int /*gauge*/) {
  pv->index_pt_delta_g = -1;
  pv->index_pt_theta_g = -1;
  pv->index_pt_shear_g = -1;
  pv->index_pt_l3_g    = -1;
  pv->index_pt_pol0_g  = -1;
  pv->index_pt_pol1_g  = -1;
  pv->index_pt_pol2_g  = -1;
  pv->index_pt_pol3_g  = -1;

  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  class_define_index(pv->index_pt_delta_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_shear_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_l3_g, _TRUE_, index_pt, pv->l_max_g - 2);
  class_define_index(pv->index_pt_pol0_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol1_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol2_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_pol3_g, _TRUE_, index_pt, pv->l_max_pol_g - 2);
}

void PhotonsSpecies::PerturbDerivs(double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;

  const double k                 = ppaw.k;
  const double k2                = k * k;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;
  const double metric_shear      = ctx.metric_shear;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double cotKgen           = ctx.cotKgen;
  const double s2_squared        = ctx.s2_squared;
  const double R                 = ctx.R;

  /* RSA active: photons are handled analytically */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[pv->index_pt_delta_g];
  const double theta_g = y[pv->index_pt_theta_g];
  const double theta_b = ctx.theta_b;

  /* delta equation (same whether TCA on or off) */
  dy[pv->index_pt_delta_g] = -4. / 3. * (theta_g + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on) {
    /* TCA: photon theta sourced by baryon theta derivative (already set by BaryonsSpecies) */
    const double delta_p_b_over_rho_b = ctx.delta_p_b_over_rho_b;
    dy[pv->index_pt_theta_g] =
        -(dy[pv->index_pt_theta_b] + a_prime_over_a * theta_b - k2 * delta_p_b_over_rho_b) / R +
        k2 * (0.25 * delta_g - s2_squared * ppw->tca_shear_g) + (1. + R) / R * metric_euler;
    return;
  }

  /* Full hierarchy (TCA off) */
  const double P0 = (y[pv->index_pt_pol0_g] + y[pv->index_pt_pol2_g] +
                     2. * s_l[2] * y[pv->index_pt_shear_g]) /
                    8.;

  dy[pv->index_pt_theta_g] = k2 * (delta_g / 4. - s2_squared * y[pv->index_pt_shear_g]) +
                             metric_euler + dkappa * (theta_b - theta_g);

  dy[pv->index_pt_shear_g] = 0.5 *
                             (8. / 15. * (theta_g + metric_shear) -
                              3. / 5. * k * s_l[3] / s_l[2] * y[pv->index_pt_l3_g] -
                              dkappa * (2. * y[pv->index_pt_shear_g] - 4. / 5. / s_l[2] * P0));

  /* l = 3 */
  int l                 = 3;
  dy[pv->index_pt_l3_g] = k / (2. * l + 1.) *
                              (l * s_l[l] * 2. * s_l[2] * y[pv->index_pt_shear_g] -
                               (l + 1.) * s_l[l + 1] * y[pv->index_pt_l3_g + 1]) -
                          dkappa * y[pv->index_pt_l3_g];

  /* l = 4 .. l_max_g - 1 */
  for (l = 4; l < pv->l_max_g; l++) {
    dy[pv->index_pt_delta_g + l] = k / (2. * l + 1.) *
                                       (l * s_l[l] * y[pv->index_pt_delta_g + l - 1] -
                                        (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_g + l + 1]) -
                                   dkappa * y[pv->index_pt_delta_g + l];
  }

  /* l = l_max_g (free-streaming truncation) */
  l                            = pv->l_max_g;
  dy[pv->index_pt_delta_g + l] = k * (s_l[l] * y[pv->index_pt_delta_g + l - 1] -
                                      (1. + l) * cotKgen * y[pv->index_pt_delta_g + l]) -
                                 dkappa * y[pv->index_pt_delta_g + l];

  /* Polarization l = 0 */
  dy[pv->index_pt_pol0_g] = -k * y[pv->index_pt_pol0_g + 1] -
                            dkappa * (y[pv->index_pt_pol0_g] - 4. * P0);

  /* Polarization l = 1 */
  dy[pv->index_pt_pol1_g] =
      k / 3. * (y[pv->index_pt_pol1_g - 1] - 2. * s_l[2] * y[pv->index_pt_pol1_g + 1]) -
      dkappa * y[pv->index_pt_pol1_g];

  /* Polarization l = 2 */
  dy[pv->index_pt_pol2_g] = k / 5. *
                                (2. * s_l[2] * y[pv->index_pt_pol2_g - 1] -
                                 3. * s_l[3] * y[pv->index_pt_pol2_g + 1]) -
                            dkappa * (y[pv->index_pt_pol2_g] - 4. / 5. * P0);

  /* Polarization l = 3 .. l_max_pol_g - 1 */
  for (l = 3; l < pv->l_max_pol_g; l++) {
    dy[pv->index_pt_pol0_g + l] = k / (2. * l + 1.) *
                                      (l * s_l[l] * y[pv->index_pt_pol0_g + l - 1] -
                                       (l + 1.) * s_l[l + 1] * y[pv->index_pt_pol0_g + l + 1]) -
                                  dkappa * y[pv->index_pt_pol0_g + l];
  }

  /* Polarization l = l_max_pol_g (truncation) */
  l                           = pv->l_max_pol_g;
  dy[pv->index_pt_pol0_g + l] = k * (s_l[l] * y[pv->index_pt_pol0_g + l - 1] -
                                     (1. + l) * cotKgen * y[pv->index_pt_pol0_g + l]) -
                                dkappa * y[pv->index_pt_pol0_g + l];
}

void PhotonsSpecies::PerturbVectorDerivs(double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector* pv     = ppw->pv;
  const double* s_l            = ppw->s_l;

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  const int gauge      = ppaw.perturbations_module->GetPerturbs()->gauge;
  const double k       = ppaw.k;
  const double k2      = k * k;
  const double cotKgen = ppw->cotKgen;
  const double ssqrt3  = sqrt(1. - 2. * pba_.K / k2);

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[pv->index_pt_delta_g];
  const double theta_g = y[pv->index_pt_theta_g];
  const double shear_g = y[pv->index_pt_shear_g];

  /* P^{(1)}: polarization combination (Eq. B.23 in arXiv:1305.3261) */
  const double P1 = -_SQRT6_ / 40. *
                    (4. / (3. * k) * theta_g + y[pv->index_pt_delta_g + 3] +
                     2. * y[pv->index_pt_pol0_g] + 10. / 7. * y[pv->index_pt_pol2_g] -
                     4. / 7. * y[pv->index_pt_pol0_g + 4]);

  if (gauge == synchronous) {
    dy[pv->index_pt_delta_g] = -4. / 3. * theta_g -
                               dkappa * (delta_g + 2. * _SQRT2_ * y[pv->index_pt_theta_b]);

    dy[pv->index_pt_theta_g] = k2 * (delta_g / 4. - s_l[2] * shear_g) -
                               dkappa * (theta_g + 4.0 / _SQRT6_ * P1) +
                               4.0 / (3.0 * _SQRT2_) * ssqrt3 * y[pv->index_pt_hv_prime];
  }
  else { /* newtonian */

    dy[pv->index_pt_delta_g] = -4. / 3. * theta_g -
                               dkappa * (delta_g + 2. * _SQRT2_ * y[pv->index_pt_theta_b]) -
                               2. * _SQRT2_ * ppw->pvecmetric[ppw->index_mt_V_prime];

    dy[pv->index_pt_theta_g] = k2 * (delta_g / 4. - s_l[2] * shear_g) -
                               dkappa * (theta_g + 4.0 / _SQRT6_ * P1);
  }

  /* photon shear (F_2/2) */
  dy[pv->index_pt_shear_g] = 4. / 15. * s_l[2] * theta_g -
                             3. / 10. * k * s_l[3] * y[pv->index_pt_shear_g + 1] - dkappa * shear_g;

  /* l = 3 */
  dy[pv->index_pt_l3_g] = k / 7. *
                              (6. * s_l[3] * shear_g - 4. * s_l[4] * y[pv->index_pt_l3_g + 1]) -
                          dkappa * y[pv->index_pt_l3_g];

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < pv->l_max_g; l++)
    dy[pv->index_pt_delta_g + l] = k / (2. * l + 1.) *
                                       (l * s_l[l] * y[pv->index_pt_delta_g + l - 1] -
                                        (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_g + l + 1]) -
                                   dkappa * y[pv->index_pt_delta_g + l];

  /* l = l_max_g (truncation) */
  {
    const int l                  = pv->l_max_g;
    dy[pv->index_pt_delta_g + l] = k * (s_l[l] * y[pv->index_pt_delta_g + l - 1] -
                                        (1. + l) * cotKgen * y[pv->index_pt_delta_g + l]) -
                                   dkappa * y[pv->index_pt_delta_g + l];
  }

  /* Polarization l = 0 (G_0) */
  dy[pv->index_pt_pol0_g] = -k * y[pv->index_pt_pol0_g + 1] -
                            dkappa * (y[pv->index_pt_pol0_g] - _SQRT6_ * P1);

  /* Polarization l = 1 .. l_max_pol_g - 1 */
  for (int l = 1; l < pv->l_max_pol_g; l++)
    dy[pv->index_pt_pol0_g + l] = k / (2. * l + 1.) *
                                      (l * s_l[l] * y[pv->index_pt_pol0_g + l - 1] -
                                       (l + 1.) * s_l[l + 1] * y[pv->index_pt_pol0_g + l + 1]) -
                                  dkappa * y[pv->index_pt_pol0_g + l];

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l                 = pv->l_max_pol_g;
    dy[pv->index_pt_pol0_g + l] = k * (s_l[l] * y[pv->index_pt_pol0_g + l - 1] -
                                       (1. + l) * cotKgen * y[pv->index_pt_pol0_g + l]) -
                                  dkappa * y[pv->index_pt_pol0_g + l];
  }
}

void PhotonsSpecies::PerturbTensorDerivs(double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector* pv     = ppw->pv;
  const double* s_l            = ppw->s_l;

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  const double k       = ppaw.k;
  const double k2      = k * k;
  const double cotKgen = ppw->cotKgen;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[pv->index_pt_delta_g];
  const double theta_g = y[pv->index_pt_theta_g];
  const double shear_g = y[pv->index_pt_shear_g];

  /* P^{(2)}: polarization combination */
  const double P2 = -1.0 / _SQRT6_ *
                    (1. / 10. * delta_g + 2. / 7. * shear_g +
                     3. / 70. * y[pv->index_pt_delta_g + 4] - 3. / 5. * y[pv->index_pt_pol0_g] +
                     6. / 7. * y[pv->index_pt_pol2_g] - 3. / 70. * y[pv->index_pt_pol0_g + 4]);

  dy[pv->index_pt_delta_g] = -4. / 3. * theta_g - dkappa * (delta_g + _SQRT6_ * P2) +
                             _SQRT6_ * y[pv->index_pt_gwdot]; /* gravitational wave source */

  dy[pv->index_pt_theta_g] = k2 * (delta_g / 4. - s_l[2] * shear_g) - dkappa * theta_g;

  dy[pv->index_pt_shear_g] = 4. / 15. * s_l[2] * theta_g -
                             3. / 10. * k * s_l[3] * y[pv->index_pt_shear_g + 1] - dkappa * shear_g;

  /* l = 3 */
  dy[pv->index_pt_l3_g] = k / 7. *
                              (6. * s_l[3] * shear_g - 4. * s_l[4] * y[pv->index_pt_l3_g + 1]) -
                          dkappa * y[pv->index_pt_l3_g];

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < pv->l_max_g; l++)
    dy[pv->index_pt_delta_g + l] = k / (2. * l + 1.) *
                                       (l * s_l[l] * y[pv->index_pt_delta_g + l - 1] -
                                        (l + 1.) * s_l[l + 1] * y[pv->index_pt_delta_g + l + 1]) -
                                   dkappa * y[pv->index_pt_delta_g + l];

  /* l = l_max_g (truncation) */
  {
    const int l                  = pv->l_max_g;
    dy[pv->index_pt_delta_g + l] = k * (s_l[l] * y[pv->index_pt_delta_g + l - 1] -
                                        (1. + l) * cotKgen * y[pv->index_pt_delta_g + l]) -
                                   dkappa * y[pv->index_pt_delta_g + l];
  }

  /* Polarization l = 0 (G_0) */
  dy[pv->index_pt_pol0_g] = -k * y[pv->index_pt_pol0_g + 1] -
                            dkappa * (y[pv->index_pt_pol0_g] - _SQRT6_ * P2);

  /* Polarization l = 1 .. l_max_pol_g - 1 */
  for (int l = 1; l < pv->l_max_pol_g; l++)
    dy[pv->index_pt_pol0_g + l] = k / (2. * l + 1.) *
                                      (l * s_l[l] * y[pv->index_pt_pol0_g + l - 1] -
                                       (l + 1.) * s_l[l + 1] * y[pv->index_pt_pol0_g + l + 1]) -
                                  dkappa * y[pv->index_pt_pol0_g + l];

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l                 = pv->l_max_pol_g;
    dy[pv->index_pt_pol0_g + l] = k * (s_l[l] * y[pv->index_pt_pol0_g + l - 1] -
                                       (1. + l) * cotKgen * y[pv->index_pt_pol0_g + l]) -
                                  dkappa * y[pv->index_pt_pol0_g + l];
  }
}

void PhotonsSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv             = ctx.ppw->pv;
  const PerturbationsModule* mod = ctx.p_mod;
  if (pv->index_pt_delta_g < 0 || pv->index_pt_theta_g < 0)
    return;

  if (ctx.index_ic == mod->index_ic_ad_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_bi_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_nid_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_niv_) {
    y[pv->index_pt_delta_g] = ctx.delta_g_ic;
    y[pv->index_pt_theta_g] = ctx.theta_g_ic;
  }
}

void PhotonsSpecies::FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const perturb_vector* pv   = ppw->pv;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // Convenient lambda to write into source table
  auto set_source = [&](int index_tp, double value) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp, ctx.index_tau, ctx.index_k, value);
  };

  double delta_g;
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    delta_g = ppw->rsa_delta_g;
  else
    delta_g = y[pv->index_pt_delta_g];

  if (p_mod->has_source_delta_g_ == _TRUE_) {
    set_source(p_mod->index_tp_delta_g_, delta_g + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);
  }

  if (p_mod->has_source_theta_g_ == _TRUE_) {
    double theta_g;
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      theta_g = y[pv->index_pt_theta_g];
    else
      theta_g = ppw->rsa_theta_g;
    set_source(p_mod->index_tp_theta_g_, theta_g + ctx.theta_shift);
  }
}

void PhotonsSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                        const PerturbationsModule& mod,
                                        enum file_format fmt,
                                        BaseSpecies::TransferColumnSection section) const {
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_g", mod.index_tp_delta_g_, true);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_g", mod.index_tp_theta_g_, true);
  }
  // camb_format: photons not written separately
}

void PhotonsSpecies::PrintVariables(PerturbColumnWriter& w,
                                    double /*tau*/,
                                    const double* y,
                                    const PerturbationsModule& mod,
                                    const perturb_workspace* ppw) const {
  double delta_g = 0., theta_g = 0., shear_g = 0.;
  double pol0_g = 0., pol1_g = 0., pol2_g = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    // scalar_ctx.k is always current here: the evolver sets it in perturb_derivs_member
    // before the print callback fires, so it equals pppaw->k at this call site.
    const double k           = ppw->scalar_ctx.k;
    const double* pvecthermo = ppw->pvecthermo;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) {
      delta_g = y[pv->index_pt_delta_g];
      theta_g = y[pv->index_pt_theta_g];
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_on) {
        const double dkappa = pvecthermo[mod.GetThermodynamicsModule()->index_th_dkappa_];
        shear_g             = ppw->tca_shear_g;
        pol0_g              = 2.5 * ppw->tca_shear_g;
        pol1_g              = 7. / 12. * 6. / 7. * k / dkappa * ppw->tca_shear_g;
        pol2_g              = 0.5 * ppw->tca_shear_g;
      }
      else {
        shear_g = y[pv->index_pt_shear_g];
        pol0_g  = y[pv->index_pt_pol0_g];
        pol1_g  = y[pv->index_pt_pol1_g];
        pol2_g  = y[pv->index_pt_pol2_g];
      }
    }
    else {
      // RSA on: use approximated values; shear/pol remain 0
      delta_g = ppw->rsa_delta_g;
      theta_g = ppw->rsa_theta_g;
    }

    // Synchronous → Newtonian gauge conversion
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_g            -= 4. * H * a * alpha;
      theta_g            += k * k * alpha;
    }
  }

  w.Add("delta_g", delta_g, true);
  w.Add("theta_g", theta_g, true);
  w.Add("shear_g", shear_g, true);
  w.Add("pol0_g", pol0_g, true);
  w.Add("pol1_g", pol1_g, true);
  w.Add("pol2_g", pol2_g, true);
}

double PhotonsSpecies::RhoPlusPShear(const perturb_vector* pv,
                                     const double* y,
                                     const double* pvecback,
                                     const perturb_workspace* ppw) const {
  const double rho_g = pvecback[index_bg_rho_];
  /* Use TCA-corrected shear_g from the pre-computed context when the shear
     perturbation is not independently evolved (TCA or RSA modes). */
  if (pv->index_pt_shear_g < 0)
    return 4. / 3. * rho_g * ppw->scalar_ctx.shear_g;
  return 4. / 3. * rho_g * y[pv->index_pt_shear_g];
}
