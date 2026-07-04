#include "photons.h"

#include "background_module.h"
#include "baryons.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

// ── Perturbation index registration ──────────────────────────────────────────

void PhotonsSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                 perturb_vector* /*pv*/,
                                                 const precision* ppr,
                                                 int& index_pt,
                                                 const perturb_workspace* ppw,
                                                 int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  /* Initialize all photon indices to sentinel -1 so that any early return
     (RSA, TCA) leaves them in a well-defined state that guards can check. */
  layout.idx_delta = layout.idx_theta = layout.idx_shear = layout.idx_l3 = -1;
  layout.idx_pol0 = layout.idx_pol1 = layout.idx_pol2 = layout.idx_pol3 = -1;

  /* Copy l_max values from precision into layout for later use in derivs/switch. */
  layout.l_max     = ppr->l_max_g;
  layout.l_max_pol = ppr->l_max_pol_g;

  /* RSA active: photons are handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;

  layout.idx_delta = index_pt++;
  layout.idx_theta = index_pt++;

  /* TCA active: only delta and theta; shear/polarization are analytic */
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  /* Full hierarchy: shear (l=2), l3..l_max_g, polarization */
  layout.idx_shear = index_pt++;

  /* l=3 .. l_max_g → (l_max_g - 2) slots */
  if (layout.l_max >= 3) {
    layout.idx_l3  = index_pt;
    index_pt      += layout.l_max - 2;
  }

  layout.idx_pol0 = index_pt++;
  layout.idx_pol1 = index_pt++;
  layout.idx_pol2 = index_pt++;

  /* pol l=3 .. l_max_pol_g → (l_max_pol_g - 2) slots */
  if (layout.l_max_pol >= 3) {
    layout.idx_pol3  = index_pt;
    index_pt        += layout.l_max_pol - 2;
  }
}

void PhotonsSpecies::RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                       perturb_vector* /*pv*/,
                                                       const precision* ppr,
                                                       int& index_pt,
                                                       const perturb_workspace* ppw,
                                                       int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = layout.idx_theta = layout.idx_shear = layout.idx_l3 = -1;
  layout.idx_pol0 = layout.idx_pol1 = layout.idx_pol2 = layout.idx_pol3 = -1;

  /* Photons read their own l_max values from ppr — no caller pre-setup. */
  layout.l_max     = ppr->l_max_g_ten;
  layout.l_max_pol = ppr->l_max_pol_g_ten;

  /* In vector mode there is no TCA shortcut: either full hierarchy or nothing. */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  layout.idx_delta  = index_pt++;
  layout.idx_theta  = index_pt++;
  layout.idx_shear  = index_pt++;
  layout.idx_l3     = index_pt;
  index_pt         += layout.l_max - 2;

  layout.idx_pol0  = index_pt++;
  layout.idx_pol1  = index_pt++;
  layout.idx_pol2  = index_pt++;
  layout.idx_pol3  = index_pt;
  index_pt        += layout.l_max_pol - 2;
}

void PhotonsSpecies::RegisterTensorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                       perturb_vector* /*pv*/,
                                                       const precision* ppr,
                                                       int& index_pt,
                                                       const perturb_workspace* ppw,
                                                       int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = layout.idx_theta = layout.idx_shear = layout.idx_l3 = -1;
  layout.idx_pol0 = layout.idx_pol1 = layout.idx_pol2 = layout.idx_pol3 = -1;

  /* Photons read their own l_max values from ppr — no caller pre-setup. */
  layout.l_max     = ppr->l_max_g_ten;
  layout.l_max_pol = ppr->l_max_pol_g_ten;

  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on)
    return;
  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on)
    return;

  layout.idx_delta  = index_pt++;
  layout.idx_theta  = index_pt++;
  layout.idx_shear  = index_pt++;
  layout.idx_l3     = index_pt;
  index_pt         += layout.l_max - 2;

  layout.idx_pol0  = index_pt++;
  layout.idx_pol1  = index_pt++;
  layout.idx_pol2  = index_pt++;
  layout.idx_pol3  = index_pt;
  index_pt        += layout.l_max_pol - 2;
}

// ── PerturbDerivs ─────────────────────────────────────────────────────────────

void PhotonsSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                   double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();

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
  if (layout.idx_delta < 0)
    return;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[layout.idx_delta];
  const double theta_g = y[layout.idx_theta];
  const double theta_b = ctx.theta_b;

  /* delta equation (same whether TCA on or off) */
  dy[layout.idx_delta] = -4. / 3. * (theta_g + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == (int) tca_on) {
    /* TCA: photon theta sourced by baryon theta derivative (already set by BaryonsSpecies).
       In TCA, perturbed recombination is off, so cb2 * delta_b is used. */
    const double delta_p_b_over_rho_b = ppw->pvecthermo[thm_->index_th_cb2_] * ctx.delta_b;
    dy[layout.idx_theta] =
        -(dy[ctx.b_idx_theta] + a_prime_over_a * theta_b - k2 * delta_p_b_over_rho_b) / R +
        k2 * (0.25 * delta_g - s2_squared * ppw->tca_shear_g) + (1. + R) / R * metric_euler;
    return;
  }

  /* Full hierarchy (TCA off) */
  const double P0 = (y[layout.idx_pol0] + y[layout.idx_pol2] + 2. * s_l[2] * y[layout.idx_shear]) /
                    8.;

  dy[layout.idx_theta] = k2 * (delta_g / 4. - s2_squared * y[layout.idx_shear]) + metric_euler +
                         dkappa * (theta_b - theta_g);

  dy[layout.idx_shear] = 0.5 * (8. / 15. * (theta_g + metric_shear) -
                                3. / 5. * k * s_l[3] / s_l[2] * y[layout.idx_l3] -
                                dkappa * (2. * y[layout.idx_shear] - 4. / 5. / s_l[2] * P0));

  /* l = 3 */
  {
    const int l       = 3;
    dy[layout.idx_l3] = k / (2. * l + 1.) *
                            (l * s_l[l] * 2. * s_l[2] * y[layout.idx_shear] -
                             (l + 1.) * s_l[l + 1] * y[layout.idx_l3 + 1]) -
                        dkappa * y[layout.idx_l3];
  }

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < layout.l_max; l++) {
    dy[layout.idx_delta + l] = k / (2. * l + 1.) *
                                   (l * s_l[l] * y[layout.idx_delta + l - 1] -
                                    (l + 1.) * s_l[l + 1] * y[layout.idx_delta + l + 1]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* l = l_max_g (free-streaming truncation) */
  {
    const int l              = layout.l_max;
    dy[layout.idx_delta + l] = k * (s_l[l] * y[layout.idx_delta + l - 1] -
                                    (1. + l) * cotKgen * y[layout.idx_delta + l]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* Polarization l = 0 */
  dy[layout.idx_pol0] = -k * y[layout.idx_pol0 + 1] - dkappa * (y[layout.idx_pol0] - 4. * P0);

  /* Polarization l = 1 */
  dy[layout.idx_pol1] = k / 3. * (y[layout.idx_pol1 - 1] - 2. * s_l[2] * y[layout.idx_pol1 + 1]) -
                        dkappa * y[layout.idx_pol1];

  /* Polarization l = 2 */
  dy[layout.idx_pol2] =
      k / 5. * (2. * s_l[2] * y[layout.idx_pol2 - 1] - 3. * s_l[3] * y[layout.idx_pol2 + 1]) -
      dkappa * (y[layout.idx_pol2] - 4. / 5. * P0);

  /* Polarization l = 3 .. l_max_pol_g - 1 */
  for (int l = 3; l < layout.l_max_pol; l++) {
    dy[layout.idx_pol0 + l] = k / (2. * l + 1.) *
                                  (l * s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (l + 1.) * s_l[l + 1] * y[layout.idx_pol0 + l + 1]) -
                              dkappa * y[layout.idx_pol0 + l];
  }

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l             = layout.l_max_pol;
    dy[layout.idx_pol0 + l] = k * (s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (1. + l) * cotKgen * y[layout.idx_pol0 + l]) -
                              dkappa * y[layout.idx_pol0 + l];
  }
}

void PhotonsSpecies::PerturbVectorDerivs(const BaseSpecies::PerturbLayout& base,
                                         double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout           = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw = ppaw.ppw;
  const double* s_l            = ppw->s_l.data();

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (layout.idx_delta < 0)
    return;

  const auto gauge     = ppaw.perturbations_module->GetPerturbs()->gauge;
  const double k       = ppaw.k;
  const double k2      = k * k;
  const double cotKgen = ppw->cotKgen;
  const double ssqrt3  = sqrt(1. - 2. * pba_.K / k2);

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[layout.idx_delta];
  const double theta_g = y[layout.idx_theta];
  const double shear_g = y[layout.idx_shear];

  /* P^{(1)}: polarization combination (Eq. B.23 in arXiv:1305.3261) */
  const double P1 = -_SQRT6_ / 40. *
                    (4. / (3. * k) * theta_g + y[layout.idx_delta + 3] + 2. * y[layout.idx_pol0] +
                     10. / 7. * y[layout.idx_pol2] - 4. / 7. * y[layout.idx_pol0 + 4]);

  /* Vector-mode baryon theta_b lives in the Baryons layout (registered in the
     module's _vectors_ branch via `b_vec_lay.idx_theta`). */
  const auto& species   = ppaw.perturbations_module->all_species_;
  const size_t b_phv_i  = species.index_of("Baryons");
  const auto& b_phv_lay = static_cast<const BaryonsSpecies::PerturbLayout&>(
      *ppw->pv->species_layouts[b_phv_i]);
  const double theta_b_vec = y[b_phv_lay.idx_theta];

  if (gauge == possible_gauges::synchronous) {
    dy[layout.idx_delta] = -4. / 3. * theta_g - dkappa * (delta_g + 2. * _SQRT2_ * theta_b_vec);

    dy[layout.idx_theta] = k2 * (delta_g / 4. - s_l[2] * shear_g) -
                           dkappa * (theta_g + 4.0 / _SQRT6_ * P1) +
                           4.0 / (3.0 * _SQRT2_) * ssqrt3 * y[ppw->pv->index_pt_hv_prime];
  }
  else { /* newtonian */
    dy[layout.idx_delta] = -4. / 3. * theta_g - dkappa * (delta_g + 2. * _SQRT2_ * theta_b_vec) -
                           2. * _SQRT2_ * ppw->pvecmetric[ppw->index_mt_V_prime];

    dy[layout.idx_theta] = k2 * (delta_g / 4. - s_l[2] * shear_g) -
                           dkappa * (theta_g + 4.0 / _SQRT6_ * P1);
  }

  /* photon shear (F_2/2) */
  dy[layout.idx_shear] = 4. / 15. * s_l[2] * theta_g -
                         3. / 10. * k * s_l[3] * y[layout.idx_shear + 1] - dkappa * shear_g;

  /* l = 3 */
  dy[layout.idx_l3] = k / 7. * (6. * s_l[3] * shear_g - 4. * s_l[4] * y[layout.idx_l3 + 1]) -
                      dkappa * y[layout.idx_l3];

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < layout.l_max; l++)
    dy[layout.idx_delta + l] = k / (2. * l + 1.) *
                                   (l * s_l[l] * y[layout.idx_delta + l - 1] -
                                    (l + 1.) * s_l[l + 1] * y[layout.idx_delta + l + 1]) -
                               dkappa * y[layout.idx_delta + l];

  /* l = l_max_g (truncation) */
  {
    const int l              = layout.l_max;
    dy[layout.idx_delta + l] = k * (s_l[l] * y[layout.idx_delta + l - 1] -
                                    (1. + l) * cotKgen * y[layout.idx_delta + l]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* Polarization l = 0 (G_0) */
  dy[layout.idx_pol0] = -k * y[layout.idx_pol0 + 1] - dkappa * (y[layout.idx_pol0] - _SQRT6_ * P1);

  /* Polarization l = 1 .. l_max_pol_g - 1 */
  for (int l = 1; l < layout.l_max_pol; l++)
    dy[layout.idx_pol0 + l] = k / (2. * l + 1.) *
                                  (l * s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (l + 1.) * s_l[l + 1] * y[layout.idx_pol0 + l + 1]) -
                              dkappa * y[layout.idx_pol0 + l];

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l             = layout.l_max_pol;
    dy[layout.idx_pol0 + l] = k * (s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (1. + l) * cotKgen * y[layout.idx_pol0 + l]) -
                              dkappa * y[layout.idx_pol0 + l];
  }
}

void PhotonsSpecies::PerturbTensorDerivs(const BaseSpecies::PerturbLayout& base,
                                         double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout           = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw = ppaw.ppw;
  const double* s_l            = ppw->s_l.data();

  /* Nothing to do when RSA or TCA is active (indices not registered). */
  if (layout.idx_delta < 0)
    return;

  const double k       = ppaw.k;
  const double k2      = k * k;
  const double cotKgen = ppw->cotKgen;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[layout.idx_delta];
  const double theta_g = y[layout.idx_theta];
  const double shear_g = y[layout.idx_shear];

  /* P^{(2)}: polarization combination */
  const double P2 = -1.0 / _SQRT6_ *
                    (1. / 10. * delta_g + 2. / 7. * shear_g + 3. / 70. * y[layout.idx_delta + 4] -
                     3. / 5. * y[layout.idx_pol0] + 6. / 7. * y[layout.idx_pol2] -
                     3. / 70. * y[layout.idx_pol0 + 4]);

  dy[layout.idx_delta] = -4. / 3. * theta_g - dkappa * (delta_g + _SQRT6_ * P2) +
                         _SQRT6_ * y[ppw->pv->index_pt_gwdot]; /* gravitational wave source */

  dy[layout.idx_theta] = k2 * (delta_g / 4. - s_l[2] * shear_g) - dkappa * theta_g;

  dy[layout.idx_shear] = 4. / 15. * s_l[2] * theta_g -
                         3. / 10. * k * s_l[3] * y[layout.idx_shear + 1] - dkappa * shear_g;

  /* l = 3 */
  dy[layout.idx_l3] = k / 7. * (6. * s_l[3] * shear_g - 4. * s_l[4] * y[layout.idx_l3 + 1]) -
                      dkappa * y[layout.idx_l3];

  /* l = 4 .. l_max_g - 1 */
  for (int l = 4; l < layout.l_max; l++)
    dy[layout.idx_delta + l] = k / (2. * l + 1.) *
                                   (l * s_l[l] * y[layout.idx_delta + l - 1] -
                                    (l + 1.) * s_l[l + 1] * y[layout.idx_delta + l + 1]) -
                               dkappa * y[layout.idx_delta + l];

  /* l = l_max_g (truncation) */
  {
    const int l              = layout.l_max;
    dy[layout.idx_delta + l] = k * (s_l[l] * y[layout.idx_delta + l - 1] -
                                    (1. + l) * cotKgen * y[layout.idx_delta + l]) -
                               dkappa * y[layout.idx_delta + l];
  }

  /* Polarization l = 0 (G_0) */
  dy[layout.idx_pol0] = -k * y[layout.idx_pol0 + 1] - dkappa * (y[layout.idx_pol0] - _SQRT6_ * P2);

  /* Polarization l = 1 .. l_max_pol_g - 1 */
  for (int l = 1; l < layout.l_max_pol; l++)
    dy[layout.idx_pol0 + l] = k / (2. * l + 1.) *
                                  (l * s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (l + 1.) * s_l[l + 1] * y[layout.idx_pol0 + l + 1]) -
                              dkappa * y[layout.idx_pol0 + l];

  /* Polarization l = l_max_pol_g (truncation) */
  {
    const int l             = layout.l_max_pol;
    dy[layout.idx_pol0 + l] = k * (s_l[l] * y[layout.idx_pol0 + l - 1] -
                                   (1. + l) * cotKgen * y[layout.idx_pol0 + l]) -
                              dkappa * y[layout.idx_pol0 + l];
  }
}

// Direct fused override: one dispatch, no inner virtual calls. Each field
// reproduces its individual function bit-for-bit (incl. the shear branches).
BaseSpecies::StressEnergyContribution PhotonsSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  const double rho_g = pvecback[index_bg_rho_];
  StressEnergyContribution se;
  se.rho = rho_g;
  se.p   = rho_g / 3.;
  if (layout.idx_delta >= 0) {
    se.delta_rho = rho_g * y[layout.idx_delta];
    se.delta_p   = rho_g * y[layout.idx_delta] / 3.;
  }
  if (layout.idx_theta >= 0)
    se.rho_plus_p_theta = 4. / 3. * rho_g * y[layout.idx_theta];
  if (layout.idx_shear >= 0)
    se.rho_plus_p_shear = 4. / 3. * rho_g * y[layout.idx_shear];  // no approximation
  else if (ppw->approx[ppw->index_ap_tca] == (int) tca_off)
    se.rho_plus_p_shear = 0.;                          // RSA: shear neglected
  else if (ppt_->gauge == possible_gauges::newtonian)  // TCA, 1st order
    se.rho_plus_p_shear = 4. / 3. * rho_g *
                          (16. / 45. / ppw->pvecthermo[thm_->index_th_dkappa_] *
                           y[layout.idx_theta]);
  // else synchronous TCA: 0 here, re-set later in perturb_einstein
  return se;
}

// ── ApplyInitialConditions ────────────────────────────────────────────────────

void PhotonsSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                            double* y,
                                            const PerturbIcContext& ctx) {
  const auto& layout             = static_cast<const PerturbLayout&>(base);
  const PerturbationsModule* mod = ctx.p_mod;
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return;

  if (ctx.index_ic == mod->index_ic_ad_) {
    y[layout.idx_delta] = ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_cdi_) {
    y[layout.idx_delta] = ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_bi_) {
    y[layout.idx_delta] = ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_nid_) {
    y[layout.idx_delta] = ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
  else if (ctx.index_ic == mod->index_ic_niv_) {
    y[layout.idx_delta] = ctx.delta_g_ic;
    y[layout.idx_theta] = ctx.theta_g_ic;
  }
}

// ── FillSources ───────────────────────────────────────────────────────────────

void PhotonsSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                 const double* y,
                                 const double* /*dy*/,
                                 PerturbSourceContext& ctx) const {
  const auto& layout         = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;

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
    delta_g = y[layout.idx_delta];

  if (index_tp_delta_ >= 0) {
    set_source(index_tp_delta_, delta_g + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);
  }

  if (index_tp_theta_ >= 0) {
    double theta_g;
    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)
      theta_g = y[layout.idx_theta];
    else
      theta_g = ppw->rsa_theta_g;
    set_source(index_tp_theta_, theta_g + ctx.theta_shift);
  }
}

// ── Newtonian-gauge transform ─────────────────────────────────────────────────

void PhotonsSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                   double* y,
                                                   const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

// ── Switch-copy hook ──────────────────────────────────────────────────────────

void PhotonsSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                   const BaseSpecies::PerturbLayout& new_base,
                                                   const double* old_y,
                                                   double* new_y,
                                                   const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);

  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;

  // delta, theta, shear, l3, pol0..pol3, all multipoles up to l_max.
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];

  if (old_l.idx_shear >= 0 && new_l.idx_shear >= 0) {
    new_y[new_l.idx_shear] = old_y[old_l.idx_shear];
    if (old_l.idx_l3 >= 0 && new_l.idx_l3 >= 0) {
      new_y[new_l.idx_l3] = old_y[old_l.idx_l3];
      for (int l = 4; l <= new_l.l_max; ++l)
        new_y[new_l.idx_delta + l] = old_y[old_l.idx_delta + l];
    }
    if (new_l.idx_pol0 >= 0 && old_l.idx_pol0 >= 0) {
      new_y[new_l.idx_pol0] = old_y[old_l.idx_pol0];
      new_y[new_l.idx_pol1] = old_y[old_l.idx_pol1];
      new_y[new_l.idx_pol2] = old_y[old_l.idx_pol2];
      if (old_l.idx_pol3 >= 0 && new_l.idx_pol3 >= 0) {
        new_y[new_l.idx_pol3] = old_y[old_l.idx_pol3];
        for (int l = 4; l <= new_l.l_max_pol; ++l)
          new_y[new_l.idx_pol0 + l] = old_y[old_l.idx_pol0 + l];
      }
    }
  }
}

// ── Output ────────────────────────────────────────────────────────────────────

void PhotonsSpecies::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void PhotonsSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                        const PerturbationsModule& mod,
                                        file_format fmt,
                                        BaseSpecies::TransferColumnSection section) const {
  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
      w.Add("d_g", index_tp_delta_, index_tp_delta_ >= 0);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
      w.Add("t_g", index_tp_theta_, index_tp_theta_ >= 0);
  }
  // file_format::camb_format: photons not written separately
}

void PhotonsSpecies::PrintVariables(PerturbColumnWriter& w,
                                    double /*tau*/,
                                    const double* y,
                                    const PerturbationsModule& mod,
                                    const perturb_workspace* ppw) const {
  double delta_g = 0., theta_g = 0., shear_g = 0.;
  double pol0_g = 0., pol1_g = 0., pol2_g = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv.get();
    // Use collection_index_ (set by SpeciesCollection::freeze) to access our layout.
    const auto& g_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
    // scalar_ctx.k is always current here: the evolver sets it in perturb_derivs_member
    // before the print callback fires, so it equals pppaw->k at this call site.
    const double k           = ppw->scalar_ctx.k;
    const double* pvecthermo = ppw->pvecthermo.data();
    const double* pvecback   = ppw->pvecback.data();
    const double* pvecmetric = ppw->pvecmetric.data();

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) {
      delta_g = y[g_lay.idx_delta];
      theta_g = y[g_lay.idx_theta];
      if (ppw->approx[ppw->index_ap_tca] == (int) tca_on) {
        const double dkappa = pvecthermo[mod.GetThermodynamicsModule()->index_th_dkappa_];
        shear_g             = ppw->tca_shear_g;
        pol0_g              = 2.5 * ppw->tca_shear_g;
        pol1_g              = 7. / 12. * 6. / 7. * k / dkappa * ppw->tca_shear_g;
        pol2_g              = 0.5 * ppw->tca_shear_g;
      }
      else {
        shear_g = y[g_lay.idx_shear];
        pol0_g  = y[g_lay.idx_pol0];
        pol1_g  = y[g_lay.idx_pol1];
        pol2_g  = y[g_lay.idx_pol2];
      }
    }
    else {
      // RSA on: use approximated values; shear/pol remain 0
      delta_g = ppw->rsa_delta_g;
      theta_g = ppw->rsa_theta_g;
    }

    // Synchronous → Newtonian gauge conversion
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == possible_gauges::synchronous) {
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

// ── MarkUsedInSources ─────────────────────────────────────────────────────────

void PhotonsSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                       const perturb_workspace* ppw,
                                       int* used_in_sources) const {
  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  /* Photon temperature l>=3 and pol l=1,3+ are not needed in source evaluation
     when both rsa and tca are off. */
  if (ppw->approx[ppw->index_ap_rsa] != (int) rsa_off)
    return;
  if (ppw->approx[ppw->index_ap_tca] != (int) tca_off)
    return;

  for (int idx = g_lay.idx_l3; idx <= g_lay.idx_delta + g_lay.l_max; ++idx)
    used_in_sources[idx] = false;

  used_in_sources[g_lay.idx_pol1] = false;
  for (int idx = g_lay.idx_pol3; idx <= g_lay.idx_pol0 + g_lay.l_max_pol; ++idx)
    used_in_sources[idx] = false;
}

// ── MarkTensorUsedInSources ────────────────────────────────────────────────────

void PhotonsSpecies::MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& base,
                                             const perturb_workspace* ppw,
                                             int* used_in_sources) const {
  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  /* In tensor mode, we only need temperature l=0,2,4 and pol l=0,2,4
     when both rsa and tca are off. */
  if (ppw->approx[ppw->index_ap_rsa] != (int) rsa_off)
    return;
  if (ppw->approx[ppw->index_ap_tca] != (int) tca_off)
    return;

  used_in_sources[g_lay.idx_theta] = false;
  used_in_sources[g_lay.idx_l3]    = false;
  for (int idx = g_lay.idx_delta + 5; idx <= g_lay.idx_delta + g_lay.l_max; ++idx)
    used_in_sources[idx] = false;

  used_in_sources[g_lay.idx_pol1] = false;
  used_in_sources[g_lay.idx_pol3] = false;
  for (int idx = g_lay.idx_pol0 + 5; idx <= g_lay.idx_pol0 + g_lay.l_max_pol; ++idx)
    used_in_sources[idx] = false;
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> PhotonsSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  result.push_back({"Photons", std::make_unique<PhotonsSpecies>(*ctx.pba)});
  return result;
}
