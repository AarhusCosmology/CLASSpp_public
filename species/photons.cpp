#include "photons.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

// ── Perturbations ─────────────────────────────────────────────────────────────

void PhotonsSpecies::RegisterPerturbationIndices(
    perturb_vector* pv, int& index_pt,
    const perturb_workspace* ppw, int /*gauge*/) {

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
  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) return;

  class_define_index(pv->index_pt_delta_g, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_g, _TRUE_, index_pt, 1);

  /* TCA active: only delta and theta; shear/polarization are analytic */
  if (ppw->approx[ppw->index_ap_tca] == (int)tca_on) return;

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

void PhotonsSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                    const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector*    pv     = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;
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
  if (ppw->approx[ppw->index_ap_rsa] == (int)rsa_on) return;

  const double dkappa = ppw->pvecthermo[
      ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[pv->index_pt_delta_g];
  const double theta_g = y[pv->index_pt_theta_g];
  const double theta_b = ctx.theta_b;

  /* delta equation (same whether TCA on or off) */
  dy[pv->index_pt_delta_g] = -4./3. * (theta_g + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == (int)tca_on) {
    /* TCA: photon theta sourced by baryon theta derivative (already set by BaryonsSpecies) */
    const double delta_p_b_over_rho_b = ctx.delta_p_b_over_rho_b;
    dy[pv->index_pt_theta_g] =
        -(dy[pv->index_pt_theta_b] + a_prime_over_a*theta_b - k2*delta_p_b_over_rho_b) / R
        + k2*(0.25*delta_g - s2_squared*ppw->tca_shear_g)
        + (1. + R)/R * metric_euler;
    return;
  }

  /* Full hierarchy (TCA off) */
  const double P0 = (y[pv->index_pt_pol0_g] + y[pv->index_pt_pol2_g]
                     + 2.*s_l[2]*y[pv->index_pt_shear_g]) / 8.;

  dy[pv->index_pt_theta_g] =
      k2*(delta_g/4. - s2_squared*y[pv->index_pt_shear_g])
      + metric_euler
      + dkappa*(theta_b - theta_g);

  dy[pv->index_pt_shear_g] =
      0.5*(8./15.*(theta_g + metric_shear)
           - 3./5.*k*s_l[3]/s_l[2]*y[pv->index_pt_l3_g]
           - dkappa*(2.*y[pv->index_pt_shear_g] - 4./5./s_l[2]*P0));

  /* l = 3 */
  int l = 3;
  dy[pv->index_pt_l3_g] =
      k/(2.*l + 1.)*(l*s_l[l]*2.*s_l[2]*y[pv->index_pt_shear_g]
                     - (l + 1.)*s_l[l + 1]*y[pv->index_pt_l3_g + 1])
      - dkappa*y[pv->index_pt_l3_g];

  /* l = 4 .. l_max_g - 1 */
  for (l = 4; l < pv->l_max_g; l++) {
    dy[pv->index_pt_delta_g + l] =
        k/(2.*l + 1.)*(l*s_l[l]*y[pv->index_pt_delta_g + l - 1]
                       - (l + 1.)*s_l[l + 1]*y[pv->index_pt_delta_g + l + 1])
        - dkappa*y[pv->index_pt_delta_g + l];
  }

  /* l = l_max_g (free-streaming truncation) */
  l = pv->l_max_g;
  dy[pv->index_pt_delta_g + l] =
      k*(s_l[l]*y[pv->index_pt_delta_g + l - 1]
         - (1. + l)*cotKgen*y[pv->index_pt_delta_g + l])
      - dkappa*y[pv->index_pt_delta_g + l];

  /* Polarization l = 0 */
  dy[pv->index_pt_pol0_g] =
      -k*y[pv->index_pt_pol0_g + 1]
      - dkappa*(y[pv->index_pt_pol0_g] - 4.*P0);

  /* Polarization l = 1 */
  dy[pv->index_pt_pol1_g] =
      k/3.*(y[pv->index_pt_pol1_g - 1] - 2.*s_l[2]*y[pv->index_pt_pol1_g + 1])
      - dkappa*y[pv->index_pt_pol1_g];

  /* Polarization l = 2 */
  dy[pv->index_pt_pol2_g] =
      k/5.*(2.*s_l[2]*y[pv->index_pt_pol2_g - 1] - 3.*s_l[3]*y[pv->index_pt_pol2_g + 1])
      - dkappa*(y[pv->index_pt_pol2_g] - 4./5.*P0);

  /* Polarization l = 3 .. l_max_pol_g - 1 */
  for (l = 3; l < pv->l_max_pol_g; l++) {
    dy[pv->index_pt_pol0_g + l] =
        k/(2.*l + 1.)*(l*s_l[l]*y[pv->index_pt_pol0_g + l - 1]
                       - (l + 1.)*s_l[l + 1]*y[pv->index_pt_pol0_g + l + 1])
        - dkappa*y[pv->index_pt_pol0_g + l];
  }

  /* Polarization l = l_max_pol_g (truncation) */
  l = pv->l_max_pol_g;
  dy[pv->index_pt_pol0_g + l] =
      k*(s_l[l]*y[pv->index_pt_pol0_g + l - 1]
         - (1. + l)*cotKgen*y[pv->index_pt_pol0_g + l])
      - dkappa*y[pv->index_pt_pol0_g + l];
}

double PhotonsSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y,
                                     const double* pvecback, const perturb_workspace* ppw) const {
  const double rho_g = pvecback[index_bg_rho_];
  /* Use TCA-corrected shear_g from the pre-computed context when the shear
     perturbation is not independently evolved (TCA or RSA modes). */
  if (pv->index_pt_shear_g < 0)
    return 4./3. * rho_g * ppw->scalar_ctx.shear_g;
  return 4./3. * rho_g * y[pv->index_pt_shear_g];
}
