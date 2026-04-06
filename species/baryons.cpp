#include "baryons.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics_module.h"

void BaryonsSpecies::RegisterPerturbationIndices(perturb_vector* pv, const precision* /*ppr*/,
                                                  int& index_pt,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {
  class_define_index(pv->index_pt_delta_b, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_b, _TRUE_, index_pt, 1);
}

void BaryonsSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                    const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector*    pv     = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  const double k2                = ppaw.k * ppaw.k;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;
  const double R                 = ctx.R;
  const double s2_squared        = ctx.s2_squared;
  const double delta_p_b_over_rho_b = ctx.delta_p_b_over_rho_b;

  const double theta_b = y[pv->index_pt_theta_b];
  /* theta_g from context: RSA-corrected photon velocity */
  const double theta_g = ctx.theta_g;
  const double delta_g = ctx.delta_g;

  const double dkappa = ppw->pvecthermo[
      ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  /* density equation */
  dy[pv->index_pt_delta_b] = -(theta_b + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == (int)tca_off) {
    /* Full equation */
    dy[pv->index_pt_theta_b] =
        -a_prime_over_a*theta_b
        + metric_euler
        + k2*delta_p_b_over_rho_b
        + R*dkappa*(theta_g - theta_b);
  } else {
    /* TCA on: tight-coupling approximation for theta_b */
    dy[pv->index_pt_theta_b] =
        (-a_prime_over_a*theta_b
         + k2*(delta_p_b_over_rho_b + R*(delta_g/4. - s2_squared*ppw->tca_shear_g))
         + R*ppw->tca_slip)
        / (1. + R)
        + metric_euler;
  }
}

double BaryonsSpecies::DeltaP(const perturb_vector* /*pv*/, const double* /*y*/,
                               const double* pvecback, const perturb_workspace* ppw) const {
  return pvecback[index_bg_rho_] * ppw->scalar_ctx.delta_p_b_over_rho_b;
}