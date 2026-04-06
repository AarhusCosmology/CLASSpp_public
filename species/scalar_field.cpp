#include "scalar_field.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "background_module.h"

double ScalarFieldSpecies::V_scf(double phi) const {
  double lambda = pba_.scf_parameters[0];
  double alpha  = pba_.scf_parameters[1];
  double A      = pba_.scf_parameters[2];
  double B      = pba_.scf_parameters[3];
  return exp(-lambda * phi) * (pow(phi - B, alpha) + A);
}

double ScalarFieldSpecies::dV_scf(double phi) const {
  double lambda = pba_.scf_parameters[0];
  double alpha  = pba_.scf_parameters[1];
  double A      = pba_.scf_parameters[2];
  double B      = pba_.scf_parameters[3];
  double Ve = exp(-lambda * phi);
  double Vp = pow(phi - B, alpha) + A;
  double dVe = -lambda * Ve;
  double dVp = alpha * pow(phi - B, alpha - 1);
  return dVe * Vp + Ve * dVp;
}

double ScalarFieldSpecies::ddV_scf(double phi) const {
  double lambda = pba_.scf_parameters[0];
  double alpha  = pba_.scf_parameters[1];
  double A      = pba_.scf_parameters[2];
  double B      = pba_.scf_parameters[3];
  double Ve   = exp(-lambda * phi);
  double Vp   = pow(phi - B, alpha) + A;
  double dVe  = -lambda * Ve;
  double dVp  = alpha * pow(phi - B, alpha - 1);
  double ddVe = lambda * lambda * Ve;
  double ddVp = alpha * (alpha - 1.) * pow(phi - B, alpha - 2);
  return ddVe * Vp + 2 * dVe * dVp + Ve * ddVp;
}

void ScalarFieldSpecies::ComputeBackground(double a_rel, const double* pvecback_B,
                                            double* pvecback) {
  const double a         = a_rel * pba_.a_today;
  const double phi       = pvecback_B[index_bi_phi_scf_];
  const double phi_prime = pvecback_B[index_bi_phi_prime_scf_];
  pvecback[index_bg_phi_scf_]       = phi;
  pvecback[index_bg_phi_prime_scf_] = phi_prime;
  pvecback[index_bg_V_scf_]         = V_scf(phi);
  pvecback[index_bg_dV_scf_]        = dV_scf(phi);
  pvecback[index_bg_ddV_scf_]       = ddV_scf(phi);
  pvecback[index_bg_rho_] = (phi_prime * phi_prime / (2.*a*a) + V_scf(phi)) / 3.;
  pvecback[index_bg_p_]   = (phi_prime * phi_prime / (2.*a*a) - V_scf(phi)) / 3.;
  pvecback[index_bg_p_prime_scf_] = 0.;
}

void ScalarFieldSpecies::BackgroundDerivs(double /*tau*/, const double* y, double* dy,
                                           const double* pvecback) {
  const double a         = pvecback[bgm_->index_bg_a_];
  const double H         = pvecback[bgm_->index_bg_H_];
  const double phi       = y[index_bi_phi_scf_];
  const double phi_prime = y[index_bi_phi_prime_scf_];
  /** phi'' + 2*a*H*phi' + a^2*dV = 0 */
  dy[index_bi_phi_scf_]       = phi_prime;
  dy[index_bi_phi_prime_scf_] = -a*(2.*H*phi_prime + a*dV_scf(phi));
}

void ScalarFieldSpecies::RegisterPerturbationIndices(perturb_vector* pv, const precision* /*ppr*/, int& index_pt,
                                                     const perturb_workspace* /*ppw*/,
                                                     int /*gauge*/) {  class_define_index(pv->index_pt_phi_scf,       _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_phi_prime_scf, _TRUE_, index_pt, 1);
}

void ScalarFieldSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                        const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector*    pv     = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  /* Use BackgroundModule's canonical pvecback slots for background quantities */
  auto bgm = ppaw.perturbations_module->GetBackgroundModule();
  const double phi_prime_bg = ppw->pvecback[bgm->index_bg_phi_prime_scf_];
  const double ddV_bg       = ppw->pvecback[bgm->index_bg_ddV_scf_];

  const double k2             = ctx.k2;
  const double a2             = ctx.a2;
  const double a_prime_over_a = ctx.a_prime_over_a;
  const double metric_continuity = ctx.metric_continuity;

  dy[pv->index_pt_phi_scf] = y[pv->index_pt_phi_prime_scf];
  dy[pv->index_pt_phi_prime_scf] =
      -2.*a_prime_over_a * y[pv->index_pt_phi_prime_scf]
      - metric_continuity * phi_prime_bg
      - (k2 + a2*ddV_bg) * y[pv->index_pt_phi_scf];
}

double ScalarFieldSpecies::Delta(const perturb_vector* pv, const double* y,
                                  const double* pvecback, const perturb_workspace* ppw) const {
  const double rho = pvecback[index_bg_rho_];
  if (rho == 0.) return 0.;
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double dV        = pvecback[index_bg_dV_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;

  double delta_rho = (1./3.) *
    (1./a2 * phi_prime * y[pv->index_pt_phi_prime_scf]
     + dV * y[pv->index_pt_phi_scf]);

  if (ppw->scalar_ctx.gauge == newtonian) {
    double psi = y[pv->index_pt_phi] - 4.5 * (a2/k2) * ppw->rho_plus_p_shear;
    delta_rho -= (1./3.) * (1./a2) * phi_prime * phi_prime * psi;
  }

  return delta_rho / rho;
}

double ScalarFieldSpecies::Theta(const perturb_vector* pv, const double* y,
                                  const double* pvecback, const perturb_workspace* ppw) const {
  const double rho = pvecback[index_bg_rho_];
  const double p   = pvecback[index_bg_p_];
  const double rho_plus_p = rho + p;
  if (rho_plus_p == 0.) return 0.;
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;
  return (1./3.) * k2 / a2 * phi_prime * y[pv->index_pt_phi_scf] / rho_plus_p;
}

double ScalarFieldSpecies::DeltaP(const perturb_vector* pv, const double* y,
                                   const double* pvecback, const perturb_workspace* ppw) const {
  const double phi_prime = pvecback[index_bg_phi_prime_scf_];
  const double dV        = pvecback[index_bg_dV_scf_];
  const double a2        = ppw->scalar_ctx.a2;
  const double k2        = ppw->scalar_ctx.k2;

  double delta_p = (1./3.) *
    (1./a2 * phi_prime * y[pv->index_pt_phi_prime_scf]
     - dV * y[pv->index_pt_phi_scf]);

  if (ppw->scalar_ctx.gauge == newtonian) {
    double psi = y[pv->index_pt_phi] - 4.5 * (a2/k2) * ppw->rho_plus_p_shear;
    delta_p -= (1./3.) * (1./a2) * phi_prime * phi_prime * psi;
  }

  return delta_p;
}
