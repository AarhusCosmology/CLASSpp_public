#include <cmath>

#include "../species/idm_dr.h"
#include "../species/idm_drmd.h"
#include "../species/idr.h"
#include "../species/idr_drmd.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "thermodynamics.h"
#include "thermodynamics_module.h"

// ── IDM_DR ─────────────────────────────────────────────────────────────────

void IDM_DRSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                                const precision* ppr,
                                                int& index_pt,
                                                const perturb_workspace* /*ppw*/,
                                                int /*gauge*/) {
  class_define_index(pv->index_pt_delta_idm_dr, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idm_dr, _TRUE_, index_pt, 1);
}

void IDM_DRSpecies::PerturbDerivs(double /*tau*/,
                                  const double* y,
                                  double* dy,
                                  const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  // idm_dr density: same as CDM continuity equation
  dy[pv->index_pt_delta_idm_dr] = -(y[pv->index_pt_theta_idm_dr] + ctx.metric_continuity);

  // idm_dr velocity: Hubble friction only (coupling to IDR added by IDM_DR_IDR_Species)
  dy[pv->index_pt_theta_idm_dr] = -ctx.a_prime_over_a * y[pv->index_pt_theta_idm_dr] +
                                  ctx.metric_euler;
}

// ── IDR ───────────────────────────────────────────────────────────────────

void IDRSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                             const precision* ppr,
                                             int& index_pt,
                                             const perturb_workspace* ppw,
                                             int /*gauge*/) {
  /* Initialize all IDR indices to sentinel -1 */
  pv->index_pt_delta_idr = -1;
  pv->index_pt_theta_idr = -1;
  pv->index_pt_shear_idr = -1;
  pv->index_pt_l3_idr    = -1;

  /* RSA active: IDR handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)
    return;

  class_define_index(pv->index_pt_delta_idr, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idr, _TRUE_, index_pt, 1);

  /* Full hierarchy: only register shear/l3 if idr_nature == idr_free_streaming
     AND (no IDM_DR OR TCA is off) */
  if (ppw->scalar_ctx.idr_nature == idr_free_streaming) {
    if (pba_.has_idm_dr == _FALSE_ ||
        ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off) {
      class_define_index(pv->index_pt_shear_idr, _TRUE_, index_pt, 1);
      pv->l_max_idr = pba_.l_max_idr;
      class_define_index(pv->index_pt_l3_idr, (pv->l_max_idr >= 3), index_pt, pv->l_max_idr - 2);
    }
  }
}

void IDRSpecies::PerturbDerivs(double /*tau*/,
                               const double* y,
                               double* dy,
                               const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)
    return;

  const double* s_l    = ppw->s_l;
  const double k       = ctx.k;
  const double cotKgen = ctx.cotKgen;

  // idr density
  dy[pv->index_pt_delta_idr] = -4. / 3. * (y[pv->index_pt_theta_idr] + ctx.metric_continuity);

  if (ctx.idr_nature == idr_free_streaming) {
    // idr velocity and hierarchy (only when shear is registered; under TCA shear is
    // not registered and AddCouplingDerivs will overwrite theta_idr)
    if (pv->index_pt_shear_idr >= 0) {
      dy[pv->index_pt_theta_idr] = ctx.k2 * (y[pv->index_pt_delta_idr] / 4. -
                                             ctx.s2_squared * y[pv->index_pt_shear_idr]) +
                                   ctx.metric_euler;

      // idr shear (l=2)
      dy[pv->index_pt_shear_idr] = 0.5 *
                                   (8. / 15. * (y[pv->index_pt_theta_idr] + ctx.metric_shear) -
                                    3. / 5. * k * s_l[3] / s_l[2] * y[pv->index_pt_shear_idr + 1]);

      // idr l=3
      int l                   = 3;
      dy[pv->index_pt_l3_idr] = k / (2. * l + 1.) *
                                (l * 2. * s_l[l] * s_l[2] * y[pv->index_pt_shear_idr] -
                                 (l + 1.) * s_l[l + 1] * y[pv->index_pt_l3_idr + 1]);

      // idr l>3
      for (l = 4; l < pv->l_max_idr; l++)
        dy[pv->index_pt_delta_idr + l] = k / (2. * l + 1) *
                                         (l * s_l[l] * y[pv->index_pt_delta_idr + l - 1] -
                                          (l + 1.) * s_l[l + 1] *
                                              y[pv->index_pt_delta_idr + l + 1]);

      // idr lmax
      l                              = pv->l_max_idr;
      dy[pv->index_pt_delta_idr + l] = k * (s_l[l] * y[pv->index_pt_delta_idr + l - 1] -
                                            (1. + l) * cotKgen * y[pv->index_pt_delta_idr + l]);
    }
    // else: TCA-on; theta_idr will be set by IDM_DR_IDR_Species::AddCouplingDerivs
  }
  else {
    // Fluid idr velocity
    dy[pv->index_pt_theta_idr] = ctx.k2 / 4. * y[pv->index_pt_delta_idr] + ctx.metric_euler;
  }
}

double IDRSpecies::RhoPlusPShear(const perturb_vector* pv,
                                 const double* y,
                                 const double* pvecback,
                                 const perturb_workspace* ppw) const {
  if (ppw->scalar_ctx.idr_nature != idr_free_streaming)
    return 0.;

  const double shear_idr = (pv->index_pt_shear_idr >= 0) ? y[pv->index_pt_shear_idr]
                                                         : ppw->tca_shear_idm_dr;

  return 4. / 3. * pvecback[index_bg_rho_] * shear_idr;
}

// ── IDM_DRMD ───────────────────────────────────────────────────────────────

void IDM_DRMDSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                                  const precision* ppr,
                                                  int& index_pt,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {
  class_define_index(pv->index_pt_delta_idm_drmd, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idm_drmd, _TRUE_, index_pt, 1);
}

void IDM_DRMDSpecies::PerturbDerivs(double /*tau*/,
                                    const double* y,
                                    double* dy,
                                    const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[pv->index_pt_delta_idm_drmd] = -(y[pv->index_pt_theta_idm_drmd] + ctx.metric_continuity);
  // Free-streaming velocity: coupling to IDR_DRMD added by IDM_DRMD_IDR_DRMD_Species
  dy[pv->index_pt_theta_idm_drmd] = -ctx.a_prime_over_a * y[pv->index_pt_theta_idm_drmd] +
                                    ctx.metric_euler;
}

// ── IDR_DRMD ───────────────────────────────────────────────────────────────

void IDR_DRMDSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                                  const precision* ppr,
                                                  int& index_pt,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {
  class_define_index(pv->index_pt_delta_idr_drmd, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idr_drmd, _TRUE_, index_pt, 1);
}

void IDR_DRMDSpecies::PerturbDerivs(double /*tau*/,
                                    const double* y,
                                    double* dy,
                                    const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[pv->index_pt_delta_idr_drmd] = -4. / 3. *
                                    (y[pv->index_pt_theta_idr_drmd] + ctx.metric_continuity);
  // Free-streaming velocity: coupling to IDM_DRMD added by IDM_DRMD_IDR_DRMD_Species
  dy[pv->index_pt_theta_idr_drmd] = 0.25 * ctx.k2 * y[pv->index_pt_delta_idr_drmd] +
                                    ctx.metric_euler;
}
