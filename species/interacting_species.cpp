#include "../species/idm_dr.h"
#include "../species/idr.h"
#include "../species/idm_drmd.h"
#include "../species/idr_drmd.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "background_module.h"
#include "thermodynamics_module.h"
#include "thermodynamics.h"
#include <cmath>

// ── IDM_DR ─────────────────────────────────────────────────────────────────

void IDM_DRSpecies::RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {
  class_define_index(pv->index_pt_delta_idm_dr, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idm_dr, _TRUE_, index_pt, 1);
}

void IDM_DRSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                   const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  auto* pth_mod = ppaw.perturbations_module->GetThermodynamicsModule().get();
  const double* pvecthermo = ppw->pvecthermo;

  // idm_dr density
  dy[pv->index_pt_delta_idm_dr] = -(y[pv->index_pt_theta_idm_dr] + ctx.metric_continuity);

  // idm_dr velocity
  if (ppw->approx[ppw->index_ap_tca_idm_dr] == (int)tca_idm_dr_off) {
    const double dmu_idm_dr = pvecthermo[pth_mod->index_th_dmu_idm_dr_];
    const double Sinv = ctx.R_idr;
    const double theta_idr = ctx.theta_idr;

    dy[pv->index_pt_theta_idm_dr] = -ctx.a_prime_over_a * y[pv->index_pt_theta_idm_dr] + ctx.metric_euler;
    dy[pv->index_pt_theta_idm_dr] -= (Sinv * dmu_idm_dr * (y[pv->index_pt_theta_idm_dr] - theta_idr)
                                      - ctx.k2 * pvecthermo[pth_mod->index_th_cidm_dr2_] * y[pv->index_pt_delta_idm_dr]);
  } else {
    const double Sinv = ctx.R_idr;
    dy[pv->index_pt_theta_idm_dr] = 1./(1. + Sinv)*(-ctx.a_prime_over_a*y[pv->index_pt_theta_idm_dr] + ctx.k2*pvecthermo[pth_mod->index_th_cidm_dr2_]*
                                                   y[pv->index_pt_delta_idm_dr] + ctx.k2*Sinv*(ctx.delta_idr/4. - ctx.tca_shear_idm_dr)) + ctx.metric_euler + Sinv/(1.+Sinv)*ctx.tca_slip_idm_dr;
  }
}

// ── IDR ───────────────────────────────────────────────────────────────────

void IDRSpecies::RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                               const perturb_workspace* ppw,
                                               int /*gauge*/) {
  /* Initialize all IDR indices to sentinel -1 */
  pv->index_pt_delta_idr = -1;
  pv->index_pt_theta_idr = -1;
  pv->index_pt_shear_idr = -1;
  pv->index_pt_l3_idr    = -1;

  /* RSA active: IDR handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_on) return;

  class_define_index(pv->index_pt_delta_idr, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idr, _TRUE_, index_pt, 1);

  /* Full hierarchy: only register shear/l3 if idr_nature == idr_free_streaming
     AND (no IDM_DR OR TCA is off) */
  if (ppw->scalar_ctx.idr_nature == idr_free_streaming) {
    if (pba_.has_idm_dr == _FALSE_ || ppw->approx[ppw->index_ap_tca_idm_dr] == (int)tca_idm_dr_off) {
      class_define_index(pv->index_pt_shear_idr, _TRUE_, index_pt, 1);
      pv->l_max_idr = pba_.l_max_idr;
      class_define_index(pv->index_pt_l3_idr, (pv->l_max_idr >= 3), index_pt, pv->l_max_idr - 2);
    }
  }
}

void IDRSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  if (ppw->approx[ppw->index_ap_rsa_idr] == (int)rsa_idr_off) {
    // idr density
    dy[pv->index_pt_delta_idr] = -4./3.*(y[pv->index_pt_theta_idr] + ctx.metric_continuity);

    auto* pth_mod = ppaw.perturbations_module->GetThermodynamicsModule().get();
    const thermo* pth = pth_mod->GetThermodynamics();
    const double* pvecthermo = ppw->pvecthermo;
    auto* ppt = ppaw.perturbations_module->GetPerturbs();

    const double dmu_idm_dr = pvecthermo[pth_mod->index_th_dmu_idm_dr_];
    const double dmu_idr = pth->b_idr/pth->a_idm_dr*pba_.Omega0_idr/pba_.Omega0_idm_dr*dmu_idm_dr;

    if ((pba_.has_idm_dr == _FALSE_) || (ppw->approx[ppw->index_ap_tca_idm_dr] == (int)tca_idm_dr_off)) {
      // idr velocity
      if (ctx.idr_nature == idr_free_streaming)
        dy[pv->index_pt_theta_idr] = ctx.k2*(y[pv->index_pt_delta_idr]/4. - ctx.s2_squared*y[pv->index_pt_shear_idr]) + ctx.metric_euler;
      else
        dy[pv->index_pt_theta_idr] = ctx.k2/4. * y[pv->index_pt_delta_idr] + ctx.metric_euler;

      if (pba_.has_idm_dr == _TRUE_)
        dy[pv->index_pt_theta_idr] += dmu_idm_dr*(ctx.theta_idm_dr - y[pv->index_pt_theta_idr]);

      if (ctx.idr_nature == idr_free_streaming) {
        const double* s_l = ppw->s_l;
        const double k = ctx.k;
        const double cotKgen = ctx.cotKgen;

        // idr shear
        int l = 2;
        dy[pv->index_pt_shear_idr] = 0.5*(8./15.*(y[pv->index_pt_theta_idr] + ctx.metric_shear) - 3./5.*k*s_l[3]/s_l[2]*y[pv->index_pt_shear_idr+1]);
        if (pba_.has_idm_dr == _TRUE_)
          dy[pv->index_pt_shear_idr] -= (ppt->alpha_idm_dr[l-2]*dmu_idm_dr + ppt->beta_idr[l-2]*dmu_idr)*y[pv->index_pt_shear_idr];

        // idr l=3
        l = 3;
        dy[pv->index_pt_l3_idr] = k/(2.*l+1.)*(l*2.*s_l[l]*s_l[2]*y[pv->index_pt_shear_idr] - (l+1.)*s_l[l+1]*y[pv->index_pt_l3_idr+1]);
        if (pba_.has_idm_dr == _TRUE_)
          dy[pv->index_pt_l3_idr] -= (ppt->alpha_idm_dr[l-2]*dmu_idm_dr + ppt->beta_idr[l-2]*dmu_idr)*y[pv->index_pt_l3_idr];

        // idr l>3
        for (l = 4; l < pv->l_max_idr; l++) {
          dy[pv->index_pt_delta_idr+l] = k/(2.*l+1)*(l*s_l[l]*y[pv->index_pt_delta_idr+l-1] - (l+1.)*s_l[l+1]*y[pv->index_pt_delta_idr+l+1]);
          if (pba_.has_idm_dr == _TRUE_)
            dy[pv->index_pt_delta_idr+l] -= (ppt->alpha_idm_dr[l-2]*dmu_idm_dr + ppt->beta_idr[l-2]*dmu_idr)*y[pv->index_pt_delta_idr+l];
        }

        // idr lmax
        l = pv->l_max_idr;
        dy[pv->index_pt_delta_idr+l] = k*(s_l[l]*y[pv->index_pt_delta_idr+l-1] - (1.+l)*cotKgen*y[pv->index_pt_delta_idr+l]);
        if (pba_.has_idm_dr == _TRUE_)
          dy[pv->index_pt_delta_idr+l] -= (ppt->alpha_idm_dr[l-2]*dmu_idm_dr + ppt->beta_idr[l-2]*dmu_idr)*y[pv->index_pt_delta_idr+l];
      }
    } else {
      const double Sinv = ctx.R_idr;
      dy[pv->index_pt_theta_idr] = 1./(1. + Sinv)*(-ctx.a_prime_over_a*ctx.theta_idm_dr
                                                   + ctx.k2*pvecthermo[pth_mod->index_th_cidm_dr2_]*ctx.delta_idm_dr
                                                   + ctx.k2*Sinv*(y[pv->index_pt_delta_idr]/4. - ctx.tca_shear_idm_dr)) + ctx.metric_euler - 1./(1.+Sinv)*ctx.tca_slip_idm_dr;
    }
  }
}

double IDRSpecies::RhoPlusPShear(const perturb_vector* pv, const double* y,
                                   const double* pvecback, const perturb_workspace* ppw) const {
  if (ppw->scalar_ctx.idr_nature != idr_free_streaming) return 0.;

  const double shear_idr = (pv->index_pt_shear_idr >= 0)
                             ? y[pv->index_pt_shear_idr]
                             : ppw->scalar_ctx.shear_idr;

  return 4. / 3. * pvecback[index_bg_rho_] * shear_idr;
}

// ── IDM_DRMD ───────────────────────────────────────────────────────────────

void IDM_DRMDSpecies::RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                                   const perturb_workspace* /*ppw*/,
                                                   int /*gauge*/) {
  class_define_index(pv->index_pt_delta_idm_drmd, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idm_drmd, _TRUE_, index_pt, 1);
}

void IDM_DRMDSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                     const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  dy[pv->index_pt_delta_idm_drmd] = -(y[pv->index_pt_theta_idm_drmd] + ctx.metric_continuity);
  dy[pv->index_pt_theta_idm_drmd] = -ctx.a_prime_over_a * y[pv->index_pt_theta_idm_drmd] + ctx.metric_euler;

  if (pba_.has_idr_drmd) {
    auto* bgm = ppaw.perturbations_module->GetBackgroundModule().get();
    const double* pvecback = ppw->pvecback;
    double Rint, csp2, Gint;
    class_call(bgm->background_idm_drmd(ctx.a, pvecback[bgm->index_bg_rho_idm_drmd_] / pvecback[bgm->index_bg_rho_idr_drmd_], &Rint, &csp2, &Gint), bgm->error_message_, bgm->error_message_);

    if (ppw->approx[ppw->index_ap_tca_idm_drmd] == (int)tca_idm_drmd_on) {
      double GdDelta = 3.*csp2*(ctx.a_prime_over_a*y[pv->index_pt_theta_idm_drmd] + ctx.k2*ctx.delta_idr_drmd/4.);
      dy[pv->index_pt_theta_idm_drmd] = 0.25*ctx.k2*ctx.delta_idr_drmd + ctx.metric_euler - GdDelta*Rint;
    } else {
      dy[pv->index_pt_theta_idm_drmd] += Gint*(ctx.theta_idr_drmd - y[pv->index_pt_theta_idm_drmd]);
    }
  }
}

// ── IDR_DRMD ───────────────────────────────────────────────────────────────

void IDR_DRMDSpecies::RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                                   const perturb_workspace* /*ppw*/,
                                                   int /*gauge*/) {
  class_define_index(pv->index_pt_delta_idr_drmd, _TRUE_, index_pt, 1);
  class_define_index(pv->index_pt_theta_idr_drmd, _TRUE_, index_pt, 1);
}

void IDR_DRMDSpecies::PerturbDerivs(double /*tau*/, const double* y, double* dy,
                                     const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx  = ppw->scalar_ctx;

  dy[pv->index_pt_delta_idr_drmd] = -4. / 3. * (y[pv->index_pt_theta_idr_drmd] + ctx.metric_continuity);
  dy[pv->index_pt_theta_idr_drmd] = 0.25 * ctx.k2 * y[pv->index_pt_delta_idr_drmd] + ctx.metric_euler;

  if (pba_.has_idm_drmd) {
    auto* bgm = ppaw.perturbations_module->GetBackgroundModule().get();
    const double* pvecback = ppw->pvecback;
    double Rint, csp2, Gint;
    class_call(bgm->background_idm_drmd(ctx.a, pvecback[bgm->index_bg_rho_idm_drmd_] / pvecback[bgm->index_bg_rho_idr_drmd_], &Rint, &csp2, &Gint), bgm->error_message_, bgm->error_message_);

    if (ppw->approx[ppw->index_ap_tca_idm_drmd] == (int)tca_idm_drmd_on) {
      double GdDelta = 3.*csp2*(ctx.a_prime_over_a*ctx.theta_idm_drmd + ctx.k2*y[pv->index_pt_delta_idr_drmd]/4.);
      dy[pv->index_pt_theta_idr_drmd] -= GdDelta*Rint;
    } else {
      dy[pv->index_pt_theta_idr_drmd] -= Gint*Rint*(y[pv->index_pt_theta_idr_drmd] - ctx.theta_idm_drmd);
    }
  }
}
