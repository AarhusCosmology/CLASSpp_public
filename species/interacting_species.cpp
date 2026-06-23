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

// ── RegisterPerturbationIndices ──────────────

void IDM_DRSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                perturb_vector* /*pv*/,
                                                const precision* /*ppr*/,
                                                int& index_pt,
                                                const perturb_workspace* /*ppw*/,
                                                int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = index_pt;
  ++index_pt;

  layout.idx_theta = index_pt;
  ++index_pt;
}

// ── PerturbDerivs ────────────────────────────

void IDM_DRSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                  double /*tau*/,
                                  const double* y,
                                  double* dy,
                                  const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  // idm_dr density: same as CDM continuity equation
  dy[layout.idx_delta] = -(y[layout.idx_theta] + ctx.metric_continuity);

  // idm_dr velocity: Hubble friction only (coupling to IDR added by IDM_DR_IDR_Species)
  dy[layout.idx_theta] = -ctx.a_prime_over_a * y[layout.idx_theta] + ctx.metric_euler;
}

// ── Stress-energy observables ──────────────────────────────────────────────

BaseSpecies::StressEnergyContribution IDM_DRSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  // IDM_DR is cold (P = 0, delta_p = 0, shear = 0)
  const double rho = pvecback[index_bg_rho_];
  return {rho,
          /*p=*/0.,
          /*delta_rho=*/rho * y[layout.idx_delta],
          /*rho_plus_p_theta=*/rho * y[layout.idx_theta],
          /*delta_p=*/0.,
          /*rho_plus_p_shear=*/0.};
}

// ── PerturbSynchronousToNewtonian ─────────────────────────────────────────

void IDM_DRSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                  double* y,
                                                  const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

// ── CopyPerturbationsAcrossSwitch ─────────────────────────────────────────

void IDM_DRSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                  const BaseSpecies::PerturbLayout& new_base,
                                                  const double* old_y,
                                                  double* new_y,
                                                  const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}

// ── IDR ───────────────────────────────────────────────────────────────────

// ── RegisterPerturbationIndices ───────────────

void IDRSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                             perturb_vector* pv,
                                             const precision* /*ppr*/,
                                             int& index_pt,
                                             const perturb_workspace* ppw,
                                             int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  /* Initialize layout slots to sentinel -1 */
  layout.idx_delta = -1;
  layout.idx_theta = -1;
  layout.idx_shear = -1;
  layout.idx_l3    = -1;
  layout.l_max     = -1;

  /* RSA active: IDR handled analytically, nothing to integrate. */
  if (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)
    return;

  layout.idx_delta = index_pt;
  ++index_pt;

  layout.idx_theta = index_pt;
  ++index_pt;

  /* Full hierarchy: only register shear/l3 if idr_nature == idr_free_streaming
     AND (no IDM_DR OR TCA is off) */
  if (ppw->scalar_ctx.idr_nature == idr_free_streaming) {
    if (!has_sibling_idm_dr_ || ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_off) {
      layout.idx_shear = index_pt;
      ++index_pt;

      layout.l_max = l_max_idr_;

      if (l_max_idr_ >= 3) {
        layout.idx_l3  = index_pt;
        index_pt      += l_max_idr_ - 2;  // l3, l4, ..., l_max
      }
    }
  }
}

// ── PerturbDerivs ────────────────────────────

void IDRSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                               double /*tau*/,
                               const double* y,
                               double* dy,
                               const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  /* RSA active (idx_delta < 0): IDR handled analytically. */
  if (layout.idx_delta < 0)
    return;

  const double* s_l    = ppw->s_l.data();
  const double k       = ctx.k;
  const double cotKgen = ctx.cotKgen;

  const int idx_delta = layout.idx_delta;
  const int idx_theta = layout.idx_theta;
  const int idx_shear = layout.idx_shear;
  const int idx_l3    = layout.idx_l3;
  const int l_max     = layout.l_max;

  // idr density
  dy[idx_delta] = -4. / 3. * (y[idx_theta] + ctx.metric_continuity);

  if (ctx.idr_nature == idr_free_streaming) {
    // idr velocity and hierarchy (only when shear is registered; under TCA shear is
    // not registered and AddCouplingDerivs will overwrite theta_idr)
    if (idx_shear >= 0) {
      dy[idx_theta] = ctx.k2 * (y[idx_delta] / 4. - ctx.s2_squared * y[idx_shear]) +
                      ctx.metric_euler;

      // idr shear (l=2)
      dy[idx_shear] = 0.5 * (8. / 15. * (y[idx_theta] + ctx.metric_shear) -
                             3. / 5. * k * s_l[3] / s_l[2] * y[idx_shear + 1]);

      // idr l=3
      int l      = 3;
      dy[idx_l3] = k / (2. * l + 1.) *
                   (l * 2. * s_l[l] * s_l[2] * y[idx_shear] -
                    (l + 1.) * s_l[l + 1] * y[idx_l3 + 1]);

      // idr l>3
      for (l = 4; l < l_max; l++)
        dy[idx_delta + l] = k / (2. * l + 1) *
                            (l * s_l[l] * y[idx_delta + l - 1] -
                             (l + 1.) * s_l[l + 1] * y[idx_delta + l + 1]);

      // idr lmax
      l                 = l_max;
      dy[idx_delta + l] = k *
                          (s_l[l] * y[idx_delta + l - 1] - (1. + l) * cotKgen * y[idx_delta + l]);
    }
    // else: TCA-on; theta_idr will be set by IDM_DR_IDR_Species::AddCouplingDerivs
  }
  else {
    // Fluid idr velocity
    dy[idx_theta] = ctx.k2 / 4. * y[idx_delta] + ctx.metric_euler;
  }
}

// ── Stress-energy observables ──────────────────────────────────────────────

BaseSpecies::StressEnergyContribution IDRSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* RSA active when idx_delta < 0: delta/theta handled analytically; return
     zero contributions to avoid double-counting. */
  const double rho = Rho(pvecback);  // pvecback[index_bg_rho_]
  const double p   = rho / 3.;       // P = rho/3 for radiation

  StressEnergyContribution sec;
  sec.rho = rho;
  sec.p   = p;

  if (layout.idx_delta >= 0) {
    sec.delta_rho = rho * y[layout.idx_delta];
    sec.delta_p   = rho * y[layout.idx_delta] / 3.;
  }
  if (layout.idx_theta >= 0) {
    sec.rho_plus_p_theta = (rho + p) * y[layout.idx_theta];
  }

  if (ppw->scalar_ctx.idr_nature == idr_free_streaming) {
    /* Use layout slot if shear is registered; fall back to TCA formula otherwise
       (matches RhoPlusPShear logic exactly). */
    const double shear_idr = (layout.idx_shear >= 0) ? y[layout.idx_shear]
                                                     : TcaShearIdr(layout, y, ppw);
    sec.rho_plus_p_shear   = 4. / 3. * rho * shear_idr;
  }
  /* else: fluid idr — rho_plus_p_shear stays 0. */

  return sec;
}

// ── PerturbSynchronousToNewtonian ──────────────────────────────────────────

void IDRSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                               double* y,
                                               const PerturbIcContext& ctx) {
  // delta/theta shift; shear/l3 gauge-invariant. idm_dr's synchronous theta is
  // already theta_ur (tight-coupling lock), so the universal += reproduces it.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

// ── CopyPerturbationsAcrossSwitch ──────────────────────────────────────────

void IDRSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                               const BaseSpecies::PerturbLayout& new_base,
                                               const double* old_y,
                                               double* new_y,
                                               const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);

  /* Guard: only copy if both layouts have slots (rsa_idr must be off in both). */
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;

  /* delta and theta are always present when rsa_idr is off. */
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];

  /* shear and l3+ are only present in free-streaming mode with tca_idm_dr off.
     Copy slot-by-slot only when both sides have them registered. */
  if (old_l.idx_shear >= 0 && new_l.idx_shear >= 0) {
    new_y[new_l.idx_shear] = old_y[old_l.idx_shear];
  }
  if (old_l.idx_l3 >= 0 && new_l.idx_l3 >= 0 && old_l.l_max >= 3 && new_l.l_max >= 3) {
    new_y[new_l.idx_l3] = old_y[old_l.idx_l3];
    for (int l = 4; l <= new_l.l_max; ++l)
      new_y[new_l.idx_delta + l] = old_y[old_l.idx_delta + l];
  }
}

// ── TcaShearIdr ────────────────────────────────────────────────────────────

double IDRSpecies::TcaShearIdr(const PerturbLayout& layout,
                               const double* y,
                               const perturb_workspace* ppw) const {
  // TCA shear is meaningful when shear_idr is NOT in the y-vector (i.e. we
  // didn't register it because TCA was on at vector-init). The layout check
  // is the canonical signal — we deliberately don't gate on
  // ppw->approx[index_ap_tca_idm_dr] because perturb_vector_init calls this
  // during the tca_on -> tca_off transition (after the approx flag has flipped
  // but with the old layout still pointing at the old slot set), and we need
  // the last TCA prediction to seed the new shear_idr integration.
  if (ppw->scalar_ctx.idr_nature != idr_free_streaming)
    return 0.;
  if (layout.idx_shear >= 0)
    return 0.;
  if (ppw->approx[ppw->index_ap_rsa_idr] != (int) rsa_idr_off)
    return 0.;
  if (!has_sibling_idm_dr_)
    return 0.;
  if (ppt_->gauge != possible_gauges::newtonian)
    return 0.;  // synchronous gauge derives shear in perturb_einstein

  return 0.5 * (8. / 15. / ppw->pvecthermo[thm_->index_th_dmu_idm_dr_] / alpha_idm_dr_[0] *
                y[layout.idx_theta]);
}

// ── IDM_DRMD ───────────────────────────────────────────────────────────────

// ── RegisterPerturbationIndices ──────────────

void IDM_DRMDSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                  perturb_vector* /*pv*/,
                                                  const precision* /*ppr*/,
                                                  int& index_pt,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = index_pt;
  ++index_pt;

  layout.idx_theta = index_pt;
  ++index_pt;
}

// ── PerturbDerivs ────────────────────────────

void IDM_DRMDSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                    double /*tau*/,
                                    const double* y,
                                    double* dy,
                                    const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[layout.idx_delta] = -(y[layout.idx_theta] + ctx.metric_continuity);
  // Free-streaming velocity: coupling to IDR_DRMD added by IDM_DRMD_IDR_DRMD_Species
  dy[layout.idx_theta] = -ctx.a_prime_over_a * y[layout.idx_theta] + ctx.metric_euler;
}

// ── Stress-energy observables ──────────────────────────────────────────────

BaseSpecies::StressEnergyContribution IDM_DRMDSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  // IDM_DRMD is cold (P = 0, delta_p = 0, shear = 0)
  const double rho = pvecback[index_bg_rho_];
  return {rho,
          /*p=*/0.,
          /*delta_rho=*/rho * y[layout.idx_delta],
          /*rho_plus_p_theta=*/rho * y[layout.idx_theta],
          /*delta_p=*/0.,
          /*rho_plus_p_shear=*/0.};
}

// ── PerturbSynchronousToNewtonian ─────────────────────────────────────────

void IDM_DRMDSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

// ── CopyPerturbationsAcrossSwitch ─────────────────────────────────────────

void IDM_DRMDSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                    const BaseSpecies::PerturbLayout& new_base,
                                                    const double* old_y,
                                                    double* new_y,
                                                    const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}

// ── IDR_DRMD ───────────────────────────────────────────────────────────────

// ── RegisterPerturbationIndices ──────────────

void IDR_DRMDSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                  perturb_vector* /*pv*/,
                                                  const precision* /*ppr*/,
                                                  int& index_pt,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.idx_delta = index_pt;
  ++index_pt;

  layout.idx_theta = index_pt;
  ++index_pt;
}

// ── PerturbDerivs ────────────────────────────

void IDR_DRMDSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                    double /*tau*/,
                                    const double* y,
                                    double* dy,
                                    const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[layout.idx_delta] = -4. / 3. * (y[layout.idx_theta] + ctx.metric_continuity);
  // Free-streaming velocity: coupling to IDM_DRMD added by IDM_DRMD_IDR_DRMD_Species
  dy[layout.idx_theta] = 0.25 * ctx.k2 * y[layout.idx_delta] + ctx.metric_euler;
}

// ── Stress-energy observables ──────────────────────────────────────────────

BaseSpecies::StressEnergyContribution IDR_DRMDSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  // IDR_DRMD: radiation (p = rho/3), no shear term
  const double rho = pvecback[index_bg_rho_];
  const double p   = rho / 3.;

  StressEnergyContribution sec;
  sec.rho = rho;
  sec.p   = p;
  if (layout.idx_delta >= 0) {
    sec.delta_rho = rho * y[layout.idx_delta];
    sec.delta_p   = rho * y[layout.idx_delta] / 3.;
  }
  if (layout.idx_theta >= 0) {
    sec.rho_plus_p_theta = (rho + p) * y[layout.idx_theta];
  }
  /* rho_plus_p_shear = 0. always (no shear slot registered) */
  return sec;
}

// ── PerturbSynchronousToNewtonian ─────────────────────────────────────────

void IDR_DRMDSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

// ── CopyPerturbationsAcrossSwitch ─────────────────────────────────────────

void IDR_DRMDSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                    const BaseSpecies::PerturbLayout& new_base,
                                                    const double* old_y,
                                                    double* new_y,
                                                    const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}
