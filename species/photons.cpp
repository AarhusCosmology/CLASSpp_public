#include "photons.h"

#include <string>

#include "background_module.h"
#include "errors.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "photons_optimal.h"
#include "photons_tam.h"
#include "thermodynamics_module.h"

/* Leading order in eps = 1/kappa' of the tensor hierarchy, in the optimal basis
 * (PhotonsOptimalSpecies::PerturbTensorDerivs).
 *
 * F_0 and G_0 are O(eps); theta is O(eps^2); shear and G_2 are O(eps^3); F_4 and
 * G_4 are O(eps^5). Keeping only the O(1) terms of the two monopole equations,
 *
 *   kappa' (F_0 + sqrt6 P_2) = sqrt6 hdot        (temperature, sourced by the GW)
 *   kappa' (G_0 - sqrt6 P_2) = 0                 (polarisation, unsourced)
 *
 * and P_2 = -(1/10 F_0 - 3/5 G_0)/sqrt6 once the higher-order multipoles drop
 * out. The second equation gives G_0 = sqrt6 P_2, which closes P_2 = -F_0/4;
 * substituting into the first gives (3/4) F_0 = sqrt6 hdot.
 *
 * Same result as class_public v3.3.0 (credits C. Pitrou), reached independently
 * here; CLASS++ carried the pre-v3.3.0 expressions until #398.
 *
 * First order. F_0 and G_0 track hdot/kappa', so they vary on the timescale over
 * which tau_c itself varies -- and through recombination tau_c climbs by orders
 * of magnitude, so kappa'^-1 F_0' is NOT negligible. Restoring the two dropped
 * derivative terms (the -4/3 theta term stays out: it is O(eps^3), as are the
 * shear and G_2 that the multipole combination omits),
 *
 *   kappa' (F_0 + sqrt6 P_2) = sqrt6 hdot - F_0'
 *   kappa' (G_0 - sqrt6 P_2) = -G_0'
 *
 * with F_0' = (4 sqrt6/3) D and G_0' = -(sqrt6/3) D evaluated on the leading
 * solution, D = d/dtau (hdot/kappa'). Solving the same 3x3 system gives the
 * coefficients below. G_0 = -F_0/4 holds only at leading order: the drift enters
 * F_0 and G_0 with -22/9 against +13/9, and that asymmetry is the correction.
 * Not in class_public. */
PhotonsSpecies::TensorTcaClosure PhotonsSpecies::TensorTightCoupling(double gwdot,
                                                                     double dkappa,
                                                                     double drift) {
  const double u = gwdot / dkappa; /* hdot/kappa' */
  const double v = drift / dkappa; /* D/kappa' */
  return {_SQRT6_ * (4. / 3. * u - 22. / 9. * v),
          _SQRT6_ * (-1. / 3. * u + 13. / 9. * v),
          -1. / 3. * u + 10. / 9. * v};
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
  if (layout.idx_theta >= 0) {
    se.rho_plus_p_theta = 4. / 3. * rho_g * y[layout.idx_theta];
  }
  if (layout.idx_shear >= 0) {
    se.rho_plus_p_shear = 4. / 3. * rho_g * y[layout.idx_shear];  // no approximation
  }
  else if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_off)) {
    se.rho_plus_p_shear = 0.;  // RSA: shear neglected
  }
  else if (ppt_->gauge == possible_gauges::newtonian) {  // TCA, 1st order
    se.rho_plus_p_shear = 4. / 3. * rho_g *
                          (16. / 45. / ppw->pvecthermo[thm_->index_th_dkappa_] *
                           y[layout.idx_theta]);
  }
  // else synchronous TCA: 0 here, re-set later in perturb_einstein
  return se;
}

// ── ApplyInitialConditions ────────────────────────────────────────────────────

void PhotonsSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                            double* y,
                                            const PerturbIcContext& ctx) {
  const auto& layout             = static_cast<const PerturbLayout&>(base);
  const PerturbationsModule* mod = ctx.p_mod;
  if (layout.idx_delta < 0 || layout.idx_theta < 0) {
    return;
  }

  // Every initial condition gives the photons the same delta_g_ic and theta_g_ic.
  if (ctx.index_ic == mod->index_ic_ad_ || ctx.index_ic == mod->index_ic_cdi_ ||
      ctx.index_ic == mod->index_ic_bi_ || ctx.index_ic == mod->index_ic_nid_ ||
      ctx.index_ic == mod->index_ic_niv_) {
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

  if (ctx.index_md != p_mod->index_md_scalars_) {
    return;
  }

  // Convenient lambda to write into source table
  auto set_source = [&](int index_tp, double value) {
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp, ctx.index_tau, ctx.index_k, value);
  };

  double delta_g;
  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_on)) {
    delta_g = ppw->rsa_delta_g;
  }
  else {
    delta_g = y[layout.idx_delta];
  }

  if (index_tp_delta_ >= 0) {
    set_source(index_tp_delta_, delta_g + 4. * ctx.a_prime_over_a * ctx.theta_over_k2);
  }

  if (index_tp_theta_ >= 0) {
    double theta_g;
    if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_off)) {
      theta_g = y[layout.idx_theta];
    }
    else {
      theta_g = ppw->rsa_theta_g;
    }
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

// ── Output ────────────────────────────────────────────────────────────────────

void PhotonsSpecies::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void PhotonsSpecies::AddTransferOutput(OutputTable& row,
                                       const double* tk,
                                       const PerturbationsModule& mod,
                                       file_format fmt,
                                       BaseSpecies::TransferColumnSection section) const {
  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers) {
      if (index_tp_delta_ >= 0) {
        row.Set("d_g", tk[index_tp_delta_]);
      }
    }
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers) {
      if (index_tp_theta_ >= 0) {
        row.Set("t_g", tk[index_tp_theta_]);
      }
    }
  }
  // file_format::camb_format: photons not written separately
}

void PhotonsSpecies::AddPerturbationOutput(OutputTable& row,
                                           const BaseSpecies::PerturbLayout& base,
                                           double /*tau*/,
                                           const double* y,
                                           const PerturbationsModule& mod,
                                           const perturb_workspace* ppw) const {
  double delta_g = 0.;
  double theta_g = 0.;
  double shear_g = 0.;

  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  // scalar_ctx.k is always current here: the evolver sets it in perturb_derivs
  // before the print callback fires, so it equals pppaw->k at this call site.
  const double k           = ppw->scalar_ctx.k;
  const double* pvecthermo = ppw->pvecthermo.data();
  const double* pvecback   = ppw->pvecback.data();
  const double* pvecmetric = ppw->pvecmetric.data();
  const double dkappa      = pvecthermo[mod.GetThermodynamicsModule()->index_th_dkappa_];

  if (ppw->approx[ppw->index_ap_rsa] == static_cast<int>(rsa_off)) {
    delta_g = y[g_lay.idx_delta];
    theta_g = y[g_lay.idx_theta];
    if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
      shear_g = ppw->tca_shear_g;
    }
    else {
      shear_g = y[g_lay.idx_shear];
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

  row.Set("delta_g", delta_g);
  row.Set("theta_g", theta_g);
  row.Set("shear_g", shear_g);
  AddPolarisationOutput(row, g_lay, y, k, dkappa, ppw);
}

// ── Scalar derivatives both bases share ──────────────────────────────────────

bool PhotonsSpecies::DensityAndVelocityDerivs(const PerturbLayout& layout,
                                              const double* y,
                                              double* dy,
                                              const perturb_parameters_and_workspace& ppaw) const {
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const double k2                = ppaw.k * ppaw.k;
  const double metric_continuity = ctx.metric_continuity;
  const double metric_euler      = ctx.metric_euler;
  const double a_prime_over_a    = ctx.a_prime_over_a;
  const double s2_squared        = ctx.s2_squared;
  const double R                 = ctx.R;

  const double dkappa =
      ppw->pvecthermo[ppaw.perturbations_module->GetThermodynamicsModule()->index_th_dkappa_];

  const double delta_g = y[layout.idx_delta];
  const double theta_g = y[layout.idx_theta];
  const double theta_b = ctx.theta_b;

  /* delta equation (same whether TCA on or off) */
  dy[layout.idx_delta] = -4. / 3. * (theta_g + metric_continuity);

  if (ppw->approx[ppw->index_ap_tca] == static_cast<int>(tca_on)) {
    /* TCA: photon theta sourced by baryon theta derivative (already set by BaryonsSpecies).
       In TCA, perturbed recombination is off, so cb2 * delta_b is used. */
    const double delta_p_b_over_rho_b = ppw->pvecthermo[thm_->index_th_cb2_] * ctx.delta_b;
    dy[layout.idx_theta] =
        -(dy[ctx.b_idx_theta] + a_prime_over_a * theta_b - k2 * delta_p_b_over_rho_b) / R +
        k2 * (0.25 * delta_g - s2_squared * ppw->tca_shear_g) + (1. + R) / R * metric_euler;
    return true;
  }

  dy[layout.idx_theta] = k2 * (delta_g / 4. - s2_squared * y[layout.idx_shear]) + metric_euler +
                         dkappa * (theta_b - theta_g);
  return false;
}

// ── Factory ───────────────────────────────────────────────────────────────────

/* `hierarchy` picks the basis of the photon hierarchy; exactly one photon species is built.
   The key is class_public's (v3.4.0), where it also switches ur; here the massless species
   keep the optimal basis, which is exact for unpolarised radiation even in curved space. */
std::vector<Named> PhotonsSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::string hierarchy = "optimal";
  if (ctx.pfc != nullptr) {
    hierarchy = ctx.pfc->get<std::string>("hierarchy").value_or("optimal");
  }
  std::vector<Named> result;
  if (hierarchy == "optimal") {
    result.push_back({"Photons", std::make_unique<PhotonsOptimalSpecies>(*ctx.pba)});
  }
  else if (hierarchy == "tam" || hierarchy == "TAM") {
    result.push_back({"Photons", std::make_unique<PhotonsTamSpecies>(*ctx.pba)});
  }
  else {
    class_stop_severe("hierarchy = '%s': it must be 'optimal' or 'tam'", hierarchy.c_str());
  }
  return result;
}
