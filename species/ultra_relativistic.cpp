#include "ultra_relativistic.h"

#include <cmath>
#include <stdexcept>

#include "background.h"
#include "background_column_writer.h"
#include "background_module.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"

UltraRelativisticSpecies::UltraRelativisticSpecies(const background& pba, double omega0_ur)
    : BaseSpecies("UR", EnergyType::Radiation), Omega0_ur_(omega0_ur), H0_(pba.H0) {}

// ── Background ────────────────────────────────────────────────────────────────

void UltraRelativisticSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_, _TRUE_, index_bg, 1);
}

void UltraRelativisticSpecies::ComputeBackground(double a_rel,
                                                 const double* /*pvecback_B*/,
                                                 double* pvecback) {
  const double a4         = a_rel * a_rel * a_rel * a_rel;
  pvecback[index_bg_rho_] = Omega0_ur_ * H0_ * H0_ / a4;
}

double UltraRelativisticSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_];
}

double UltraRelativisticSpecies::P(const double* pvecback) const {
  return pvecback[index_bg_rho_] / 3.;
}

double UltraRelativisticSpecies::DpDloga(const double* pvecback) const {
  return -4. / 3. * pvecback[index_bg_rho_];
}

void UltraRelativisticSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_ur", 0.);
}

void UltraRelativisticSpecies::WriteBackgroundData(const double* pvecback,
                                                   BackgroundColumnWriter& w) const {
  w.Add("(.)rho_ur", Rho(pvecback));
}

// ── Perturbations ─────────────────────────────────────────────────────────────

void UltraRelativisticSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                           perturb_vector* /*pv*/,
                                                           const precision* ppr,
                                                           int& index_pt,
                                                           const perturb_workspace* ppw,
                                                           int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);

  layout.l_max         = ppr->l_max_ur; /* moved from perturb_vector_init pre-loop setup */
  const int l_max_full = layout.l_max;

  /* If radiation streaming approximation is on, no UR variables are integrated.
     Sentinel -1 signals RSA to Delta/Theta/DeltaP/RhoPlusPShear. */
  if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_on) {
    layout.idx_delta = -1;
    layout.idx_theta = -1;
    layout.idx_shear = -1;
    layout.idx_l3    = -1;
    layout.l_max     = -1;
    return;
  }

  layout.idx_delta = index_pt++;
  layout.idx_theta = index_pt++;
  layout.idx_shear = index_pt++;

  if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off) {
    /* Full hierarchy: l_max set by caller */
    layout.l_max = l_max_full;
    if (l_max_full >= 3) {
      layout.idx_l3  = index_pt;
      index_pt      += l_max_full - 2;  // l3, l4, ..., l_max
    }
  }
  else {
    /* UFA active: only delta, theta, shear are integrated.
       l_max in the layout is set to 2 to indicate UFA truncation. */
    layout.l_max  = 2;
    layout.idx_l3 = -1;
  }
}

void UltraRelativisticSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                             double tau,
                                             const double* y,
                                             double* dy,
                                             const perturb_parameters_and_workspace& ppaw) {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;
  const double k                  = ppaw.k;
  const double k2                 = k * k;
  const double cotKgen            = ctx.cotKgen;
  const double s2_squared         = ctx.s2_squared;
  const double a_prime_over_a     = ctx.a_prime_over_a;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double metric_ufa_class   = ctx.metric_ufa_class;

  const perturbs* ppt   = ppaw.perturbations_module->GetPerturbs();
  const background* pba = ppaw.perturbations_module->GetBackground();
  const precision* ppr  = ppaw.perturbations_module->GetPrecision();

  /* RSA active: UR quantities handled analytically, nothing to integrate */
  if (layout.idx_delta < 0)
    return;

  int l;
  const double three_ceff2 = ppt->three_ceff2_ur;
  const double three_cvis2 = ppt->three_cvis2_ur;

  const int idx_delta = layout.idx_delta;
  const int idx_theta = layout.idx_theta;
  const int idx_shear = layout.idx_shear;
  const int idx_l3    = layout.idx_l3;
  const int l_max     = layout.l_max;

  /* ── density ─────────────────────────────────────────────────────────── */
  dy[idx_delta] = -4. / 3. * (y[idx_theta] + metric_continuity) +
                  (1. - three_ceff2) * a_prime_over_a *
                      (y[idx_delta] + 4. * a_prime_over_a * y[idx_theta] / k2);

  /* ── velocity ────────────────────────────────────────────────────────── */
  dy[idx_theta] = k2 * (three_ceff2 * y[idx_delta] / 4. - s2_squared * y[idx_shear]) +
                  metric_euler - (1. - three_ceff2) * a_prime_over_a * y[idx_theta];

  if (ppw->approx[ppw->index_ap_ufa] == (int) ufa_off) {
    /* ── exact shear ──────────────────────────────────────────────────── */
    dy[idx_shear] = 0.5 * (8. / 15. * (y[idx_theta] + metric_shear) -
                           3. / 5. * k * s_l[3] / s_l[2] * y[idx_shear + 1] -
                           (1. - three_cvis2) * 8. / 15. * (y[idx_theta] + metric_shear));

    /* ── l = 3 ────────────────────────────────────────────────────────── */
    l          = 3;
    dy[idx_l3] = k / (2. * l + 1.) *
                 (l * 2. * s_l[l] * s_l[2] * y[idx_shear] - (l + 1.) * s_l[l + 1] * y[idx_l3 + 1]);

    /* ── l = 4 .. l_max - 1 ──────────────────────────────────────────── */
    for (l = 4; l < l_max; l++) {
      dy[idx_delta + l] = k / (2. * l + 1) *
                          (l * s_l[l] * y[idx_delta + l - 1] -
                           (l + 1.) * s_l[l + 1] * y[idx_delta + l + 1]);
    }

    /* ── l = l_max (truncation) ──────────────────────────────────────── */
    l                 = l_max;
    dy[idx_delta + l] = k * (s_l[l] * y[idx_delta + l - 1] - (1. + l) * cotKgen * y[idx_delta + l]);

    /* ── optional self-interaction (G_eff_ur) damping ─────────────────── */
    if (ppt->G_eff_ur > 0.) {
      const double a = ctx.a;
      double taudot  = pow(a, -4) * pow(pow(4. / 11., 1. / 3.) * pba->T_cmb * _k_B_, 5) *
                       pow(ppt->G_eff_ur / (1e12 * _eV_ * _eV_), 2) * (2. * _PI_ / _h_P_) / _c_ *
                       _Mpc_over_m_;
      taudot         = std::min(taudot, a_prime_over_a * 1e9);
      if (taudot > 0.) {
        const double alpha_RTA[5] = {0.40, 0.43, 0.46, 0.47, 0.48};
        for (l = 2; l <= l_max; l++) {
          int alpha_index    = std::min(4, l - 2);
          dy[idx_delta + l] -= alpha_RTA[alpha_index] * taudot * y[idx_delta + l];
        }
      }
    }
  }
  else {
    /* ── UFA: fluid approximation for shear only ──────────────────────── */
    const int method = ppr->ur_fluid_approximation;

    if (method == (int) ufa_mb) {
      dy[idx_shear] = -3. / tau * y[idx_shear] + 2. / 3. * (y[idx_theta] + metric_shear);
    }
    else if (method == (int) ufa_hu) {
      dy[idx_shear] = -3. * a_prime_over_a * y[idx_shear] + 2. / 3. * (y[idx_theta] + metric_shear);
    }
    else { /* ufa_CLASS (default) */
      dy[idx_shear] = -3. / tau * y[idx_shear] + 2. / 3. * (y[idx_theta] + metric_ufa_class);
    }

    /* optional G_eff_ur damping on shear */
    if (ppt->G_eff_ur > 0.) {
      const double a  = ctx.a;
      double taudot   = pow(a, -4) * pow(pow(4. / 11., 1. / 3.) * pba->T_cmb * _k_B_, 5) *
                        pow(ppt->G_eff_ur / (1e12 * _eV_ * _eV_), 2) * (2. * _PI_ / _h_P_) / _c_ *
                        _Mpc_over_m_;
      taudot          = std::min(taudot, a_prime_over_a * 1e9);
      dy[idx_shear]  -= 0.40 * taudot * y[idx_shear];
    }
  }
}

double UltraRelativisticSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                                       const perturb_vector* /*pv*/,
                                       const double* y,
                                       const double* /*pvecback*/,
                                       const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* RSA active (idx_delta < 0): UR delta is handled by perturb_rsa_delta_and_theta,
     which adds rsa_delta_ur to ppw->delta_rho directly. Return 0 to avoid double-counting. */
  return (layout.idx_delta >= 0) ? y[layout.idx_delta] : 0.;
}

double UltraRelativisticSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                                       const perturb_vector* /*pv*/,
                                       const double* y,
                                       const double* /*pvecback*/,
                                       const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* RSA active (idx_theta < 0): UR theta is handled by perturb_rsa_delta_and_theta. Return 0. */
  return (layout.idx_theta >= 0) ? y[layout.idx_theta] : 0.;
}

double UltraRelativisticSpecies::DeltaP(const BaseSpecies::PerturbLayout& base,
                                        const perturb_vector* /*pv*/,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* delta_p = c_s^2 * delta_rho = (1/3) * rho * delta */
  return (layout.idx_delta >= 0) ? Rho(pvecback) * y[layout.idx_delta] / 3. : 0.;
}

double UltraRelativisticSpecies::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                               const perturb_vector* /*pv*/,
                                               const double* y,
                                               const double* pvecback,
                                               const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  /* (rho + p) * sigma = 4/3 * rho_ur * shear_ur */
  return (layout.idx_shear >= 0) ? 4. / 3. * Rho(pvecback) * y[layout.idx_shear] : 0.;
}

void UltraRelativisticSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                                      double* y,
                                                      const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_delta >= 0)
    y[layout.idx_delta] = ctx.delta_ur;
  if (layout.idx_theta >= 0)
    y[layout.idx_theta] = ctx.theta_ur;
  if (layout.idx_shear >= 0)
    y[layout.idx_shear] = ctx.shear_ur;
  if (layout.idx_l3 >= 0)
    y[layout.idx_l3] = ctx.l3_ur;
}

void UltraRelativisticSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                                           const double* y,
                                           const double* /*dy*/,
                                           PerturbSourceContext& ctx) {
  const auto& layout         = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;

  const double a_prime_over_a = ctx.a_prime_over_a;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_ur ───────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_ur_ == _TRUE_) {
    const double delta_ur = (layout.idx_delta >= 0) ? y[layout.idx_delta] : ppw->rsa_delta_ur;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_ur_,
                          ctx.index_tau,
                          ctx.index_k,
                          delta_ur +
                              4. * a_prime_over_a * ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_ur ───────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_ur_ == _TRUE_) {
    const double theta_ur = (layout.idx_theta >= 0) ? y[layout.idx_theta] : ppw->rsa_theta_ur;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_ur_,
                          ctx.index_tau,
                          ctx.index_k,
                          theta_ur + ctx.theta_shift);  // N-body gauge correction
  }
}

// ── Newtonian-gauge transform ─────────────────────────────────────────────────

void UltraRelativisticSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                             double* y,
                                                             const PerturbIcContext& ctx) {
  // delta/theta shift; shear (idx_shear) and l3 (idx_l3) are gauge-invariant.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback, ctx);
}

void UltraRelativisticSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y,
    double* new_y,
    const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);

  /* Guard: only copy if both layouts have slots registered (not RSA mode) */
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;

  /* delta, theta, shear always present in both layouts */
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
  new_y[new_l.idx_shear] = old_y[old_l.idx_shear];

  /* l3+ only present when ufa is off (full hierarchy on both sides) */
  if (old_l.l_max >= 3 && new_l.l_max >= 3 && old_l.idx_l3 >= 0 && new_l.idx_l3 >= 0) {
    new_y[new_l.idx_l3] = old_y[old_l.idx_l3];
    for (int l = 4; l <= new_l.l_max; ++l)
      new_y[new_l.idx_delta + l] = old_y[old_l.idx_delta + l];
  }
}

void UltraRelativisticSpecies::WriteOutputColumns(
    PerturbColumnWriter& w,
    const PerturbationsModule& mod,
    enum file_format fmt,
    BaseSpecies::TransferColumnSection section) const {
  const background* pba = mod.GetBackground();
  if (fmt == class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
      w.Add("d_ur", mod.index_tp_delta_ur_, _TRUE_);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
      w.Add("t_ur", mod.index_tp_theta_ur_, _TRUE_);
  }
  else if (fmt == camb_format) {
    if (section != TransferColumnSection::velocity)
      w.Add("-T_ur/k2", mod.index_tp_delta_ur_, _TRUE_);
  }
}

void UltraRelativisticSpecies::PrintVariables(PerturbColumnWriter& w,
                                              double /*tau*/,
                                              const double* y,
                                              const PerturbationsModule& mod,
                                              const perturb_workspace* ppw) const {
  double delta_ur = 0., theta_ur = 0., shear_ur = 0.;

  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    const double k           = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback;
    const double* pvecmetric = ppw->pvecmetric;

    if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) {
      const auto& lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
      delta_ur        = y[lay.idx_delta];
      theta_ur        = y[lay.idx_theta];
      shear_ur        = y[lay.idx_shear];
    }
    else {
      delta_ur = ppw->rsa_delta_ur;
      theta_ur = ppw->rsa_theta_ur;
      shear_ur = 0.;
    }

    // Synchronous → Newtonian gauge conversion (relativistic: factor 4)
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_ur           -= 4. * H * a * alpha;
      theta_ur           += k * k * alpha;
    }
  }

  w.Add("delta_ur", delta_ur, true);
  w.Add("theta_ur", theta_ur, true);
  w.Add("shear_ur", shear_ur, true);
}

// ── MarkUsedInSources ─────────────────────────────────────────────────────────

void UltraRelativisticSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                                 const perturb_workspace* ppw,
                                                 int* used_in_sources) const {
  const auto& ur_lay = static_cast<const PerturbLayout&>(base);
  /* idx_l3 < 0 when RSA is on, UFA is on, or we are in tensor/vector mode
     (UR has no tensor slots). index_ap_ufa is only defined in scalar mode,
     so we must not access it when idx_l3 < 0. */
  if (ur_lay.idx_l3 < 0)
    return;
  /* UR l>=3 multipoles not needed in sources when both rsa and ufa are off.
     (The idx_l3 >= 0 guard above already implies rsa_off and ufa_off, but we
     check explicitly for clarity and to mirror the legacy module condition.) */
  if (ppw->approx[ppw->index_ap_rsa] != (int) rsa_off)
    return;
  if (ppw->approx[ppw->index_ap_ufa] != (int) ufa_off)
    return;

  for (int idx = ur_lay.idx_l3; idx <= ur_lay.idx_delta + ur_lay.l_max; ++idx)
    used_in_sources[idx] = _FALSE_;
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> UltraRelativisticSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // Read the UR density from the input file, supporting four equivalent forms:
  //   N_ur   (or deprecated N_eff) — effective neutrino number (massless)
  //   Omega_ur                     — physical density parameter
  //   omega_ur                     — reduced density (omega = Omega * h^2)
  // At most one of the three forms may be specified; if none is given, the
  // standard default N_ur = 3.046 is used.

  // (a) Try to read N_ur (and deprecated N_eff alias).
  double n_ur    = 0.;
  double n_eff   = 0.;
  bool flag_nur  = ctx.pfc->read_double("N_ur", n_ur);
  bool flag_neff = ctx.pfc->read_double("N_eff", n_eff);

  // N_eff is a deprecated synonym for N_ur — only one may be present.
  if (flag_nur && flag_neff)
    throw std::invalid_argument(
        "In input file, you can only enter one of N_eff (deprecated syntax) or N_ur "
        "(up-to-date syntax), since they both describe the same, i.e. the contribution "
        "of ultra-relativistic species to the effective neutrino number");
  // Normalise: treat N_eff as N_ur going forward.
  if (flag_neff) {
    n_ur      = n_eff;
    flag_nur  = true;
    flag_neff = false;
  }

  // (b) Try to read Omega_ur and omega_ur.
  double omega_ur_big   = 0.;  // Omega_ur
  double omega_ur_small = 0.;  // omega_ur (= Omega_ur * h^2)
  bool flag_Omega       = ctx.pfc->read_double("Omega_ur", omega_ur_big);
  bool flag_omega       = ctx.pfc->read_double("omega_ur", omega_ur_small);

  // At most one of the three forms may be specified.
  const int n_specified = (flag_nur ? 1 : 0) + (flag_Omega ? 1 : 0) + (flag_omega ? 1 : 0);
  if (n_specified >= 2)
    throw std::invalid_argument(
        "In input file, you can only enter one of N_ur (or deprecated N_eff), Omega_ur, or "
        "omega_ur, choose one");

  // (c) Compute Omega0_ur from whichever form was provided.
  double omega0_ur      = 0.;
  const double omega0_g = ctx.pba->Omega0_g;

  if (n_specified == 0) {
    // Default: standard 3.046 massless neutrinos.
    omega0_ur = 3.046 * 7. / 8. * pow(4. / 11., 4. / 3.) * omega0_g;
  }
  else if (flag_nur) {
    omega0_ur = n_ur * 7. / 8. * pow(4. / 11., 4. / 3.) * omega0_g;
  }
  else if (flag_Omega) {
    omega0_ur = omega_ur_big;
  }
  else {
    omega0_ur = omega_ur_small / ctx.pba->h / ctx.pba->h;
  }

  if (omega0_ur == 0.)
    return result;

  result.push_back({"UR", std::make_unique<UltraRelativisticSpecies>(*ctx.pba, omega0_ur)});
  return result;
}
