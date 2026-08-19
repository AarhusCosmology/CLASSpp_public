#include "dcdm.h"

#include "background.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

DCDMSpecies::DCDMSpecies(const background& pba,
                         double omega0_dcdmdr,
                         double Gamma_dcdm,
                         double Omega_ini_dcdm)
    : BaseSpecies("DCDM", EnergyType::Matter), pba_(pba), Omega0_dcdmdr_(omega0_dcdmdr),
      Gamma_dcdm_(Gamma_dcdm), Omega_ini_dcdm_(Omega_ini_dcdm) {}

void DCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_dcdm_, true, index_bg, 1);
  index_bg_rho_ = index_bg_rho_dcdm_;
}

void DCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_rho_dcdm_, true, index_bi, 1);
}

void DCDMSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  ctx.pvecback_integration[index_bi_rho_dcdm_] = Omega_ini_dcdm_ * std::pow(pba_.H0, 2) *
                                                 std::pow(1.0 / ctx.a_ini, 3);
}

void DCDMSpecies::ComputeBackground(double /*a*/, const double* pvecback_B, double* pvecback) {
  pvecback[index_bg_rho_dcdm_] = pvecback_B[index_bi_rho_dcdm_];
}

void DCDMSpecies::BackgroundDerivs(double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const double* pvecback) {
  const double a   = pvecback[bgm_->index_bg_a_];
  const double H   = pvecback[bgm_->index_bg_H_];
  const double rho = y[index_bi_rho_dcdm_];
  /** rho' = -a*(3H + Gamma) * rho */
  dy[index_bi_rho_dcdm_] = -a * (3. * H + Gamma_dcdm_) * rho;
}

double DCDMSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_dcdm_];
}

double DCDMSpecies::P(const double* /*pvecback*/) const {
  return 0.;
}

double DCDMSpecies::RhoPrimeOverRho(const double* pvecback, double a_prime_over_a) const {
  // ρ̇/ρ = -a(3H + Γ) = -3ℋ - aΓ  (same physics as BackgroundDerivs).
  const double a = pvecback[bgm_->index_bg_a_];
  return -3. * a_prime_over_a - a * Gamma_dcdm_;
}

void DCDMSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                              perturb_vector* pv,
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

void DCDMSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                double /*tau*/,
                                const double* y,
                                double* dy,
                                const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  dy[layout.idx_delta] = -(y[layout.idx_theta] + ctx.metric_continuity) -
                         ctx.a * Gamma_dcdm_ / ctx.k2 * ctx.metric_euler;
  dy[layout.idx_theta] = -ctx.a_prime_over_a * y[layout.idx_theta] + ctx.metric_euler;
}

BaseSpecies::StressEnergyContribution DCDMSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho              = pvecback[index_bg_rho_dcdm_];
  se.p                = 0.;
  se.delta_rho        = pvecback[index_bg_rho_] * y[layout.idx_delta];
  se.rho_plus_p_theta = pvecback[index_bg_rho_] * y[layout.idx_theta];
  se.delta_p          = 0.;
  se.rho_plus_p_shear = 0.;
  return se;
}

void DCDMSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                         double* y,
                                         const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.idx_delta >= 0)
    y[layout.idx_delta] = 3. / 4. * ctx.delta_g_ic;
}

void DCDMSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                double* y,
                                                const PerturbIcContext& ctx) {
  // delta += (ρ̇/ρ)·α with ρ̇/ρ = -3ℋ - aΓ (decay-aware RhoPrimeOverRho), theta += k²·α.
  const auto& l = static_cast<const PerturbLayout&>(base);
  ApplyFluidLikeNewtonianShift(y, l.idx_delta, l.idx_theta, ctx.ppw->pvecback.data(), ctx);
}

void DCDMSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                const BaseSpecies::PerturbLayout& new_base,
                                                const double* old_y,
                                                double* new_y,
                                                const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta >= 0 && new_l.idx_delta >= 0)
    new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  if (old_l.idx_theta >= 0 && new_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}

void DCDMSpecies::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_, ctx.wants_velocity, index_tp, 1);
}

void DCDMSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                              const double* y,
                              const double* /*dy*/,
                              PerturbSourceContext& ctx) const {
  PerturbationsModule* p_mod = ctx.p_mod;
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const auto& layout     = static_cast<const PerturbLayout&>(base);
  const double* pvecback = ctx.ppw->pvecback.data();

  if (index_tp_delta_ >= 0) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_delta_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[layout.idx_delta] - RhoPrimeOverRho(pvecback, ctx.a_prime_over_a) *
                                                    ctx.theta_over_k2);  // N-body gauge correction
  }
  if (index_tp_theta_ >= 0) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          index_tp_theta_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[layout.idx_theta] + ctx.theta_shift);  // N-body gauge correction
  }
}

void DCDMSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                     const PerturbationsModule& mod,
                                     file_format fmt,
                                     BaseSpecies::TransferColumnSection section) const {
  if (fmt == file_format::class_format) {
    const perturbs* ppt = mod.GetPerturbs();
    if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
      w.Add("d_dcdm", index_tp_delta_, index_tp_delta_ >= 0);
    if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
      w.Add("t_dcdm", index_tp_theta_, index_tp_theta_ >= 0);
  }
}

void DCDMSpecies::PrintVariables(PerturbColumnWriter& w,
                                 const BaseSpecies::PerturbLayout* base,
                                 double /*tau*/,
                                 const double* y,
                                 const PerturbationsModule& mod,
                                 const perturb_workspace* ppw) const {
  double delta_dcdm = 0., theta_dcdm = 0.;

  if (!w.IsTitleMode()) {
    const auto& layout       = static_cast<const PerturbLayout&>(*base);
    const double k           = ppw->scalar_ctx.k;
    const double* pvecback   = ppw->pvecback.data();
    const double* pvecmetric = ppw->pvecmetric.data();

    delta_dcdm = y[layout.idx_delta];
    theta_dcdm = y[layout.idx_theta];

    // Synchronous → Newtonian gauge conversion for DCDM (matter + decay term)
    const perturbs* ppt = mod.GetPerturbs();
    if (ppt->gauge == possible_gauges::synchronous) {
      const double alpha  = pvecmetric[ppw->index_mt_alpha];
      const double H      = pvecback[mod.GetBackgroundModule()->index_bg_H_];
      const double a      = pvecback[mod.GetBackgroundModule()->index_bg_a_];
      delta_dcdm         += alpha * (-a * Gamma_dcdm_ - 3. * a * H);
      theta_dcdm         += k * k * alpha;
    }
  }

  w.Add("delta_dcdm", delta_dcdm, true);
  w.Add("theta_dcdm", theta_dcdm, true);
}
