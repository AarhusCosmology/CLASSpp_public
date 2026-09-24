#include "nede_species.h"

#include <cmath>
#include <stdexcept>

#include "background.h"
#include "background_module.h"
#include "bisection.h"
#include "errors.h"
#include "parser.h"
#include "perturbations.h"
#include "perturbations_module.h"
#include "species_collection.h"

NEDESpecies::NEDESpecies(const background& pba,
                         double f_NEDE,
                         double m_NEDE,
                         double w_NEDE,
                         double H_star,
                         double a_star)
    : BaseSpecies("NEDE", EnergyType::DarkEnergy), pba_(pba), m_NEDE_(m_NEDE), w_NEDE_(w_NEDE),
      rho_star_(f_NEDE * H_star * H_star), a_star_(a_star) {}

double NEDESpecies::RhoAt(double a) const {
  return (a < a_star_) ? rho_star_ : rho_star_ * pow(a_star_ / a, 3. * (1. + w_NEDE_));
}

double NEDESpecies::GetOmega0() const {
  return RhoAt(1.) / (pba_.H0 * pba_.H0);
}

std::optional<double> NEDESpecies::GetParam(const std::string& name) const {
  if (name == "z_star")
    return 1. / a_star_ - 1.;
  return std::nullopt;
}

// ── Background ─────────────────────────────────────────────────────────────────

void NEDESpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_, true, index_bg, 1);
  class_define_index(index_bg_phi_prime_, true, index_bg, 1);
}

void NEDESpecies::RegisterIntegrationIndices(int& index_bi) {
  class_define_index(index_bi_phi_, true, index_bi, 1);
  class_define_index(index_bi_phi_prime_, true, index_bi, 1);
}

void NEDESpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  /* Frozen by Hubble friction; the amplitude is arbitrary, see the class comment. */
  ctx.pvecback_integration[index_bi_phi_]       = 1.;
  ctx.pvecback_integration[index_bi_phi_prime_] = 0.;
}

void NEDESpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  pvecback[index_bg_rho_]       = RhoAt(a);
  pvecback[index_bg_phi_prime_] = pvecback_B[index_bi_phi_prime_];
}

void NEDESpecies::BackgroundDerivs(double /*tau*/,
                                   const double* y,
                                   double* dy,
                                   const double* pvecback) {
  const double a = pvecback[bgm_->index_bg_a_];
  /* After the transition nothing reads the trigger, which would otherwise oscillate at
     m = H_* / H_over_m for the rest of the integration. */
  if (a >= a_star_) {
    dy[index_bi_phi_]       = 0.;
    dy[index_bi_phi_prime_] = 0.;
    return;
  }
  /* phi'' + 2 aH phi' + a^2 m^2 phi = 0 */
  dy[index_bi_phi_]       = y[index_bi_phi_prime_];
  dy[index_bi_phi_prime_] = -2. * a * pvecback[bgm_->index_bg_H_] * y[index_bi_phi_prime_] -
                            a * a * m_NEDE_ * m_NEDE_ * y[index_bi_phi_];
}

double NEDESpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_];
}

double NEDESpecies::P(const double* pvecback) const {
  return WAt(pvecback[bgm_->index_bg_a_]) * pvecback[index_bg_rho_];
}

double NEDESpecies::PPrime(double a,
                           double H,
                           const double* /*pvecback_B*/,
                           const double* pvecback) const {
  /* p = w rho with w constant on either side, and rho' = -3 aH (1+w) rho. */
  const double w = WAt(a);
  return -3. * a * H * (1. + w) * w * pvecback[index_bg_rho_];
}

void NEDESpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_NEDE", 0.);
}

void NEDESpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_NEDE", pvecback[index_bg_rho_]);
}

// ── Perturbations ──────────────────────────────────────────────────────────────

int NEDESpecies::ApproximationRegimeAt(int /*which*/, double /*k*/, const double* pvecback) const {
  return (pvecback[bgm_->index_bg_a_] < a_star_) ? kBeforeTransition : kAfterTransition;
}

void NEDESpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                              perturb_vector* /*pv*/,
                                              const precision* /*ppr*/,
                                              int& index_pt,
                                              const perturb_workspace* ppw,
                                              int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  if (ApproximationRegime(ppw) == kBeforeTransition) {
    class_define_index(layout.idx_delta_phi, true, index_pt, 1);
    class_define_index(layout.idx_delta_phi_prime, true, index_pt, 1);
  }
  else {
    class_define_index(layout.idx_delta, true, index_pt, 1);
    class_define_index(layout.idx_theta, true, index_pt, 1);
  }
}

void NEDESpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                double /*tau*/,
                                const double* y,
                                double* dy,
                                const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout              = static_cast<const PerturbLayout&>(base);
  const PerturbScalarContext& ctx = ppaw.ppw->scalar_ctx;

  if (layout.idx_delta_phi >= 0) {
    /* The trigger, sourced by the metric as a scalar field is in synchronous gauge
       (metric_continuity = h'/2). */
    dy[layout.idx_delta_phi] = y[layout.idx_delta_phi_prime];
    dy[layout.idx_delta_phi_prime] =
        -2. * ctx.a_prime_over_a * y[layout.idx_delta_phi_prime] -
        ctx.metric_continuity * ppaw.ppw->pvecback[index_bg_phi_prime_] -
        (ctx.k2 + ctx.a2 * m_NEDE_ * m_NEDE_) * y[layout.idx_delta_phi];
    return;
  }

  /* A fluid with constant w and c_s^2 = w, so the non-adiabatic terms vanish. */
  const double w       = w_NEDE_;
  dy[layout.idx_delta] = -(1. + w) * (y[layout.idx_theta] + ctx.metric_continuity);
  dy[layout.idx_theta] = -(1. - 3. * w) * ctx.a_prime_over_a * y[layout.idx_theta] +
                         w * ctx.k2 / (1. + w) * y[layout.idx_delta] + ctx.metric_euler;
}

void NEDESpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                const BaseSpecies::PerturbLayout& new_base,
                                                const double* old_y,
                                                double* new_y,
                                                const PerturbSwitchContext& ctx) const {
  const auto& old_layout = static_cast<const PerturbLayout&>(old_base);
  const auto& new_layout = static_cast<const PerturbLayout&>(new_base);
  if (old_layout.idx_delta_phi < 0 || new_layout.idx_delta < 0) {
    /* Not NEDE's switch: carry the fluid across. */
    if (old_layout.idx_delta >= 0 && new_layout.idx_delta >= 0) {
      new_y[new_layout.idx_delta] = old_y[old_layout.idx_delta];
      new_y[new_layout.idx_theta] = old_y[old_layout.idx_theta];
    }
    else if (old_layout.idx_delta_phi >= 0 && new_layout.idx_delta_phi >= 0) {
      new_y[new_layout.idx_delta_phi]       = old_y[old_layout.idx_delta_phi];
      new_y[new_layout.idx_delta_phi_prime] = old_y[old_layout.idx_delta_phi_prime];
    }
    return;
  }
  /* Junction conditions of 2006.06686: the transition surface is delayed by
     delta_phi/phi' in conformal time. */
  const double delay          = old_y[old_layout.idx_delta_phi] / ctx.pvecback[index_bg_phi_prime_];
  new_y[new_layout.idx_delta] = -3. * (1. + w_NEDE_) * ctx.a * ctx.pvecback[bgm_->index_bg_H_] *
                                delay;
  new_y[new_layout.idx_theta] = ctx.k * ctx.k * delay;
}

BaseSpecies::StressEnergyContribution NEDESpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = Rho(pvecback);
  se.p   = P(pvecback);
  /* Before the transition NEDE is a cosmological constant and the trigger a clock. */
  if (layout.idx_delta >= 0) {
    se.delta_rho        = se.rho * y[layout.idx_delta];
    se.rho_plus_p_theta = (se.rho + se.p) * y[layout.idx_theta];
    se.delta_p          = w_NEDE_ * se.delta_rho;
  }
  return se;
}

// ── Factory ────────────────────────────────────────────────────────────────────

std::vector<Named> NEDESpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  const auto f_NEDE = ctx.pfc->get<double>("f_NEDE");
  if (!f_NEDE)
    return result;

  /* m_NEDE is in 1/Mpc, as TriggerCLASS's trigger mass: the H0 Olympics prior of
     log10 m_NEDE in [1.3, 3.3] is labelled eV, but 20-2000/Mpc is what puts the
     transition before recombination. */
  const auto m_NEDE = ctx.pfc->get<double>("m_NEDE");
  class_test_severe(!m_NEDE, "NEDE needs 'm_NEDE', the trigger mass in 1/Mpc");
  const double w_NEDE   = ctx.pfc->get_or("w_NEDE", 2. / 3.);
  const double H_over_m = ctx.pfc->get_or("H_over_m_NEDE", 0.2);
  class_test(*f_NEDE <= 0. || *f_NEDE >= 1., "f_NEDE = %g must lie in (0, 1)", *f_NEDE);
  class_test(*m_NEDE <= 0., "m_NEDE = %g must be positive", *m_NEDE);
  class_test(w_NEDE < 0. || w_NEDE > 1., "w_NEDE = %g must lie in [0, 1]", w_NEDE);
  class_test(H_over_m <= 0., "H_over_m_NEDE = %g must be positive", H_over_m);
  if (auto gauge = ctx.pfc->get<std::string>("gauge"))
    class_test_severe(gauge->find("new") != std::string::npos ||
                          gauge->find("New") != std::string::npos,
                      "NEDE is implemented in synchronous gauge only");

  /* The transition is where H = H_*, i.e. where everything but NEDE contributes
     H_*^2 (1 - f_NEDE) to H^2 = rho_tot - K/a^2. The species are all built by now except
     the closure species, which is dark energy and negligible there. */
  if (ctx.all_species == nullptr)
    throw std::logic_error("NEDE needs the species collection (internal: all_species missing)");
  const double H0     = ctx.pba->H0;
  const double H_star = H_over_m * *m_NEDE;
  const double target = H_star * H_star * (1. - *f_NEDE) / (H0 * H0);
  auto H2_others      = [&](double a) {  // in units of H0^2
    double sum = ctx.pba->Omega0_k / (a * a);
    for (const auto& sp : *ctx.all_species)
      sum += sp->BackgroundDensityOverH0Sq(a, H0);
    return sum;
  };
  const double a_ini = ctx.ppr->a_ini_over_a_today_default;
  class_test(H2_others(1.) >= target || H2_others(a_ini) <= target,
             "NEDE: H never equals H_* = %g/Mpc between a = %g and today",
             H_star,
             a_ini);
  const double ln_a_star = bisect_value(log(a_ini), 0., 1e-12, [&](double ln_a) {
    return H2_others(exp(ln_a)) < target;
  });

  result.push_back(
      {"NEDE",
       std::make_unique<NEDESpecies>(*ctx.pba, *f_NEDE, *m_NEDE, w_NEDE, H_star, exp(ln_a_star))});
  return result;
}
