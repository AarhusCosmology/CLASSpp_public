#include "dcdm_wdm_species.h"

#include <algorithm>
#include <cctype>
#include <cmath>

#include "background_module.h"
#include "errors.h"
#include "perturbations_module.h"

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

DCDM_WDM_Species::DCDM_WDM_Species(std::unique_ptr<WdmDecayProductSpecies> wdm_arg,
                                   const background* pba,
                                   const BackgroundModule* bgm,
                                   double omega0_combined,
                                   double Omega_ini_dcdm)
    : CompositeSpecies(wdm_arg->name(), BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm) {
  auto dcdm =
      std::make_unique<DCDMSpecies>(*pba, omega0_combined, wdm_arg->Gamma(), Omega_ini_dcdm);
  dcdm_ = dcdm.get();
  wdm_  = wdm_arg.get();
  // Child order matters: DCDM first so its rho is in pvecback before the
  // daughter's injection columns are written (ComputeBackground below).
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(wdm_arg));
  scratch_J_.assign(wdm_->q_size(), 0.);
  scratch_dJ_.assign(wdm_->q_size(), 0.);
}

void DCDM_WDM_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

// ─────────────────────────────────────────────────────────────────────────────
// Background wiring
// ─────────────────────────────────────────────────────────────────────────────

void DCDM_WDM_Species::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  CompositeSpecies::ComputeBackground(a, pvecback_B, pvecback);
  // Injection columns need rho_dcdm, so they are written after the child loop.
  wdm_->FillInjection(a, dcdm_->Rho(pvecback), pvecback + wdm_->bg_inj_index(), nullptr);
}

void DCDM_WDM_Species::BackgroundDerivs(double tau,
                                        const double* y,
                                        double* dy,
                                        const double* pvecback) {
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);  // dcdm: rho' (dilution+decay)

  const double a = pvecback[bgm_->index_bg_a_];
  wdm_->FillInjection(a, dcdm_->Rho(pvecback), scratch_J_.data(), scratch_dJ_.data());
  const int N    = wdm_->q_size();
  const int bi_f = wdm_->bi_f_index();
  const int bi_g = wdm_->bi_dfdlnq_index();
  for (int i = 0; i < N; ++i) {
    dy[bi_f + i] = scratch_J_[i];
    dy[bi_g + i] = scratch_dJ_[i];
  }
}

void DCDM_WDM_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm_" + name(), 0.);
  wdm_->WriteBackgroundColumnTitles(w);
}

void DCDM_WDM_Species::WriteBackgroundData(const double* pvecback,
                                           BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm_" + name(), dcdm_->Rho(pvecback));
  wdm_->WriteBackgroundData(pvecback, w);
}

// ─────────────────────────────────────────────────────────────────────────────
// Closure + shooter hooks
// ─────────────────────────────────────────────────────────────────────────────

double DCDM_WDM_Species::GetOmega0() const {
  if (wdm_->Omega_combined_pending().has_value())
    return *wdm_->Omega_combined_pending();
  return CompositeSpecies::GetOmega0();  // pre-shooting discovery fallback
}

std::vector<ShootingTarget> DCDM_WDM_Species::GetShootingTargets() const {
  if (needs_shooting_)
    return {shooting_target_};
  return {};
}

void DCDM_WDM_Species::ComputeShootingGuess(const SpeciesBuildContext& /*ctx*/,
                                            std::vector<double>& guess,
                                            std::vector<double>& dxdy) const {
  const std::string& t = shooting_target_.target_name;
  const double tv      = shooting_target_.target_value;
  if (t == name() + ".Omega_dcdmwdm_fixedpoint") {
    // Initial mode: seed the combined-today reserve with Omega_ini (matter
    // scaling would make them equal; kicks only redshift the kinetic part).
    guess.push_back(tv);
    dxdy.push_back(1.0);
  }
  else {
    // Combined mode: shoot Omega_ini; for slow decays combined-today ~ Omega_ini.
    // For fast decays (Gamma >> H0) the kinetic fraction has redshifted away and
    // combined-today -> eps * Omega_ini, so seed accordingly.
    const double gamma_over_H0 = wdm_->Gamma() / pba_->H0;
    const double eps           = std::max(wdm_->epsilon_retention(), 1e-3);
    const double seed          = (gamma_over_H0 > 1.) ? tv / eps : tv;
    guess.push_back(seed);
    dxdy.push_back(seed / std::max(tv, 1e-30));
  }
}

double DCDM_WDM_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                 const ShootingTarget& target) const {
  const double* bg      = ctx.bg_today;
  const double H0       = ctx.pba->H0;
  const double combined = (dcdm_->Rho(bg) + wdm_->Rho(bg)) / (H0 * H0);
  if (target.target_name == name() + ".Omega_dcdmwdm_fixedpoint") {
    // Initial mode: drive the reserved Omega_dcdmwdm (= GetOmega0()) to the
    // integrated combined density.
    return -combined + GetOmega0();
  }
  // Combined mode: drive the integrated combined to the requested target.
  return combined - target.target_value;
}

// ─────────────────────────────────────────────────────────────────────────────
// Factory
// ─────────────────────────────────────────────────────────────────────────────

std::vector<Named> DCDM_WDM_Species::CreateAll(const SpeciesBuildContext& ctx) {
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  std::vector<Named> result;
  if (instances.empty())
    return result;

  // Synchronous gauge only (the l=1 injection source and the gauge transform of
  // the injected distribution are not implemented; the reference paper works in
  // synchronous gauge). Mirror the Type3 guard: any gauge value starting with
  // "new" (case-insensitive) is Newtonian.
  if (auto gauge = ctx.pfc->get<std::string>("gauge")) {
    std::string g = *gauge;
    std::transform(g.begin(), g.end(), g.begin(), [](unsigned char c) { return std::tolower(c); });
    class_test_severe(g.rfind("new", 0) == 0,
                      "dcdm_wdm supports the synchronous gauge only (arXiv:2606.14849 "
                      "implementation)");
  }

  result.reserve(instances.size());
  for (const auto& nm : instances) {
    (void) ctx.pfc->get<std::string>(nm + ".type");  // mark consumed

    auto wdm =
        std::make_unique<WdmDecayProductSpecies>(ctx.pfc, nm, *ctx.ncdm_settings, ctx.pba, ctx.bgm);

    const bool has_initial  = wdm->Omega_ini_pending().has_value();
    const bool has_combined = wdm->Omega_combined_pending().has_value();
    const double Omega_ini  = wdm->Omega_ini_pending().value_or(0.);
    double omega0_combined  = wdm->Omega_combined_pending().value_or(0.);

    ShootingTarget target{};
    bool needs_shooting    = false;
    const bool in_shooting = ctx.pfc->is_shooting;
    if (has_initial) {
      // Mode (a): fixed-point shoot of the combined-today reserve. During a
      // shooting iteration the shot value is already present (has_combined).
      target         = {nm + ".Omega_dcdmwdm_fixedpoint",
                        nm + ".Omega_dcdmwdm",
                        has_combined ? omega0_combined : Omega_ini};
      needs_shooting = true;
      if (!has_combined)
        omega0_combined = Omega_ini;  // discovery-build seed for the reserve
    }
    else if (has_combined && !in_shooting) {
      // Mode (b): shoot the parent's initial abundance to hit combined-today.
      target         = {nm + ".Omega_dcdmwdm", nm + ".Omega_ini", omega0_combined};
      needs_shooting = true;
    }

    // Mode (b) discovery build needs an Omega_ini seed for the DCDM child.
    double Omega_ini_dcdm = Omega_ini;
    if (!has_initial) {
      std::vector<double> g, d;
      auto tmp = std::make_unique<
          DCDM_WDM_Species>(std::make_unique<WdmDecayProductSpecies>(ctx.pfc,
                                                                     nm,
                                                                     *ctx.ncdm_settings,
                                                                     ctx.pba,
                                                                     ctx.bgm),
                            ctx.pba,
                            ctx.bgm,
                            omega0_combined,
                            /*Omega_ini_dcdm=*/omega0_combined);
      tmp->shooting_target_ = target;
      tmp->needs_shooting_  = true;
      tmp->ComputeShootingGuess(ctx, g, d);
      Omega_ini_dcdm = g[0];
    }

    auto composite = std::make_unique<DCDM_WDM_Species>(std::move(wdm),
                                                        ctx.pba,
                                                        ctx.bgm,
                                                        omega0_combined,
                                                        Omega_ini_dcdm);
    // Seed the reserve for the discovery build in mode (a): GetOmega0() falls
    // back to the child sum (= omega0_combined via the DCDM child) until the
    // shooter pins <inst>.Omega_dcdmwdm.
    composite->shooting_target_ = target;
    composite->needs_shooting_  = needs_shooting;
    result.push_back({nm, std::move(composite)});
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Perturbation coupling: injection source on the daughter's psi_0 bins
// ─────────────────────────────────────────────────────────────────────────────

void DCDM_WDM_Species::AddCouplingDerivs(double /*tau*/,
                                         const double* y,
                                         double* dy,
                                         const perturb_parameters_and_workspace& ppaw) const {
  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector* pv     = ppw->pv.get();

  const auto& dcdm_lay = dcdm_layout(*pv->species_layouts[collection_index_]);
  const auto& wdm_lay  = wdm_layout(*pv->species_layouts[collection_index_]);
  if (dcdm_lay.idx_delta < 0 || wdm_lay.q_size <= 0 || wdm_lay.index_per_q.empty())
    return;

  const double* pvecback  = ppw->pvecback.data();
  const double delta_dcdm = y[dcdm_lay.idx_delta];
  const int inj0          = wdm_->bg_inj_index();
  for (int iq = 0; iq < wdm_lay.q_size; ++iq)
    dy[wdm_lay.index_per_q[iq]] += pvecback[inj0 + iq] * delta_dcdm;
}
