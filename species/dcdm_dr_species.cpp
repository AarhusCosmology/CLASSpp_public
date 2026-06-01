#include "dcdm_dr_species.h"

#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "perturbations_module.h"

DCDM_DR_Species::DCDM_DR_Species(const background* pba,
                                 const BackgroundModule* bgm,
                                 double omega0_dcdmdr,
                                 double Gamma_dcdm,
                                 double Omega_ini_dcdm)
    : CompositeSpecies("DCDM_DR", BaseSpecies::EnergyType::Other), pba_(pba), bgm_(bgm) {
  auto dcdm  = std::make_unique<DCDMSpecies>(*pba, omega0_dcdmdr, Gamma_dcdm, Omega_ini_dcdm);
  auto dr_sp = std::make_unique<DarkRadiationSpecies>("DR", pba, bgm);
  dcdm_      = dcdm.get();
  dr_sp_     = dr_sp.get();
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(dr_sp));
}

void DCDM_DR_Species::SetBackgroundModule(const BackgroundModule* bgm) {
  bgm_ = bgm;
  CompositeSpecies::SetBackgroundModule(bgm);
}

void DCDM_DR_Species::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm", 0.);
  w.Add("(.)rho_dr", 0.);
}

void DCDM_DR_Species::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_dcdm", dcdm().Rho(pvecback));
  w.Add("(.)rho_dr", dr().Rho(pvecback));
}

void DCDM_DR_Species::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  // Initialize children first (DCDM)
  CompositeSpecies::SetBackgroundInitialConditions(ctx);

  // Then add the DR initial condition from DCDM decay.  This method only runs
  // when a DCDM_DR_Species exists, so the legacy pba->has_dcdm guard is gone.
  // Total radiation Omega0 at a_ini, recovered from the module-wide rho_rad
  // (rho_rad = Omega_rad * H0^2 / a_rel^4). This includes every radiation-like
  // species (photons, UR, IDR, ...) consistently and needs no cross-species lookup.
  const double Omega_rad    = ctx.rho_rad * std::pow(ctx.a_rel, 4) / std::pow(pba_->H0, 2);
  const double rho_dcdm_ini = ctx.pvecback_integration[dcdm_->bi_rho_index()];
  double f = 1. / 3. * std::pow(ctx.a_rel, 6) * rho_dcdm_ini * dcdm_->Gamma_dcdm() /
             std::pow(pba_->H0, 3) / std::sqrt(Omega_rad);
  ctx.pvecback_integration[dr_sp_->bi_rho_index()] = f * std::pow(pba_->H0, 2) /
                                                     std::pow(ctx.a_rel, 4);
}

void DCDM_DR_Species::BackgroundDerivs(double tau,
                                       const double* y,
                                       double* dy,
                                       const double* pvecback) {
  // Children handle their own dilution terms
  CompositeSpecies::BackgroundDerivs(tau, y, dy, pvecback);

  // DCDM->DR decay source
  const double a              = pvecback[bgm_->index_bg_a_];
  dy[dr_sp_->bi_rho_index()] += a * dcdm_->Gamma_dcdm() * dcdm_->Rho(pvecback);
}

// ── Perturbation layout-based overrides ───────────────────────────────────────

void DCDM_DR_Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                  perturb_vector* pv,
                                                  const precision* ppr,
                                                  int& index_pt,
                                                  const perturb_workspace* ppw,
                                                  int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  dcdm_->RegisterPerturbationIndices(my.dcdm, pv, ppr, index_pt, ppw, gauge);
  dr_sp_->RegisterPerturbationIndices(my.dr, pv, ppr, index_pt, ppw, gauge);
}

void DCDM_DR_Species::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                    double tau,
                                    const double* y,
                                    double* dy,
                                    const perturb_parameters_and_workspace& ppaw) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dcdm_->PerturbDerivs(my.dcdm, tau, y, dy, ppaw);
  dr_sp_->PerturbDerivs(my.dr, tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);
}

void DCDM_DR_Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                             double* y,
                                             const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dcdm_->ApplyInitialConditions(my.dcdm, y, ctx);
  dr_sp_->ApplyInitialConditions(my.dr, y, ctx);
}

void DCDM_DR_Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                    double* y,
                                                    const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dcdm_->PerturbSynchronousToNewtonian(my.dcdm, y, ctx);  // matter+decay via RhoDotOverRho
  const double* pvecback  = ctx.ppw->pvecback;
  const double rho_dr     = dr_sp_->Rho(pvecback);
  const double decay_corr = (rho_dr > 0.)
                                ? ctx.a * dcdm_->Gamma_dcdm() * dcdm_->Rho(pvecback) / rho_dr
                                : 0.;
  dr_sp_->PerturbNewtonianReseed(my.dr, y, ctx, decay_corr);
}

double DCDM_DR_Species::Delta(const BaseSpecies::PerturbLayout& base,
                              const perturb_vector* pv,
                              const double* y,
                              const double* pvecback,
                              const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rho_d   = dcdm_->Rho(pvecback);
  const double rho_dr  = dr_sp_->Rho(pvecback);
  const double rho_tot = rho_d + rho_dr;
  if (rho_tot <= 0.)
    return 0.;
  return (rho_d * dcdm_->Delta(my.dcdm, pv, y, pvecback, ppw) +
          rho_dr * dr_sp_->Delta(my.dr, pv, y, pvecback, ppw)) /
         rho_tot;
}

double DCDM_DR_Species::Theta(const BaseSpecies::PerturbLayout& base,
                              const perturb_vector* pv,
                              const double* y,
                              const double* pvecback,
                              const perturb_workspace* ppw) const {
  const auto& my       = static_cast<const PerturbLayout&>(base);
  const double rpp_d   = dcdm_->Rho(pvecback) + dcdm_->P(pvecback);
  const double rpp_dr  = dr_sp_->Rho(pvecback) + dr_sp_->P(pvecback);
  const double rpp_tot = rpp_d + rpp_dr;
  if (rpp_tot <= 0.)
    return 0.;
  return (rpp_d * dcdm_->Theta(my.dcdm, pv, y, pvecback, ppw) +
          rpp_dr * dr_sp_->Theta(my.dr, pv, y, pvecback, ppw)) /
         rpp_tot;
}

double DCDM_DR_Species::MatterRhoDelta(const perturb_vector* pv,
                                       const double* y,
                                       const double* pvecback,
                                       const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return dcdm_->IsMatterSpecies()
             ? dcdm_->Rho(pvecback) * dcdm_->Delta(my.dcdm, pv, y, pvecback, ppw)
             : 0.;
}

double DCDM_DR_Species::MatterRhoPlusPTheta(const perturb_vector* pv,
                                            const double* y,
                                            const double* pvecback,
                                            const perturb_workspace* ppw) const {
  if (collection_index_ >= pv->species_layouts.size())
    return 0.;
  const auto& my = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);
  return dcdm_->IsMatterSpecies() ? (dcdm_->Rho(pvecback) + dcdm_->P(pvecback)) *
                                        dcdm_->Theta(my.dcdm, pv, y, pvecback, ppw)
                                  : 0.;
}

double DCDM_DR_Species::DeltaP(const BaseSpecies::PerturbLayout& base,
                               const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  return dcdm_->DeltaP(my.dcdm, pv, y, pvecback, ppw) + dr_sp_->DeltaP(my.dr, pv, y, pvecback, ppw);
}

double DCDM_DR_Species::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                      const perturb_vector* pv,
                                      const double* y,
                                      const double* pvecback,
                                      const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  return dcdm_->RhoPlusPShear(my.dcdm, pv, y, pvecback, ppw) +
         dr_sp_->RhoPlusPShear(my.dr, pv, y, pvecback, ppw);
}

void DCDM_DR_Species::FillSources(const double* y,
                                  const double* /*dy*/,
                                  PerturbSourceContext& ctx) {
  PerturbationsModule* p_mod = ctx.p_mod;
  perturb_workspace* ppw     = ctx.ppw;
  const perturb_vector* pv   = ppw->pv;
  const double* pvecback     = ppw->pvecback;

  const double a_prime_over_a = ctx.a_prime_over_a;
  const double a2_rel         = ctx.a2_rel;

  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  const auto& my_lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);

  // ── delta_dcdm ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_delta_dcdm_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_dcdm_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[my_lay.dcdm.idx_delta] +
                              (3. * a_prime_over_a + ctx.a_rel * dcdm_->Gamma_dcdm()) *
                                  ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_dcdm ─────────────────────────────────────────────────────────────
  if (p_mod->has_source_theta_dcdm_ == _TRUE_) {
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_dcdm_,
                          ctx.index_tau,
                          ctx.index_k,
                          y[my_lay.dcdm.idx_theta] + ctx.theta_shift);  // N-body gauge correction
  }

  // ── delta_dr (this channel's slot) ───────────────────────────────────────────
  // r_dr == 0 when the channel carries no DR (e.g. Gamma_dcdm == 0); write 0 then,
  // matching DarkRadiationSpecies::Delta's rho_dr <= 0 convention (avoids 0/0).
  if (p_mod->has_source_delta_dr_ == _TRUE_) {
    const double r_dr = (a2_rel / pba_->H0) * (a2_rel / pba_->H0) * dr_sp_->Rho(pvecback);
    const double src  = (r_dr > 0.)
                            ? y[my_lay.dr.idx_F0] / r_dr +
                                  4. * a_prime_over_a * ctx.theta_over_k2  // N-body gauge corr.
                            : 0.;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_delta_dr_ + dr_sp_->source_slot(),
                          ctx.index_tau,
                          ctx.index_k,
                          src);
  }

  // ── theta_dr (this channel's slot) ───────────────────────────────────────────
  if (p_mod->has_source_theta_dr_ == _TRUE_) {
    const double r_dr = (a2_rel / pba_->H0) * (a2_rel / pba_->H0) * dr_sp_->Rho(pvecback);
    const double src  = (r_dr > 0.) ? 3. / 4. * ctx.k * y[my_lay.dr.idx_F0 + 1] / r_dr +
                                          ctx.theta_shift  // N-body gauge correction
                                    : 0.;
    p_mod->SetSourceValue(ctx.index_md,
                          ctx.index_ic,
                          p_mod->index_tp_theta_dr_ + dr_sp_->source_slot(),
                          ctx.index_tau,
                          ctx.index_k,
                          src);
  }
}

void DCDM_DR_Species::WriteOutputColumns(PerturbColumnWriter& w,
                                         const PerturbationsModule& mod,
                                         enum file_format fmt,
                                         BaseSpecies::TransferColumnSection section) const {
  if (fmt != class_format)
    return;
  const perturbs* ppt = mod.GetPerturbs();
  const int slot      = dr_sp_->source_slot();
  if (section != TransferColumnSection::velocity && ppt->has_density_transfers == _TRUE_)
    w.Add("d_" + dr_sp_->name(), mod.index_tp_delta_dr_ + slot, mod.has_source_delta_dr_ == _TRUE_);
  if (section != TransferColumnSection::density && ppt->has_velocity_transfers == _TRUE_)
    w.Add("t_" + dr_sp_->name(), mod.index_tp_theta_dr_ + slot, mod.has_source_theta_dr_ == _TRUE_);
}

void DCDM_DR_Species::AddCouplingDerivs(double /*tau*/,
                                        const double* y,
                                        double* dy,
                                        const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  const auto& my_cd_lay = static_cast<const PerturbLayout&>(
      *pv->species_layouts[collection_index_]);
  if (my_cd_lay.dcdm.idx_delta < 0 || my_cd_lay.dr.idx_F0 < 0)
    return;

  const double* pvecback = ppw->pvecback;
  const double a         = ctx.a;
  const double k         = ctx.k;

  const int base         = my_cd_lay.dr.idx_F0;
  const double rprime_dr = dcdm_->Gamma_dcdm() * dcdm_->Rho(pvecback) * std::pow(a, 5) /
                           (pba_->H0 * pba_->H0);

  const double delta_dcdm = y[my_cd_lay.dcdm.idx_delta];
  const double theta_dcdm = y[my_cd_lay.dcdm.idx_theta];

  // Add DCDM source to this channel's DR l=0 and l=1
  dy[base + 0] += rprime_dr * (delta_dcdm + ctx.metric_euler / (k * k));
  dy[base + 1] += rprime_dr / k * theta_dcdm;
}

void DCDM_DR_Species::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                    const BaseSpecies::PerturbLayout& new_base,
                                                    const double* old_y,
                                                    double* new_y,
                                                    const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  dcdm_->CopyPerturbationsAcrossSwitch(old_l.dcdm, new_l.dcdm, old_y, new_y, ctx);
  dr_sp_->CopyPerturbationsAcrossSwitch(old_l.dr, new_l.dr, old_y, new_y, ctx);
}

// ── Shooter hooks ─────────────────────────────────────────────────────────────

std::vector<ShootingTarget> DCDM_DR_Species::GetShootingTargets() const {
  if (needs_shooting_)
    return {shooting_target_};
  return {};
}

void DCDM_DR_Species::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                           std::vector<double>& guess,
                                           std::vector<double>& dxdy) const {
  const background& ba = *ctx.pba;
  const std::string& t = shooting_target_.target_name;
  const double tv      = shooting_target_.target_value;
  // For Omega_M we need CDM / IDM_DR / DCDM_DR contributions, sourced from the
  // SpeciesOmegaBudget that InputModule fills before species construction.
  const SpeciesOmegaBudget* budget = ctx.omega_budget;
  const double Omega0_cdm_eff      = budget ? budget->cdm.value_or(0.) : 0.;
  const double Omega0_idm_dr_eff   = budget ? budget->idm_dr.value_or(0.) : 0.;
  const double Omega0_dcdmdr_eff   = budget ? budget->dcdmdr.value_or(0.) : 0.;
  if (t == "Omega_dcdmdr" || t == "omega_dcdmdr") {
    double Omega_M = Omega0_cdm_eff + Omega0_idm_dr_eff + Omega0_dcdmdr_eff + ba.Omega0_b;
    double gamma = dcdm_->Gamma_dcdm() / ba.H0, a_decay = 1.0;
    if (gamma > 1)
      a_decay = pow(1 + (gamma * gamma - 1.) / Omega_M, -1. / 3.);
    double tgt = (t == "omega_dcdmdr") ? tv / ba.h / ba.h : tv;
    guess.push_back(tgt / a_decay);
    dxdy.push_back((tgt / a_decay) / tv);
  }
  else {  // Omega_ini_dcdm / omega_ini_dcdm
    double Omega0_dcdmdr_ini = (t == "omega_ini_dcdm") ? tv / (ba.h * ba.h) : tv;
    double Omega_M           = Omega0_cdm_eff + Omega0_idm_dr_eff + Omega0_dcdmdr_ini + ba.Omega0_b;
    double gamma = dcdm_->Gamma_dcdm() / ba.H0, a_decay = 1.0;
    if (gamma > 1)
      a_decay = pow(1 + (gamma * gamma - 1.) / Omega_M, -1. / 3.);
    guess.push_back(tv * a_decay);
    double d = a_decay;
    if (gamma > 100)
      d *= gamma / 100;
    dxdy.push_back(d);
  }
}

double DCDM_DR_Species::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                                const ShootingTarget& target) const {
  const background& ba  = *ctx.pba;
  const double* bg      = ctx.bg_today;
  const double H0       = ba.H0;
  const double rho_dcdm = dcdm_->Rho(bg);
  // DCDM_DR composite always owns a DR sub-species, so the legacy has_dr guard is gone.
  const double rho_dr  = dr_sp_->Rho(bg);
  const std::string& t = target.target_name;
  if (t == "Omega_dcdmdr") {
    return (rho_dcdm + rho_dr) / (H0 * H0) - target.target_value;
  }
  else if (t == "omega_dcdmdr") {
    return (rho_dcdm + rho_dr) / (H0 * H0) - target.target_value / ba.h / ba.h;
  }
  // Omega_ini_dcdm / omega_ini_dcdm: target is the ini, we vary the today density.
  return -(rho_dcdm + rho_dr) / (H0 * H0) + dcdm_->GetOmega0();
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> DCDM_DR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // Existence detection: dcdm_dr exists iff the budget says so (parser block in
  // ReadCoupledOmegaBudget filled .dcdmdr from Omega_dcdmdr/omega_dcdmdr) OR a
  // shooter target (Omega_ini_dcdm/omega_ini_dcdm) is present in pfc.
  const double omega0_dcdmdr_from_budget = (ctx.omega_budget ? ctx.omega_budget->dcdmdr.value_or(0.)
                                                             : 0.);
  double omega0_dcdmdr                   = omega0_dcdmdr_from_budget;
  {
    double tmp;
    bool exists = (omega0_dcdmdr != 0.) || ctx.pfc->read_double("Omega_ini_dcdm", tmp) ||
                  ctx.pfc->read_double("omega_ini_dcdm", tmp);
    if (!exists)
      return result;
  }

  // ── Parse DCDM physics params (Gamma_dcdm + Omega_ini_dcdm) ────────────────
  // Gamma_dcdm is read in km/(s·Mpc) and converted to Mpc^-1 (same as the
  // original input_module code path).  Default 0 keeps the legacy semantics.
  double Gamma_dcdm = 0.;
  if (ctx.pfc->read_double("Gamma_dcdm", Gamma_dcdm))
    Gamma_dcdm *= 1.e3 / _c_;

  // Omega_ini_dcdm: direct (Omega_ini_dcdm) XOR rescaled (omega_ini_dcdm = Omega*h^2).
  // The two are mutually exclusive (replicates ReadCoupledOmegaBudget's class_test).
  double Omega_ini_dcdm     = 0.;
  double omega_ini_dcdm_val = 0.;
  bool has_Omega_ini        = ctx.pfc->read_double("Omega_ini_dcdm", Omega_ini_dcdm);
  bool has_omega_ini        = ctx.pfc->read_double("omega_ini_dcdm", omega_ini_dcdm_val);
  if (has_Omega_ini && has_omega_ini)
    throw std::invalid_argument(
        "In input file, you can only enter one of Omega_ini_dcdm or omega_ini_dcdm, choose one");
  if (has_omega_ini)
    Omega_ini_dcdm = omega_ini_dcdm_val / (ctx.pba->h * ctx.pba->h);

  // Detect guess-driven construction: target key present but direct unknown absent.
  // Mapping: Omega_dcdmdr/omega_dcdmdr -> unknown Omega_ini_dcdm
  //          Omega_ini_dcdm/omega_ini_dcdm -> unknown Omega_dcdmdr
  struct TargetSpec {
    const char* target_name;
    const char* unknown_param;
  };
  static const TargetSpec kTargets[] = {
      {"Omega_dcdmdr", "Omega_ini_dcdm"},
      {"omega_dcdmdr", "Omega_ini_dcdm"},
      {"Omega_ini_dcdm", "Omega_dcdmdr"},
      {"omega_ini_dcdm", "omega_dcdmdr"},
  };

  ShootingTarget shooting_target{};
  bool needs_shooting = false;

  for (const auto& ts : kTargets) {
    double target_val    = 0.;
    double unknown_val   = 0.;
    bool target_present  = ctx.pfc->read_double(ts.target_name, target_val);
    bool unknown_present = ctx.pfc->read_double(ts.unknown_param, unknown_val);
    if (target_present && !unknown_present) {
      shooting_target = {ts.target_name, ts.unknown_param, target_val};
      needs_shooting  = true;

      // Set the target's own value for use by the discovery build's
      // ComputeShootingGuess (ini-form targets were not parsed above when
      // Omega0_dcdmdr was 0 in the budget — pfc reads above only set
      // Omega_ini_dcdm when the user gave it directly).
      const std::string tname(ts.target_name);
      if (tname == "Omega_ini_dcdm") {
        Omega_ini_dcdm = target_val;
      }
      else if (tname == "omega_ini_dcdm") {
        Omega_ini_dcdm = target_val / (ctx.pba->h * ctx.pba->h);
      }
      // Omega_dcdmdr/omega_dcdmdr: already accounted for in omega0_dcdmdr above.

      // Build a temporary composite to compute the shooting guess; discard it
      // once the guess has been applied. The final composite is built below
      // with the post-mutation omega0_dcdmdr / Omega_ini_dcdm.
      auto tmp_composite              = std::make_unique<DCDM_DR_Species>(ctx.pba,
                                                                          ctx.bgm,
                                                                          omega0_dcdmdr,
                                                                          Gamma_dcdm,
                                                                          Omega_ini_dcdm);
      tmp_composite->shooting_target_ = shooting_target;
      tmp_composite->needs_shooting_  = true;

      std::vector<double> g, d;
      tmp_composite->ComputeShootingGuess(ctx, g, d);

      // Seed the guess into the UNKNOWN so this (discovery) build is valid.
      // DoShooting later writes the *resolved* unknown back into the fc; the
      // final build then re-parses it from there.
      const std::string uname(ts.unknown_param);
      if (uname == "Omega_ini_dcdm") {
        Omega_ini_dcdm = g[0];
      }
      else if (uname == "Omega_dcdmdr") {
        omega0_dcdmdr = g[0];
      }
      else if (uname == "omega_dcdmdr") {
        omega0_dcdmdr = g[0] / (ctx.pba->h * ctx.pba->h);
      }
      break;
    }
  }

  // Construct the final composite with the post-mutation omega0_dcdmdr and
  // Omega_ini_dcdm so the DCDM child stores the correct existence-time values.
  auto composite              = std::make_unique<DCDM_DR_Species>(ctx.pba,
                                                                  ctx.bgm,
                                                                  omega0_dcdmdr,
                                                                  Gamma_dcdm,
                                                                  Omega_ini_dcdm);
  composite->shooting_target_ = shooting_target;
  composite->needs_shooting_  = needs_shooting;

  result.push_back({"DCDM_DR", std::move(composite)});
  return result;
}
