#include "dncdm_species.h"

#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

#include "arrays.h"
#include "background_module.h"
#include "errors.h"
#include "perturbations_module.h"
#include "species/species_input.h"

namespace {

constexpr const char* kLegacyDecayDrKeys[] = {
    "N_ncdm_decay_dr",
    "m_ncdm_decay_dr",
    "T_ncdm_decay_dr",
    "ksi_ncdm_decay_dr",
    "deg_ncdm_decay_dr",
    "quadrature_strategy_ncdm_decay_dr",
    "N_momentum_bins_ncdm_decay_dr",
    "maximum_q_ncdm_decay_dr",
    "Gamma_ncdm_decay_dr",
    "log10Gamma_ncdm_decay_dr",
    "lifetime_ncdm_decay_dr",
    "log10lifetime_ncdm_decay_dr",
    "Omega_dncdmdr",
    "omega_dncdmdr",
    "Omega_ini_dncdm",
    "omega_ini_dncdm",
    "Neff_ini_dncdm",
};

void RejectLegacyDecayDrKeys(FileContent& pfc) {
  for (const char* key : kLegacyDecayDrKeys) {
    class_test_severe(pfc.get<std::string>(key).has_value(),
                      "'%s' is no longer supported. Use dot-syntax: '<instance>.<dot-name> = "
                      "...' with '<instance>.type = ncdm_decay_dr'.",
                      key);
  }
}

void ApplyDncdmInitialClosure(DNCDMSpecies& sp, const SpeciesBuildContext& ctx) {
  const auto& Omega_ini = sp.Omega_ini_pending();
  const auto& Neff_ini  = sp.Neff_ini_pending();
  if (!Omega_ini.has_value() && !Neff_ini.has_value()) {
    return;
  }

  const double a_ini_seed = 1e-14;
  const double a_ini      = sp.GetIni(a_ini_seed, ctx.ncdm_settings->tol_ncdm);
  const double z_ini      = 1.0 / a_ini - 1.0;
  const double H0         = ctx.pba ? ctx.pba->H0 : ctx.ncdm_settings->h * 1.e5 / _c_;

  if (Omega_ini.has_value()) {
    sp.SetDeg_from_Omega_ini(z_ini, H0, *Omega_ini);
  }
  else {
    // Neff_ini → Omega_ini using photon energy density at the same epoch.
    const double Omega0_g        = ctx.pba ? ctx.pba->Omega0_g : 0.;
    const double Omega_ini_value = *Neff_ini * 7. / 8. * std::pow(4. / 11., 4. / 3.) * Omega0_g;
    sp.SetDeg_from_Omega_ini(z_ini, H0, Omega_ini_value);
  }
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

DNCDMSpecies::DNCDMSpecies(FileContent* pfc,
                           const std::string& instance_name,
                           const NcdmSettings& settings,
                           const background* pba,
                           const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings), pba_(pba) {
  bgm_ = bgm;

  SpeciesInput input(pfc, instance_name);
  const int strategy = input.get_or("quadrature_strategy", static_cast<int>(qm_auto));
  if (strategy != qm_auto) {
    if (auto bg_bins = input.get<int>("momenta_bins_bg")) {
      class_test_severe(*bg_bins != q_size(),
                        "species '%s': momenta_bins_bg must match momenta_bins for "
                        "ncdm_decay_dr; its background distribution is evolved on the "
                        "perturbation momentum grid",
                        instance_name.c_str());
    }
  }

  // Read decay rate (exactly one of these must be specified)
  auto Gamma_value         = input.get<double>("Gamma");
  auto log10Gamma_value    = input.get<double>("log10Gamma");
  auto lifetime_value      = input.get<double>("lifetime");
  auto log10lifetime_value = input.get<double>("log10lifetime");

  int n_provided = (int) Gamma_value.has_value() + (int) log10Gamma_value.has_value() +
                   (int) lifetime_value.has_value() + (int) log10lifetime_value.has_value();
  class_test_severe(n_provided != 1,
                    "species '%s': specify exactly one of Gamma, log10Gamma, lifetime, "
                    "log10lifetime",
                    instance_name.c_str());

  double Gamma_raw = 0.;
  if (Gamma_value) {
    Gamma_raw = *Gamma_value;
  }
  else if (log10Gamma_value) {
    Gamma_raw = std::pow(10., *log10Gamma_value);
  }
  else if (lifetime_value) {
    Gamma_raw = 1. / *lifetime_value / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  else {  // log10lifetime
    double lifetime = std::pow(10., *log10lifetime_value);
    Gamma_raw       = 1. / lifetime / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  Gamma_ = Gamma_raw * (1.e3 / _c_);

  // Read DNCDM-specific deg / Omega_ini values from the same instance.
  // These were previously bulk-loaded in CreateAll from CSV lists.
  // Semantics: at most one of {deg, Omega_ini, omega_ini, Neff_ini, Omega_dncdmdr}
  // (resolved later by ConstructSpecies / DoShooting). Bare today-matter Omega/omega is
  // rejected: it leaves the decay radiation out of the flatness budget.
  class_test_severe(
      input.get<double>("Omega").has_value() || input.get<double>("omega").has_value(),
      "species '%s': a decaying species (ncdm_decay_dr) cannot be normalized by today-matter "
      "'Omega'/'omega' (the decay radiation would be left out of the flatness budget). Use "
      "'Omega_ini'/'omega_ini'/'Neff_ini' to pin the initial abundance, or "
      "'Omega_dncdmdr'/'omega_dncdmdr' to pin the combined matter+radiation density today.",
      instance_name.c_str());

  auto deg_opt           = input.get<double>("deg");
  auto Omega_ini_opt     = input.get<double>("Omega_ini");
  auto omega_ini_opt     = input.get<double>("omega_ini");
  auto Neff_ini_opt      = input.get<double>("Neff_ini");
  auto Omega_dncdmdr_opt = input.get<double>("Omega_dncdmdr");
  auto omega_dncdmdr_opt = input.get<double>("omega_dncdmdr");

  double deg_local           = deg_opt.value_or(0.);
  double Omega_ini_local     = Omega_ini_opt.value_or(0.);
  double omega_ini_local     = omega_ini_opt.value_or(0.);
  double Neff_ini_local      = Neff_ini_opt.value_or(0.);
  double Omega_dncdmdr_local = Omega_dncdmdr_opt.value_or(0.);
  double omega_dncdmdr_local = omega_dncdmdr_opt.value_or(0.);
  bool has_deg               = deg_opt.has_value();
  bool has_Omega_ini         = Omega_ini_opt.has_value();
  bool has_omega_ini         = omega_ini_opt.has_value();
  bool has_Neff_ini          = Neff_ini_opt.has_value();
  bool has_Omega_dncdmdr     = Omega_dncdmdr_opt.has_value();
  bool has_omega_dncdmdr     = omega_dncdmdr_opt.has_value();
  class_test_severe(has_Omega_dncdmdr && has_omega_dncdmdr,
                    "species '%s': specify exactly one of Omega_dncdmdr, omega_dncdmdr",
                    instance_name.c_str());
  if (has_omega_dncdmdr) {
    Omega_dncdmdr_local = omega_dncdmdr_local / settings.h / settings.h;
    has_Omega_dncdmdr   = true;
  }

  // Omega_dncdmdr is mutually exclusive with Omega_ini/omega_ini/Neff_ini — EXCEPT inside a
  // shooting build, where DoShooting legitimately writes <flavor>.Omega_dncdmdr as the
  // initial-mode fixed-point unknown (or .deg for combined mode) alongside the user's key.
  const bool in_shooting = (pfc != nullptr) && pfc->is_shooting;
  class_test_severe(has_Omega_dncdmdr && (has_Omega_ini || has_omega_ini || has_Neff_ini) &&
                        !in_shooting,
                    "species '%s': Omega_dncdmdr conflicts with Omega_ini/omega_ini/Neff_ini",
                    instance_name.c_str());

  // Count mutually-exclusive target specs (skip inside a shooting build: the shooter adds
  // Omega_dncdmdr alongside the user's initial-abundance key).
  if (!in_shooting) {
    int n_deg_options = (int) has_deg + (int) has_Omega_ini + (int) has_omega_ini +
                        (int) has_Neff_ini + (int) has_Omega_dncdmdr;
    class_test_severe(n_deg_options > 1,
                      "species '%s': specify at most one of deg, Omega_ini, omega_ini, Neff_ini, "
                      "Omega_dncdmdr",
                      instance_name.c_str());
  }

  // Stash Omega_dncdmdr_pending_ whenever Omega_dncdmdr is present: combined mode (user input),
  // or the initial-mode fixed-point unknown the shooter wrote. GetOmega0() on the DNCDM_DR
  // composite reads this to reserve the combined sector density.
  if (has_Omega_dncdmdr)
    Omega_dncdmdr_pending_ = Omega_dncdmdr_local;

  // Stash the initial-abundance target; deg is set later by DNCDMSpecies::CreateAll's
  // ApplyDncdmInitialClosure (a_ini-driven). In a shooting iteration this co-occurs with the
  // shot Omega_dncdmdr above — the initial key drives deg, Omega_dncdmdr only the closure reserve.
  const bool has_initial = has_Omega_ini || has_omega_ini || has_Neff_ini;
  if (has_Omega_ini)
    Omega_ini_pending_ = Omega_ini_local;  // already in Omega units
  else if (has_omega_ini)
    Omega_ini_pending_ = omega_ini_local / settings.h / settings.h;
  else if (has_Neff_ini)
    Neff_ini_pending_ = Neff_ini_local;

  if (has_deg) {
    // deg wins: used directly (also covers the combined-mode shooting iteration where
    // DoShooting wrote .deg).
    SetDegAndFactor(deg_local);
  }
  else if (has_Omega_dncdmdr && !has_initial) {
    // Combined mode: set a guessed deg via DegGuessFromOmegaToday.
    // DegGuessFromOmegaToday calls GetDeg() — deg_ is initialized to 1. by NCDMBaseSpecies,
    // so the per-unit-degeneracy calculation is well-defined without a prior SetDegAndFactor.
    if (pba) {
      // Build context needed for GetIni; assemble a minimal one from available data.
      SpeciesBuildContext guess_ctx{nullptr, pba, nullptr, &settings, nullptr};
      auto [g, d] = DegGuessFromOmegaToday(guess_ctx, Omega_dncdmdr_local);
      SetDegAndFactor(g);
    }
    // If pba is null (early-construction pass), leave deg_ at 1. and let the shooting
    // dispatch provide the correct value when DoShooting writes <instance>.deg.
  }
  // else: initial mode — deg comes from ApplyDncdmInitialClosure (stashed above).

  // Compute dq_[i] = w_bg_[i] / f0(q_bg_[i]) (matches existing semantics)
  dq_ = ComputeDq();
}

// ─────────────────────────────────────────────────────────────────────────────
// CreateAll factory
// ─────────────────────────────────────────────────────────────────────────────

std::vector<DNCDMSpecies::Named> DNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyDecayDrKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    (void) ctx.pfc->get<std::string>(name + ".type");  // mark consumed

    auto sp = std::make_unique<DNCDMSpecies>(ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    ApplyDncdmInitialClosure(*sp, ctx);
    result.push_back({name, std::move(sp)});
  }
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// DegGuessFromOmegaToday — initial-guess helper for Omega_dncdmdr shooting
// ─────────────────────────────────────────────────────────────────────────────

std::pair<double, double> DNCDMSpecies::DegGuessFromOmegaToday(const SpeciesBuildContext& ctx,
                                                               double Omega_target) const {
  const background& ba = *ctx.pba;
  // a_ini: earliest scale factor at which the NCDM distribution is non-relativistic enough
  // to integrate numerically.  Mirror input_module.cpp:3764-3769.
  double a_ini         = ctx.ppr ? ctx.ppr->a_ini_over_a_today_default : 1e-14;
  const double tol_ini = ctx.ppr ? ctx.ppr->tol_ncdm_initial_w : ctx.ncdm_settings->tol_ncdm;
  a_ini                = GetIni(a_ini, tol_ini);
  const double z_ini   = 1.0 / a_ini - 1.0;

  double rho_actual;
  ComputeMomenta(z_ini, nullptr, &rho_actual, nullptr, nullptr, nullptr);

  // rho per unit degeneracy at a_ini, converted to Omega units scaled to a_ini.
  const double rho_deg1   = (GetDeg() > 0.) ? rho_actual / GetDeg() : 0.;
  const double Omega_deg1 = rho_deg1 * std::pow(a_ini, 4.0) / (ba.H0 * ba.H0);

  double Omega_ini = Omega_target;
  if (Gamma() / ba.H0 > 1.0) {
    // Approximately fully decayed today: compute correction factor.
    const double a_nr         = 3.15 / GetMass();
    const double Omega0_ur    = (ctx.all_species && ctx.all_species->count("UR"))
                                    ? ctx.all_species->at("UR")->GetOmega0()
                                    : 0.;
    const double k_rad        = std::sqrt(2.0 * ba.H0 * std::sqrt(ba.Omega0_g + Omega0_ur));
    const double t_nr         = std::pow(a_nr / k_rad, 2.0);
    const double x            = Gamma() * t_nr;
    const double experfcsqrtx = (x < 20.) ? std::exp(x) * std::erfc(std::sqrt(x))
                                          : 1.0 / std::sqrt(x * _PI_);
    Omega_ini                 = std::sqrt(2.0) * a_nr * std::sqrt(Gamma()) * Omega_target / k_rad /
                                (2.0 * std::sqrt(x) + std::sqrt(_PI_) * experfcsqrtx);
  }

  const double guess = (Omega_deg1 > 0.) ? Omega_ini / Omega_deg1 : 0.;
  const double dxdy  = (Omega_target != 0.) ? guess / Omega_target : 0.;
  return {guess, dxdy};
}

// ─────────────────────────────────────────────────────────────────────────────
// GetRescaledParameters — uses dq_ (not standard w_ weights)
// ─────────────────────────────────────────────────────────────────────────────

std::tuple<double, double> DNCDMSpecies::GetRescaledParameters(double a,
                                                               const double* lnf_array) const {
  double rho_scaled      = 0.;
  double p_scaled        = 0.;
  double pseudo_p_scaled = 0.;

  const double lnN = GetRescalingFactor(lnf_array);
  for (int index_q = 0; index_q < q_size(); index_q++) {
    double dq      = dq_[index_q];
    double q       = q_[index_q];
    double lnf     = lnf_array[index_q];
    double epsilon = std::sqrt(q * q + a * a * M_ * M_);

    rho_scaled      += dq * q * q * epsilon * std::exp(lnN + lnf);
    p_scaled        += dq * std::pow(q, 4) / 3. / epsilon * std::exp(lnN + lnf);
    pseudo_p_scaled += dq * std::pow(q * q / epsilon, 3) / 3. * std::exp(lnN + lnf);
  }
  if (rho_scaled == 0.)
    return {0., 0.};
  double pseudo_p_over_p = pseudo_p_scaled / p_scaled;
  double w               = p_scaled / rho_scaled;
  return {w, pseudo_p_over_p};
}

// ── Background ─────────────────────────────────────────────────────────────

void DNCDMSpecies::RegisterBackgroundIndices(int& index_bg) {
  index_bg_number_   = index_bg++;
  index_bg_rho_      = index_bg++;
  index_bg_p_        = index_bg++;
  index_bg_pseudo_p_ = index_bg++;

  index_bg_lnf_decay_dr1_   = index_bg;
  index_bg                 += q_size();
  index_bg_dlnfdlnq_decay_  = index_bg;
  index_bg                 += q_size();
  index_bg_dlnfdlnq_sep_    = index_bg;
  index_bg                 += q_size();
}

void DNCDMSpecies::RegisterIntegrationIndices(int& index_bi) {
  if (collision_owned_) {
    // Single f-per-bin variable; the ln f + separate-dlnfdlnq pair is dropped.
    index_bi_f_parent_  = index_bi;
    index_bi           += q_size();
    return;
  }
  index_bi_lnf_decay_dr1_            = index_bi;
  index_bi                          += q_size();
  index_bi_dlnfdlnq_separate_decay_  = index_bi;
  index_bi                          += q_size();
}

void DNCDMSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  double* pvecback_integration = ctx.pvecback_integration;
  if (collision_owned_) {
    // Seed f = w_/dq_ directly (the ln f + separate-dlnfdlnq seeds are dropped),
    // in the integrator's rescaled units (see kFScale).
    for (int index_q = 0; index_q < q_size(); index_q++)
      pvecback_integration[index_bi_f_parent_ + index_q] = kFScale * w_[index_q] / dq_[index_q];
    return;
  }
  for (int index_q = 0; index_q < q_size(); index_q++) {
    double q  = q_[index_q];
    double f0 = w_[index_q] / dq_[index_q];

    pvecback_integration[index_bi_lnf_decay_dr1_ + index_q] = std::log(f0);
    pvecback_integration[index_bi_dlnfdlnq_separate_decay_ + index_q] =
        -q * std::exp(q) / (std::exp(q) + 1);  // Fermi-Dirac
  }
}

void DNCDMSpecies::ComputeBackground(double a, const double* pvecback_B, double* pvecback) {
  double z       = 1. / a - 1.;
  const int q_sz = q_size();

  if (collision_owned_) {
    // Read f from the single f-variable, floor it, take ln f and spline it vs
    // ln q for dlnf/dlnq — publishing the SAME background columns the decay-only
    // path fills, so the perturbation surface is unchanged. The separate-dlnfdlnq
    // column has no integration variable here, so it is written 0.
    std::vector<double> lnf_dlnf(2 * q_sz);
    std::vector<double> ddlnf(q_sz);
    std::vector<double> lnq(q_sz);
    for (int i = 0; i < q_sz; i++) {
      // Back to physical units here, once: everything downstream of this line --
      // the published ln f column, w_bg_, and every moment built from them -- is in
      // ordinary occupation units and needs no knowledge of kFScale.
      const double f_raw = pvecback_B[index_bi_f_parent_ + i] / kFScale;
      // SOFT floor, max(f,0) + eps rather than max(f, eps).
      //
      // The soft part: a hard floor is a kink, and this column is cubic-splined in
      // tau by the perturbation module, so a solution that reaches the floor and
      // later climbs off it puts a STEP in the spline's input -- a step measured in
      // DECADES (a f = 1e-100 -> 1e-16 recovery is a jump of 193 in ln f), and the
      // ringing that follows is what turned P(k) into NaN before db149ad3. f + eps
      // is monotone and smooth, asymptotes to the same ln(eps) plateau, and agrees
      // with ln f to a part in 1e59 anywhere this species carries energy.
      //
      // The max(f,0): f + eps ALONE is not positivity-safe. max(f, eps) was, and
      // dropping it was a regression -- at Gamma=1e9 the integrator drives f
      // negative by more than eps in some bins and log() then returns NaN (353 of
      // them in the published parent columns). Clamping the negative excursion to
      // zero keeps the plateau and the smoothness while making the logarithm total.
      // Note the old hard floor did not FIX those excursions, it hid them.
      const double f = std::max(f_raw, 0.) + kFParentFloor;
      if (f_raw < 0. && bgm_ != nullptr && bgm_->StoringBackgroundTable() && !negative_f_rows_) {
        negative_f_rows_ = true;
        fprintf(stderr,
                "WARNING: species '%s': the decaying parent's occupation went NEGATIVE on a "
                "stored background row (bin %d, a=%g, f=%g). It is clamped to zero so the "
                "published ln f stays finite, but a negative occupation is unphysical and the "
                "background is not to be trusted where it happens.\n",
                name().c_str(),
                i,
                a,
                f_raw);
      }
      // Only the floor BINDING on a STORED row matters. Trial states dip below it
      // and are then thrown away, so testing them reports a configuration as broken
      // when the accepted solution stayed 75 decades clear of the floor.
      if (bgm_ != nullptr && bgm_->StoringBackgroundTable()) {
        if (f_raw < kFParentFloor) {
          floor_touched_ = true;
        }
        else if (floor_touched_ && !floor_left_) {
          floor_left_ = true;
          fprintf(stderr,
                  "WARNING: species '%s': the decaying parent's occupation went below "
                  "kFParentFloor (%g) on stored rows and was then repopulated by inverse "
                  "decays (bin %d, a=%g). The published ln f column is pinned at the floor "
                  "over that stretch, so the recovery's amplitude is set by the floor rather "
                  "than by the solution. The floor is too high for this configuration.\n",
                  name().c_str(),
                  kFParentFloor,
                  i,
                  a);
        }
      }
      lnf_dlnf[i] = std::log(f);
      lnq[i]      = std::log(q_[i]);
      w_bg_[i]    = f * dq_[i];
    }
    array_spline_table_lines(lnq.data(),
                             q_sz,
                             lnf_dlnf.data(),
                             1,
                             ddlnf.data(),
                             _SPLINE_EST_DERIV_);
    array_derive_spline(lnq.data(), q_sz, lnf_dlnf.data(), ddlnf.data(), 1, 0, q_sz);
    for (int i = 0; i < q_sz; i++) {
      pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf[i];
      pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf[q_sz + i];
      pvecback[index_bg_dlnfdlnq_sep_ + i]   = 0.;
    }
    double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
    ComputeMomenta(z, &number_ncdm, &rho_ncdm, &p_ncdm, nullptr, &pseudo_p_ncdm);
    pvecback[index_bg_number_]   = number_ncdm;
    pvecback[index_bg_rho_]      = rho_ncdm;
    pvecback[index_bg_p_]        = p_ncdm;
    pvecback[index_bg_pseudo_p_] = pseudo_p_ncdm;
    return;
  }

  std::vector<double> lnf_dlnf_array(2 * q_sz);
  std::vector<double> ddlnf_array(q_sz);
  std::vector<double> lnq(q_sz);

  // ln f is the integration variable, and it stays finite and smooth however far
  // the parent decays: the loss term is d ln f / d tau = -a^2 M Gamma / epsilon,
  // whose integral is ~ -Gamma * t in the non-relativistic regime, i.e. -1.4e6 by
  // today at Gamma = 1e8 (km/s/Mpc) -- large, but nowhere near the range of a
  // double, and analytic in q. So there is nothing to guard against here:
  //   * the spline of ln f vs ln q is well posed at any magnitude (only the
  //     DIFFERENCES across bins enter d ln f / d ln q, and they stay O(1..1e4));
  //   * w_bg_ = exp(ln f) * dq underflowing to exactly 0 is the CORRECT statement
  //     that the parent has decayed away -- rho, p and n then vanish, the composite
  //     stops sourcing decay radiation, and every consumer of a vanishing parent
  //     density is already guarded (StressEnergy/FillSources test rho > 0, and the
  //     perturbation moments divide out the exp(lnN) rescaling from
  //     GetRescalingFactor precisely so that delta = delta rho / rho survives it).
  //
  // This replaces a `ln f <= -460 (in ANY bin) => substitute the pristine
  // Fermi-Dirac f0 in EVERY bin` fallback inherited from the C implementation.
  // That fallback was the reason the integrated (dr_representation = integrated)
  // path was unusable above Gamma ~ 3e4 (the Gamma at which Gamma*t_0 first
  // reaches 460, i.e. 460*c/(1000*t_0) in the km/s/Mpc input units): it
  // RESURRECTED an already extinct parent at its full undecayed density, and the
  // composite's decay source a*Gamma*M*n then went on injecting decay radiation
  // for the rest of the run. Measured at Gamma = 1e5, m = 0.3 eV:
  //   * evolver = ndf15 (the default): no error is raised at all. The parent
  //     climbs back from rho/rho_tot = 1e-20 to 4e-2 at a = 0.48 and the decay
  //     radiation reaches 94% of the universe; the age comes out 5.5 Gyr instead
  //     of 13.85 and P(k) is wrong by a factor 29 (482 at Gamma = 1e6). Silent.
  //   * evolver = rkdp45: the same discontinuity throws the explicit stepper,
  //     rho_dr overshoots to -4.6e-6 and the run aborts on rho_crit <= 0.
  // It also left index_bg_dlnfdlnq_sep_ unwritten, publishing a stale row.
  for (int i = 0; i < q_sz; i++) {
    lnf_dlnf_array[i] = pvecback_B[index_bi_lnf_decay_dr1_ + i];
    lnq[i]            = std::log(q_[i]);
  }
  array_spline_table_lines(lnq.data(),
                           q_sz,
                           lnf_dlnf_array.data(),
                           1,
                           ddlnf_array.data(),
                           _SPLINE_EST_DERIV_);
  array_derive_spline(lnq.data(), q_sz, lnf_dlnf_array.data(), ddlnf_array.data(), 1, 0, q_sz);
  for (int i = 0; i < q_sz; i++) {
    pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf_array[i];
    pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf_array[q_sz + i];
    pvecback[index_bg_dlnfdlnq_sep_ + i]   = pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
    w_bg_[i]                               = std::exp(lnf_dlnf_array[i]) * dq_[i];
  }

  double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
  ComputeMomenta(z, &number_ncdm, &rho_ncdm, &p_ncdm, nullptr, &pseudo_p_ncdm);
  pvecback[index_bg_number_]   = number_ncdm;
  pvecback[index_bg_rho_]      = rho_ncdm;
  pvecback[index_bg_p_]        = p_ncdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_ncdm;
}

void DNCDMSpecies::BackgroundDerivs(double /*tau*/,
                                    const double* /*y*/,
                                    double* dy,
                                    const double* pvecback) {
  if (collision_owned_)
    return;  // the composite assigns dy[index_bi_f_parent_+i] from the kernel
  const double a      = pvecback[bgm_->index_bg_a_];
  const double M_ncdm = M_;
  const double Gamma  = Gamma_;
  for (int i = 0; i < q_size(); ++i) {
    const double q                            = q_[i];
    const double epsilon                      = std::sqrt(q * q + a * a * M_ncdm * M_ncdm);
    dy[index_bi_lnf_decay_dr1_ + i]           = -a * a * M_ncdm * Gamma / epsilon;
    dy[index_bi_dlnfdlnq_separate_decay_ + i] = a * a * M_ncdm * Gamma * q * q /
                                                std::pow(epsilon, 3);
  }
}

// ── Background column output ───────────────────────────────────────────────

void DNCDMSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  const std::string& nm = name();
  w.Add("(.)number_" + nm, 0.);
  w.Add("(.)rho_" + nm, 0.);
  w.Add("(.)p_" + nm, 0.);
  for (int i = 0; i < q_size(); i++) {
    const std::string suffix = "_" + nm + "[" + std::to_string(i) + "]";
    w.Add("lnf" + suffix, 0.);
    w.Add("dlnfdlnq" + suffix, 0.);
    w.Add("dlnfdlnq_separate" + suffix, 0.);
  }
}

void DNCDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  const std::string& nm = name();
  w.Add("(.)number_" + nm, pvecback[bg_number_index()]);
  w.Add("(.)rho_" + nm, Rho(pvecback));
  w.Add("(.)p_" + nm, P(pvecback));
  const int bg_lnf_idx          = bg_lnf_index();
  const int bg_dlnfdlnq_idx     = bg_dlnfdlnq_index();
  const int bg_dlnfdlnq_sep_idx = bg_dlnfdlnq_sep_index();
  for (int i = 0; i < q_size(); i++) {
    const std::string suffix = "_" + nm + "[" + std::to_string(i) + "]";
    w.Add("lnf" + suffix, pvecback[bg_lnf_idx + i]);
    w.Add("dlnfdlnq" + suffix, pvecback[bg_dlnfdlnq_idx + i]);
    w.Add("dlnfdlnq_separate" + suffix, pvecback[bg_dlnfdlnq_sep_idx + i]);
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DNCDMSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                               perturb_vector* pv,
                                               const precision* ppr,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);

  index_pt_psi0_ = index_pt;

  // DNCDM has no fluid approximation — always full Boltzmann hierarchy.
  layout.l_max  = ppr->l_max_ncdm;
  layout.q_size = q_size();

  layout.index_per_q.clear();
  layout.index_per_q.reserve(layout.q_size);
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  index_pt += layout.total_size();
}

void DNCDMSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                 double /*tau*/,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l.data();
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback.data();
  const double M_ncdm    = M_;
  const int lmax         = layout.l_max;

  // Metric driver. The two representations differ ONLY in this one quantity:
  //   Psi-space (decay-only): the dimensionless d ln f-bar / d ln q, from pvecback.
  //   F-space (collision-owned): the SAME hierarchy multiplied through by f-bar, so
  //   the driver becomes d f-bar / d ln q = f-bar * d ln f-bar / d ln q while the
  //   streaming terms keep their form (f-bar*(qk/eps)Psi_{l+-1} = (qk/eps)F_{l+-1}).
  //   Reconstructing that product is safe for the parent -- unlike the daughters,
  //   whose f-bar sits on a positivity floor where d ln f-bar / d ln q is
  //   meaningless -- because the parent integrates its occupation directly and its
  //   published ln f is monotone.
  // Two loops rather than a branch inside one: the decay-only body must stay
  // exactly the instruction sequence it has always been (a branch in there moves
  // P(k) in the last digit under -ffast-math), and a shared scratch buffer would be
  // a data race, since perturb_derivs runs concurrently across k-modes.
  if (collision_owned_) {
    for (int iq = 0; iq < layout.q_size; ++iq) {
      const double q          = q_[iq];
      const double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + iq] *
                                std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]);

      const double epsilon        = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
      const double qk_div_epsilon = k * q / epsilon;
      const int idx               = layout.index_per_q[iq];

      dy[idx]     = -qk_div_epsilon * y[idx + 1] + metric_continuity * dlnf0_dlnq / 3.;
      dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) -
                    epsilon * metric_euler / (3. * q * k) * dlnf0_dlnq;
      dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                    s_l[2] * metric_shear * 2. / 15. * dlnf0_dlnq;
      for (int l = 3; l < lmax; ++l)
        dy[idx + l] = qk_div_epsilon / (2. * l + 1.) *
                      (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
      dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] -
                       (1. + lmax) * k * cotKgen * y[idx + lmax];
    }
    return;
  }

  for (int iq = 0; iq < layout.q_size; ++iq) {
    const double q    = q_[iq];
    double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + iq];

    const double epsilon        = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
    const double qk_div_epsilon = k * q / epsilon;
    const int idx               = layout.index_per_q[iq];

    dy[idx]     = -qk_div_epsilon * y[idx + 1] + metric_continuity * dlnf0_dlnq / 3.;
    dy[idx + 1] = qk_div_epsilon / 3. * (y[idx] - 2. * s_l[2] * y[idx + 2]) -
                  epsilon * metric_euler / (3. * q * k) * dlnf0_dlnq;
    dy[idx + 2] = qk_div_epsilon / 5. * (2. * s_l[2] * y[idx + 1] - 3. * s_l[3] * y[idx + 3]) -
                  s_l[2] * metric_shear * 2. / 15. * dlnf0_dlnq;
    for (int l = 3; l < lmax; ++l)
      dy[idx + l] = qk_div_epsilon / (2. * l + 1.) *
                    (l * s_l[l] * y[idx + l - 1] - (l + 1.) * s_l[l + 1] * y[idx + l + 1]);
    dy[idx + lmax] = qk_div_epsilon * y[idx + lmax - 1] - (1. + lmax) * k * cotKgen * y[idx + lmax];
  }
}

void DNCDMSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                          double* y,
                                          const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const double* pvecback = ctx.ppw->pvecback.data();
  for (int index_q = 0; index_q < layout.q_size; ++index_q) {
    const int idx        = layout.index_per_q[index_q];
    const double q       = q_[index_q];
    const double epsilon = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);
    // F-space seeds F = f̄·Ψ_adiabatic, i.e. the same expressions with ∂ln f̄/∂ln q
    // swapped for ∂f̄/∂ln q (see PerturbDerivs). An already-empty bin then seeds
    // F = 0 by itself — the right statement, and one Ψ-space cannot make.
    const double dlnf0_dlnq = collision_owned_
                                  ? pvecback[index_bg_dlnfdlnq_decay_ + index_q] *
                                        std::exp(pvecback[index_bg_lnf_decay_dr1_ + index_q])
                                  : pvecback[index_bg_dlnfdlnq_decay_ + index_q];
    const int lmax          = layout.l_max;

    y[idx + 0] = -0.25 * ctx.delta_ur * dlnf0_dlnq;
    if (lmax >= 1)
      y[idx + 1] = -epsilon / 3. / q / ctx.k * ctx.theta_ur * dlnf0_dlnq;
    if (lmax >= 2)
      y[idx + 2] = -0.5 * ctx.shear_ur * dlnf0_dlnq;
    if (lmax >= 3)
      y[idx + 3] = -0.25 * ctx.l3_ur * dlnf0_dlnq;
  }
}

// Fused stress-energy: avoids calling RescaledPerturbations three times.
// DeltaRho/RhoPlusPTheta/RhoPlusPShear share one RescaledPerturbations call;
// DeltaP uses its own loop (different weights: exp(lnf)*dq vs rescaled),
// matching the original per-term expressions and operand order exactly.
BaseSpecies::StressEnergyContribution DNCDMSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho = pvecback[index_bg_rho_];
  se.p   = pvecback[index_bg_p_];
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return se;

  // DeltaRho, RhoPlusPTheta, RhoPlusPShear via RescaledPerturbations (one call).
  const double a      = ppw->scalar_ctx.a;
  const double k      = ppw->scalar_ctx.k;
  auto [d, t, s]      = RescaledPerturbations(layout, a, k, ppw);
  se.delta_rho        = Rho(pvecback) * d;
  se.rho_plus_p_theta = (Rho(pvecback) + P(pvecback)) * t;
  se.rho_plus_p_shear = (Rho(pvecback) + P(pvecback)) * s;

  // DeltaP: independent loop with exp(lnf)*dq weights (matches DeltaP method).
  // Both representations integrate the same δf; only where it comes from differs
  // (f̄·Ψ in Ψ-space, the state variable itself in F-space), so the F-space weight
  // simply drops the exp(lnf). Split into two loops rather than branching inside
  // one: the decay-only body must stay exactly the instruction sequence it was.
  double delta_p_ncdm = 0.0;
  if (collision_owned_) {
    for (int iq = 0; iq < layout.q_size; ++iq) {
      const double q        = q_[iq];
      const double epsilon  = std::sqrt(q * q + std::pow(M_ * a, 2));
      delta_p_ncdm         += q * q * q * q / epsilon * dq_[iq] * y[layout.index_per_q[iq]];
    }
  }
  else {
    for (int iq = 0; iq < layout.q_size; ++iq) {
      double w0             = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * dq_[iq];
      const double q        = q_[iq];
      const double epsilon  = std::sqrt(q * q + std::pow(M_ * a, 2));
      delta_p_ncdm         += q * q * q * q / epsilon * w0 * y[layout.index_per_q[iq]];
    }
  }
  const double fac = factor_ * std::pow(1. / a, 4);
  se.delta_p       = delta_p_ncdm * fac / 3.;

  return se;
}

// Transfer sources (#309 slots, allocated by NCDMBaseSpecies but until now unfilled).
//
// Deliberately spelled like DrPsdSpecies::FillSources rather than NCDMSpecies::
// FillSources: the two are algebraically identical under the base RhoDotOverRho
// (-3ℋ(ρ+p)/ρ vs 3ℋ(1+w), the same number), but the parent and its two daughters are
// summed against each other downstream to form the sector's total δ, and a sector sum
// is only meaningful if its three terms carry the SAME gauge convention. Should the
// decay sink ever be folded into RhoDotOverRho, this spelling picks it up and the
// NCDM one silently would not.
//
// StressEnergy already resolves the representation split: it returns δρ from
// RescaledPerturbations, which handles F-space (collision-owned) and Ψ-space alike and
// carries the exp(lnN) rescaling that keeps δ = δρ/ρ meaningful after f̄ has decayed
// below underflow. Nothing here needs to know which mode is active.
void DNCDMSpecies::FillSources(const BaseSpecies::PerturbLayout& layout,
                               const double* /*y*/,
                               const double* /*dy*/,
                               PerturbSourceContext& ctx) const {
  if (ctx.index_md != ctx.p_mod->index_md_scalars_)
    return;
  if (index_tp_delta_ < 0 && index_tp_theta_ < 0)
    return;

  perturb_workspace* ppw   = ctx.ppw;
  const double* pvecback   = ppw->pvecback.data();
  const perturb_vector* pv = ppw->pv.get();
  const auto se            = StressEnergy(layout, pv, pv->y.data(), pvecback, ppw);

  if (index_tp_delta_ >= 0) {
    // ρ_H falls to ~1e-97 of ρ_tot at high Γ, so guard the ratio rather than assume
    // it: both sums vanish together and the quotient is only well posed while ρ > 0.
    const double src = (se.rho > 0.)
                           ? se.delta_rho / se.rho -
                                 RhoDotOverRho(pvecback, ctx.a_prime_over_a) * ctx.theta_over_k2
                           : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md,
                              ctx.index_ic,
                              index_tp_delta_,
                              ctx.index_tau,
                              ctx.index_k,
                              src);
  }
  if (index_tp_theta_ >= 0) {
    const double rho_plus_p = se.rho + se.p;
    const double src = (rho_plus_p > 0.) ? se.rho_plus_p_theta / rho_plus_p + ctx.theta_shift : 0.;
    ctx.p_mod->SetSourceValue(ctx.index_md,
                              ctx.index_ic,
                              index_tp_theta_,
                              ctx.index_tau,
                              ctx.index_k,
                              src);
  }
}

void DNCDMSpecies::WriteOutputColumns(PerturbColumnWriter& w,
                                      const PerturbationsModule& mod,
                                      file_format fmt,
                                      BaseSpecies::TransferColumnSection section) const {
  if (fmt != file_format::class_format)
    return;
  const perturbs* ppt = mod.GetPerturbs();
  if (section != TransferColumnSection::velocity && ppt->has_density_transfers)
    w.Add("d_ncdm_" + name(), index_tp_delta_, index_tp_delta_ >= 0);
  if (section != TransferColumnSection::density && ppt->has_velocity_transfers)
    w.Add("t_ncdm_" + name(), index_tp_theta_, index_tp_theta_ >= 0);
}

std::tuple<double, double, double> DNCDMSpecies::RescaledPerturbations(
    const NCDMBaseSpecies::PerturbLayout& layout,
    double a,
    double k,
    const perturb_workspace* ppw) const {
  double rho_scaled              = 0.;
  double rho_plus_p_scaled       = 0.;
  double rho_delta_scaled        = 0.;
  double rho_plus_p_theta_scaled = 0.;
  double rho_plus_p_shear_scaled = 0.;

  const double* lnf_array = ppw->pvecback.data() + bg_lnf_index();

  if (collision_owned_) {
    // F-space. Same moments, but the state variable already IS δf, so the
    // perturbation weight loses its exp(lnf) and keeps only the rescaling.
    //
    // That rescaling exists so a f̄ which has decayed below the underflow threshold
    // still yields a meaningful δ = δρ/ρ: both sums carry exp(lnN) and it cancels in
    // the ratio. In Ψ-space it only ever appears as exp(lnN + lnf) and is bounded by
    // construction; here the numerator needs the bare exp(lnN), which would overflow
    // once f̄_max drops below ~e^-309. Clamp it — the factor is common to both sums,
    // so capping costs headroom at the small end and nothing else, and 700 still
    // leaves the peak bin ~e^400 of room.
    const double lnN = std::min(GetRescalingFactor(lnf_array), 700.);
    const double eN  = std::exp(lnN);
    for (int index_q = 0; index_q < q_size(); index_q++) {
      const int index_pt   = layout.index_per_q[index_q];
      const double dq_val  = dq()[index_q];
      const double lnf     = lnf_array[index_q];
      const double q       = q_[index_q];
      const double q2      = q * q;
      const double epsilon = std::sqrt(q2 + a * a * M_ * M_);

      rho_scaled              += dq_val * q2 * epsilon * std::exp(lnN + lnf);
      rho_plus_p_scaled       += dq_val * q2 * (epsilon + q2 / 3. / epsilon) * std::exp(lnN + lnf);
      rho_delta_scaled        += dq_val * q2 * epsilon * eN * ppw->pv->y[index_pt];
      rho_plus_p_theta_scaled += dq_val * q2 * q * eN * ppw->pv->y[index_pt + 1];
      rho_plus_p_shear_scaled += dq_val * q2 * q2 / epsilon * eN * ppw->pv->y[index_pt + 2];
    }
  }
  else {
    const double lnN = GetRescalingFactor(lnf_array);
    for (int index_q = 0; index_q < q_size(); index_q++) {
      const int index_pt   = layout.index_per_q[index_q];
      const double dq_val  = dq()[index_q];
      const double lnf     = lnf_array[index_q];
      const double q       = q_[index_q];
      const double q2      = q * q;
      const double epsilon = std::sqrt(q2 + a * a * M_ * M_);

      rho_scaled              += dq_val * q2 * epsilon * std::exp(lnN + lnf);
      rho_plus_p_scaled       += dq_val * q2 * (epsilon + q2 / 3. / epsilon) * std::exp(lnN + lnf);
      rho_delta_scaled        += dq_val * q2 * epsilon * std::exp(lnN + lnf) * ppw->pv->y[index_pt];
      rho_plus_p_theta_scaled += dq_val * q2 * q * std::exp(lnN + lnf) * ppw->pv->y[index_pt + 1];
      rho_plus_p_shear_scaled += dq_val * q2 * q2 / epsilon * std::exp(lnN + lnf) *
                                 ppw->pv->y[index_pt + 2];
    }
  }
  rho_plus_p_theta_scaled *= k;
  rho_plus_p_shear_scaled *= 2. / 3.;

  const double delta = rho_delta_scaled / rho_scaled;
  const double theta = rho_plus_p_theta_scaled / rho_plus_p_scaled;
  const double shear = rho_plus_p_shear_scaled / rho_plus_p_scaled;

  return {delta, theta, shear};
}
