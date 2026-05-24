#include "dncdm_species.h"

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "arrays.h"
#include "background_column_writer.h"
#include "background_module.h"
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
    std::string unused;
    if (pfc.read_string(key, unused)) {
      throw std::invalid_argument(
          std::string("'") + key + "' is no longer supported. Use dot-syntax: " +
          "'<instance>.<dot-name> = ...' with '<instance>.type = ncdm_decay_dr'.");
    }
  }
}

void ApplyDncdmInitialClosure(DNCDMSpecies& sp, const SpeciesBuildContext& ctx) {
  const auto& Omega_ini = sp.Omega_ini_pending();
  const auto& Neff_ini  = sp.Neff_ini_pending();
  if (!Omega_ini.has_value() && !Neff_ini.has_value()) {
    return;
  }

  const double a_ini_seed = ctx.pba ? ctx.pba->a_today * 1e-14 : 1e-14;
  const double a_today    = ctx.pba ? ctx.pba->a_today : 1.;
  const double a_ini      = sp.GetIni(a_ini_seed, a_today, ctx.ncdm_settings->tol_ncdm);
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

  // Read decay rate (exactly one of these must be specified)
  double Gamma_value         = 0.;
  double log10Gamma_value    = 0.;
  double lifetime_value      = 0.;
  double log10lifetime_value = 0.;
  bool has_G                 = input.read_double("Gamma", Gamma_value);
  bool has_lG                = input.read_double("log10Gamma", log10Gamma_value);
  bool has_lt                = input.read_double("lifetime", lifetime_value);
  bool has_llt               = input.read_double("log10lifetime", log10lifetime_value);

  int n_provided = (int) has_G + (int) has_lG + (int) has_lt + (int) has_llt;
  if (n_provided != 1) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': specify exactly one of Gamma, log10Gamma, lifetime, log10lifetime");
  }

  double Gamma_raw = 0.;
  if (has_G) {
    Gamma_raw = Gamma_value;
  }
  else if (has_lG) {
    Gamma_raw = std::pow(10., log10Gamma_value);
  }
  else if (has_lt) {
    Gamma_raw = 1. / lifetime_value / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  else {  // has_llt
    double lifetime = std::pow(10., log10lifetime_value);
    Gamma_raw       = 1. / lifetime / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  Gamma_ = Gamma_raw * (1.e3 / _c_);

  // Read DNCDM-specific Omega / deg / Omega_ini values from the same instance.
  // These were previously bulk-loaded in CreateAll from CSV lists.
  // Semantics: at most one of {Omega, omega} and at most one of
  // {deg, Omega_ini, omega_ini, Neff_ini} (the latter group resolved later
  // by ConstructSpecies for shooting).
  double Omega_local = 0., omega_local = 0.;
  bool has_Omega = input.read_double("Omega", Omega_local);
  bool has_omega = input.read_double("omega", omega_local);
  if (has_Omega && has_omega) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of Omega, omega");
  }
  if (has_omega) {
    Omega_local = omega_local / settings.h / settings.h;
    has_Omega   = true;
  }
  if (has_Omega) {
    SetOmega0(Omega_local, settings.h);
  }

  double deg_local           = 0.;
  double Omega_ini_local     = 0.;
  double omega_ini_local     = 0.;
  double Neff_ini_local      = 0.;
  double Omega_dncdmdr_local = 0.;
  double omega_dncdmdr_local = 0.;
  bool has_deg               = input.read_double("deg", deg_local);
  bool has_Omega_ini         = input.read_double("Omega_ini", Omega_ini_local);
  bool has_omega_ini         = input.read_double("omega_ini", omega_ini_local);
  bool has_Neff_ini          = input.read_double("Neff_ini", Neff_ini_local);
  bool has_Omega_dncdmdr     = input.read_double("Omega_dncdmdr", Omega_dncdmdr_local);
  bool has_omega_dncdmdr     = input.read_double("omega_dncdmdr", omega_dncdmdr_local);
  if (has_Omega_dncdmdr && has_omega_dncdmdr) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of Omega_dncdmdr, omega_dncdmdr");
  }
  if (has_omega_dncdmdr) {
    Omega_dncdmdr_local = omega_dncdmdr_local / settings.h / settings.h;
    has_Omega_dncdmdr   = true;
  }

  // Omega_dncdmdr is mutually exclusive with Omega_ini/omega_ini/Neff_ini/Omega/omega.
  // deg may co-occur with Omega_dncdmdr (shooting-iteration case: DoShooting wrote deg).
  if (has_Omega_dncdmdr && (has_Omega_ini || has_omega_ini || has_Neff_ini || has_Omega)) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': Omega_dncdmdr conflicts with Omega_ini/omega_ini/Neff_ini/Omega/omega");
  }

  // Count mutually-exclusive target specs (Omega_dncdmdr counts as one; deg does not).
  int n_deg_options = (int) has_Omega_ini + (int) has_omega_ini + (int) has_Neff_ini +
                      (int) has_Omega_dncdmdr;
  if (n_deg_options > 1) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': specify at most one of Omega_ini, omega_ini, Neff_ini, Omega_dncdmdr");
  }

  if (has_deg) {
    // deg wins: used directly (also covers shooting-iteration case where DoShooting set deg).
    SetDegAndFactor(deg_local);
  }
  else if (has_Omega_dncdmdr) {
    // Today-density target: stash the target and set a guessed deg via DegGuessFromOmegaToday.
    // DegGuessFromOmegaToday calls GetDeg() — deg_ is initialized to 1. by NCDMBaseSpecies,
    // so the per-unit-degeneracy calculation is well-defined without a prior SetDegAndFactor.
    Omega_dncdmdr_pending_ = Omega_dncdmdr_local;
    if (pba) {
      // Build context needed for GetIni; assemble a minimal one from available data.
      SpeciesBuildContext guess_ctx{nullptr, pba, nullptr, &settings, nullptr};
      auto [g, d] = DegGuessFromOmegaToday(guess_ctx, Omega_dncdmdr_local);
      SetDegAndFactor(g);
    }
    // If pba is null (early-construction pass), leave deg_ at 1. and let the shooting
    // dispatch provide the correct value when DoShooting writes <instance>.deg.
  }
  else {
    // The Omega_ini / omega_ini / Neff_ini cases require an a_ini-driven shooting step.
    // Preserve their values by stashing them on the species so that ConstructSpecies'
    // closure pass can apply them.
    if (has_Omega_ini) {
      Omega_ini_pending_ = Omega_ini_local;  // already in Omega units
    }
    else if (has_omega_ini) {
      Omega_ini_pending_ = omega_ini_local / settings.h / settings.h;
    }
    else if (has_Neff_ini) {
      Neff_ini_pending_ = Neff_ini_local;
    }
  }

  // Compute dq_[i] = w_bg_[i] / f0(q_bg_[i]) (matches existing semantics)
  dq_ = ComputeDq();
}

// ─────────────────────────────────────────────────────────────────────────────
// CreateAll factory
// ─────────────────────────────────────────────────────────────────────────────

std::vector<DNCDMSpecies::Named> DNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyDecayDrKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_decay_dr");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);  // mark consumed

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
  double a_ini = ctx.ppr ? ctx.ppr->a_ini_over_a_today_default * ba.a_today : ba.a_today * 1e-14;
  const double tol_ini = ctx.ppr ? ctx.ppr->tol_ncdm_initial_w : ctx.ncdm_settings->tol_ncdm;
  a_ini                = GetIni(a_ini, ba.a_today, tol_ini);
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
    const double k_rad        = std::sqrt(2.0 * ba.H0 * std::sqrt(ba.Omega0_g + ba.Omega0_ur));
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
  index_bi_lnf_decay_dr1_            = index_bi;
  index_bi                          += q_size();
  index_bi_dlnfdlnq_separate_decay_  = index_bi;
  index_bi                          += q_size();
}

void DNCDMSpecies::SetBackgroundInitialConditions(double /*a_rel*/, double* pvecback_integration) {
  for (int index_q = 0; index_q < q_size(); index_q++) {
    double q  = q_[index_q];
    double f0 = w_[index_q] / dq_[index_q];

    pvecback_integration[index_bi_lnf_decay_dr1_ + index_q] = std::log(f0);
    pvecback_integration[index_bi_dlnfdlnq_separate_decay_ + index_q] =
        -q * std::exp(q) / (std::exp(q) + 1);  // Fermi-Dirac
  }
}

void DNCDMSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  double z       = 1. / a_rel - 1.;
  const int q_sz = q_size();

  std::vector<double> lnf_dlnf_array(2 * q_sz);
  std::vector<double> ddlnf_array(q_sz);
  std::vector<double> lnq(q_sz);

  if (pba_->background_method == bgevo_evolver) {
    bool has_problem = false;
    for (int i = 0; i < q_sz; i++) {
      if (pvecback_B[index_bi_lnf_decay_dr1_ + i] <= -460) {
        has_problem = true;
        break;
      }
      lnf_dlnf_array[i] = pvecback_B[index_bi_lnf_decay_dr1_ + i];
      lnq[i]            = std::log(q_[i]);
    }
    if (has_problem) {
      for (int i = 0; i < q_sz; i++) {
        double q                               = q_[i];
        double lnf_fd                          = -std::log(1. + std::exp(q));
        double dlnfdlnq_fd                     = -q * std::exp(q) / (1. + std::exp(q));
        pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_fd;
        pvecback[index_bg_dlnfdlnq_decay_ + i] = dlnfdlnq_fd;
        w_bg_[i]                               = std::exp(lnf_fd) * dq_[i];
      }
    }
    else {
      class_call(array_spline_table_lines(lnq.data(),
                                          q_sz,
                                          lnf_dlnf_array.data(),
                                          1,
                                          ddlnf_array.data(),
                                          _SPLINE_EST_DERIV_,
                                          bgm_->error_message_),
                 bgm_->error_message_,
                 bgm_->error_message_);
      class_call(array_derive_spline(lnq.data(),
                                     q_sz,
                                     lnf_dlnf_array.data(),
                                     ddlnf_array.data(),
                                     1,
                                     0,
                                     q_sz,
                                     bgm_->error_message_),
                 bgm_->error_message_,
                 bgm_->error_message_);
      for (int i = 0; i < q_sz; i++) {
        pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf_array[i];
        pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf_array[q_sz + i];
        pvecback[index_bg_dlnfdlnq_sep_ + i]   = pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
        w_bg_[i]                               = std::exp(lnf_dlnf_array[i]) * dq_[i];
      }
    }
  }
  else {
    for (int i = 0; i < q_sz; i++) {
      lnf_dlnf_array[i] = pvecback_B[index_bi_lnf_decay_dr1_ + i];
      lnq[i]            = std::log(q_[i]);
    }
    class_call(array_spline_table_lines(lnq.data(),
                                        q_sz,
                                        lnf_dlnf_array.data(),
                                        1,
                                        ddlnf_array.data(),
                                        _SPLINE_EST_DERIV_,
                                        bgm_->error_message_),
               bgm_->error_message_,
               bgm_->error_message_);
    class_call(array_derive_spline(lnq.data(),
                                   q_sz,
                                   lnf_dlnf_array.data(),
                                   ddlnf_array.data(),
                                   1,
                                   0,
                                   q_sz,
                                   bgm_->error_message_),
               bgm_->error_message_,
               bgm_->error_message_);
    for (int i = 0; i < q_sz; i++) {
      pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf_array[i];
      pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf_array[q_sz + i];
      pvecback[index_bg_dlnfdlnq_sep_ + i]   = pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
      w_bg_[i]                               = std::exp(lnf_dlnf_array[i]) * dq_[i];
    }
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
                                 const perturb_parameters_and_workspace& ppaw) {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const double* s_l               = ppw->s_l;
  const double k                  = ctx.k;
  const double a2                 = ctx.a2;
  const double metric_continuity  = ctx.metric_continuity;
  const double metric_euler       = ctx.metric_euler;
  const double metric_shear       = ctx.metric_shear;
  const double cotKgen            = ctx.cotKgen;

  const double* pvecback = ppw->pvecback;
  const double M_ncdm    = M_;
  const int lmax         = layout.l_max;

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

void DNCDMSpecies::PerturbDerivs(double /*tau*/,
                                 const double* /*y*/,
                                 double* /*dy*/,
                                 const perturb_parameters_and_workspace& /*ppaw*/) {
  // Legacy non-layout PerturbDerivs is unreachable: the composite
  // (DNCDM_DR_Species) calls the layout-based variant directly with my.dncdm.
  throw std::logic_error(
      "DNCDMSpecies::PerturbDerivs(tau,...) is unreachable — composite must pass layout");
}

void DNCDMSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                          double* y,
                                          const PerturbIcContext& ctx) {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return;

  const double* pvecback = ctx.ppw->pvecback;
  for (int index_q = 0; index_q < layout.q_size; ++index_q) {
    const int idx           = layout.index_per_q[index_q];
    const double q          = q_[index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);
    const double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + index_q];
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

// ── Integrated observables (layout-based) ──────────────────────────────────

double DNCDMSpecies::Delta(const BaseSpecies::PerturbLayout& base,
                           const perturb_vector* /*pv*/,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return 0.;
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(layout, a, k, ppw);
  return d;
}

double DNCDMSpecies::Theta(const BaseSpecies::PerturbLayout& base,
                           const perturb_vector* /*pv*/,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return 0.;
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(layout, a, k, ppw);
  return t;
}

double DNCDMSpecies::DeltaP(const BaseSpecies::PerturbLayout& base,
                            const perturb_vector* /*pv*/,
                            const double* y,
                            const double* pvecback,
                            const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return 0.;

  const double a      = ppw->scalar_ctx.a;
  double delta_p_ncdm = 0.0;
  for (int iq = 0; iq < layout.q_size; ++iq) {
    double w0             = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * dq_[iq];
    const double q        = q_[iq];
    const double epsilon  = std::sqrt(q * q + std::pow(M_ * a, 2));
    delta_p_ncdm         += q * q * q * q / epsilon * w0 * y[layout.index_per_q[iq]];
  }
  const double fac = factor_ * std::pow(pba_->a_today / a, 4);
  return delta_p_ncdm * fac / 3.;
}

double DNCDMSpecies::RhoPlusPShear(const BaseSpecies::PerturbLayout& base,
                                   const perturb_vector* /*pv*/,
                                   const double* /*y*/,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  if (layout.q_size <= 0 || layout.index_per_q.empty())
    return 0.;
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(layout, a, k, ppw);
  return (Rho(pvecback) + P(pvecback)) * s;
}

// ── Integrated observables (legacy — composite always passes layout via the
// layout-based overrides, so these legacy overloads are unreachable on the
// real call paths and exist only to satisfy the pure-virtual contract from
// BaseSpecies). The composite (DNCDM_DR_Species) overrides the layout-based
// Delta/Theta/etc. and forwards `my.dncdm` to the layout-based DNCDMSpecies
// overloads above. CompositeSpecies's legacy fallback chains through here
// only when no caller knows the layout — no such call site exists today. ──

double DNCDMSpecies::Delta(const perturb_vector* /*pv*/,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  throw std::logic_error("DNCDMSpecies::Delta(pv,...) is unreachable — composite must pass layout");
}

double DNCDMSpecies::Theta(const perturb_vector* /*pv*/,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* /*ppw*/) const {
  throw std::logic_error("DNCDMSpecies::Theta(pv,...) is unreachable — composite must pass layout");
}

double DNCDMSpecies::DeltaP(const perturb_vector* /*pv*/,
                            const double* /*y*/,
                            const double* /*pvecback*/,
                            const perturb_workspace* /*ppw*/) const {
  throw std::logic_error(
      "DNCDMSpecies::DeltaP(pv,...) is unreachable — composite must pass layout");
}

double DNCDMSpecies::RhoPlusPShear(const perturb_vector* /*pv*/,
                                   const double* /*y*/,
                                   const double* /*pvecback*/,
                                   const perturb_workspace* /*ppw*/) const {
  throw std::logic_error(
      "DNCDMSpecies::RhoPlusPShear(pv,...) is unreachable — composite must pass layout");
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

  const double* lnf_array = ppw->pvecback + bg_lnf_index();
  const double lnN        = GetRescalingFactor(lnf_array);

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
  rho_plus_p_theta_scaled *= k;
  rho_plus_p_shear_scaled *= 2. / 3.;

  const double delta = rho_delta_scaled / rho_scaled;
  const double theta = rho_plus_p_theta_scaled / rho_plus_p_scaled;
  const double shear = rho_plus_p_shear_scaled / rho_plus_p_scaled;

  return {delta, theta, shear};
}
