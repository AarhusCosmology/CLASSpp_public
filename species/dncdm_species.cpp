#include "dncdm_species.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "arrays.h"
#include "background_column_writer.h"
#include "background_module.h"
#include "perturbations_module.h"
#include "species/species_input.h"

namespace {

int readDoubleList(
    FileContent* pfc, const char* name, std::vector<double>& values, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_list_of_doubles(name, values) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}

// ── Dot-syntax → legacy key mapping for ncdm_decay_dr species ──────────────
// Each dot-syntax instance field (under nu*.…) maps to the corresponding legacy
// CSV-style key which DNCDMSpecies / NCDMBaseSpecies already knows how to read.
//   "m"                   -> "m_ncdm_decay_dr"           // mass in eV
//   "T"                   -> "T_ncdm_decay_dr"           // temperature ratio
//   "ksi"                 -> "ksi_ncdm_decay_dr"         // degeneracy parameter
//   "deg"                 -> "deg_ncdm_decay_dr"         // phase-space degeneracy
//   "quadrature_strategy" -> "quadrature_strategy_ncdm_decay_dr"
//   "momenta_bins"        -> "N_momentum_bins_ncdm_decay_dr"
//   "max_q"               -> "maximum_q_ncdm_decay_dr"
//   "Gamma"               -> "Gamma_ncdm_decay_dr"       // decay rate
//   "log10Gamma"          -> "log10Gamma_ncdm_decay_dr"
//   "lifetime"            -> "lifetime_ncdm_decay_dr"    // (singularised)
//   "log10lifetime"       -> "log10lifetime_ncdm_decay_dr"
//   "Omega"               -> "Omega_dncdmdr"             // today density
//   "omega"               -> "omega_dncdmdr"             // today density * h^2
//   "Omega_ini"           -> "Omega_ini_dncdm"
//   "omega_ini"           -> "omega_ini_dncdm"
//   "Neff_ini"            -> "Neff_ini_dncdm"
struct FieldMap {
  const char* dot;
  const char* legacy;
  const char* default_value;
};
constexpr FieldMap kDecayDrFields[] = {
    {"m", "m_ncdm_decay_dr", "1.0"},
    {"T", "T_ncdm_decay_dr", "0.71611"},
    {"ksi", "ksi_ncdm_decay_dr", "0"},
    {"deg", "deg_ncdm_decay_dr", "1.0"},
    {"quadrature_strategy", "quadrature_strategy_ncdm_decay_dr", "0"},
    {"momenta_bins", "N_momentum_bins_ncdm_decay_dr", "5"},
    {"max_q", "maximum_q_ncdm_decay_dr", "15"},
    {"Gamma", "Gamma_ncdm_decay_dr", "0"},
    {"log10Gamma", "log10Gamma_ncdm_decay_dr", "0"},
    {"lifetime", "lifetime_ncdm_decay_dr", "0"},
    {"log10lifetime", "log10lifetime_ncdm_decay_dr", "0"},
    {"Omega", "Omega_dncdmdr", "0"},
    {"omega", "omega_dncdmdr", "0"},
    {"Omega_ini", "Omega_ini_dncdm", "0"},
    {"omega_ini", "omega_ini_dncdm", "0"},
    {"Neff_ini", "Neff_ini_dncdm", "0"},
};

std::vector<std::string> synthesise_decay_dr_ncdm_flat_keys(
    FileContent& pfc, const std::vector<std::string>& instances) {
  const int n_fields = static_cast<int>(sizeof(kDecayDrFields) / sizeof(FieldMap));
  (void) CollectInstanceFieldValues(&pfc, instances, "type");
  std::vector<std::vector<std::string>> values;
  values.reserve(n_fields);
  for (int f = 0; f < n_fields; ++f) {
    values.push_back(CollectInstanceFieldValues(&pfc, instances, kDecayDrFields[f].dot));
  }

  pfc.set("N_ncdm_decay_dr", std::to_string(instances.size()));
  for (int f = 0; f < n_fields; ++f) {
    if (!AnyInstanceFieldValue(values[f]))
      continue;
    pfc.set(kDecayDrFields[f].legacy, CsvWithDefaults(values[f], kDecayDrFields[f].default_value));
  }
  return instances;
}

bool has_unconsumed_dot_type_keys(const FileContent& pfc,
                                  const std::vector<std::string>& instances) {
  for (const auto& instance : instances) {
    if (!pfc.was_read(instance + ".type"))
      return true;
  }
  return false;
}

}  // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────

DNCDMSpecies::DNCDMSpecies(FileContent* pfc,
                           int ncdm_index,
                           int dncdm_index,
                           const NcdmSettings& settings,
                           const background* pba,
                           const BackgroundModule* bgm)
    : NCDMBaseSpecies("DNCDM_" + std::to_string(dncdm_index),
                      EnergyType::Other,
                      pfc,
                      dncdm_index,
                      settings,
                      true /*is_decay_dr*/),
      ncdm_id_(ncdm_index), pba_(pba) {
  bgm_ = bgm;

  // Read Gamma for this decay_dr species.
  // Supports: Gamma_ncdm_decay_dr, log10Gamma_ncdm_decay_dr,
  //           lifetime_ncdm_decay_dr, log10lifetime_ncdm_decay_dr
  char errmsg[2048];
  int Gamma_read, log10Gamma_read, lifetime_read, log10lifetime_read;
  std::vector<double> Gamma_list, log10Gamma_list, lifetime_list, log10lifetime_list;

  readDoubleList(pfc, "Gamma_ncdm_decay_dr", Gamma_list, &Gamma_read, errmsg);
  readDoubleList(pfc, "log10Gamma_ncdm_decay_dr", log10Gamma_list, &log10Gamma_read, errmsg);
  readDoubleList(pfc, "lifetime_ncdm_decay_dr", lifetime_list, &lifetime_read, errmsg);
  readDoubleList(pfc,
                 "log10lifetime_ncdm_decay_dr",
                 log10lifetime_list,
                 &log10lifetime_read,
                 errmsg);

  // Exactly one decay-rate parameter must be provided
  int n_provided = (Gamma_read == _TRUE_) + (log10Gamma_read == _TRUE_) +
                   (lifetime_read == _TRUE_) + (log10lifetime_read == _TRUE_);
  if (n_provided > 1)
    throw std::invalid_argument(
        "DNCDM species: specify exactly one of Gamma_ncdm_decay_dr, "
        "log10Gamma_ncdm_decay_dr, lifetime_ncdm_decay_dr, log10lifetime_ncdm_decay_dr");

  auto check_size = [&](const char* name, const std::vector<double>& list) {
    if (dncdm_index >= static_cast<int>(list.size()))
      throw std::invalid_argument(std::string(name) +
                                  " has fewer entries than expected for DNCDM species index " +
                                  std::to_string(dncdm_index));
  };

  double Gamma_raw = 0.;
  if (Gamma_read == _TRUE_) {
    check_size("Gamma_ncdm_decay_dr", Gamma_list);
    Gamma_raw = Gamma_list[dncdm_index];
  }
  else if (log10Gamma_read == _TRUE_) {
    check_size("log10Gamma_ncdm_decay_dr", log10Gamma_list);
    Gamma_raw = std::pow(10., log10Gamma_list[dncdm_index]);
  }
  else if (lifetime_read == _TRUE_) {
    check_size("lifetime_ncdm_decay_dr", lifetime_list);
    Gamma_raw = 1. / lifetime_list[dncdm_index] / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  else if (log10lifetime_read == _TRUE_) {
    check_size("log10lifetime_ncdm_decay_dr", log10lifetime_list);
    double lifetime = std::pow(10., log10lifetime_list[dncdm_index]);
    Gamma_raw       = 1. / lifetime / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  // Convert to Mpc^-1
  Gamma_ = Gamma_raw * (1.e3 / _c_);

  // Compute dq_[i] = w_bg_[i] / f0(q_bg_[i])
  // (w_bg_ is already populated by NCDMBaseSpecies::InitQuadrature)
  dq_ = ComputeDq();
}

// ─────────────────────────────────────────────────────────────────────────────
// CreateAll factory
// ─────────────────────────────────────────────────────────────────────────────

std::vector<DNCDMSpecies::Named> DNCDMSpecies::CreateAll(FileContent* pfc,
                                                         const NcdmSettings& settings,
                                                         const background* pba,
                                                         const BackgroundModule* bgm) {
  std::vector<Named> result;
  ErrorMsg errmsg_buf;
  char* errmsg = errmsg_buf;

  const std::vector<std::string> dot_instances = pfc->instances_with("type", "ncdm_decay_dr");

  // Find N_ncdm_standard and N_ncdm_decay_dr
  int N_ncdm_standard = 0, N_ncdm_decay_dr = 0;
  int flag1, flag2, flag3;

  int int1, int2;
  parser_read_int(pfc, "N_ncdm_standard", &int1, &flag1, errmsg);
  parser_read_int(pfc, "N_ncdm", &int2, &flag2, errmsg);
  if (flag1 == _TRUE_)
    N_ncdm_standard = int1;
  else if (flag2 == _TRUE_)
    N_ncdm_standard = int2;

  parser_read_int(pfc, "N_ncdm_decay_dr", &N_ncdm_decay_dr, &flag3, errmsg);

  const bool has_legacy = (flag3 == _TRUE_);
  const bool has_dot    = !dot_instances.empty();
  if (has_legacy && has_dot && has_unconsumed_dot_type_keys(*pfc, dot_instances)) {
    throw std::invalid_argument(
        "cannot mix legacy N_ncdm_decay_dr with dot-syntax "
        "'*.type = ncdm_decay_dr'; use one or the other");
  }

  std::vector<std::string> keys;
  if (has_dot) {
    keys            = synthesise_decay_dr_ncdm_flat_keys(*pfc, dot_instances);
    N_ncdm_decay_dr = static_cast<int>(dot_instances.size());
  }
  else {
    keys.reserve(N_ncdm_decay_dr);
    for (int n = 0; n < N_ncdm_decay_dr; ++n) {
      keys.push_back("ncdm__" + std::to_string(N_ncdm_standard + n + 1));
    }
  }

  if (N_ncdm_decay_dr <= 0)
    return result;

  // Construct species: decay_dr indices are [N_ncdm_standard, N_ncdm_standard + N_ncdm_decay_dr)
  std::vector<std::unique_ptr<DNCDMSpecies>> species_list;
  species_list.reserve(N_ncdm_decay_dr);
  for (int i = 0; i < N_ncdm_decay_dr; i++) {
    int ncdm_index = N_ncdm_standard + i;
    species_list.push_back(std::make_unique<DNCDMSpecies>(pfc, ncdm_index, i, settings, pba, bgm));
  }

  // ── Read DNCDM-specific Omega/deg parameters ──────────────────────────────
  std::vector<double> Omega_dncdmdr_list, omega_dncdmdr_list, deg_list;
  std::vector<double> Omega_ini_dncdm_list, omega_ini_dncdm_list, Neff_ini_dncdm_list;
  int flag_omega1, flag_omega2, flag_deg, flag_Omega_ini, flag_omega_ini, flag_Neff_ini;

  auto readDoubleListLocal = [&](const char* name, std::vector<double>& values, int* found) {
    try {
      *found = pfc->read_list_of_doubles(name, values) ? _TRUE_ : _FALSE_;
    }
    catch (...) {
      *found = _FALSE_;
    }
    return _SUCCESS_;
  };

  readDoubleListLocal("Omega_dncdmdr", Omega_dncdmdr_list, &flag_omega1);
  readDoubleListLocal("omega_dncdmdr", omega_dncdmdr_list, &flag_omega2);
  readDoubleListLocal("deg_ncdm_decay_dr", deg_list, &flag_deg);
  readDoubleListLocal("Omega_ini_dncdm", Omega_ini_dncdm_list, &flag_Omega_ini);
  readDoubleListLocal("omega_ini_dncdm", omega_ini_dncdm_list, &flag_omega_ini);
  readDoubleListLocal("Neff_ini_dncdm", Neff_ini_dncdm_list, &flag_Neff_ini);

  if ((flag_omega1 == _TRUE_) && (flag_omega2 == _TRUE_)) {
    throw std::invalid_argument(
        "In input file, you can only enter one of Omega_dncdmdr or omega_dncdmdr, choose one");
  }
  if ((flag_deg == _TRUE_) && (flag_Omega_ini == _TRUE_)) {
    throw std::invalid_argument(
        "In input file, you can only enter one of deg_ncdm_decay_dr or Omega_ini_dncdm, choose "
        "one");
  }

  const double h = settings.h;

  auto check_list = [&](const char* name, const std::vector<double>& list) {
    if (static_cast<int>(list.size()) < N_ncdm_decay_dr)
      throw std::invalid_argument(
          std::string(name) + " has " + std::to_string(list.size()) +
          " entries but N_ncdm_decay_dr = " + std::to_string(N_ncdm_decay_dr));
  };

  if (flag_omega1 == _TRUE_) {
    check_list("Omega_dncdmdr", Omega_dncdmdr_list);
    for (int n = 0; n < N_ncdm_decay_dr; n++) {
      species_list[n]->SetOmega0(Omega_dncdmdr_list[n], h);
    }
  }
  else if (flag_omega2 == _TRUE_) {
    check_list("omega_dncdmdr", omega_dncdmdr_list);
    for (int n = 0; n < N_ncdm_decay_dr; n++) {
      species_list[n]->SetOmega0(omega_dncdmdr_list[n] / h / h, h);
    }
  }

  if (flag_deg == _TRUE_) {
    check_list("deg_ncdm_decay_dr", deg_list);
    for (int n = 0; n < N_ncdm_decay_dr; n++) {
      species_list[n]->SetDegAndFactor(deg_list[n]);
    }
  }
  else if (flag_Omega_ini == _TRUE_) {
    check_list("Omega_ini_dncdm", Omega_ini_dncdm_list);
    double a_ini = species_list[0]->GetIni(pba ? pba->a_today * 1e-14 : 1e-14,
                                           pba ? pba->a_today : 1.,
                                           settings.tol_ncdm);
    double z_ini = 1.0 / a_ini - 1.0;
    double H0    = pba ? pba->H0 : settings.h * 1.e5 / _c_;
    for (int n = 0; n < N_ncdm_decay_dr; n++) {
      species_list[n]->SetDeg_from_Omega_ini(z_ini, H0, Omega_ini_dncdm_list[n]);
    }
  }
  else if (flag_omega_ini == _TRUE_) {
    check_list("omega_ini_dncdm", omega_ini_dncdm_list);
    double a_ini = species_list[0]->GetIni(pba ? pba->a_today * 1e-14 : 1e-14,
                                           pba ? pba->a_today : 1.,
                                           settings.tol_ncdm);
    double z_ini = 1.0 / a_ini - 1.0;
    double H0    = pba ? pba->H0 : settings.h * 1.e5 / _c_;
    for (int n = 0; n < N_ncdm_decay_dr; n++) {
      species_list[n]->SetDeg_from_Omega_ini(z_ini, H0, omega_ini_dncdm_list[n] / h / h);
    }
  }
  else if (flag_Neff_ini == _TRUE_) {
    check_list("Neff_ini_dncdm", Neff_ini_dncdm_list);
    double a_ini    = species_list[0]->GetIni(pba ? pba->a_today * 1e-14 : 1e-14,
                                              pba ? pba->a_today : 1.,
                                              settings.tol_ncdm);
    double z_ini    = 1.0 / a_ini - 1.0;
    double H0       = pba ? pba->H0 : settings.h * 1.e5 / _c_;
    double Omega0_g = pba ? pba->Omega0_g : 0.;
    for (int n = 0; n < N_ncdm_decay_dr; n++) {
      double Omega_ini = Neff_ini_dncdm_list[n] * 7. / 8. * std::pow(4. / 11., 4. / 3.) * Omega0_g;
      species_list[n]->SetDeg_from_Omega_ini(z_ini, H0, Omega_ini);
    }
  }

  result.reserve(N_ncdm_decay_dr);
  for (int n = 0; n < N_ncdm_decay_dr; ++n) {
    result.push_back(Named{keys[n], std::move(species_list[n])});
  }
  return result;
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
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, 0.);
  for (int i = 0; i < q_size(); i++) {
    snprintf(tmp, 40, "lnf_dncdm[%d][%d]", ncdm_id_, i);
    w.Add(tmp, 0.);
    snprintf(tmp, 40, "dlnfdlnq_dncdm[%d][%d]", ncdm_id_, i);
    w.Add(tmp, 0.);
    snprintf(tmp, 40, "dlnfdlnq_separate_dncdm[%d][%d]", ncdm_id_, i);
    w.Add(tmp, 0.);
  }
}

void DNCDMSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  char tmp[40];
  snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
  w.Add(tmp, pvecback[bg_number_index()]);
  snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
  w.Add(tmp, Rho(pvecback));
  snprintf(tmp, 40, "(.)p_ncdm[%d]", ncdm_id_);
  w.Add(tmp, P(pvecback));
  const int bg_lnf_idx          = bg_lnf_index();
  const int bg_dlnfdlnq_idx     = bg_dlnfdlnq_index();
  const int bg_dlnfdlnq_sep_idx = bg_dlnfdlnq_sep_index();
  for (int i = 0; i < q_size(); i++) {
    snprintf(tmp, 40, "lnf_dncdm[%d][%d]", ncdm_id_, i);
    w.Add(tmp, pvecback[bg_lnf_idx + i]);
    snprintf(tmp, 40, "dlnfdlnq_dncdm[%d][%d]", ncdm_id_, i);
    w.Add(tmp, pvecback[bg_dlnfdlnq_idx + i]);
    snprintf(tmp, 40, "dlnfdlnq_separate_dncdm[%d][%d]", ncdm_id_, i);
    w.Add(tmp, pvecback[bg_dlnfdlnq_sep_idx + i]);
  }
}

// ── Perturbations ──────────────────────────────────────────────────────────

void DNCDMSpecies::RegisterPerturbationIndices(perturb_vector* pv,
                                               const precision* ppr,
                                               int& index_pt,
                                               const perturb_workspace* /*ppw*/,
                                               int /*gauge*/) {
  index_pt_psi0_ = index_pt;

  pv->l_max_ncdm[ncdm_id_]  = ppr->l_max_ncdm;
  pv->q_size_ncdm[ncdm_id_] = q_size();

  pv->index_ncdm_[ncdm_id_].clear();
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq)
    pv->index_ncdm_[ncdm_id_].push_back(index_pt + iq * (pv->l_max_ncdm[ncdm_id_] + 1));
  index_pt += (pv->l_max_ncdm[ncdm_id_] + 1) * pv->q_size_ncdm[ncdm_id_];
}

void DNCDMSpecies::PerturbDerivs(double /*tau*/,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) {
  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
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

  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    const double q    = q_[iq];
    double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + iq];

    const double epsilon        = std::sqrt(q * q + a2 * M_ncdm * M_ncdm);
    const double qk_div_epsilon = k * q / epsilon;
    const int idx               = pv->index_ncdm_.at(ncdm_id_)[iq];
    const int lmax              = pv->l_max_ncdm[ncdm_id_];

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

void DNCDMSpecies::ApplyInitialConditions(double* y, const PerturbIcContext& ctx) {
  perturb_vector* pv = ctx.ppw->pv;
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return;

  const double* pvecback = ctx.ppw->pvecback;
  for (int index_q = 0; index_q < pv->q_size_ncdm[ncdm_id_]; ++index_q) {
    const int idx           = pv->index_ncdm_[ncdm_id_][index_q];
    const double q          = q_[index_q];
    const double epsilon    = std::sqrt(q * q + ctx.a * ctx.a * M_ * M_);
    const double dlnf0_dlnq = pvecback[index_bg_dlnfdlnq_decay_ + index_q];
    const int lmax          = pv->l_max_ncdm[ncdm_id_];

    y[idx + 0] = -0.25 * ctx.delta_ur * dlnf0_dlnq;
    if (lmax >= 1)
      y[idx + 1] = -epsilon / 3. / q / ctx.k * ctx.theta_ur * dlnf0_dlnq;
    if (lmax >= 2)
      y[idx + 2] = -0.5 * ctx.shear_ur * dlnf0_dlnq;
    if (lmax >= 3)
      y[idx + 3] = -0.25 * ctx.l3_ur * dlnf0_dlnq;
  }
}

// ── Integrated observables ──────────────────────────────────────────────────

double DNCDMSpecies::Delta(const perturb_vector* pv,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(a, k, ppw);
  return d;
}

double DNCDMSpecies::Theta(const perturb_vector* pv,
                           const double* /*y*/,
                           const double* /*pvecback*/,
                           const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(a, k, ppw);
  return t;
}

double DNCDMSpecies::DeltaP(const perturb_vector* pv,
                            const double* y,
                            const double* pvecback,
                            const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;

  const double a      = ppw->scalar_ctx.a;
  double delta_p_ncdm = 0.0;
  for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
    double w0             = std::exp(pvecback[index_bg_lnf_decay_dr1_ + iq]) * dq_[iq];
    const double q        = q_[iq];
    const double epsilon  = std::sqrt(q * q + std::pow(M_ * a, 2));
    delta_p_ncdm         += q * q * q * q / epsilon * w0 * y[pv->index_ncdm_.at(ncdm_id_)[iq]];
  }
  const double fac = factor_ * std::pow(pba_->a_today / a, 4);
  return delta_p_ncdm * fac / 3.;
}

double DNCDMSpecies::RhoPlusPShear(const perturb_vector* pv,
                                   const double* /*y*/,
                                   const double* pvecback,
                                   const perturb_workspace* ppw) const {
  if (pv->index_ncdm_.at(ncdm_id_).empty())
    return 0.;
  const double a = ppw->scalar_ctx.a;
  const double k = ppw->scalar_ctx.k;
  auto [d, t, s] = RescaledPerturbations(a, k, ppw);
  return (Rho(pvecback) + P(pvecback)) * s;
}

std::tuple<double, double, double> DNCDMSpecies::RescaledPerturbations(
    double a, double k, const perturb_workspace* ppw) const {
  double rho_scaled              = 0.;
  double rho_plus_p_scaled       = 0.;
  double rho_delta_scaled        = 0.;
  double rho_plus_p_theta_scaled = 0.;
  double rho_plus_p_shear_scaled = 0.;

  const double* lnf_array = ppw->pvecback + bg_lnf_index();
  const double lnN        = GetRescalingFactor(lnf_array);

  for (int index_q = 0; index_q < q_size(); index_q++) {
    const int index_pt   = ppw->pv->index_ncdm_[ncdm_id()][index_q];
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
