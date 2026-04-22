#include "ncdm_interacting_species.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "perturbations_module.h"
#include "species/species_input.h"

namespace {

// Field mapping: dot-syntax key (under the instance) -> legacy CSV key.
// Legacy keys are produced by NCDMBaseSpecies::ReadParameters using suffix
// "_interacting" (e.g. base "m_ncdm" -> "m_ncdm_interacting"); NCDM-interacting
// extras (G_eff / log10G_eff) are read directly in the constructor.
struct FieldMap {
  const char* dot;
  const char* legacy;
  const char* default_value;
};
constexpr FieldMap kInteractingFields[] = {
    {"m", "m_ncdm_interacting", "0"},
    {"T", "T_ncdm_interacting", "0.71611"},
    {"deg", "deg_ncdm_interacting", "1.0"},
    {"ksi", "ksi_ncdm_interacting", "0"},
    {"Omega", "Omega_ncdm_interacting", "0"},
    {"omega", "omega_ncdm_interacting", "0"},
    {"psd_parameters", "ncdm_psd_parameters_interacting", "0"},
    {"psd_filename", "ncdm_psd_filenames_interacting", nullptr},
    {"use_psd_file", "use_ncdm_psd_files_interacting", "0"},
    {"quadrature_strategy", "quadrature_strategy_ncdm_interacting", "0"},
    {"momenta_bins", "N_momentum_bins_ncdm_interacting", "5"},
    {"max_q", "maximum_q_ncdm_interacting", "15"},
    {"G_eff", "G_eff_ncdm_interacting", "0"},
    {"log10G_eff", "log10G_eff_ncdm_interacting", "0"},
};

std::vector<std::string> synthesise_self_interacting_ncdm_flat_keys(
    FileContent& pfc, const std::vector<std::string>& instances) {
  const int n_fields = static_cast<int>(sizeof(kInteractingFields) / sizeof(FieldMap));
  (void) CollectInstanceFieldValues(&pfc, instances, "type");
  std::vector<std::vector<std::string>> values;
  values.reserve(n_fields);
  for (int f = 0; f < n_fields; ++f) {
    values.push_back(CollectInstanceFieldValues(&pfc, instances, kInteractingFields[f].dot));
  }

  pfc.set("N_ncdm_interacting", std::to_string(instances.size()));
  const auto g_eff_values      = CollectInstanceFieldValues(&pfc, instances, "G_eff");
  const auto log10g_eff_values = CollectInstanceFieldValues(&pfc, instances, "log10G_eff");
  if (AnyInstanceFieldValue(g_eff_values) && AnyInstanceFieldValue(log10g_eff_values)) {
    throw std::invalid_argument(
        "cannot mix G_eff and log10G_eff representations across dot-syntax self-interacting "
        "NCDM species");
  }

  for (int f = 0; f < n_fields; ++f) {
    if (!AnyInstanceFieldValue(values[f]))
      continue;

    if (std::string(kInteractingFields[f].dot) == "psd_filename") {
      const auto use_psd_file_values = CollectInstanceFieldValues(&pfc, instances, "use_psd_file");
      const std::string csv          = CsvForPsdFilenames(use_psd_file_values,
                                                          values[f],
                                                          "use_psd_file",
                                                          "psd_filename",
                                                          kInteractingFields[f].legacy);
      if (!csv.empty()) {
        pfc.set(kInteractingFields[f].legacy, csv);
      }
      continue;
    }

    pfc.set(kInteractingFields[f].legacy,
            CsvWithDefaults(values[f], kInteractingFields[f].default_value));
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

// ── Constructor ─────────────────────────────────────────────────────────────
NCDMInteractingSpecies::NCDMInteractingSpecies(FileContent* pfc,
                                               int species_index,
                                               const NcdmSettings& settings,
                                               const background* pba,
                                               const BackgroundModule* bgm)
    : NCDMSpecies(pfc, species_index, settings, pba, bgm, "_interacting") {
  char errmsg[2048];  // Local error message buffer to bypass private member access

  std::vector<double> G_eff_list;
  std::vector<double> log10G_eff_list;

  bool flag_G    = pfc->read_list_of_doubles("G_eff_ncdm_interacting", G_eff_list);
  bool flag_logG = pfc->read_list_of_doubles("log10G_eff_ncdm_interacting", log10G_eff_list);

  class_test(flag_G && flag_logG,
             errmsg,
             "In input file, you cannot enter both log10G_eff_ncdm_interacting and "
             "G_eff_ncdm_interacting, choose one");

  if (flag_G) {
    if (species_index < static_cast<int>(G_eff_list.size())) {
      G_eff_ = G_eff_list[species_index];
    }
  }
  else if (flag_logG) {
    if (species_index < static_cast<int>(log10G_eff_list.size())) {
      G_eff_ = std::pow(10.0, log10G_eff_list[species_index]);
    }
  }
}

// ── CreateAll factory ───────────────────────────────────────────────────────
std::vector<NCDMInteractingSpecies::Named> NCDMInteractingSpecies::CreateAll(
    FileContent* pfc,
    const NcdmSettings& settings,
    const background* pba,
    const BackgroundModule* bgm) {
  std::vector<Named> result;
  const auto dot_instances = pfc->instances_with("type", "ncdm_self_interacting");

  int N_ncdm_interacting = 0;
  bool flag              = pfc->read_int("N_ncdm_interacting", N_ncdm_interacting);
  const bool has_legacy  = flag;
  const bool has_dot     = !dot_instances.empty();
  if (has_legacy && has_dot && has_unconsumed_dot_type_keys(*pfc, dot_instances)) {
    throw std::invalid_argument(
        "cannot mix N_ncdm_interacting with dot-syntax "
        "'*.type = ncdm_self_interacting'; use one or the other");
  }

  std::vector<std::string> keys;
  int N = 0;
  if (has_dot) {
    keys = synthesise_self_interacting_ncdm_flat_keys(*pfc, dot_instances);
    N    = static_cast<int>(dot_instances.size());
  }
  else if (has_legacy) {
    int N_std = 0, N_dr = 0;
    pfc->read_int("N_ncdm_standard", N_std);
    pfc->read_int("N_ncdm_decay_dr", N_dr);
    keys.reserve(N_ncdm_interacting);
    for (int n = 0; n < N_ncdm_interacting; ++n) {
      keys.push_back("ncdm__" + std::to_string(N_std + N_dr + n + 1));
    }
    N = N_ncdm_interacting;
  }

  for (int n = 0; n < N; ++n) {
    auto sp = std::make_unique<NCDMInteractingSpecies>(pfc, n, settings, pba, bgm);
    result.push_back(Named{keys[n], std::move(sp)});
  }
  return result;
}

// ── Perturbations ──────────────────────────────────────────────────────────
void NCDMInteractingSpecies::PerturbDerivs(double tau,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) {
  // 1. Compute standard free-streaming derivatives
  NCDMSpecies::PerturbDerivs(tau, y, dy, ppaw);

  // If there are no interactions, we are done
  if (G_eff_ <= 0.) {
    return;
  }

  const perturb_workspace* ppw    = ppaw.ppw;
  const perturb_vector* pv        = ppw->pv;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;

  // 2. Extract necessary cosmological variables
  const double a              = std::sqrt(ctx.a2);
  const double a_prime_over_a = ctx.a_prime_over_a;

  // 3. Compute the collision rate (taudot)
  double taudot = std::pow(a, -4) * std::pow(std::pow(4. / 11., 1. / 3.) * T_cmb_ * _k_B_, 5) *
                  std::pow(G_eff_ / (1e12 * _eV_ * _eV_), 2) * (2. * _PI_ / _h_P_) / _c_ *
                  _Mpc_over_m_;
  double alpha_RTA[5]{0.40, 0.43, 0.46, 0.47, 0.48};

  // Cap taudot to avoid stiff differential equations crashing the integrator
  taudot = std::min(taudot, a_prime_over_a * 1e9);

  // 4. Apply the collision terms
  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);

  if (fa_on) {
    // --- Fluid Approximation ---
    const int idx  = pv->index_ncdm_.at(ncdm_id_)[0];
    dy[idx + 2]   -= 0.40 * taudot * y[idx + 2];
  }
  else {
    // --- Exact Boltzmann Hierarchy ---
    const int lmax = pv->l_max_ncdm[ncdm_id_];

    for (int iq = 0; iq < pv->q_size_ncdm[ncdm_id_]; ++iq) {
      const int idx = pv->index_ncdm_.at(ncdm_id_)[iq];

      for (int l = 2; l <= lmax; l++) {
        int alpha_index  = std::min(4, l - 2);
        double alpha     = alpha_RTA[alpha_index];
        dy[idx + l]     -= alpha * taudot * y[idx + l];
      }
    }
  }
}
