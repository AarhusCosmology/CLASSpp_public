#include "ncdm_interacting_species.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "perturbations_module.h"
#include "species/species_input.h"

namespace {

constexpr const char* kLegacyInteractingKeys[] = {
    "N_ncdm_interacting",
    "G_eff_ncdm_interacting",
    "log10G_eff_ncdm_interacting",
    "quadrature_strategy_ncdm_interacting",
    "N_momentum_bins_ncdm_interacting",
    "maximum_q_ncdm_interacting",
};

constexpr double kAlphaRtaCoefficients[5]{0.40, 0.43, 0.46, 0.47, 0.48};
constexpr double kCorrectedAlpha2            = (7 * _PI_ * (3949.0 - 400 * _PI_ * _PI_)) / 60.;
constexpr double kExactHierarchyL2Correction = 64.0;

void RejectLegacyInteractingKeys(FileContent& pfc) {
  for (const char* key : kLegacyInteractingKeys) {
    if (pfc.get<std::string>(key).has_value()) {
      throw std::invalid_argument(
          std::string("'") + key + "' is no longer supported. Use dot-syntax: " +
          "'<instance>.<dot-name> = ...' with '<instance>.type = ncdm_self_interacting'.");
    }
  }
}

}  // namespace

// ── Constructor ─────────────────────────────────────────────────────────────
// New input path: read per-instance dot-syntax parameters
NCDMInteractingSpecies::NCDMInteractingSpecies(FileContent* pfc,
                                               const std::string& instance_name,
                                               const NcdmSettings& settings,
                                               const background* pba,
                                               const BackgroundModule* bgm)
    : NCDMSpecies(pfc, instance_name, settings, pba, bgm) {
  SpeciesInput input(pfc, instance_name);

  auto G_eff_value      = input.get<double>("G_eff");
  auto log10G_eff_value = input.get<double>("log10G_eff");

  if (G_eff_value && log10G_eff_value) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of G_eff or log10G_eff");
  }
  if (G_eff_value) {
    G_eff_ = *G_eff_value;
  }
  else if (log10G_eff_value) {
    G_eff_ = std::pow(10.0, *log10G_eff_value);
  }

  std::string use_alpha_correction_string = input.get_or<std::string>("use_alpha_correction",
                                                                      "false");

  if (use_alpha_correction_string == "True" || use_alpha_correction_string == "true" ||
      use_alpha_correction_string == "yes") {
    use_alpha_correction_ = _TRUE_;
  }
  else if (use_alpha_correction_string == "False" || use_alpha_correction_string == "false" ||
           use_alpha_correction_string == "no") {
    use_alpha_correction_ = _FALSE_;
  }
  else {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': Please specify use_alpha_correction as either True/true/yes or False/false/no");
  }
}

// ── CreateAll factory ───────────────────────────────────────────────────────
std::vector<Named> NCDMInteractingSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyInteractingKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_self_interacting");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    (void) ctx.pfc->get<std::string>(name + ".type");  // mark consumed

    auto sp = std::make_unique<NCDMInteractingSpecies>(ctx.pfc,
                                                       name,
                                                       *ctx.ncdm_settings,
                                                       ctx.pba,
                                                       ctx.bgm);
    result.push_back({name, std::move(sp)});
  }
  return result;
}

double NCDMInteractingSpecies::FitIntegralOfl(double z) const {
  const double qs[6]{0.0,
                     265.1039577689589,
                     22154.177207573077,
                     687011.2661763201,
                     8041726.168424849,
                     21122401.83355329};

  const double ps[6]{5529600.0,
                     -656640000.0 + 5529600.0 * qs[1],
                     64159257600.0 - 656640000.0 * qs[1] + 5529600.0 * qs[2],
                     17894987354.00044,
                     -9680376946.240694,
                     1623310849.867793};

  double numerator   = ps[0] + ps[1] * z + ps[2] * std::pow(z, 2) + ps[3] * std::pow(z, 3) +
                       ps[4] * std::pow(z, 4) + ps[5] * std::pow(z, 5);
  double denominator = qs[0] + qs[1] * z + qs[2] * std::pow(z, 2) + qs[3] * std::pow(z, 3) +
                       qs[4] * std::pow(z, 4) + qs[5] * std::pow(z, 5);

  return std::pow(z, 3) * numerator / denominator;
}

// ── Perturbations ──────────────────────────────────────────────────────────
void NCDMInteractingSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                                           double tau,
                                           const double* y,
                                           double* dy,
                                           const perturb_parameters_and_workspace& ppaw) const {
  // 1. Compute standard free-streaming derivatives via layout-based NCDMSpecies
  NCDMSpecies::PerturbDerivs(layout, tau, y, dy, ppaw);

  // If there are no interactions, we are done
  if (G_eff_ <= 0.) {
    return;
  }

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const auto& ncdm_layout         = static_cast<const NCDMBaseSpecies::PerturbLayout&>(layout);

  // 2. Extract necessary cosmological variables
  const double a              = std::sqrt(ctx.a2);
  const double a_prime_over_a = ctx.a_prime_over_a;

  // 3. Compute the collision rate (taudot)
  double taudot = std::pow(a, -4) * std::pow(std::pow(4. / 11., 1. / 3.) * T_cmb_ * _k_B_, 5) *
                  std::pow(G_eff_ / (1e12 * _eV_ * _eV_), 2) * (2. * _PI_ / _h_P_) / _c_ *
                  _Mpc_over_m_;

  // Cap taudot to avoid stiff differential equations crashing the integrator
  taudot = std::min(taudot, a_prime_over_a * 1e9);

  // 4. Apply the collision terms
  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);

  if (fa_on) {
    // --- Fluid Approximation ---
    const int idx         = ncdm_layout.index_per_q[0];
    const double alpha_2  = use_alpha_correction_ ? kCorrectedAlpha2 : kAlphaRtaCoefficients[0];
    dy[idx + 2]          -= alpha_2 * taudot * y[idx + 2];
  }
  else {
    // --- Exact Boltzmann Hierarchy ---
    const int lmax = ncdm_layout.l_max;

    for (int iq = 0; iq < ncdm_layout.q_size; ++iq) {
      const int idx = ncdm_layout.index_per_q[iq];

      for (int l = 2; l <= lmax; l++) {
        if (use_alpha_correction_) {
          double z_correction = std::pow(l + 0.5, -2);
          double I_l          = FitIntegralOfl(z_correction);
          double N            = 7.0 * std::pow(_PI_, 4) / 720.0;
          const double corr_l = (l == 2) ? kExactHierarchyL2Correction : 0.0;
          double alpha        = N / (6.0 * std::pow(2.0 * _PI_, 3)) * (800.0 - I_l + corr_l);

          dy[idx + l] -= alpha * taudot * y[idx + l];
        }
        else {
          int alpha_index  = std::min(4, l - 2);
          double alpha     = kAlphaRtaCoefficients[alpha_index];
          dy[idx + l]     -= alpha * taudot * y[idx + l];
        }
      }
    }
  }
}
