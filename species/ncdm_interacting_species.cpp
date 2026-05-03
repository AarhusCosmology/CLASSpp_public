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

void RejectLegacyInteractingKeys(FileContent& pfc) {
  for (const char* key : kLegacyInteractingKeys) {
    std::string unused;
    if (pfc.read_string(key, unused)) {
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

  double G_eff_value      = 0.;
  double log10G_eff_value = 0.;
  bool has_G              = input.read_double("G_eff", G_eff_value);
  bool has_lG             = input.read_double("log10G_eff", log10G_eff_value);

  if (has_G && has_lG) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of G_eff or log10G_eff");
  }
  if (has_G) {
    G_eff_ = G_eff_value;
  }
  else if (has_lG) {
    G_eff_ = std::pow(10.0, log10G_eff_value);
  }
}

// ── CreateAll factory ───────────────────────────────────────────────────────
std::vector<Named> NCDMInteractingSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyInteractingKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_self_interacting");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);  // mark consumed

    auto sp = std::make_unique<NCDMInteractingSpecies>(ctx.pfc,
                                                       name,
                                                       *ctx.ncdm_settings,
                                                       ctx.pba,
                                                       ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    result.push_back({name, std::move(sp)});
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
