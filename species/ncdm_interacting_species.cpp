#include "ncdm_interacting_species.h"

#include <algorithm>
#include <cmath>

#include "perturbations_module.h"

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
std::vector<std::unique_ptr<NCDMInteractingSpecies>> NCDMInteractingSpecies::CreateAll(
    FileContent* pfc,
    const NcdmSettings& settings,
    const background* pba,
    const BackgroundModule* bgm) {
  std::vector<std::unique_ptr<NCDMInteractingSpecies>> result;
  char errmsg[2048];
  int N_ncdm_interacting = 0;

  if (pfc->read_int("N_ncdm_interacting", N_ncdm_interacting) && N_ncdm_interacting > 0) {
    for (int n = 0; n < N_ncdm_interacting; ++n) {
      result.push_back(std::make_unique<NCDMInteractingSpecies>(pfc, n, settings, pba, bgm));
    }
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