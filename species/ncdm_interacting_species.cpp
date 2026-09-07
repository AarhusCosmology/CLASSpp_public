#include "ncdm_interacting_species.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "errors.h"
#include "perturbations_module.h"
#include "species/ncdm_family.h"
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
    class_test_severe(pfc.get<std::string>(key).has_value(),
                      "'%s' is no longer supported. Use dot-syntax: '<instance>.<dot-name> = "
                      "...' with '<instance>.type = ncdm_self_interacting'.",
                      key);
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

  class_test_severe(G_eff_value && log10G_eff_value,
                    "species '%s': specify exactly one of G_eff or log10G_eff",
                    instance_name.c_str());
  if (G_eff_value) {
    G_eff_ = *G_eff_value;
  }
  else if (log10G_eff_value) {
    G_eff_ = std::pow(10.0, *log10G_eff_value);
  }

  /* A NEGATIVE coupling is rejected rather than quietly meaning "off". G_eff
     enters the rate squared, so -x and +x are indistinguishable there, while
     every gate in this class reads `G_eff_ <= 0.` as "no interaction" -- a typed
     minus sign would otherwise silently produce standard neutrinos and an
     explicit evolver instead of the model asked for. class_test (not the severe
     variant) because this is a physics value a sampler may vary: the point is
     rejected and the chain survives. */
  class_test(G_eff_ < 0.,
             "species '%s': G_eff must be positive (got %g); it enters the collision rate "
             "squared, so a negative value cannot mean anything. Omit it to switch the "
             "interaction off.",
             instance_name.c_str(),
             G_eff_);

  std::string use_alpha_correction_string = input.get_or<std::string>("use_alpha_correction",
                                                                      "false");

  if (use_alpha_correction_string == "True" || use_alpha_correction_string == "true" ||
      use_alpha_correction_string == "yes") {
    use_alpha_correction_ = true;
  }
  else if (use_alpha_correction_string == "False" || use_alpha_correction_string == "false" ||
           use_alpha_correction_string == "no") {
    use_alpha_correction_ = false;
  }
  else {
    class_stop_severe(
        "species '%s': Please specify use_alpha_correction as either True/true/yes or "
        "False/false/no",
        instance_name.c_str());
  }
}

// ── CreateAll factory ───────────────────────────────────────────────────────
std::vector<Named> NCDMInteractingSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyInteractingKeys(*ctx.pfc);
  return CreateAllNcdmInstances<NCDMInteractingSpecies>(ctx);
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
/* The self-interaction enters the hierarchy as a pure RELAXATION: every
   multipole l >= 2 is damped by its own coefficient times a common rate, and
   nothing else. Both the RHS and the exponential evolver's Jacobian DIAGONAL
   are built from the two helpers below rather than from two copies of the same
   algebra, so the diagonal cannot drift out of step with the term it
   linearises -- the failure mode BaseSpecies::PerturbDerivsDiagonal warns
   about, where a wrong sign degrades the step with nothing reporting it. */

/* G_eff_ >= 0 is a class invariant (the constructor rejects a negative value), so
   squaring it here cannot turn an "interaction off" state into a positive rate.
   Callers still test `G_eff_ <= 0.` first, which is what actually switches the
   collision term off; this function assumes the interaction is on. */
double NCDMInteractingSpecies::CollisionRate(double a, double a_prime_over_a) const {
  const double taudot = std::pow(a, -4) *
                        std::pow(std::pow(4. / 11., 1. / 3.) * T_cmb_ * _k_B_, 5) *
                        std::pow(G_eff_ / (1e12 * _eV_ * _eV_), 2) * (2. * _PI_ / _h_P_) / _c_ *
                        _Mpc_over_m_;

  // Cap taudot to avoid stiff differential equations crashing the integrator
  return std::min(taudot, a_prime_over_a * 1e9);
}

double NCDMInteractingSpecies::AlphaFluid() const {
  return use_alpha_correction_ ? kCorrectedAlpha2 : kAlphaRtaCoefficients[0];
}

double NCDMInteractingSpecies::AlphaHierarchy(int l) const {
  if (!use_alpha_correction_)
    return kAlphaRtaCoefficients[std::min(4, l - 2)];

  const double z_correction = std::pow(l + 0.5, -2);
  const double I_l          = FitIntegralOfl(z_correction);
  const double N            = 7.0 * std::pow(_PI_, 4) / 720.0;
  const double corr_l       = (l == 2) ? kExactHierarchyL2Correction : 0.0;
  return N / (6.0 * std::pow(2.0 * _PI_, 3)) * (800.0 - I_l + corr_l);
}

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

  const double taudot = CollisionRate(std::sqrt(ctx.a2), ctx.a_prime_over_a);

  if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on) {
    // --- Fluid Approximation ---
    const int idx  = ncdm_layout.index_per_q[0];
    dy[idx + 2]   -= AlphaFluid() * taudot * y[idx + 2];
  }
  else {
    // --- Exact Boltzmann Hierarchy ---
    const int lmax = ncdm_layout.l_max;
    for (int iq = 0; iq < ncdm_layout.q_size; ++iq) {
      const int idx = ncdm_layout.index_per_q[iq];
      for (int l = 2; l <= lmax; l++)
        dy[idx + l] -= AlphaHierarchy(l) * taudot * y[idx + l];
    }
  }
}

/* d(dy_i)/d(y_i) of the collision term above. The free-streaming part
   NCDMSpecies contributes is deliberately NOT included: its only self-terms are
   the l = lmax truncation and the q-drift, both of order the Hubble rate, so
   carrying them explicitly costs nothing and leaving them out keeps the
   diagonal to the one piece that is actually stiff (taudot reaches 1e7 * aH at
   a = 1e-6 for log10 G_eff = -1.35, and is capped at 1e9 * aH). */
void NCDMInteractingSpecies::PerturbDerivsDiagonal(
    const BaseSpecies::PerturbLayout& layout,
    double /*tau*/,
    const double* /*y*/,
    double* diag,
    const perturb_parameters_and_workspace& ppaw) const {
  if (G_eff_ <= 0.) {
    return;
  }

  const perturb_workspace* ppw    = ppaw.ppw;
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const auto& ncdm_layout         = static_cast<const NCDMBaseSpecies::PerturbLayout&>(layout);

  const double taudot = CollisionRate(std::sqrt(ctx.a2), ctx.a_prime_over_a);

  if (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on) {
    const int idx  = ncdm_layout.index_per_q[0];
    diag[idx + 2] -= AlphaFluid() * taudot;
  }
  else {
    const int lmax = ncdm_layout.l_max;
    for (int iq = 0; iq < ncdm_layout.q_size; ++iq) {
      const int idx = ncdm_layout.index_per_q[iq];
      for (int l = 2; l <= lmax; l++)
        diag[idx + l] -= AlphaHierarchy(l) * taudot;
    }
  }
}
