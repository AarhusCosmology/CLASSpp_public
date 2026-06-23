#pragma once

#include <memory>

#include "base_species.h"
#include "species_build_context.h"

struct background;

/**
 * Cosmological constant (Lambda) species.
 *
 * Background:  rho_lambda = Omega0_lambda * H0^2  (constant),
 *              p_lambda = -rho_lambda  (w = -1),  dp/dloga = 0.
 * Perturbations: none (Lambda is homogeneous).
 */
class LambdaSpecies : public BaseSpecies {
 public:
  explicit LambdaSpecies(const background& pba, double omega0_lambda);

  double GetOmega0() const override {
    return Omega0_lambda_;
  }

  // ── Perturbation Layout ────────────────────────────────────────────────────
  // Lambda has no perturbation slots; the layout is empty.
  struct PerturbLayout : BaseSpecies::PerturbLayout {};

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  // ── Background ─────────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // ── Perturbations (Lambda is homogeneous: everything vanishes) ────────────
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& /*layout*/,
                                        const perturb_vector* /*pv*/,
                                        const double* /*y*/,
                                        const double* pvecback,
                                        const perturb_workspace* /*ppw*/) const override {
    StressEnergyContribution se;
    se.rho = Rho(pvecback);
    se.p   = P(pvecback);
    // All perturbation fields are zero (Lambda is homogeneous).
    return se;
  }

  void PerturbDerivs(const BaseSpecies::PerturbLayout& /*layout*/,
                     double /*tau*/,
                     const double* /*y*/,
                     double* /*dy*/,
                     const perturb_parameters_and_workspace& /*ppaw*/) const override {}

 private:
  double Omega0_lambda_;
  double H0_;

  int index_bg_rho_lambda_ = -1;
};
