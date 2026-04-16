#pragma once

#include "base_species.h"

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
  explicit LambdaSpecies(const background& pba);

  // ── Background ─────────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  /** Lambda has no perturbations; does nothing. */
  void RegisterPerturbationIndices(perturb_vector* /*pv*/,
                                   const precision* /*ppr*/,
                                   int& /*index_pt*/,
                                   const perturb_workspace* /*ppw*/,
                                   int /*gauge*/) override {}
  void PerturbDerivs(double /*tau*/,
                     const double* /*y*/,
                     double* /*dy*/,
                     const perturb_parameters_and_workspace& /*ppaw*/) override {}
  double Delta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }
  double Theta(const perturb_vector* /*pv*/,
               const double* /*y*/,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }
  double DeltaP(const perturb_vector* /*pv*/,
                const double* /*y*/,
                const double* /*pvecback*/,
                const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }
  double RhoPlusPShear(const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

 private:
  double Omega0_lambda_;
  double H0_;

  int index_bg_rho_lambda_ = -1;
};
