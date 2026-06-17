#include "lambda.h"

#include "background.h" /* for class_define_index, _TRUE_ */

LambdaSpecies::LambdaSpecies(const background& pba, double omega0_lambda)
    : BaseSpecies("Lambda", EnergyType::DarkEnergy), Omega0_lambda_(omega0_lambda), H0_(pba.H0) {}

void LambdaSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_lambda_, _TRUE_, index_bg, 1);
  index_bg_rho_ = index_bg_rho_lambda_;
}

void LambdaSpecies::ComputeBackground(double /*a_rel*/,
                                      const double* /*pvecback_B*/,
                                      double* pvecback) {
  pvecback[index_bg_rho_lambda_] = Omega0_lambda_ * H0_ * H0_;
}

double LambdaSpecies::Rho(const double* pvecback) const {
  return pvecback[index_bg_rho_lambda_];
}

double LambdaSpecies::P(const double* pvecback) const {
  return -pvecback[index_bg_rho_lambda_];
}

void LambdaSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_lambda", 0.);
}

void LambdaSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_lambda", pvecback[index_bg_rho_lambda_]);
}

// ── Factory ───────────────────────────────────────────────────────────────────

std::vector<Named> LambdaSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // Lambda only exists as the closure species. ConstructSpecies's Pass 2 sets
  // omega0_closure_override before calling this factory. In all other passes
  // (Pass 1, shooting-guess construction) the override is absent and Lambda is
  // absent.
  if (!ctx.omega0_closure_override.has_value())
    return result;

  const double omega0_lambda = *ctx.omega0_closure_override;
  if (omega0_lambda == 0.)
    return result;

  result.push_back({"Lambda", std::make_unique<LambdaSpecies>(*ctx.pba, omega0_lambda)});
  return result;
}
