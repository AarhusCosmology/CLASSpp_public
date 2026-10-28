#include "lambda.h"

#include "background.h" /* for class_define_index, true */

LambdaSpecies::LambdaSpecies(const background& pba, double omega0_lambda)
    : BaseSpecies("Lambda", EnergyType::DarkEnergy), Omega0_lambda_(omega0_lambda), H0_(pba.H0) {}

void LambdaSpecies::RegisterBackgroundIndices(int& index_bg) {
  class_define_index(index_bg_rho_lambda_, true, index_bg, 1);
  index_bg_rho_ = index_bg_rho_lambda_;
}

void LambdaSpecies::ComputeBackground(double /*a*/,
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

  // As the closure species (Omega_Lambda not given), Lambda is built in ConstructSpecies's
  // Pass 2 from omega0_closure_override. Given in the input, it is a fixed component built in
  // Pass 1 like any other, and another species closes the budget: the fluid, or the scalar
  // field with Omega_scf < 0. Until 2026-10-27 Lambda existed only as the closure species, so
  // Omega_Lambda = 0.3 with a fluid gave a fluid of 0.69 and no Lambda, without a word.
  const double omega0_lambda = ctx.omega0_closure_override.has_value()
                                   ? *ctx.omega0_closure_override
                                   : ctx.pfc->get_or("Omega_Lambda", 0.);
  if (omega0_lambda == 0.) {
    return result;
  }

  result.push_back({"Lambda", std::make_unique<LambdaSpecies>(*ctx.pba, omega0_lambda)});
  return result;
}
