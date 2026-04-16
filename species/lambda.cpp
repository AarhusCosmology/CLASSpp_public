#include "lambda.h"

#include "background.h" /* for class_define_index, _TRUE_ */
#include "background_column_writer.h"

LambdaSpecies::LambdaSpecies(const background& pba)
    : BaseSpecies("Lambda", EnergyType::DarkEnergy), Omega0_lambda_(pba.Omega0_lambda),
      H0_(pba.H0) {}

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

double LambdaSpecies::DpDloga(const double* /*pvecback*/) const {
  return 0.;
}

void LambdaSpecies::WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const {
  w.Add("(.)rho_lambda", 0.);
}

void LambdaSpecies::WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const {
  w.Add("(.)rho_lambda", pvecback[index_bg_rho_lambda_]);
}
