#pragma once
#include "../source/background_column_writer.h"
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

/** Baryons: rho ~ a^{-3}. Two perturbation variables: delta_b, theta_b. */
class BaryonsSpecies : public BaseSpecies {
 public:
  explicit BaryonsSpecies(const background& pba)
      : BaseSpecies("Baryons", EnergyType::Matter), pba_(pba) {}

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_b * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* /*pvecback*/) const override {
    return 0.;
  }
  double DpDloga(const double* /*pvecback*/) const override {
    return 0.;
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_b", 0.);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_b", pvecback[index_bg_rho_]);
  }

  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  void FillSources(const double* y, const double* dy, PerturbSourceContext& ctx) override;
  void ApplyInitialConditions(double* y, const PerturbIcContext& ctx) override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      enum file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  double Delta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return y[pv->index_pt_delta_b];
  }
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return y[pv->index_pt_theta_b];
  }
  /** Baryon pressure perturbation: rho_b * (delta_p_b / rho_b) from pre-computed context. */
  double DeltaP(const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* ppw) const override;
  double RhoPlusPShear(const perturb_vector* /*pv*/,
                       const double* /*y*/,
                       const double* /*pvecback*/,
                       const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }

 private:
  const background& pba_;
};
