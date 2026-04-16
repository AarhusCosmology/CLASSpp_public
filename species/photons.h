#pragma once
#include "../source/background_column_writer.h"
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

/** Photons: rho ~ a^{-4}. Boltzmann hierarchy with l_max_g multipoles. */
class PhotonsSpecies : public BaseSpecies {
 public:
  explicit PhotonsSpecies(const background& pba)
      : BaseSpecies("Photons", EnergyType::Radiation), pba_(pba) {}

  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_rho_ = index_bg++;
  }

  void ComputeBackground(double a_rel, const double* /*pvecback_B*/, double* pvecback) override {
    pvecback[index_bg_rho_] = pba_.Omega0_g * pba_.H0 * pba_.H0 / (a_rel * a_rel * a_rel * a_rel);
  }

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_rho_] / 3.;
  }
  double DpDloga(const double* pvecback) const override {
    return -4. / 3. * pvecback[index_bg_rho_];
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", 0.);
  }
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override {
    w.Add("(.)rho_g", pvecback[index_bg_rho_]);
  }

  void RegisterPerturbationIndices(perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;

  void RegisterVectorPerturbationIndices(perturb_vector* pv,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  void RegisterTensorPerturbationIndices(perturb_vector* pv,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;

  void PerturbDerivs(double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  void PerturbVectorDerivs(double tau,
                           const double* y,
                           double* dy,
                           const perturb_parameters_and_workspace& ppaw) override;

  void PerturbTensorDerivs(double tau,
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

  /** RSA/TCA active when pv->index_pt_delta_g == -1 (sentinel set by RegisterPerturbationIndices). */
  double Delta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return (pv->index_pt_delta_g >= 0) ? y[pv->index_pt_delta_g] : 0.;
  }
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* /*pvecback*/,
               const perturb_workspace* /*ppw*/) const override {
    return (pv->index_pt_theta_g >= 0) ? y[pv->index_pt_theta_g] : 0.;
  }
  /** δp_g = δρ_g / 3 = ρ_g * δ_g / 3. Returns 0 when RSA is active. */
  double DeltaP(const perturb_vector* pv,
                const double* y,
                const double* pvecback,
                const perturb_workspace* /*ppw*/) const override {
    return (pv->index_pt_delta_g >= 0) ? pvecback[index_bg_rho_] * y[pv->index_pt_delta_g] / 3.
                                       : 0.;
  }
  /**
   * (rho+p)*shear for photons. Uses the TCA-corrected shear_g stored in
   * ppw->scalar_ctx when the shear perturbation index is not evolved.
   */
  double RhoPlusPShear(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

 private:
  const background& pba_;
};
