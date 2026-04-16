#pragma once
#include "../species/base_species.h"
#include "background.h"
#include "perturbations.h"

class BackgroundModule;

/**
 * Scalar field: phi and phi' integrated via ODE.
 * Potential V(phi) = exp(-lambda*phi) * ((phi-B)^alpha + A).
 */
class ScalarFieldSpecies : public BaseSpecies {
 public:
  explicit ScalarFieldSpecies(const background& pba)
      : BaseSpecies("ScalarField", EnergyType::Other), pba_(pba) {}

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override {
    bgm_ = bgm;
  }
  void RegisterBackgroundIndices(int& index_bg) override {
    index_bg_phi_scf_       = index_bg++;
    index_bg_phi_prime_scf_ = index_bg++;
    index_bg_V_scf_         = index_bg++;
    index_bg_dV_scf_        = index_bg++;
    index_bg_ddV_scf_       = index_bg++;
    index_bg_rho_           = index_bg++;
    index_bg_p_             = index_bg++;
    index_bg_p_prime_scf_   = index_bg++;
  }

  void RegisterIntegrationIndices(int& index_bi) override {
    index_bi_phi_scf_       = index_bi++;
    index_bi_phi_prime_scf_ = index_bi++;
  }

  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;

  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  double Rho(const double* pvecback) const override {
    return pvecback[index_bg_rho_];
  }
  double P(const double* pvecback) const override {
    return pvecback[index_bg_p_];
  }
  double DpDloga(const double* /*pvecback*/) const override {
    return 0.;
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  // ── Perturbations ────────────────────────────────────────────────────────
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

  /**
   * Gauge-dependent fractional density perturbation delta_rho_scf / rho_scf.
   * In Newtonian gauge, includes the metric perturbation psi computed from
   * the accumulated rho_plus_p_shear available in ppw.
   */
  double Delta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Theta(const perturb_vector* pv,
               const double* y,
               const double* pvecback,
               const perturb_workspace* ppw) const override;
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

  /** Compute p'_scf, write it into pvecback, and return its value.
   *  Called by BackgroundModule after ComputeBackground to update p_tot_prime. */
  double ComputePPrimeAndWrite(double a, double* pvecback) const;

  int bi_phi_index() const {
    return index_bi_phi_scf_;
  }
  int bi_phi_prime_index() const {
    return index_bi_phi_prime_scf_;
  }

 private:
  double V_scf(double phi) const;
  double dV_scf(double phi) const;
  double ddV_scf(double phi) const;

  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
  int index_bg_phi_scf_        = -1;
  int index_bg_phi_prime_scf_  = -1;
  int index_bg_V_scf_          = -1;
  int index_bg_dV_scf_         = -1;
  int index_bg_ddV_scf_        = -1;
  int index_bg_p_prime_scf_    = -1;
  int index_bi_phi_scf_        = -1;
  int index_bi_phi_prime_scf_  = -1;
};
