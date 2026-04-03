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
  void SetBackgroundModule(const BackgroundModule* bgm) override { bgm_ = bgm; }
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

  void ComputeBackground(double a_rel, const double* pvecback_B,
                          double* pvecback) override;

  void BackgroundDerivs(double tau, const double* y, double* dy,
                         const double* pvecback) override;

  double Rho(const double* pvecback) const override { return pvecback[index_bg_rho_]; }
  double P(const double* pvecback) const override   { return pvecback[index_bg_p_]; }
  double DpDloga(const double* /*pvecback*/) const override { return 0.; }

  // ── Perturbations ────────────────────────────────────────────────────────
  void RegisterPerturbationIndices(perturb_vector* pv, int& index_pt,
                                   const perturb_workspace* ppw, int gauge) override;

  void PerturbDerivs(double tau, const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;

  /**
   * @note These are approximations for informational use. ScalarField stress-energy
   *       is not yet dispatched via the species interface in perturb_total_stress_energy;
   *       the full expressions (using phi', dV, metric terms) remain in the original
   *       perturbations module code. Migrate to species dispatch in a future PR.
   */
  double Delta(const perturb_vector* pv, const double* y, const double* pvecback) const override;
  double Theta(const perturb_vector* pv, const double* y, const double* pvecback) const override;
  double DeltaP(const perturb_vector* pv, const double* y, const double* pvecback) const override;
  double RhoPlusPShear(const perturb_vector* /*pv*/, const double* /*y*/, const double* /*pvecback*/) const override { return 0.; }

  int bi_phi_index()       const { return index_bi_phi_scf_; }
  int bi_phi_prime_index() const { return index_bi_phi_prime_scf_; }

private:
  double V_scf(double phi) const;
  double dV_scf(double phi) const;
  double ddV_scf(double phi) const;

  const background& pba_;
  const BackgroundModule* bgm_ = nullptr;
  int index_bg_phi_scf_       = -1;
  int index_bg_phi_prime_scf_ = -1;
  int index_bg_V_scf_         = -1;
  int index_bg_dV_scf_        = -1;
  int index_bg_ddV_scf_       = -1;
  int index_bg_p_prime_scf_   = -1;
  int index_bi_phi_scf_       = -1;
  int index_bi_phi_prime_scf_ = -1;
};
