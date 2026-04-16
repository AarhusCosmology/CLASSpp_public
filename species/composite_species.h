#pragma once
#include <memory>
#include <string>
#include <vector>

#include "base_species.h"

/**
 * CompositeSpecies: a BaseSpecies that owns N child species and acts as a
 * single all_species_ entry for a physically coupled sector.
 *
 * Background methods sum over children.
 * PerturbDerivs runs a two-phase dispatch:
 *   1. Each child's PerturbDerivs (free-streaming terms)
 *   2. AddCouplingDerivs (coupling terms — override in concrete subclasses)
 *
 * Delta returns the rho-weighted average so that
 *   Rho() * Delta() == sum_i(rho_i * delta_i)
 * Theta returns the (rho+p)-weighted average so that
 *   (Rho()+P()) * Theta() == sum_i((rho_i+p_i) * theta_i)
 * DeltaP and RhoPlusPShear return the direct sum over children.
 * These are what the Einstein equations need.
 */
class CompositeSpecies : public BaseSpecies {
 public:
  CompositeSpecies(std::string name, EnergyType energy_type)
      : BaseSpecies(std::move(name), energy_type) {}

  // ── Registration ────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
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

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) override;
  void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double DpDloga(const double* pvecback) const override;
  double FreestreamingRho(const double* pvecback) const override;

  // ── Perturbations ────────────────────────────────────────────────────────
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

  // Weighted averages: see class doc-comment above for weighting details.
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
  double RhoPlusPShear(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const override;

 protected:
  std::vector<std::unique_ptr<BaseSpecies>> children_;

  /**
   * Override in concrete subclasses to add coupling terms to dy after
   * all children have written their free-streaming contributions.
   * Default: no-op.
   */
  virtual void AddCouplingDerivs(double tau,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw);
};
