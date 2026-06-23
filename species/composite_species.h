#pragma once
#include <memory>
#include <string>
#include <vector>

#include "base_species.h"

/**
 * CompositeSpecies: a BaseSpecies that owns N child species and acts as a
 * single all_species_ entry for a physically coupled sector.
 *
 * Background methods sum over children. The perturbation interface (layout
 * signatures) is implemented by each concrete composite: its PerturbLayout
 * nests one sub-layout per child, and every method forwards the matching
 * sub-layout to the child, e.g. child->Delta(my.child_layout, ...).
 *
 * Conventions concrete composites must follow so the Einstein equations
 * get the right totals:
 *   Delta  = rho-weighted average:      Rho() * Delta() == sum_i(rho_i * delta_i)
 *   Theta  = (rho+p)-weighted average:  (Rho()+P()) * Theta() == sum_i((rho_i+p_i) * theta_i)
 *   delta_p / rho_plus_p_shear = summed inside StressEnergy.
 * PerturbDerivs runs a two-phase dispatch:
 *   1. Each child's PerturbDerivs (free-streaming terms)
 *   2. AddCouplingDerivs (coupling terms — override in concrete subclasses)
 */
class CompositeSpecies : public BaseSpecies {
 public:
  CompositeSpecies(std::string name, EnergyType energy_type)
      : BaseSpecies(std::move(name), energy_type) {}

  // ── Registration ────────────────────────────────────────────────────────
  void RegisterBackgroundIndices(int& index_bg) override;
  void RegisterIntegrationIndices(int& index_bi) override;
  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;
  // ── Omega0 (closure) ────────────────────────────────────────────────────
  /** Sums GetOmega0() over all children. */
  double GetOmega0() const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->GetOmega0();
    return sum;
  }

  /** Sums GetRadiationOmega0() over all children (dark-radiation children
   *  contribute their Omega0; matter children contribute 0). */
  double GetRadiationOmega0() const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->GetRadiationOmega0();
    return sum;
  }

  /** Sums DarkRadiationRhoToday() over all children. */
  double DarkRadiationRhoToday(const double* pvecback_integration) const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->DarkRadiationRhoToday(pvecback_integration);
    return sum;
  }

  /** Sums NeutrinoOmega0() over children (matter NCDM children contribute). */
  double NeutrinoOmega0() const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->NeutrinoOmega0();
    return sum;
  }

  /** Sums NeffContribution() over children. */
  double NeffContribution(double z) const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->NeffContribution(z);
    return sum;
  }

  void PrintNeffInfo() const override {
    for (const auto& c : children_)
      c->PrintNeffInfo();
  }
  void PrintMassInfo() const override {
    for (const auto& c : children_)
      c->PrintMassInfo();
  }

  double TensorMasslessRelativisticRho(const double* pvecback) const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->TensorMasslessRelativisticRho(pvecback);
    return sum;
  }
  void CheckUltraRelativisticAtIc(const double* pvecback, double tol) const override {
    for (const auto& c : children_)
      c->CheckUltraRelativisticAtIc(pvecback, tol);
  }
  bool IsUltraRelativisticAtIc(const double* pvecback, double tol) const override {
    for (const auto& c : children_)
      if (!c->IsUltraRelativisticAtIc(pvecback, tol))
        return false;
    return true;
  }
  void WarnIfTooHeavyForHalofit(double m_ev_threshold) const override {
    for (const auto& c : children_)
      c->WarnIfTooHeavyForHalofit(m_ev_threshold);
  }

  /** Threads a_proposed through all children (e.g. a wrapped DNCDM child
   *  may pull the earliest integration start earlier). */
  double BackgroundAIni(double a_proposed, double tol) const override {
    for (const auto& c : children_)
      a_proposed = c->BackgroundAIni(a_proposed, tol);
    return a_proposed;
  }

  // ── Background ──────────────────────────────────────────────────────────
  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void SetThermodynamicsModule(const ThermodynamicsModule* thm) override;
  void SetPerturbs(const perturbs* ppt) override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;
  double Rho(const double* pvecback) const override;
  double P(const double* pvecback) const override;
  double PPrime(double a,
                double H,
                const double* pvecback_B,
                const double* pvecback) const override;
  double FreestreamingRho(const double* pvecback) const override;

  // ── Matter tally ────────────────────────────────────────────────────────
  bool ClustersAsMatter() const override;
  bool IsColdMatterSpecies() const override;

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
                                 const perturb_parameters_and_workspace& ppaw) const;
};
