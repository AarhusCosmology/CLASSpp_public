#pragma once
#include <string_view>
#include <vector>

#include "background.h"
#include "cdm.h"
#include "composite_species.h"
#include "scalar_field.h"
#include "species_build_context.h"

// ── Type-3 momentum-transfer coupling, as pure free functions (arXiv:1604.04222,
// synchronous gauge). Factored out of AddCouplingDerivs / StressEnergy so the
// physics is unit-testable without a perturb_parameters_and_workspace, and so the
// same closed form feeds both the derivs and the stress-energy paths. Conventions:
// Zbar = -phi_prime_bg/a;  D = 3*rho_cdm*a2 - 2*beta*Zbar^2. Every term carries an
// explicit `beta` factor, so both return exactly 0 at beta = 0.

// phi-KG momentum source (added to dy[scf.idx_phi_prime]).
double Type3CouplingDeltaPhiPrime(double beta, double phi_prime_bg, double theta_cdm);

// CDM-Euler momentum-transfer source (added to dy[cdm.idx_theta]).
double Type3CouplingDeltaThetaCdm(double beta,
                                  double k2,
                                  double a_prime_over_a,
                                  double rho_cdm,
                                  double a2,
                                  double Zbar,
                                  double dV,
                                  double phi,
                                  double phi_prime,
                                  double phi_prime_bg,
                                  double theta_cdm);

/**
 * Type3Species: composite for the Pourtsidou-Tram Type-3 pure-momentum-transfer
 * coupling between cold dark matter and a quintessence scalar field
 * (arXiv:1604.04222). Children handle free-streaming; AddCouplingDerivs (Task 5)
 * adds the synchronous-gauge momentum exchange, and StressEnergy adds the
 * -(2 beta/3) Zbar^2 theta_cdm cross-term. Synchronous gauge only.
 */
class Type3Species : public CompositeSpecies {
 public:
  static constexpr std::string_view kTypeName = "cdm_scf_momentum";

  enum ChildIndex { kCdm = 0, kScf = 1 };  // children_ order, set in the ctor
  static const CDMSpecies::PerturbLayout& cdm_layout(const CompositeSpecies::PerturbLayout& my) {
    return static_cast<const CDMSpecies::PerturbLayout&>(*my.child_layouts[kCdm]);
  }
  static const ScalarFieldSpecies::PerturbLayout& scf_layout(
      const CompositeSpecies::PerturbLayout& my) {
    return static_cast<const ScalarFieldSpecies::PerturbLayout&>(*my.child_layouts[kScf]);
  }

  Type3Species(const background& pba, double omega0_cdm, std::unique_ptr<ScalarFieldSpecies> scf);

  CDMSpecies& cdm() {
    return *cdm_;
  }
  ScalarFieldSpecies& scf() {
    return *scf_;
  }
  const CDMSpecies& cdm() const {
    return *cdm_;
  }
  const ScalarFieldSpecies& scf() const {
    return *scf_;
  }
  double beta() const {
    return scf_->beta();
  }

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  std::vector<ShootingTarget> GetShootingTargets() const override {
    return scf_->GetShootingTargets();
  }
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override {
    scf_->ComputeShootingGuess(ctx, guess, dxdy);
  }
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override {
    return scf_->ComputeShootingResidual(ctx, target);
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // Generic child sum (CompositeSpecies::StressEnergy) PLUS the coupling's
  // -(2 beta/3) Zbar^2 theta_cdm contribution to the total (rho+p)theta. Used
  // un-tallied at the IC delta_tot loop; the hot-path total goes through
  // DelegateTally below, which adds the same cross-term.
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;

  // Hot-path total stress-energy for composites bypasses StressEnergy (it sums
  // each child's StressEnergy via TallyStressEnergy). Override so the coupling
  // cross-term reaches ppw->rho_plus_p_theta during the ODE, not only at ICs.
  void DelegateTally(const BaseSpecies::PerturbLayout& layout,
                     const perturb_vector* pv,
                     const double* y,
                     const double* pvecback,
                     const perturb_workspace* ppw,
                     StressEnergyContribution& total,
                     StressEnergyContribution& total_cold,
                     StressEnergyContribution& total_warm) const override;

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  // The -(2 beta/3) Zbar^2 theta_cdm coupling contribution to (rho+p)theta, shared
  // by StressEnergy (ICs) and DelegateTally (hot path). Returns 0 at beta = 0.
  double CrossTermRhoPlusPTheta(const BaseSpecies::PerturbLayout& layout,
                                const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const;

  CDMSpecies* cdm_         = nullptr;
  ScalarFieldSpecies* scf_ = nullptr;
};
