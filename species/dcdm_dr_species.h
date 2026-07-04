#pragma once
#include <string_view>

#include "background.h"
#include "composite_species.h"
#include "dark_radiation_species.h"
#include "dcdm.h"
#include "species_build_context.h"

class BackgroundModule;

/**
 * DCDM_DR_Species: composite for Decaying Cold Dark Matter + its decay radiation.
 *
 * EnergyType::Other so the background loop splits rho correctly:
 *   rho_m += rho - 3p  (= rho_dcdm)
 *   rho_r += 3p        (= rho_dr)
 *
 * BackgroundDerivs override adds the DCDM->DR decay source on top of children.
 * AddCouplingDerivs adds the DR Boltzmann l=0,1 coupling source terms from DCDM.
 */
class DCDM_DR_Species : public CompositeSpecies {
 public:
  static constexpr std::string_view kTypeName = "dcdm_dr";

  // Uses the generic CompositeSpecies::PerturbLayout (one owning sub-layout per
  // child in children_ order): the base-class child loops — including the
  // TallyStressEnergy/DelegateTally hot path — index child_layouts, so every
  // composite layout must carry it (#358). Typed views below.
  enum ChildIndex { kDcdm = 0, kDr = 1 };  // children_ order, set in the ctor
  static const DCDMSpecies::PerturbLayout& dcdm_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const DCDMSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kDcdm]);
  }
  static const DarkRadiationSpecies::PerturbLayout& dr_layout(
      const BaseSpecies::PerturbLayout& my) {
    return static_cast<const DarkRadiationSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kDr]);
  }

  DCDM_DR_Species(const background* pba,
                  const BackgroundModule* bgm,
                  double omega0_dcdmdr,
                  double Gamma_dcdm,
                  double Omega_ini_dcdm);

  void SetBackgroundModule(const BackgroundModule* bgm) override;
  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;

  // Override to add DCDM->DR decay source after children
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  // Registration, derivs (children + AddCouplingDerivs), ICs, stress-energy and
  // approximation-switch copies all use the generic CompositeSpecies child
  // loops. Only composite-specific logic is overridden.
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;

  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;

  void WriteOutputColumns(
      PerturbColumnWriter& writer,
      const PerturbationsModule& mod,
      file_format fmt,
      TransferColumnSection section = TransferColumnSection::all) const override;

  // Typed accessors so callers can capture child indices
  DCDMSpecies& dcdm() {
    return *dcdm_;
  }
  DarkRadiationSpecies& dr() {
    return *dr_sp_;
  }
  const DCDMSpecies& dcdm() const {
    return *dcdm_;
  }
  const DarkRadiationSpecies& dr() const {
    return *dr_sp_;
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  DCDMSpecies* dcdm_           = nullptr;  // non-owning pointer into children_
  DarkRadiationSpecies* dr_sp_ = nullptr;  // non-owning pointer into children_
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;

  ShootingTarget shooting_target_{};  // unknown_param empty => no shooting target
  bool needs_shooting_ = false;       // true iff the direct unknown was absent (we guessed)
};
