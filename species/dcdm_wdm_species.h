#pragma once
#include <memory>
#include <string_view>
#include <vector>

#include "background.h"
#include "composite_species.h"
#include "dcdm.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"
#include "wdm_decay_product.h"

class BackgroundModule;

/**
 * DCDM_WDM_Species: composite for cold dark matter decaying to two identical
 * massive daughters (arXiv:2606.14849). Children: DCDMSpecies parent (index
 * kDcdm) + WdmDecayProductSpecies daughter (kWdm).
 *
 * EnergyType::Other so the background loop splits rho correctly:
 *   rho_m += rho - 3p, rho_r += 3p.
 *
 * Background wiring: ComputeBackground writes the daughter's injection columns
 * after the child loop (needs rho_dcdm); BackgroundDerivs adds df_i/dtau = J_i
 * and dg_i/dtau = dJ_i/dlnq. AddCouplingDerivs adds the perturbation injection
 * source J_i * delta_dcdm to the daughter's psi_0 bins.
 *
 * Synchronous gauge only (guard in CreateAll). The l=1 parent-velocity source
 * is omitted: theta_dcdm == 0 identically in synchronous gauge.
 */
class DCDM_WDM_Species : public CompositeSpecies {
 public:
  static constexpr std::string_view kTypeName = "dcdm_wdm";

  // Uses the generic CompositeSpecies::PerturbLayout (one owning sub-layout per
  // child in children_ order) — required by the TallyStressEnergy/DelegateTally
  // hot path (#358). Typed views below.
  enum ChildIndex { kDcdm = 0, kWdm = 1 };
  static const DCDMSpecies::PerturbLayout& dcdm_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const DCDMSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kDcdm]);
  }
  static const NCDMBaseSpecies::PerturbLayout& wdm_layout(const BaseSpecies::PerturbLayout& my) {
    return static_cast<const NCDMBaseSpecies::PerturbLayout&>(
        *static_cast<const CompositeSpecies::PerturbLayout&>(my).child_layouts[kWdm]);
  }

  /** Takes ownership of a pre-built daughter (from CreateAll). */
  DCDM_WDM_Species(std::unique_ptr<WdmDecayProductSpecies> wdm,
                   const background* pba,
                   const BackgroundModule* bgm,
                   double omega0_combined,
                   double Omega_ini_dcdm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  void SetBackgroundModule(const BackgroundModule* bgm) override;

  // Injection wiring (see class doc)
  void ComputeBackground(double a, const double* pvecback_B, double* pvecback) override;
  void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) override;

  void WriteBackgroundColumnTitles(BackgroundColumnWriter& w) const override;
  void WriteBackgroundData(const double* pvecback, BackgroundColumnWriter& w) const override;

  /** Combined sector density today (pinned/shot Omega_dcdmwdm), mirroring
   *  DNCDM_DR_Species::GetOmega0: at closure time the daughter's density is
   *  not yet integrated, so the child sum alone would under-reserve. */
  double GetOmega0() const override;

  // ── Shooter hooks ──────────────────────────────────────────────────────────
  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

  // ── Perturbations ──────────────────────────────────────────────────────────
  // Registration, derivs (children + AddCouplingDerivs), ICs, stress-energy and

  DCDMSpecies& dcdm() {
    return *dcdm_;
  }
  const DCDMSpecies& dcdm() const {
    return *dcdm_;
  }
  WdmDecayProductSpecies& wdm() {
    return *wdm_;
  }
  const WdmDecayProductSpecies& wdm() const {
    return *wdm_;
  }

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  DCDMSpecies* dcdm_           = nullptr;  // non-owning pointers into children_
  WdmDecayProductSpecies* wdm_ = nullptr;
  const background* pba_;
  const BackgroundModule* bgm_ = nullptr;

  ShootingTarget shooting_target_{};
  bool needs_shooting_ = false;

  // Scratch for BackgroundDerivs (background integration is single-threaded).
  mutable std::vector<double> scratch_J_, scratch_dJ_;
};
