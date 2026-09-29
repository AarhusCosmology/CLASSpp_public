#pragma once

#include <string>
#include <vector>

#include "composite_species.h"
#include "dncdm_species.h"
#include "species/shooting_target.h"
#include "species/species_build_context.h"

/**
 * What the three ncdm_decay_dr representations (integrated, psd, proxy) share: a decaying
 * NCDM parent and its decay radiation, and the flatness budget that has to reserve their
 * density today before the background can say what it is.
 *
 * Three modes, one per normalisation key:
 *
 *   combined  (Omega_dncdmdr)        shoot deg until the sector holds the requested density;
 *   initial   (Omega_ini, Neff_ini)  deg follows from the initial abundance, and
 *   deg       (deg, or nothing)      the reserve is a fixed point: reserve = Omega_sector(today).
 *
 * The fixed point's unknown is the internal key `<name>.Omega_dncdmdr_reserve`, written only
 * by the shooter. It feeds Omega_Lambda alone, and Omega_Lambda feeds the sector only through
 * H(a) late on, so the map contracts strongly and its seed decides only how many backgrounds
 * it costs, never the answer.
 * docs/superpowers/specs/2026-09-29-shooting-reuse-and-dncdm-budget-design.md (C)
 */
class DNCDMSector : public CompositeSpecies {
 public:
  /** The reserve once the shooter has written it, the requested density in combined mode.
   *  Only a discovery build, before any shooting, falls back to the children's sum. */
  double GetOmega0() const override;

  std::vector<ShootingTarget> GetShootingTargets() const override;
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override;
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override;

 protected:
  explicit DNCDMSector(DNCDMSpecies* parent);

  /** The daughters' initial populations as free radiation today, a^4 rho / H0^2. 0 when they
   *  start empty. */
  virtual double DaughtersInitialRadiationOmega0(double H0) const = 0;

  /** Whether inverse decays hold the parent at chemical equilibrium (the high-Gamma seed). */
  virtual bool HasInverseDecays() const = 0;

 private:
  /** The sector's density today if nothing decayed: the parent's quadrature plus its
   *  daughters' initial populations as radiation. Exact at Gamma = 0. */
  double StableOmega0(double H0) const;

  /** The whole sector's initial energy, redshifted to today as radiation. */
  double RadiationOmega0(double H0) const;

  /** Seed for the reserve (see the class comment and the spec's seed table). */
  double ReserveSeed(double H0) const;

  bool CombinedMode() const;

  DNCDMSpecies* sector_parent_;
};
