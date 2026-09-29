#include "dncdm_sector.h"

#include <cmath>

#include "background.h"

namespace {

// With inverse decays the parent is held at chemical equilibrium until it goes
// non-relativistic, so the sector ends with a universal comoving energy: 1.233 times its
// initial one, independent of Gamma above H(a_nr) and of the mass (psd on converged grids,
// m = 0.06 and 0.3 eV; the proxy's 16-bin grid gives 1.2235-1.2254 for Gamma12 >= 1e-4).
// A seed only: the fixed point checks it. Spec (C), seed table.
constexpr double kEquilibriumEnergyGain = 1.233;

// Deep enough in the radiation era that (M a / q)^2 ~ 1e-14: a^4 rho is the massless limit.
constexpr double kRadiationEraA = 1e-10;

}  // namespace

DNCDMSector::DNCDMSector(DNCDMSpecies* parent)
    : CompositeSpecies(parent->name(), BaseSpecies::EnergyType::Other), sector_parent_(parent) {}

bool DNCDMSector::CombinedMode() const {
  return sector_parent_->Omega_dncdmdr_pending().has_value();
}

double DNCDMSector::GetOmega0() const {
  if (const auto& reserve = sector_parent_->Omega_dncdmdr_reserve())
    return *reserve;
  if (CombinedMode())
    return *sector_parent_->Omega_dncdmdr_pending();
  return CompositeSpecies::GetOmega0();
}

std::vector<ShootingTarget> DNCDMSector::GetShootingTargets() const {
  if (CombinedMode()) {
    // target_name = the input key; unknown_param = the fc key DoShooting varies.
    return {{name() + ".Omega_dncdmdr", name() + ".deg", *sector_parent_->Omega_dncdmdr_pending()}};
  }
  // Initial and deg modes: the reserve is a fixed point, so there is no target value.
  return {{name() + ".Omega_dncdmdr_fixedpoint", name() + ".Omega_dncdmdr_reserve", 0.}};
}

void DNCDMSector::ComputeShootingGuess(const SpeciesBuildContext& ctx,
                                       std::vector<double>& guess,
                                       std::vector<double>& dxdy) const {
  const double H0 = ctx.pba->H0;
  if (CombinedMode()) {
    // Below H0 the sector today is barely decayed: the parent's stable density, which deg
    // scales, plus its daughters' initial populations, which deg does not.
    double target = *sector_parent_->Omega_dncdmdr_pending();
    if (sector_parent_->Gamma() < H0)
      target -= DaughtersInitialRadiationOmega0(H0);
    auto [g, d] = sector_parent_->DegGuessFromOmegaToday(ctx, target);
    guess.push_back(g);
    dxdy.push_back(d);
    return;
  }
  // The residual is reserve - Omega_sector(reserve), with d Omega_sector / d reserve ~ 0.
  guess.push_back(ReserveSeed(H0));
  dxdy.push_back(1.);
}

double DNCDMSector::ComputeShootingResidual(const ShootingResidualContext& ctx,
                                            const ShootingTarget& target) const {
  const double H0     = ctx.pba->H0;
  const double sector = Rho(ctx.bg_today) / (H0 * H0);
  if (target.target_name == name() + ".Omega_dncdmdr_fixedpoint")
    return GetOmega0() - sector;
  return sector - target.target_value;
}

double DNCDMSector::StableOmega0(double H0) const {
  return sector_parent_->BackgroundDensityOverH0Sq(1., H0) + DaughtersInitialRadiationOmega0(H0);
}

double DNCDMSector::RadiationOmega0(double H0) const {
  const double a4 = std::pow(kRadiationEraA, 4);
  return a4 * sector_parent_->BackgroundDensityOverH0Sq(kRadiationEraA, H0) +
         DaughtersInitialRadiationOmega0(H0);
}

double DNCDMSector::ReserveSeed(double H0) const {
  // Below H0 the parent has barely decayed by today; above it, it has gone, and the sector is
  // radiation carrying its initial energy (plus the rest mass inverse decays convert).
  if (sector_parent_->Gamma() < H0)
    return StableOmega0(H0);
  return (HasInverseDecays() ? kEquilibriumEnergyGain : 1.) * RadiationOmega0(H0);
}
