#ifndef NCDM_INTERACTING_SPECIES_H
#define NCDM_INTERACTING_SPECIES_H

#include <string_view>

#include "ncdm_species.h"
#include "species/species_build_context.h"

class NCDMInteractingSpecies : public NCDMSpecies {
 public:
  /** With an interaction switched on this sector is stiff and a plain explicit
   *  method does not merely slow down, it fails to complete: rkdp45 and tsit5
   *  both time out at every tolerance from 3e-2 to 1e-6, because the limit is
   *  the stability bound on alpha_l * taudot, not accuracy. With G_eff <= 0 the
   *  collision term is switched off entirely and this is an ordinary NCDM
   *  species, so it inherits the parent's opt-in.
   *  See BaseSpecies::SupportsExplicitPerturbationEvolver(). */
  bool SupportsExplicitPerturbationEvolver() const override {
    return G_eff_ <= 0. && NCDMSpecies::SupportsExplicitPerturbationEvolver();
  }

  /** The collision term is a PURE DIAGONAL relaxation, so `etd` integrates the
   *  whole of this species' stiffness exactly. Measured against ndf15 at matched
   *  accuracy: 1.8-2.5x faster at log10 G_eff = -1.35, 2.0-2.7x at -3.5.
   *  Conditional on the interaction actually being on -- at G_eff <= 0 there is
   *  no diagonal to integrate and etd degenerates to RK4, which LOSES to
   *  rkdp45/tsit5 at tight tolerance on ordinary content.
   *  See BaseSpecies::PrefersExponentialPerturbationEvolver(). */
  bool PrefersExponentialPerturbationEvolver() const override {
    return G_eff_ > 0.;
  }
  static constexpr std::string_view kTypeName = "ncdm_self_interacting";

  // New input path: parameters read from PFC under <instance_name>.<field>
  NCDMInteractingSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm);

  // Factory method to read N_ncdm_interacting and create instances
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  double FitIntegralOfl(double z) const;

  // Layout-based PerturbDerivs: runs NCDMSpecies hierarchy + collision terms.
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;

  /** Diagonal of the collision term for the exponential (etd) evolver. The
   *  relaxation is exactly diagonal, so this is the whole stiff part. */
  void PerturbDerivsDiagonal(const BaseSpecies::PerturbLayout& layout,
                             double tau,
                             const double* y,
                             double* diag,
                             const perturb_parameters_and_workspace& ppaw) const override;

  /** Conformal-time relaxation rate [1/Mpc], capped exactly as the RHS caps it.
   *  Public so the finite-difference test can pin the diagonal against it. */
  double CollisionRate(double a, double a_prime_over_a) const;
  /** Damping coefficient for multipole l >= 2 of the exact hierarchy. */
  double AlphaHierarchy(int l) const;
  /** Damping coefficient for l = 2 under the ncdm fluid approximation. */
  double AlphaFluid() const;

  double GetGeff() const {
    return G_eff_;
  }

 private:
  double G_eff_              = 0.0;
  bool use_alpha_correction_ = false;
};

#endif  // NCDM_INTERACTING_SPECIES_H