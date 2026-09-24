#pragma once

#include "errors.h"
#include "species_collection.h"

struct perturb_workspace;

/**
 * A species that modifies gravity itself: the one place CLASS otherwise assumes
 * General Relativity.
 *
 * GR is assumed in exactly two places: the Friedmann equations in
 * BackgroundModule::background_functions, and the scalar Einstein equations in
 * PerturbationsModule::perturb_einstein. A species that changes them (the
 * Jordan-frame scalar of a scalar-tensor theory, say) implements this interface
 * next to BaseSpecies. Each module resolves the role once, as it resolves the PPF
 * fluid, and calls these hooks exactly where it would otherwise use GR. With no
 * such species the modules run their GR lines unchanged. At most one species may
 * modify gravity.
 *
 * Design: docs/superpowers/specs/2026-09-24-scalar-tensor-gravity-design.md
 */
class ModifiedGravity {
 public:
  virtual ~ModifiedGravity()                         = default;
  ModifiedGravity(const ModifiedGravity&)            = delete;
  ModifiedGravity& operator=(const ModifiedGravity&) = delete;
  ModifiedGravity(ModifiedGravity&&)                 = delete;
  ModifiedGravity& operator=(ModifiedGravity&&)      = delete;

  /**
   * Close the Friedmann equations. Called in background_functions after every
   * species has computed its background, with the running totals of all of them;
   * the modified-gravity species itself has contributed nothing yet.
   *
   * On return the species' own density and pressure are the effective ones that
   * make GR's expressions exact, H = sqrt(rho_tot - K/a^2) and
   * H' = -3/2 (rho_tot + p_tot) a + K/a, and they have been added to rho_tot and
   * p_tot. rho_r and rho_m, the radiation- and matter-like parts of rho_tot, are
   * updated to match.
   */
  virtual void CloseFriedmann(double a,
                              double K,
                              const double* pvecback_B,
                              double* pvecback,
                              double& rho_tot,
                              double& p_tot,
                              double& rho_r,
                              double& rho_m) const = 0;

  /* The scalar Einstein equations in synchronous gauge, each replacing the GR line
     of perturb_einstein it stands next to, and called at the same point, so it sees
     the same running totals in ppw (delta_rho, rho_plus_p_theta, delta_p,
     rho_plus_p_shear) and the metric already computed (h' for h'', alpha for
     alpha'). */

  /** h' from the 00 (energy) constraint. */
  virtual double HPrime(double k, const double* y, const perturb_workspace* ppw) const = 0;
  /** eta' from the 0i (momentum) constraint. */
  virtual double EtaPrime(double k, const double* y, const perturb_workspace* ppw) const = 0;
  /** h'' from the trace of the ij equation. */
  virtual double HPrimePrime(double k, const double* y, const perturb_workspace* ppw) const = 0;
  /** alpha' from the traceless part of the ij equation. */
  virtual double AlphaPrime(double k, const double* y, const perturb_workspace* ppw) const = 0;

 protected:
  ModifiedGravity() = default;
};

/** The species that modifies gravity, or nullptr. Each module resolves the role once
 *  through this, so the two agree on what they accept: at most one such species, which
 *  is a statement about the input (severe), not about a point in parameter space. */
inline const ModifiedGravity* FindModifiedGravity(const SpeciesCollection& all_species) {
  const ModifiedGravity* found = nullptr;
  for (const auto& [name, sp] : all_species) {
    if (const auto* gravity = dynamic_cast<const ModifiedGravity*>(sp.get())) {
      class_test_severe(found != nullptr, "at most one species may modify gravity");
      found = gravity;
    }
  }
  return found;
}
