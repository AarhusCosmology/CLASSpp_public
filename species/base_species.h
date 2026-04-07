#pragma once

#include <string>
#include <memory>

// Forward declarations to avoid circular includes
struct background;
struct precision;
struct perturb_vector;
struct perturb_workspace;
struct perturb_parameters_and_workspace;

#include "perturb_context.h"

class BackgroundModule;  // forward declaration

/**
 * Abstract base class for all cosmological species.
 *
 * Each species:
 *  - claims slots in the background vector (pvecback) and integration vector (y)
 *    by implementing Register*Indices()
 *  - computes its background density/pressure in ComputeBackground()
 *  - optionally provides ODE contributions in BackgroundDerivs()
 *  - provides perturbation equations in PerturbDerivs()
 *  - exposes Delta/Theta/DeltaP/RhoPlusPShear for the Einstein equations
 *
 * The map in BaseModule is: const std::map<std::string, std::unique_ptr<BaseSpecies>>.
 * Use .at("CDM") – never operator[] – to maintain const correctness.
 */
class BaseSpecies {
public:
  /** Classification used by background_functions() to accumulate rho_r, rho_m, etc. */
  enum class EnergyType { Radiation, Matter, DarkEnergy, Other };

  virtual ~BaseSpecies() = default;
  BaseSpecies(const BaseSpecies&) = delete;
  BaseSpecies& operator=(const BaseSpecies&) = delete;

  const std::string& name() const { return name_; }
  EnergyType energy_type() const { return energy_type_; }

  /**
   * Called by BackgroundModule after construction to provide access to its
   * indices (index_bg_a_, index_bg_H_, etc.) and methods (dV_scf, etc.).
   * Species that need it override this; default is no-op.
   */
  virtual void SetBackgroundModule(const BackgroundModule* /*bgm*/) {}


  // ── Background ────────────────────────────────────────────────────────────

  /**
   * Claim consecutive slots in pvecback. Called once during background_indices().
   * Implementation must assign index_bg_rho_ (and any other needed indices) and
   * increment index_bg accordingly.
   */
  virtual void RegisterBackgroundIndices(int& index_bg) = 0;

  /**
   * Index of this species' density in pvecback. Valid after RegisterBackgroundIndices().
   * Returns -1 if species is absent (has_* flag was false).
   */
  int bg_rho_index() const { return index_bg_rho_; }

  /**
   * Index of this species' pressure in pvecback, or -1 if p is not stored separately.
   * Most species compute p analytically from rho; NCDM stores it.
   */
  int bg_p_index() const { return index_bg_p_; }

  /**
   * Claim slots in the ODE integration vector y. Called during background_indices().
   * Default: no integrated variables (analytic species do not need this).
   */
  virtual void RegisterIntegrationIndices(int& index_bi) {}

  /**
   * Compute species background quantities at relative scale factor a_rel = a/a_today.
   * pvecback_B is the current ODE integration vector.
   * Write density (and any other owned quantities) into pvecback.
   */
  virtual void ComputeBackground(double a_rel,
                                  const double* pvecback_B,
                                  double* pvecback) = 0;

  /**
   * Contribute to dy/dtau for species with ODE-integrated background variables.
   * pvecback already contains all evaluated background quantities.
   * Default: nothing (species with analytic rho(a) don't override this).
   */
  virtual void BackgroundDerivs(double tau,
                                  const double* y,
                                  double* dy,
                                  const double* pvecback) {}

  /** Energy density at current background state. */
  virtual double Rho(const double* pvecback) const = 0;

  /** Pressure at current background state. */
  virtual double P(const double* pvecback) const = 0;

  /**
   * d(p)/d(ln a): used by BackgroundModule to compute H' and dp_tot_prime.
   * For radiation (rho ~ a^-4): dp/dloga = -4/3 * rho.
   * For matter (p=0): 0.
   * For Lambda (p=-rho=-const): 0.
   */
  virtual double DpDloga(const double* pvecback) const = 0;

  /**
   * Returns true if this species' PerturbDerivs must run AFTER all other
   * species in a second pass. Used for PPF fluid (FluidSpecies).
   */
  virtual bool RequiresDeferredPerturbDerivs() const { return false; }

  /**
   * Returns true if this species' ComputeBackground must be deferred.
   * Used for FluidSpecies which needs w_fld evaluated before it can run.
   */
  virtual bool RequiresDeferredBackground() const { return false; }

  // ── Perturbations ─────────────────────────────────────────────────────────

  /**
   * Claim consecutive slots in the perturbation integration vector (scalar mode).
   * Also writes the species' indices into the corresponding fields of pv
   * (which other code still uses for e.g. initial conditions and sources).
   * @param ppw   Current workspace; use approx[] flags and curvature info.
   * @param gauge Cast of enum possible_gauges (0=newtonian, 1=synchronous).
   */
  virtual void RegisterPerturbationIndices(perturb_vector* pv, const precision* ppr,
                                            int& index_pt,
                                            const perturb_workspace* ppw,
                                            int gauge) = 0;

  /** Claim perturbation slots for vector modes. Default: no-op (most species are scalar-only). */
  virtual void RegisterVectorPerturbationIndices(perturb_vector* /*pv*/, int& /*index_pt*/,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {}

  /** Claim perturbation slots for tensor modes. Default: no-op (most species are scalar-only). */
  virtual void RegisterTensorPerturbationIndices(perturb_vector* /*pv*/, int& /*index_pt*/,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {}

  /**
   * Contribute to dy for the scalar perturbation ODE at conformal time tau.
   * The PerturbScalarContext inside ppaw->ppw has pre-computed metric terms,
   * cross-species state (delta_g, theta_g, theta_b, ...), and approximation flags.
   */
  virtual void PerturbDerivs(double tau,
                               const double* y,
                               double* dy,
                               const perturb_parameters_and_workspace& ppaw) = 0;

  /** Contribute to dy for the vector perturbation ODE. Default: no-op. */
  virtual void PerturbVectorDerivs(double /*tau*/, const double* /*y*/, double* /*dy*/,
                                    const perturb_parameters_and_workspace& /*ppaw*/) {}

  /** Contribute to dy for the tensor perturbation ODE. Default: no-op. */
  virtual void PerturbTensorDerivs(double /*tau*/, const double* /*y*/, double* /*dy*/,
                                    const perturb_parameters_and_workspace& /*ppaw*/) {}

  /**
   * Fractional density perturbation delta = delta_rho / rho.
   * @param pv       Per-thread perturbation vector; read index_pt_* from here (NOT from species members).
   * @param y        Current ODE state vector (ppw->pv->y).
   * @param pvecback Per-thread background vector (ppw->pvecback).
   * @param ppw      Per-thread workspace; provides scalar_ctx, accumulated stress-energy,
   *                 pvecthermo, and approximation flags for species that need them.
   */
  virtual double Delta(const perturb_vector* pv, const double* y,
                       const double* pvecback, const perturb_workspace* ppw) const = 0;

  /** Velocity divergence theta. */
  virtual double Theta(const perturb_vector* pv, const double* y,
                       const double* pvecback, const perturb_workspace* ppw) const = 0;

  /** Pressure perturbation delta_p. */
  virtual double DeltaP(const perturb_vector* pv, const double* y,
                        const double* pvecback, const perturb_workspace* ppw) const = 0;

  /** (rho + p) * sigma: anisotropic stress contribution to Einstein equations. */
  virtual double RhoPlusPShear(const perturb_vector* pv, const double* y,
                               const double* pvecback, const perturb_workspace* ppw) const = 0;

protected:
  BaseSpecies(std::string name, EnergyType energy_type)
    : name_(std::move(name)), energy_type_(energy_type) {}

  std::string name_;
  EnergyType energy_type_;

  // Set by RegisterBackgroundIndices(); -1 means "not registered / species absent"
  int index_bg_rho_ = -1;
  int index_bg_p_   = -1;  // only set by species that store p separately (e.g. NCDM)
};
