#pragma once

#include <memory>
#include <string>

// Forward declarations to avoid circular includes
struct background;
struct precision;
struct perturb_vector;
struct perturb_workspace;
struct perturb_parameters_and_workspace;

#include "perturb_source_context.h"

class BackgroundModule;        // forward declaration
class BackgroundColumnWriter;  // forward declaration
class ThermodynamicsModule;    // forward declaration
struct perturbs;               // forward declaration

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
 * The collection in BaseModule is a const SpeciesCollection&. Use .at("CDM")
 * (throws if absent) or .find("CDM") (returns pointer-or-nullptr); iterate
 * with `for (auto& sp : all_species_)` and use sp.key / sp->name() as needed.
 */
class BaseSpecies {
 public:
  /** Classification used by background_functions() to accumulate rho_r, rho_m, etc. */
  enum class EnergyType { Radiation, Matter, DarkEnergy, Other };
  enum class TransferColumnSection { all, density, velocity };

  virtual ~BaseSpecies()                     = default;
  BaseSpecies(const BaseSpecies&)            = delete;
  BaseSpecies& operator=(const BaseSpecies&) = delete;

  const std::string& name() const {
    return name_;
  }
  EnergyType energy_type() const {
    return energy_type_;
  }

  /**
   * Called by BackgroundModule after construction to provide access to its
   * indices (index_bg_a_, index_bg_H_, etc.) and methods (dV_scf, etc.).
   * Species that need it override this; default is no-op.
   */
  virtual void SetBackgroundModule(const BackgroundModule* /*bgm*/) {}

  /**
   * Called by PerturbationsModule during construction to provide access to
   * the thermodynamics module (pvecthermo indices, etc.). Default: no-op.
   */
  virtual void SetThermodynamicsModule(const ThermodynamicsModule* /*thm*/) {}

  /**
   * Called by PerturbationsModule during construction to provide access to
   * the perturbs struct (alpha_idm_dr, gauge, etc.). Default: no-op.
   */
  virtual void SetPerturbs(const perturbs* /*ppt*/) {}

  // ── Background ────────────────────────────────────────────────────────────

  /**
   * Claim consecutive slots in pvecback. Called once during background_indices().
   * Implementation must assign index_bg_rho_ (and any other needed indices) and
   * increment index_bg accordingly.
   */
  virtual void RegisterBackgroundIndices(int& index_bg) = 0;

  /**
   * Claim slots in the ODE integration vector y. Called during background_indices().
   * Default: no integrated variables (analytic species do not need this).
   */
  virtual void RegisterIntegrationIndices(int& index_bi) {}

  /**
   * Set initial conditions for ODE-integrated background variables.
   * @param a_rel                Initial relative scale factor.
   * @param pvecback_integration ODE integration vector to be populated.
   */
  virtual void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {}

  /**
   * Compute species background quantities at relative scale factor a_rel = a/a_today.
   * pvecback_B is the current ODE integration vector.
   * Write density (and any other owned quantities) into pvecback.
   */
  virtual void ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) = 0;

  /**
   * Contribute to dy/dtau for species with ODE-integrated background variables.
   * pvecback already contains all evaluated background quantities.
   * Default: nothing (species with analytic rho(a) don't override this).
   */
  virtual void BackgroundDerivs(double tau, const double* y, double* dy, const double* pvecback) {}

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

  // ── Background output ─────────────────────────────────────────────────────

  /**
   * True for species that are free-streaming and massless at IC time
   * (deep radiation domination). Used to accumulate rho_nu / fracnu in
   * perturbation IC setup. Static at construction time.
   * Default: false.
   */
  virtual bool IsFreestreaming() const {
    return false;
  }

  /**
   * Free-streaming radiation density contributed by this species at the current
   * background state. Plain species derive this from IsFreestreaming(); composites
   * override it to sum over their children.
   */
  virtual double FreestreamingRho(const double* pvecback) const {
    return IsFreestreaming() ? Rho(pvecback) : 0.;
  }

  /**
   * Returns true if this species is present (has registered its background
   * indices). For top-level species use all_species_.count(); IsPresent() is
   * for sub-components of composites where all_species_.count() is unavailable.
   */
  bool IsPresent() const {
    return index_bg_rho_ >= 0;
  }

  /**
   * Write this species' background column titles into writer.
   * Called by BackgroundModule::background_output_titles().
   * Default: no-op (species with no background output need not override).
   */
  virtual void WriteBackgroundColumnTitles(BackgroundColumnWriter& /*writer*/) const {}

  /**
   * Write this species' background data values into writer.
   * pvecback is the full background vector at the current time step.
   * Called by BackgroundModule::background_output_data().
   * Default: no-op.
   */
  virtual void WriteBackgroundData(const double* /*pvecback*/,
                                   BackgroundColumnWriter& /*writer*/) const {}

  /**
   * Returns true if this species' PerturbDerivs must run AFTER all other
   * species in a second pass. Used for PPF fluid (FluidSpecies).
   */
  virtual bool RequiresDeferredPerturbDerivs() const {
    return false;
  }

  /**
   * Returns true if this species' ComputeBackground must be deferred.
   * Used for FluidSpecies which needs w_fld evaluated before it can run.
   */
  virtual bool RequiresDeferredBackground() const {
    return false;
  }

  // ── Perturbations ─────────────────────────────────────────────────────────

  /**
   * Claim consecutive slots in the perturbation integration vector (scalar mode).
   * Also writes the species' indices into the corresponding fields of pv
   * (which other code still uses for e.g. initial conditions and sources).
   * @param ppw   Current workspace; use approx[] flags and curvature info.
   * @param gauge Cast of enum possible_gauges (0=newtonian, 1=synchronous).
   */
  virtual void RegisterPerturbationIndices(perturb_vector* pv,
                                           const precision* ppr,
                                           int& index_pt,
                                           const perturb_workspace* ppw,
                                           int gauge) = 0;

  /** Claim perturbation slots for vector modes. Default: no-op (most species are scalar-only). */
  virtual void RegisterVectorPerturbationIndices(perturb_vector* /*pv*/,
                                                 int& /*index_pt*/,
                                                 const perturb_workspace* /*ppw*/,
                                                 int /*gauge*/) {}

  /** Claim perturbation slots for tensor modes. Default: no-op (most species are scalar-only). */
  virtual void RegisterTensorPerturbationIndices(perturb_vector* /*pv*/,
                                                 int& /*index_pt*/,
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
  virtual void PerturbVectorDerivs(double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const perturb_parameters_and_workspace& /*ppaw*/) {}

  /** Contribute to dy for the tensor perturbation ODE. Default: no-op. */
  virtual void PerturbTensorDerivs(double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const perturb_parameters_and_workspace& /*ppaw*/) {}

  /**
   * Fractional density perturbation delta = delta_rho / rho.
   * @param pv       Per-thread perturbation vector; read index_pt_* from here (NOT from species members).
   * @param y        Current ODE state vector (ppw->pv->y).
   * @param pvecback Per-thread background vector (ppw->pvecback).
   * @param ppw      Per-thread workspace; provides scalar_ctx, accumulated stress-energy,
   *                 pvecthermo, and approximation flags for species that need them.
   */
  virtual double Delta(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const = 0;

  /** Velocity divergence theta. */
  virtual double Theta(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const = 0;

  /** Pressure perturbation delta_p. */
  virtual double DeltaP(const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const = 0;

  /** (rho + p) * sigma: anisotropic stress contribution to Einstein equations. */
  virtual double RhoPlusPShear(const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const = 0;

  // ── Stage 1: Output ──────────────────────────────────────────────────────

  /**
   * Append this species' output columns (delta, theta) to writer.
   * Called from perturb_output_titles (title mode) and perturb_output_data (data mode).
   * In class_format, the module may call this separately for density and velocity
   * sections to preserve the historical transfer-column ordering.
   * Use writer.Add(title, mod.index_tp_XXX_, active) — writer looks up tk[tp_index].
   */
  virtual void WriteOutputColumns(
      PerturbColumnWriter& /*writer*/,
      const PerturbationsModule& /*mod*/,
      enum file_format /*fmt*/,
      TransferColumnSection /*section*/ = TransferColumnSection::all) const {}

  /**
   * Append this species' time-evolution variables to writer.
   * Called from perturb_prepare_k_output (title mode, tau=0/y=nullptr/ppw=nullptr)
   * and perturb_print_variables_member (data mode).
   * Guard computation behind: if (!writer.IsTitleMode()) { ... compute ... }
   * Use writer.Add(title, computed_value, active).
   */
  virtual void PrintVariables(PerturbColumnWriter& /*writer*/,
                              double /*tau*/,
                              const double* /*y*/,
                              const PerturbationsModule& /*mod*/,
                              const perturb_workspace* /*ppw*/) const {}

  // ── Stage 2: Source filling ───────────────────────────────────────────────

  /**
   * Write per-k, per-tau source values for this species into the source table.
   * Called in perturb_sources_member() after metric/temperature sources are set.
   * Use ctx.p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, ctx.p_mod->index_tp_XXX_,
   *                               ctx.index_tau, ctx.index_k, value).
   * All addressing (index_md, index_ic, index_k, index_tau) is in ctx.
   */
  virtual void FillSources(const double* /*y*/,
                           const double* /*dy*/,
                           PerturbSourceContext& /*ctx*/) {}

  // ── Stage 3: Initial conditions ───────────────────────────────────────────

  /**
   * Write synchronous-gauge initial conditions for this species into y[].
   * y[] == ppw->pv->y. Use ctx.ppw->pv->index_pt_* for index lookup.
   * The module pre-computes ctx.delta_g_ic, ctx.delta_ur, etc.; species use them.
   * Guard each IC type: if (ctx.index_ic == ctx.p_mod->index_ic_ad_) { ... }
   * Newtonian gauge transformation is handled by the module after this loop.
   */
  virtual void ApplyInitialConditions(double* /*y*/, const PerturbIcContext& /*ctx*/) {}

  // ── Matter tally ──────────────────────────────────────────────────────────

  /**
   * True iff this sector participates in the matter tally (delta_m, theta_m,
   * and the delta_cb / theta_cb passes that cap out in perturb_total_stress_energy).
   * Default: based on EnergyType::Matter. Overrides: IDM_DR returns false
   * (preserves a pre-existing asymmetry in the matter tally — see follow-up issue).
   * Composites return true iff any child is matter.
   */
  virtual bool IsMatterSpecies() const {
    return energy_type_ == EnergyType::Matter;
  }

  /** Rho contribution to the matter tally. Default: Rho() if IsMatterSpecies, else 0. */
  virtual double MatterRho(const double* pvecback) const {
    return IsMatterSpecies() ? Rho(pvecback) : 0.;
  }

  /** Rho * Delta contribution to delta_rho_m. Default: Rho*Delta if IsMatterSpecies, else 0. */
  virtual double MatterRhoDelta(const perturb_vector* pv,
                                const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const {
    return IsMatterSpecies() ? Rho(pvecback) * Delta(pv, y, pvecback, ppw) : 0.;
  }

  /** (Rho+P) * Theta contribution to rho_plus_p_theta_m. */
  virtual double MatterRhoPlusPTheta(const perturb_vector* pv,
                                     const double* y,
                                     const double* pvecback,
                                     const perturb_workspace* ppw) const {
    return IsMatterSpecies() ? (Rho(pvecback) + P(pvecback)) * Theta(pv, y, pvecback, ppw) : 0.;
  }

  /** Rho+P contribution to rho_plus_p_m. */
  virtual double MatterRhoPlusP(const double* pvecback) const {
    return IsMatterSpecies() ? Rho(pvecback) + P(pvecback) : 0.;
  }

  /**
   * True iff this sector participates in the cold-matter tally (delta_cb,
   * theta_cb). Cold matter excludes NCDM/DNCDM (which are "warm" matter) but
   * otherwise mirrors IsMatterSpecies. Default: IsMatterSpecies. Overrides:
   * NCDMBaseSpecies returns false; composites scan children.
   */
  virtual bool IsColdMatterSpecies() const {
    return IsMatterSpecies();
  }

 protected:
  BaseSpecies(std::string name, EnergyType energy_type)
      : name_(std::move(name)), energy_type_(energy_type) {}

  std::string name_;
  EnergyType energy_type_;

  // Set by RegisterBackgroundIndices(); -1 means "not registered / species absent"
  int index_bg_rho_ = -1;
  int index_bg_p_   = -1;  // only set by species that store p separately (e.g. NCDM)
};
