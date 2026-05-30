#pragma once

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "shooting_target.h"

// Forward declarations to avoid circular includes
struct background;
struct SpeciesBuildContext;
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

  // ── Perturbation layout (per-pv per-species index storage) ───────────────
  /**
   * Polymorphic layout base. Each concrete species defines its own subclass
   * holding the perturbation slot indices/sizes it needs. Layouts are owned by
   * perturb_vector::species_layouts (parallel to all_species_, lex-key order).
   * Per-thread isolation comes for free because each thread allocates its own pv.
   */
  struct PerturbLayout {
    virtual ~PerturbLayout() = default;
  };

  /**
   * Produce a fresh layout for this species. Called once per pv during
   * perturb_vector_init, before Register*PerturbationIndices.
   * Default: empty base PerturbLayout (suitable for species with no
   * perturbation slots, or as a placeholder during migration).
   */
  virtual std::unique_ptr<PerturbLayout> CreatePerturbLayout() const {
    return std::make_unique<PerturbLayout>();
  }

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

  /** Generic species-param accessor for the Python wrapper layer and
   *  cross-module readers that don't want to downcast.  Each species
   *  may opt in by overriding this and returning the named param;
   *  default is std::nullopt. */
  virtual std::optional<double> GetParam(const std::string& /*name*/) const {
    return std::nullopt;
  }

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

  // ── Per-species layout signatures ───────────────────────────────────────

  virtual void RegisterPerturbationIndices(PerturbLayout& /*layout*/,
                                           perturb_vector* /*pv*/,
                                           const precision* /*ppr*/,
                                           int& /*index_pt*/,
                                           const perturb_workspace* /*ppw*/,
                                           int /*gauge*/) {}

  virtual void RegisterVectorPerturbationIndices(PerturbLayout& /*layout*/,
                                                 perturb_vector* /*pv*/,
                                                 const precision* /*ppr*/,
                                                 int& /*index_pt*/,
                                                 const perturb_workspace* /*ppw*/,
                                                 int /*gauge*/) {}

  virtual void RegisterTensorPerturbationIndices(PerturbLayout& /*layout*/,
                                                 perturb_vector* /*pv*/,
                                                 const precision* /*ppr*/,
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

  virtual void PerturbDerivs(const PerturbLayout& /*layout*/,
                             double tau,
                             const double* y,
                             double* dy,
                             const perturb_parameters_and_workspace& ppaw) {
    PerturbDerivs(tau, y, dy, ppaw);
  }

  /** Contribute to dy for the vector perturbation ODE. Default: no-op. */
  virtual void PerturbVectorDerivs(double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const perturb_parameters_and_workspace& /*ppaw*/) {}

  virtual void PerturbVectorDerivs(const PerturbLayout& /*layout*/,
                                   double tau,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) {
    PerturbVectorDerivs(tau, y, dy, ppaw);
  }

  /** Contribute to dy for the tensor perturbation ODE. Default: no-op. */
  virtual void PerturbTensorDerivs(double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const perturb_parameters_and_workspace& /*ppaw*/) {}

  virtual void PerturbTensorDerivs(const PerturbLayout& /*layout*/,
                                   double tau,
                                   const double* y,
                                   double* dy,
                                   const perturb_parameters_and_workspace& ppaw) {
    PerturbTensorDerivs(tau, y, dy, ppaw);
  }

  /** Write this species' tensor-mode output column titles. Default: no-op. */
  virtual void WriteTensorOutputColumnTitles(char* /*tensor_titles*/) const {}

  /** Contribute to the gravitational-wave source term (anisotropic stress) for tensor mode.
   *  Called after perturb_workspace::gw_source is reset to zero.  Default: no-op.
   *  The layout is this species' PerturbLayout slot (ppw->pv->species_layouts[i]). */
  virtual void ContributeTensorGwSource(const PerturbLayout& /*layout*/,
                                        double /*a*/,
                                        double /*a_today*/,
                                        const double* /*y*/,
                                        perturb_workspace* /*ppw*/) const {}

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

  virtual double Delta(const PerturbLayout& /*layout*/,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const {
    return Delta(pv, y, pvecback, ppw);
  }

  /** Velocity divergence theta. */
  virtual double Theta(const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const = 0;

  virtual double Theta(const PerturbLayout& /*layout*/,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const {
    return Theta(pv, y, pvecback, ppw);
  }

  /** Pressure perturbation delta_p. */
  virtual double DeltaP(const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const = 0;

  virtual double DeltaP(const PerturbLayout& /*layout*/,
                        const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const {
    return DeltaP(pv, y, pvecback, ppw);
  }

  /** (rho + p) * sigma: anisotropic stress contribution to Einstein equations. */
  virtual double RhoPlusPShear(const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const = 0;

  virtual double RhoPlusPShear(const PerturbLayout& /*layout*/,
                               const perturb_vector* pv,
                               const double* y,
                               const double* pvecback,
                               const perturb_workspace* ppw) const {
    return RhoPlusPShear(pv, y, pvecback, ppw);
  }

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

  virtual void FillSources(const PerturbLayout& /*layout*/,
                           const double* y,
                           const double* dy,
                           PerturbSourceContext& ctx) {
    FillSources(y, dy, ctx);
  }

  // ── Stage 3: Initial conditions ───────────────────────────────────────────

  /**
   * Write synchronous-gauge initial conditions for this species into y[].
   * y[] == ppw->pv->y. Use ctx.ppw->pv->index_pt_* for index lookup.
   * The module pre-computes ctx.delta_g_ic, ctx.delta_ur, etc.; species use them.
   * Guard each IC type: if (ctx.index_ic == ctx.p_mod->index_ic_ad_) { ... }
   * Newtonian gauge transformation is handled by the module after this loop.
   */
  virtual void ApplyInitialConditions(double* /*y*/, const PerturbIcContext& /*ctx*/) {}

  virtual void ApplyInitialConditions(const PerturbLayout& /*layout*/,
                                      double* y,
                                      const PerturbIcContext& ctx) {
    ApplyInitialConditions(y, ctx);
  }

  /**
   * Copy perturbation state from one layout to another across an approximation switch.
   * Called when the perturbation vector is reallocated (e.g., NCDM FA collapse from
   * full Boltzmann hierarchy to fluid variables delta/theta/shear).
   * Default: no-op. Most species have no approximation switches.
   */
  virtual void CopyPerturbationsAcrossSwitch(const PerturbLayout& /*old_layout*/,
                                             const PerturbLayout& /*new_layout*/,
                                             const double* /*old_y*/,
                                             double* /*new_y*/,
                                             const PerturbSwitchContext& /*ctx*/) const {}

  // ── Synchronous → Newtonian gauge transformation ──────────────────────────
  /**
   * Log conformal-time derivative of the background density, ρ̇/ρ (ρ̇ ≡ dρ̄/dτ).
   * Single source of truth for the density-contrast gauge shift; mirrors the
   * continuity equation used in BackgroundDerivs. Default: ρ̇/ρ = -3ℋ(Rho+P)/Rho.
   * Species with extra source/sink terms (e.g. decay) override this.
   * @param a_prime_over_a  ℋ = a'/a (base does not store H/a indices).
   */
  virtual double RhoDotOverRho(const double* pvecback, double a_prime_over_a) const {
    return -3. * a_prime_over_a * (Rho(pvecback) + P(pvecback)) / Rho(pvecback);
  }

  /**
   * Transform this species' own perturbation variables from synchronous to
   * Newtonian gauge, in place in y[]. Called by the module after alpha is known
   * (ctx.alpha / ctx.alpha_prime filled). The synchronous IC is already in y[]
   * (ApplyInitialConditions runs in both gauges), so this is a pure shift/re-seed.
   * Default: no-op (species with no perturbed variables, e.g. cosmological constant).
   */
  virtual void PerturbSynchronousToNewtonian(const PerturbLayout& /*layout*/,
                                             double* /*y*/,
                                             const PerturbIcContext& /*ctx*/) {}

  // ── Shooter / root-finding ────────────────────────────────────────────────
  // Reported after construction; non-empty iff this species guessed its unknown
  // (target-form input set, direct unknown absent). Drives target detection, the
  // unknown vector, and the lazy-shooting trigger.
  virtual std::vector<ShootingTarget> GetShootingTargets() const {
    return {};
  }
  // Initial guess + Jacobian seed for each unknown, same order as GetShootingTargets.
  // Also called by CreateAll to obtain the value to build from when the unknown is absent.
  virtual void ComputeShootingGuess(const SpeciesBuildContext& /*ctx*/,
                                    std::vector<double>& /*guess*/,
                                    std::vector<double>& /*dxdy*/) const {}
  // computed-today minus the requested target. The authoritative target is supplied by
  // the shooter (collected once at discovery), so the species never re-derives it from
  // the file content — in iteration builds DoShooting has set the unknown, which for some
  // species (e.g. DCDM_DR's overlapping Omega_dcdmdr/Omega_ini_dcdm) makes the user's
  // target indistinguishable by key presence.
  virtual double ComputeShootingResidual(const ShootingResidualContext& /*ctx*/,
                                         const ShootingTarget& /*target*/) const {
    return 0.;
  }

  /**
   * Mark which source slots this species writes during FillSources.
   * Called to build the used_in_sources bitmask/index list.
   * ppw is supplied so species can consult their approximation flags.
   * Default: no-op. Sources-writing species override to flag their tp indices.
   */
  virtual void MarkUsedInSources(const PerturbLayout& /*layout*/,
                                 const perturb_workspace* /*ppw*/,
                                 int* /*used_in_sources*/) const {}

  virtual void MarkVectorUsedInSources(const PerturbLayout& /*layout*/,
                                       const perturb_workspace* /*ppw*/,
                                       int* /*used_in_sources*/) const {}

  virtual void MarkTensorUsedInSources(const PerturbLayout& /*layout*/,
                                       const perturb_workspace* /*ppw*/,
                                       int* /*used_in_sources*/) const {}

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

  /**
   * Returns the species' contribution to Omega0 (used for budget closure during
   * construction). For composites, this defaults to summing children. For
   * decay-product species starting at zero, override returning 0.
   */
  virtual double GetOmega0() const = 0;

  /** Rho contribution to the matter tally. Default: Rho() if IsMatterSpecies, else 0. */
  virtual double MatterRho(const double* pvecback) const {
    return IsMatterSpecies() ? Rho(pvecback) : 0.;
  }

  /** Rho * Delta contribution to delta_rho_m. Default: Rho*Delta if IsMatterSpecies, else 0.
      Prefers the layout-based Delta variant when a species_layouts slot is owned
      (i.e. for first-class species registered via SpeciesCollection::freeze).
      Composite children with no collection_index_ fall back to the legacy variant.
      Defined out-of-line in base_species.cpp because it dereferences
      perturb_vector, which is forward-declared here. */
  virtual double MatterRhoDelta(const perturb_vector* pv,
                                const double* y,
                                const double* pvecback,
                                const perturb_workspace* ppw) const;

  /** (Rho+P) * Theta contribution to rho_plus_p_theta_m. Out-of-line; see above. */
  virtual double MatterRhoPlusPTheta(const perturb_vector* pv,
                                     const double* y,
                                     const double* pvecback,
                                     const perturb_workspace* ppw) const;

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

  /**
   * Universal fluid-like gauge shift: delta += (ρ̇/ρ)·alpha ; theta += k²·alpha
   * (shear, l3 and higher moments are gauge-invariant). Fluid-like species call
   * this from PerturbSynchronousToNewtonian, supplying their own delta/theta slots.
   */
  void ApplyFluidLikeNewtonianShift(double* y,
                                    int idx_delta,
                                    int idx_theta,
                                    const double* pvecback,
                                    const PerturbIcContext& ctx) const {
    if (idx_delta >= 0)
      y[idx_delta] += RhoDotOverRho(pvecback, ctx.a_prime_over_a) * ctx.alpha;
    if (idx_theta >= 0)
      y[idx_theta] += ctx.k * ctx.k * ctx.alpha;
  }

  std::string name_;
  EnergyType energy_type_;

  // Set by SpeciesCollection::freeze(); used by PrintVariables to look up layout.
  std::size_t collection_index_ = SIZE_MAX;

  // Set by RegisterBackgroundIndices(); -1 means "not registered / species absent"
  int index_bg_rho_ = -1;
  int index_bg_p_   = -1;  // only set by species that store p separately (e.g. NCDM)

  friend class SpeciesCollection;
};
