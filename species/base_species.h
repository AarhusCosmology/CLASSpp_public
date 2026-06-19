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

#include "background_ic_context.h"
#include "common.h"  // class_store_columntitle, class_store_double, _TRUE_, _FALSE_
#include "perturb_source_context.h"

class BackgroundModule;      // forward declaration
class ThermodynamicsModule;  // forward declaration

// ── BackgroundColumnWriter ───────────────────────────────────────────────────
/**
 * Thin helper for background output — analogous to PerturbColumnWriter.
 * Construct in title mode (titles != nullptr) or data mode (dataptr != nullptr).
 * Call Add() once per column in both modes; the writer handles the branching.
 */
class BackgroundColumnWriter {
 public:
  explicit BackgroundColumnWriter(std::string& titles) : titles_(&titles) {}
  BackgroundColumnWriter(double* dataptr, int& storeidx)
      : dataptr_(dataptr), storeidx_(&storeidx) {}

  bool IsTitleMode() const {
    return titles_ != nullptr;
  }

  void Add(const char* title, double value, bool condition = true) {
    if (titles_) {
      class_store_columntitle(*titles_, title, condition ? _TRUE_ : _FALSE_);
    }
    else if (dataptr_) {
      class_store_double(dataptr_, value, condition ? _TRUE_ : _FALSE_, (*storeidx_));
    }
  }
  void Add(const std::string& title, double value, bool condition = true) {
    Add(title.c_str(), value, condition);
  }

 private:
  std::string* titles_ = nullptr;
  double* dataptr_     = nullptr;
  int* storeidx_       = nullptr;
};
struct perturbs;  // forward declaration

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
   * indices (index_bg_a_, index_bg_H_, etc.) and background state.
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
   * @param ctx Background initial-condition context (a_rel, a_ini, rho_rad,
   *            and the integration vector to fill).
   */
  virtual void SetBackgroundInitialConditions(const BackgroundICContext& ctx) {}

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
   * Post-Friedmann. This species' conformal-time pressure derivative p'
   * (a' / a = a*H). pvecback_B is the ODE state (IN); pvecback is fully
   * populated through H. Default 0 (matter / Lambda).
   */
  virtual double PPrime(double a,
                        double H,
                        const double* pvecback_B,
                        const double* pvecback) const {
    return 0.;
  }

  /**
   * Post-Friedmann. Write this species' H-dependent owned output slots.
   * H passed explicitly so the hook need not re-read it. Default no-op.
   */
  virtual void FinalizeBackground(double a, double H, const double* pvecback_B, double* pvecback) {}

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
   * Called once after the full background table is built. Lets a species run
   * table-scope analysis over its own columns. Default no-op.
   */
  virtual void ProcessBackgroundTable(const double* /*background_table*/,
                                      int /*n_rows*/,
                                      int /*row_stride*/,
                                      const double* /*z_table*/) {}

  /**
   * Returns true if this species' PerturbDerivs must run AFTER all other
   * species in a second pass. Used for PPF fluid (FluidSpecies).
   */
  virtual bool RequiresDeferredPerturbDerivs() const {
    return false;
  }

  // ── Perturbations ─────────────────────────────────────────────────────────

  // ── Per-species layout signatures ───────────────────────────────────────

  /**
   * Register this species' transfer-source (index_tp) slots. Set once in
   * perturb_indices_of_perturbs() before the k-loop; mirrors
   * RegisterBackgroundIndices (k-independent, cached in a plain species member).
   * Takes no PerturbLayout — that absence marks it as non-per-k state.
   * Default: no-op (species that emit no transfer functions).
   */
  virtual void RegisterTransferSourceIndices(int& /*index_tp*/,
                                             const SourceRequestContext& /*ctx*/) {}

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
   * The layout is this species' PerturbLayout slot (ppw->pv->species_layouts[i]);
   * composites receive their own layout and pass nested sub-layouts to children.
   */
  virtual void PerturbDerivs(const PerturbLayout& layout,
                             double tau,
                             const double* y,
                             double* dy,
                             const perturb_parameters_and_workspace& ppaw) = 0;

  /** Contribute to dy for the vector perturbation ODE. Default: no-op. */
  virtual void PerturbVectorDerivs(const PerturbLayout& /*layout*/,
                                   double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const perturb_parameters_and_workspace& /*ppaw*/) {}

  /** Contribute to dy for the tensor perturbation ODE. Default: no-op. */
  virtual void PerturbTensorDerivs(const PerturbLayout& /*layout*/,
                                   double /*tau*/,
                                   const double* /*y*/,
                                   double* /*dy*/,
                                   const perturb_parameters_and_workspace& /*ppaw*/) {}

  /** Write this species' tensor-mode output column titles. Default: no-op. */
  virtual void WriteTensorOutputColumnTitles(std::string& /*tensor_titles*/) const {}

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
   * @param layout   This species' PerturbLayout slot (ppw->pv->species_layouts[i]);
   *                 composites receive their own layout and pass nested
   *                 sub-layouts to children.
   * @param pv       Per-thread perturbation vector.
   * @param y        Current ODE state vector (ppw->pv->y).
   * @param pvecback Per-thread background vector (ppw->pvecback).
   * @param ppw      Per-thread workspace; provides scalar_ctx, accumulated stress-energy,
   *                 pvecthermo, and approximation flags for species that need them.
   */
  virtual double Delta(const PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const = 0;

  /** Velocity divergence theta. */
  virtual double Theta(const PerturbLayout& layout,
                       const perturb_vector* pv,
                       const double* y,
                       const double* pvecback,
                       const perturb_workspace* ppw) const = 0;

  /** Pressure perturbation delta_p. */
  virtual double DeltaP(const PerturbLayout& layout,
                        const perturb_vector* pv,
                        const double* y,
                        const double* pvecback,
                        const perturb_workspace* ppw) const = 0;

  /** (rho + p) * sigma: anisotropic stress contribution to Einstein equations. */
  virtual double RhoPlusPShear(const PerturbLayout& layout,
                               const perturb_vector* pv,
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
      file_format /*fmt*/,
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
   * Default: no-op.
   */
  virtual void FillSources(const PerturbLayout& /*layout*/,
                           const double* /*y*/,
                           const double* /*dy*/,
                           PerturbSourceContext& /*ctx*/) {}

  // ── Stage 3: Initial conditions ───────────────────────────────────────────

  /**
   * Write synchronous-gauge initial conditions for this species into y[].
   * y[] == ppw->pv->y. The module pre-computes ctx.delta_g_ic, ctx.delta_ur, etc.
   * Guard each IC type: if (ctx.index_ic == ctx.p_mod->index_ic_ad_) { ... }
   * Newtonian gauge transformation is handled by the module after this loop.
   * Default: no-op.
   */
  virtual void ApplyInitialConditions(const PerturbLayout& /*layout*/,
                                      double* /*y*/,
                                      const PerturbIcContext& /*ctx*/) {}

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

  /** Relativistic (radiation-like) Omega0 contribution at early times.
   *  Default 0; ultra-relativistic and interacting-dark-radiation species override. */
  virtual double GetRadiationOmega0() const {
    return 0.;
  }

  // ── NCDM/DR-family hooks (#308) ───────────────────────────────────────────
  // Neutral defaults let modules loop over all_species_ and call these directly
  // instead of downcasting; composites forward/sum over their children.

  /** Dark-radiation energy density today, read from the integrated background
   *  vector. Default 0; DarkRadiationSpecies returns its own; composites sum. */
  virtual double DarkRadiationRhoToday(const double* /*pvecback_integration*/) const {
    return 0.;
  }

  /** This species' contribution to Omega0 of the neutrino/NCDM sector (fnu).
   *  Default 0; NCDM returns GetOmega0(); composites sum children. */
  virtual double NeutrinoOmega0() const {
    return 0.;
  }

  /** Contribution to N_eff at redshift z. Default 0; NCDM returns GetNeff(z);
   *  composites sum children. */
  virtual double NeffContribution(double /*z*/) const {
    return 0.;
  }

  /** Verbose N_eff line (background_verbose). Default no-op; NCDM prints;
   *  composites forward to children. */
  virtual void PrintNeffInfo() const {}
  /** Verbose mass/omega line (background_verbose). Default no-op; NCDM prints;
   *  composites forward to children. */
  virtual void PrintMassInfo() const {}

  /** Contribution to the tensor-mode "relativistic" density in the massless
   *  approximation. Default 0; NCDM returns 3*P; composites sum. */
  virtual double TensorMasslessRelativisticRho(const double* /*pvecback*/) const {
    return 0.;
  }

  /** Initial-time ultra-relativistic check (throws via class_test on failure).
   *  Default no-op; NCDM checks |w-1/3|; composites forward to children. */
  virtual void CheckUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const {}

  /** Initial-time ultra-relativistic predicate. Default true; NCDM checks
   *  |w-1/3| <= tol; composites AND over children. */
  virtual bool IsUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const {
    return true;
  }

  /** Warn (stdout) if this species is too heavy for Halofit/HMcode. The policy
   *  threshold is passed by the nonlinear module. Default no-op; NCDM compares
   *  its mass; composites forward. */
  virtual void WarnIfTooHeavyForHalofit(double /*m_ev_threshold*/) const {}

  /** Earliest a/a_today this species needs integration to start from.
   *  Default: returns a_proposed unchanged. NCDM species may pull it earlier. */
  virtual double BackgroundAIni(double a_proposed, double /*a_today*/, double /*tol*/) const {
    return a_proposed;
  }

  /** Rho contribution to the matter tally. Default: Rho() if IsMatterSpecies, else 0. */
  virtual double MatterRho(const double* pvecback) const {
    return IsMatterSpecies() ? Rho(pvecback) : 0.;
  }

  /** Rho * Delta contribution to delta_rho_m. Default: Rho*Delta if IsMatterSpecies, else 0.
      Uses the species' own layout slot (pv->species_layouts[collection_index_]); only
      valid for first-class species registered via SpeciesCollection::freeze.
      Composites override to tally matter children via their nested sub-layouts.
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
