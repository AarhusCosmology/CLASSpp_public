# Remove species-specific physics parameters from pba — design

**Date:** 2026-05-29
**Issue:** #275 (final closing PR)
**Predecessor:** PR #276 (`275-remove-has-flags` — removed `pba->has_<species>` and `pba->Omega0_<species>`)

## Goal

Move every species-specific physics parameter off the `background` struct (`pba`) and onto the corresponding species class.  After this PR, `pba` retains only universal numerical params (H0, T_cmb, Omega0_g, Omega0_b, Omega0_k, h, K, sgnK, a_today, background_method, info/inter flags, verbose).  Each species owns its physics state via private member fields with per-field public getters; modules and the Cython wrapper read those state values exclusively through the species accessors.

## Scope

Twenty species-specific fields, single PR:

| Species | Fields |
|---|---|
| Fluid | `fluid_equation_of_state`, `w0_fld`, `wa_fld`, `cs2_fld`, `Omega_EDE`, `use_ppf`, `c_gamma_over_c_fld` |
| ScalarField | `attractor_ic_scf`, `phi_ini_scf`, `phi_prime_ini_scf`, `scf_parameters`, `scf_tuning_index` |
| DCDM | `Gamma_dcdm`, `Omega_ini_dcdm` |
| IDR | `T_idr`, `l_max_idr` |
| IDM_DRMD | `f_idm_drmd`, `G_over_aH_drmd`, `delta_Neff_drmd`, `z_stop` |

Out of scope:
- `pba->fluid_equation_of_state` enum declaration itself stays in `background.h` (or moves to `fluid.h` — see Open question below).
- The `ClosureSpecies` enum, info/inter flags, verbose, H0/h/T_cmb/Omega0_g/b/k stay on pba (universal numeric state).

## Per-species pattern

Each affected species class gains:

1. **Private member fields**, one per physics param, with the same defaults the field had on pba.
2. **Public const-method getters**, one per field (`w0_fld()`, `wa_fld()`, …) returning the stored value.
3. **An expanded constructor** that takes one positional parameter per physics field (in addition to the existing `(background&, Omega0_<species>)` head).  Verbose, but called only from `CreateAll`, so the cost is contained.
4. **`CreateAll(ctx)`** parses each field from `ctx.pfc` using the same key-name/format logic that lived in `input_module::input_read_parameters`, then constructs the species with the parsed values.

Worked example — Fluid:

```cpp
// fluid.h
class FluidSpecies : public BaseSpecies {
 public:
  FluidSpecies(const background& pba,
               double omega0_fld,
               equation_of_state fluid_eos,
               double w0_fld,
               double wa_fld,
               double cs2_fld,
               double Omega_EDE,
               short  use_ppf,
               double c_gamma_over_c_fld);

  // Accessors — replace pba->w0_fld etc. at all read sites
  equation_of_state fluid_eos() const { return fluid_eos_; }
  double w0_fld() const { return w0_fld_; }
  double wa_fld() const { return wa_fld_; }
  double cs2_fld() const { return cs2_fld_; }
  double Omega_EDE() const { return Omega_EDE_; }
  short  use_ppf() const { return use_ppf_; }
  double c_gamma_over_c_fld() const { return c_gamma_over_c_fld_; }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

 private:
  // existing fields (Omega0_fld_, bgm_ ref, indices, …)
  equation_of_state fluid_eos_ = CLP;
  double w0_fld_              = -1.;
  double wa_fld_              = 0.;
  double cs2_fld_             = 1.;
  double Omega_EDE_           = 0.;
  short  use_ppf_             = _TRUE_;
  double c_gamma_over_c_fld_  = 0.4;
};
```

`CreateAll` does the parse (mirroring the deleted parser block in `input_module.cpp`):

```cpp
std::vector<Named> FluidSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  // Existing Omega_fld / closure-override logic stays.
  double omega0_fld = /* unchanged: closure override OR parser_read_double("Omega_fld") */;
  if (omega0_fld == 0.)
    return {};

  // Parse physics params (logic copied verbatim from input_module).
  equation_of_state fluid_eos = CLP;
  std::string s;
  if (ctx.pfc->read_string("fluid_equation_of_state", s)) {
    if (s.find("EDE") != std::string::npos || s.find("ede") != std::string::npos)
      fluid_eos = EDE;
    /* … */
  }

  double w0_fld = -1., wa_fld = 0., cs2_fld = 1., Omega_EDE = 0.;
  ctx.pfc->read_double("w0_fld", w0_fld);
  ctx.pfc->read_double("wa_fld", wa_fld);
  ctx.pfc->read_double("cs2_fld", cs2_fld);
  if (fluid_eos == EDE) ctx.pfc->read_double("Omega_EDE", Omega_EDE);

  short use_ppf = _TRUE_;
  if (ctx.pfc->read_string("use_ppf", s))
    use_ppf = (s.find("y") != std::string::npos) ? _TRUE_ : _FALSE_;
  double c_gamma_over_c_fld = 0.4;
  if (use_ppf == _TRUE_) ctx.pfc->read_double("c_gamma_over_c_fld", c_gamma_over_c_fld);

  return {{"Fluid", std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld, fluid_eos,
                                                    w0_fld, wa_fld, cs2_fld, Omega_EDE,
                                                    use_ppf, c_gamma_over_c_fld)}};
}
```

`InputModule::input_read_parameters` deletes the matching parser block.

The same shape applies to:

- **ScalarField** — adds `attractor_ic_scf_`, `phi_ini_scf_`, `phi_prime_ini_scf_`, `scf_parameters_` (`std::vector<double>`), `scf_tuning_index_`.  `scf_parameters` parser handling moves over with `readDoubleList`.

- **DCDMSpecies (sub-species)** — owns `Gamma_dcdm_` and `Omega_ini_dcdm_` directly (they're particle properties of DCDM, and `dcdm.cpp` already reads them as `pba_.Gamma_dcdm` / `pba_.Omega_ini_dcdm` from its own `pba` reference; switching to private fields is a one-for-one rename).  The `DCDM_DR_Species` composite's constructor takes the two fields and forwards to `DCDMSpecies`'s constructor; the composite reads them as `dcdm_->Gamma_dcdm()` at the call sites in `dcdm_dr_species.cpp` (DR initial-condition seed, decay coupling, Newtonian re-seed).  `Omega_ini_dcdm` is the existing shooter target; the shooter pipeline already routes through `ctx.omega_budget`, and `DoShooting` writes the resolved unknown back into `fc`, so the next iteration's `CreateAll` reads the updated value and stores it on `DCDMSpecies` directly.

- **IDRSpecies (sub-species)** — owns `T_idr_` and `l_max_idr_` directly (`interacting_species.cpp:140-144` already reads `pba_.l_max_idr` from IDR's own `pba` reference).  `T_idr` is currently derived in `ReadCoupledOmegaBudget` (from `N_idr`/`N_dg`/`xi_idr`); after this PR the derivation produces a local, the `IDM_DR_IDR_Species` composite's constructor takes both `T_idr` and `l_max_idr` and forwards to `IDRSpecies`'s constructor.  Cross-module reads (`thermodynamics_module.cpp` `cidm_dr2_` / `Tidm_dr_` computations, 10 sites) go via `bgm_->all_species_.find("IDM_DR_IDR")` → cast to composite → `.idr().T_idr()`.

- **IDM_DRMD_IDR_DRMD_Species (composite)** — `f_idm_drmd_`, `G_over_aH_drmd_`, `delta_Neff_drmd_`, `z_stop_` live on the composite, since they're interaction-channel params governing the IDM_DRMD ↔ IDR_DRMD relationship rather than properties of either sub-species individually.  `background_module.cpp` reads (Gamma0_drmd_ etc.) migrate to `comp.f_idm_drmd()` / `comp.z_stop()` / etc. via the existing composite cast.

## Coupling: `ReadCoupledOmegaBudget`

`ReadCoupledOmegaBudget` currently reads `pba->T_idr`, `pba->f_idm_drmd`, `pba->delta_Neff_drmd`, `pba->Gamma_dcdm`, `pba->Omega_ini_dcdm` to compute Omega derivations and CDM subtractions.  After this PR those fields are gone, so the budget step re-parses each needed key inline from `pfc`:

```cpp
// inside ReadCoupledOmegaBudget
double T_idr_local = 0.;
double xi_idr;
if (ctx.pfc->read_double("xi_idr", xi_idr))
  T_idr_local = xi_idr * pba->T_cmb;
// (and the N_idr / N_dg variants)

omega_budget_.idr = stat_f_idr * std::pow(T_idr_local / pba->T_cmb, 4.) * pba->Omega0_g;
```

The same `parser_read_double` calls run again in `IDR_Species::CreateAll` (for storage on the species).  One duplicate parse per key — `FileContent::read_double` is a `std::map` lookup; cost is negligible.

This keeps:
- The budget struct (`SpeciesOmegaBudget`) **Omega-only** — it does not grow into a physics-params container.
- Each species's `CreateAll` self-contained — it parses every field it stores, with no dependency on what the budget chose to read.

## Cross-module read migration

Same pattern as PR #276: every existing `pba->X` read outside the owning species migrates to a species lookup, typically `bgm_->all_species_.find("<Key>")` followed by the accessor call.  Touch points:

| File | Count | Notes |
|---|---|---|
| `species/<self>.cpp` | 23 (Fluid), 14 (SCF), 13 (DCDM), 3 (IDR) | Self-reads — become `field_` / `field()` calls. |
| `source/input_module.cpp` | 49 (across species) | Parser blocks deleted entirely. Closure decision still reads pfc flags (unchanged). |
| `source/background_module.cpp` | 4 (Fluid), 36 (SCF), 2 (DCDM), 6 (DRMD) | Mostly verbose budget print + V_scf evaluation + scf IC seeding.  Migrate to species lookup. |
| `source/perturbations_module.cpp` | 5 (Fluid) | Same lookup pattern. |
| `source/thermodynamics_module.cpp` | 10 (IDR — for `a_idm_dr` coupling math) | Same lookup pattern. |
| `source/nonlinear_module.cpp` | 1 (halofit `pba->wa_fld != 0.` check) | Look up Fluid via `all_species_`. |

## Cython wrapper

`classy.pyx` exposes only `pba->T_idr` directly (line 1481, computing `xi_idr`).  Add `BackgroundModule::GetTIdr()` (returns `T_idr` of the IDR sub-species inside `IDM_DR_IDR_Species`, or 0 if absent) — `generate_wrapper.py` auto-picks-up the new method, and classy.pyx calls it via the same `deref(bam).GetTIdr()` pattern PR #276 introduced for `GetOmega0Species`.

No other classy.pyx changes.  No other wrapper-visible accessors need to be added.

## Defaults

Each new species private field carries an inline default matching the value pba currently declares (e.g. `double w0_fld_ = -1.;`).  `InputModule::input_default_params` did not write any of these (we cleared it in PR #276); the default-cosmology path therefore just gets the in-class defaults for free.

## Verification

Identical harness to PR #276:

1. `make -j8 class` — clean build.
2. `pip install .` — Cython regen succeeds; `cclassy.pxd` auto-generated.
3. `TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py` → 84/84 pass.
4. Before merge: `TEST_LEVEL=2 python -m pytest -q -m test_scenario python/test_class.py` → 1260/1260 pass.

No new tests required — the migration is value-preserving and the existing cross-gauge / classyref reference comparisons catch behavioral drift.

## Open questions

- **Where does the `equation_of_state` enum live?**  Currently declared at `source/background.h:21`.  Two options: (a) leave it there (it's a tiny enum and `background.h` is widely included), or (b) move it to `species/fluid.h` next to `FluidSpecies`.  Either is fine; (a) minimises diff churn, (b) tightens encapsulation.  Pick during implementation.

- **`l_max_idr` is a precision-like int, not really a physics param.**  Currently set by `pba->l_max_idr = ppr->l_max_idr` in `input_read_parameters`.  Move it to `IDR_Species` as a private field — keeps the "everything species-specific lives on the species" invariant.

## What this leaves on pba

After this PR, `background` contains only:

```
double H0;                           // derived from h
double Omega0_g;                     // photons
double T_cmb = 2.7255;
double Omega0_b;                     // baryons
double Omega0_k = 0.;                // curvature
double h = 0.67556;
double K = 0.;
int    sgnK = 0;
double a_today = 1.;
enum ClosureSpecies closure_species; // closure dispatch
enum background_evolution_method background_method;
short short_info / normal_info / long_info;
short inter_normal / inter_closeby;
short background_verbose;
```

`pba` becomes the universal numerical/global state; species own everything species-specific.  Issue #275 closes.
