# Unified species `CreateAll` loop and deferred Omega closure

**Status:** design approved, awaiting implementation plan
**Date:** 2026-04-25
**Repo:** AarhusCosmology/CLASSpp

## Problem

`InputModule::input_read_parameters` currently builds full NCDM / DNCDM / NCDMInteracting species *temporarily* (lines 1087–1105 of `source/input_module.cpp`) just to read off their `Omega0`, `N_ncdm`, and `N_decay_dr`, then throws them away. The same species are then constructed for real in `ConstructSpecies`. Two construction paths for the same objects, one of which exists only to feed mid-parser closure logic.

`ConstructSpecies` itself is a hand-written sequence of `if (pba->has_X) all_species_.insert(...)` blocks — the species roster is wired into module code rather than declared in the species subsystem. Adding a new species requires editing `source/input_module.cpp`.

## Goal

Construct every species exactly once, after parameter reading is complete, by iterating a uniform factory list owned by the species subsystem. Adding a new species should not require touching `source/`.

## Non-goals

- Removing `pba->N_ncdm` / `pba->N_decay_dr` from the ~40 module/species use sites. The "build species once" cleanup is in scope; eliminating these counters is a larger follow-up effort that touches every module.
- Changing physics, output, or any user-visible behaviour. The acceptance bar is bit-identical output on reference inputs.

## Architecture

### Two phases, cleanly separated

**Phase 1 — `input_read_parameters` (parsing only, no species built):**

- Reads every input parameter into `pba` / `ppr` / `ppt` exactly as today, including the IDM_DRMD redistribution at lines 1147–1186 and other derivations.
- Drops the temp-NCDM block 1087–1105 entirely. `N_ncdm`, `N_decay_dr`, and `Omega0_ncdm_tot` are no longer set during parameter reading.
- Keeps the validation tests at lines 1196–1200 and 1256–1260 (they read flags from `pfc`, no species needed).
- Replaces the closure assignment at 1213–1243 with closure *detection only*: sets `pba->closure_species` to one of `{None, Lambda, Fluid, ScalarField}` based on the same flag logic. Does not compute the closure value or maintain `Omega_tot`.
- Drops the `Omega_tot += …` accumulator everywhere — no longer needed.

**Phase 2 — `ConstructSpecies` (uniform construction):**

1. Iterate `kAllSpeciesFactories` (see below). For every entry whose `name` ≠ `pba->closure_species`, call its `create_all`, insert results into `all_species_`, and accumulate `Σ Omega0` from `species->GetOmega0()` for each returned entry.
2. After the loop, if `pba->closure_species != None`, write the closure value to the matching field on `pba` (e.g. `pba->Omega0_lambda = 1 − Omega_k − Σ`) and emit the verbose message.
3. Run the closure species' `create_all` exactly the same way and insert.
4. Set `pba->N_ncdm`, `pba->N_decay_dr`, `pba->Omega0_ncdm_tot` from the actually-built species (sum / count).
5. `all_species_.freeze()`.

The closure species' factory reads `pba` like every other factory — it has no special path. Closure orchestration lives entirely in `ConstructSpecies`.

### Factory list lives in `species/`

A new header `species/all_species.h` is the single place edited when adding a species:

```cpp
#pragma once
#include "species/species_build_context.h"
#include "species/photons.h"
#include "species/baryons.h"
#include "species/cdm.h"
// ... one include per species class

struct SpeciesFactoryEntry {
  std::string_view  name;
  std::vector<Named> (*create_all)(const SpeciesBuildContext&);
};

inline constexpr std::array kAllSpeciesFactories = {{
  {"Photons",            &PhotonsSpecies::CreateAll},
  {"Baryons",            &BaryonsSpecies::CreateAll},
  {"CDM",                &CDMSpecies::CreateAll},
  {"UR",                 &UltraRelativisticSpecies::CreateAll},
  {"DCDM_DR",            &DCDM_DR_Species::CreateAll},
  {"NCDM",               &NCDMSpecies::CreateAll},
  {"DNCDM_DR",           &DNCDM_DR_Species::CreateAll},
  {"NCDMInt",            &NCDMInteractingSpecies::CreateAll},
  {"IDM_DR_IDR",         &IDM_DR_IDR_Species::CreateAll},
  {"IDM_DRMD_IDR_DRMD",  &IDM_DRMD_IDR_DRMD_Species::CreateAll},
  {"Lambda",             &LambdaSpecies::CreateAll},
  {"Fluid",              &FluidSpecies::CreateAll},
  {"ScalarField",        &ScalarFieldSpecies::CreateAll},
}};
```

`source/input_module.cpp`'s `ConstructSpecies` does not name any species — it just iterates this constant. Adding a species: write the class, register one row here. No module code changes.

`DNCDMSpecies::CreateAll` is **not** in this list — it's an internal helper used by `DNCDM_DR_Species::CreateAll` to build the children of a composite. Same as today, just signature-migrated to the new context.

### Shared infrastructure (`species/species_build_context.h`)

```cpp
struct SpeciesBuildContext {
  FileContent*            pfc;
  const background*       pba;
  const precision*        ppr;
  const NcdmSettings&     ncdm_settings;
  const BackgroundModule* bgm;   // nullptr at construction time
};

struct Named {
  std::string                  key;
  std::unique_ptr<BaseSpecies> species;
};
```

The unified `Named` replaces the per-class nested `Named` types in the NCDM family (which currently hold concrete-typed unique_ptrs). Inside the NCDM-family `CreateAll` bodies, construction can still use `unique_ptr<NCDMSpecies>` etc.; conversion to `unique_ptr<BaseSpecies>` happens implicitly at return.

### `background` struct gains one field

```cpp
enum class ClosureSpecies { None, Lambda, Fluid, ScalarField };
ClosureSpecies closure_species;
```

Set once in `input_read_parameters`. The orchestrator in `ConstructSpecies` matches it against factory entry names by string, with a small enum-to-string mapping co-located with the field.

### Per-species `CreateAll` migration

| Species | Today's `ConstructSpecies` line | New `CreateAll` body |
|---|---|---|
| Photons | always inserted | always returns 1 entry |
| Baryons | always inserted | always returns 1 entry |
| CDM | `if (has_cdm)` | guard on `has_cdm`; returns 1 or 0 |
| Lambda | `if (has_lambda)` | guard on `has_lambda`; reads `pba->Omega0_lambda` (orchestrator may have just written it) |
| UR | `if (has_ur)` | guard on `has_ur` |
| Fluid | `if (has_fld)` | guard on `has_fld` |
| DCDM_DR | `if (has_dcdm)` | guard on `has_dcdm`; passes `(ctx.pba, ctx.bgm)` |
| ScalarField | `if (has_scf)` | guard on `has_scf` |
| IDM_DR_IDR | `if (has_idm_dr \|\| has_idr)` | guard on the same disjunction |
| IDM_DRMD_IDR_DRMD | `if (has_idm_drmd \|\| has_idr_drmd)` | guard on the same disjunction |
| NCDM | already has `CreateAll` | signature → `(const SpeciesBuildContext&)`, body unchanged |
| DNCDM_DR | already has `CreateAll` | signature changes; internal call to `DNCDMSpecies::CreateAll` also takes context |
| NCDMInteracting | already has `CreateAll` | signature changes |

The `has_*` guard moves from the `ConstructSpecies` call site into the species' own `CreateAll` body. A factory returns 0 entries when its species isn't present, 1 for the simple cases, and N for the NCDM family.

## Treatment of `N_ncdm`, `N_decay_dr`, `Omega0_ncdm_tot`

These three values are still read after construction by `thermodynamics_module`, `nonlinear_module`, perturbations / background / output modules, and the species composites — so they must remain populated on `pba`. They are simply set in a different place.

**In scope:**

- The temp block at 1087–1105 is deleted. After `CreateAll` runs in `ConstructSpecies`, populate `pba->N_ncdm`, `pba->N_decay_dr`, `pba->Omega0_ncdm_tot` from the actually-built species. This becomes the only assignment site.
- Verify and resolve the use sites *inside* `input_read_parameters` and `input_default_params` that read these fields: the precision test at `input_module.cpp:2843` (depends on `N_ncdm > 0`), the derivations in `input_default_params` at `input_module.cpp:3021` (`has_ncdm` from `Omega0_ncdm_tot`) and `:3030` (`n_dncdm_dr` from `N_decay_dr`), and the use at `:3129`. Likely resolution: move those checks/derivations into `ConstructSpecies` post-loop where the values are available, or — if shooting iterations require lightweight per-iteration access (see below) — provide static helpers.
- Verbose prints whose only purpose was to report a counter we're trying to remove → drop the print rather than keep the counter alive.

**Shooting interaction (must be resolved at implementation time):**

`input_init()` runs recursively during shooting. `input_default_params` (line ~3092) reads `Omega0_ncdm_tot` and `N_decay_dr` to derive `has_*` flags from prior-iteration state. Removing the temp block from `input_read_parameters` may leave shooting iterations without these values. Two acceptable resolutions:

1. Confirm shooting does not actually depend on these values during iteration (only on the converged final state, which goes through the full `ConstructSpecies`).
2. Provide static lightweight helpers `NCDMSpecies::TotalOmega0FromInput(pfc, settings, pba)` and `…::CountFromInput(pfc)` that compute from `pfc` without constructing full objects, and call them from inside the shooting loop only.

Option 2 falls slightly short of "zero per-iteration species computation" but keeps the post-shooting flow clean. The implementation plan must investigate and choose.

**Out of scope (follow-up):**

Removing `pba->N_ncdm` / `pba->N_decay_dr` from the ~40 module / species use sites by rewriting them as loops over `all_species_` with virtual dispatch. This is the same direction as the existing "no species-picking in module code" principle, but it is a multi-PR refactor across `perturbations_module`, `background_module`, `output_module`, `dark_radiation_species`, `dcdm_dr_species`. Outside this spec.

## Verification

This is a structural refactor with no physics change. Acceptance bar is bit-identical (or numerically-identical-to-tolerance) output on reference inputs.

- **Reference-output diff.** Run before/after on at minimum: LCDM, ΛCDM+NCDM (1 species), ΛCDM+multiple NCDM, DCDM_DR, DNCDM_DR, NCDM_interacting, IDM_DR_IDR, IDM_DRMD_IDR_DRMD, fluid-as-closure, scalar-field-as-closure (Omega_scf < 0). Diff background, thermodynamics, and Cl outputs.
- **Shooting tests.** At least one shooting target (e.g. `theta_s`). Verify convergence and that converged parameters match pre-refactor values.
- **Closure correctness.** The verbose print of `Omega_lambda` (or fld / scf) at construction time should match the pre-refactor value exactly for every reference input.
- **Build systems.** Compile cleanly under all three: `Makefile`, `setup.py`, `CLASS.xcodeproj`. The new headers `species/all_species.h` and `species/species_build_context.h` must be registered in each.

## Files touched

**New:**
- `species/all_species.h`
- `species/species_build_context.h`

**Modified (species side):**
- Every species `.h` / `.cpp` to add or migrate a static `CreateAll(const SpeciesBuildContext&)`. Trivial wrapper around the existing constructor for the non-NCDM species; signature change for the NCDM family.

**Modified (source side):**
- `include/background.h` — add `closure_species` enum field.
- `source/input_module.cpp` — `input_read_parameters` loses the temp block, the `Omega_tot` accumulator, and the closure-value computation; gains closure-species detection. `ConstructSpecies` reduces to a generic loop over `kAllSpeciesFactories`.
- `source/input_module.h` — minor signature changes if any.

**Modified (build):**
- `Makefile`, `setup.py`, `CLASS.xcodeproj/project.pbxproj` — register the new headers.
