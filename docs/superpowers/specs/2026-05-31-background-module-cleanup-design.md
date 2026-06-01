# Background Module Cleanup: drop old RK solver, full index ownership, shared bisection

## Motivation

`background_module.cpp` still carries several artifacts that predate or sit
awkwardly beside the species-as-plugins refactor:

- A second background integrator (`background_solve()`, the old RK path) selected
  by a runtime switch, even though `background_solve_evolver()` is always the
  right choice. It keeps ~210 lines of dead-weight code alive, plus a `growTable`
  usage, an enum, an input read, and a guard in `dncdm_species.cpp`.
- `BackgroundColumnWriter` lives in its own `.h/.cpp` file pair, inconsistent with
  its analog `PerturbColumnWriter`, which rides along in a shared species header.
- Index registration and initial-condition setup reach into the integration
  vector for specific species (DCDM/Fluid/ScalarField) by name and `static_cast`,
  caching four `index_bi_*_` members that are all provably redundant. This
  violates the project rule that module code should loop + dispatch, never pick
  species by type. `background_initial_conditions()` is the worst offender, and
  `V_scf`/`w_fld` physics is duplicated between the module and the species.
- Manual bisection loops are duplicated across four modules.

This is a single bundled PR with four logically-ordered commits. It is a large
diff but cohesive; it is reviewed commit-by-commit.

## Guiding constraints

- **Verification, not bit-identity.** Verify with ~0.1% tolerance on Cl spectra,
  handling Cl^TE zero-crossings; never blind max-rel-diff. Regenerate any
  baseline only after confirming TT/TE drift is acceptable. (See memory:
  `feedback_no_bit_identical_requirement`, `feedback_vectorization_reduction_drift`.)
- **No species-type picking in modules.** After Part C, `background_module.cpp`
  must not name DCDM / Fluid / ScalarField / IDR sub-species to reach into the
  integration vector or set initial conditions. (See memory:
  `feedback_no_species_picking_in_modules`.)
- Sole authors on this codebase; no cross-PR coordination needed. #278
  (`perturb_vector_init`) is in a different module and does not collide.

---

## Part A — Item 3: delete the old RK solver (first commit)

**Goal:** `background_solve_evolver()` becomes the only integrator.

- Delete `BackgroundModule::background_solve()` (`source/background_module.cpp`
  ~902–1113, ~212 lines) and its `growTable` usage within that function.
- **`growTable` stays.** `thermodynamics_module.cpp` and `tools/growTable.cpp` /
  `include/growTable.h` are untouched.
- Replace the `switch (pba->background_method)` (~621–632) with a direct call to
  `background_solve_evolver()`.
- Remove `enum background_evolution_method { bgevo_rk, bgevo_evolver }` and the
  `background_method` field from `source/background.h` (27, 108). Remove the
  `class_read_int("background_method", …)` in `source/input_module.cpp:730`.
- Collapse `species/dncdm_species.cpp:349`: `if (pba_->background_method ==
  bgevo_evolver)` becomes unconditional (drop the `else` branch's special-casing,
  keeping the evolver path's body).
- Remove the `background_solve()` declaration in `source/background_module.h:132`
  and update the file-header doc comment (~41) that describes it.

**Verification:** build + run the standard test suite; spectra must match the
pre-change baseline within tolerance (the default method was already
`bgevo_evolver`, so output should be unchanged).

---

## Part B — Item 1: relocate `BackgroundColumnWriter` (second commit)

**Goal:** match the `PerturbColumnWriter` pattern; eliminate the standalone file
pair.

- `BackgroundColumnWriter` cannot be a `BackgroundModule` method — species call it
  via `WriteBackgroundColumnTitles(BackgroundColumnWriter&)` /
  `WriteBackgroundData(...)` (declared in `species/base_species.h:191,199`).
- Move the class definition (currently `source/background_column_writer.h`) into
  the species-facing header where the background write-hooks are declared
  (alongside `base_species.h`, mirroring how `PerturbColumnWriter` lives in
  `species/perturb_source_context.h`). Move the `Add` implementation accordingly
  (inline in the header, as it is a trivial branch — or into an existing species
  `.cpp` if a translation unit is preferred).
- Delete `source/background_column_writer.h` and
  `source/background_column_writer.cpp`; update the build (Makefile / setup.py
  source lists) and any `#include "background_column_writer.h"`.
- Pure move, no behavior change.

**Verification:** build; output byte-for-byte unchanged (no logic touched).

---

## Part C — Item 2: full index ownership + IC migration (third commit)

**Goal:** `background_module.cpp` never names a specific dynamical species to
address the integration vector or set initial conditions.

### C1 — single registration loop, drop cached indices
- Replace the special-cased block (`source/background_module.cpp:838–868`) with
  one loop:
  ```cpp
  for (auto& [name, sp] : all_species_)
    sp->RegisterIntegrationIndices(index_bi);   // all land in the {B} block
  ```
  Order within `{B}` is irrelevant; each species owns its own offsets.
- Delete the four members `index_bi_rho_dcdm_`, `index_bi_rho_fld_`,
  `index_bi_phi_scf_`, `index_bi_phi_prime_scf_` (`source/background_module.h:
  181–185`).
- At each former read site, use the existing species accessor at point of use:
  - `dcdm().bi_rho_index()` for the dcdm density (the evolver code at ~1274
    already `dynamic_cast`s to `DCDM_DR_Species` right beside the read).
  - `Fluid` `bi_rho_index()` (already on `fluid.h:224`).
  - `ScalarField` `bi_phi_index()` / `bi_phi_prime_index()` (already on
    `scalar_field.h:202,205`).

### C2 — initial conditions via a context struct
- Introduce `BackgroundICContext` (mirroring `PerturbSourceContext`), holding at
  least: `double a_rel; double rho_rad; double* pvecback_integration;
  BackgroundModule* mod;` (extend only as the migration proves necessary).
- Change the virtual `SetBackgroundInitialConditions` (base + the 5 existing
  overrides: `composite_species`, `dcdm`, `dcdm_dr_species`, `dncdm_species`,
  `dncdm_dr_species`) to take the context.
- Add `Fluid::SetBackgroundInitialConditions` override: compute `integral_fld`
  via the species' own `ComputeWFld` and set its `bi_rho_index()` slot (logic
  currently at `background_module.cpp:1467–1485`).
- Add `ScalarField::SetBackgroundInitialConditions` override: use
  `ctx.rho_rad` and the species' own `V_scf` to set `bi_phi_index()` /
  `bi_phi_prime_index()` slots (logic currently at 1495–1524, including the
  `attractor_ic_scf` branch and the `isfinite` `class_test`).
- The module's `background_initial_conditions` then reduces to: set `a`, compute
  `rho_rad`, build the context, run the single dispatch loop
  (`SetBackgroundInitialConditions(ctx)`), then set the universal `{C}` quantities
  (time, tau, rs, D, D').

### C3 — relocate `V_scf`/`w_fld` physics onto the species
- The scalar-field potential (`V_scf`/`dV_scf`/`ddV_scf`) and the fluid `w(a)`
  belong on `ScalarFieldSpecies` / `FluidSpecies`. `ScalarFieldSpecies::V_scf`
  and `FluidSpecies::ComputeWFld` already exist.
- Remove the module-side `V_scf`/`dV_scf`/`ddV_scf`
  (`source/background_module.cpp:1905+`, decls `background_module.h:23–25`) and
  repoint `background_functions` (the scalar-field EOM hot path) at the species
  versions.
- Keep `BackgroundModule::background_w_fld` only as the thin
  delegating wrapper external callers (HyRec, etc.) still need — it already
  delegates to `FluidSpecies::ComputeWFld` with a no-Fluid fallback.
- **Hot-path caution:** `background_functions` runs in the integration inner
  loop. Confirm the species-method indirection does not regress the NCDM/scf hot
  path (watch for auto-vectorization reduction drift per
  `feedback_vectorization_reduction_drift`). Verify TT/TE < 0.1% before
  regenerating baseline.

### C4 — `Omega_rad` dispatch, dead-code removal, `GetNcdmSpecies`
- Replace the radiation-density special-casts (UR / IDM_DR_IDR /
  IDM_DRMD_IDR_DRMD at 1444–1455) with a virtual returning each species'
  relativistic `Omega0` contribution (default 0; photons/UR/idr return theirs),
  summed in a loop.
- **Delete the dead block** at 1457–1461 (`rho_ncdm_rel_tot = 0; rho_rad += 0`).
- Convert the IC-side `GetNcdmSpecies` use (`GetIni`, 1436–1438) to a dispatch:
  each species reports its required initial `a` (default = the passed value),
  module takes the min.
- **`GetNcdmSpecies` free function:** chase every remaining use. The public
  n-NCDM accessor API (`numberOfNCDMSpecies`, deg / mass / q-size / q at
  ~2142–2163, and verbose printing at 547–604) likely has external callers
  (output module, etc.). Remove the free function only if *every* use can be
  expressed as dispatch; otherwise leave that public API in place and remove only
  the IC-side use. Decide based on what the dependency chase finds — do not break
  the public NCDM accessor contract.

**Verification:** full suite across the model matrix that exercises DCDM_DR,
Fluid (CLP w0/wa), ScalarField, and NCDM. Spectra within ~0.1%; handle Cl^TE
zero-crossings. Regenerate baseline only after confirming acceptable drift.

---

## Part D — Item 4: shared bisection helper (fourth commit)

**Goal:** one small helper replaces the duplicated manual bisections; apply at
all clear-win sites.

- Add a header-only templated helper in `tools/` (e.g. `tools/bisection.h`),
  taking a predicate lambda and a stop condition. The known sites come in two
  shapes:
  - **Value-based refinement:** loop while `(hi - lo) > tol`, `mid = 0.5*(lo+hi)`,
    move `lo`/`hi` by a predicate. (`background_find_equality` second loop 1615–
    1632; `thermodynamics_module.cpp:2400`; `primordial_module.cpp:2213`;
    `perturbations_module.cpp:1363/1392, 2636/2689, 3035/3049`.)
  - **Integer index bracketing:** loop while `(hi - lo) > 1`,
    `mid = (lo+hi)/2`, move by predicate. (`background_find_equality` first loop
    1595–1605.)
- Provide whatever minimal surface covers both (one helper with a stop predicate,
  or two small overloads). Convert each site where the helper reduces lines
  without hurting clarity; leave any site where inlining is genuinely clearer
  (note which, and why).

**Verification:** build + full suite; spectra within tolerance (bisection results
must converge to the same roots).

---

## Out of scope

- Removing `growTable` (still used by thermodynamics).
- Touching the public n-NCDM accessor API beyond what the dependency chase shows
  is safe.
- Any perturbations-module work overlapping #278.

## Commit / verification sequence

1. Part A → build + suite.
2. Part B → build (no-logic move).
3. Part C → full model-matrix suite; regenerate baseline only after confirming
   drift within tolerance.
4. Part D → build + suite.

Each commit is independently buildable and testable.
