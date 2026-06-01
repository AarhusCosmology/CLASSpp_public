# Background Module Cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the dead RK background integrator, relocate `BackgroundColumnWriter`, give species full ownership of their integration indices and initial conditions, and replace duplicated manual bisections with one shared helper.

**Architecture:** Four logically-ordered commits on one branch. Part A deletes the old solver (shrinks the file first). Part B is a pure file move. Part C is the substantive refactor: one registration loop, cached `index_bi_*_` members deleted, IC math pushed into species via a `BackgroundICContext`, `V_scf`/`w_fld` physics relocated onto the species, radiation-`Omega0` and earliest-`a` turned into dispatch. Part D adds a header-only `bisect` helper and applies it everywhere it's a clear win.

**Tech Stack:** C++17, GNU Make build (`make class -j`), Python verification harness (`python/test_class.py` via `pytest -m test_scenario`, and `test/scenarios/compare_tol.py`, RTOL=1e-3, zero-crossing aware).

---

## Verification model (read first)

This codebase has **no unit-test framework for module internals**. Correctness of a
refactor is established by **characterization**: the scenario spectra must not move
beyond ~0.1% (`compare_tol.py`, RTOL=1e-3, which already handles Cl^TE zero-crossings).
The only thing we unit-test directly is the new `bisect` helper (Part D), via the
existing `test-parser`-style standalone-executable pattern.

**Local build + smoke loop used throughout:**

```bash
make class -j          # builds ./class
./class explanatory.ini   # smoke test: must exit 0
```

**Baseline snapshot (capture ONCE, before Part C — the only behavior-touching part):**

```bash
make class -j
mkdir -p /tmp/bgcleanup_baseline
for ini in test/scenarios/gauge_lcdm.ini test/scenarios/gauge_dcdm.ini \
           test/scenarios/gauge_fluid.ini test/scenarios/gauge_scf.ini \
           test/scenarios/gauge_ncdm.ini test/scenarios/gauge_idmdr.ini \
           test/scenarios/dncdm_dr.ini test/scenarios/idm_dr_full.ini \
           test/scenarios/idm_drmd_full.ini; do
  ./class "$ini"
done
cp output/*.dat /tmp/bgcleanup_baseline/   # adjust if a scenario writes elsewhere
```

**Drift check (run after each behavior-touching task in Part C):**

```bash
make class -j
for ini in <same list as above>; do ./class "$ini"; done
python test/scenarios/compare_tol.py /tmp/bgcleanup_baseline output
# Expected: every line "OK ...". Any "FAIL" => investigate before continuing.
```

Parts A, B, D must produce **identical** output (no physics changed); Part C must
stay within tolerance. Regenerate any committed reference baseline only after the
full `pytest -m test_scenario` run passes against the `classyref` build.

---

## Part A — Item 3: delete the old RK solver

### Task A1: Delete `background_solve()` and its `growTable` usage

**Files:**
- Modify: `source/background_module.cpp` (delete ~902–1113)
- Modify: `source/background_module.h:132` (remove declaration)

- [ ] **Step 1: Capture the exact span.** Confirm the function boundaries:

```bash
grep -n "int BackgroundModule::background_solve()" source/background_module.cpp
# then find the matching closing brace before "int BackgroundModule::background_solve_evolver()"
grep -n "int BackgroundModule::background_solve_evolver()" source/background_module.cpp
```

Expected: `background_solve()` starts at ~902, `background_solve_evolver()` starts at ~1199. Delete everything from the `int BackgroundModule::background_solve() {` line through its closing `}` (the line immediately before the doc-comment block that precedes `background_solve_evolver`). This removes the `growTable gTable;` declaration and all `gt_init/gt_add/gt_getPtr/gt_free` calls **inside this function only**.

- [ ] **Step 2: Remove the declaration.** In `source/background_module.h`, delete the line:

```cpp
  int background_solve();
```

- [ ] **Step 3: Verify `growTable` is still referenced elsewhere (must NOT be fully removed).**

Run: `grep -rn "growTable\|gt_init" source/thermodynamics_module.cpp`
Expected: still present (thermodynamics keeps it). Do not touch `tools/growTable.cpp` or `include/growTable.h`.

- [ ] **Step 4: Defer build** — the file won't compile until A2 removes the `switch` case calling `background_solve()`. Proceed directly to A2.

### Task A2: Replace the switch with a direct evolver call; drop the enum + input read

**Files:**
- Modify: `source/background_module.cpp:619–632`
- Modify: `source/background.h:27` (enum) and `:108` (field)
- Modify: `source/input_module.cpp:730`

- [ ] **Step 1: Replace the switch.** In `source/background_module.cpp`, replace:

```cpp
  /** - this function integrates the background over time, allocates
      and fills the background table */
  switch (pba->background_method) {
    case (bgevo_rk):
      class_call(background_solve(), error_message_, error_message_);
      break;
    case (bgevo_evolver):
      class_call(background_solve_evolver(), error_message_, error_message_);
      break;
    default:
      printf(
          "Invalid background method selected. Please set it to 0 or 1 or omit it from your "
          "input.\n");
  }
```

with:

```cpp
  /** - this function integrates the background over time, allocates
      and fills the background table */
  class_call(background_solve_evolver(), error_message_, error_message_);
```

- [ ] **Step 2: Remove the enum and field.** In `source/background.h`, delete:

```cpp
enum background_evolution_method { bgevo_rk, bgevo_evolver };
```

and the field (line ~108):

```cpp
  enum background_evolution_method background_method = bgevo_evolver;
```

- [ ] **Step 3: Remove the input read.** In `source/input_module.cpp`, delete:

```cpp
  class_read_int("background_method", pba->background_method);
```

- [ ] **Step 4: Update the file-header doc comment.** In `source/background_module.cpp` near line 41, remove the bullet describing `background_solve()` (the line beginning `* - background_solve() integrates the quantities {B} and {C} with`). Leave the surrounding comment intact.

- [ ] **Step 5: Defer build** — `dncdm_species.cpp` still references `bgevo_evolver`; fix in A3, then build.

### Task A3: Collapse the `background_method` branch in `dncdm_species.cpp`

**Files:**
- Modify: `species/dncdm_species.cpp:349–427`

- [ ] **Step 1: Make the evolver branch unconditional.** The function `DNCDMSpecies::ComputeBackground` currently has `if (pba_->background_method == bgevo_evolver) { <evolver body> } else { <RK body> }`. The RK `else` branch (lines ~397–427) is now dead. Remove the `if (pba_->background_method == bgevo_evolver) {` line and its matching `else { ... }` block, leaving the evolver body (the `has_problem` logic, lines ~350–395) as the unconditional body. Concretely, after the edit the structure is:

```cpp
void DNCDMSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  double z       = 1. / a_rel - 1.;
  const int q_sz = q_size();

  std::vector<double> lnf_dlnf_array(2 * q_sz);
  std::vector<double> ddlnf_array(q_sz);
  std::vector<double> lnq(q_sz);

  bool has_problem = false;
  for (int i = 0; i < q_sz; i++) {
    if (pvecback_B[index_bi_lnf_decay_dr1_ + i] <= -460) {
      has_problem = true;
      break;
    }
    lnf_dlnf_array[i] = pvecback_B[index_bi_lnf_decay_dr1_ + i];
    lnq[i]            = std::log(q_[i]);
  }
  if (has_problem) {
    for (int i = 0; i < q_sz; i++) {
      double q                               = q_[i];
      double lnf_fd                          = -std::log(1. + std::exp(q));
      double dlnfdlnq_fd                     = -q * std::exp(q) / (1. + std::exp(q));
      pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_fd;
      pvecback[index_bg_dlnfdlnq_decay_ + i] = dlnfdlnq_fd;
      w_bg_[i]                               = std::exp(lnf_fd) * dq_[i];
    }
  }
  else {
    class_call(array_spline_table_lines(lnq.data(), q_sz, lnf_dlnf_array.data(), 1,
                                        ddlnf_array.data(), _SPLINE_EST_DERIV_,
                                        bgm_->error_message_),
               bgm_->error_message_, bgm_->error_message_);
    class_call(array_derive_spline(lnq.data(), q_sz, lnf_dlnf_array.data(),
                                   ddlnf_array.data(), 1, 0, q_sz, bgm_->error_message_),
               bgm_->error_message_, bgm_->error_message_);
    for (int i = 0; i < q_sz; i++) {
      pvecback[index_bg_lnf_decay_dr1_ + i]  = lnf_dlnf_array[i];
      pvecback[index_bg_dlnfdlnq_decay_ + i] = lnf_dlnf_array[q_sz + i];
      pvecback[index_bg_dlnfdlnq_sep_ + i]   = pvecback_B[index_bi_dlnfdlnq_separate_decay_ + i];
      w_bg_[i]                               = std::exp(lnf_dlnf_array[i]) * dq_[i];
    }
  }

  double number_ncdm, rho_ncdm, p_ncdm, pseudo_p_ncdm;
  ComputeMomenta(z, &number_ncdm, &rho_ncdm, &p_ncdm, nullptr, &pseudo_p_ncdm);
  pvecback[index_bg_number_]   = number_ncdm;
  pvecback[index_bg_rho_]      = rho_ncdm;
  pvecback[index_bg_p_]        = p_ncdm;
  pvecback[index_bg_pseudo_p_] = pseudo_p_ncdm;
}
```

(Keep the original multi-line `class_call(...)` formatting from the file; the above is condensed for the plan. The behavior is: the former evolver branch, de-indented one level, with the dead RK `else` removed.)

- [ ] **Step 2: Confirm no other `bgevo_`/`background_method` references remain.**

Run: `grep -rn "bgevo_\|background_method" source/ species/ include/ | grep -v "\.worktrees\|/build/"`
Expected: **no matches.**

- [ ] **Step 3: Build and smoke-test.**

Run: `make class -j && ./class explanatory.ini`
Expected: builds clean; `./class` exits 0.

- [ ] **Step 4: Spectra unchanged.** The default was already `bgevo_evolver`, so output must be identical. Run two representative scenarios and diff against a quick pre-A snapshot if available; otherwise rely on the full `pytest -m test_scenario` at PR time.

Run: `./class test/scenarios/gauge_dcdm.ini && ./class test/scenarios/dncdm_dr.ini`
Expected: both exit 0.

- [ ] **Step 5: Commit.**

```bash
git add source/background_module.cpp source/background_module.h source/background.h \
        source/input_module.cpp species/dncdm_species.cpp
git commit -m "Remove old RK background_solve(); evolver is the only integrator

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Part B — Item 1: relocate `BackgroundColumnWriter`

### Task B1: Move the class into the shared species-facing header

**Files:**
- Modify: `species/base_species.h` (add the class near the background write-hook decls, ~line 185–200)
- Delete: `source/background_column_writer.h`, `source/background_column_writer.cpp`
- Modify: `Makefile:84`, `setup.py:70`
- Modify: any file with `#include "background_column_writer.h"`

- [ ] **Step 1: Find the includes and the hook declarations.**

Run: `grep -rln "background_column_writer.h" source/ species/ | grep -v "\.worktrees\|/build/"`
Run: `grep -n "WriteBackgroundColumnTitles\|class BackgroundColumnWriter" species/base_species.h source/background_module.cpp`
Expected includes: at least `species/base_species.h` and `source/background_module.cpp`.

- [ ] **Step 2: Insert the class into `species/base_species.h`** immediately above the first use of `BackgroundColumnWriter` (the `WriteBackgroundColumnTitles` virtual at ~191). Mirror how `PerturbColumnWriter` lives inline in `species/perturb_source_context.h`. Paste the full class with the `Add` body inlined:

```cpp
// ── BackgroundColumnWriter ───────────────────────────────────────────────────
/**
 * Thin helper for background output — analogous to PerturbColumnWriter.
 * Construct in title mode (titles != nullptr) or data mode (dataptr != nullptr).
 * Call Add() once per column in both modes; the writer handles the branching.
 */
class BackgroundColumnWriter {
 public:
  explicit BackgroundColumnWriter(char* titles) : titles_(titles) {}
  BackgroundColumnWriter(double* dataptr, int& storeidx)
      : dataptr_(dataptr), storeidx_(&storeidx) {}

  bool IsTitleMode() const { return titles_ != nullptr; }

  void Add(const char* title, double value, bool condition = true) {
    if (titles_) {
      class_store_columntitle(titles_, title, condition ? _TRUE_ : _FALSE_);
    }
    else if (dataptr_) {
      class_store_double(dataptr_, value, condition ? _TRUE_ : _FALSE_, (*storeidx_));
    }
  }
  void Add(const std::string& title, double value, bool condition = true) {
    Add(title.c_str(), value, condition);
  }

 private:
  char* titles_    = nullptr;
  double* dataptr_ = nullptr;
  int* storeidx_   = nullptr;
};
```

Ensure `base_species.h` includes what the body needs (`<string>`, and `common.h` for `class_store_*`, `_TRUE_`, `_FALSE_`). Check existing includes first:

Run: `grep -n "#include" species/base_species.h | head`
Add `#include <string>` and/or `#include "common.h"` only if not already present (transitively).

- [ ] **Step 3: Remove the old include lines.** In every file found in Step 1, delete `#include "background_column_writer.h"`. Files that use `BackgroundColumnWriter` and already include `base_species.h` need nothing further; for any that used the writer without including `base_species.h`, add `#include "base_species.h"`.

- [ ] **Step 4: Delete the old files.**

```bash
git rm source/background_column_writer.h source/background_column_writer.cpp
```

- [ ] **Step 5: Update the build.**
  - `Makefile:84`: remove `background_column_writer.opp ` from the `SOURCE = ...` list.
  - `setup.py:70`: remove the `'background_column_writer.cpp',` entry from the sources list.

- [ ] **Step 6: Build and smoke-test.**

Run: `make class -j && ./class explanatory.ini`
Expected: builds clean; exits 0. Background output (`*_background.dat`) byte-identical to before (no logic changed).

- [ ] **Step 7: Commit.**

```bash
git add species/base_species.h Makefile setup.py
git commit -m "Fold BackgroundColumnWriter into base_species.h, matching PerturbColumnWriter

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Part C — Item 2: full index ownership + IC migration

> **Capture the baseline snapshot now** (see Verification model) before any task in this Part.

### Task C1: Single registration loop; delete the four cached `index_bi_*_` members

**Files:**
- Modify: `source/background_module.cpp:838–868` (registration), `:1274` (dcdm read)
- Modify: `source/background_module.h:181–185` (member decls)

- [ ] **Step 1: Replace the special-cased registration block.** In `source/background_module.cpp`, replace lines ~838–868 (the `index_bi_rho_dcdm_ = -1; if (...DCDM_DR...) {...}` block, the "all other species" loop that skips DCDM_DR/Fluid/ScalarField, and the Fluid/ScalarField blocks) with one loop:

```cpp
  /* -> integration indices for all species (each species owns its own offsets) */
  for (auto& [name, sp] : all_species_) {
    sp->RegisterIntegrationIndices(index_bi);
  }
```

- [ ] **Step 2: Delete the cached members.** In `source/background_module.h`, remove:

```cpp
  int index_bi_rho_dcdm_; /**< {B} dcdm density */
  int index_bi_rho_fld_;       /**< {B} fluid density */
  int index_bi_phi_scf_;       /**< {B} scalar field value */
  int index_bi_phi_prime_scf_; /**< {B} scalar field derivative wrt conformal time */
```

- [ ] **Step 3: Repoint the dcdm density read.** At `source/background_module.cpp:1274` (inside `background_solve_evolver`), replace:

```cpp
  if (all_species_.count("DCDM_DR")) {
    Omega0_dcdm_   = pvecback_integration[index_bi_rho_dcdm_] / pba->H0 / pba->H0;
    auto& dcdm_dr  = dynamic_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
    Omega0_dr_    += pvecback_integration[dcdm_dr.dr().bi_rho_index()] / pba->H0 / pba->H0;
  }
```

with (read the dcdm density via the same composite, no cached member):

```cpp
  if (all_species_.count("DCDM_DR")) {
    auto& dcdm_dr  = dynamic_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
    Omega0_dcdm_   = pvecback_integration[dcdm_dr.dcdm().bi_rho_index()] / pba->H0 / pba->H0;
    Omega0_dr_    += pvecback_integration[dcdm_dr.dr().bi_rho_index()] / pba->H0 / pba->H0;
  }
```

(Confirm `DCDM_DR_Species` exposes `dcdm()` returning the `DCDMSpecies` with `bi_rho_index()`: `grep -n "dcdm()\|dr()" species/dcdm_dr_species.h`.)

- [ ] **Step 4: The IC function (`background_initial_conditions`) still references `index_bi_rho_fld_`, `index_bi_phi_scf_`, `index_bi_phi_prime_scf_`.** Those are migrated in C2; until then the file won't compile. Proceed to C2 before building.

### Task C2: Introduce `BackgroundICContext`; push Fluid & ScalarField IC into species

**Files:**
- Create/modify: a header for `BackgroundICContext` (add to `species/base_species.h` near the other context usages, or a small `species/background_ic_context.h` — match where `PerturbSourceContext` lives; it is in `species/perturb_source_context.h`, so create `species/background_ic_context.h`)
- Modify: `species/base_species.h:126` (virtual signature)
- Modify: `species/composite_species.h`, `species/dcdm.h`, `species/dcdm_dr_species.h`, `species/dncdm_species.h`, `species/dncdm_dr_species.h` (existing overrides: signature change)
- Modify: `species/fluid.h`/`.cpp`, `species/scalar_field.h`/`.cpp` (new overrides)
- Modify: `source/background_module.cpp:1420–1524` (IC function)

- [ ] **Step 1: Define the context.** Create `species/background_ic_context.h`:

```cpp
#pragma once

class BackgroundModule;

/**
 * Context passed to BaseSpecies::SetBackgroundInitialConditions().
 * a_rel    : a / a_today at the initial time.
 * rho_rad  : total radiation density at the initial time (units of H0^2),
 *            needed by ScalarField attractor ICs.
 * pvecback_integration : the {B}/{C} integration vector to fill.
 * mod      : back-pointer for module helpers (error_message, etc.).
 */
struct BackgroundICContext {
  double a_rel = 0.;   // a / a_today at the initial time
  double a_ini = 0.;   // absolute initial a (= a_rel * a_today); the value the
                       // module stores at index_bi_a_, exposed so species need
                       // not know the {B} layout
  double rho_rad = 0.;
  double* pvecback_integration = nullptr;
  BackgroundModule* mod = nullptr;
};
```

- [ ] **Step 2: Change the virtual signature.** In `species/base_species.h`, change:

```cpp
  virtual void SetBackgroundInitialConditions(double a_rel, double* pvecback_integration) {}
```

to:

```cpp
  virtual void SetBackgroundInitialConditions(const BackgroundICContext& ctx) {}
```

Add `#include "background_ic_context.h"` to `base_species.h`.

- [ ] **Step 3: Update existing overrides.** For each of `composite_species`, `dcdm`, `dcdm_dr_species`, `dncdm_species`, `dncdm_dr_species`: change the declaration and definition signature from `(double a_rel, double* pvecback_integration)` to `(const BackgroundICContext& ctx)`, and inside each body replace `a_rel` → `ctx.a_rel` and `pvecback_integration` → `ctx.pvecback_integration`. Find them:

Run: `grep -rn "SetBackgroundInitialConditions" species/*.h species/*.cpp | grep -v "\.worktrees"`
For each definition, update the body's two identifiers. (Composite forwards to children — update its forwarding call to pass `ctx` through.)

- [ ] **Step 4: Add the Fluid override.** In `species/fluid.h`, declare:

```cpp
  void SetBackgroundInitialConditions(const BackgroundICContext& ctx) override;
```

In `species/fluid.cpp`, define it, moving the logic currently at `background_module.cpp:1467–1485`. The fluid computes its own `w(a)` integral via its existing `ComputeWFld`:

```cpp
void FluidSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  /* rho_fld today (units of H0^2) */
  const double rho_fld_today = GetOmega0() * pow(ctx.mod->pba->H0, 2);

  /* integrate rho_fld(a) from a_ini to a_0 to get rho_fld(a_ini) */
  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(ctx.a_rel * ctx.mod->pba->a_today, &w_fld, &dw_over_da_fld, &integral_fld);

  ctx.pvecback_integration[bi_rho_index()] = rho_fld_today * exp(integral_fld);
}
```

Confirm `FluidSpecies::ComputeWFld` signature and whether it takes `a` (absolute) or `a_rel`: `grep -n "ComputeWFld" species/fluid.h species/fluid.cpp`. The module previously called `background_w_fld(a, ...)` with `a = a_rel * a_today` (absolute `a`); match that. If `ctx.mod->pba` is not accessible, expose what's needed via an accessor or add `H0`/`a_today` to the context — prefer adding them to `BackgroundICContext` if the back-pointer is awkward. (Decide during implementation; keep the context minimal but sufficient.)

- [ ] **Step 5: Add the ScalarField override.** In `species/scalar_field.h`, declare the same override. In `species/scalar_field.cpp`, define it, moving the logic from `background_module.cpp:1495–1524` and using the species' own `V_scf` and `ctx.rho_rad`:

```cpp
void ScalarFieldSpecies::SetBackgroundInitialConditions(const BackgroundICContext& ctx) {
  double* y = ctx.pvecback_integration;
  const double scf_lambda = scf_parameters()[0];
  if (attractor_ic_scf() == _TRUE_) {
    y[bi_phi_index()] = -1. / scf_lambda *
                        log(ctx.rho_rad * 4. / (3 * pow(scf_lambda, 2) - 12)) *
                        phi_ini_scf();
    if (3. * pow(scf_lambda, 2) - 12. < 0) {
      y[bi_phi_index()] = 1. / scf_lambda;  // no attractor solution; avoid nan
      if (ctx.mod->pba->background_verbose > 0)
        printf(" No attractor IC for lambda = %.3e ! \n ", scf_lambda);
    }
    y[bi_phi_prime_index()] =
        2 * ctx.a_ini *  // was pvecback_integration[index_bi_a_]; see note below
        sqrt(V_scf(y[bi_phi_index()])) * phi_prime_ini_scf();
  }
  else {
    printf("Not using attractor initial conditions\n");
    y[bi_phi_index()]       = phi_ini_scf();
    y[bi_phi_prime_index()] = phi_prime_ini_scf();
  }
}
```

**Note on `index_bi_a_`:** the original used `pvecback_integration[index_bi_a_]` (the scale factor slot, owned by the module, not the species). `index_bi_a_` is `0` (first registered, `class_define_index(index_bi_a_, _TRUE_, index_bi, 1)` at `background_module.cpp:836`). Pass the initial `a` (absolute) through the context instead of reaching for `index_bi_a_`: add `double a_ini = 0.;` to `BackgroundICContext` (set to `ctx.a_rel * a_today`, which equals the value the module stored at `index_bi_a_`), and use `ctx.a_ini` here. This keeps the species from depending on the module's `{B}` layout. The `class_test(isfinite(...))` from lines 1518–1523 can move here too (use `ctx.mod->error_message_`) or stay in the module after the dispatch — keep it in the species for cohesion.

- [ ] **Step 6: Rewrite the module IC function.** In `source/background_module.cpp`, `background_initial_conditions` becomes (replacing the Fluid block 1467–1485 and ScalarField block 1495–1524 with a single dispatch, and threading `rho_rad`):

```cpp
  /* Set initial conditions for all species (each owns its own ODE variables) */
  BackgroundICContext ic;
  ic.a_rel   = a / pba->a_today;
  ic.a_ini   = a;
  ic.rho_rad = rho_rad;             // computed above (see C4 for the dispatch)
  ic.pvecback_integration = pvecback_integration;
  ic.mod     = this;
  for (auto& [name, sp] : all_species_) {
    sp->SetBackgroundInitialConditions(ic);
  }
```

Keep the existing `index_bi_a_`, `index_bi_time_`, `index_bi_tau_`, `index_bi_rs_`, `index_bi_D_`, `index_bi_D_prime_` assignments (the `{C}` quantities and the scale factor) in the module — those are module-owned, not species-owned. The `SetBackgroundInitialConditions` dispatch loop that already existed at 1463–1465 is now the single per-species IC entry point; remove the old separate call if duplicated.

- [ ] **Step 7: Add `#include "background_ic_context.h"`** where `BackgroundModule` constructs the context (`background_module.cpp` — likely transitively via species headers; add explicitly to be safe).

- [ ] **Step 8: Build, smoke, drift-check.**

Run: `make class -j && ./class explanatory.ini`
Then run the drift check (Verification model) over the gauge_* scenarios.
Expected: builds clean; `compare_tol.py` prints all `OK`. Pay special attention to `gauge_fluid.ini` and `gauge_scf.ini` — those exercise the migrated IC. If FAIL, the most likely culprit is a units mismatch in `a` (absolute vs `a_rel`) in the Fluid/ScalarField IC — re-derive against the original lines.

- [ ] **Step 9: Commit.**

```bash
git add species/background_ic_context.h species/base_species.h species/fluid.* \
        species/scalar_field.* species/composite_species.* species/dcdm.* \
        species/dcdm_dr_species.* species/dncdm_species.* species/dncdm_dr_species.* \
        source/background_module.cpp source/background_module.h
git commit -m "Species own integration indices + initial conditions via BackgroundICContext

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task C3: Relocate `V_scf`/`w_fld` physics onto the species

**Files:**
- Modify: `source/background_module.cpp` (`V_scf`/`dV_scf`/`ddV_scf` at ~1905+, `background_functions` scalar-field EOM), `source/background_module.h:23–25`
- Modify: `species/scalar_field.*` (ensure potential is the single source of truth)

- [ ] **Step 1: Map all module-side callers of `V_scf`/`dV_scf`/`ddV_scf`.**

Run: `grep -n "V_scf\|dV_scf\|ddV_scf\|V_e_scf\|V_p_scf" source/background_module.cpp source/background_module.h`
Expected: definitions at ~1905+ and call sites inside `background_functions` (scalar-field EOM) and possibly `background_solve_evolver`.

- [ ] **Step 2: Confirm the species potential matches.** Compare `BackgroundModule::V_scf` (`V_e_scf(phi)*V_p_scf(phi)`) against `ScalarFieldSpecies::V_scf` (`species/scalar_field.cpp:27`). They must compute the identical potential. If the module also defines `V_e_scf`/`V_p_scf`, confirm the species has equivalents (`grep -n "V_e_scf\|V_p_scf" species/scalar_field.*`). If the species lacks any helper, port it there first so the species is the single source of truth.

- [ ] **Step 3: Repoint `background_functions` at the species potential.** In each scalar-field EOM call site, replace `V_scf(phi)` / `dV_scf(phi)` / `ddV_scf(phi)` with a call through the ScalarField species, e.g.:

```cpp
  auto& scf = static_cast<ScalarFieldSpecies&>(*all_species_.at("ScalarField"));
  // ... scf.V_scf(phi), scf.dV_scf(phi), scf.ddV_scf(phi)
```

(This block is already guarded by `all_species_.count("ScalarField")`/the scf index path; fetch the species once outside the inner arithmetic. **Hot path** — do not call `all_species_.at()` per-iteration; hoist it.)

- [ ] **Step 4: Delete the module-side potential.** Remove `BackgroundModule::V_scf`, `dV_scf`, `ddV_scf` definitions (~1905+) and their declarations (`background_module.h:23–25`). If `V_e_scf`/`V_p_scf` are now unused in the module, remove them too (check with grep first).

- [ ] **Step 5: Keep `background_w_fld` as the delegating wrapper.** Do **not** remove `BackgroundModule::background_w_fld` — it already delegates to `FluidSpecies::ComputeWFld` with a no-Fluid fallback and external callers (HyRec/thermo) use it. Verify external callers:

Run: `grep -rn "background_w_fld" source/ | grep -v background_module.cpp`
Expected: at least one external caller — leave the wrapper intact.

- [ ] **Step 6: Build, smoke, drift-check (hot-path sensitive).**

Run: `make class -j && ./class explanatory.ini && ./class test/scenarios/gauge_scf.ini`
Then the full drift check. **Watch for reduction-drift** (memory: `feedback_vectorization_reduction_drift`): if `gauge_scf` drifts slightly but TT/TE < 0.1% via `compare_tol.py`, that is acceptable — do not chase bit-identity.
Expected: `compare_tol.py` all `OK`.

- [ ] **Step 7: Commit.**

```bash
git add source/background_module.cpp source/background_module.h species/scalar_field.*
git commit -m "Move scalar-field potential onto ScalarFieldSpecies; drop module duplicates

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task C4: `Omega_rad` dispatch, dead-code removal, `GetNcdmSpecies` chase

**Files:**
- Modify: `species/base_species.h` (new virtual), `species/photons.*`, `species/ultra_relativistic.*`, the IDR sub-species (`idm_dr_idr_species`, `idm_drmd_idr_drmd_species`)
- Modify: `source/background_module.cpp:1436–1461` (IC radiation block + dead block), and `GetNcdmSpecies` uses

- [ ] **Step 1: Add a relativistic-`Omega0` virtual.** In `species/base_species.h`:

```cpp
  /** Relativistic (radiation-like) Omega0 contribution at early times.
   *  Default 0; photons / UR / interacting dark radiation override. */
  virtual double GetRadiationOmega0() const { return 0.; }
```

- [ ] **Step 2: Override it where the module special-cased.** The module summed `pba->Omega0_g` + UR + idr (IDM_DR_IDR) + idr_drmd (IDM_DRMD_IDR_DRMD). Add overrides:
  - `PhotonsSpecies::GetRadiationOmega0()` → returns its `Omega0` (photons). Confirm how photons expose `Omega0_g`; if photons already return it via `GetOmega0()`, the override can be `return GetOmega0();`.
  - `UltraRelativisticSpecies::GetRadiationOmega0()` → `return GetOmega0();`
  - For composites `IDM_DR_IDR_Species` / `IDM_DRMD_IDR_DRMD_Species`: the radiation piece is the `idr()` / `idr_drmd()` sub-species. The cleanest is for the **sub-species** (the IDR dark-radiation component) to override `GetRadiationOmega0() { return GetOmega0(); }`, and the composite to sum its children's `GetRadiationOmega0()` (composites already iterate children). Verify how the composite aggregates: `grep -n "GetOmega0\|GetRadiationOmega0\|children" species/composite_species.* species/idm_dr_idr_species.* species/idm_drmd_idr_drmd_species.*`.

- [ ] **Step 3: Replace the module radiation sum.** In `source/background_module.cpp`, replace lines ~1444–1461:

```cpp
  double Omega_rad = pba->Omega0_g;
  if (all_species_.count("UR"))
    Omega_rad += all_species_.at("UR")->GetOmega0();
  if (all_species_.count("IDM_DR_IDR")) { ... }
  if (all_species_.count("IDM_DRMD_IDR_DRMD")) { ... }
  double rho_rad = Omega_rad * pow(pba->H0, 2) / pow(a / pba->a_today, 4);
  if (!GetNcdmSpecies(all_species_).empty()) {
    double rho_ncdm_rel_tot  = 0.;
    rho_rad                 += rho_ncdm_rel_tot;
  }
```

with:

```cpp
  double Omega_rad = pba->Omega0_g;
  for (auto& [name, sp] : all_species_)
    Omega_rad += sp->GetRadiationOmega0();
  double rho_rad = Omega_rad * pow(pba->H0, 2) / pow(a / pba->a_today, 4);
```

This **deletes the dead `rho_ncdm_rel_tot = 0` block** (it added zero). If `pba->Omega0_g` is itself owned by the Photons species, prefer dropping the seed and letting Photons' `GetRadiationOmega0()` supply it — confirm there isn't double-counting (Photons override + `pba->Omega0_g` seed). Pick exactly one source for the photon term.

- [ ] **Step 4: Chase `GetNcdmSpecies`.** Enumerate every use:

Run: `grep -n "GetNcdmSpecies" source/background_module.cpp`
Categorize:
  - **IC earliest-`a` use (1436–1438):** convert to dispatch. Add a virtual to `base_species.h`:
    ```cpp
    /** Earliest a/a_today this species needs integration to start from.
     *  Default: returns the proposed a unchanged. NCDM species may pull it earlier. */
    virtual double BackgroundAIni(double a_proposed, double a_today, double tol) const { return a_proposed; }
    ```
    Override in the NCDM base species to wrap the existing `GetIni`. Replace the module loop with:
    ```cpp
    for (auto& [name, sp] : all_species_)
      a = std::min(a, sp->BackgroundAIni(a, pba->a_today, ppr->tol_ncdm_initial_w));
    ```
    (Confirm `GetIni` returns a *smaller* `a` to start earlier; `std::min` keeps the earliest. Check the sign/semantics of `GetIni` at its definition before finalizing.)
  - **Verbose mass printing (547–604), n-NCDM accessor API (2142–2163):** these back a public API (`numberOfNCDMSpecies`, `GetNCDMMass`, q-size/q). **Leave them and the `GetNcdmSpecies` free function in place** unless every caller can be expressed as dispatch without breaking the public contract. Document in the commit message which uses were converted and which were intentionally kept.

- [ ] **Step 5: Build, smoke, drift-check across the full model matrix.**

Run: `make class -j && ./class explanatory.ini`
Then the drift check over **all** gauge_* scenarios plus `idm_dr_full.ini`, `idm_drmd_full.ini`, `ncdm_*` scenarios (these exercise the radiation dispatch and earliest-`a`).
Expected: `compare_tol.py` all `OK`. `gauge_idmdr.ini`, `idm_dr_full.ini`, `idm_drmd_full.ini`, and `gauge_ncdm.ini` are the sensitive ones.

- [ ] **Step 6: Confirm the module no longer names dynamical species for indices/ICs.**

Run: `grep -n "static_cast<.*Species\|dynamic_cast<.*Species\|all_species_.at(\"Fluid\")\|all_species_.at(\"ScalarField\")\|all_species_.at(\"UR\")\|IDM_DR_IDR\|IDM_DRMD" source/background_module.cpp`
Expected: remaining hits should be only the legitimate ones (e.g., the `DCDM_DR` density read in C1, the `background_w_fld` wrapper, the `numberOfNCDMSpecies` public API). No species-typed reach-ins remain in `background_initial_conditions` or the registration block.

- [ ] **Step 7: Commit.**

```bash
git add species/base_species.h species/photons.* species/ultra_relativistic.* \
        species/idm_dr_idr_species.* species/idm_drmd_idr_drmd_species.* \
        species/ncdm_base_species.* source/background_module.cpp
git commit -m "Radiation Omega0 and earliest-a via species dispatch; drop dead ncdm-rel block

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Part D — Item 4: shared bisection helper

### Task D1: Write the helper with a failing unit test

**Files:**
- Create: `tools/bisection.h`
- Create: `tools/bisection_test.cpp`
- Modify: `Makefile` (add a `test-bisection` target mirroring `test-parser` at line 112)

- [ ] **Step 1: Write the failing test.** Create `tools/bisection_test.cpp`:

```cpp
#include "bisection.h"
#include <cassert>
#include <cmath>
#include <cstdio>

int main() {
  // Value-based: find root of f(x)=x-2 in [0,10] to tol 1e-9.
  double root = bisect_value(0.0, 10.0, 1e-9,
                             [](double x) { return x > 2.0; });  // predicate: f(x) crossed?
  assert(std::fabs(root - 2.0) < 1e-6);

  // Integer-bracket: smallest index where table[i] >= 5 in 0..9 (table[i]=i).
  int idx = bisect_index(0, 9, [](int i) { return i >= 5; });
  assert(idx == 5);

  std::printf("bisection tests passed\n");
  return 0;
}
```

- [ ] **Step 2: Add the build target.** In `Makefile`, after the `test-parser` rule (~112–113), add:

```make
test-bisection: common.opp
	$(CXX) $(OPTFLAG) $(CXXFLAG) -Iinclude -Itools -Isource -I. tools/bisection_test.cpp -o test-bisection $(LIBRARIES)
```

- [ ] **Step 3: Run the test to verify it fails (no header yet).**

Run: `make test-bisection && ./test-bisection`
Expected: FAIL — compile error, `bisection.h` not found / `bisect_value` undeclared.

- [ ] **Step 4: Write the helper.** Create `tools/bisection.h`:

```cpp
#pragma once

/**
 * Minimal bisection helpers, header-only, replacing the hand-rolled loops
 * scattered across the modules.
 *
 * bisect_value: continuous bracketing. `lo`,`hi` bracket the root; `pred(mid)`
 *   returns true when `mid` is on the `hi` side of the root. Iterates until
 *   (hi - lo) <= tol; returns the final midpoint.
 *
 * bisect_index: integer bracketing. `lo`,`hi` are indices with pred(lo)==false,
 *   pred(hi)==true; returns the smallest index in (lo,hi] where pred flips true,
 *   matching the `(hi - lo) > 1` table-bracketing loops.
 */
template <typename Pred>
double bisect_value(double lo, double hi, double tol, Pred pred) {
  double mid = 0.5 * (lo + hi);
  while ((hi - lo) > tol) {
    mid = 0.5 * (lo + hi);
    if (pred(mid)) hi = mid;
    else           lo = mid;
  }
  return mid;
}

template <typename Pred>
int bisect_index(int lo, int hi, Pred pred) {
  while ((hi - lo) > 1) {
    int mid = (lo + hi) / 2;
    if (pred(mid)) hi = mid;
    else           lo = mid;
  }
  return hi;
}
```

- [ ] **Step 5: Run the test to verify it passes.**

Run: `make test-bisection && ./test-bisection`
Expected: prints `bisection tests passed`, exits 0.

- [ ] **Step 6: Commit.**

```bash
git add tools/bisection.h tools/bisection_test.cpp Makefile
git commit -m "Add header-only bisect_value/bisect_index helpers with unit test

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task D2: Apply the helper to `background_find_equality`

**Files:**
- Modify: `source/background_module.cpp:1587–1645`

- [ ] **Step 1: Add the include.** At the top of `source/background_module.cpp`, add `#include "bisection.h"` among the tools includes.

- [ ] **Step 2: Replace the integer bracketing loop (1595–1605).** Replace:

```cpp
  while ((index_tau_plus - index_tau_minus) > 1) {
    index_tau_mid = (int) (0.5 * (index_tau_plus + index_tau_minus));
    Omega_m_over_Omega_r = background_table_[index_tau_mid * bg_size_ + index_bg_Omega_m_] /
                           background_table_[index_tau_mid * bg_size_ + index_bg_Omega_r_];
    if (Omega_m_over_Omega_r > 1)
      index_tau_plus = index_tau_mid;
    else
      index_tau_minus = index_tau_mid;
  }
```

with:

```cpp
  index_tau_plus = bisect_index(index_tau_minus, index_tau_plus, [&](int i) {
    return background_table_[i * bg_size_ + index_bg_Omega_m_] /
               background_table_[i * bg_size_ + index_bg_Omega_r_] > 1.;
  });
  index_tau_minus = index_tau_plus - 1;
```

(`bisect_index` returns the upper bracket; `index_tau_minus`/`index_tau_plus` then frame the root for the refinement step exactly as before.)

- [ ] **Step 3: Replace the value refinement loop (1615–1632).** Replace the `while ((tau_plus - tau_minus) > ppr->tol_tau_eq) { ... }` block with:

```cpp
  tau_mid = bisect_value(tau_minus, tau_plus, ppr->tol_tau_eq, [&](double tau) {
    class_call(background_at_tau(tau, pba->long_info, pba->inter_closeby,
                                 &index_tau_minus, pvecback.data()),
               error_message_, error_message_);
    return pvecback[index_bg_Omega_m_] / pvecback[index_bg_Omega_r_] > 1.;
  });
```

**Caveat:** the lambda must propagate `class_call` errors. Since `class_call` returns from the enclosing function on error, it cannot live inside a lambda that returns `bool`. Handle this by either (a) keeping an error flag captured by reference and checking after, or (b) leaving this particular loop inline if the error-handling macro makes the helper awkward. **If (b), document it**: the value helper fits the clean loops (thermo/primordial/perturbations) better than ones with `class_call` inside. Decide per the actual macro expansion; do not force the helper where it harms clarity or correctness.

- [ ] **Step 4: Build and drift-check.**

Run: `make class -j && ./class explanatory.ini`
Then drift check (equality redshift appears in verbose output and affects nothing spectral if unchanged). Expected: `compare_tol.py` all `OK`.

- [ ] **Step 5: Commit.**

```bash
git add source/background_module.cpp
git commit -m "Use bisect_index/bisect_value in background_find_equality

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task D3: Apply the helper to the other clear-win sites

**Files:**
- Modify: `source/thermodynamics_module.cpp:2400`, `source/primordial_module.cpp:2213`, `source/perturbations_module.cpp` (1363/1392, 2636/2689, 3035/3049)

- [ ] **Step 1: Triage each site.** For each location, read the surrounding loop and classify:
  - Pure value bisection with a side-effect-free predicate and a `(hi-lo) > tol` stop → **convert** to `bisect_value`.
  - Loop containing `class_call`/early-return inside the predicate, or a non-standard stop/update → **leave inline** and note why in the commit message.

Run for each: `grep -n "0.5 \* (.*+.*)\|while ((.*-.*) >" source/thermodynamics_module.cpp source/primordial_module.cpp source/perturbations_module.cpp`

- [ ] **Step 2: Convert the clear wins.** Add `#include "bisection.h"` to each module touched. For a representative value site, e.g. `thermodynamics_module.cpp:2400`:

```cpp
  // before:
  // while ((tau_sup - tau_inf) > tau_reionization_ * ppr->reionization_optical_depth_tol) {
  //   double z_mid = 0.5 * (z_sup + z_inf);
  //   ... compute optical depth at z_mid ...
  //   if (<too much depth>) z_inf = z_mid; else z_sup = z_mid;  // (match original direction!)
  // }
  // after: express the same predicate and tol via bisect_value, preserving the
  // exact comparison direction and the variable that is the convergence target.
```

**Critical:** preserve each loop's exact predicate direction and convergence variable. Several of these bisect on `z` while the stop is on `tau`; if the stop variable differs from the bisected variable, the helper's `(hi-lo) > tol` on the bisected variable is **not** equivalent — **leave those inline**. Only convert sites where the bisected variable and the stop variable are the same.

- [ ] **Step 3: Build and full drift-check.**

Run: `make class -j && ./class explanatory.ini`
Then drift check over the full scenario list (thermo/reionization and primordial shooting affect spectra directly).
Expected: `compare_tol.py` all `OK`.

- [ ] **Step 4: Commit.**

```bash
git add source/thermodynamics_module.cpp source/primordial_module.cpp source/perturbations_module.cpp
git commit -m "Use shared bisect helper at clear-win sites; note sites left inline

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Final verification (before opening the PR)

- [ ] **Full build + smoke:** `make class -j && ./class explanatory.ini`
- [ ] **Scenario suite vs reference build** (the canonical PR gate):

```bash
cd python
COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py
```

Expected: all scenarios pass within tolerance. This requires the `classyref` reference build per `.github/workflows/test_on_pull_request.yml`; if running locally, build `master` as `classyref` first (see that workflow's "make reference" step).

- [ ] **Confirm the cleanup goals are met:**
  - `grep -rn "bgevo_\|background_method\|background_solve\b" source/ species/ include/ | grep -v "\.worktrees\|/build/"` → no matches.
  - `ls source/background_column_writer.*` → no such files.
  - `grep -n "index_bi_rho_dcdm_\|index_bi_rho_fld_\|index_bi_phi_scf_\|index_bi_phi_prime_scf_" source/background_module.*` → no matches.
  - `git grep -c "0.5 \* (.*plus.*minus" source/background_module.cpp` → reduced.

- [ ] **Open the PR** referencing this plan and the design spec, with a body listing the four commits and the scenario-comparison evidence.

---

## Self-review notes (coverage vs spec)

- Spec Part A → Tasks A1–A3 (solver delete, switch/enum/input, dncdm collapse). `growTable` retention asserted in A1 Step 3. ✓
- Spec Part B → Task B1 (move into `base_species.h`, delete files, build lists). ✓
- Spec Part C1 (loop + drop members) → Task C1. ✓
- Spec Part C2 (IC context, Fluid/ScalarField overrides) → Task C2. ✓
- Spec Part C3 (`V_scf`/`w_fld` relocation, hot-path caution) → Task C3. ✓
- Spec Part C4 (`Omega_rad` dispatch, dead block, `GetNcdmSpecies` chase) → Task C4. ✓
- Spec Part D (helper + all clear-win sites) → Tasks D1–D3, with explicit guidance to leave `class_call`-bearing and mismatched-variable loops inline. ✓

**Known implementation-time decisions flagged in tasks (not placeholders — genuine forks to resolve against the code):** exact contents of `BackgroundICContext` (back-pointer vs. copied `H0`/`a_today`); whether `pba->Omega0_g` seed vs. a Photons `GetRadiationOmega0()` override (avoid double count); which D3 sites are safe to convert; whether the `background_find_equality` value loop keeps `class_call` inline.
