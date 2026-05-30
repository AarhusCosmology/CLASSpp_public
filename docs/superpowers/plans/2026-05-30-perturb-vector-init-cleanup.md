# Clean up `perturb_vector_init` — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the three accidental complexities in `PerturbationsModule::perturb_vector_init` (`source/perturbations_module.cpp:3269-3747`) — precision tests on the hot path, duplicated registration calls, and a long sequence of per-species `if (all_species_.count("X"))` blocks — with one dispatch loop per mode plus a single upstream test block. Behavior shifts the integration `y`-vector layout for DNCDM_DR setups; Cl spectra stay within ~0.1% per project tolerance.

**Architecture:**
1. Every `class_test(ppr->l_max_X < 4, …)` moves once into `InputModule::input_read_parameters`, unconditional.
2. The legacy no-layout `RegisterPerturbationIndices(perturb_vector*, …)` overload is deleted; the layout overload is the only registration path. `DNCDM_DR_Species` gains a single-position override matching the existing `DCDM_DR_Species` pattern, so its DR + DNCDM children register contiguously at its lex position.
3. Each species reads its own `ppr->l_max_*` inside the register call instead of having the module pre-populate the layout.
4. `MarkUsedInSources` gains a `ppw` argument and splits into a `Register{,Vector,Tensor}`-mirroring trio; each approximation-conditional source mask folds into the owning species.

**Tech Stack:** C++17, CLASSpp Makefile (`make -j8 class`), pytest scenario harness (`python -m pytest -m test_scenario python/test_class.py`), `test/scenarios/compare_tol.py` for Cl-level diffing at ~0.1%.

**Design spec:** `docs/superpowers/specs/2026-05-30-perturb-vector-init-cleanup-design.md` (head `7b42c7b2`).

**Verification model:**
- Per-task: build clean + `TEST_LEVEL=1 python -m pytest -m test_scenario python/test_class.py` passes (the harness checks "runs without crashing + produces correct-shape output"; small Cl drifts don't fail it).
- Tight numerical verification happens twice: **once before any code change** (Task 0 captures master baselines into `test/scenarios/_baseline_master/`), and **at the end of Task 5** (after the y-vector reorder), comparing the new outputs against the captured baseline with `compare_tol.py` at 0.1%.
- Final task runs `TEST_LEVEL=3` (1260 scenarios).

**Working directory:** `/Users/au192734/Projects/class_claude/.claude/worktrees/remove-species-physics-params` (branch `worktree-remove-species-physics-params`).

---

## File Structure

| File | Change |
|---|---|
| `source/input_module.cpp` | Add 8 unconditional `class_test` calls before `return _SUCCESS_;` at end of `input_read_parameters` (line ~2722). |
| `source/perturbations_module.cpp` | Replace per-species scaffolding in `perturb_vector_init` (lines 3269-3747) with dispatch loops; delete in-place precision tests, drop the duplicated `RegisterPerturbationIndices(ppv, …)` calls, delete per-species `MarkUsedInSources` masks. |
| `species/base_species.h` | Remove the three `Register*(perturb_vector*, …)` virtuals; expand `MarkUsedInSources` signature and add `MarkVectorUsedInSources`, `MarkTensorUsedInSources` no-op defaults. |
| `species/composite_species.{h,cpp}` | Remove the three `Register*(perturb_vector*, …)` forwarders. |
| `species/photons.{h,cpp}` | Remove the no-layout stubs. Implement scalar + tensor `MarkUsedInSources` / `MarkTensorUsedInSources` body. |
| `species/baryons.{h,cpp}` | Remove no-layout stub. Implement `RegisterVectorPerturbationIndices(layout, …)`. |
| `species/cdm.{h,cpp}`, `lambda.h`, `fluid.{h,cpp}`, `scalar_field.{h,cpp}`, `dcdm.{h,cpp}`, `idm_dr.{h,cpp}`, `idr.{h,cpp}`, `idm_drmd.{h,cpp}`, `idr_drmd.{h,cpp}`, `dark_radiation_species.{h,cpp}`, `dncdm_species.{h,cpp}` | Remove the no-layout stubs. |
| `species/ultra_relativistic.{h,cpp}` | Remove no-layout stub. Read `ppr->l_max_ur` inside `RegisterPerturbationIndices`. Implement scalar `MarkUsedInSources` body. |
| `species/ncdm_species.{h,cpp}` | Remove no-layout stub. Move tensor pre-loop (`l_max = ppr->l_max_ncdm; q_size = q_size()`) into `RegisterTensorPerturbationIndices`. Implement scalar `MarkUsedInSources` body. |
| `species/idm_dr_idr_species.{h,cpp}` | Implement scalar `MarkUsedInSources` body. |
| `species/dncdm_dr_species.{h,cpp}` | Add `RegisterPerturbationIndices(layout, …)` override (DR child first, DNCDM second). Drop the "intentionally not overridden" header comment. |

---

## Task 0: Baseline + survey

Confirm the worktree starts clean, the test harness works, and capture master Cl baselines for the y-vector-reorder verification later.

**Files:** none modified. Creates `test/scenarios/_baseline_master/` (gitignored — directory not committed).

- [ ] **Step 1: Verify clean working tree**

```bash
git status --short
```

Expected: empty output (head is `7b42c7b2`, spec commit).

- [ ] **Step 2: Verify build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: last line is the `g++ -O3 ... -o class` link, no errors.

- [ ] **Step 3: Install classy**

```bash
pip install . 2>&1 | tail -3
```

Expected: `Successfully installed classy-community-...`.

- [ ] **Step 4: Run scenario regression**

```bash
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -5
```

Expected: `84 passed, 98 deselected` (some seconds; passes).

- [ ] **Step 5: Capture master Cl baselines for the DNCDM_DR and adjacent scenarios**

```bash
mkdir -p test/scenarios/_baseline_master
for ini in test/scenarios/gauge_lcdm.ini \
           test/scenarios/gauge_ncdm.ini \
           test/scenarios/gauge_dcdm.ini \
           test/scenarios/gauge_idmdr.ini \
           test/scenarios/gauge_fluid.ini \
           test/scenarios/gauge_scf.ini \
           test/scenarios/dncdm_dr.ini \
           test/scenarios/ncdm_dncdm_idmdr_combined.ini \
           test/scenarios/ncdm_multi_unsorted.ini \
           test/scenarios/ncdm_self_interacting.ini \
           test/scenarios/idm_dr_full.ini \
           test/scenarios/idm_drmd_full.ini ; do
  ./class "$ini" 2>&1 | tail -2
done
# Move every .dat produced this run into the baseline dir. The .ini files
# direct their output via root= settings; collect the resulting cl/pk files.
mv output/*.dat test/scenarios/_baseline_master/ 2>/dev/null || true
ls test/scenarios/_baseline_master/ | head -10
```

Expected: a non-empty listing of `*_cl_lensed.dat` / `*_cl.dat` / `*_pk.dat` files. **Do not commit this directory** — it lives across the branch only.

- [ ] **Step 6: No commit** — baseline capture only.

---

## Task 1: Move precision tests upstream

Add unconditional `class_test` calls to `InputModule::input_read_parameters` and delete the gated per-block tests inside `perturb_vector_init`.

**Files:**
- Modify: `source/input_module.cpp` (insertion at end of `input_read_parameters`, line ~2722)
- Modify: `source/perturbations_module.cpp` (delete lines 3294-3336, 3492-3500, 3533-3541, 3571-3576)

- [ ] **Step 1: Add the test block at end of `input_read_parameters`**

In `source/input_module.cpp`, immediately before the final `return _SUCCESS_;` of `input_read_parameters` (line ~2722), add:

```cpp
  /* ── Precision-consistency tests for perturbation-hierarchy l_max values.
     Moved here from perturb_vector_init: these checks depend only on ppr
     (plus ppt->idr_nature for the IDR test), so they belong at input-parse
     time, not inside the per-(k, approximation) hot path.  Tests run
     unconditionally — a too-low l_max is a user-config error whether or not
     the species ends up active. */
  class_test(ppr->l_max_g < 4, errmsg,
             "ppr->l_max_g should be at least 4, i.e. we must integrate at least over photon "
             "density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_pol_g < 4, errmsg,
             "ppr->l_max_pol_g should be at least 4");
  class_test(ppr->l_max_ur < 4, errmsg,
             "ppr->l_max_ur should be at least 4, i.e. we must integrate at least over "
             "neutrino/relic density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_dr < 4, errmsg,
             "ppr->l_max_dr should be at least 4, i.e. we must integrate at least over "
             "neutrino/relic density, velocity, shear, third and fourth momentum");
  class_test((ppr->l_max_idr < 4) && (ppt->idr_nature == idr_free_streaming), errmsg,
             "ppr->l_max_idr should be at least 4, i.e. we must integrate at least over "
             "interacting dark radiation density, velocity, shear, third and fourth momentum");
  class_test(ppr->l_max_g_ten < 4, errmsg,
             "ppr->l_max_g_ten should be at least 4, i.e. we must integrate at least over photon "
             "density, velocity, shear, third momentum");
  class_test(ppr->l_max_pol_g_ten < 4, errmsg,
             "ppr->l_max_pol_g_ten should be at least 4");
  class_test(ppr->l_max_ncdm < 4, errmsg,
             "ppr->l_max_ncdm=%d should be at least 4, i.e. we must integrate at least "
             "over first four momenta of non-cold dark matter perturbed phase-space "
             "distribution", ppr->l_max_ncdm);

```

(`errmsg` is the local error buffer already in scope in `input_read_parameters`.)

- [ ] **Step 2: Delete the gated tests inside `perturb_vector_init`**

In `source/perturbations_module.cpp`, delete:

- Lines 3294-3336 (the `_scalars_` precision tests: `l_max_g`, `l_max_pol_g`, the `has_any_dr_species` scan + `l_max_dr` test, `l_max_ur` test, IDR `l_max_idr` test).
- Lines 3492-3500 (the `_vectors_` `l_max_g_ten` and `l_max_pol_g_ten` tests).
- Lines 3533-3541 (the `_tensors_` `l_max_g_ten` and `l_max_pol_g_ten` tests).
- Lines 3571-3576 (the `evolve_tensor_ncdm_` `l_max_ncdm` test inside the NCDM tensor block).

After deletion, the `if (_scalars_) {` block opens directly to the photon registration block, and similarly for vectors/tensors.

- [ ] **Step 3: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 4: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 5: Commit**

```bash
git add source/input_module.cpp source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
refactor: move perturbation l_max tests to input

Moves the eight class_test(ppr->l_max_X < 4, ...) checks out of
perturb_vector_init (where they ran once per (k, approx-scheme) and were
hand-gated on species presence with a dynamic_cast scan for DR-emitting
species) into a single unconditional block at the end of input_read_
parameters. Tests now fire once per cosmology at config-parse time.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Drop the duplicate `RegisterPerturbationIndices(ppv, …)` call sites

The second call after each layout-based registration is a no-op stub on every species (and a child-forwarder on `CompositeSpecies` that calls the same stubs). Remove these calls from the module first so the API removal in Task 3 builds cleanly.

**Files:**
- Modify: `source/perturbations_module.cpp` (lines 3344, 3353, 3362, 3371, 3380, 3389, 3419, 3428, 3446, 3465, 3472, 3519, 3551)

- [ ] **Step 1: Delete each duplicate no-layout call**

Inside `perturb_vector_init`, delete each line that calls the no-layout overload after a layout call. The pattern is (current line numbers in parens):

```cpp
sp->RegisterPerturbationIndices(layout, ppv, ppr, index_pt, ppw, ppt->gauge);
sp->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);   // ← DELETE
```

Specifically delete the no-layout call from:

- Photons scalar block (line 3344)
- Baryons scalar block (line 3353)
- CDM scalar block (line 3362)
- IDM_DR_IDR scalar block (line 3371)
- IDM_DRMD_IDR_DRMD scalar block (line 3380)
- DCDM_DR scalar block (line 3389)
- Fluid scalar block (line 3419)
- ScalarField scalar block (line 3428)
- UR scalar block (line 3446)
- NCDM dispatch loop inside `HasNcdm(...)` (line 3465: `n->RegisterPerturbationIndices(ppv, …)`)
- NCDM dispatch loop, DNCDM_DR branch (line 3472: `composite->dncdm().RegisterPerturbationIndices(ppv, …)`)
- Photons vector block (line 3519)
- Photons tensor block (line 3551)

Keep every layout-based call.

- [ ] **Step 2: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build (the no-layout overloads still exist as `= 0` pure virtuals in `BaseSpecies` plus `{}` stubs everywhere; removing call sites doesn't break the API).

- [ ] **Step 3: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed` (the deleted calls were no-ops, so behavior is unchanged).

- [ ] **Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
refactor: drop dead RegisterPerturbationIndices(ppv,...) call sites

The second of each pair in perturb_vector_init invoked the no-layout
overload, which has been an empty stub on every species since the layout
migration. CompositeSpecies's implementation just iterates children's
stubs. Remove the 13 call sites so the overload itself can be deleted
in the next step.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Delete the legacy no-layout `Register*` overloads

Remove the three virtuals from `BaseSpecies`, the `CompositeSpecies` forwarders, and every species-level stub override.

**Files:**
- Modify: `species/base_species.h` (lines 227-243)
- Modify: `species/composite_species.h` and `species/composite_species.cpp` (3 method declarations + bodies)
- Modify: each species header that overrides any of the three (≈20 files; full list in File Structure above)

- [ ] **Step 1: Remove from `BaseSpecies`**

In `species/base_species.h`, delete:

```cpp
  virtual void RegisterPerturbationIndices(perturb_vector* pv,
                                           const precision* ppr,
                                           int& index_pt,
                                           const perturb_workspace* ppw,
                                           int gauge) = 0;

  virtual void RegisterVectorPerturbationIndices(perturb_vector* /*pv*/,
                                                 int& /*index_pt*/,
                                                 const perturb_workspace* /*ppw*/,
                                                 int /*gauge*/) {}

  virtual void RegisterTensorPerturbationIndices(perturb_vector* /*pv*/,
                                                 int& /*index_pt*/,
                                                 const perturb_workspace* /*ppw*/,
                                                 int /*gauge*/) {}
```

(approx lines 227-243). The layout-based overloads declared further down (lines 249-266 in master) stay.

- [ ] **Step 2: Remove from `CompositeSpecies`**

In `species/composite_species.h`, delete the three method declarations matching:

```cpp
void RegisterPerturbationIndices(perturb_vector* pv, ...) override;
void RegisterVectorPerturbationIndices(perturb_vector* pv, ...) override;
void RegisterTensorPerturbationIndices(perturb_vector* pv, ...) override;
```

In `species/composite_species.cpp`, delete the matching bodies (lines 13-36 in current file). The other `CompositeSpecies` methods (`RegisterBackgroundIndices`, `RegisterIntegrationIndices`, all the perturb-derivs forwarders, etc.) stay untouched.

- [ ] **Step 3: Remove the stub overrides from every species header**

Each of these files has an override that looks like:

```cpp
void RegisterPerturbationIndices(perturb_vector* /*pv*/,
                                 const precision* /*ppr*/,
                                 int& /*index_pt*/,
                                 const perturb_workspace* /*ppw*/,
                                 int /*gauge*/) override {}
```

(plus similar `Vector` / `Tensor` variants in `photons.h`). Delete every such block from:

- `species/photons.h` (lines 69-73, 82-85, plus the tensor variant if present)
- `species/baryons.h`
- `species/cdm.h` (lines 49-53)
- `species/lambda.h` (line 46)
- `species/fluid.h`
- `species/scalar_field.h`
- `species/dcdm.h` (lines 65-69)
- `species/dark_radiation_species.h` (lines 71-…)
- `species/dncdm_species.h` (lines 84-…)
- `species/ncdm_species.h` (lines 55-…)
- `species/ultra_relativistic.h` (lines 62-…)
- `species/idm_dr.h` (lines 60-…)
- `species/idr.h` (lines 116-…)
- `species/idm_drmd.h` (lines 59-…)
- `species/idr_drmd.h` (lines 64-…)

(Composites like `dcdm_dr_species.h`, `idm_dr_idr_species.h`, `idm_drmd_idr_drmd_species.h`, `dncdm_dr_species.h` inherit from `CompositeSpecies` and have no direct override of the no-layout method, so nothing to delete from them — confirm by `grep "perturb_vector\* /\*pv\*/" species/<file>.h`.)

To find every file that needs editing:

```bash
grep -rln "RegisterPerturbationIndices(perturb_vector\* /\*pv\*/" species/
grep -rln "RegisterVectorPerturbationIndices(perturb_vector\* /\*pv\*/" species/
grep -rln "RegisterTensorPerturbationIndices(perturb_vector\* /\*pv\*/" species/
```

Delete the stub block in each match.

- [ ] **Step 4: Build**

```bash
make -j8 class 2>&1 | tail -10
```

Expected: clean build. (If anything still calls the old API, the compiler complains; that's the verification.)

- [ ] **Step 5: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 6: Commit**

```bash
git add species/ source/
git commit -m "$(cat <<'EOF'
refactor: delete legacy no-layout RegisterPerturbationIndices overload

The (perturb_vector*, ppr, index_pt, ppw, gauge) overload of
Register{,Vector,Tensor}PerturbationIndices was the transitional API
during the layout migration. Every species now stubs it to {} and
CompositeSpecies just forwards to those stubs. Drop the pure virtuals
from BaseSpecies, the forwarders from CompositeSpecies, and every
stub override.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: DNCDM_DR composite registers its own children

Today `perturb_vector_init` splits `DNCDM_DR_Species` across two locations (DR children with photons-side, DNCDM children with NCDM-side) because the composite has no `RegisterPerturbationIndices` override. Match the pattern `DCDM_DR_Species` already uses.

**Files:**
- Modify: `species/dncdm_dr_species.h`
- Modify: `species/dncdm_dr_species.cpp`

- [ ] **Step 1: Declare the override**

In `species/dncdm_dr_species.h`, in the public section, add (and delete the existing `NOTE: RegisterPerturbationIndices is intentionally NOT overridden here.` comment if it sits next to the declarations):

```cpp
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
```

- [ ] **Step 2: Implement**

In `species/dncdm_dr_species.cpp`, add at the bottom of the file (or alongside other methods):

```cpp
void DNCDM_DR_Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                   perturb_vector* pv,
                                                   const precision* ppr,
                                                   int& index_pt,
                                                   const perturb_workspace* ppw,
                                                   int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  /* DR child first, DNCDM child second — matches legacy registration order
     (DR was registered with the photon-side block; DNCDM with the NCDM-side
     block).  Inside one composite, DR slots come before DNCDM slots in pv->y. */
  dr_sp_->RegisterPerturbationIndices(my.dr,    pv, ppr, index_pt, ppw, gauge);
  dncdm_->RegisterPerturbationIndices(my.dncdm, pv, ppr, index_pt, ppw, gauge);
}
```

(Private member names confirmed: `DNCDMSpecies* dncdm_` and `DarkRadiationSpecies* dr_sp_` — `species/dncdm_dr_species.h:153-154`. Sub-layout names confirmed: `PerturbLayout::dncdm` and `PerturbLayout::dr` — `species/dncdm_dr_species.h:22-23`.)

- [ ] **Step 3: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build. (This task doesn't change behavior yet — the module still routes DNCDM_DR via the manual split; this override just becomes available.)

- [ ] **Step 4: Regression**

```bash
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 5: Commit**

```bash
git add species/dncdm_dr_species.h species/dncdm_dr_species.cpp
git commit -m "$(cat <<'EOF'
refactor: DNCDM_DR_Species overrides RegisterPerturbationIndices

Matches the existing DCDM_DR_Species pattern: composite's register method
dispatches to its DR child (first) and DNCDM child (second) using the
per-child sub-layouts. The module-side manual split disappears in the
next commit; this commit just makes the override available.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Scalar registration → single dispatch loop (y-vector reorder lands here)

Replace the long per-species `if (all_species_.count("X"))` chain in `_scalars_` with one loop over `all_species_`. UR pre-loop layout setup moves into the species. **This is where the integration y-vector reorders for DNCDM_DR setups.**

**Files:**
- Modify: `source/perturbations_module.cpp` (lines 3338-3476 → replaced with the loop)
- Modify: `species/ultra_relativistic.h` and `species/ultra_relativistic.cpp` (read `ppr->l_max_ur` inside the register call)

- [ ] **Step 1: UR reads its own `ppr->l_max_ur`**

In `species/ultra_relativistic.cpp`, inside `UltraRelativisticSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base, …)`, add at the top (after the static_cast):

```cpp
  auto& layout = static_cast<PerturbLayout&>(base);
  layout.l_max = ppr->l_max_ur;   /* moved from perturb_vector_init pre-loop setup */
```

(if `layout.l_max` is already set elsewhere in this function, replace that assignment with the line above.)

- [ ] **Step 2: Replace the scalar registration chain**

In `source/perturbations_module.cpp`, in `perturb_vector_init`, replace lines 3338-3476 (everything from the `/* photons (...) */` comment up to the end of the `HasNcdm(...) { … }` block) with:

```cpp
    /* Per-species perturbation registration — single dispatch loop.
       Each species reads its own ppr->l_max_* values and fills its layout.
       Composites (IDM_DR_IDR, IDM_DRMD_IDR_DRMD, DCDM_DR, DNCDM_DR) register
       all children at their lex position in all_species_. */
    for (size_t i = 0; i < all_species_.size(); ++i) {
      all_species_[i]->RegisterPerturbationIndices(
          *ppv->species_layouts[i], ppv, ppr, index_pt, ppw, ppt->gauge);
    }

    /* perturbed_recombination — owned by the module (extension of the
       baryon-photon system, not a standalone species). Indices defined
       once TCA is off. */
    if ((ppt->has_perturbed_recombination == _TRUE_) &&
        (ppw->approx[ppw->index_ap_tca] == (int) tca_off)) {
      class_define_index(ppv->index_pt_perturbed_recombination_delta_temp, _TRUE_, index_pt, 1);
      class_define_index(ppv->index_pt_perturbed_recombination_delta_chi, _TRUE_, index_pt, 1);
    }
```

The metric perturbations block immediately after (`class_define_index(ppv->index_pt_eta, …)` and `class_define_index(ppv->index_pt_phi, …)`, originally at lines 3478-3487) stays untouched.

- [ ] **Step 3: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 4: Quick scenario sanity**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed` (the harness doesn't enforce tight Cl matching, so the reorder passes here).

- [ ] **Step 5: Tight Cl comparison vs. master baseline (the reorder verification)**

```bash
mkdir -p /tmp/perturb_cleanup_new
for ini in test/scenarios/gauge_lcdm.ini \
           test/scenarios/gauge_ncdm.ini \
           test/scenarios/gauge_dcdm.ini \
           test/scenarios/gauge_idmdr.ini \
           test/scenarios/gauge_fluid.ini \
           test/scenarios/gauge_scf.ini \
           test/scenarios/dncdm_dr.ini \
           test/scenarios/ncdm_dncdm_idmdr_combined.ini \
           test/scenarios/ncdm_multi_unsorted.ini \
           test/scenarios/ncdm_self_interacting.ini \
           test/scenarios/idm_dr_full.ini \
           test/scenarios/idm_drmd_full.ini ; do
  ./class "$ini" 2>&1 | tail -1
done
mv output/*.dat /tmp/perturb_cleanup_new/ 2>/dev/null || true
python test/scenarios/compare_tol.py test/scenarios/_baseline_master /tmp/perturb_cleanup_new "*.dat" 2>&1 | tail -20
```

Expected: every file marked `OK` with `worst_vs_colpeak < 1e-3`. If any file `FAIL`s, inspect the worst column. The reorder should not shift Cl spectra meaningfully — only the y-vector ordering changes, and the physics is the same. A genuine numerical shift > 0.1% means a logic bug (likely a forgotten ppr read in a species's register method).

- [ ] **Step 6: Commit**

```bash
git add source/perturbations_module.cpp species/ultra_relativistic.h species/ultra_relativistic.cpp
git commit -m "$(cat <<'EOF'
refactor: scalar perturb registration → single dispatch loop

Replaces the per-species if (all_species_.count("X")) chain in
perturb_vector_init's scalar block with one loop over all_species_.
UltraRelativisticSpecies reads ppr->l_max_ur inside its own register
method instead of having the module pre-populate layout.l_max.

The y-vector layout shifts for setups with DNCDM_DR composites
(per-composite DR+DNCDM children are now contiguous instead of being
split across the photon-side and NCDM-side blocks). Verified within
0.1% Cl tolerance via test/scenarios/compare_tol.py against master
baseline; project tolerance per
feedback_no_bit_identical_requirement.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Vector registration → single dispatch loop

Add `BaryonsSpecies::RegisterVectorPerturbationIndices` (currently inline in the module) and collapse the vector block to one loop. The vector/tensor base virtuals also need to gain a `const precision* ppr` argument so species can read their own `ppr->l_max_g_ten` etc. — Step 1 makes this signature change atomically.

**Files:**
- Modify: `species/base_species.h` (vector + tensor layout-form virtuals gain `const precision* ppr`)
- Modify: `species/photons.h`, `species/photons.cpp` (match new sig, fold the layout-l_max writes into the species body)
- Modify: `species/baryons.h`, `species/baryons.cpp` (new vector method)
- Modify: `source/perturbations_module.cpp` (current lines 3502-3520)

- [ ] **Step 1: Add `const precision* ppr` to the vector + tensor layout virtuals**

In `species/base_species.h`, change the two layout-form virtuals (lines 256-266):

```cpp
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
```

In `species/photons.h` and `species/photons.cpp`, change both Photons overrides to match the new signature (add `const precision* ppr` between `perturb_vector*` and `int& index_pt`). Inside each Photons override body, replace any `layout.l_max = …` / `layout.l_max_pol = …` setup with reads from ppr:

```cpp
  /* Photons read their own l_max values from ppr — no caller pre-setup. */
  layout.l_max     = ppr->l_max_g_ten;
  layout.l_max_pol = ppr->l_max_pol_g_ten;
```

(Both vector and tensor variants use `l_max_g_ten` and `l_max_pol_g_ten` — same as the current module pre-setup lines 3516-3517 and 3548-3549.)

The existing module call sites for vector (line 3518) and tensor (line 3550) must also get a `ppr` argument inserted (`sp->RegisterVectorPerturbationIndices(layout, ppv, ppr, index_pt, ppw, ppt->gauge);`) so the file builds at this step. The full loop-collapse happens in Step 4.

Build to check:

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build. Behavior unchanged so far.

- [ ] **Step 2: Declare `BaryonsSpecies::RegisterVectorPerturbationIndices`**

In `species/baryons.h`, in the public section near the other `Register*` declarations, add:

```cpp
  void RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                         perturb_vector* pv,
                                         const precision* ppr,
                                         int& index_pt,
                                         const perturb_workspace* ppw,
                                         int gauge) override;
```

- [ ] **Step 3: Implement it**

In `species/baryons.cpp`, add:

```cpp
void BaryonsSpecies::RegisterVectorPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                       perturb_vector* /*pv*/,
                                                       const precision* /*ppr*/,
                                                       int& index_pt,
                                                       const perturb_workspace* /*ppw*/,
                                                       int /*gauge*/) {
  auto& layout      = static_cast<PerturbLayout&>(base);
  layout.idx_theta  = index_pt++;   /* v_b^{(1)} */
}
```

- [ ] **Step 4: Replace the vector block with a dispatch loop**

In `source/perturbations_module.cpp`, in `perturb_vector_init`'s `if (_vectors_) { … }` block, replace lines 3502-3520 (everything from the baryon `idx_theta` inline write through the closing brace of the Photons block) with:

```cpp
    /* Per-species vector-mode registration. */
    for (size_t i = 0; i < all_species_.size(); ++i) {
      all_species_[i]->RegisterVectorPerturbationIndices(
          *ppv->species_layouts[i], ppv, ppr, index_pt, ppw, ppt->gauge);
    }
```

- [ ] **Step 5: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 6: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 7: Commit**

```bash
git add species/baryons.h species/baryons.cpp species/photons.h species/photons.cpp \
        species/base_species.h source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
refactor: vector perturb registration → single dispatch loop

Adds BaryonsSpecies::RegisterVectorPerturbationIndices for the v_b^(1)
slot (previously inline in perturb_vector_init). Vector + tensor
Register* base virtuals gain the precision* argument so each species
can read its own ppr->l_max_g_ten / ppr->l_max_pol_g_ten. Vector block
collapses to a single loop over all_species_.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Tensor registration → single dispatch loop (drop NCDM downcast)

Move the tensor NCDM pre-loop layout setup (`l_max = ppr->l_max_ncdm; q_size = nsp->q_size();`) into `NCDMSpecies::RegisterTensorPerturbationIndices` and drop the `dynamic_cast<NCDMSpecies*>` scan. The tensor UR hierarchy (`evolve_tensor_ur_`) stays module-owned.

**Files:**
- Modify: `species/ncdm_species.h` and `species/ncdm_species.cpp`
- Modify: `source/perturbations_module.cpp` (current lines 3569-3595)

- [ ] **Step 1: NCDM reads its own `ppr->l_max_ncdm` + `q_size`**

In `species/ncdm_species.cpp`, inside `NCDMSpecies::RegisterTensorPerturbationIndices(layout, …)`, at the top after the static_cast, add:

```cpp
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);
  layout.l_max  = ppr->l_max_ncdm;
  layout.q_size = q_size();
```

If `NCDMSpecies::RegisterTensorPerturbationIndices` does not currently exist, declare and implement it (likely it does — the tensor dispatch loop at module line 3586-3594 calls into a virtual). If only `RegisterTensorPerturbationIndices(perturb_vector*, …)` exists (the no-layout legacy now removed), add the layout-taking version with the assignments above plus the existing per-q tensor-slot reservation logic (lift it from the legacy implementation if present).

- [ ] **Step 2: Replace the tensor NCDM block with a clean loop**

In `source/perturbations_module.cpp`, in `perturb_vector_init`'s `if (_tensors_) { … if (evolve_tensor_ncdm_ == _TRUE_) { … } … }` block, replace the entire body of `if (evolve_tensor_ncdm_ == _TRUE_) { … }` (current lines 3577-3594, which contained both the `dynamic_cast<NCDMSpecies*>` pre-loop and the dispatch loop) with:

```cpp
    if (evolve_tensor_ncdm_ == _TRUE_) {
      for (size_t i = 0; i < all_species_.size(); ++i) {
        all_species_[i]->RegisterTensorPerturbationIndices(
            *ppv->species_layouts[i], ppv, ppr, index_pt, ppw, ppt->gauge);
      }
    }
```

The Photons tensor block immediately above (current lines 3543-3552) becomes part of this loop too — move the `layout.l_max = ppr->l_max_g_ten; layout.l_max_pol = ppr->l_max_pol_g_ten;` lines into `PhotonsSpecies::RegisterTensorPerturbationIndices` so the loop can subsume the manual call. Delete the manual Photons tensor block.

The `evolve_tensor_ur_` block (current lines 3554-3567) and the `gw` / `gwdot` `class_define_index` (current lines 3600-3604) stay untouched — they're module-owned, not species.

- [ ] **Step 3: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 4: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 5: Tight Cl comparison vs. master baseline (tensor channels)**

```bash
mkdir -p /tmp/perturb_cleanup_tensor
# any .ini that exercises a B-mode / tensor-NCDM path
for ini in test/scenarios/gauge_lcdm.ini \
           test/scenarios/gauge_ncdm.ini ; do
  ./class "$ini" 2>&1 | tail -1
done
mv output/*.dat /tmp/perturb_cleanup_tensor/ 2>/dev/null || true
python test/scenarios/compare_tol.py test/scenarios/_baseline_master /tmp/perturb_cleanup_tensor "*.dat" 2>&1 | tail -10
```

Expected: every file `OK`, `worst_vs_colpeak < 1e-3`.

- [ ] **Step 6: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp species/photons.cpp \
        source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
refactor: tensor perturb registration → single dispatch loop

NCDMSpecies::RegisterTensorPerturbationIndices now sets its own l_max
and q_size from ppr; PhotonsSpecies::RegisterTensorPerturbationIndices
reads l_max_g_ten / l_max_pol_g_ten. The module's tensor block reduces
to one loop over all_species_, dropping the dynamic_cast<NCDMSpecies*>
pre-loop. evolve_tensor_ur_ and gw/gwdot stay module-owned (non-species).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 8: Expand `MarkUsedInSources` signature + add Vector / Tensor variants

Interface-only change: pass `perturb_workspace*` to the scalar method, add `MarkVectorUsedInSources` and `MarkTensorUsedInSources` as no-op defaults. No behavior change yet — module still owns the approximation-conditional masking.

**Files:**
- Modify: `species/base_species.h`
- Modify: every species header that already overrides `MarkUsedInSources` (find via `grep -rln "MarkUsedInSources" species/`)

- [ ] **Step 1: Update `BaseSpecies` declarations**

In `species/base_species.h`, change the existing `MarkUsedInSources` virtual to:

```cpp
  /* Mark perturbation slots NOT needed for source evaluation as _FALSE_.
     ppw is supplied so species can consult their approximation flags. */
  virtual void MarkUsedInSources(const PerturbLayout& /*layout*/,
                                 const perturb_workspace* /*ppw*/,
                                 int* /*used_in_sources*/) const {}

  virtual void MarkVectorUsedInSources(const PerturbLayout& /*layout*/,
                                       const perturb_workspace* /*ppw*/,
                                       int* /*used_in_sources*/) const {}

  virtual void MarkTensorUsedInSources(const PerturbLayout& /*layout*/,
                                       const perturb_workspace* /*ppw*/,
                                       int* /*used_in_sources*/) const {}
```

- [ ] **Step 2: Update existing overrides to the new signature**

For every species file that currently overrides `MarkUsedInSources(const PerturbLayout&, int*)`, change the signature to add `const perturb_workspace* ppw` as the second argument. Identify them:

```bash
grep -rln "MarkUsedInSources" species/
```

For each match, update both the declaration and the body to the three-arg form. If a species's existing body uses no approximation flags, it can ignore the new `ppw` argument (mark as `/*ppw*/`).

- [ ] **Step 3: Update the module's existing call site**

In `source/perturbations_module.cpp`, change the existing dispatch loop (currently lines 3628-3630):

```cpp
  for (size_t i = 0; i < all_species_.size(); ++i) {
    all_species_[i]->MarkUsedInSources(*ppv->species_layouts[i], ppv->used_in_sources);
  }
```

to:

```cpp
  for (size_t i = 0; i < all_species_.size(); ++i) {
    all_species_[i]->MarkUsedInSources(
        *ppv->species_layouts[i], ppw, ppv->used_in_sources);
  }
```

- [ ] **Step 4: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 5: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed` (no behavior change — signature-only commit).

- [ ] **Step 6: Commit**

```bash
git add species/ source/
git commit -m "$(cat <<'EOF'
refactor: expand MarkUsedInSources, add Vector/Tensor variants

Scalar MarkUsedInSources now takes (layout, ppw, used_in_sources) so
species can consult approximation flags. Adds Mark{Vector,Tensor}
UsedInSources as no-op defaults, mirroring Register{,Vector,Tensor}
PerturbationIndices.

Interface-only commit: module still owns all approximation-conditional
source masks; species methods get bodies in the next commits.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 9: Move scalar source masks into species

Migrate the four approximation-conditional masking blocks in `perturb_vector_init`'s `if (_scalars_)` source-mask section into the owning species's `MarkUsedInSources` body.

**Files:**
- Modify: `species/photons.h`, `species/photons.cpp` (scalar mask)
- Modify: `species/ultra_relativistic.h`, `species/ultra_relativistic.cpp` (UR scalar mask)
- Modify: `species/idm_dr_idr_species.h`, `species/idm_dr_idr_species.cpp` (IDR scalar mask)
- Modify: `species/ncdm_species.h`, `species/ncdm_species.cpp` (NCDM scalar mask)
- Modify: `source/perturbations_module.cpp` (delete current lines 3636-3713)

- [ ] **Step 1: Photons scalar mask**

In `species/photons.cpp`, implement (and declare in `species/photons.h`):

```cpp
void PhotonsSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                       const perturb_workspace* ppw,
                                       int* used_in_sources) const {
  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  /* In scalar mode, photon l>=3 (temperature) and pol1, pol l>=3 are not
     needed in source evaluation when both rsa and tca are off. */
  if (ppw->approx[ppw->index_ap_rsa] != (int) rsa_off) return;
  if (ppw->approx[ppw->index_ap_tca] != (int) tca_off) return;

  for (int idx = g_lay.idx_l3; idx <= g_lay.idx_delta + g_lay.l_max; ++idx)
    used_in_sources[idx] = _FALSE_;

  used_in_sources[g_lay.idx_pol1] = _FALSE_;
  for (int idx = g_lay.idx_pol3; idx <= g_lay.idx_pol0 + g_lay.l_max_pol; ++idx)
    used_in_sources[idx] = _FALSE_;
}
```

- [ ] **Step 2: UR scalar mask**

In `species/ultra_relativistic.cpp`, implement (and declare):

```cpp
void UltraRelativisticSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                                 const perturb_workspace* ppw,
                                                 int* used_in_sources) const {
  const auto& ur_lay = static_cast<const PerturbLayout&>(base);
  /* UR l>=3 not needed in sources when both rsa and ufa are off. */
  if (ppw->approx[ppw->index_ap_rsa] != (int) rsa_off) return;
  if (ppw->approx[ppw->index_ap_ufa] != (int) ufa_off) return;

  for (int idx = ur_lay.idx_l3; idx <= ur_lay.idx_delta + ur_lay.l_max; ++idx)
    used_in_sources[idx] = _FALSE_;
}
```

- [ ] **Step 3: IDM_DR_IDR scalar mask**

The composite needs to read `ppt->idr_nature`. It inherits from `CompositeSpecies` (which already forwards `SetPerturbs` to children) but doesn't currently store its own `ppt_`. Override `SetPerturbs` to capture it while still calling the base for child forwarding.

In `species/idm_dr_idr_species.h`, in the public section, add:

```cpp
  void SetPerturbs(const perturbs* ppt) override {
    ppt_ = ppt;
    CompositeSpecies::SetPerturbs(ppt);
  }

  void MarkUsedInSources(const BaseSpecies::PerturbLayout& layout,
                         const perturb_workspace* ppw,
                         int* used_in_sources) const override;
```

In the private section:

```cpp
  const perturbs* ppt_ = nullptr;
```

In `species/idm_dr_idr_species.cpp`, implement:

```cpp
void IDM_DR_IDR_Species::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                           const perturb_workspace* ppw,
                                           int* used_in_sources) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  /* IDR l>=3 not needed in sources when rsa_idr is off, idr is
     free-streaming, and tca_idm_dr is off. */
  if (ppw->approx[ppw->index_ap_rsa_idr] != (int) rsa_idr_off) return;
  if (ppt_->idr_nature != idr_free_streaming) return;
  if (ppw->approx[ppw->index_ap_tca_idm_dr] != (int) tca_idm_dr_off) return;

  for (int idx = my.idr.idx_l3; idx <= my.idr.idx_delta + my.idr.l_max; ++idx)
    used_in_sources[idx] = _FALSE_;
}
```

- [ ] **Step 4: NCDM scalar mask**

In `species/ncdm_species.cpp`, implement (and declare). NCDM masks per-q `l>2`:

```cpp
void NCDMSpecies::MarkUsedInSources(const BaseSpecies::PerturbLayout& base,
                                    const perturb_workspace* ppw,
                                    int* used_in_sources) const {
  const auto& ncdm_lay = static_cast<const NCDMBaseSpecies::PerturbLayout&>(base);
  /* NCDM multipoles l>2 not needed in sources when ncdmfa is off. */
  if (ppw->approx[ppw->index_ap_ncdmfa] != (int) ncdmfa_off) return;

  for (int q = 0; q < ncdm_lay.q_size; ++q) {
    int idx = ncdm_lay.index_per_q[q];
    for (int l = 0; l <= ncdm_lay.l_max; ++l) {
      if (l > 2) used_in_sources[idx] = _FALSE_;
      ++idx;
    }
  }
}
```

- [ ] **Step 5: Delete the module-side bespoke masks**

In `source/perturbations_module.cpp`, in `perturb_vector_init`'s `if (_scalars_) { … }` source-mask section, delete the four blocks (current lines 3636-3713):

- The Photons scalar masking block (lines 3636-3656)
- The UR scalar masking block (lines 3658-3670)
- The IDM_DR_IDR scalar masking block (lines 3672-3692)
- The NCDM scalar masking block (lines 3694-3713)

The species's `MarkUsedInSources` is already called by the existing dispatch loop (Task 8's surviving loop at line 3628).

- [ ] **Step 6: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 7: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 8: Tight Cl comparison**

```bash
mkdir -p /tmp/perturb_cleanup_mask
for ini in test/scenarios/gauge_lcdm.ini \
           test/scenarios/gauge_ncdm.ini \
           test/scenarios/gauge_idmdr.ini \
           test/scenarios/dncdm_dr.ini ; do
  ./class "$ini" 2>&1 | tail -1
done
mv output/*.dat /tmp/perturb_cleanup_mask/ 2>/dev/null || true
python test/scenarios/compare_tol.py test/scenarios/_baseline_master /tmp/perturb_cleanup_mask "*.dat" 2>&1 | tail -10
```

Expected: `OK` with `worst_vs_colpeak < 1e-3`. (Masking only affects sampling-of-source efficiency, not physics, so any drift here is concerning.)

- [ ] **Step 9: Commit**

```bash
git add species/ source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
refactor: scalar source masks → per-species MarkUsedInSources

PhotonsSpecies, UltraRelativisticSpecies, IDM_DR_IDR_Species, and
NCDMSpecies each own their approximation-conditional source masking
for scalar mode. Module-side bespoke blocks in perturb_vector_init
gone; the existing per-species dispatch loop handles everything.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 10: Move tensor source mask into species + add tensor mark dispatch loop

The only tensor-mode approximation-conditional mask is Photons-specific (current lines 3722-3741). Move it into `PhotonsSpecies::MarkTensorUsedInSources` and add the matching dispatch loop in the module.

**Files:**
- Modify: `species/photons.h`, `species/photons.cpp`
- Modify: `source/perturbations_module.cpp` (delete current lines 3716-3743, then add a tensor dispatch loop)

- [ ] **Step 1: Photons tensor mask**

In `species/photons.cpp`, implement (and declare in `photons.h`):

```cpp
void PhotonsSpecies::MarkTensorUsedInSources(const BaseSpecies::PerturbLayout& base,
                                             const perturb_workspace* ppw,
                                             int* used_in_sources) const {
  const auto& g_lay = static_cast<const PerturbLayout&>(base);
  /* In tensor mode, we only need temperature l=0,2,4 and pol l=0,2,4
     when both rsa and tca are off. */
  if (ppw->approx[ppw->index_ap_rsa] != (int) rsa_off) return;
  if (ppw->approx[ppw->index_ap_tca] != (int) tca_off) return;

  used_in_sources[g_lay.idx_theta] = _FALSE_;
  used_in_sources[g_lay.idx_l3]    = _FALSE_;
  for (int idx = g_lay.idx_delta + 5; idx <= g_lay.idx_delta + g_lay.l_max; ++idx)
    used_in_sources[idx] = _FALSE_;

  used_in_sources[g_lay.idx_pol1] = _FALSE_;
  used_in_sources[g_lay.idx_pol3] = _FALSE_;
  for (int idx = g_lay.idx_pol0 + 5; idx <= g_lay.idx_pol0 + g_lay.l_max_pol; ++idx)
    used_in_sources[idx] = _FALSE_;
}
```

- [ ] **Step 2: Replace the module tensor mask block with a dispatch loop**

In `source/perturbations_module.cpp`, in `perturb_vector_init`'s `if (_tensors_) { … }` source-mask section, replace the bespoke Photons tensor block (current lines 3716-3743, the entire `if (ppw->approx[ppw->index_ap_rsa] == (int) rsa_off)` nested block) with:

```cpp
  if (_tensors_) {
    for (size_t i = 0; i < all_species_.size(); ++i) {
      all_species_[i]->MarkTensorUsedInSources(
          *ppv->species_layouts[i], ppw, ppv->used_in_sources);
    }
    /* gw is a metric slot, not a species — module-owned. */
    ppv->used_in_sources[ppv->index_pt_gw] = _FALSE_;
  }
```

(Keep the existing `gw` line; just relocate it inside the simplified block.)

- [ ] **Step 3: Build**

```bash
make -j8 class 2>&1 | tail -3
```

Expected: clean build.

- [ ] **Step 4: Regression**

```bash
pip install . 2>&1 | tail -3 && \
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -3
```

Expected: `84 passed`.

- [ ] **Step 5: Commit**

```bash
git add species/photons.h species/photons.cpp source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
refactor: tensor source mask → PhotonsSpecies::MarkTensorUsedInSources

Photons owns its tensor-mode mask (theta, l3, l5..lmax, pol1, pol3,
pol5..lmax_pol when rsa_off + tca_off). Module's tensor source-mask
block becomes a single dispatch loop over all_species_ plus the
standalone gw line (metric slot, not a species).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 11: Final regression at TEST_LEVEL=3

Run the full 1260-scenario suite plus a full Cl comparison vs. master before declaring done.

**Files:** none modified.

- [ ] **Step 1: Full pytest regression**

```bash
TEST_LEVEL=3 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -5
```

Expected: `1260 passed`.

- [ ] **Step 2: Full Cl comparison vs. master baseline (every captured `.ini`)**

```bash
mkdir -p /tmp/perturb_cleanup_final
for ini in test/scenarios/gauge_lcdm.ini \
           test/scenarios/gauge_ncdm.ini \
           test/scenarios/gauge_dcdm.ini \
           test/scenarios/gauge_idmdr.ini \
           test/scenarios/gauge_fluid.ini \
           test/scenarios/gauge_scf.ini \
           test/scenarios/dncdm_dr.ini \
           test/scenarios/ncdm_dncdm_idmdr_combined.ini \
           test/scenarios/ncdm_multi_unsorted.ini \
           test/scenarios/ncdm_self_interacting.ini \
           test/scenarios/idm_dr_full.ini \
           test/scenarios/idm_drmd_full.ini ; do
  ./class "$ini" 2>&1 | tail -1
done
mv output/*.dat /tmp/perturb_cleanup_final/ 2>/dev/null || true
python test/scenarios/compare_tol.py test/scenarios/_baseline_master /tmp/perturb_cleanup_final "*.dat"
```

Expected: every file `OK`, `worst_vs_colpeak < 1e-3`.

- [ ] **Step 3: Clean up the local baseline dir** (it was never committed)

```bash
rm -rf test/scenarios/_baseline_master /tmp/perturb_cleanup_*
```

- [ ] **Step 4: Survey the final diff**

```bash
git log --oneline 7616ffe9..HEAD
git diff --stat 7616ffe9..HEAD
```

Expected: ~10 commits, mostly subtractive (lines removed > added). The big remover is `source/perturbations_module.cpp`.

- [ ] **Step 5: No commit** — verification only.
