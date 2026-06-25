# `has_ncdm` Collection-Property Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-call `HasNcdm(all_species_)` O(n) double-`dynamic_cast` scan (19 sites in `perturbations_module`) and the now-redundant `ResolvedSpecies::has_ncdm` cache with a single boolean computed once in `SpeciesCollection::freeze()` and read via `all_species_.has_ncdm()`.

**Architecture:** `SpeciesCollection` is frozen after construction and already caches derived hot-path state in `freeze()` (`photons_`/`baryons_` pointers + indices). Add a `bool has_ncdm_`, computed in the existing `freeze()` index loop with the exact predicate used today (`dynamic_cast<NCDMBaseSpecies*> || dynamic_cast<DNCDM_DR_Species*>`), exposed via `bool has_ncdm() const`. Migrate every `perturbations_module` caller to it, delete the `HasNcdm` free function, and remove `ResolvedSpecies::has_ncdm` (its three per-step readers point at `all_species_.has_ncdm()`, itself an inlinable cached-bool read).

**Tech Stack:** C++17, CMake (build tree `build_prof/`, binaries land in repo root as `./class`), `python/test_class.py` regression suite, `test/scenarios/compare_tol.py` scale-relative comparator, clang-format 22.1.3.

## Global Constraints

- **Verification gate: scale-relative ≤0.1%, NOT byte-identical.** Under fast-math-by-default (#338), logically-identical changes can drift at ULP level from codegen reassociation, so byte-identity is not a reliable gate in this repo (per `feedback_no_bit_identical_requirement` and the fluid-ppf / PR #341 experience). Every code task reruns the baseline scenarios and compares against the Task-0 master baseline with `python3 test/scenarios/compare_tol.py /tmp/has_ncdm_golden <new_dir> '*.dat'` (RTOL=1e-3, zero-crossing-aware): **every file must report `OK`** (n_exceed=0). `all_species_.has_ncdm()` returns the *same boolean* as `HasNcdm(all_species_)`, so differences should be zero or ULP-level; a real `FAIL`/`MISSING`/`SHAPE` is a logic error (e.g. a typo in the moved predicate) — STOP and investigate.
- **Final correctness gate:** `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` full `python/test_class.py` suite, 0 failures vs master classyref. This is the comprehensive gate covering the DNCDM_DR-composite branch and both gauges (the fast compare_tol scenarios do not exercise DNCDM_DR perturbations).
- **Out of scope — do NOT change `source/output_module.cpp:668`.** It uses the same per-species `dynamic_cast` pair but as a `count` feeding `class_test(n_ncdm_family > 1, ...)` (CAMB-format single-ncdm constraint). That is a `> 1` check, not presence; a boolean `has_ncdm()` cannot express it. Leave it untouched.
- **Leave `perturbations_module.cpp` ncdm `#include`s in place** — `NCDMBaseSpecies`/`DNCDM_DR_Species` are still used elsewhere in the file (e.g. lines 575, 5006, 5627).
- **Never `git add -A`** in this repo (build artifacts leak in); always `git add` explicit paths.
- **Subagents work in the shared current working tree on branch `has-ncdm-collection-property`.** Do NOT `git checkout`/`switch`/`reset`/`branch` (`feedback_subagents_no_git_checkout`). Only edit, build, `git add <explicit paths>`, `git commit`.
- **clang-format 22.1.3** every touched file before the final commit.
- No perf benchmark is required — this is a pure dedup; the per-step perf win was already captured by the resolved-species-handles programme.
- Line numbers below are as of the current branch base and **will shift** as edits land — match on the quoted code, not the line number.

---

### Task 0: Capture the master golden baselines

**Files:** none modified (baseline capture only).

**Interfaces:**
- Produces: rootless scenario inis `/tmp/has_ncdm_base.ini` + `/tmp/has_ncdm_lcdm.ini` and baseline outputs in `/tmp/has_ncdm_golden/`, used by every later task's compare_tol gate. The current branch HEAD (spec + plan files only, no code change) is code-identical to master, so building it produces master output.

- [ ] **Step 1: Build the current tree (sanity + baseline binary)**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: builds `./class` with no errors (last line a link step / `Built target class`).

- [ ] **Step 2: Prepare rootless scenario inis (reused by every task)**

`base_2018...ini` has `N_ncdm = 1` (true path; full TT + pk). `gauge_lcdm.ini` has no ncdm (false path). Strip any `root=` so each run can pin its own output dir:
```bash
grep -vE '^\s*root\s*=' base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini > /tmp/has_ncdm_base.ini
grep -vE '^\s*root\s*=' test/scenarios/gauge_lcdm.ini                  > /tmp/has_ncdm_lcdm.ini
```

- [ ] **Step 3: Capture the master baseline outputs**

```bash
GOLD=/tmp/has_ncdm_golden; rm -rf "$GOLD"; mkdir -p "$GOLD"
{ cat /tmp/has_ncdm_base.ini; printf '\nroot = %s/base_\n' "$GOLD"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
{ cat /tmp/has_ncdm_lcdm.ini; printf '\nroot = %s/lcdm_\n' "$GOLD"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
ls -l "$GOLD"/*.dat
```
Expected: non-empty baseline files exist — `base_cl_lensed.dat`, `base_pk.dat`, `lcdm_cl.dat`, `lcdm_pk.dat` (and possibly additional `base_*`/`lcdm_*` `.dat`).

- [ ] **Step 4: DNCDM_DR composite smoke (construction-path coverage)**

`dncdm_dr.ini` builds a `DNCDM_DR_Species` composite (the second `dynamic_cast` branch in `freeze()`) but requests no Cl/Pk output. Confirm it runs cleanly now so any later breakage is attributable.
```bash
./class test/scenarios/dncdm_dr.ini >/tmp/has_ncdm_dncdm_dr.log 2>&1; echo "exit=$?"
```
Expected: `exit=0`.

- [ ] **Step 5: No commit** (nothing changed).

---

### Task 1: Add `SpeciesCollection::has_ncdm()` cached in `freeze()`

**Files:**
- Modify: `species/species_collection.h` (private member + accessor)
- Modify: `species/species_collection.cpp` (includes + `freeze()` computation)

**Interfaces:**
- Produces, for Task 2: `bool SpeciesCollection::has_ncdm() const` — true iff any species is an `NCDMBaseSpecies` or a `DNCDM_DR_Species`. Valid only after `freeze()` (asserted).

- [ ] **Step 1: Add the accessor to the header**

In `species/species_collection.h`, immediately after the `baryons_index()` accessor:
```cpp
  std::size_t baryons_index() const {
    assert(frozen_);
    return baryons_index_;
  }
```
add:
```cpp
  /** True iff the collection holds any NCDM-family species (an NCDMBaseSpecies,
   *  or the DNCDM_DR_Species composite). Valid only after freeze(). */
  bool has_ncdm() const {
    assert(frozen_);
    return has_ncdm_;
  }
```

- [ ] **Step 2: Add the private member**

In `species/species_collection.h`, in the private members, after:
```cpp
  std::size_t baryons_index_ = 0;
```
add:
```cpp
  bool has_ncdm_ = false;
```

- [ ] **Step 3: Add the ncdm includes to the .cpp**

In `species/species_collection.cpp`, replace the include block:
```cpp
#include "species_collection.h"

#include <algorithm>
#include <stdexcept>
```
with:
```cpp
#include "species_collection.h"

#include <algorithm>
#include <stdexcept>

#include "dncdm_dr_species.h"
#include "ncdm_base_species.h"
```

- [ ] **Step 4: Compute `has_ncdm_` in the existing `freeze()` loop**

In `species/species_collection.cpp`, in `freeze()`, replace the photons/baryons index loop:
```cpp
  for (std::size_t i = 0; i < species_.size(); ++i) {
    if (species_[i].key == "Photons")
      photons_index_ = i;
    if (species_[i].key == "Baryons")
      baryons_index_ = i;
  }
```
with:
```cpp
  for (std::size_t i = 0; i < species_.size(); ++i) {
    if (species_[i].key == "Photons")
      photons_index_ = i;
    if (species_[i].key == "Baryons")
      baryons_index_ = i;
    if (dynamic_cast<NCDMBaseSpecies*>(species_[i].get()) ||
        dynamic_cast<DNCDM_DR_Species*>(species_[i].get()))
      has_ncdm_ = true;
  }
```

- [ ] **Step 5: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: clean build. (`has_ncdm_` is populated but not yet read — fine; members do not warn when unused.)

- [ ] **Step 6: Verify ≤0.1% (accessor unused → output must match baseline)**

```bash
NEW=/tmp/has_ncdm_new; rm -rf "$NEW"; mkdir -p "$NEW"
{ cat /tmp/has_ncdm_base.ini; printf '\nroot = %s/base_\n' "$NEW"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
{ cat /tmp/has_ncdm_lcdm.ini; printf '\nroot = %s/lcdm_\n' "$NEW"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
python3 test/scenarios/compare_tol.py /tmp/has_ncdm_golden "$NEW" '*.dat'
```
Expected: every line `OK` (n_exceed=0); no `FAIL`/`MISSING`/`SHAPE`.

- [ ] **Step 7: Commit**

```bash
git add species/species_collection.h species/species_collection.cpp
git commit -m "species: cache has_ncdm() on SpeciesCollection (computed in freeze)"
```

---

### Task 2: Migrate `perturbations_module` callers; delete `HasNcdm` + `ResolvedSpecies::has_ncdm`

**Files:**
- Modify: `source/perturbations_module.h` (remove `ResolvedSpecies::has_ncdm` field)
- Modify: `source/perturbations_module.cpp` (remove free function + `resolved_` assignment; repoint readers; migrate 19 call sites)

**Interfaces:**
- Consumes: `bool SpeciesCollection::has_ncdm() const` from Task 1.

**Ordering matters** — do the deletions before the global replacements so no edit lands on a line about to be deleted.

- [ ] **Step 1: Remove the `has_ncdm` field from `ResolvedSpecies`**

In `source/perturbations_module.h`, delete the line:
```cpp
  bool has_ncdm                             = false;            // ncdmfa scheme present
```

- [ ] **Step 2: Remove the `resolved_.has_ncdm` assignment in `ResolveSpecies()`**

In `source/perturbations_module.cpp`, delete the line:
```cpp
  resolved_.has_ncdm  = HasNcdm(all_species_);
```

- [ ] **Step 3: Repoint the three per-step readers**

In `source/perturbations_module.cpp`, run a global replace of `resolved_.has_ncdm` → `all_species_.has_ncdm()`. This updates exactly the three remaining readers (in `perturb_approximations`, `perturb_timescale_member`, `perturb_total_stress_energy`), e.g.:
```cpp
    if (resolved_.has_ncdm) {
```
→
```cpp
    if (all_species_.has_ncdm()) {
```
and the `||`-guard in `perturb_timescale_member`:
```cpp
    if ((ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) || (resolved_.has_ncdm))
```
→
```cpp
    if ((ppw->approx[ppw->index_ap_rsa] == (int) rsa_off) || (all_species_.has_ncdm()))
```
Then confirm none remain:
```bash
grep -n 'resolved_.has_ncdm' source/perturbations_module.cpp
```
Expected: no output.

- [ ] **Step 4: Delete the `HasNcdm` free function and its anonymous namespace**

In `source/perturbations_module.cpp`, delete the entire block (it is the only member of the anonymous namespace):
```cpp
namespace {

/** True iff any NCDM-family species (NCDM, NCDMInteracting, or the DNCDMSpecies
 *  child of a DNCDM_DR_Species composite) is present in all_species_. */
bool HasNcdm(const SpeciesCollection& all_species) {
  for (const auto& sp : all_species) {
    if (dynamic_cast<NCDMBaseSpecies*>(sp.get()) || dynamic_cast<DNCDM_DR_Species*>(sp.get()))
      return true;
  }
  return false;
}

}  // namespace
```

- [ ] **Step 5: Migrate the 19 remaining call sites**

In `source/perturbations_module.cpp`, run a global replace of `HasNcdm(all_species_)` → `all_species_.has_ncdm()`. These are all boolean presence checks; the surrounding expressions are unchanged. Representative sites (only the token changes):
```cpp
  if (HasNcdm(all_species_)) {
```
→
```cpp
  if (all_species_.has_ncdm()) {
```
```cpp
        if ((all_species_.count("UR")) || HasNcdm(all_species_))
```
→
```cpp
        if ((all_species_.count("UR")) || all_species_.has_ncdm())
```
```cpp
    class_define_index(ppw->index_ap_ncdmfa, HasNcdm(all_species_), index_ap, 1);
```
→
```cpp
    class_define_index(ppw->index_ap_ncdmfa, all_species_.has_ncdm(), index_ap, 1);
```

- [ ] **Step 6: Confirm no `HasNcdm` reference remains**

```bash
grep -n 'HasNcdm' source/perturbations_module.cpp
```
Expected: no output.

- [ ] **Step 7: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: clean build.

- [ ] **Step 8: Verify ≤0.1% on both paths**

```bash
NEW=/tmp/has_ncdm_new; rm -rf "$NEW"; mkdir -p "$NEW"
{ cat /tmp/has_ncdm_base.ini; printf '\nroot = %s/base_\n' "$NEW"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
{ cat /tmp/has_ncdm_lcdm.ini; printf '\nroot = %s/lcdm_\n' "$NEW"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
python3 test/scenarios/compare_tol.py /tmp/has_ncdm_golden "$NEW" '*.dat'
```
Expected: every line `OK` (n_exceed=0).

- [ ] **Step 9: DNCDM_DR composite smoke (second predicate branch still detected)**

```bash
./class test/scenarios/dncdm_dr.ini >/tmp/has_ncdm_dncdm_dr_after.log 2>&1; echo "exit=$?"
```
Expected: `exit=0` (matches the Task 0 smoke).

- [ ] **Step 10: Commit**

```bash
git add source/perturbations_module.h source/perturbations_module.cpp
git commit -m "perturb: use all_species_.has_ncdm(); drop HasNcdm + ResolvedSpecies::has_ncdm"
```

---

### Task 3: Format, full regression gate, finalize

**Files:**
- Modify: any touched file (formatting only)

- [ ] **Step 1: clang-format every touched file**

```bash
clang-format --version   # expect 22.1.3
clang-format -i species/species_collection.h species/species_collection.cpp \
  source/perturbations_module.h source/perturbations_module.cpp
git diff --stat
```
Expected: only whitespace/formatting changes, if any.

- [ ] **Step 2: Rebuild and re-confirm ≤0.1%**

```bash
cmake --build build_prof -j 2>&1 | tail -3
NEW=/tmp/has_ncdm_new; rm -rf "$NEW"; mkdir -p "$NEW"
{ cat /tmp/has_ncdm_base.ini; printf '\nroot = %s/base_\n' "$NEW"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
{ cat /tmp/has_ncdm_lcdm.ini; printf '\nroot = %s/lcdm_\n' "$NEW"; } > /tmp/has_ncdm_run.ini; ./class /tmp/has_ncdm_run.ini >/dev/null
python3 test/scenarios/compare_tol.py /tmp/has_ncdm_golden "$NEW" '*.dat'
```
Expected: clean build; every line `OK`.

- [ ] **Step 3: Full regression suite vs master classyref (comprehensive gate)**

This exercises the DNCDM_DR-composite perturbation path, both gauges, and ncdm scenarios the fast compare_tol checks do not cover.
```bash
TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python3 -m pytest python/test_class.py -q 2>&1 | tail -20
```
Expected: all pass, 0 failures vs the current master classyref. If classyref is stale relative to master, regenerate it per the project's classyref procedure, then re-run.

- [ ] **Step 4: Commit formatting (if any)**

```bash
git add species/species_collection.h species/species_collection.cpp \
  source/perturbations_module.h source/perturbations_module.cpp
git commit -m "species/perturb: clang-format has_ncdm property changes" || echo "nothing to format-commit"
```

- [ ] **Step 5: Open the PR**

Summarize: `has_ncdm` now a cached property of `SpeciesCollection` (computed in `freeze()`), 19 `perturbations_module` call sites migrated, `HasNcdm` free function and `ResolvedSpecies::has_ncdm` removed, ≤0.1% on base_2018 + gauge_lcdm and green on the `TEST_LEVEL=2` suite. Note `output_module.cpp`'s `n_ncdm_family > 1` count check was intentionally left unchanged. Link the spec `docs/superpowers/specs/2026-06-25-has-ncdm-collection-property-design.md`.

---

## Notes for the implementer

- The moved `freeze()` predicate is the *exact* two-`dynamic_cast` test from the old `HasNcdm`; a real compare_tol `FAIL` means a transcription error — debug before continuing (see `superpowers:systematic-debugging`).
- The gate is scale-relative ≤0.1% (`compare_tol.py`), NOT byte-identical — ULP-level drift under ffast-math is expected and passes; only a meaningfully-off column fails.
- Do not touch `source/output_module.cpp` — its `n_ncdm_family > 1` count is a different question (Global Constraints).
- Do not remove the ncdm `#include`s from `perturbations_module.cpp`; `NCDMBaseSpecies`/`DNCDM_DR_Species` remain in use there.
