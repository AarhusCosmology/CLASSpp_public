# Unified Species `CreateAll` Loop and Deferred Omega Closure — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the current mid-parser temp-NCDM-construction + hand-written `if (has_X)` block in `InputModule::ConstructSpecies` with a uniform `CreateAll` loop driven by a factory list owned by the species subsystem, and defer Omega closure to after the loop.

**Architecture:** Every species class exposes a static `CreateAll(const SpeciesBuildContext&) -> std::vector<Named>`. A constexpr factory list in `species/all_species.h` is the single source of truth for the species roster. `ConstructSpecies` iterates that list, sums `GetOmega0()` over the produced sectors, computes the closure value if any (`pba->Omega0_lambda = 1 − Omega_k − Σ`, or fld / scf), and runs the closure species' `CreateAll` last. `input_read_parameters` only *detects* which species is the closure; it no longer builds species or computes the closure value.

**Tech Stack:** C++17, `make class -j` for build, pytest regression suite in `python/test_class.py` for verification. No unit-test framework for C++ — verification is bit-identical reference comparison at the end.

**Spec:** [docs/superpowers/specs/2026-04-25-unified-species-createall-design.md](../specs/2026-04-25-unified-species-createall-design.md)

---

## Conventions

- After every code change: `make class -j 2>&1 | tail -10` and `./class explanatory.ini 2>&1 | tail -3`. Both must succeed.
- Fast smoke at checkpoints: `./class explanatory.ini 2>&1 | tail -3`.
- Reference gate (used at the major checkpoints in Tasks 19, 21, 23, and 26): `cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30`.
- Every new `.cpp`/`.h` must land in **all three** build systems — `Makefile`, `setup.py`, **and** `CLASS.xcodeproj/project.pbxproj`. The user maintains Xcode manually but Make and setup.py must work. Two new headers in this plan: `species/species_build_context.h` and `species/all_species.h` (the latter replaces the existing `species/all_species_list.h`).
- Commits follow the existing repo style (short prefix like `species:`, `input:`, `background:` matching the touched area). Use `git -c commit.gpgsign=false commit ...` if gpg signing is not set up.
- Existing nested per-class `Named` structs (currently in `NCDMSpecies`, `DNCDMSpecies`, `DNCDM_DR_Species`, `NCDMInteractingSpecies`) get replaced by a single shared `Named` in `species/species_build_context.h`. Only `DNCDMSpecies` keeps a concrete-typed return (an internal helper for `DNCDM_DR_Species`, not on the factory list).

---

## Task 1: Create branch and capture baseline reference outputs

**Files:** none (git + reference capture).

- [ ] **Step 1: Create branch**

```bash
cd /Users/au192734/Projects/class_claude
git checkout -b unified-species-createall master
```

- [ ] **Step 2: Baseline check**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build; `All parameters and all quantities computed successfully`. If not, master is broken — stop and investigate.

- [ ] **Step 3: Capture baseline reference outputs**

```bash
mkdir -p /tmp/baseline_pre
cp -r output/. /tmp/baseline_pre/explanatory_dump 2>/dev/null || true
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

Expected: reference suite passes on master before any code change. If anything fails on master, fix or note before proceeding — we need a green baseline to diff against.

- [ ] **Step 4: Commit (no code changes; just a marker)**

Skip — no commit needed, this is verification only.

---

## Task 2: Create `species/species_build_context.h`

**Files:**
- Create: `species/species_build_context.h`

- [ ] **Step 1: Write the new header**

Create `species/species_build_context.h`:

```cpp
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "base_species.h"

class FileContent;
struct background;
struct precision;
class BackgroundModule;
struct NcdmSettings;

/**
 * Inputs every species' static CreateAll factory needs.
 * Bundled into one struct to keep factory signatures uniform.
 */
struct SpeciesBuildContext {
  FileContent*            pfc;
  const background*       pba;
  const precision*        ppr;
  const NcdmSettings&     ncdm_settings;
  const BackgroundModule* bgm;   // nullptr at species-construction time
};

/**
 * One species sector produced by a CreateAll factory.
 * `key` is the SpeciesCollection insertion key.
 */
struct Named {
  std::string                  key;
  std::unique_ptr<BaseSpecies> species;
};
```

- [ ] **Step 2: Forward-declare `NcdmSettings` correctly**

`NcdmSettings` is currently a `struct` defined in `species/ncdm_base_species.h`. Confirm that the forward declaration above (`struct NcdmSettings`) matches by checking the file:

```bash
grep -n "struct NcdmSettings\|class NcdmSettings" species/ncdm_base_species.h
```

If the result shows `struct NcdmSettings`, leave the forward declaration as-is. If it shows `class NcdmSettings`, change the forward declaration to `class`.

- [ ] **Step 3: Add to all three build systems**

In `Makefile`, find the species headers line and add `species/species_build_context.h`. Most species headers in the Makefile are in the `EXTERNAL` headers list — locate it:

```bash
grep -n "species_build_context\|all_species_list\|species/" Makefile | head -20
```

If species headers are not explicitly listed in Makefile, add nothing (Makefile uses dependency generation). Otherwise add the entry.

In `setup.py`, locate the species-headers list:

```bash
grep -n "species/" setup.py | head -20
```

Add `'species/species_build_context.h'` if there's a header list; otherwise no change.

In `CLASS.xcodeproj/project.pbxproj`, the user will add this manually after the rest of the plan converges. The plan does not touch the pbxproj for new headers; flag in the final commit message that the pbxproj needs the new header.

- [ ] **Step 4: Build check**

```bash
make class -j 2>&1 | tail -5
```

Expected: clean build (the new header is not yet included anywhere).

- [ ] **Step 5: Commit**

```bash
git add species/species_build_context.h Makefile setup.py
git commit -m "species: add SpeciesBuildContext and shared Named struct"
```

---

## Task 3: Add `ClosureSpecies` enum to `background`

**Files:**
- Modify: `source/background.h`

- [ ] **Step 1: Find the right insertion point**

```bash
grep -n "Omega0_lambda\|Omega0_fld\|Omega0_scf" source/background.h | head -5
```

Note the line numbers near `Omega0_lambda` etc.

- [ ] **Step 2: Add the enum and field**

In `source/background.h`, near the `Omega0_lambda` / `Omega0_fld` / `Omega0_scf` block, add:

```cpp
  /** Which species (if any) is the budget-closure species — set by
   *  input_read_parameters, consumed by InputModule::ConstructSpecies. */
  enum class ClosureSpecies { None, Lambda, Fluid, ScalarField };
  ClosureSpecies closure_species = ClosureSpecies::None;
```

If `background` is a plain `struct` (it currently is), the in-class default initializer is fine in C++17.

- [ ] **Step 3: Build check**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, explanatory still runs.

- [ ] **Step 4: Commit**

```bash
git add source/background.h
git commit -m "background: add closure_species enum field"
```

---

## Task 4: Replace `species/all_species_list.h` with `species/all_species.h` (types-only, no factory list yet)

**Files:**
- Create: `species/all_species.h`
- Delete: `species/all_species_list.h`
- Modify: `source/input_module.cpp` (one include line)

- [ ] **Step 1: Create `species/all_species.h`**

Initial content (factory list filled in later in Task 19):

```cpp
#pragma once
/**
 * Aggregator: include every species class, and (added in Task 19) declare
 * the constexpr factory list iterated by InputModule::ConstructSpecies.
 *
 * Adding a new species: write the class, include its header here, add one
 * row to kAllSpeciesFactories below.
 */
#include "baryons.h"
#include "base_species.h"
#include "cdm.h"
#include "dark_radiation_species.h"
#include "dcdm.h"
#include "dcdm_dr_species.h"
#include "dncdm_decay_radiation_species.h"
#include "dncdm_dr_species.h"
#include "fluid.h"
#include "idm_dr.h"
#include "idm_dr_idr_species.h"
#include "idm_drmd.h"
#include "idm_drmd_idr_drmd_species.h"
#include "idr.h"
#include "idr_drmd.h"
#include "lambda.h"
#include "ncdm_interacting_species.h"
#include "ncdm_species.h"
#include "photons.h"
#include "scalar_field.h"
#include "species_build_context.h"
#include "ultra_relativistic.h"
```

- [ ] **Step 2: Update the one include site**

In `source/input_module.cpp`, find:

```cpp
#include "../species/all_species_list.h"
```

Replace with:

```cpp
#include "../species/all_species.h"
```

- [ ] **Step 3: Delete the old file**

```bash
git rm species/all_species_list.h
```

- [ ] **Step 4: Build systems**

In `Makefile`, if `all_species_list.h` is referenced by name, replace with `all_species.h`. Otherwise no change.

In `setup.py`, same — replace `all_species_list.h` with `all_species.h` if listed.

In `CLASS.xcodeproj/project.pbxproj` — flag for the user to update manually (rename `all_species_list.h` → `all_species.h`).

- [ ] **Step 5: Build check**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, explanatory still runs.

- [ ] **Step 6: Commit**

```bash
git add species/all_species.h source/input_module.cpp Makefile setup.py
git commit -m "species: rename all_species_list.h to all_species.h, add build_context include"
```

---

## Task 5: Add `CreateAll` to `PhotonsSpecies` (always returns 1 entry)

**Files:**
- Modify: `species/photons.h`
- Modify: `species/photons.cpp`

- [ ] **Step 1: Add the static method to the header**

In `species/photons.h`, add `#include "species_build_context.h"` near the top. In the `public:` section of `PhotonsSpecies`, add:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Add the implementation**

In `species/photons.cpp`, add at the bottom of the file:

```cpp
std::vector<Named> PhotonsSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  result.push_back({"Photons", std::make_unique<PhotonsSpecies>(*ctx.pba)});
  return result;
}
```

If `<memory>` and `<vector>` aren't already included transitively, add `#include <memory>` and `#include <vector>` at the top.

- [ ] **Step 3: Build check**

```bash
make class -j 2>&1 | tail -5
```

Expected: clean build (new method is not yet called).

- [ ] **Step 4: Commit**

```bash
git add species/photons.h species/photons.cpp
git commit -m "species: add CreateAll factory to PhotonsSpecies"
```

---

## Task 6: Add `CreateAll` to `BaryonsSpecies` (always returns 1 entry)

**Files:**
- Modify: `species/baryons.h`
- Modify: `species/baryons.cpp`

- [ ] **Step 1: Header**

In `species/baryons.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

In `species/baryons.cpp`:

```cpp
std::vector<Named> BaryonsSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  result.push_back({"Baryons", std::make_unique<BaryonsSpecies>(*ctx.pba)});
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/baryons.h species/baryons.cpp
git commit -m "species: add CreateAll factory to BaryonsSpecies"
```

---

## Task 7: Add `CreateAll` to `CDMSpecies` (guarded by `has_cdm`)

**Files:**
- Modify: `species/cdm.h`
- Modify: `species/cdm.cpp`

- [ ] **Step 1: Header**

In `species/cdm.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

In `species/cdm.cpp`:

```cpp
std::vector<Named> CDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_cdm == _TRUE_) {
    result.push_back({"CDM", std::make_unique<CDMSpecies>(*ctx.pba)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/cdm.h species/cdm.cpp
git commit -m "species: add CreateAll factory to CDMSpecies"
```

---

## Task 8: Add `CreateAll` to `LambdaSpecies` (guarded by `has_lambda`)

**Files:**
- Modify: `species/lambda.h`
- Modify: `species/lambda.cpp` (if it exists; if not, add to header)

- [ ] **Step 1: Check if `.cpp` exists**

```bash
ls species/lambda.cpp 2>/dev/null && echo "exists" || echo "header-only"
```

If header-only, the implementation can go inline in the header (in an `inline` definition). For consistency with the other species, create `species/lambda.cpp` if not present.

- [ ] **Step 2: Header**

In `species/lambda.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 3: Implementation**

In `species/lambda.cpp` (create if absent, with appropriate `#include "lambda.h"`):

```cpp
std::vector<Named> LambdaSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_lambda == _TRUE_) {
    result.push_back({"Lambda", std::make_unique<LambdaSpecies>(*ctx.pba)});
  }
  return result;
}
```

If `lambda.cpp` did not exist, register it in `Makefile` (under `SPECIES_OPP` or equivalent) and `setup.py` (species source list). Flag for Xcode manual add.

- [ ] **Step 4: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/lambda.h species/lambda.cpp Makefile setup.py
git commit -m "species: add CreateAll factory to LambdaSpecies"
```

---

## Task 9: Add `CreateAll` to `UltraRelativisticSpecies` (guarded by `has_ur`)

**Files:**
- Modify: `species/ultra_relativistic.h`
- Modify: `species/ultra_relativistic.cpp` (create if absent, same as Task 8)

- [ ] **Step 1: Header**

In `species/ultra_relativistic.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

```cpp
std::vector<Named> UltraRelativisticSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_ur == _TRUE_) {
    result.push_back({"UR", std::make_unique<UltraRelativisticSpecies>(*ctx.pba)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/ultra_relativistic.h species/ultra_relativistic.cpp
git commit -m "species: add CreateAll factory to UltraRelativisticSpecies"
```

---

## Task 10: Add `CreateAll` to `FluidSpecies` (guarded by `has_fld`)

**Files:**
- Modify: `species/fluid.h`
- Modify: `species/fluid.cpp` (create if absent)

- [ ] **Step 1: Header**

In `species/fluid.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

```cpp
std::vector<Named> FluidSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_fld == _TRUE_) {
    result.push_back({"Fluid", std::make_unique<FluidSpecies>(*ctx.pba)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/fluid.h species/fluid.cpp
git commit -m "species: add CreateAll factory to FluidSpecies"
```

---

## Task 11: Add `CreateAll` to `DCDM_DR_Species` (guarded by `has_dcdm`, takes `pba` + `bgm`)

**Files:**
- Modify: `species/dcdm_dr_species.h`
- Modify: `species/dcdm_dr_species.cpp`

- [ ] **Step 1: Header**

In `species/dcdm_dr_species.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

In `species/dcdm_dr_species.cpp`:

```cpp
std::vector<Named> DCDM_DR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_dcdm == _TRUE_) {
    result.push_back({"DCDM_DR", std::make_unique<DCDM_DR_Species>(ctx.pba, ctx.bgm)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/dcdm_dr_species.h species/dcdm_dr_species.cpp
git commit -m "species: add CreateAll factory to DCDM_DR_Species"
```

---

## Task 12: Add `CreateAll` to `ScalarFieldSpecies` (guarded by `has_scf`)

**Files:**
- Modify: `species/scalar_field.h`
- Modify: `species/scalar_field.cpp` (create if absent)

- [ ] **Step 1: Header**

In `species/scalar_field.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

```cpp
std::vector<Named> ScalarFieldSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_scf == _TRUE_) {
    result.push_back({"ScalarField", std::make_unique<ScalarFieldSpecies>(*ctx.pba)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/scalar_field.h species/scalar_field.cpp
git commit -m "species: add CreateAll factory to ScalarFieldSpecies"
```

---

## Task 13: Add `CreateAll` to `IDM_DR_IDR_Species` (guarded by `has_idm_dr || has_idr`)

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`

- [ ] **Step 1: Header**

In `species/idm_dr_idr_species.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

```cpp
std::vector<Named> IDM_DR_IDR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_idm_dr == _TRUE_ || ctx.pba->has_idr == _TRUE_) {
    result.push_back({"IDM_DR_IDR", std::make_unique<IDM_DR_IDR_Species>(*ctx.pba)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp
git commit -m "species: add CreateAll factory to IDM_DR_IDR_Species"
```

---

## Task 14: Add `CreateAll` to `IDM_DRMD_IDR_DRMD_Species` (guarded by `has_idm_drmd || has_idr_drmd`)

**Files:**
- Modify: `species/idm_drmd_idr_drmd_species.h`
- Modify: `species/idm_drmd_idr_drmd_species.cpp`

- [ ] **Step 1: Header**

In `species/idm_drmd_idr_drmd_species.h`, add `#include "species_build_context.h"`. In the `public:` section:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Implementation**

```cpp
std::vector<Named> IDM_DRMD_IDR_DRMD_Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  if (ctx.pba->has_idm_drmd == _TRUE_ || ctx.pba->has_idr_drmd == _TRUE_) {
    result.push_back(
        {"IDM_DRMD_IDR_DRMD", std::make_unique<IDM_DRMD_IDR_DRMD_Species>(*ctx.pba)});
  }
  return result;
}
```

- [ ] **Step 3: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp
git commit -m "species: add CreateAll factory to IDM_DRMD_IDR_DRMD_Species"
```

---

## Task 15: Migrate `NCDMSpecies::CreateAll` signature

The existing `NCDMSpecies::CreateAll(FileContent*, const NcdmSettings&, const background*, const BackgroundModule*)` becomes `CreateAll(const SpeciesBuildContext&)`. The nested `NCDMSpecies::Named` type goes away — uses the global `Named` from `species_build_context.h`.

**Files:**
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp`
- Modify: `source/input_module.cpp` (two call sites: line ~232 in `ConstructSpecies`, line ~1089 in the temp block — both updated together so the build stays green)

- [ ] **Step 1: Update header**

In `species/ncdm_species.h`:
- Add `#include "species_build_context.h"`.
- Remove the nested `struct Named { std::string key; std::unique_ptr<NCDMSpecies> species; };` declaration.
- Change the `CreateAll` declaration to:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Update implementation**

In `species/ncdm_species.cpp`, change the `CreateAll` definition. The body stays the same; only the signature, the parameter references, and the `Named` construction change. Replace `pfc` → `ctx.pfc`, `settings` → `ctx.ncdm_settings`, `pba` → `ctx.pba`, `bgm` → `ctx.bgm`. The `result.push_back(Named{keys[n], std::move(sp)});` line continues to work — `Named` is now the shared type, and `unique_ptr<NCDMSpecies>` implicitly converts to `unique_ptr<BaseSpecies>`.

```cpp
std::vector<Named> NCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  const std::vector<std::string> dot_instances =
      ctx.pfc->instances_with("type", "ncdm_standard");

  int int1 = 0, int2 = 0, flag1 = 0, flag2 = 0;
  char errmsg[2048];
  class_call(parser_read_int(ctx.pfc, "N_ncdm_standard", &int1, &flag1, errmsg),
             errmsg, errmsg);
  class_call(parser_read_int(ctx.pfc, "N_ncdm", &int2, &flag2, errmsg),
             errmsg, errmsg);
  if (flag1 == _TRUE_ && flag2 == _TRUE_) {
    throw std::invalid_argument(
        "In input file, you can only enter one of N_ncdm_standard and N_ncdm, choose one");
  }
  int N_ncdm_standard = 0;
  if (flag1 == _TRUE_)
    N_ncdm_standard = int1;
  else if (flag2 == _TRUE_)
    N_ncdm_standard = int2;

  const bool has_dot    = !dot_instances.empty();
  const bool has_legacy = (flag1 == _TRUE_ || flag2 == _TRUE_);
  if (has_dot && has_legacy && has_unconsumed_dot_type_keys(*ctx.pfc, dot_instances)) {
    throw std::invalid_argument(
        "cannot mix legacy N_ncdm/N_ncdm_standard with dot-syntax "
        "'*.type = ncdm_standard'; use one or the other");
  }

  std::vector<std::string> keys;
  if (has_dot) {
    keys            = synthesise_standard_ncdm_flat_keys(*ctx.pfc, dot_instances);
    N_ncdm_standard = static_cast<int>(dot_instances.size());
  }
  else {
    keys.reserve(N_ncdm_standard);
    for (int n = 0; n < N_ncdm_standard; ++n) {
      keys.push_back("ncdm__" + std::to_string(n + 1));
    }
  }

  for (int n = 0; n < N_ncdm_standard; ++n) {
    auto sp = std::make_unique<NCDMSpecies>(ctx.pfc, n, ctx.ncdm_settings, ctx.pba, ctx.bgm);
    result.push_back({keys[n], std::move(sp)});
  }
  return result;
}
```

- [ ] **Step 3: Update both call sites in `input_module.cpp`**

Two existing calls to `NCDMSpecies::CreateAll(pfc, settings, pba, bgm)`:

1. In `InputModule::ConstructSpecies` (around line 232). Build a `SpeciesBuildContext`:

```cpp
    NcdmSettings ncdm_settings_for_species;
    ncdm_settings_for_species.h           = pba->h;
    ncdm_settings_for_species.T_cmb       = pba->T_cmb;
    ncdm_settings_for_species.tol_ncdm    = precision_.tol_ncdm;
    ncdm_settings_for_species.tol_ncdm_bg = precision_.tol_ncdm_bg;
    ncdm_settings_for_species.tol_M_ncdm  = precision_.tol_M_ncdm;

    SpeciesBuildContext ctx{
        /*pfc=*/&file_content_,
        /*pba=*/pba,
        /*ppr=*/&precision_,
        /*ncdm_settings=*/ncdm_settings_for_species,
        /*bgm=*/nullptr,
    };

    auto ncdm_list = NCDMSpecies::CreateAll(ctx);
    for (auto& e : ncdm_list) {
      all_species_.insert(e.key, std::move(e.species));
    }
```

(Note: `Named` no longer has typed `species` — it holds `unique_ptr<BaseSpecies>`. The code already moves into `all_species_.insert`, which takes `unique_ptr<BaseSpecies>`; this just works.)

2. In the temp block at lines 1087–1105: also build a `SpeciesBuildContext` and call `NCDMSpecies::CreateAll(ctx)`. The temp block survives Tasks 15–18 (signatures only); it's deleted later in Task 23.

- [ ] **Step 4: Build + smoke**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, explanatory still runs.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp source/input_module.cpp
git commit -m "ncdm: migrate NCDMSpecies::CreateAll to SpeciesBuildContext"
```

---

## Task 16: Migrate `DNCDMSpecies::CreateAll` signature

`DNCDMSpecies::CreateAll` is an internal helper — it's called by `DNCDM_DR_Species::CreateAll` to build the children of a composite. It is **not** in the factory list. Its return type stays a concrete `vector<DNCDMSpecies::Named>` because the caller needs the concrete type to build the composite. Only the input changes: `(FileContent*, NcdmSettings, background*, BackgroundModule*)` → `(const SpeciesBuildContext&)`.

**Files:**
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`
- Modify: `species/dncdm_dr_species.cpp` (the one caller)

- [ ] **Step 1: Update header**

In `species/dncdm_species.h`, add `#include "species_build_context.h"`. The nested `DNCDMSpecies::Named` stays — it holds `unique_ptr<DNCDMSpecies>` (concrete) and is consumed only by `DNCDM_DR_Species::CreateAll`. Update the `CreateAll` declaration:

```cpp
  struct Named {
    std::string                   key;
    std::unique_ptr<DNCDMSpecies> species;
  };
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Update implementation**

In `species/dncdm_species.cpp`, change the function signature and replace internal references (`pfc` → `ctx.pfc`, etc.). Body otherwise unchanged.

- [ ] **Step 3: Update the caller in `dncdm_dr_species.cpp`**

Find the line `auto dncdm_vec = DNCDMSpecies::CreateAll(pfc, settings, pba, bgm);`. Build a `SpeciesBuildContext` from the existing arguments:

```cpp
  SpeciesBuildContext ctx{pfc, pba, /*ppr=*/nullptr, settings, bgm};
  auto dncdm_vec = DNCDMSpecies::CreateAll(ctx);
```

(`ppr` is unused inside `DNCDMSpecies::CreateAll` so passing `nullptr` is fine. Verify with `grep "ppr\|precision" species/dncdm_species.cpp` before committing — if `ppr` IS used, add the precision pointer to the context construction.)

- [ ] **Step 4: Build + commit**

```bash
make class -j 2>&1 | tail -5
git add species/dncdm_species.h species/dncdm_species.cpp species/dncdm_dr_species.cpp
git commit -m "ncdm: migrate DNCDMSpecies::CreateAll to SpeciesBuildContext (internal helper)"
```

---

## Task 17: Migrate `DNCDM_DR_Species::CreateAll` signature

**Files:**
- Modify: `species/dncdm_dr_species.h`
- Modify: `species/dncdm_dr_species.cpp`
- Modify: `source/input_module.cpp` (two call sites)

- [ ] **Step 1: Update header**

In `species/dncdm_dr_species.h`, add `#include "species_build_context.h"`. Remove the nested `Named` struct (replaced by global `Named`). Update declaration:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Update implementation**

In `species/dncdm_dr_species.cpp`, change the signature. The current body has no `has_*` guard because `DNCDMSpecies::CreateAll` already returns an empty vector when no DNCDM is requested. So no guard is needed at this level either:

```cpp
std::vector<Named> DNCDM_DR_Species::CreateAll(const SpeciesBuildContext& ctx) {
  auto dncdm_vec = DNCDMSpecies::CreateAll(ctx);
  std::vector<Named> result;
  result.reserve(dncdm_vec.size());
  for (auto& e : dncdm_vec) {
    result.push_back(
        {e.key,
         std::make_unique<DNCDM_DR_Species>(std::move(e.species), ctx.pba, ctx.bgm)});
  }
  return result;
}
```

`e.species` is `unique_ptr<DNCDMSpecies>` (from `DNCDMSpecies::Named`, the concrete-typed internal helper that survived Task 16). The `make_unique<DNCDM_DR_Species>` call wraps it into the composite, which is then pushed into a global `Named` (from `species_build_context.h`) holding `unique_ptr<BaseSpecies>`.

- [ ] **Step 3: Update both call sites in `input_module.cpp`**

The same `SpeciesBuildContext ctx{...}` already constructed in Task 15 is reused. Replace `DNCDM_DR_Species::CreateAll(pfc, settings, pba, bgm)` (or its equivalent) with `DNCDM_DR_Species::CreateAll(ctx)` at both sites.

- [ ] **Step 4: Build + smoke + commit**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
git add species/dncdm_dr_species.h species/dncdm_dr_species.cpp source/input_module.cpp
git commit -m "ncdm: migrate DNCDM_DR_Species::CreateAll to SpeciesBuildContext"
```

---

## Task 18: Migrate `NCDMInteractingSpecies::CreateAll` signature

**Files:**
- Modify: `species/ncdm_interacting_species.h`
- Modify: `species/ncdm_interacting_species.cpp`
- Modify: `source/input_module.cpp` (two call sites)

- [ ] **Step 1: Update header**

Same pattern as Task 15: add `#include "species_build_context.h"`, drop the nested `Named`, change declaration:

```cpp
  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);
```

- [ ] **Step 2: Update implementation**

In `species/ncdm_interacting_species.cpp`, change the signature; replace `pfc`/`settings`/`pba`/`bgm` references with `ctx.*`. The existing body already handles the "no NCDM_interacting requested" case (returns empty when both `dot_instances` is empty and `flag` is false), so no extra guard is needed:

```cpp
std::vector<Named> NCDMInteractingSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;
  const auto dot_instances = ctx.pfc->instances_with("type", "ncdm_self_interacting");

  int N_ncdm_interacting = 0;
  bool flag              = ctx.pfc->read_int("N_ncdm_interacting", N_ncdm_interacting);
  const bool has_legacy  = flag;
  const bool has_dot     = !dot_instances.empty();
  if (has_legacy && has_dot && has_unconsumed_dot_type_keys(*ctx.pfc, dot_instances)) {
    throw std::invalid_argument(
        "cannot mix N_ncdm_interacting with dot-syntax "
        "'*.type = ncdm_self_interacting'; use one or the other");
  }

  std::vector<std::string> keys;
  int N = 0;
  if (has_dot) {
    keys = synthesise_self_interacting_ncdm_flat_keys(*ctx.pfc, dot_instances);
    N    = static_cast<int>(dot_instances.size());
  }
  else if (has_legacy) {
    int N_std = 0, N_dr = 0;
    ctx.pfc->read_int("N_ncdm_standard", N_std);
    ctx.pfc->read_int("N_ncdm_decay_dr", N_dr);
    keys.reserve(N_ncdm_interacting);
    for (int n = 0; n < N_ncdm_interacting; ++n) {
      keys.push_back("ncdm__" + std::to_string(N_std + N_dr + n + 1));
    }
    N = N_ncdm_interacting;
  }

  for (int n = 0; n < N; ++n) {
    auto sp = std::make_unique<NCDMInteractingSpecies>(
        ctx.pfc, n, ctx.ncdm_settings, ctx.pba, ctx.bgm);
    result.push_back({keys[n], std::move(sp)});
  }
  return result;
}
```

- [ ] **Step 3: Update both call sites in `input_module.cpp`**

Replace `NCDMInteractingSpecies::CreateAll(pfc, settings, pba, bgm)` with `NCDMInteractingSpecies::CreateAll(ctx)` at both sites (one in `ConstructSpecies`, one in the temp block — temp block gets deleted in Task 23).

- [ ] **Step 4: Build + smoke + commit**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
git add species/ncdm_interacting_species.h species/ncdm_interacting_species.cpp source/input_module.cpp
git commit -m "ncdm: migrate NCDMInteractingSpecies::CreateAll to SpeciesBuildContext"
```

---

## Task 19: Populate `kAllSpeciesFactories` in `species/all_species.h`

**Files:**
- Modify: `species/all_species.h`

- [ ] **Step 1: Add the factory list**

Append to `species/all_species.h` (after the existing includes):

```cpp
#include <array>
#include <string_view>

struct SpeciesFactoryEntry {
  std::string_view  name;
  std::vector<Named> (*create_all)(const SpeciesBuildContext&);
};

inline constexpr std::array kAllSpeciesFactories = {
  SpeciesFactoryEntry{"Photons",            &PhotonsSpecies::CreateAll},
  SpeciesFactoryEntry{"Baryons",            &BaryonsSpecies::CreateAll},
  SpeciesFactoryEntry{"CDM",                &CDMSpecies::CreateAll},
  SpeciesFactoryEntry{"UR",                 &UltraRelativisticSpecies::CreateAll},
  SpeciesFactoryEntry{"DCDM_DR",            &DCDM_DR_Species::CreateAll},
  SpeciesFactoryEntry{"NCDM",               &NCDMSpecies::CreateAll},
  SpeciesFactoryEntry{"DNCDM_DR",           &DNCDM_DR_Species::CreateAll},
  SpeciesFactoryEntry{"NCDMInt",            &NCDMInteractingSpecies::CreateAll},
  SpeciesFactoryEntry{"IDM_DR_IDR",         &IDM_DR_IDR_Species::CreateAll},
  SpeciesFactoryEntry{"IDM_DRMD_IDR_DRMD",  &IDM_DRMD_IDR_DRMD_Species::CreateAll},
  SpeciesFactoryEntry{"Lambda",             &LambdaSpecies::CreateAll},
  SpeciesFactoryEntry{"Fluid",              &FluidSpecies::CreateAll},
  SpeciesFactoryEntry{"ScalarField",        &ScalarFieldSpecies::CreateAll},
};

/** Maps the closure-species enum to the matching factory entry name. */
inline constexpr std::string_view ClosureSpeciesName(background::ClosureSpecies cs) {
  switch (cs) {
    case background::ClosureSpecies::Lambda:      return "Lambda";
    case background::ClosureSpecies::Fluid:       return "Fluid";
    case background::ClosureSpecies::ScalarField: return "ScalarField";
    case background::ClosureSpecies::None:        return "";
  }
  return "";
}
```

(If `background::ClosureSpecies` is defined as a free enum in the global namespace rather than nested in `background`, drop the `background::` prefix.)

- [ ] **Step 2: Build check**

```bash
make class -j 2>&1 | tail -5
```

Expected: clean build (the factory list is unused so far).

- [ ] **Step 3: Commit**

```bash
git add species/all_species.h
git commit -m "species: declare kAllSpeciesFactories factory list"
```

---

## Task 20: Switch `ConstructSpecies` to a uniform factory-loop

This task replaces the body of `InputModule::ConstructSpecies` (currently `source/input_module.cpp` lines ~205–265) with a generic loop. **Closure handling is unchanged in this task** — `pba->Omega0_lambda` (or fld / scf) is still computed in `input_read_parameters` and is set on `pba` before `ConstructSpecies` runs, so the closure species' `CreateAll` reads the correct value. Closure detection moves in Task 21.

**Files:**
- Modify: `source/input_module.cpp` (`ConstructSpecies` body)

- [ ] **Step 1: Rewrite `ConstructSpecies`**

Replace the entire body of `ConstructSpecies` with:

```cpp
void InputModule::ConstructSpecies() {
  const background* pba = &background_;

  NcdmSettings ncdm_settings;
  ncdm_settings.h           = pba->h;
  ncdm_settings.T_cmb       = pba->T_cmb;
  ncdm_settings.tol_ncdm    = precision_.tol_ncdm;
  ncdm_settings.tol_ncdm_bg = precision_.tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = precision_.tol_M_ncdm;

  const SpeciesBuildContext ctx{
      /*pfc=*/&file_content_,
      /*pba=*/pba,
      /*ppr=*/&precision_,
      /*ncdm_settings=*/ncdm_settings,
      /*bgm=*/nullptr,
  };

  for (const auto& entry : kAllSpeciesFactories) {
    auto produced = entry.create_all(ctx);
    for (auto& e : produced) {
      all_species_.insert(e.key, std::move(e.species));
    }
  }

  all_species_.freeze();
}
```

- [ ] **Step 2: Build + smoke**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, explanatory still runs.

- [ ] **Step 3: Reference gate**

```bash
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

Expected: all reference scenarios still pass. If any fail, investigate before proceeding.

- [ ] **Step 4: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: switch ConstructSpecies to uniform factory-loop"
```

---

## Task 21: Move closure-value computation from `input_read_parameters` to `ConstructSpecies`

This is the task that actually changes when closure happens. After this task, `input_read_parameters` only *records* `pba->closure_species`; `ConstructSpecies` builds non-closure species, sums their Omega0, computes the closure value, and runs the closure species' factory last.

**Files:**
- Modify: `source/input_module.cpp` (closure block at lines ~1190–1245, and `ConstructSpecies`)

- [ ] **Step 1: Replace closure-value-computation in `input_read_parameters` with closure-species detection**

Find lines ~1212–1243 in `source/input_module.cpp` (the `/* Step 1 */` and `/* Step 2 */` blocks). Keep `Step 1` as-is (it sets `pba->Omega0_X = paramN` for any directly-specified value and accumulates into `Omega_tot` — that accumulation will be deleted in Task 24, leave it for now). Replace `Step 2` (the `if (flag1 == _FALSE_) { pba->Omega0_lambda = 1. - …; }` etc. cascade) with closure-species recording only:

```cpp
  /* Step 2: record which species (if any) is the closure species. */
  if (flag1 == _FALSE_) {
    pba->closure_species = background::ClosureSpecies::Lambda;
  }
  else if (flag2 == _FALSE_) {
    pba->closure_species = background::ClosureSpecies::Fluid;
  }
  else if ((flag3 == _TRUE_) && (param3 < 0.)) {
    pba->closure_species = background::ClosureSpecies::ScalarField;
  }
  else {
    pba->closure_species = background::ClosureSpecies::None;
  }
```

The verbose-print at line ~1230 (`" -> matched budget equations by adjusting Omega_Lambda = %e\n"`) — move it to `ConstructSpecies` (after the closure value is computed, see Step 2 below).

- [ ] **Step 2: Update `ConstructSpecies` to handle closure**

Replace the loop in `ConstructSpecies` from Task 20 with a closure-aware version:

```cpp
void InputModule::ConstructSpecies() {
  background* pba = &background_;   // need non-const to write closure Omega

  NcdmSettings ncdm_settings;
  ncdm_settings.h           = pba->h;
  ncdm_settings.T_cmb       = pba->T_cmb;
  ncdm_settings.tol_ncdm    = precision_.tol_ncdm;
  ncdm_settings.tol_ncdm_bg = precision_.tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = precision_.tol_M_ncdm;

  const SpeciesBuildContext ctx{
      /*pfc=*/&file_content_,
      /*pba=*/pba,
      /*ppr=*/&precision_,
      /*ncdm_settings=*/ncdm_settings,
      /*bgm=*/nullptr,
  };

  const std::string_view closure_name = ClosureSpeciesName(pba->closure_species);

  // Pass 1: build every non-closure species, accumulate their Omega0.
  double omega0_sum = 0.0;
  for (const auto& entry : kAllSpeciesFactories) {
    if (entry.name == closure_name) {
      continue;
    }
    auto produced = entry.create_all(ctx);
    for (auto& e : produced) {
      omega0_sum += e.species->GetOmega0();
      all_species_.insert(e.key, std::move(e.species));
    }
  }

  // Compute closure value and write to pba so the closure species' factory
  // reads it like every other species reads its own Omega0_X.
  if (pba->closure_species != background::ClosureSpecies::None) {
    const double closure_value = 1.0 - pba->Omega0_k - omega0_sum;
    switch (pba->closure_species) {
      case background::ClosureSpecies::Lambda:
        pba->Omega0_lambda = closure_value;
        if (input_verbose_ > 0) {
          printf(" -> matched budget equations by adjusting Omega_Lambda = %e\n",
                 pba->Omega0_lambda);
        }
        break;
      case background::ClosureSpecies::Fluid:
        pba->Omega0_fld = closure_value;
        if (input_verbose_ > 0) {
          printf(" -> matched budget equations by adjusting Omega_fld = %e\n",
                 pba->Omega0_fld);
        }
        break;
      case background::ClosureSpecies::ScalarField:
        pba->Omega0_scf = closure_value;
        if (input_verbose_ > 0) {
          printf(" -> matched budget equations by adjusting Omega_scf = %e\n",
                 pba->Omega0_scf);
        }
        break;
      case background::ClosureSpecies::None:
        break;
    }

    // Pass 2: build the closure species now that pba->Omega0_X is set.
    for (const auto& entry : kAllSpeciesFactories) {
      if (entry.name != closure_name) {
        continue;
      }
      auto produced = entry.create_all(ctx);
      for (auto& e : produced) {
        all_species_.insert(e.key, std::move(e.species));
      }
      break;
    }
  }

  all_species_.freeze();
}
```

`input_verbose_` may need to be looked up — confirm by `grep -n "input_verbose" source/input_module.h source/input_module.cpp`. If it's stored as a member, use it; otherwise capture it from `precision_` or pass through the context.

- [ ] **Step 3: Build + smoke**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, explanatory still runs.

- [ ] **Step 4: Reference gate**

```bash
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

Expected: all reference scenarios pass.

- [ ] **Step 5: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: defer Omega closure to ConstructSpecies"
```

---

## Task 22: Investigate the shooting interaction

Before deleting the temp-NCDM block, confirm whether shooting iterations actually depend on the per-iteration values of `pba->Omega0_ncdm_tot`, `pba->N_ncdm`, `pba->N_decay_dr`. Two possible outcomes drive Task 23.

**Files:** none (investigation).

- [ ] **Step 1: Trace shooting flow**

Read `source/input_module.cpp` around lines 478–500 (the shooting trigger) and the shooting workspace's `input_try_unknown_parameters` / `input_fzerofun_1d`. Identify whether shooting iterations:
- (a) re-enter `input_init` recursively, in which case the temp-NCDM block runs each iteration and `pba->Omega0_ncdm_tot` is set per-iteration; or
- (b) only check a single residual that doesn't depend on `Omega0_ncdm_tot` / `N_ncdm` / `N_decay_dr`.

Run a single shooting input under instrumentation: add a temporary `printf("shooting iter: Omega0_ncdm_tot = %e\n", pba->Omega0_ncdm_tot);` at the top of `input_default_params` (where `Omega0_ncdm_tot` is read at line ~3021) and run `./class <some_shooting_input>.ini`. Note whether the value changes between iterations.

- [ ] **Step 2: Decide**

If the value does **not** change meaningfully between iterations (i.e. shooting doesn't depend on it), Task 23 simply deletes the temp block.

If the value **does** change between iterations (shooting depends on per-iteration NCDM Omega), Task 23 must additionally provide static helpers `NCDMSpecies::TotalOmega0FromInput(const SpeciesBuildContext&)`, `DNCDMSpecies::TotalOmega0FromInput(const SpeciesBuildContext&)`, `NCDMInteractingSpecies::TotalOmega0FromInput(const SpeciesBuildContext&)` (and corresponding count helpers) that compute the totals from `pfc` without instantiating species. The shooting iteration calls these instead of the deleted block; the post-shooting final state still uses the unified `ConstructSpecies` loop.

Document the choice in the commit message of Task 23.

- [ ] **Step 3: Remove instrumentation**

Delete the temporary `printf` from Step 1.

- [ ] **Step 4: No commit**

Investigation only — no code lands.

---

## Task 23: Delete the temp-NCDM block; populate `N_ncdm` / `N_decay_dr` / `Omega0_ncdm_tot` from built species

This is the destructive task. After it lands, the temp block at `source/input_module.cpp` lines ~1087–1105 is gone. The three counter fields are populated in `ConstructSpecies` from the actually-built species.

**Files:**
- Modify: `source/input_module.cpp` (delete temp block; update `ConstructSpecies`; potentially relocate the `N_ncdm > 0` test at line ~2843)
- Possibly modify: `species/ncdm_species.{h,cpp}`, `species/dncdm_species.{h,cpp}`, `species/ncdm_interacting_species.{h,cpp}` to add static `TotalOmega0FromInput` / `CountFromInput` helpers (only if Task 22 concluded shooting requires them)

- [ ] **Step 1: Delete the temp block**

In `source/input_module.cpp`, delete the entire braced block lines ~1087–1105 starting with `// Build NCDM species temporarily to get N_ncdm, Omega0_ncdm_tot, and DNCDM count`. Also delete line 1108: `Omega_tot += pba->Omega0_ncdm_tot;`.

Keep lines 1080–1086 (the `NcdmSettings ncdm_settings; …` setup) for now; this block will be evaluated for deletion in Task 24 (it's no longer needed at this point in `input_read_parameters` since species are no longer built here, but leaving it for one extra task keeps the diff focused).

- [ ] **Step 2: Populate the three counters in `ConstructSpecies` post-loop**

Inside `ConstructSpecies`, after Pass 1 (and before or after Pass 2 — choose after, when all species are built), iterate `all_species_` and tally:

```cpp
  // Tally NCDM-family counters from the actually-built species. These remain
  // load-bearing for downstream modules (thermo, nonlinear, etc.); removing
  // them from those modules is out of scope for this refactor.
  pba->N_ncdm          = 0;
  pba->N_decay_dr      = (pba->Omega0_dcdmdr > 0 ? 1 : 0);
  pba->Omega0_ncdm_tot = 0.0;
  for (const auto& [name, sp] : all_species_) {
    // Match by SpeciesCollection key prefix — these are the NCDM-family sectors
    // that contribute to N_ncdm. NCDMSpecies keys are "ncdm__N", DNCDM_DR_Species
    // keys are "dncdm_dr__N", NCDMInteractingSpecies keys are "ncdm_interacting__N"
    // (verify exact key formats by reading each CreateAll body).
    if (name.rfind("ncdm__", 0) == 0
        || name.rfind("dncdm_dr__", 0) == 0
        || name.rfind("ncdm_interacting__", 0) == 0) {
      pba->N_ncdm += 1;
      pba->Omega0_ncdm_tot += sp->GetOmega0();
    }
    if (name.rfind("dncdm_dr__", 0) == 0) {
      pba->N_decay_dr += 1;
    }
  }
```

This counter-tallying is type-introspection-by-name-prefix, which is a code smell (key-prefix matching is fragile). It is acceptable here as a *transient* arrangement: the eventual goal (out of scope) is to remove these counters from module code entirely. Note this in a code comment.

If you find this prefix-matching too brittle, an alternative is to add a virtual `bool IsNcdmFamily() const { return false; }` to `BaseSpecies`, override in NCDM/DNCDM_DR/NCDMInteracting to return true, and a virtual `bool IsDecayDr() const` similarly. Pick one approach and document the choice in the commit. The spec considers either acceptable; key-prefix is a lighter touch.

- [ ] **Step 3: Resolve the `N_ncdm > 0` test at line ~2843 in `input_read_parameters`**

The test at lines 2843–2855 (precision-parameter validation depending on `N_ncdm > 0`) cannot stay in `input_read_parameters` because `N_ncdm` is no longer set there. Move this test to the end of `ConstructSpecies` (after the loop, where `N_ncdm` is populated):

```cpp
  // Precision-parameter validation that depends on N_ncdm — moved here from
  // input_read_parameters when N_ncdm became a post-construction value.
  if (pba->N_ncdm > 0) {
    if (precision_.ncdm_fluid_trigger_tau_over_tau_k
        == precision_.radiation_streaming_trigger_tau_over_tau_k) {
      throw std::invalid_argument(
          "please choose different values for precision parameters "
          "ncdm_fluid_trigger_tau_over_tau_k and radiation_streaming_trigger_tau_over_tau_k, "
          "in order to avoid switching two approximation schemes at the same time");
    }
    if (precision_.ncdm_fluid_trigger_tau_over_tau_k
        == precision_.ur_fluid_trigger_tau_over_tau_k) {
      throw std::invalid_argument(
          "please choose different values for precision parameters "
          "ncdm_fluid_trigger_tau_over_tau_k and ur_fluid_trigger_tau_over_tau_k, in order to "
          "avoid switching two approximation schemes at the same time");
    }
  }
```

Read the original block at lines 2843–2855 for the exact error messages and any additional tests in that span; copy verbatim.

- [ ] **Step 4: Resolve `input_default_params` readers (lines ~3021, ~3030, ~3129)**

`input_default_params` is called at line 627 of `input_read_parameters`. It reads `pba->Omega0_ncdm_tot` (3021) and `pba->N_decay_dr` (3030, 3129). Per Task 22:
- If shooting does not depend on these per-iteration → these reads simply use whatever value is currently on `pba` (zero on the first call; possibly stale on later calls). Verify by re-reading those code paths what they do when the value is zero. If they correctly do nothing when zero, leave them alone.
- If shooting *does* depend on these per-iteration → introduce the static helpers from Task 22 Step 2 and replace the reads accordingly. Or: extract that section of `input_default_params` and call it from `ConstructSpecies` post-loop instead.

Document the chosen path in the commit message.

- [ ] **Step 5: Build + smoke**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, explanatory still runs.

- [ ] **Step 6: Reference gate (full)**

```bash
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -50
cd ..
```

Expected: all reference scenarios pass, including any NCDM, DNCDM_DR, NCDMInteracting, IDM_DRMD, and shooting-based scenarios. If any fail, investigate before proceeding.

- [ ] **Step 7: Commit**

```bash
git add source/input_module.cpp species/  # any species/ changes from helper additions
git commit -m "input: delete temp-NCDM block; populate N_ncdm/Omega0_ncdm_tot in ConstructSpecies"
```

---

## Task 24: Cleanup — drop `Omega_tot` accumulator and dead closure code

`Omega_tot` is no longer used now that closure happens in `ConstructSpecies`. Every `Omega_tot += …` line in `input_read_parameters` is dead.

**Files:**
- Modify: `source/input_module.cpp`

- [ ] **Step 1: Find and delete every `Omega_tot` reference in `input_read_parameters`**

```bash
grep -n "Omega_tot" source/input_module.cpp
```

Expected: lines like 614 (the declaration), 740, 753, 819, 865, 881, 913, 914, 918, 923, 1108 (already deleted in Task 23), 1187, 1188, 1215, 1219, 1223, 1228, 1234, 1240, 1251.

Delete the declaration at line 614 (`double Omega_tot;`) and every `Omega_tot += …` and `Omega_tot -= …` line in `input_read_parameters`. Also delete the closure-detection block's `Step 1` accumulators (`Omega_tot += pba->Omega0_lambda;` and similar at ~1215/1219/1223) — they're no longer needed since `ConstructSpecies` sums Omega0 from species.

The closure-detection block (Task 21) still needs the `flag1`/`flag2`/`flag3`/`param1`/`param2`/`param3` to identify the closure species; keep those `parser_read_double` calls and the validation tests at lines 1196–1200 and 1256–1260.

The commented-out fprintf at lines 1245–1252 can be deleted entirely (dead code).

- [ ] **Step 2: Delete the now-unused NCDM settings block at lines ~1080–1086**

This block was kept in Task 23 for diff focus. Now delete it:

```cpp
  /** - non-cold relics (ncdm) */
  NcdmSettings ncdm_settings;
  ncdm_settings.h           = pba->h;
  ncdm_settings.T_cmb       = pba->T_cmb;
  ncdm_settings.tol_ncdm    = ppr->tol_ncdm;
  ncdm_settings.tol_ncdm_bg = ppr->tol_ncdm_bg;
  ncdm_settings.tol_M_ncdm  = ppr->tol_M_ncdm;
```

(`ncdm_settings` is constructed inside `ConstructSpecies` instead.)

- [ ] **Step 3: Drop redundant verbose prints**

Per the design: any verbose-print in `input_read_parameters` whose only purpose was to report a value we no longer compute → drop the print rather than keep dead state. Likely candidates: any `printf` in the deleted block neighbourhood that printed `N_ncdm` or `Omega0_ncdm_tot` summaries. If found, delete the print.

- [ ] **Step 4: Build + smoke + reference gate**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

Expected: all reference scenarios pass.

- [ ] **Step 5: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: drop Omega_tot accumulator and dead closure code"
```

---

## Task 25: Build-system finalization

The Makefile and setup.py have already been touched incrementally as files were added. This task is a final consolidation pass plus the explicit Xcode-flag in the PR description.

**Files:**
- Possibly modify: `Makefile`, `setup.py`

- [ ] **Step 1: Verify all new files are in Makefile**

```bash
grep -c "species_build_context\.h\|all_species\.h" Makefile
```

If 0, the Makefile uses dependency generation and explicit headers aren't listed — that's fine. Otherwise, confirm both new headers are present.

- [ ] **Step 2: Verify all new species `.cpp` files are listed**

If Tasks 8 / 9 / 10 / 12 created new `.cpp` files (`lambda.cpp`, `ultra_relativistic.cpp`, `fluid.cpp`, `scalar_field.cpp`):

```bash
grep -E "lambda\.o|ultra_relativistic\.o|fluid\.o|scalar_field\.o" Makefile
grep -E "lambda\.cpp|ultra_relativistic\.cpp|fluid\.cpp|scalar_field\.cpp" setup.py
```

Expected: every new `.cpp` is registered in both files. Add any missing entries.

- [ ] **Step 3: Clean rebuild**

```bash
make clean && make class -j 2>&1 | tail -10
```

Expected: clean build from scratch.

- [ ] **Step 4: Note Xcode-update requirement in commit message and PR description**

The user maintains `CLASS.xcodeproj/project.pbxproj` manually. Final commit message in Task 26 should list:
- New header `species/species_build_context.h`
- Renamed file `species/all_species_list.h` → `species/all_species.h`
- Any new `.cpp` files added in Tasks 8 / 9 / 10 / 12 (whichever existed prior to this plan)

so the user can update Xcode in one pass.

- [ ] **Step 5: Commit (only if anything changed)**

```bash
git add Makefile setup.py
git commit -m "build: register new species infrastructure files"
```

If nothing changed in this task, skip the commit.

---

## Task 26: Final regression run and PR

**Files:** none (verification + git ops).

- [ ] **Step 1: Full reference gate**

```bash
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -50
cd ..
```

Expected: all reference scenarios pass.

- [ ] **Step 2: Run the targeted reference inputs from the spec**

Run each of these and confirm they complete without error:

```bash
./class explanatory.ini 2>&1 | tail -3
./class test/dotsyntax_ncdm.ini 2>&1 | tail -3
./class test/dotsyntax_ncdm_mixed.ini 2>&1 | tail -3
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini 2>&1 | tail -3
```

Each should end with `All parameters and all quantities computed successfully`.

- [ ] **Step 3: Closure correctness spot-check**

Pick one fluid-as-closure and one scalar-field-as-closure scenario from the test inputs (or construct one). Run with `input_verbose = 1` and confirm the printed `Omega_fld` / `Omega_scf` matches what master prints for the same input.

- [ ] **Step 4: Push and open PR**

```bash
git push -u origin unified-species-createall
gh pr create --title "Unified species CreateAll loop and deferred Omega closure" --body "$(cat <<'EOF'
## Summary

- Replaces the temp-NCDM construction inside `input_read_parameters` (lines 1087–1105) with a single uniform `CreateAll` loop in `ConstructSpecies`, driven by a constexpr factory list in `species/all_species.h` (renamed from `all_species_list.h`).
- Defers Omega closure to after species are built. `input_read_parameters` now only records `pba->closure_species`; `ConstructSpecies` sums `Omega0` over non-closure species, computes the closure value, writes it to `pba`, and runs the closure species' `CreateAll` last.
- Adds a static `CreateAll(const SpeciesBuildContext&)` to every species class. Adding a new species: write the class, include its header in `species/all_species.h`, add one row to `kAllSpeciesFactories`. No source/-side changes needed.
- Keeps `pba->N_ncdm`, `pba->N_decay_dr`, `pba->Omega0_ncdm_tot` populated (load-bearing for thermo / nonlinear / etc.) — populated in `ConstructSpecies` post-loop. Removing them entirely from module code is out of scope.

## Xcode

The user must update `CLASS.xcodeproj/project.pbxproj` to:
- Add `species/species_build_context.h`.
- Rename `species/all_species_list.h` → `species/all_species.h`.
- Add any new `.cpp` files created during this plan (lambda / ultra_relativistic / fluid / scalar_field — whichever didn't already have a `.cpp`).

## Test plan

- [x] Clean `make class -j` build
- [x] `./class explanatory.ini` succeeds
- [x] Full `python/test_class.py` reference gate passes (`TEST_LEVEL=2 COMPARE_OUTPUT_REF=1`)
- [x] Targeted scenarios pass: NCDM, DNCDM_DR, NCDMInteracting, IDM_DR_IDR, IDM_DRMD_IDR_DRMD, fluid-as-closure, scalar-field-as-closure, shooting

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 5: Done**

Done. The PR awaits review. Mention in the PR conversation that the Xcode project file needs the manual updates listed in the description.
