# Composite Child Iteration + Per-Child Matter Tally — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the perturbation matter tally correct for composites that straddle the matter/non-matter partition (Type3 = CDM + scalar field) by tallying per child, while giving `CompositeSpecies` generic child-iteration infrastructure the whole composite pipeline can reuse.

**Architecture:** Three layered changes. (A) `CompositeSpecies` gains an owning `child_layouts` vector aligned with `children_` and generic forwarding defaults, so `Type3Species` sheds its boilerplate. (B) A non-virtual `BaseSpecies::TallyStressEnergy` fuses stress-energy accumulation and the matter tally into one call that invokes the virtual `StressEnergy` once; composites delegate it per child via a virtual `DelegateTally` seam — one indirect call for plain species, no double-eval, `StressEnergy` stays pure. (C) The `ρ−3P`/`δρ−3δP` matter proxy is retired; the per-child gate already excludes radiation, so each clustering species contributes its actual `ρ`/`δρ`/`(ρ+P)θ`.

**Tech Stack:** C++17, CMake (single-config Release build in `build/`), CTest C++ unit tests (`species/*_test.cpp` linking `classpp`), Python integration tests (`python/test_class.py`, classyref reference workflow).

**Design doc:** `docs/superpowers/specs/2026-06-27-composite-child-iteration-and-matter-tally-design.md`

**Sequencing (part of the Type3 effort).** This plan resolves finding **I-2** of the Type-3 momentum-transfer work (`docs/superpowers/plans/2026-06-26-type3-momentum-transfer.md`; matter+DE composite polluting the cold `delta_m` → P(k) ~10× low). In the combined SDD ledger (`.superpowers/sdd/progress.md`) its tasks run **after** Type3 Task 4 and **before** Type3 Task 5 (the coupling): this plan's Task 1 → ledger **Task 5**, Task 2 → **Task 6**, Task 3 → **Task 7** (resolves I-2); the coupling and Fig-2 validation are ledger Tasks 8-9; this plan's Task 4 (verification) is **moved to the very end** as ledger **Task 10**.

## Global Constraints

- **Verification tolerance:** ≤0.1% on TT and mP(k); never require bit-identical output; handle Cl^TE zero-crossings (never blind max-rel-diff). Pure-FP-reordering drift at ULP level is accepted.
- **Git staging:** never `git add -A` / `git add .` in this repo (in-source CMake/Xcode artifacts get swept in). Stage explicit paths only; use `git add -f` for gitignored test `.ini` files.
- **`cclassy.pxd` is auto-generated** from C++ headers at build time — never hand-edit it. New public `BaseSpecies` methods are picked up by the wrapper generator on rebuild.
- **C++-only repo:** plain C++ headers, no `extern "C"`, no `#ifdef __cplusplus`.
- **Type3 is synchronous-gauge only** (already guarded in `Type3Species::CreateAll`); do not add Newtonian support.
- **Build (BUILD_DIR = `build/cmake`):** re-configure after editing `CMakeLists.txt` with `cmake -S . -B build/cmake`; build a target with `cmake --build build/cmake --target <t> --parallel`. Run unit tests with `ctest --test-dir build/cmake -R <t> --output-on-failure` (or `./build/cmake/<t>`, working dir = project root). The CLI shim `make class` builds `class` to the **repo root** (`./class <ini>`). Wrapper rebuild for Python tests: `make classy-pip-dev`.

---

## File Structure

| File | Responsibility | Change |
|------|----------------|--------|
| `species/base_species.h` | Species interface | Add cached `clusters_as_matter_`/`is_cold_`/`delegates_tally_`, `FinalizeMatterClassification()`, public cached accessors (Task 1); inline `TallyStressEnergy` + virtual `DelegateTally` (Task 3) |
| `species/composite_species.h` | Composite base interface | `FinalizeMatterClassification` override (Task 1); `PerturbLayout` w/ `child_layouts`, generic `CreatePerturbLayout` + 7 forwarder decls (Task 2); `DelegateTally` decl + ctor sets `delegates_tally_` (Task 3) |
| `species/composite_species.cpp` | Composite base impl | `FinalizeMatterClassification` (Task 1); generic forwarders (Task 2); `DelegateTally` (Task 3) |
| `species/type3_species.h` | Type3 composite interface | Drop `PerturbLayout`/`CreatePerturbLayout` + 8 forwarder decls; add `ChildIndex` enum + cast helpers (Task 2) |
| `species/type3_species.cpp` | Type3 composite impl | Delete generic-now forwarders; rewrite `PrintVariables` to use cast helpers (Task 2) |
| `species/species_collection.cpp` | Collection lifecycle | Call `FinalizeMatterClassification()` in `freeze()` (Task 1) |
| `source/perturbations.h` | Module workspace/views | Drop the two `ActiveSpecies` bools (Task 3) |
| `source/perturbations_module.cpp` | RHS hot path | Drop bools at 3 `active_species` construction sites; rewrite the scalar tally loop (Task 3) |
| `species/composite_classification_test.cpp` | NEW unit test | cached-bool invariant (Task 1) |
| `species/composite_layout_test.cpp` | NEW unit test | `child_layouts` alignment (Task 2) |
| `CMakeLists.txt` | Build | Register the two new test executables (Tasks 1, 2) |

---

## Task 1: Cached matter classification + finalize pass

Adds two cached booleans per species (`clusters_as_matter_`, `is_cold_`), stamped once after construction by `SpeciesCollection::freeze()` (composites recurse into children). Pure addition — no behavior change. The hot-path tally (Task 3) reads these instead of calling the virtuals per RHS step.

**Files:**
- Modify: `species/base_species.h` (protected members + public accessors + virtual `FinalizeMatterClassification`)
- Modify: `species/composite_species.h:131-133` (declare `FinalizeMatterClassification` override)
- Modify: `species/composite_species.cpp` (implement recursion)
- Modify: `species/species_collection.cpp:33-34` (call in `freeze()`)
- Create: `species/composite_classification_test.cpp`
- Modify: `CMakeLists.txt:214-219`

**Interfaces:**
- Produces: `bool BaseSpecies::ClustersAsMatterCached() const`, `bool BaseSpecies::IsColdCached() const`, `virtual void BaseSpecies::FinalizeMatterClassification()`, protected `bool clusters_as_matter_`, `bool is_cold_`.

- [ ] **Step 1: Write the failing test**

Create `species/composite_classification_test.cpp`:

```cpp
#include <cassert>
#include <cstdio>

#include "background.h"
#include "cdm.h"

int main() {
  background pba{};
  pba.H0 = 1e-4;

  // Cold dark matter: clusters as matter AND is cold. Cached values must match
  // the virtual predicates after FinalizeMatterClassification().
  CDMSpecies cdm(pba, 0.12, /*coupled=*/false);
  cdm.FinalizeMatterClassification();
  assert(cdm.ClustersAsMatterCached() == cdm.ClustersAsMatter());
  assert(cdm.IsColdCached() == cdm.IsColdMatterSpecies());
  assert(cdm.ClustersAsMatterCached() == true);
  assert(cdm.IsColdCached() == true);

  std::printf("composite classification test passed\n");
  return 0;
}
```

- [ ] **Step 2: Register the test and run it to verify it fails to build**

In `CMakeLists.txt`, add the executable after line 217 (`test-cdm-coupled`):

```cmake
  add_executable(test-composite-classification species/composite_classification_test.cpp)
```

And add `test-composite-classification` to the `foreach(... IN ITEMS ...)` list at line 219.

Run:
```bash
cmake -S . -B build/cmake && cmake --build build/cmake --target test-composite-classification 2>&1 | tail -20
```
Expected: FAIL — `'ClustersAsMatterCached' is not a member of 'BaseSpecies'` and `'FinalizeMatterClassification' was not declared`.

- [ ] **Step 3: Add the members, accessors, and finalize method to `base_species.h`**

In the protected block of `BaseSpecies` (near `species/base_species.h:606-614`, alongside `collection_index_`), add:

```cpp
  // Cached matter-tally classification, stamped once by FinalizeMatterClassification()
  // (driven by SpeciesCollection::freeze(); composites recurse into children). Read on
  // the stress-energy hot path, so it must be a plain member, not a per-step virtual.
  bool clusters_as_matter_ = false;
  bool is_cold_            = false;
```

In the public section (e.g. just before the `ClustersAsMatter()` doc-comment at line 559), add:

```cpp
  /** Cached counterparts of ClustersAsMatter()/IsColdMatterSpecies(), valid after
   *  FinalizeMatterClassification(). The hot-path tally reads these. */
  bool ClustersAsMatterCached() const { return clusters_as_matter_; }
  bool IsColdCached() const { return is_cold_; }

  /** Stamp the cached classification from the (virtual) predicates. Called once by
   *  SpeciesCollection::freeze() after every species and its children are built.
   *  CompositeSpecies overrides to also recurse into children. */
  virtual void FinalizeMatterClassification() {
    clusters_as_matter_ = ClustersAsMatter();
    is_cold_            = IsColdMatterSpecies();
  }
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```bash
cmake --build build/cmake --target test-composite-classification && ./build/cmake/test-composite-classification
```
Expected: PASS — `composite classification test passed`.

- [ ] **Step 5: Add the composite recursion override**

In `species/composite_species.h`, in the `// ── Matter tally ──` block (after line 133), declare:

```cpp
  void FinalizeMatterClassification() override;
```

In `species/composite_species.cpp` (after the `IsColdMatterSpecies` impl, ~line 108), add:

```cpp
void CompositeSpecies::FinalizeMatterClassification() {
  BaseSpecies::FinalizeMatterClassification();   // stamp this composite's own (unused but valid) bits
  for (auto& child : children_)
    child->FinalizeMatterClassification();        // children's bits drive DelegateTally (Task 3)
}
```

- [ ] **Step 6: Wire the finalize pass into `freeze()`**

In `species/species_collection.cpp`, immediately after the `collection_index_` stamping loop (after line 34), add:

```cpp
  /* Stamp the cached matter classification once, now that every species (and every
     composite child) is fully constructed. Composites recurse into their children. */
  for (std::size_t i = 0; i < species_.size(); ++i)
    species_[i].species->FinalizeMatterClassification();
```

- [ ] **Step 7: Build everything and run the full unit-test suite**

Run:
```bash
cmake --build build/cmake --parallel 2>&1 | tail -5 && ctest --test-dir build/cmake --output-on-failure
```
Expected: build succeeds; all tests PASS (including `test-composite-classification`).

- [ ] **Step 8: Commit**

```bash
git add species/base_species.h species/composite_species.h species/composite_species.cpp \
        species/species_collection.cpp species/composite_classification_test.cpp CMakeLists.txt
git commit -m "species: cache matter classification, stamped at freeze (recursing composites)"
```

---

## Task 2: Generic composite child iteration (Part A)

Give `CompositeSpecies` an owning `child_layouts` vector aligned with `children_`, a generic `CreatePerturbLayout`, and generic forwarding defaults; migrate `Type3Species` to drop its identical boilerplate. **Pure refactor — no behavior change.** Verified earlier: DCDM_DR / DNCDM_DR / IDM_DR_IDR all override every one of these 7 methods, so the new generic defaults are reached only by Type3.

**Files:**
- Modify: `species/composite_species.h` (`PerturbLayout` + `CreatePerturbLayout` + 7 forwarder decls)
- Modify: `species/composite_species.cpp` (7 forwarder impls + `CreatePerturbLayout`)
- Modify: `species/type3_species.h` (delete `PerturbLayout`/`CreatePerturbLayout` + 8 overrides; add `ChildIndex` + cast helpers)
- Modify: `species/type3_species.cpp` (delete the now-generic method bodies; rewrite `PrintVariables`)
- Create: `species/composite_layout_test.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: `BaseSpecies::PerturbLayout`, `BaseSpecies::CreatePerturbLayout()`, the per-method signatures already on `BaseSpecies`.
- Produces: `struct CompositeSpecies::PerturbLayout : BaseSpecies::PerturbLayout { std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> child_layouts; }`; generic `CompositeSpecies::{CreatePerturbLayout, RegisterPerturbationIndices, ApplyInitialConditions, StressEnergy, PerturbDerivs, FillSources, PerturbSynchronousToNewtonian, CopyPerturbationsAcrossSwitch}`; `Type3Species::ChildIndex{kCdm=0,kScf=1}` + `cdm_layout()`/`scf_layout()`.

- [ ] **Step 1: Write the failing test**

Create `species/composite_layout_test.cpp`:

```cpp
#include <cassert>
#include <cstdio>
#include <memory>
#include <vector>

#include "background.h"
#include "cdm.h"
#include "composite_species.h"
#include "scalar_field.h"
#include "type3_species.h"

int main() {
  background pba{};
  pba.H0 = 1e-4;
  const std::vector<double> params = {1.22, 0.0, 0.0, 0.0};
  auto scf = std::make_unique<ScalarFieldSpecies>(pba, /*omega0=*/0.7, params, /*tuning=*/0,
                                                  /*attractor=*/true, /*phi_ini=*/1.,
                                                  /*phi_prime_ini=*/1.,
                                                  DefaultScalarFieldPotential(), /*beta=*/-2.0);
  Type3Species t3(pba, /*omega0_cdm=*/0.12, std::move(scf));

  auto layout = t3.CreatePerturbLayout();
  auto& comp  = static_cast<CompositeSpecies::PerturbLayout&>(*layout);

  // Generic CreatePerturbLayout builds one sub-layout per child, in children_ order.
  assert(comp.child_layouts.size() == 2);
  assert(dynamic_cast<CDMSpecies::PerturbLayout*>(comp.child_layouts[0].get()) != nullptr);
  assert(dynamic_cast<ScalarFieldSpecies::PerturbLayout*>(comp.child_layouts[1].get()) != nullptr);

  std::printf("composite layout test passed\n");
  return 0;
}
```

- [ ] **Step 2: Register and run to verify it fails**

In `CMakeLists.txt` add after the Task-1 executable:
```cmake
  add_executable(test-composite-layout species/composite_layout_test.cpp)
```
and add `test-composite-layout` to the `foreach(... IN ITEMS ...)` list.

Run:
```bash
cmake -S . -B build/cmake && cmake --build build/cmake --target test-composite-layout 2>&1 | tail -20
```
Expected: FAIL — `CompositeSpecies::PerturbLayout` has no member `child_layouts` (it doesn't exist yet).

- [ ] **Step 3: Add `CompositeSpecies::PerturbLayout` + generic `CreatePerturbLayout`**

In `species/composite_species.h`, inside `class CompositeSpecies`, public section (after the ctor at line 29), add:

```cpp
  /** Generic composite layout: one owning sub-layout per child, aligned 1:1 with
   *  children_ (construction order). Concrete composites that need a typed child
   *  field cast child_layouts[k] to the concrete child layout. */
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> child_layouts;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override;
```

And declare the 7 generic forwarders (in the registration / perturbation section):

```cpp
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout, perturb_vector* pv,
                                   const precision* ppr, int& index_pt,
                                   const perturb_workspace* ppw, int gauge) override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout, double* y,
                              const PerturbIcContext& ctx) override;
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv, const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout, double tau, const double* y,
                     double* dy, const perturb_parameters_and_workspace& ppaw) const override;
  void FillSources(const BaseSpecies::PerturbLayout& layout, const double* y, const double* dy,
                   PerturbSourceContext& ctx) const override;
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout, double* y,
                                     const PerturbIcContext& ctx) override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y, double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 4: Implement the generic forwarders in `composite_species.cpp`**

At the top of `species/composite_species.cpp`, after the existing include, add the headers whose types are dereferenced/constructed here:

```cpp
#include "perturbations.h"  // perturb_vector, perturb_workspace, contexts, precision
```

Then append the implementations:

```cpp
std::unique_ptr<BaseSpecies::PerturbLayout> CompositeSpecies::CreatePerturbLayout() const {
  auto l = std::make_unique<PerturbLayout>();
  l->child_layouts.reserve(children_.size());
  for (const auto& c : children_)
    l->child_layouts.push_back(c->CreatePerturbLayout());
  return l;
}

void CompositeSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                   perturb_vector* pv, const precision* ppr,
                                                   int& index_pt, const perturb_workspace* ppw,
                                                   int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->RegisterPerturbationIndices(*my.child_layouts[i], pv, ppr, index_pt, ppw, gauge);
}

void CompositeSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base, double* y,
                                              const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->ApplyInitialConditions(*my.child_layouts[i], y, ctx);
}

BaseSpecies::StressEnergyContribution CompositeSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base, const perturb_vector* pv, const double* y,
    const double* pvecback, const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  for (size_t i = 0; i < children_.size(); ++i)
    se += children_[i]->StressEnergy(*my.child_layouts[i], pv, y, pvecback, ppw);
  return se;
}

void CompositeSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base, double tau,
                                     const double* y, double* dy,
                                     const perturb_parameters_and_workspace& ppaw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->PerturbDerivs(*my.child_layouts[i], tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);   // two-phase contract: children first, then coupling
}

void CompositeSpecies::FillSources(const BaseSpecies::PerturbLayout& base, const double* y,
                                   const double* dy, PerturbSourceContext& ctx) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->FillSources(*my.child_layouts[i], y, dy, ctx);
}

void CompositeSpecies::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                     double* y, const PerturbIcContext& ctx) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->PerturbSynchronousToNewtonian(*my.child_layouts[i], y, ctx);
}

void CompositeSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                     const BaseSpecies::PerturbLayout& new_base,
                                                     const double* old_y, double* new_y,
                                                     const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->CopyPerturbationsAcrossSwitch(*old_l.child_layouts[i], *new_l.child_layouts[i],
                                                old_y, new_y, ctx);
}
```

- [ ] **Step 5: Strip Type3's boilerplate and add cast helpers (`type3_species.h`)**

In `species/type3_species.h`:

Delete the nested `struct PerturbLayout { CDMSpecies::PerturbLayout cdm; ScalarFieldSpecies::PerturbLayout scf; };` (lines 22-25) and the `CreatePerturbLayout` override (lines 27-29). Replace with the child-index enum and typed cast helpers:

```cpp
  enum ChildIndex { kCdm = 0, kScf = 1 };  // children_ order, set in the ctor
  static const CDMSpecies::PerturbLayout& cdm_layout(const CompositeSpecies::PerturbLayout& my) {
    return static_cast<const CDMSpecies::PerturbLayout&>(*my.child_layouts[kCdm]);
  }
  static const ScalarFieldSpecies::PerturbLayout& scf_layout(const CompositeSpecies::PerturbLayout& my) {
    return static_cast<const ScalarFieldSpecies::PerturbLayout&>(*my.child_layouts[kScf]);
  }
```

Delete these now-generic override declarations (they fall through to `CompositeSpecies`):
`RegisterPerturbationIndices` (44-49), `PerturbDerivs` (50-54), `ApplyInitialConditions` (55-57), `StressEnergy` (58-62), `FillSources` (63-66), `RegisterTransferSourceIndices` (67 — `CompositeSpecies` already forwards it), `PerturbSynchronousToNewtonian` (68-70), `CopyPerturbationsAcrossSwitch` (71-75).

Keep: `WriteBackgroundColumnTitles/Data`, `WriteOutputColumns`, `PrintVariables`, the shooting forwarders, `kTypeName`, `CreateAll`, `AddCouplingDerivs`, and the `cdm_`/`scf_` typed species pointers.

- [ ] **Step 6: Delete the now-generic bodies and fix `PrintVariables` (`type3_species.cpp`)**

In `species/type3_species.cpp`, delete the bodies of `RegisterPerturbationIndices` (31-40), `PerturbDerivs` (42-51), `ApplyInitialConditions` (53-59), `StressEnergy` (61-71), `FillSources` (73-80), `RegisterTransferSourceIndices` (82-85), `PerturbSynchronousToNewtonian` (87-93), `CopyPerturbationsAcrossSwitch` (95-104).

In `PrintVariables` (still Type3-owned), replace the layout access at lines 137-139:

```cpp
    const auto& my_lay  = static_cast<const CompositeSpecies::PerturbLayout&>(*pv->species_layouts[collection_index_]);
    const auto& cdm_lay = cdm_layout(my_lay);
    const auto& scf_lay = scf_layout(my_lay);
```

(The rest of `PrintVariables` is unchanged — it reads `cdm_lay.idx_delta`, `scf_lay.idx_phi`, etc.)

- [ ] **Step 7: Build and run the layout test + full suite**

Run:
```bash
cmake --build build/cmake --parallel 2>&1 | tail -5 && ./build/cmake/test-composite-layout && ctest --test-dir build/cmake --output-on-failure
```
Expected: build succeeds; `composite layout test passed`; all tests PASS.

- [ ] **Step 8: Integration smoke — refactor is behavior-preserving**

Build `class` and run an existing scalar-field scenario; it must still run cleanly (this is a pure refactor):
```bash
make class && ./class test/scenarios/gauge_scf.ini 2>&1 | tail -15
```
Expected: completes without error, writes output as before. (Full numeric A/B vs master happens in Task 4.)

- [ ] **Step 9: Commit**

```bash
git add species/composite_species.h species/composite_species.cpp species/type3_species.h \
        species/type3_species.cpp species/composite_layout_test.cpp CMakeLists.txt
git commit -m "composite: generic child_layouts + forwarders; Type3 sheds boilerplate"
```

---

## Task 3: Per-child matter tally + honest density (Parts B + C)

Add the non-virtual `TallyStressEnergy` wrapper + virtual `DelegateTally` seam, switch the module's scalar stress-energy loop to use them, and drop the `ρ−3P`/`δρ−3δP` proxy in favor of each clustering species' actual `ρ`/`δρ`/`(ρ+P)θ`. **This is the behavior-change task:** Type3's matter tally becomes CDM-only (fix); warm NCDM/DNCDM `delta_m`/`theta_m` shift to the full density (intentional). `delta_cb`/`theta_cb` and cold-only cosmologies are unchanged (cold species are pressureless).

**Files:**
- Modify: `species/base_species.h` (`delegates_tally_` already added in Task 1's protected block — add the `TallyStressEnergy` inline + `DelegateTally` virtual)
- Modify: `species/composite_species.h` (`DelegateTally` decl; ctor sets `delegates_tally_`)
- Modify: `species/composite_species.cpp` (`DelegateTally` impl; ctor body)
- Modify: `source/perturbations.h:260-266` (drop the two `ActiveSpecies` bools)
- Modify: `source/perturbations_module.cpp` (3 `active_species` construction sites: 2993-3000, 3054-3057, 3087-3090; the scalar tally loop: 4833-4894)

**Interfaces:**
- Consumes: `BaseSpecies::clusters_as_matter_`/`is_cold_`/`delegates_tally_` (Task 1), `CompositeSpecies::PerturbLayout::child_layouts` (Task 2), `BaseSpecies::StressEnergyContribution` + `operator+=`.
- Produces: `void BaseSpecies::TallyStressEnergy(const PerturbLayout&, const perturb_vector*, const double* y, const double* pvecback, const perturb_workspace*, StressEnergyContribution& total, StressEnergyContribution& total_cold, StressEnergyContribution& total_warm) const` (non-virtual, inline); `virtual void BaseSpecies::DelegateTally(...same args...) const`; `CompositeSpecies::DelegateTally` override.

- [ ] **Step 1: Add `TallyStressEnergy` (inline) + `DelegateTally` to `base_species.h`**

In `species/base_species.h`, immediately after the `StressEnergy` pure-virtual declaration (after line 356), add. **Note:** the body must only *pass* `ppw` through — never dereference it — so this stays in `base_species.h` with `perturb_workspace` forward-declared (do NOT include `perturbations.h`):

```cpp
  // ── Stress-energy + matter tally (hot path) ───────────────────────────────
  // Non-virtual wrapper, inlined at the loop call site (e.species is BaseSpecies*).
  // Plain species: one virtual call (StressEnergy) + accumulate. Composites set
  // delegates_tally_ and forward the WHOLE call per child via DelegateTally, so each
  // child is evaluated exactly once. The three accumulators are owned by the module
  // loop; StressEnergyContribution::operator+= is reused as the bucket combiner.
  // ppw is only passed through to StressEnergy (never dereferenced), so this compiles
  // with perturb_workspace forward-declared.
  void TallyStressEnergy(const PerturbLayout& layout, const perturb_vector* pv,
                         const double* y, const double* pvecback, const perturb_workspace* ppw,
                         StressEnergyContribution& total,
                         StressEnergyContribution& total_cold,
                         StressEnergyContribution& total_warm) const {
    if (delegates_tally_) {
      DelegateTally(layout, pv, y, pvecback, ppw, total, total_cold, total_warm);
      return;
    }
    const StressEnergyContribution se = StressEnergy(layout, pv, y, pvecback, ppw);
    total += se;
    // Actual ρ/δρ/(ρ+P)θ — no ρ−3P proxy. Radiation never reaches here (it does not
    // cluster as matter), so there is nothing to "zero out".
    if (clusters_as_matter_)
      (is_cold_ ? total_cold : total_warm) += se;
  }

  // Seam: only composites (delegates_tally_ == true) override this; plain species
  // never reach the base body.
  virtual void DelegateTally(const PerturbLayout& /*layout*/, const perturb_vector* /*pv*/,
                             const double* /*y*/, const double* /*pvecback*/,
                             const perturb_workspace* /*ppw*/,
                             StressEnergyContribution& /*total*/,
                             StressEnergyContribution& /*total_cold*/,
                             StressEnergyContribution& /*total_warm*/) const {}
```

Also add the `delegates_tally_` member to the protected block (next to `clusters_as_matter_` from Task 1):

```cpp
  bool delegates_tally_ = false;  // true for composites (set in CompositeSpecies ctor)
```

- [ ] **Step 2: Composite overrides `DelegateTally` and sets the flag**

In `species/composite_species.h`, declare in the public section:

```cpp
  void DelegateTally(const BaseSpecies::PerturbLayout& layout, const perturb_vector* pv,
                     const double* y, const double* pvecback, const perturb_workspace* ppw,
                     StressEnergyContribution& total, StressEnergyContribution& total_cold,
                     StressEnergyContribution& total_warm) const override;
```

Change the ctor (line 28-29) to set the flag:

```cpp
  CompositeSpecies(std::string name, EnergyType energy_type)
      : BaseSpecies(std::move(name), energy_type) {
    delegates_tally_ = true;
  }
```

In `species/composite_species.cpp`, add:

```cpp
void CompositeSpecies::DelegateTally(const BaseSpecies::PerturbLayout& base, const perturb_vector* pv,
                                     const double* y, const double* pvecback,
                                     const perturb_workspace* ppw,
                                     StressEnergyContribution& total,
                                     StressEnergyContribution& total_cold,
                                     StressEnergyContribution& total_warm) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->TallyStressEnergy(*my.child_layouts[i], pv, y, pvecback, ppw,
                                    total, total_cold, total_warm);
}
```

- [ ] **Step 3: Drop the two booleans from `ActiveSpecies`**

In `source/perturbations.h:260-266`, change to:

```cpp
  struct ActiveSpecies {
    const BaseSpecies* species;
    const BaseSpecies::PerturbLayout* layout;
  };
```

- [ ] **Step 4: Update the three `active_species` construction sites**

In `source/perturbations_module.cpp`:

Scalar (lines 2993-2996):
```cpp
          perturb_vector::ActiveSpecies active{entry.get(),
                                               ppv->species_layouts[i].get()};
```

Vector (lines 3054-3057):
```cpp
          ppv->active_species.push_back({entry.get(),
                                         ppv->species_layouts[i].get()});
```

Tensor (lines 3087-3090):
```cpp
          ppv->active_species.push_back({entry.get(),
                                         ppv->species_layouts[i].get()});
```

Also update the two derivs structured bindings that name the removed fields (`perturbations_module.cpp:5871, 5931, 5947`) from `for (const auto& [species, layout, clusters_as_matter, is_cold] : pv->active_species)` to `for (const auto& [species, layout] : pv->active_species)`.

- [ ] **Step 5: Rewrite the scalar tally loop**

In `source/perturbations_module.cpp`, replace the block from line 4833 (`ppw->delta_rho = 0.;`) through line 4894 (the closing of the `_scalars_` matter-output section, i.e. the old `struct Tally`, the single-pass loop, and the post-loop `delta_cb`/`theta_cb`/`delta_m`/`theta_m`) with:

```cpp
    StressEnergyContribution total, total_cold, total_warm;   // module-owned accumulators

    /* Single pass: each species (or composite, via DelegateTally over its children)
       accumulates its stress-energy into `total` and folds its matter contribution
       into the cold/warm buckets (per-child, using the species' cached classification).
       A canonical scalar field's Newtonian StressEnergy reconstructs psi from the
       running shear in ppw->rho_plus_p_shear (active_species orders scalar fields last,
       see perturb_vector_init), so we keep that field synced from `total` after each
       species. See docs/superpowers/specs/2026-06-27-composite-child-iteration-and-matter-tally-design.md */
    ppw->rho_plus_p_shear = 0.;
    for (const auto& e : pv->active_species) {
      e.species->TallyStressEnergy(*e.layout, pv, y, pb, ppw, total, total_cold, total_warm);
      ppw->rho_plus_p_shear = total.rho_plus_p_shear;
    }

    ppw->delta_rho        = total.delta_rho;
    ppw->rho_plus_p_theta = total.rho_plus_p_theta;
    ppw->delta_p          = total.delta_p;
    ppw->rho_plus_p_tot   = total.rho + total.p;
    // ppw->rho_plus_p_shear already holds total.rho_plus_p_shear (last sync above).

    if (has_source_delta_m_ && has_source_delta_cb_)
      ppw->delta_cb = total_cold.delta_rho / total_cold.rho;
    if ((has_source_delta_m_ || has_source_theta_m_) &&
        (has_source_delta_cb_ || has_source_theta_cb_))
      ppw->theta_cb = total_cold.rho_plus_p_theta / (total_cold.rho + total_cold.p);

    if (auto* f = ppf_fluid()) {
      f->ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
      ppw->delta_rho        += ppw->delta_rho_fld;
      ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
      ppw->delta_p          += ppw->delta_p_fld;
      ppw->rho_plus_p_tot   += f->Rho(pb) + f->P(pb);
    }

    if (has_source_delta_m_)
      ppw->delta_m = (total_cold.delta_rho + total_warm.delta_rho) /
                     (total_cold.rho + total_warm.rho);
    if (has_source_delta_m_ || has_source_theta_m_)
      ppw->theta_m = (total_cold.rho_plus_p_theta + total_warm.rho_plus_p_theta) /
                     ((total_cold.rho + total_cold.p) + (total_warm.rho + total_warm.p));
```

Keep the `scalar_ctx` population block above it (lines 4826-4831) untouched.

- [ ] **Step 6: Build and run the full unit-test suite**

Run:
```bash
cmake --build build/cmake --parallel 2>&1 | tail -10 && ctest --test-dir build/cmake --output-on-failure
```
Expected: build succeeds; all unit tests PASS.

- [ ] **Step 7: Integration — cold cosmology is byte-stable (modulo ULP)**

Build `class`, run a cold scenario, and diff against a master baseline captured before this branch's changes. Generate the master baseline once from a clean `master` checkout into `python/baseline_ref/` per the classyref workflow ([[reference_classyref_testing]]), then:

```bash
make class
./class test/scenarios/gauge_lcdm.ini
python3 test/scenarios/compare_tol.py <master_lcdm_output> <new_lcdm_output> --rtol 1e-3
```
Expected: TT and mP(k) agree to ≤0.1% (cold species are pressureless, so Part C is a no-op for them).

- [ ] **Step 8: Integration — Type3 fix + NCDM shift are as designed**

Create a Type3 scenario from the existing scalar-field one plus the coupling and matter output. Copy `test/scenarios/gauge_scf.ini` to `test/scenarios/type3_scf_veta.ini`, set `gauge = synchronous`, add `scf_veta = 0.2`, and ensure `output = mPk,mTk` with `z_pk = 0`. Then:

```bash
git add -f test/scenarios/type3_scf_veta.ini
./class test/scenarios/type3_scf_veta.ini 2>&1 | tail -10
```
Expected: completes; the cold matter tally now excludes the scalar field. Sanity-check that `delta_cb` at z=0 equals the CDM child's `δρ_cdm/ρ_cdm` (the scalar field no longer pollutes it) by inspecting the `mTk` output columns.

For the intended NCDM change, run `gauge_ncdm.ini` and confirm TT is ≤0.1% vs master (the tally does not feed the ODE) while `delta_m`/mP(k) shift in the expected direction (full δρ vs pressure-subtracted; ~`3w·f_ν`):
```bash
./class test/scenarios/gauge_ncdm.ini
python3 test/scenarios/compare_tol.py <master_ncdm_output> <new_ncdm_output> --rtol 1e-3 || echo "mPk diff expected for NCDM (intentional)"
```

- [ ] **Step 9: Commit**

```bash
git add species/base_species.h species/composite_species.h species/composite_species.cpp \
        source/perturbations.h source/perturbations_module.cpp
git add -f test/scenarios/type3_scf_veta.ini
git commit -m "perturbations: per-child matter tally via TallyStressEnergy; drop rho-3P proxy"
```

---

## Task 4: Verification sweep & reference regeneration

Run the full reference test suite, confirm the cold/no-change invariants and the intended warm/Type3 changes, and regenerate the committed references/goldens for the intentional output changes.

**Files:**
- Modify (regenerate): `python/baseline_ref/`, transfer/background goldens as needed.

- [ ] **Step 1: Full unit-test suite (Debug, asserts on)**

```bash
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=Debug && cmake --build build/debug --parallel
ctest --test-dir build/debug --output-on-failure
```
Expected: all PASS, no assertion failures (catches `child_layouts`/`children_` misalignment and any null layout).

- [ ] **Step 2: Reference comparison vs a fresh `master` build (`classyref`)**

Build the `classyref` module from a clean `master` worktree, rebuild this branch's wrapper, then run the comparison suite with reference comparison enabled ([[reference_classyref_testing]]):

```bash
make classy-pip-dev   # rebuild the classy wrapper for this branch (see build gotcha in progress.md)
cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest test_class.py -q 2>&1 | tail -40
```
Expected: cold/radiation/PPF/tensor scenarios within ≤0.1%; only the **warm NCDM/DNCDM** matter-output scenarios flag differences (intentional). Inspect any flagged plot in `python/faulty_figs/` and confirm it is a matter-power shift, not a CMB regression.

- [ ] **Step 3: Confirm the intended changes are matter-only, not CMB**

For each flagged warm scenario, confirm TT is within ≤0.1% (the tally feeds only transfer/P(k) sources) while mP(k) shifts. If any TT exceeds 0.1%, STOP — that indicates the totals (not just the tally) changed; re-check the `ppw->rho_plus_p_shear` sync in Task 3 Step 5.

- [ ] **Step 4: Regenerate goldens for the intentional output changes**

Regenerate the transfer/background column goldens and the `classyref` baseline for the changed warm-species scenarios:

```bash
python3 python/gen_transfer_golden.py
python3 python/gen_background_golden.py
```
Review the diff to confirm only warm-species `delta_m`/`theta_m`/mP(k) columns moved, then stage explicitly.

- [ ] **Step 5: Commit the regenerated references**

```bash
git add python/baseline_ref python/transfer_golden python/background_golden
git commit -m "test: regenerate references for honest warm-species matter tally"
```

- [ ] **Step 6: (Optional) Re-profile the hot path**

Confirm the plain-species hot path is unchanged (one indirect call per entry):
```bash
make class_profiled && ./class_profiled base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini
```
Expected: `perturb_total_stress_energy` cost in line with the pre-change profile (no extra virtual per plain species).

---

## Self-Review

**Spec coverage:** Part A → Task 2; Part B (`TallyStressEnergy`/`DelegateTally`, cached bools, `delegates_tally_`) → Tasks 1+3; Part C (drop `ρ−3P`) → Task 3 Step 5; `ActiveSpecies` shrink → Task 3 Steps 3-4; scalar-field shear read-back → Task 3 Step 5 (the `ppw->rho_plus_p_shear` sync); ordering note (Type3 synchronous-only) → no code, covered by the synchronous-only guard; verification/regeneration → Task 4. All spec sections map to a task.

**Placeholder scan:** no TBD/TODO; every code step shows full code; commands have expected outcomes. The two integration baselines (`<master_*_output>`) are explicitly produced via the documented classyref workflow rather than left vague.

**Type consistency:** `TallyStressEnergy`/`DelegateTally` use the identical 8-parameter signature (`layout, pv, y, pvecback, ppw, total, total_cold, total_warm`) in `base_species.h`, `composite_species.h/.cpp`, and the module loop. `CompositeSpecies::PerturbLayout::child_layouts` is referenced consistently in Task 2 (forwarders) and Task 3 (`DelegateTally`). `ChildIndex{kCdm,kScf}` + `cdm_layout()/scf_layout()` defined in Task 2 and used in `PrintVariables`. Cached accessors `ClustersAsMatterCached()`/`IsColdCached()` defined and used in Task 1.
