# Encapsulate Species Perturbation Layouts (and Remove `ncdm_id`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `perturb_vector`'s hardcoded per-species fields (and the integer `ncdm_id` flat arrays) with polymorphic per-species `PerturbLayout` subclasses owned by `perturb_vector`. Each species fully owns its perturbation layout; `perturb_vector` is reduced to a container.

**Architecture:** Each species defines a `PerturbLayout` subclass holding its indices/sizes. `perturb_vector` holds `std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> species_layouts` parallel to `all_species_` (lex-key order). Per-thread isolation comes for free because each OpenMP thread allocates its own `pv`. Composites use nested composition. After all species migrate, `ncdm_id` and the legacy flat arrays are deleted.

**Tech Stack:** C++17, `make class -j` build, `pytest` in `python/`.

**Reference:** Design spec at `docs/superpowers/specs/2026-05-07-encapsulate-species-layouts-design.md`.

---

## Strategy

Tasks fall into nine phases:

- **Phase 0** (Task 1): branch hygiene + scenario / baseline confirmation.
- **Phase A** (Tasks 2-9): infrastructure — `PerturbLayout` base, `species_layouts` on `perturb_vector`, `SpeciesCollection::operator[]` and `index_of`, layout-aware Register and read virtuals on `BaseSpecies`, module-side dual-call dispatch, Lambda pattern check.
- **Phase B** (Tasks 10-23): per-species migration in dual-write mode. One species per task. Each task adds the layout type, populates it, migrates readers.
- **Phase C** (Tasks 24-26): bookkeeping hooks — wire `CopyPerturbationsAcrossSwitch` + `MarkUsedInSources` into `perturb_vector_init`; delete the 8 NCDM switch-copy blocks and analogous UR/IDR/IDM_DRMD/Photons blocks; delete tensor-mode NCDM emission blocks.
- **Phase D** (Tasks 27-29): cleanup of legacy `pv->*` fields once all species have migrated.
- **Phase E** (Tasks 30-34): `ncdm_id` removal — fields, `pba->N_ncdm/Omega0_ncdm_tot/N_decay_dr`, `ncdm_species_sorted_`, `index_pt_psi0_ncdm1`, `RescaledNCDMPerturbations` signature.
- **Phase F** (Task 35): rewrite NCDM-family output column titles to instance-name based.
- **Phase G** (Task 36): rename DNCDM_DR composite + DR child to use the dncdm instance name.
- **Phase H** (Tasks 37-40): shooter hooks + per-species migration of DCDM_DR, DNCDM_DR, ScalarField shooting cases.
- **Phase I** (Task 41): final audit + open PR.

Total: 41 tasks.

---

## Regression cycle (`<regression-check>`)

After every behaviour-touching task, run this exact sequence. **Do not skip steps.** The previous attempt failed because stale `.dat` files masked crashes:

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/remove-ncdm-id

# 1. Clean output of stale scenario .dat files (CRITICAL).
rm -f output/scen_*.dat

# 2. Build clean.
make class -j 2>&1 | tail -10

# 3. Run each scenario; explicitly check exit code per scenario.
set -e
for f in test/scenarios/ncdm_single.ini \
         test/scenarios/ncdm_multi_unsorted.ini \
         test/scenarios/ncdm_self_interacting.ini \
         test/scenarios/dncdm_dr.ini \
         test/scenarios/ncdm_dncdm_idmdr_combined.ini; do
  echo "===== $f ====="
  ./class "$f" || { echo "FAIL: $f exited non-zero"; exit 1; }
done
set +e

# 4. Smoke test explanatory.ini too.
./class explanatory.ini > /tmp/exp.log 2>&1 || { tail -30 /tmp/exp.log; exit 1; }

# 5. Diff against baseline.
python3 - <<'PY'
import pathlib, numpy as np, sys
base = pathlib.Path("/tmp/ncdm_baseline"); new = pathlib.Path("output")
ok = True
for f in sorted(base.glob("*.dat")):
    n = new / f.name
    if not n.exists(): print(f"MISSING: {f.name}"); ok = False; continue
    a = np.loadtxt(f, comments="#"); b = np.loadtxt(n, comments="#")
    if a.shape != b.shape: print(f"SHAPE: {f.name} base={a.shape} new={b.shape}"); ok = False; continue
    rel = np.where(np.abs(a) > 0, np.abs((a - b) / a), np.abs(a - b))
    rmax = float(rel.max()) if rel.size else 0.0
    print(f"{f.name}: max_rel_diff={rmax:.2e}")
    if rmax > 1e-10: ok = False
sys.exit(0 if ok else 1)
PY
```

Expected: all 23 baseline `.dat` files report `max_rel_diff=0.00e+00` (bit-identical) for all dual-write tasks. Tasks that knowingly break bit-identicality (column-title rename in Phase F, NCDM ordering in `ncdm_multi_unsorted` once flat arrays are gone) document the expected diff in the task description.

Default thread count (do NOT set `OMP_NUM_THREADS=1` for regression). Multi-thread coverage is the whole point.

---

## Task 1: Branch hygiene, scenario verification, baseline check

**Files:** none (verification only)

- [ ] **Step 1: Confirm branch and HEAD**

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/remove-ncdm-id
git branch --show-current      # expect: 267-remove-ncdm-id
git log --oneline -3
# Expected last 3 commits:
#   <sha> Spec: encapsulate species perturbation layouts (and remove ncdm_id)
#   2f64cffe test: scenarios for ncdm_id removal regression checking
#   07077010 Plan: remove ncdm_id from CLASSpp
```

- [ ] **Step 2: Confirm test scenarios are present**

```bash
ls test/scenarios/
```

Expected: 5 files — `dncdm_dr.ini`, `ncdm_dncdm_idmdr_combined.ini`, `ncdm_multi_unsorted.ini`, `ncdm_self_interacting.ini`, `ncdm_single.ini`.

- [ ] **Step 3: Confirm baseline directory**

```bash
ls /tmp/ncdm_baseline/ | wc -l
```

Expected: ~23 (5 logs + 18 .dat). If empty (e.g., on a fresh machine), regenerate by checking out master, building, and running each scenario:

```bash
cd /Users/au192734/Projects/class_claude
git branch --show-current  # expect: master
mkdir -p test/scenarios
cp .claude/worktrees/remove-ncdm-id/test/scenarios/*.ini test/scenarios/
test -x ./class || make class -j
mkdir -p /tmp/ncdm_baseline output
for ini in test/scenarios/ncdm_single.ini test/scenarios/ncdm_multi_unsorted.ini \
           test/scenarios/dncdm_dr.ini test/scenarios/ncdm_self_interacting.ini \
           test/scenarios/ncdm_dncdm_idmdr_combined.ini; do
  ./class "$ini" 2>&1 | tee "/tmp/ncdm_baseline/$(basename $ini .ini).log"
done
cp output/scen_*.dat /tmp/ncdm_baseline/
rm -rf test/scenarios/ output/scen_*.dat
```

- [ ] **Step 4: Run regression cycle from current HEAD as a smoke test**

Run the full `<regression-check>` block above. Expected: all scenarios pass with `max_rel_diff=0.00e+00`. This proves the worktree is in a clean state.

- [ ] **Step 5: No commit (verification only)**

---

## Task 2: Add `BaseSpecies::PerturbLayout` base type and `CreatePerturbLayout` virtual

**Files:**
- Modify: `species/base_species.h`

The base `PerturbLayout` is empty with a virtual destructor. `CreatePerturbLayout()` is a non-pure virtual returning a default-constructed layout — concrete species override it during migration. This means the build stays green even before any concrete species defines its own layout.

- [ ] **Step 1: Add `PerturbLayout` base struct and `CreatePerturbLayout` virtual to `BaseSpecies`**

Edit `species/base_species.h`. Locate the start of `class BaseSpecies` (currently around the `public:` after the constructor doc comment). Add this **at the very top of the public section**, before any other member:

```cpp
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
```

You may also need `#include <memory>` at the top of `species/base_species.h` if it's not already included. Check first.

- [ ] **Step 2: Build to confirm header changes compile**

```bash
make class -j 2>&1 | tail -10
```

Expected: clean build (no warnings/errors from your edit). `species/base_species.h` is included by every concrete species, so the build exercises every translation unit.

- [ ] **Step 3: Run `<regression-check>`**

Expected: bit-identical (`max_rel_diff=0.00e+00`) — the new infrastructure is unused.

- [ ] **Step 4: Commit**

```bash
git add species/base_species.h
git commit -m "$(cat <<'EOF'
species: add PerturbLayout base + CreatePerturbLayout virtual

Polymorphic base type for per-species perturbation index storage. Default
CreatePerturbLayout returns an empty base PerturbLayout; concrete species
will override during migration to per-species layouts.

Spec: docs/superpowers/specs/2026-05-07-encapsulate-species-layouts-design.md
EOF
)"
```

---

## Task 3: Add `SpeciesCollection::operator[](size_t)` for index-based dispatch

**Files:**
- Modify: `species/species_collection.h`

- [ ] **Step 1: Add `operator[](size_t)` to `SpeciesCollection`**

Edit `species/species_collection.h`. After the existing `size()` method (around line 126-128), add:

```cpp
  /** Index-based access for parallel-vector dispatch (e.g. with
   *  perturb_vector::species_layouts). i must be < size(). */
  BaseSpecies* operator[](std::size_t i) {
    assert(i < species_.size());
    return species_[i].species.get();
  }
  const BaseSpecies* operator[](std::size_t i) const {
    assert(i < species_.size());
    return species_[i].species.get();
  }
```

- [ ] **Step 2: Build**

```bash
make class -j 2>&1 | tail -5
```

Expected: clean.

- [ ] **Step 3: Run `<regression-check>`**

Bit-identical expected.

- [ ] **Step 4: Commit**

```bash
git add species/species_collection.h
git commit -m "species_collection: add operator[](size_t) for parallel-vector dispatch"
```

---

## Task 4: Add `species_layouts` to `perturb_vector` and populate it during `perturb_vector_init`

**Files:**
- Modify: `source/perturbations.h`
- Modify: `source/perturbations_module.cpp`

This task adds the new `species_layouts` storage on `perturb_vector` and populates it in `perturb_vector_init` parallel to all the existing legacy field initialisation. No species's behaviour changes yet — every species's `CreatePerturbLayout` still returns an empty base `PerturbLayout` (from Task 2).

- [ ] **Step 1: Add `species_layouts` to `perturb_vector`**

Edit `source/perturbations.h`. Inside `struct perturb_vector` (line 257-334), add **at the very top** (above `index_pt_delta_g`):

```cpp
  // Per-species perturbation layouts, parallel to all_species_ (lex-key order).
  // Each thread allocates its own pv, so this storage is per-thread.
  std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> species_layouts;
```

You'll need `#include <memory>` and `#include <vector>` if not already present, and `#include "../species/base_species.h"` for `BaseSpecies::PerturbLayout`. Check existing includes first; `perturbations.h` may already include some via transitive includes.

- [ ] **Step 2: Populate `species_layouts` at the top of `perturb_vector_init`**

Edit `source/perturbations_module.cpp`. Find `perturb_vector_init` (line ~3200, look for `int PerturbationsModule::perturb_vector_init(`). Immediately after `struct perturb_vector* ppv = new perturb_vector();`, add:

```cpp
  // Create per-species layouts, one per all_species_ entry, in iteration order.
  ppv->species_layouts.reserve(all_species_.size());
  for (size_t i = 0; i < all_species_.size(); ++i) {
    ppv->species_layouts.push_back(all_species_[i]->CreatePerturbLayout());
  }
```

(`operator[]` from Task 3 makes this clean.)

- [ ] **Step 3: Build**

```bash
make class -j 2>&1 | tail -10
```

Expected: clean. The new vector is populated but unused.

- [ ] **Step 4: Run `<regression-check>`**

Expected: bit-identical.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations.h source/perturbations_module.cpp
git commit -m "perturbations: add species_layouts on perturb_vector, populate in perturb_vector_init"
```

---

## Task 5: Add layout-based `Register*PerturbationIndices` virtuals on `BaseSpecies`

**Files:**
- Modify: `species/base_species.h`

This task adds the new layout-based virtuals as no-op defaults, *alongside* the existing `(perturb_vector*, ...)` virtuals. Both signatures coexist until Phase D deletes the old one. Migrated species override the new (layout-based) signature and reduce the legacy signature to a no-op; un-migrated species keep using the legacy signature; the module call site (Task 7) calls BOTH overloads, so combined effect is "migrated does the work via new; un-migrated does the work via legacy".

- [ ] **Step 1: Add new virtuals**

Edit `species/base_species.h`. Find each existing `Register{Scalar,Vector,Tensor}PerturbationIndices(perturb_vector*, ...)` virtual. Beside each, add a parallel virtual taking `BaseSpecies::PerturbLayout&` and the same `perturb_vector*` (the latter is needed for dual-write during migration; deleted in Phase D):

```cpp
  // ── New per-species layout signatures (in-progress migration) ────────────
  // Migrated species override these; un-migrated species ignore them and the
  // module call site passes through to the legacy (perturb_vector*) overload.

  virtual void RegisterPerturbationIndices(PerturbLayout& /*layout*/,
                                            perturb_vector* /*pv*/,
                                            const precision* /*ppr*/,
                                            int& /*index_pt*/,
                                            const perturb_workspace* /*ppw*/,
                                            int /*gauge*/) {}

  virtual void RegisterVectorPerturbationIndices(PerturbLayout& /*layout*/,
                                                  perturb_vector* /*pv*/,
                                                  int& /*index_pt*/,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {}

  virtual void RegisterTensorPerturbationIndices(PerturbLayout& /*layout*/,
                                                  perturb_vector* /*pv*/,
                                                  int& /*index_pt*/,
                                                  const perturb_workspace* /*ppw*/,
                                                  int /*gauge*/) {}
```

- [ ] **Step 2: Build, regression**

`make class -j` and `<regression-check>`. Expected: bit-identical (new virtuals unused).

- [ ] **Step 3: Commit**

```bash
git add species/base_species.h
git commit -m "species: add layout-based Register*PerturbationIndices virtuals (parallel to legacy)"
```

---

## Task 6: Add `SpeciesCollection::index_of(key)` helper

**Files:**
- Modify: `species/species_collection.h`

Used by `perturb_vector_init` to map a named species to its slot in `all_species_`/`species_layouts`.

- [ ] **Step 1: Add `index_of`**

Edit `species/species_collection.h`. After the existing `count(key)` / `at(key)` / `find(key)` cluster, add:

```cpp
  /** Lex-sorted index of an inserted species, or size() if absent. */
  std::size_t index_of(const std::string& key) const {
    auto it = std::lower_bound(species_.begin(), species_.end(), key,
                               [](const Entry& e, const std::string& k) { return e.key < k; });
    if (it == species_.end() || it->key != key) return species_.size();
    return static_cast<std::size_t>(it - species_.begin());
  }
```

You'll need `#include <algorithm>` if not already.

- [ ] **Step 2: Build + regression**

`<regression-check>`. Expected: bit-identical (helper unused).

- [ ] **Step 3: Commit**

```bash
git add species/species_collection.h
git commit -m "species_collection: add index_of(key) helper"
```

---

## Task 7: Module-side dual-call dispatch for `Register*PerturbationIndices`

**Files:**
- Modify: `source/perturbations_module.cpp`

Change every call site in `perturb_vector_init` that calls `sp->Register*PerturbationIndices(ppv, ...)` to call BOTH the new (layout) signature and the legacy signature. Each species defines exactly ONE of the two with real work and the other as no-op (override-and-empty). The combined effect: `index_pt` advances exactly once.

- [ ] **Step 1: Replace each Register* call with dual-call pattern**

Find every call:
```bash
grep -n "RegisterPerturbationIndices\|RegisterTensorPerturbationIndices\|RegisterVectorPerturbationIndices" source/perturbations_module.cpp | head -30
```

For each named-species call site like `all_species_.at("CDM")->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);`, replace with:
```cpp
{
  const size_t i = all_species_.index_of("CDM");
  if (i < all_species_.size()) {
    auto& layout = *ppv->species_layouts[i];
    all_species_[i]->RegisterPerturbationIndices(layout, ppv, ppr, index_pt, ppw, ppt->gauge);  // new (no-op until species migrated)
    all_species_[i]->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);          // legacy (no-op once species migrated)
  }
}
```

Apply uniformly to every named-species call site. There are ~10-15 such sites (Photons, Baryons, CDM, DCDM, DCDM_DR, Fluid, ScalarField, UR, IDR, IDM_DR_IDR, IDM_DRMD_IDR_DRMD, plus the NCDM-family loop and the tensor-mode block).

For NCDM-family inside `if (!ncdm_species_sorted_.empty()) { ... ncdm_reg_vec ... }`: replace the `ncdm_reg_vec` machinery with iteration over `all_species_` filtered by `dynamic_cast<NCDMBaseSpecies*>`:
```cpp
for (size_t i = 0; i < all_species_.size(); ++i) {
  if (auto* ncdm_sp = dynamic_cast<NCDMBaseSpecies*>(all_species_[i])) {
    auto& layout = *ppv->species_layouts[i];
    ncdm_sp->RegisterPerturbationIndices(layout, ppv, ppr, index_pt, ppw, ppt->gauge);
    ncdm_sp->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);
  }
}
```

(The `if (!ncdm_species_sorted_.empty())` outer guard can keep using `ncdm_species_sorted_` for now — it's deleted in Task 32.)

- [ ] **Step 2: Build + regression**

`<regression-check>`. Expected: bit-identical. Both overloads run; only legacy does work (no species migrated yet); index_pt advances once per species.

- [ ] **Step 3: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: dual-call Register* dispatch (layout + legacy)"
```

---

## Task 8: Migrate `Lambda` species (no-perturbation pattern establishment)

**Files:**
- Modify: `species/lambda.h`

Lambda (cosmological constant) has no perturbation slots — no `RegisterPerturbationIndices` override exists. This task adds the empty layout subclass. Subsequent species (CDM, Fluid, …) follow the same pattern but with non-empty layouts.

- [ ] **Step 1: Define `LambdaPerturbLayout` and override `CreatePerturbLayout`**

Edit `species/lambda.h`. Inside the class's public section:
```cpp
  // Lambda has no perturbation slots; the layout is empty.
  struct PerturbLayout : BaseSpecies::PerturbLayout {};

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }
```

(The nested `PerturbLayout` shadows the base name — fine since each species's layout is independent. This pattern repeats for every species.)

- [ ] **Step 2: Build + regression**

`<regression-check>`. Bit-identical.

- [ ] **Step 3: Commit**

```bash
git add species/lambda.h
git commit -m "lambda: add empty PerturbLayout"
```

---

## Task 9: Add layout-based variants for read virtuals

**Files:**
- Modify: `species/base_species.h`
- Modify: `species/perturb_source_context.h`

Add layout-aware overloads for `PerturbDerivs`, `PerturbVectorDerivs`, `PerturbTensorDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`, `CopyPerturbationsAcrossSwitch`, `MarkUsedInSources`. Defaults forward to the legacy (non-layout) variant where one exists, so un-migrated species keep working unchanged. Also defines `PerturbSwitchContext`.

- [ ] **Step 1: Define `PerturbSwitchContext`**

Edit `species/perturb_source_context.h`. After the existing `PerturbIcContext`, add:
```cpp
// ── PerturbSwitchContext ─────────────────────────────────────────────────────
struct PerturbSwitchContext {
  double k        = 0.;
  double a        = 0.;
  double a_today  = 1.;
  const double* pvecback = nullptr;
};
```

- [ ] **Step 2: Add layout-based virtuals on `BaseSpecies`**

Edit `species/base_species.h`. Place each new virtual immediately after the legacy one. For `Delta`:
```cpp
  virtual double Delta(const PerturbLayout& /*layout*/, const perturb_vector* pv,
                       const double* y, const double* pvecback,
                       const perturb_workspace* ppw) const {
    return Delta(pv, y, pvecback, ppw);  // forward to legacy by default
  }
```

Repeat for `Theta`, `DeltaP`, `RhoPlusPShear` (return double; same shape).

For `PerturbDerivs`:
```cpp
  virtual void PerturbDerivs(const PerturbLayout& /*layout*/, double tau,
                             const double* y, double* dy,
                             const perturb_parameters_and_workspace& ppaw) {
    PerturbDerivs(tau, y, dy, ppaw);  // forward to legacy by default
  }
```

Same shape for `PerturbVectorDerivs`, `PerturbTensorDerivs`.

For `ApplyInitialConditions(double* y, const PerturbIcContext& ctx)`:
```cpp
  virtual void ApplyInitialConditions(const PerturbLayout& /*layout*/, double* y,
                                       const PerturbIcContext& ctx) {
    ApplyInitialConditions(y, ctx);
  }
```

For `FillSources(const double* y, const double* dy, PerturbSourceContext& ctx)`:
```cpp
  virtual void FillSources(const PerturbLayout& /*layout*/, const double* y,
                           const double* dy, PerturbSourceContext& ctx) {
    FillSources(y, dy, ctx);
  }
```

`CopyPerturbationsAcrossSwitch` and `MarkUsedInSources` are NEW (no legacy):
```cpp
  virtual void CopyPerturbationsAcrossSwitch(
      const PerturbLayout& /*old_layout*/, const PerturbLayout& /*new_layout*/,
      const double* /*old_y*/, double* /*new_y*/,
      const PerturbSwitchContext& /*ctx*/) const {}
  virtual void MarkUsedInSources(
      const PerturbLayout& /*layout*/, int* /*used_in_sources*/) const {}
```

Add `#include "perturb_source_context.h"` if not present.

- [ ] **Step 3: Build + regression**

`<regression-check>`. Bit-identical (new virtuals forward to legacy by default).

- [ ] **Step 4: Commit**

```bash
git add species/base_species.h species/perturb_source_context.h
git commit -m "species: add layout-based variants for read virtuals + PerturbSwitchContext"
```

---

## Task 10: Migrate `FluidSpecies` to per-species layout

**Files:**
- Modify: `species/fluid.h`
- Modify: `species/fluid.cpp`
- Modify: `source/perturbations_module.cpp` (generic-loop dispatches)

Establishes the full migration pattern for a species with non-trivial perturbation slots. Subsequent species (CDM, DCDM, ...) follow the same shape.

- [ ] **Step 1: Inspect existing Fluid usage**

```bash
grep -n "index_pt_delta_fld\|index_pt_theta_fld\|index_pt_Gamma_fld" species/fluid.{h,cpp} source/perturbations_module.cpp
```

Expected fields: `pv->index_pt_delta_fld`, `pv->index_pt_theta_fld`, `pv->index_pt_Gamma_fld` (PPF dynamical variable). Read sites are inside `species/fluid.cpp` (PerturbDerivs, Delta, Theta, etc.) and possibly the existing `RegisterPerturbationIndices(perturb_vector*, ...)` in `species/fluid.cpp`. Module code only reads these via species methods.

- [ ] **Step 2: Define `FluidSpecies::PerturbLayout`**

Edit `species/fluid.h`. Add inside the class's public section:
```cpp
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    int idx_delta = -1;
    int idx_theta = -1;
    int idx_Gamma = -1;  // PPF dynamical variable (only populated when use_ppf == true)
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }
```

- [ ] **Step 3: Override new-signature `RegisterPerturbationIndices`; reduce legacy to no-op**

Edit `species/fluid.h`, declare both overrides:
```cpp
  // Layout signature: does the work (and dual-writes legacy fields).
  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                    perturb_vector* pv, const precision* ppr,
                                    int& index_pt, const perturb_workspace* ppw,
                                    int gauge) override;
  // Legacy signature: no-op (Fluid is migrated).
  void RegisterPerturbationIndices(perturb_vector* /*pv*/, const precision* /*ppr*/,
                                    int& /*index_pt*/, const perturb_workspace* /*ppw*/,
                                    int /*gauge*/) override {}
```

Edit `species/fluid.cpp`. Find the existing `RegisterPerturbationIndices(perturb_vector*, ...)` body and port its logic into the new method. Dual-write to BOTH `layout.idx_*` AND `pv->index_pt_*_fld`:

```cpp
void FluidSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                perturb_vector* pv,
                                                const precision* /*ppr*/,
                                                int& index_pt,
                                                const perturb_workspace* /*ppw*/,
                                                int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  if (use_ppf_) {
    layout.idx_Gamma       = index_pt;
    pv->index_pt_Gamma_fld = index_pt;
    ++index_pt;
  } else {
    layout.idx_delta       = index_pt;
    pv->index_pt_delta_fld = index_pt;
    ++index_pt;
    layout.idx_theta       = index_pt;
    pv->index_pt_theta_fld = index_pt;
    ++index_pt;
  }
}
```

(Adapt the exact PPF / standard-fluid branching from existing code — verify `use_ppf_` is the actual member name with `grep use_ppf_ species/fluid.{h,cpp}`.)

Delete the body of the legacy `RegisterPerturbationIndices(perturb_vector*, ...)` — it's now a header-only `{}` no-op declared above.

- [ ] **Step 4: Override layout-based read virtuals; reduce legacy variants to no-op**

For each of `PerturbDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`: declare both signatures in `species/fluid.h` (layout-aware override + legacy no-op):
```cpp
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout, double tau,
                     const double* y, double* dy,
                     const perturb_parameters_and_workspace& ppaw) override;
  void PerturbDerivs(double /*tau*/, const double* /*y*/, double* /*dy*/,
                     const perturb_parameters_and_workspace& /*ppaw*/) override {}

  double Delta(const BaseSpecies::PerturbLayout& layout, const perturb_vector* pv,
               const double* y, const double* pvecback,
               const perturb_workspace* ppw) const override;
  double Delta(const perturb_vector* /*pv*/, const double* /*y*/,
               const double* /*pvecback*/, const perturb_workspace* /*ppw*/) const override {
    return 0.;
  }
  // ... same shape for Theta, DeltaP, RhoPlusPShear (all return double)
  // ... and:
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout, double* y,
                               const PerturbIcContext& ctx) override;
  void ApplyInitialConditions(double* /*y*/, const PerturbIcContext& /*ctx*/) override {}

  void FillSources(const BaseSpecies::PerturbLayout& layout, const double* y,
                   const double* dy, PerturbSourceContext& ctx) override;
  void FillSources(const double* /*y*/, const double* /*dy*/,
                   PerturbSourceContext& /*ctx*/) override {}
```

In `species/fluid.cpp`, port each method body:
```cpp
void FluidSpecies::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                  double tau, const double* y, double* dy,
                                  const perturb_parameters_and_workspace& ppaw) {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  // ... port the existing body, replacing every pv->index_pt_delta_fld with
  //     layout.idx_delta, pv->index_pt_theta_fld with layout.idx_theta,
  //     pv->index_pt_Gamma_fld with layout.idx_Gamma.
}
```

Apply the same substitution pattern to `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`.

- [ ] **Step 5: Convert module-side generic dispatch loops to index-based**

Find generic loops in `source/perturbations_module.cpp` that iterate `all_species_` calling `PerturbDerivs`, `Delta`, etc.:
```bash
grep -nE 'for.*\[.*sp\].*all_species_.*\{[^}]*sp->(PerturbDerivs|Delta|Theta|DeltaP|RhoPlusPShear|ApplyInitialConditions|FillSources)' source/perturbations_module.cpp
```

For each match, convert to index-based:
```cpp
// Before:
for (auto& [name, sp] : all_species_) sp->PerturbDerivs(tau, y, dy, *pppaw);

// After:
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->PerturbDerivs(*ppw->pv->species_layouts[i], tau, y, dy, *pppaw);
}
```

For un-migrated species the layout-based default forwards to legacy (Task 9), so behaviour is unchanged. For migrated species (Fluid after this task) the layout override does the work; legacy is no-op.

There are roughly 5-10 such generic loops. Convert all of them.

- [ ] **Step 6: Build + regression**

`<regression-check>`. Bit-identical.

- [ ] **Step 7: Commit**

```bash
git add species/fluid.h species/fluid.cpp source/perturbations_module.cpp
git commit -m "fluid: migrate to per-species PerturbLayout; module dispatch via species_layouts"
```

---

## Task 11: Migrate `CDMSpecies` (single-instance, simple)

**Files:**
- Modify: `species/cdm.h`
- Modify: `species/cdm.cpp`

CDM has `pv->index_pt_delta_cdm` (and `pv->index_pt_theta_cdm` only in newtonian gauge — synchronous gauge uses metric).

- [ ] **Step 1: Define `CDMSpecies::PerturbLayout`**

```cpp
// species/cdm.h, in public section
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1;
  int idx_theta = -1;  // newtonian gauge only
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 2: Override layout-based Register, reduce legacy to no-op**

Port the body of the existing `RegisterPerturbationIndices(perturb_vector*, ...)`. Dual-write `layout.idx_*` AND `pv->index_pt_*`. Make the legacy overload a no-op.

- [ ] **Step 3: Override layout-based read virtuals (PerturbDerivs, Delta, Theta, DeltaP, RhoPlusPShear, ApplyInitialConditions, FillSources)**

Port each body to use `layout.idx_delta` / `idx_theta` instead of `pv->index_pt_delta_cdm` / `theta_cdm`. Reduce legacy overloads to no-op (or forward-to-zero for return-doubles).

- [ ] **Step 4: Build + regression**

Bit-identical.

- [ ] **Step 5: Commit**

```bash
git add species/cdm.h species/cdm.cpp
git commit -m "cdm: migrate to per-species PerturbLayout"
```

---

## Task 12: Migrate `DCDMSpecies`

**Files:**
- Modify: `species/dcdm.h`
- Modify: `species/dcdm.cpp`

DCDM has `pv->index_pt_delta_dcdm`, `pv->index_pt_theta_dcdm`. Same shape as CDM (single density + velocity).

- [ ] **Step 1: Define `DCDMSpecies::PerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1;
  int idx_theta = -1;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 2-4: Mirror Task 11 — Register, read methods, regression**

Port `pv->index_pt_delta_dcdm/theta_dcdm` reads to `layout.idx_delta/theta`.

- [ ] **Step 5: Commit**

```bash
git add species/dcdm.h species/dcdm.cpp
git commit -m "dcdm: migrate to per-species PerturbLayout"
```

---

## Task 13: Migrate `DRSpecies` (Dark Radiation)

**Files:**
- Modify: `species/dark_radiation_species.h`
- Modify: `species/dark_radiation_species.cpp`

DR has multipole hierarchy: `pv->index_pt_F0_dr_sum`, `pv->index_pt_F0_dr_species`, `l_max_dr`, `l_max_dr_col`. Inspect first:

```bash
grep -n "index_pt_F0_dr\|l_max_dr" species/dark_radiation_species.{h,cpp}
```

- [ ] **Step 1: Define `DRPerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_F0_sum     = -1;
  int idx_F0_species = -1;
  int l_max          = -1;
  int l_max_col      = -1;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 2-4: Standard migration**

Port reads. Keep dual-write to `pv->index_pt_F0_dr_*` and `pv->l_max_dr_*`.

- [ ] **Step 5: Commit**

```bash
git add species/dark_radiation_species.h species/dark_radiation_species.cpp
git commit -m "dr: migrate to per-species PerturbLayout"
```

---

## Task 14: Migrate `ScalarFieldSpecies`

**Files:**
- Modify: `species/scalar_field.h`
- Modify: `species/scalar_field.cpp`

ScalarField has `pv->index_pt_phi_scf`, `pv->index_pt_phi_prime_scf`.

- [ ] **Step 1: Define `ScalarFieldPerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_phi       = -1;
  int idx_phi_prime = -1;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 2-4: Standard migration**

Port `pv->index_pt_phi_scf` → `layout.idx_phi`, `pv->index_pt_phi_prime_scf` → `layout.idx_phi_prime`.

- [ ] **Step 5: Commit**

```bash
git add species/scalar_field.h species/scalar_field.cpp
git commit -m "scalar_field: migrate to per-species PerturbLayout"
```

---

## Task 15: Migrate `URSpecies` (Ultra-Relativistic)

**Files:**
- Modify: `species/ultra_relativistic.cpp`
- Modify: `species/ultra_relativistic.h` (if it exists; otherwise inline header inside .cpp)

UR has `pv->index_pt_delta_ur/theta_ur/shear_ur/l3_ur` and `pv->l_max_ur`. Plus `ufa` switch behaviour.

- [ ] **Step 1: Locate the UR class and headers**

```bash
ls species/ultra_relativistic.* species/ur.* species/all_species.h
grep -n "URSpecies\|class.*UR" species/*.h species/*.cpp
```

- [ ] **Step 2: Define `URPerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1;
  int idx_theta = -1;
  int idx_shear = -1;
  int idx_l3    = -1;
  int l_max     = -1;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 3: Override layout-based Register, dual-write**

Port the existing `Register*PerturbationIndices(perturb_vector*, ...)` body. Set `layout.idx_*` AND `pv->index_pt_*_ur` (and `pv->l_max_ur`). The `ufa_on` branch collapses the hierarchy similarly to how NCDM's fa_on does — port that branch correctly into the layout sets (e.g., when `ufa_on`, `layout.l_max = 2`).

- [ ] **Step 4: Override read virtuals**

Port `PerturbDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`.

- [ ] **Step 5: Override `CopyPerturbationsAcrossSwitch` for UR's ufa transition**

Currently the UR ufa switch-copy block in `perturb_vector_init` is the slot-by-slot ppv->y[delta_ur] = ppw->pv->y[delta_ur] copy. Move that into:

```cpp
void URSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y, double* new_y,
    const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0) return;
  // Slot-by-slot copy of valid slots.
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
  new_y[new_l.idx_shear] = old_y[old_l.idx_shear];
  // l3 only present when ufa is off (full hierarchy):
  if (old_l.l_max >= 3 && new_l.l_max >= 3) {
    new_y[new_l.idx_l3] = old_y[old_l.idx_l3];
    for (int l = 4; l <= new_l.l_max; ++l)
      new_y[new_l.idx_delta + l] = old_y[old_l.idx_delta + l];
  }
}
```

- [ ] **Step 6: Build + regression**

Bit-identical (UR ufa switch-copy is now done by the species hook AND the legacy module-side block; both run, both write the same slots, so y values are unchanged. The legacy block is removed in Phase C.)

- [ ] **Step 7: Commit**

```bash
git add species/ultra_relativistic.{h,cpp}
git commit -m "ur: migrate to per-species PerturbLayout incl. ufa switch hook"
```

---

## Task 16: Migrate `IDRSpecies`

**Files:**
- Modify: `species/idr.h` (or wherever IDR class is defined)

IDR has `pv->index_pt_delta_idr/theta_idr/shear_idr/l3_idr`, `pv->l_max_idr`, plus the `idr_nature` (free-streaming vs fluid) gating which fields are used.

- [ ] **Step 1: Locate the IDR class**

```bash
grep -rn "class.*IDR.*BaseSpecies\|IDRSpecies\b" species/ | head -5
```

- [ ] **Step 2: Define `IDRPerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1;
  int idx_theta = -1;
  int idx_shear = -1;  // free-streaming only
  int idx_l3    = -1;  // free-streaming only
  int l_max     = -1;  // free-streaming only
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 3-5: Standard migration with switch hooks**

Port Register, reads, CopyPerturbationsAcrossSwitch (slot-by-slot for IDR's tca/rsa transitions).

- [ ] **Step 6: Commit**

```bash
git add species/idr.h <species/idr.cpp if exists>
git commit -m "idr: migrate to per-species PerturbLayout incl. tca/rsa switch hook"
```

---

## Task 17: Migrate `IDM_DRSpecies` and `IDM_DRMDSpecies` (interacting DM species)

**Files:**
- Modify: `species/idm_dr.h`
- Modify: `species/idm_drmd.h` (and `.cpp` files if any)

Both have the shape `pv->index_pt_delta_*/theta_*` for their interacting species.

- [ ] **Steps 1-5**: Standard pattern. Each gets its own `PerturbLayout`. CopyPerturbationsAcrossSwitch handles tca_idm_dr / tca_idm_drmd transitions (slot-by-slot).

```cpp
// species/idm_dr.h
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1;
  int idx_theta = -1;
};
// species/idm_drmd.h — same shape.
```

- [ ] **Step 6: Commit**

```bash
git add species/idm_dr.h species/idm_drmd.h <any .cpp files>
git commit -m "idm_dr, idm_drmd: migrate to per-species PerturbLayout"
```

---

## Task 18: Migrate `BaryonsSpecies`

**Files:**
- Modify: `species/baryons.h`
- Modify: `species/baryons.cpp`

Baryons have `pv->index_pt_delta_b`, `pv->index_pt_theta_b`. Tensor mode adds nothing (baryons participate in scalar only via tight coupling with photons).

- [ ] **Step 1: Define `BaryonsPerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int idx_delta = -1;
  int idx_theta = -1;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Steps 2-5: Standard migration**

Port `pv->index_pt_delta_b/theta_b` to `layout.idx_delta/theta`. Update `Baryons::PerturbDerivs`, `Delta`, `Theta`, `ApplyInitialConditions`, `FillSources` to take and use the layout.

`SpeciesCollection::baryons()` still returns the species reference (not a layout), so module-side hot-path access stays via `all_species_.baryons().Delta(layout, ...)`. The dispatch site needs to look up the layout via `all_species_.index_of("Baryons")`.

- [ ] **Step 6: Commit**

```bash
git add species/baryons.h species/baryons.cpp source/perturbations_module.cpp
git commit -m "baryons: migrate to per-species PerturbLayout"
```

---

## Task 19: Migrate `PhotonsSpecies` (largest single species)

**Files:**
- Modify: `species/photons.h`
- Modify: `species/photons.cpp`
- Modify: `source/perturbations_module.cpp`

Photons have many fields:
- Scalar: `pv->index_pt_delta_g/theta_g/shear_g/l3_g`, `l_max_g`, `pv->index_pt_pol{0,1,2,3}_g`, `l_max_pol_g`
- Tensor: subset of the above

Plus the tight-coupling approximation has its own switch behaviour.

- [ ] **Step 1: Define `PhotonsPerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  // Scalar / tensor — the same struct holds whichever mode this pv is for.
  int idx_delta = -1, idx_theta = -1, idx_shear = -1, idx_l3 = -1;
  int idx_pol0  = -1, idx_pol1  = -1, idx_pol2  = -1, idx_pol3 = -1;
  int l_max     = -1;
  int l_max_pol = -1;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 2: Override `RegisterPerturbationIndices` (scalar) and `RegisterTensorPerturbationIndices`**

Port the existing bodies (in `species/photons.cpp` for scalar, also separately for tensor). Dual-write `layout.idx_*` AND `pv->index_pt_*_g`. Both Register methods write into the SAME `PhotonsPerturbLayout` instance (one pv has one mode, so no conflict).

- [ ] **Step 3: Override read virtuals**

Port `PerturbDerivs`, `PerturbTensorDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`. All use `layout.idx_*` instead of `pv->index_pt_*_g`.

The tight-coupling approximation logic in `PerturbDerivs` (lots of conditional branches based on `tca_on/off`) ports without structural change — just substitution.

- [ ] **Step 4: Override `CopyPerturbationsAcrossSwitch`**

The 6 non-NCDM module-side switch-copy blocks include photon copies for tca_off, rsa_on, etc. Move these into `Photons::CopyPerturbationsAcrossSwitch`:

```cpp
void PhotonsSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y, double* new_y,
    const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0) return;
  // delta, theta, shear, l3, pol0..pol3, all multipoles up to l_max.
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
  new_y[new_l.idx_shear] = old_y[old_l.idx_shear];
  new_y[new_l.idx_l3]    = old_y[old_l.idx_l3];
  for (int l = 4; l <= new_l.l_max; ++l)
    new_y[new_l.idx_delta + l] = old_y[old_l.idx_delta + l];
  if (new_l.idx_pol0 >= 0) {
    new_y[new_l.idx_pol0] = old_y[old_l.idx_pol0];
    new_y[new_l.idx_pol1] = old_y[old_l.idx_pol1];
    new_y[new_l.idx_pol2] = old_y[old_l.idx_pol2];
    new_y[new_l.idx_pol3] = old_y[old_l.idx_pol3];
    for (int l = 4; l <= new_l.l_max_pol; ++l)
      new_y[new_l.idx_pol0 + l] = old_y[old_l.idx_pol0 + l];
  }
}
```

- [ ] **Step 5: Module-side dispatch update**

Find the photon-specific calls in `source/perturbations_module.cpp`:
```bash
grep -n "all_species_.photons()" source/perturbations_module.cpp
```

Each `photons().PerturbDerivs(tau, y, dy, *pppaw)` becomes `photons().PerturbDerivs(*ppw->pv->species_layouts[photons_idx], tau, y, dy, *pppaw)`. Cache `photons_idx = all_species_.index_of("Photons")` once at the top of the function.

- [ ] **Step 6: Build + regression**

Bit-identical.

- [ ] **Step 7: Commit**

```bash
git add species/photons.h species/photons.cpp source/perturbations_module.cpp
git commit -m "photons: migrate to per-species PerturbLayout (scalar+tensor) incl. switch hook"
```

---

## Task 20: Migrate `NCDMBaseSpecies` family (NCDM, DNCDM, NCDMInteracting)

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_species.{h,cpp}`
- Modify: `species/dncdm_species.{h,cpp}`
- Modify: `species/ncdm_interacting_species.{h,cpp}`
- Modify: `source/perturbations_module.cpp`

This is the big one. NCDM family currently uses `pv->index_ncdm_[ncdm_id_][q]`, `pv->l_max_ncdm[ncdm_id_]`, `pv->q_size_ncdm[ncdm_id_]`. Migrate to per-species layout.

- [ ] **Step 1: Define `NCDMPerturbLayout` on `NCDMBaseSpecies`**

Edit `species/ncdm_base_species.h`. Add inside the class:

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  int l_max  = -1;
  int q_size = -1;
  std::vector<int> index_per_q;  // absolute offsets into pv->y, one per momentum bin

  int total_size() const { return q_size * (l_max + 1); }
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

(Both NCDMSpecies and DNCDMSpecies inherit, but they need their own `CreatePerturbLayout` — actually no, since the layout struct is the same shape for all NCDM-family species, the base's CreatePerturbLayout suffices. NCDMInteractingSpecies inherits via NCDMSpecies's RegisterPerturbationIndices.)

- [ ] **Step 2: Override layout-based `RegisterPerturbationIndices` on NCDMSpecies and DNCDMSpecies**

In `species/ncdm_species.cpp`:
```cpp
void NCDMSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                perturb_vector* pv,
                                                const precision* ppr,
                                                int& index_pt,
                                                const perturb_workspace* ppw,
                                                int /*gauge*/) {
  if (!pba_->has_ncdm) return;
  auto& layout = static_cast<NCDMBaseSpecies::PerturbLayout&>(base);

  if (ncdm_id_ == 0) {
    pv->index_pt_psi0_ncdm1 = index_pt;  // dual-write (deleted in Phase E)
  }
  index_pt_psi0_ = index_pt;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);
  if (fa_on) {
    layout.l_max  = 2;
    layout.q_size = 1;
  } else {
    layout.l_max  = ppr->l_max_ncdm;
    layout.q_size = q_size();
  }
  layout.index_per_q.clear();
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  // Dual-write to legacy pv arrays (deleted in Phase D).
  pv->l_max_ncdm[ncdm_id_]  = layout.l_max;
  pv->q_size_ncdm[ncdm_id_] = layout.q_size;
  pv->index_ncdm_[ncdm_id_] = layout.index_per_q;

  index_pt += layout.total_size();
}
```

In `species/dncdm_species.cpp`: same structure, but no `ncdmfa` branch (DNCDM doesn't participate in fluid approximation), and no `if (ncdm_id_ == 0)` block.

- [ ] **Step 3: Override `RegisterTensorPerturbationIndices` on NCDMBaseSpecies**

Tensor mode for NCDM uses full hierarchy (no fluid approximation). Implementation lives on the BASE so both NCDM and DNCDM inherit:

```cpp
void NCDMBaseSpecies::RegisterTensorPerturbationIndices(
    BaseSpecies::PerturbLayout& base, perturb_vector* pv,
    int& index_pt, const perturb_workspace* /*ppw*/, int /*gauge*/) {
  auto& layout = static_cast<PerturbLayout&>(base);
  layout.l_max  = pt_l_max_ncdm_from_ppr_;  // captured during scalar Register
  layout.q_size = q_size();
  layout.index_per_q.clear();
  for (int iq = 0; iq < layout.q_size; ++iq)
    layout.index_per_q.push_back(index_pt + iq * (layout.l_max + 1));

  pv->l_max_ncdm[ncdm_id_]  = layout.l_max;
  pv->q_size_ncdm[ncdm_id_] = layout.q_size;
  pv->index_ncdm_[ncdm_id_] = layout.index_per_q;

  index_pt += layout.total_size();
}
```

(`pt_l_max_ncdm_from_ppr_` is captured during the scalar Register, when ppr is in scope.)

- [ ] **Step 4: Override read virtuals on NCDMSpecies and DNCDMSpecies**

For each of `PerturbDerivs`, `PerturbTensorDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear`, `ApplyInitialConditions`, `FillSources`: add layout-based override; reduce legacy to no-op/0.

Inside each method, downcast and use `layout.q_size`, `layout.l_max`, `layout.index_per_q[iq]` instead of `pv->q_size_ncdm[ncdm_id_]`, `pv->l_max_ncdm[ncdm_id_]`, `pv->index_ncdm_[ncdm_id_][iq]`.

`NCDMInteractingSpecies` inherits from NCDMSpecies, so it doesn't need its own overrides — its layout is the inherited NCDMBaseSpecies::PerturbLayout populated by the base's RegisterPerturbationIndices.

- [ ] **Step 5: Override `CopyPerturbationsAcrossSwitch` on NCDMSpecies (FA-collapse case)**

This is the non-trivial switch case. From the spec, port the FA-collapse logic from `perturbations_module.cpp:4344-4383`:

```cpp
void NCDMSpecies::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y, double* new_y,
    const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const NCDMBaseSpecies::PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const NCDMBaseSpecies::PerturbLayout&>(new_base);
  if (old_l.q_size < 0 || new_l.q_size < 0) return;

  const bool fa_was_on = (old_l.q_size == 1 && old_l.l_max == 2);
  const bool fa_is_on  = (new_l.q_size == 1 && new_l.l_max == 2);

  if (fa_was_on == fa_is_on) {
    // Slot-by-slot copy.
    for (int iq = 0; iq < new_l.q_size; ++iq) {
      const int new_base = new_l.index_per_q[iq];
      const int old_base = old_l.index_per_q[iq];
      for (int l = 0; l <= new_l.l_max; ++l)
        new_y[new_base + l] = old_y[old_base + l];
    }
    return;
  }

  if (!fa_was_on && fa_is_on) {
    // Collapse old hierarchy → (delta, theta, shear).
    const double a            = ctx.a;
    const double k            = ctx.k;
    const double a_today      = ctx.a_today;
    const double rho          = Rho(ctx.pvecback);
    const double p            = P(ctx.pvecback);
    const double rho_plus_p   = rho + p;
    const double M_local      = M();
    const double factor_local = factor() * pow(a_today / a, 4);

    const int idx_new = new_l.index_per_q[0];
    for (int l = 0; l <= 2; ++l) new_y[idx_new + l] = 0.0;

    double delta = 0., theta = 0., shear = 0.;
    for (int iq = 0; iq < old_l.q_size; ++iq) {
      const int idx_old    = old_l.index_per_q[iq];
      const double w0      = w()[iq];
      const double q       = this->q()[iq];
      const double epsilon = sqrt(q * q + a * a * M_local * M_local);

      delta += w0 * q * q * epsilon * old_y[idx_old];
      theta += w0 * q * q * q * old_y[idx_old + 1];
      shear += w0 * q * q * q * q / epsilon * old_y[idx_old + 2];
    }
    delta *= factor_local / rho;
    theta *= k * factor_local / rho_plus_p;
    shear *= 2. / 3. * factor_local / rho_plus_p;

    new_y[idx_new]     = delta;
    new_y[idx_new + 1] = theta;
    new_y[idx_new + 2] = shear;
    return;
  }

  // fa_on → fa_off: unreachable in CLASS (once fa is on, stays on).
  // Defensive fallback: leave new_y as zero-initialized — will warn during
  // integrator step-size if ever hit. NOT an assert (assert(false) crashes
  // multi-threaded scenarios as we discovered in the previous attempt).
}
```

- [ ] **Step 6: Override `MarkUsedInSources` on `NCDMBaseSpecies`**

```cpp
void NCDMBaseSpecies::MarkUsedInSources(
    const BaseSpecies::PerturbLayout& base, int* used_in_sources) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  if (layout.q_size < 0) return;
  // NCDM family: multipoles l > 2 do NOT enter source functions.
  for (int iq = 0; iq < layout.q_size; ++iq) {
    const int b = layout.index_per_q[iq];
    for (int l = 0; l <= layout.l_max; ++l)
      if (l > 2) used_in_sources[b + l] = _FALSE_;
  }
}
```

- [ ] **Step 7: Module-side dispatch update for NCDM dispatch**

In `source/perturbations_module.cpp`, every site iterating `ncdm_species_sorted_` for evolution / sources / tensor changes to iterate `all_species_` filtered by `dynamic_cast<NCDMBaseSpecies*>` and using the species's index in `species_layouts`.

- [ ] **Step 8: Build + regression**

Bit-identical.

- [ ] **Step 9: Commit**

```bash
git add species/ncdm_base_species.{h,cpp} species/ncdm_species.{h,cpp} \
        species/dncdm_species.{h,cpp} species/ncdm_interacting_species.{h,cpp} \
        source/perturbations_module.cpp
git commit -m "ncdm-family: migrate to per-species PerturbLayout incl. fa-collapse hook"
```

---

## Task 21: Migrate `DCDM_DR_Species` composite (nested composition)

**Files:**
- Modify: `species/dcdm_dr_species.h`
- Modify: `species/dcdm_dr_species.cpp`

DCDM_DR is a composite owning DCDMSpecies + DRSpecies children. Children's layouts live INSIDE the composite's layout.

- [ ] **Step 1: Define `DCDM_DR_PerturbLayout`**

```cpp
// species/dcdm_dr_species.h
struct PerturbLayout : BaseSpecies::PerturbLayout {
  DCDMSpecies::PerturbLayout dcdm;
  DRSpecies::PerturbLayout   dr;
};

std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
  return std::make_unique<PerturbLayout>();
}
```

- [ ] **Step 2: Override Register, delegating to children**

```cpp
void DCDM_DR_Species::RegisterPerturbationIndices(
    BaseSpecies::PerturbLayout& base, perturb_vector* pv,
    const precision* ppr, int& index_pt,
    const perturb_workspace* ppw, int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  dcdm_->RegisterPerturbationIndices(my.dcdm, pv, ppr, index_pt, ppw, gauge);
  dr_->RegisterPerturbationIndices(my.dr,   pv, ppr, index_pt, ppw, gauge);
}
```

(`dcdm_` and `dr_` are existing raw-pointer children members.)

- [ ] **Step 3: Override read virtuals, delegating**

```cpp
void DCDM_DR_Species::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                    double tau, const double* y, double* dy,
                                    const perturb_parameters_and_workspace& ppaw) {
  const auto& my = static_cast<const PerturbLayout&>(base);
  dcdm_->PerturbDerivs(my.dcdm, tau, y, dy, ppaw);
  dr_->PerturbDerivs(my.dr,    tau, y, dy, ppaw);
}
```

Same shape for `Delta` (sum?), `Theta`, etc. — depends on existing CompositeSpecies semantics. Inspect the existing `CompositeSpecies` base to see how it aggregates.

- [ ] **Step 4: Build + regression**

Bit-identical.

- [ ] **Step 5: Commit**

```bash
git add species/dcdm_dr_species.{h,cpp}
git commit -m "dcdm_dr: migrate composite to nested PerturbLayout composition"
```

---

## Task 22: Migrate `DNCDM_DR_Species` composite

**Files:**
- Modify: `species/dncdm_dr_species.{h,cpp}`
- Modify: `species/dncdm_decay_radiation_species.{h,cpp}`

Same pattern as Task 21. The composite owns `DNCDMSpecies` + `DNCDM_DecayRadiationSpecies` children.

- [ ] **Step 1: Define `DNCDM_DR_PerturbLayout`**

```cpp
struct PerturbLayout : BaseSpecies::PerturbLayout {
  NCDMBaseSpecies::PerturbLayout dncdm;  // DNCDM uses NCDMBaseSpecies::PerturbLayout
  DRSpecies::PerturbLayout       dr;
};
```

- [ ] **Steps 2-5: Same as Task 21**

```bash
git add species/dncdm_dr_species.{h,cpp} species/dncdm_decay_radiation_species.{h,cpp}
git commit -m "dncdm_dr: migrate composite to nested PerturbLayout composition"
```

---

## Task 23: Migrate `IDM_DR_IDR_Species` and `IDM_DRMD_IDR_DRMD_Species` composites

**Files:**
- Modify: `species/idm_dr_idr_species.{h,cpp}`
- Modify: `species/idm_drmd_idr_drmd_species.{h,cpp}`

Same composite pattern.

- [ ] **Steps 1-5**: Define each composite's nested layout, override Register and read methods, delegating to children's layouts.

```bash
git add species/idm_dr_idr_species.{h,cpp} species/idm_drmd_idr_drmd_species.{h,cpp}
git commit -m "idm_dr_idr, idm_drmd_idr_drmd: migrate composites to nested PerturbLayout"
```

---

## Phase C — Bookkeeping hooks

## Task 24: Wire `CopyPerturbationsAcrossSwitch` into `perturb_vector_init`

**Files:**
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 1: Insert dispatch loop before each existing module-side switch-copy block**

In `perturb_vector_init`, the `if (pa_old != NULL)` branch contains 8 NCDM-shaped switch-copy blocks plus analogous UR/IDR/IDM_DRMD/Photons blocks. **Don't delete them yet** — first add the species dispatch loop. The species hooks dual-write the same slots, so behaviour is unchanged (each slot is written twice with the same value).

Find the start of the `else` branch (the "switching approximation" branch) where copy logic begins. Add at the top:

```cpp
// Dispatch per-species CopyPerturbationsAcrossSwitch. Each migrated species
// performs its own slot-by-slot copy or non-trivial migration (NCDM fa-collapse).
{
  PerturbSwitchContext switch_ctx;
  switch_ctx.k        = k;
  switch_ctx.a        = ppw->pvecback[background_module_->index_bg_a_];
  switch_ctx.a_today  = pba->a_today;
  switch_ctx.pvecback = ppw->pvecback;
  for (size_t i = 0; i < all_species_.size(); ++i) {
    all_species_[i]->CopyPerturbationsAcrossSwitch(
        *ppw->pv->species_layouts[i],
        *ppv->species_layouts[i],
        ppw->pv->y, ppv->y,
        switch_ctx);
  }
}
```

- [ ] **Step 2: Wire `MarkUsedInSources`**

After the `for (index_pt = 0; index_pt < ppv->pt_size; index_pt++) ppv->used_in_sources[index_pt] = _TRUE_;` default-init loop, add:

```cpp
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->MarkUsedInSources(*ppv->species_layouts[i], ppv->used_in_sources);
}
```

- [ ] **Step 3: Build + regression**

Bit-identical (every slot is written by both the legacy block and the new species hook; same value).

- [ ] **Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: wire CopyPerturbationsAcrossSwitch + MarkUsedInSources hooks"
```

---

## Task 25: Delete the 8 NCDM switch-copy blocks + analogous UR/IDR/IDM_DRMD/Photons blocks

**Files:**
- Modify: `source/perturbations_module.cpp`

These are the legacy module-side blocks now redundant after Task 24's species hooks.

- [ ] **Step 1: Delete the 8 NCDM-shaped blocks**

Locate each block (search for `index_ncdm_[n]` inside `if (!ncdm_species_sorted_.empty())` guards in the switch branch). There are 8 — at the various transition cases (tca-off, rsa-on, ufa-off, rsa_idr-off, tca_idr-off, NCDMFA-on FA-collapse, tca_idm_drmd-off DRMD branch, tensor mode).

Delete each entire block (the `if (!ncdm_species_sorted_.empty())` guard plus the inner loop).

- [ ] **Step 2: Delete the analogous UR/IDR/IDM_DRMD/Photons blocks**

Each species's `CopyPerturbationsAcrossSwitch` override (Tasks 15, 16, 17, 19) handles these. Find and delete the slot-by-slot copy blocks for UR (`pv->index_pt_delta_ur` etc.), IDR (`pv->index_pt_delta_idr`), IDM_DRMD (`pv->index_pt_delta_idm_drmd`), Photons (`pv->index_pt_delta_g`). Also delete the NCDM `used_in_sources` masking block at lines 3571-3589 (now done by `MarkUsedInSources` hook).

- [ ] **Step 3: Build + regression**

Still bit-identical (the species hooks are doing the work).

- [ ] **Step 4: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: delete legacy switch-copy blocks (replaced by species hooks)"
```

---

## Task 26: Delete tensor-mode NCDM emission blocks (replaced by species hooks)

**Files:**
- Modify: `source/perturbations_module.cpp`

The tensor PerturbDerivs block (lines ~7290), tensor GW source block (lines ~6238), tensor source title emission (lines ~2824), and tensor registration block (lines ~3456) were all migrated to species hooks. Delete the module-side bodies.

- [ ] **Step 1: Replace each block with a per-species dispatch loop**

Pattern for `perturb_derivs_member` tensor block:
```cpp
if (evolve_tensor_ncdm_ == _TRUE_) {
  for (size_t i = 0; i < all_species_.size(); ++i)
    all_species_[i]->PerturbTensorDerivs(*ppw->pv->species_layouts[i], tau, y, dy, *pppaw);
}
```

For tensor GW source (line ~6238):
```cpp
if (evolve_tensor_ncdm_ == _TRUE_) {
  for (size_t i = 0; i < all_species_.size(); ++i)
    all_species_[i]->ContributeTensorGwSource(a, pba->a_today, y, ppw, *ppw->pv->species_layouts[i]);
}
```

For tensor source title emission (line ~2824):
```cpp
if (evolve_tensor_ncdm_ == _TRUE_) {
  for (size_t i = 0; i < all_species_.size(); ++i)
    all_species_[i]->WriteTensorOutputColumnTitles(tensor_titles_);
}
```

For tensor registration (line ~3456):
```cpp
if (evolve_tensor_ncdm_ == _TRUE_) {
  ppv->index_pt_psi0_ncdm1 = index_pt;  // dual-write (Phase E will delete)
  for (size_t i = 0; i < all_species_.size(); ++i) {
    auto& layout = *ppv->species_layouts[i];
    all_species_[i]->RegisterTensorPerturbationIndices(layout, ppv, index_pt, ppw, ppt->gauge);
  }
}
```

- [ ] **Step 2: Build + regression**

Bit-identical.

- [ ] **Step 3: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: tensor NCDM emission via species hooks (delete module blocks)"
```

---

## Phase D — Cleanup of legacy fields

## Task 27: Delete migrated bare fields from `perturb_vector` (group 1: simple species)

**Files:**
- Modify: `source/perturbations.h`
- Modify: `species/*.cpp` (delete dual-writes)
- Modify: `source/perturbations_module.cpp` (delete remaining stale references, if any)

Delete from `struct perturb_vector`:
- `index_pt_delta_g/theta_g/shear_g/l3_g`, `l_max_g`, `index_pt_pol{0,1,2,3}_g`, `l_max_pol_g` (Photons)
- `index_pt_delta_b/theta_b` (Baryons)
- `index_pt_delta_cdm/theta_cdm` (CDM)
- `index_pt_delta_dcdm/theta_dcdm` (DCDM)
- `index_pt_delta_fld/theta_fld/Gamma_fld` (Fluid)
- `index_pt_phi_scf/phi_prime_scf` (ScalarField)

For each field, also delete the dual-write in the species's `RegisterPerturbationIndices` body. Search-and-delete:
```bash
grep -rn "pv->index_pt_delta_g\b\|pv->index_pt_theta_g\b" source/ species/
# ... ensure all reads are gone (only writes remain — those are the dual-writes to delete).
```

- [ ] **Step 1: Verify no remaining reads**

For each field to delete, verify only dual-writes remain (in the species's Register method):
```bash
for field in index_pt_delta_g index_pt_theta_g index_pt_shear_g index_pt_l3_g \
             l_max_g index_pt_pol0_g index_pt_pol1_g index_pt_pol2_g index_pt_pol3_g \
             l_max_pol_g index_pt_delta_b index_pt_theta_b index_pt_delta_cdm \
             index_pt_theta_cdm index_pt_delta_dcdm index_pt_theta_dcdm \
             index_pt_delta_fld index_pt_theta_fld index_pt_Gamma_fld \
             index_pt_phi_scf index_pt_phi_prime_scf; do
  echo "=== $field ==="
  grep -rn "$field\b" source/ species/
done
```

Expected: each field has at most ONE remaining hit per species (the dual-write line in Register*PerturbationIndices). If reads remain, fix them before deleting.

- [ ] **Step 2: Delete dual-writes**

For each field, in the species's Register method, remove the `pv->index_pt_X = index_pt;` line.

- [ ] **Step 3: Delete fields from `struct perturb_vector`**

Edit `source/perturbations.h`. Remove the field declarations.

- [ ] **Step 4: Build + regression**

Bit-identical.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations.h source/perturbations_module.cpp species/*.cpp species/*.h
git commit -m "perturb_vector: delete migrated bare fields (Photons, Baryons, CDM, DCDM, Fluid, ScalarField)"
```

---

## Task 28: Delete migrated bare fields (group 2: UR, IDR, IDM_DR, IDM_DRMD)

**Files:**
- Modify: `source/perturbations.h`
- Modify: `species/*.cpp` (delete dual-writes)

Delete:
- `index_pt_delta_ur/theta_ur/shear_ur/l3_ur`, `l_max_ur`
- `index_pt_delta_idr/theta_idr/shear_idr/l3_idr`, `l_max_idr`
- `index_pt_delta_idr_drmd/theta_idr_drmd`
- `index_pt_delta_idm_dr/theta_idm_dr`
- `index_pt_delta_idm_drmd/theta_idm_drmd`

- [ ] **Steps 1-5: Same pattern as Task 27**

```bash
git commit -m "perturb_vector: delete migrated bare fields (UR, IDR, IDM_DR, IDM_DRMD)"
```

---

## Task 29: Delete migrated bare fields (group 3: NCDM-family flat arrays)

**Files:**
- Modify: `source/perturbations.h`
- Modify: `species/*.cpp` (delete dual-writes)

Delete:
- `index_pt_psi0_ncdm1`
- `N_ncdm`
- `l_max_ncdm` (raw pointer)
- `q_size_ncdm` (raw pointer)
- `index_ncdm_` (`std::map<int, std::vector<int>>`)
- `l_max_ncdm_storage`
- `q_size_ncdm_storage`

Plus delete the manual storage allocation at `perturbations_module.cpp:3357-3360`, `3452-3455`, etc. and all dual-writes in NCDMSpecies/DNCDMSpecies.

- [ ] **Steps 1-5: Standard cleanup**

```bash
git commit -m "perturb_vector: delete migrated NCDM flat arrays"
```

---

## Phase E — `ncdm_id` removal

## Task 30: Delete `ncdm_id_` field, `SetNcdmId`, and `ncdm_id()` accessor

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_species.{h,cpp}`
- Modify: `species/dncdm_species.{h,cpp}`
- Modify: `species/dncdm_dr_species.{h,cpp}`
- Modify: `species/species_build_context.h`

- [ ] **Step 1: Verify no remaining reads of `ncdm_id_` outside the species classes**

```bash
grep -rn "ncdm_id_\b\|->ncdm_id()\|->SetNcdmId\|ncdm_id_next" source/ species/ include/
```

Each remaining hit must be one of:
- The deletion targets (the field declaration, the SetNcdmId method, the ncdm_id() accessor).
- The `ctx.ncdm_id_next` counter in SpeciesBuildContext.

If anything else remains (e.g., a module-side dispatch using ncdm_id), fix it first.

- [ ] **Step 2: Delete the field, methods, and counter**

In `species/ncdm_base_species.h`, remove:
- `int ncdm_id_ = -1;` field
- `virtual void SetNcdmId(int id) = 0;` declaration
- `int ncdm_id() const { return ncdm_id_; }` accessor

In `species/ncdm_species.h`, `species/dncdm_species.h`, `species/dncdm_dr_species.h`: remove the `SetNcdmId` overrides.

In `species/ncdm_species.cpp`: in the `CreateAll` factory, remove `sp->SetNcdmId((*ctx.ncdm_id_next)++);`.

In `species/species_build_context.h`: remove `int* ncdm_id_next;` field.

In `source/cosmology.cpp` (or wherever SpeciesBuildContext is constructed): remove the `ncdm_id_next` initialization.

- [ ] **Step 3: Build + regression**

Bit-identical (the int was only used by deleted code paths).

- [ ] **Step 4: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_species.{h,cpp} species/dncdm_species.{h,cpp} \
        species/dncdm_dr_species.{h,cpp} species/species_build_context.h source/cosmology.cpp
git commit -m "ncdm: delete ncdm_id_ field, SetNcdmId, and ncdm_id_next counter"
```

---

## Task 31: Delete `pba->N_ncdm`

**Files:**
- Modify: `source/background.h`
- Modify: 5 consumer sites (find via grep)

- [ ] **Step 1: Find all consumers**

```bash
grep -rn "pba->N_ncdm\b\|pba_->N_ncdm\b" source/ species/
```

- [ ] **Step 2: Replace with dynamic_cast filter sums**

For each consumer, replace `pba->N_ncdm` with:
```cpp
[&]{ size_t n = 0; for (size_t i = 0; i < all_species_.size(); ++i) {
        if (dynamic_cast<NCDMBaseSpecies*>(all_species_[i])) ++n;
     } return n; }()
```

Or extract a helper function `CountNCDMSpecies(const SpeciesCollection&)` to avoid repeating. Mark each substitution with `// TODO(architecture): unify with Omega0_b/Omega0_cdm tally`.

- [ ] **Step 3: Delete `N_ncdm` from `struct background`**

Edit `source/background.h`. Remove `int N_ncdm;`.

- [ ] **Step 4: Build + regression**

Bit-identical.

- [ ] **Step 5: Commit**

```bash
git add source/background.h source/perturbations_module.cpp source/cosmology.cpp species/*.cpp
git commit -m "background: delete pba->N_ncdm (replace 5 consumers with NCDM filter)"
```

---

## Task 32: Delete `pba->Omega0_ncdm_tot` and `pba->N_decay_dr`

**Files:**
- Modify: `source/background.h`
- Modify: consumer sites

- [ ] **Step 1: Find all consumers**

```bash
grep -rn "pba->Omega0_ncdm_tot\|pba->N_decay_dr" source/ species/
```

- [ ] **Step 2: Replace with dynamic_cast filter sums**

Pattern same as Task 31. `Omega0_ncdm_tot` becomes a sum of `dynamic_cast<NCDMBaseSpecies*>(sp)->GetOmega0()`. `N_decay_dr` becomes a count of DNCDM_DR composites in `all_species_`. Mark each `// TODO(architecture)`.

- [ ] **Step 3: Delete from background struct**

Remove `Omega0_ncdm_tot` and `N_decay_dr` declarations.

- [ ] **Step 4: Build + regression**

Bit-identical.

- [ ] **Step 5: Commit**

```bash
git add source/background.h source/cosmology.cpp source/perturbations_module.cpp species/*.cpp
git commit -m "background: delete pba->Omega0_ncdm_tot, pba->N_decay_dr"
```

---

## Task 33: Delete `PerturbationsModule::ncdm_species_sorted_`

**Files:**
- Modify: `source/perturbations_module.h`
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 1: Verify no remaining uses**

```bash
grep -rn "ncdm_species_sorted_\b" source/ species/
```

Each remaining hit must be the declaration, construction, or destruction. Replace any remaining iteration:
- `if (!ncdm_species_sorted_.empty())` → `if (std::any_of(all_species_.begin(), all_species_.end(), [](const auto& e){ return dynamic_cast<NCDMBaseSpecies*>(e.species.get()) != nullptr; }))`
- `for (auto* sp : ncdm_species_sorted_)` → index loop with dynamic_cast filter.

- [ ] **Step 2: Delete the field and its construction**

Edit `source/perturbations_module.h`: remove `std::vector<NCDMBaseSpecies*> ncdm_species_sorted_;` (or whatever its declared type is).

Edit `source/perturbations_module.cpp`: find where `ncdm_species_sorted_` is populated (likely in `perturbations_init` or constructor) and delete.

- [ ] **Step 3: Build + regression**

Bit-identical.

- [ ] **Step 4: Commit**

```bash
git add source/perturbations_module.h source/perturbations_module.cpp
git commit -m "perturbations: delete ncdm_species_sorted_ (replaced by all_species_ + dynamic_cast)"
```

---

## Task 34: Refactor `RescaledNCDMPerturbations(int n_ncdm, ...)` to take a species pointer

**Files:**
- Modify: `source/perturbations_module.h`
- Modify: `source/perturbations_module.cpp`
- Modify: any callers

- [ ] **Step 1: Find current signature and consumers**

```bash
grep -rn "RescaledNCDMPerturbations" source/ species/
```

Currently:
```cpp
std::tuple<double, double, double> RescaledNCDMPerturbations(int n_ncdm, double a, double k, perturb_workspace* ppw);
```

- [ ] **Step 2: Rewrite to take species pointer**

```cpp
// In source/perturbations_module.h:
std::tuple<double, double, double> RescaledNCDMPerturbations(
    BaseSpecies* species, double a, double k, perturb_workspace* ppw);

// In source/perturbations_module.cpp:
std::tuple<double, double, double> PerturbationsModule::RescaledNCDMPerturbations(
    BaseSpecies* species, double a, double k, perturb_workspace* ppw) {
  if (auto* composite = dynamic_cast<DNCDM_DR_Species*>(species)) {
    return composite->dncdm().RescaledPerturbations(a, k, ppw);
  }
  throw std::runtime_error("RescaledNCDMPerturbations: invalid species");
}
```

Update each caller to pass the species pointer they already have in scope.

- [ ] **Step 3: Build + regression**

Bit-identical.

- [ ] **Step 4: Commit**

```bash
git add source/perturbations_module.h source/perturbations_module.cpp
git commit -m "perturbations: RescaledNCDMPerturbations(int) -> (BaseSpecies*)"
```

---

## Phase F — Output column titles

## Task 35: Rewrite NCDM-family output column titles to instance-name based

**Files:**
- Modify: `species/ncdm_species.cpp`
- Modify: `species/dncdm_species.cpp`
- Modify: `source/perturbations_module.cpp` (for tensor-title emission)

- [ ] **Step 1: Find existing title emissions**

```bash
grep -rn 'snprintf.*ncdm\[\|"ncdm\[%d\]\|class_store_columntitle.*ncdm\[' species/ source/
```

These are the titles using `[%d]` formatting with `ncdm_id_` or `n` (the legacy id). Examples:
- `(.)number_ncdm[%d]`
- `(.)rho_ncdm[%d]`
- `(.)p_ncdm[%d]`
- `(.)pseudo_p_ncdm[%d]`
- `(.)dlnfdlnq_ncdm[%d][q]`
- `delta_ncdm[%d]`, `theta_ncdm[%d]`, `shear_ncdm[%d]` (tensor)
- `d_ncdm[%d]`, `t_ncdm[%d]` (transfer)

- [ ] **Step 2: Replace `ncdm[%d]` with `<instance_name>`**

For each emission site, replace the format-string code with the species's instance name.

For `(.)rho_ncdm[%d]` in `NCDMSpecies::WriteBackgroundColumnTitles`:
```cpp
// Before:
snprintf(tmp, 40, "(.)rho_ncdm[%d]", ncdm_id_);
// After:
snprintf(tmp, 40, "(.)rho_%s", instance_name_.c_str());
```

(`instance_name_` is the existing field on `NCDMBaseSpecies`. Verify with `grep -n "instance_name_" species/ncdm_base_species.h`.)

Apply this pattern to every `ncdm[%d]` title in `WriteBackgroundColumnTitles`, `WriteOutputColumns`, `PrintVariables`, and the tensor title emission in `WriteTensorOutputColumnTitles`.

For the camb_format aggregate column at `WriteOutputColumns` (`-T_ncdm/k2`): keep as-is — it's an aggregate not indexed by id.

- [ ] **Step 3: Build + regression**

**This task knowingly breaks bit-identicality** for the column-title comparison. The numerical values at each column are unchanged but the `# 1: tau` header line differs. Verify by:
- `head -50 output/scen_ncdm_single_cl.dat` shows the new title lines.
- `head -50 /tmp/ncdm_baseline/scen_ncdm_single_cl.dat` shows the old `[0]` style.
- The numerical data rows match (same column ordering for single NCDM).

For multi-NCDM scenarios: the column ORDER changes too because lex order (`nu_a < nu_b`) differs from construction order (`nu_b` first). Verify by comparing column titles between baseline and new — `nu_a` should now precede `nu_b` regardless of construction order.

The diff helper from `<regression-check>` will report large `max_rel_diff` for the multi-unsorted scenario. **Inspect manually**: find the column titles, confirm `nu_a` precedes `nu_b`, confirm numerical values are physical.

- [ ] **Step 4: Commit**

```bash
git add species/ncdm_species.cpp species/dncdm_species.cpp source/perturbations_module.cpp
git commit -m "ncdm: rewrite output column titles to use instance name (BREAKING)"
```

---

## Phase G — Composite renames

## Task 36: Rename `DNCDM_DR_<id>` composite key to `<dncdm_instance>`

**Files:**
- Modify: `species/dncdm_dr_species.cpp`
- Modify: any string-literal lookups via `all_species_.at("DNCDM_DR_...")`

- [ ] **Step 1: Find current key construction**

```bash
grep -n 'DNCDM_DR_\|"DNCDM_DR_"' species/dncdm_dr_species.cpp source/perturbations_module.cpp source/input_module.cpp source/cosmology.cpp species/*.h
```

- [ ] **Step 2: Update key construction**

Edit `species/dncdm_dr_species.cpp`. Find the constructor:
```cpp
// Before:
DNCDM_DR_Species::DNCDM_DR_Species(...)
    : CompositeSpecies("DNCDM_DR_" + std::to_string(dncdm_arg->ncdm_id()), ...) ...
// After:
DNCDM_DR_Species::DNCDM_DR_Species(...)
    : CompositeSpecies(dncdm_arg->instance_name(), ...) ...
```

(`instance_name()` accessor on NCDMBaseSpecies — add if missing.) The DR child's name becomes `<dncdm_instance_name>_DR`.

- [ ] **Step 3: Update string-literal lookups**

Find any `all_species_.at("DNCDM_DR_...")` or `count("DNCDM_DR_...")` literal lookups. After this rename they will fail. Replace with iterating `all_species_` and dynamic_cast filtering for `DNCDM_DR_Species*`.

- [ ] **Step 4: Build + regression**

Bit-identical (the rename only affects internal keys; output-column titles already use instance names from Task 35; numerical computations are unchanged).

- [ ] **Step 5: Commit**

```bash
git add species/dncdm_dr_species.{h,cpp} source/perturbations_module.cpp source/input_module.cpp species/ncdm_base_species.h
git commit -m "dncdm_dr: composite key uses dncdm instance name (was DNCDM_DR_<id>)"
```

---

## Phase H — Shooter

## Task 37: Add shooter hooks to `BaseSpecies` and define `ShootingContext`

**Files:**
- Create: `species/shooting_context.h`
- Modify: `species/base_species.h`

- [ ] **Step 1: Define `ShootingContext`**

```cpp
// species/shooting_context.h
#pragma once
#include "background.h"

class FileContent;

struct ShootingContext {
  FileContent* pfc                = nullptr;
  const background* pba           = nullptr;
  // Currently-fitted parameter values during root finding:
  const double* current_values    = nullptr;
  // (Indices of this species's slots are stored on the species itself,
  //  populated during RegisterShootingIndices.)
};
```

- [ ] **Step 2: Add hooks to `BaseSpecies`**

Edit `species/base_species.h`. Add:

```cpp
  // ── Shooter / root-finding ────────────────────────────────────────────────
  // Single-threaded; runs inside InputModule before evolution. Each species
  // claims as many slots as it needs in the parameter vector.
  virtual void RegisterShootingIndices(int& /*index_sh*/,
                                        FileContent* /*pfc*/,
                                        const background& /*pba*/) {}
  virtual void ComputeShootingResidual(const ShootingContext& /*ctx*/,
                                        double* /*residual*/) const {}
  virtual void ComputeShootingGuess(const ShootingContext& /*ctx*/,
                                     double* /*guess*/) const {}
```

Include `shooting_context.h` (or forward-declare).

- [ ] **Step 3: Build + regression**

Bit-identical (hooks are no-op defaults).

- [ ] **Step 4: Commit**

```bash
git add species/shooting_context.h species/base_species.h
git commit -m "species: add shooter hooks (RegisterShootingIndices etc.) and ShootingContext"
```

---

## Task 38: Move `DNCDM_DR_Species` shooting cases into the species

**Files:**
- Modify: `species/dncdm_dr_species.{h,cpp}`
- Modify: `source/input_module.cpp` (delete the moved code)

- [ ] **Step 1: Find existing DNCDM_DR shooting code in InputModule**

```bash
grep -n "DNCDM_DR\|dncdm_id\|target_values.*dncdm" source/input_module.cpp | head -20
```

The block is around lines 3625-3845 of `source/input_module.cpp` (Omega_dncdm_decay_dr fitting).

- [ ] **Step 2: Add overrides on DNCDM_DR_Species**

```cpp
// species/dncdm_dr_species.h
void RegisterShootingIndices(int& index_sh, FileContent* pfc,
                              const background& pba) override;
void ComputeShootingResidual(const ShootingContext& ctx,
                              double* residual) const override;
void ComputeShootingGuess(const ShootingContext& ctx,
                           double* guess) const override;

// species/dncdm_dr_species.cpp
// Port the body of the DNCDM_DR shooting cases verbatim, replacing
// `target_values[idx + dncdm_id]` with `ctx.current_values[shooter_slot_]`.
```

- [ ] **Step 3: Delete the moved code from InputModule**

Remove the DNCDM_DR-specific blocks from `source/input_module.cpp`. The InputModule shooter dispatch will iterate `all_species_` calling `ComputeShootingResidual` / `ComputeShootingGuess` instead — that wiring is in Task 40.

- [ ] **Step 4: Build + regression**

Bit-identical (the shooter-side dispatch hasn't changed yet — the DNCDM_DR overrides are unused until Task 40).

- [ ] **Step 5: Commit**

```bash
git add species/dncdm_dr_species.{h,cpp} source/input_module.cpp
git commit -m "dncdm_dr: move shooting cases into species (RegisterShootingIndices + ...)"
```

---

## Task 39: Move `DCDM_DR_Species` and `ScalarFieldSpecies` shooting cases

**Files:**
- Modify: `species/dcdm_dr_species.{h,cpp}`
- Modify: `species/scalar_field.{h,cpp}`
- Modify: `source/input_module.cpp`

Same pattern as Task 38, but for DCDM_DR (Omega_dcdmdr fitting) and ScalarField (phi_ini, phi_prime_ini).

- [ ] **Steps 1-4**: Add overrides; delete moved code; build + regression.

```bash
git commit -m "dcdm_dr, scalar_field: move shooting cases into species"
```

---

## Task 40: Refactor `InputModule` shooter dispatch

**Files:**
- Modify: `source/input_module.cpp`

The remaining `theta_s` cosmological-level case stays in InputModule. All species-owned cases dispatch via `RegisterShootingIndices` / `ComputeShootingResidual` / `ComputeShootingGuess`.

- [ ] **Step 1: Replace shooter index allocation with species dispatch**

Find the shooter index-allocation code (where `target_values_size`, `target_values`, etc. are sized). Replace species-specific allocation with a per-species loop:

```cpp
int index_sh = 0;
// Cosmological-level theta_s slot first:
const int idx_theta_s = index_sh++;
// Per-species shooter slots:
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->RegisterShootingIndices(index_sh, pfc, pba);
}
const int total_shooter_slots = index_sh;
```

- [ ] **Step 2: Replace residual/guess computation with species dispatch**

For each iteration of the root-finder:
```cpp
ShootingContext ctx{ pfc, &pba, current_values };
// theta_s residual computed by InputModule (cosmological).
ComputeThetaSResidual(...);  // existing module-level
// Per-species residuals / guesses:
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->ComputeShootingResidual(ctx, residual);
}
```

(Same for ComputeShootingGuess.)

- [ ] **Step 3: Build + regression**

Bit-identical (same root-finder behavior, different code paths to compute it).

- [ ] **Step 4: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: shooter dispatch via species hooks (theta_s remains module-level)"
```

---

## Phase I — Final audit + PR

## Task 41: Final grep audit, full pytest, open PR

**Files:** none (verification + PR)

- [ ] **Step 1: Grep audit**

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/remove-ncdm-id

# ncdm_id should be completely gone:
grep -rn "ncdm_id" source/ species/ include/ | grep -v "^Binary"

# ncdm_species_sorted_ should be gone:
grep -rn "ncdm_species_sorted_" source/ species/

# index_pt_psi0_ncdm1 should be gone:
grep -rn "index_pt_psi0_ncdm1" source/ species/

# pba->N_ncdm/Omega0_ncdm_tot/N_decay_dr should be gone:
grep -rn "pba->N_ncdm\|pba->Omega0_ncdm_tot\|pba->N_decay_dr" source/ species/

# Hardcoded perturb_vector species fields should be gone:
grep -rn "pv->index_pt_phi\b\|pv->index_pt_delta_g\b\|pv->index_pt_delta_b\b" source/ species/
```

All five greps must return ZERO matches.

- [ ] **Step 2: Run regression cycle one more time**

`<regression-check>`. The Phase F column-title rename is the only knowingly-divergent task; manually inspect affected `.dat` headers.

- [ ] **Step 3: Run full pytest**

```bash
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py
```

Must pass. If any scenario fails, fix before PR.

- [ ] **Step 4: Push branch and open PR**

```bash
git push -u origin 267-remove-ncdm-id
gh pr create \
  --title "Encapsulate species perturbation layouts; remove ncdm_id" \
  --body "$(cat <<'EOF'
Resolves #267.

## Summary
- Per-species `PerturbLayout` polymorphic subclasses replace hardcoded `pv->index_pt_*` fields and the NCDM flat arrays. `perturb_vector` holds them in `species_layouts` parallel to `all_species_`.
- Per-thread isolation comes from each thread's own `pv` allocation inside the OpenMP region.
- `ncdm_id` integer is fully removed (along with `pba->N_ncdm`, `pba->Omega0_ncdm_tot`, `pba->N_decay_dr`, `ncdm_species_sorted_`, `index_pt_psi0_ncdm1`, `RescaledNCDMPerturbations(int n_ncdm,…)` → `(BaseSpecies*,…)`).
- DNCDM_DR composite key renamed from `DNCDM_DR_<id>` to `<dncdm_instance_name>` (and DR child to `<dncdm_instance_name>_DR`).
- Output column titles rename from `(.)rho_ncdm[<id>]` etc. to instance-name based — **breaking change for downstream parsers**.
- Shooter slot allocation moves from hardcoded InputModule cases to per-species `RegisterShootingIndices` hooks.

## Architecture
After this work, every species fully owns its perturbation layout; multi-instance support exists at the architecture level (input parser still single-instance for non-NCDM, follow-up).

## Test plan
- [x] `make class -j` clean build, no new warnings
- [x] `./class explanatory.ini` runs successfully
- [x] All five regression scenarios in `test/scenarios/` pass with default thread count
- [x] `TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py` passes
- [x] Numerical regression diff against master baseline within `1e-10` for non-renamed columns
- [x] grep -rn "ncdm_id" source/ species/ include/ → 0 matches

## Migration notes
- Downstream parsers matching `\\bncdm\\[\\d+\\]\\b` in column titles must update to match `<instance_name>` (e.g., `nu1`, `dncdm1`).
- C-API consumers reading `pba->N_ncdm` etc. must compute the count from species.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 5: Note the PR URL** — return it to the user.

---

## Plan summary

41 tasks total:

- **Phase 0** (Task 1): branch / scenario / baseline verification.
- **Phase A** (Tasks 2-9): infrastructure — `PerturbLayout` base + `CreatePerturbLayout` virtual; `SpeciesCollection::operator[]` + `index_of`; `species_layouts` on `perturb_vector`; layout-aware Register virtuals; module-side dual-call dispatch; Lambda migration (pattern check); layout-aware read virtuals + `PerturbSwitchContext`.
- **Phase B** (Tasks 10-23): per-species migrations — Fluid, CDM, DCDM, DR, ScalarField, UR, IDR, IDM_DR/IDM_DRMD, Baryons, Photons, NCDM-family (NCDM+DNCDM+NCDMInteracting), DCDM_DR, DNCDM_DR, IDM_DR_IDR+IDM_DRMD_IDR_DRMD composites.
- **Phase C** (Tasks 24-26): wire bookkeeping hooks + delete legacy switch-copy + tensor blocks.
- **Phase D** (Tasks 27-29): delete migrated bare fields from `perturb_vector` (3 groups).
- **Phase E** (Tasks 30-34): `ncdm_id`, `pba->N_ncdm/Omega0_ncdm_tot/N_decay_dr`, `ncdm_species_sorted_`, `RescaledNCDMPerturbations` refactor.
- **Phase F** (Task 35): output column title rename.
- **Phase G** (Task 36): DNCDM_DR composite-key rename.
- **Phase H** (Tasks 37-40): shooter hooks.
- **Phase I** (Task 41): final audit + PR.

## Spec coverage check

| Spec section | Implementing task(s) |
|---|---|
| Goal 1 (per-species layouts) | Tasks 2, 8, 10-23 |
| Goal 2 (`species_layouts` on pv, parallel-vector dispatch) | Tasks 3, 4, 7, 24, 26 |
| Goal 3 (all species methods take layout) | Task 9, Tasks 10-23 |
| Goal 4 (delete hardcoded `pv->*` fields) | Tasks 27-29 |
| Goal 5 (delete `ncdm_id` and family) | Tasks 30-34 |
| Goal 6 (instance-name column titles + DNCDM_DR rename) | Tasks 35, 36 |
| Goal 7 (shooter hooks) | Tasks 37-40 |
| Goal 8 (`RescaledNCDMPerturbations(BaseSpecies*)`) | Task 34 |
| Composite nested composition | Tasks 21, 22, 23 |
| Bookkeeping hooks (Stash removed; Copy + Mark) | Task 9 (virtuals), Task 24 (wiring), Task 25 (delete legacy); NCDM override in Task 20, UR/IDR/Photons overrides in 15/16/19 |
| Regression methodology (clean output, exit code, multi-thread) | `<regression-check>` defined at top, referenced in every behaviour-touching task |
| Risk: per-thread storage on pv | Task 4 (architecture), regression-check enforces multi-thread default |
| Risk: column-title rename ABI break | Task 35 documents and tests; PR body has migration note |
| Risk: composite layout-routing mistakes | Tasks 21-23; `dncdm_dr` and `combined` regression scenarios cover them |
