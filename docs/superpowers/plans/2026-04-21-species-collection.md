# SpeciesCollection Refactor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace `std::map<std::string, std::unique_ptr<BaseSpecies>> all_species_` with a purpose-built `SpeciesCollection` backed by a sorted vector, and cache direct pointers for the always-present `Photons`/`Baryons` species.

**Architecture:** Introduce a `SpeciesCollection` value type owning `std::vector<Entry>` where `Entry = { std::string key; std::unique_ptr<BaseSpecies> species; }`, exposing the subset of the map API that is actually used (`count`, `at`, iteration, `empty`, `size`) plus new accessors (`find` returning pointer-or-nullptr, `photons()`, `baryons()`). Construction-time mutation is confined to `InputModule::ConstructSpecies()`; after `freeze()`, the container is sorted by key and the hot-path pointers are resolved. This is an API-compatible refactor for the ~464 keyed-lookup sites (`.count("X")`, `.at("X")`) — only the ~43 iteration sites using structured bindings (`for (auto& [name, sp] : all_species_)`) need mechanical changes. The deeper cleanup of those keyed lookups into polymorphic calls is out of scope and already tracked in `docs/superpowers/plans/2026-04-15-remove-has-guards.md`.

**Naming / key policy (issue #246):** GitHub issue AarhusCosmology/CLASSpp#246 plans to reinterpret `BaseSpecies::name_` as the *class/type identifier* (e.g. two `massive_ncdm` instances share `name_`), with unique instance identifiers coming from user input (`nu1`, `ncdm__1`). `SpeciesCollection` therefore stores the lookup key explicitly alongside each species rather than deriving it from `sp->name()`. Today the caller passes the old map key (`"Photons"`, `"CDM"`, `"IDM_DR_IDR"`, ...) — these still match `sp->name()` for singletons. After #246 lands the caller simply passes the user-chosen identifier; no change to `SpeciesCollection` is required.

**Tech Stack:** C++17, existing CLASSpp build (Makefile / Xcode project), `class default.ini` + bbn output as the integration test.

---

## Preconditions

- Work in a dedicated git worktree or feature branch. Suggested branch name: `species-collection`.
- Establish a reference run before touching code:
  ```bash
  make class -j && ./class default.ini
  ```
  Save the generated `output/` files (or at minimum `*cl.dat`, `*pk.dat`, `*background.dat`) to a scratch directory, e.g. `/tmp/class_ref/`. Byte-for-byte equality with these files after each task is the regression test.
- Do NOT rebase onto or integrate with `refactor_perturbations.diff` during this refactor — it is a separate in-progress change that also touches `all_species_` and will need a manual merge afterwards. Flag this to the user at the end.

---

## File Structure

**Create:**
- `species/species_collection.h` — declares `SpeciesCollection` class
- `species/species_collection.cpp` — implements it

**Modify:**
- `source/input_module.h` — change `all_species_` type; include new header
- `source/input_module.cpp` — switch insertion from `map[name] = ...` to `.insert(...)`; add `freeze()` call at end of `ConstructSpecies()`; migrate 4 iteration sites
- `source/base_module.h` — change the reference type of `all_species_`; drop `<map>` include if unused
- `source/background_module.cpp` — migrate ~6 iteration sites
- `source/perturbations_module.cpp` — migrate ~32 iteration sites; convert hot-path `all_species_.at("Photons")`/`.at("Baryons")` calls in `Rho(...)` accumulators to `all_species_.photons()`/`baryons()`
- `source/nonlinear_module.cpp` — migrate 1 iteration site
- `species/base_species.h` — update the docstring that still references `std::map`
- Xcode project / Makefile — register the two new source/header files if the build system does not auto-discover them

---

## Task 1: Create `SpeciesCollection` skeleton (header)

**Files:**
- Create: `species/species_collection.h`

- [ ] **Step 1: Write the header**

```cpp
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "base_species.h"

/**
 * Ordered, named collection of cosmological species.
 *
 * Mutation is restricted to the construction phase (InputModule::ConstructSpecies).
 * Each species is inserted with an explicit lookup key — today this is the
 * hard-coded class identifier (e.g. "Photons", "CDM", "IDM_DR_IDR"); under the
 * object-based input syntax planned in issue #246 it will be the user-supplied
 * instance name (e.g. "nu1", "ncdm__1"). The key is NOT derived from
 * BaseSpecies::name() because post-#246 several instances may share the same
 * type name (e.g. two "massive_ncdm" instances).
 *
 * Call freeze() once all species have been inserted; after that the container
 * is sorted by key and the hot-path accessors photons()/baryons() are valid.
 *
 * Iteration yields `Entry&`, where Entry exposes .key and .species. The common
 * form
 *     for (auto& sp : all_species_) { sp->Foo(); }
 * is supported via an Entry that acts like a unique_ptr<BaseSpecies>& (operator->
 * and operator* forward to the owned species). Use `sp.key` when the lookup key
 * is needed; use `sp->name()` only when the BaseSpecies type name is wanted
 * (these differ once #246 lands).
 */
class SpeciesCollection {
 public:
  using Ptr = std::unique_ptr<BaseSpecies>;

  struct Entry {
    std::string key;
    Ptr         species;

    BaseSpecies*       operator->()       { return species.get(); }
    const BaseSpecies* operator->() const { return species.get(); }
    BaseSpecies&       operator*()        { return *species; }
    const BaseSpecies& operator*()  const { return *species; }
    BaseSpecies*       get()              { return species.get(); }
    const BaseSpecies* get()        const { return species.get(); }
    explicit operator bool() const { return static_cast<bool>(species); }
  };

  using Container = std::vector<Entry>;

  SpeciesCollection()                                    = default;
  SpeciesCollection(const SpeciesCollection&)            = delete;
  SpeciesCollection& operator=(const SpeciesCollection&) = delete;
  SpeciesCollection(SpeciesCollection&&)                 = default;
  SpeciesCollection& operator=(SpeciesCollection&&)      = default;

  // ── Construction-phase mutation ─────────────────────────────────────────
  /** Insert a species under an explicit lookup key. Must be called before
   *  freeze(). Duplicate keys abort in debug builds. */
  void insert(std::string key, Ptr species);

  /** Sort by name and resolve cached Photons/Baryons pointers.
   *  Must be called exactly once, after all insert() calls. */
  void freeze();

  // ── Keyed lookup (replaces std::map API) ────────────────────────────────
  /** Returns 1 if a species with this name is present, else 0.
   *  Kept as size_t/bool-compatible to match std::map::count call sites. */
  std::size_t count(const std::string& name) const;

  /** Access by name; throws std::out_of_range if absent. */
  BaseSpecies&       at(const std::string& name);
  const BaseSpecies& at(const std::string& name) const;

  /** Pointer-or-nullptr lookup (no exception). */
  BaseSpecies*       find(const std::string& name);
  const BaseSpecies* find(const std::string& name) const;

  // ── Hot-path cached accessors ───────────────────────────────────────────
  /** Always-present species; valid only after freeze(). */
  BaseSpecies&       photons()        { return *photons_; }
  const BaseSpecies& photons() const  { return *photons_; }
  BaseSpecies&       baryons()        { return *baryons_; }
  const BaseSpecies& baryons() const  { return *baryons_; }

  // ── Iteration / size ────────────────────────────────────────────────────
  Container::iterator       begin()        { return species_.begin(); }
  Container::iterator       end()          { return species_.end(); }
  Container::const_iterator begin()  const { return species_.begin(); }
  Container::const_iterator end()    const { return species_.end(); }

  bool        empty() const { return species_.empty(); }
  std::size_t size()  const { return species_.size(); }

 private:
  Container    species_;
  BaseSpecies* photons_ = nullptr;
  BaseSpecies* baryons_ = nullptr;
  bool         frozen_  = false;
};
```

- [ ] **Step 2: Commit**

```bash
git add species/species_collection.h
git commit -m "species: add SpeciesCollection header"
```

---

## Task 2: Implement `SpeciesCollection`

**Files:**
- Create: `species/species_collection.cpp`

- [ ] **Step 1: Write the implementation**

```cpp
#include "species_collection.h"

#include <algorithm>
#include <cassert>
#include <stdexcept>

void SpeciesCollection::insert(std::string key, Ptr species) {
  assert(!frozen_ && "SpeciesCollection: insert() after freeze()");
  assert(species && "SpeciesCollection: null species");
  for (const auto& e : species_) {
    assert(e.key != key && "SpeciesCollection: duplicate key");
  }
  species_.push_back(Entry{std::move(key), std::move(species)});
}

void SpeciesCollection::freeze() {
  assert(!frozen_ && "SpeciesCollection: freeze() called twice");
  std::sort(species_.begin(), species_.end(),
            [](const Entry& a, const Entry& b) { return a.key < b.key; });
  photons_ = find("Photons");
  baryons_ = find("Baryons");
  assert(photons_ && "SpeciesCollection: Photons species missing at freeze()");
  assert(baryons_ && "SpeciesCollection: Baryons species missing at freeze()");
  frozen_ = true;
}

namespace {
// Linear search is fine: species counts are ~5–15.
template <class Vec>
auto find_impl(Vec& v, const std::string& key)
    -> decltype(v[0].species.get()) {
  for (auto& e : v)
    if (e.key == key) return e.species.get();
  return nullptr;
}
}  // namespace

std::size_t SpeciesCollection::count(const std::string& key) const {
  return find(key) ? 1u : 0u;
}

BaseSpecies* SpeciesCollection::find(const std::string& key) {
  return find_impl(species_, key);
}
const BaseSpecies* SpeciesCollection::find(const std::string& key) const {
  return find_impl(species_, key);
}

BaseSpecies& SpeciesCollection::at(const std::string& key) {
  if (auto* p = find(key)) return *p;
  throw std::out_of_range("SpeciesCollection::at: no species with key '" + key + "'");
}
const BaseSpecies& SpeciesCollection::at(const std::string& key) const {
  if (auto* p = find(key)) return *p;
  throw std::out_of_range("SpeciesCollection::at: no species with key '" + key + "'");
}
```

- [ ] **Step 2: Register the files in the build**

Check whether `species/*.cpp` is auto-globbed by the Makefile:
```bash
grep -E 'species/.*\.(cpp|o)' Makefile
```
- If globbed: nothing to do.
- If listed explicitly: add `species/species_collection.o` (and the `.cpp` source) next to the other species entries.

For Xcode:
```bash
grep -c 'species_collection' CLASS.xcodeproj/project.pbxproj
```
If 0, add the files via Xcode GUI (or defer until the user runs once — the Makefile build is primary).

- [ ] **Step 3: Compile standalone (new file only)**

```bash
make class -j 2>&1 | tail -20
```
Expected: build succeeds; no linker call-sites yet so the new file's symbols are unused but present.

- [ ] **Step 4: Commit**

```bash
git add species/species_collection.cpp Makefile  # (and pbxproj if changed)
git commit -m "species: implement SpeciesCollection (sorted-vector backing store)"
```

---

## Task 3: Switch `InputModule` to own a `SpeciesCollection`

**Files:**
- Modify: `source/input_module.h:57-58`
- Modify: `source/input_module.cpp` (insertion sites at lines 207–265; iteration sites at 3653, 3711, 3802, 3852; `.at(...)` sites at 3593, 3607, 3623, 3632)

- [ ] **Step 1: Change the member type**

In `source/input_module.h`, replace:

```cpp
  std::map<std::string, std::unique_ptr<BaseSpecies>> all_species_;
```

with:

```cpp
  SpeciesCollection all_species_;
```

Add the include near the other species include (top of file):
```cpp
#include "../species/species_collection.h"
```
and remove `#include <map>` from this header if nothing else needs it (grep first).

- [ ] **Step 2: Convert insertion sites in `ConstructSpecies()`**

In `source/input_module.cpp`, every line of the form

```cpp
all_species_["Name"] = std::make_unique<SomeSpecies>(...);
```

becomes

```cpp
all_species_.insert("Name", std::make_unique<SomeSpecies>(...));
```

For the move-in variants at lines ~234, 242, 250:
```cpp
all_species_[name] = std::move(sp);
all_species_[name] = std::make_unique<DNCDM_DR_Species>(std::move(dncdm_sp), pba, nullptr);
all_species_[name] = std::move(int_sp);
```
become
```cpp
all_species_.insert(name, std::move(sp));
all_species_.insert(name, std::make_unique<DNCDM_DR_Species>(std::move(dncdm_sp), pba, nullptr));
all_species_.insert(name, std::move(int_sp));
```

The first argument is the lookup key — keep the exact string that was being used as the map key today. (Under issue #246 this will later be replaced by user-supplied instance names, but that is not in scope here.)

- [ ] **Step 3: Call `freeze()` at the end of `ConstructSpecies()`**

Append `all_species_.freeze();` as the last line of `InputModule::ConstructSpecies()`.

- [ ] **Step 4: Migrate input-module iteration sites**

At lines 3653, 3711, 3802, 3852, replace

```cpp
for (auto& [key, sp] : all_species_) {          // or input_module->all_species_ / bam->all_species_
  ... sp->Foo() ...
}
```

with

```cpp
for (auto& sp : all_species_) {
  ... sp->Foo() ...
}
```

If any body uses `key`, substitute `sp->name()`. Run:
```bash
grep -n 'key' source/input_module.cpp | head
```
to eyeball each migrated loop body.

- [ ] **Step 5: Build**

```bash
make class -j 2>&1 | tail -40
```
Expected: builds cleanly. Remaining `all_species_.count(...)`, `all_species_.at(...)` calls elsewhere in the project are API-compatible and still compile.

- [ ] **Step 6: Regression test**

```bash
./class default.ini && diff -r output/ /tmp/class_ref/
```
Expected: identical output.

- [ ] **Step 7: Commit**

```bash
git add source/input_module.h source/input_module.cpp
git commit -m "input: own all_species_ as SpeciesCollection"
```

---

## Task 4: Update `BaseModule`'s reference type

**Files:**
- Modify: `source/base_module.h:43`

- [ ] **Step 1: Update the reference and include**

Replace
```cpp
#include <map>
...
  const std::map<std::string, std::unique_ptr<BaseSpecies>>& all_species_;
```
with
```cpp
#include "../species/species_collection.h"
...
  const SpeciesCollection& all_species_;
```
Remove `#include <map>` if nothing else in the header needs it.

- [ ] **Step 2: Update the docstring on `BaseSpecies` that still mentions `std::map`**

In `species/base_species.h` (~line 29), replace the paragraph:
```
 * The map in BaseModule is: const std::map<std::string, std::unique_ptr<BaseSpecies>>.
 * Use .at("CDM") – never operator[] – to maintain const correctness.
```
with:
```
 * The collection in BaseModule is a const SpeciesCollection&. Use .at("CDM")
 * (throws if absent) or .find("CDM") (returns pointer-or-nullptr) for keyed
 * lookup; iterate with `for (auto& sp : all_species_)` and use sp->name() if
 * the name is needed.
```

Update the analogous `Use .at() – never operator[]` comment on `BaseModule::all_species_` in `source/base_module.h` to drop the "keyed by name" phrasing (the collection is still named-lookup-capable but is no longer a map).

- [ ] **Step 3: Build**

```bash
make class -j 2>&1 | tail -40
```
Expected: success — all `.count()`/`.at()` call sites across the project resolve against the new type.

- [ ] **Step 4: Regression test**

```bash
./class default.ini && diff -r output/ /tmp/class_ref/
```
Expected: identical output.

- [ ] **Step 5: Commit**

```bash
git add source/base_module.h species/base_species.h
git commit -m "base_module: expose all_species_ as const SpeciesCollection&"
```

---

## Task 5: Migrate iteration sites in `background_module.cpp`

**Files:**
- Modify: `source/background_module.cpp` at lines 759, 775, 898, 1499, 1699, 1734 (survey with grep first — line numbers shift).

- [ ] **Step 1: Locate all iteration sites**

```bash
grep -n 'for (auto& \[name, sp\] : all_species_' source/background_module.cpp
grep -n 'for (auto& \[key, sp\] : all_species_'  source/background_module.cpp
```

- [ ] **Step 2: Migrate each site**

For each hit: replace
```cpp
for (auto& [name, sp] : all_species_) {
  ... uses of `name` and `sp` ...
}
```
with
```cpp
for (auto& sp : all_species_) {
  ... uses of `sp.key` and `sp` ...
}
```

`sp` is now a `SpeciesCollection::Entry&`; `sp->Foo()` and `*sp` still forward to the owned `BaseSpecies`. If the original loop read `name` to compare against hard-coded strings (e.g. `if (name == "DCDM_DR")`), that is a lookup-key comparison, so use `sp.key`. Only use `sp->name()` if the code genuinely wants the species type identifier (post-#246 these diverge).

Example (line ~898): `if (name == "DCDM_DR") continue;` becomes `if (sp.key == "DCDM_DR") continue;`.
Example (line ~1499): no use of `name` — just drop the binding.

- [ ] **Step 3: Build**

```bash
make class -j 2>&1 | tail -20
```
Expected: clean build.

- [ ] **Step 4: Regression test**

```bash
./class default.ini && diff -r output/ /tmp/class_ref/
```
Expected: identical output.

- [ ] **Step 5: Commit**

```bash
git add source/background_module.cpp
git commit -m "background_module: migrate species iteration to SpeciesCollection"
```

---

## Task 6: Migrate iteration sites in `perturbations_module.cpp`

**Files:**
- Modify: `source/perturbations_module.cpp` (~32 sites).

- [ ] **Step 1: Locate**

```bash
grep -n 'for (auto& \[name, sp\] : all_species_' source/perturbations_module.cpp
```

- [ ] **Step 2: Migrate each site** using the same `[name, sp]` → `sp` + `sp->name()` pattern as Task 5.

Pay attention to any site that uses `name` in a conditional body — substitute `sp->name()`. Watch for lambdas capturing `name` (unlikely based on the survey, but grep the body of each loop).

- [ ] **Step 3: Build**

```bash
make class -j 2>&1 | tail -20
```

- [ ] **Step 4: Regression test**

```bash
./class default.ini && diff -r output/ /tmp/class_ref/
```
Expected: identical output.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations_module: migrate species iteration to SpeciesCollection"
```

---

## Task 7: Migrate the remaining iteration site in `nonlinear_module.cpp`

**Files:**
- Modify: `source/nonlinear_module.cpp:1008`.

- [ ] **Step 1: Migrate**

Replace
```cpp
for (auto& [name, sp] : all_species_) {
  auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
  ...
}
```
with
```cpp
for (auto& sp : all_species_) {
  auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
  ...
}
```

- [ ] **Step 2: Build + regression test**

```bash
make class -j && ./class default.ini && diff -r output/ /tmp/class_ref/
```
Expected: identical output.

- [ ] **Step 3: Verify no iteration sites remain**

```bash
grep -rn 'for (auto& \[.*, sp\] : .*all_species_' source/ species/
```
Expected: no hits.

- [ ] **Step 4: Commit**

```bash
git add source/nonlinear_module.cpp
git commit -m "nonlinear_module: migrate species iteration to SpeciesCollection"
```

---

## Task 8: Wire up hot-path `photons()` / `baryons()` accessors

**Files:**
- Modify: `source/perturbations_module.cpp` — the ~39 sites using `all_species_.at("Photons")` / `all_species_.at("Baryons")`.
- Modify: `source/background_module.cpp` — any analogous sites (grep first).

The goal is a cheap mechanical swap: `all_species_.at("Photons")` → `all_species_.photons()` (same for Baryons). The returned reference type is identical to dereferencing `unique_ptr::operator*`, so `.Rho(...)` etc. calls work unchanged.

- [ ] **Step 1: Locate**

```bash
grep -n 'all_species_\.at("Photons")' source/perturbations_module.cpp source/background_module.cpp source/thermodynamics_module.cpp
grep -n 'all_species_\.at("Baryons")' source/perturbations_module.cpp source/background_module.cpp source/thermodynamics_module.cpp
```

- [ ] **Step 2: Rewrite each site**

For expressions of the form `all_species_.at("Photons")->Rho(...)` (which return via `operator->` on `unique_ptr`), the new form is `all_species_.photons().Rho(...)` — note the `.` not `->` because `photons()` returns a reference, not a pointer.

For expressions already in reference form like `*all_species_.at("Photons")` or `const auto& PH = all_species_.at("Photons");` (which is actually a `unique_ptr&` in the old code — keep an eye out), rewrite to `all_species_.photons()` and `const auto& PH = all_species_.photons();` respectively.

When uncertain about a particular site, compile — the type system will catch `->` vs `.` mistakes immediately.

- [ ] **Step 3: Build**

```bash
make class -j 2>&1 | tail -40
```
Expected: clean build. Fix any `->` vs `.` mismatches flagged by the compiler.

- [ ] **Step 4: Regression test**

```bash
./class default.ini && diff -r output/ /tmp/class_ref/
```
Expected: identical output.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp source/background_module.cpp source/thermodynamics_module.cpp
git commit -m "modules: use SpeciesCollection::photons()/baryons() on hot paths"
```

---

## Task 9: Broader regression sweep

**Files:** none modified.

- [ ] **Step 1: Run a non-default input file**

Pick one of the bundled `.ini` files exercising more species:
```bash
for ini in explanatory.ini idm_dr.ini idm_drmd_idr_drmd.ini dcdm_dr.ini scalar_field.ini; do
  test -f $ini && (./class $ini > /tmp/class_$ini.log 2>&1 && echo OK $ini || echo FAIL $ini)
done
```
Expected: all OK (file-not-found is acceptable; failures are not).

- [ ] **Step 2: Compare outputs against a reference run on master**

From a clean checkout of master in a sibling worktree, run the same inis and `diff -r` each output directory against this branch. Expected: identical.

- [ ] **Step 3: Sanity-grep for stragglers**

```bash
grep -rn 'std::map<std::string, std::unique_ptr<BaseSpecies>>' source/ species/
```
Expected: no hits.

```bash
grep -rn 'for (auto& \[.*, sp\] : .*all_species_' source/ species/
```
Expected: no hits.

- [ ] **Step 4: Final commit (docs if anything drifted)**

Only if something needed adjusting:
```bash
git add -u
git commit -m "species-collection: docs + straggler cleanup"
```

---

## Out of Scope / Follow-ups

- The ~232 `.count("SpeciesName")` guards and ~70 `.at("SpeciesName")` calls paired with them remain after this refactor. They still work via the new API, but they are the species-specific-logic-in-main-code smell the user flagged. Their removal is already tracked in `docs/superpowers/plans/2026-04-15-remove-has-guards.md` — do not expand this plan.
- `refactor_perturbations.diff` in the repo root touches many `all_species_` iteration sites and will need a manual merge after this plan lands. Flag to the user at the end; do not try to apply it.
- A future enhancement could expose cached pointers for other always-present species (e.g. `lambda()` when present) or replace `find`/`at` with a lightweight `std::string_view`-keyed API. Not needed now.
