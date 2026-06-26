# Dot-notation for single-instance legacy species — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let the single-instance legacy species (photons, baryons, cdm, ur, lambda, fluid) accept dot-notation (`bary.type=baryons`, `bary.Omega=0.05`) by translating it to legacy keys (`Omega_b=0.05`) in an input-preprocessing pass, producing byte-identical output.

**Architecture:** Two parts in one PR. **Part A** (output-neutral refactor) gives each species a canonical `kTypeName` string constant and makes `kAllSpeciesFactories` and the `instances_with("type", …)` call sites source the type string from it — one source of truth, owned by the species. **Part B** adds `TranslateSingleInstanceDotSyntax(FileContent*)` in `species/species_input.cpp`, a module-owned table keyed by `kTypeName` that rewrites recognised `<name>.<dot-field>` entries to legacy keys, with single-instance and conflict guards, called once before all input reads.

**Tech Stack:** C++17, CMake + ctest (plain `assert`-based test executables linked against the `classpp` static lib), Python/Cython wrapper via scikit-build-core, scenario `.ini` files in `test/scenarios/`.

## Global Constraints

- Spec: `docs/superpowers/specs/2026-06-26-dot-notation-single-instance-species-design.md`.
- **Output bar: byte-identical (`cmp`-clean).** This is a translation + rename, not arithmetic. Part A scenarios must be `cmp`-identical to master; a dot-form `.ini` must be `cmp`-identical to the equivalent legacy-form `.ini`.
- `kTypeName` is the factory/registry/`type`-input identifier. It is **distinct from** the species' runtime `name_` (display/output name, e.g. `BaseSpecies("UR", …)`). **Never touch `name_`** — output columns derive from it.
- `kTypeName` values (the complete set, = the 14 `kAllSpeciesFactories` entries): `photons`, `baryons`, `cdm`, `ur`, `dcdm_dr`, `ncdm_standard`, `ncdm_greybody`, `ncdm_decay_dr`, `ncdm_self_interacting`, `idm_dr_idr`, `idm_drmd_idr_drmd`, `lambda`, `fluid`, `scalar_field`.
- `GreyBodyNCDMSpecies` and `NCDMInteractingSpecies` inherit `NCDMSpecies`; each **must declare its own** `kTypeName` (otherwise `Derived::kTypeName` silently resolves to the inherited `"ncdm_standard"`).
- ScalarField gets a `kTypeName` (Part A) but is **excluded** from Part B's translation table.
- Never hand-edit `cclassy.pxd` (auto-generated). Never `git add -A` (the tree has unrelated untracked artifacts) — stage explicit paths.
- Branch already exists: `single-instance-dot-syntax` (the spec commit is on it).
- Build: `cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j`. Tests: `ctest --test-dir build --output-on-failure`. The `./class` binary lands in the repo root.

---

### Task 1: Baseline

Establish the byte-identical reference. The branch currently differs from master only by the spec doc, so the baseline equals master's output.

**Files:** none modified.

- [ ] **Step 1: Confirm branch + clean build**

Run:
```bash
cd /Users/au192734/Projects/class_claude
git branch --show-current        # expect: single-instance-dot-syntax
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```
Expected: builds clean; `./class` exists in the repo root.

- [ ] **Step 2: Record baseline outputs for the representative scenarios**

Run:
```bash
mkdir -p /tmp/dotsyn/base
for s in gauge_lcdm gauge_fluid gauge_ncdm gauge_dcdm dncdm_dr; do
  rm -f output/*
  ./class test/scenarios/$s.ini || echo "RUN FAILED: $s"
  mkdir -p /tmp/dotsyn/base/$s && mv output/* /tmp/dotsyn/base/$s/ 2>/dev/null
done
ls /tmp/dotsyn/base/*/
```
Expected: each directory holds `.dat` files. (These cover photons/baryons/cdm/lambda + ur, fluid, ncdm, dcdm composite, dncdm composite — i.e. every type-discovered factory plus the single-instance species.)

- [ ] **Step 3: Record baseline ctest**

Run:
```bash
ctest --test-dir build --output-on-failure
```
Expected: `test-parser`, `test-bisection`, `test-photons` all pass.

---

### Task 2: Part A — `kTypeName` constants + factory/type-string unification

Output-neutral. Add a `kTypeName` to every species, route the factory table, the four `instances_with("type", …)` literals, the `input_module` synthesis call, `ClosureSpeciesName`, and `is_ncdm_family` through it. Test-first via a static registry-consistency test, then verify `cmp`-identical scenarios.

**Files:**
- Create: `species/species_type_name_test.cpp`
- Modify: `CMakeLists.txt:210-216` (register the new test)
- Modify (add `kTypeName`): `species/photons.h`, `species/baryons.h`, `species/cdm.h`, `species/ultra_relativistic.h`, `species/dcdm_dr_species.h`, `species/ncdm_species.h`, `species/greybody_ncdm_species.h`, `species/dncdm_species.h`, `species/ncdm_interacting_species.h`, `species/idm_dr_idr_species.h`, `species/idm_drmd_idr_drmd_species.h`, `species/lambda.h`, `species/fluid.h`, `species/scalar_field.h`
- Modify (route literals): `species/all_species.h:40-70`, `species/ncdm_species.cpp:205`, `species/greybody_ncdm_species.cpp:298`, `species/dncdm_species.cpp:226`, `species/ncdm_interacting_species.cpp:85`, `source/input_module.cpp:259-260,594`

**Interfaces:**
- Produces: `static constexpr std::string_view <Species>::kTypeName` on all 14 factory species, with the values listed in Global Constraints. Consumed by Task 3's translation table and by `kAllSpeciesFactories`.

- [ ] **Step 1: Write the failing static registry test**

Create `species/species_type_name_test.cpp`:
```cpp
// Pins each species' canonical type string and the kAllSpeciesFactories
// registry to the agreed set, so a rename or a missing override is caught.
#include "species/all_species.h"

#include <cassert>
#include <set>
#include <string>

int main() {
  // Per-species constants (the inherited-shadowing cases matter: GreyBody and
  // NCDMInt derive from NCDMSpecies and must declare their own kTypeName).
  assert(PhotonsSpecies::kTypeName == "photons");
  assert(BaryonsSpecies::kTypeName == "baryons");
  assert(CDMSpecies::kTypeName == "cdm");
  assert(UltraRelativisticSpecies::kTypeName == "ur");
  assert(DCDM_DR_Species::kTypeName == "dcdm_dr");
  assert(NCDMSpecies::kTypeName == "ncdm_standard");
  assert(GreyBodyNCDMSpecies::kTypeName == "ncdm_greybody");
  assert(DNCDMSpecies::kTypeName == "ncdm_decay_dr");
  assert(NCDMInteractingSpecies::kTypeName == "ncdm_self_interacting");
  assert(IDM_DR_IDR_Species::kTypeName == "idm_dr_idr");
  assert(IDM_DRMD_IDR_DRMD_Species::kTypeName == "idm_drmd_idr_drmd");
  assert(LambdaSpecies::kTypeName == "lambda");
  assert(FluidSpecies::kTypeName == "fluid");
  assert(ScalarFieldSpecies::kTypeName == "scalar_field");

  // The factory registry carries exactly these type strings.
  std::set<std::string> names;
  for (const auto& e : kAllSpeciesFactories) {
    names.insert(std::string(e.name));
  }
  const std::set<std::string> expected = {
      "photons", "baryons", "cdm", "ur", "dcdm_dr", "ncdm_standard",
      "ncdm_greybody", "ncdm_decay_dr", "ncdm_self_interacting", "idm_dr_idr",
      "idm_drmd_idr_drmd", "lambda", "fluid", "scalar_field",
  };
  assert(names == expected);

  // ClosureSpeciesName routes through the same strings.
  assert(ClosureSpeciesName(ClosureSpecies::Lambda) == "lambda");
  assert(ClosureSpeciesName(ClosureSpecies::Fluid) == "fluid");
  assert(ClosureSpeciesName(ClosureSpecies::ScalarField) == "scalar_field");
  return 0;
}
```

- [ ] **Step 2: Register the test in CMake**

In `CMakeLists.txt`, add the executable after line 212 and the name to the `foreach` ITEMS list (line 213):
```cmake
  add_executable(test-parser tools/parser_test.cpp)
  add_executable(test-bisection tools/bisection_test.cpp)
  add_executable(test-photons species/photons_formula_test.cpp)
  add_executable(test-species-types species/species_type_name_test.cpp)
  foreach(_t IN ITEMS test-parser test-bisection test-photons test-species-types)
```

- [ ] **Step 3: Verify it fails to compile (kTypeName undefined)**

Run:
```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target test-species-types -j 2>&1 | tail -15
```
Expected: compile error — `kTypeName is not a member of PhotonsSpecies` (or similar). This confirms the test exercises the new constants.

- [ ] **Step 4: Add the `kTypeName` constant to each species header**

In the `public:` section of each class, add the line (ensure the header has `#include <string_view>` — add it near the other includes if absent). Exact constant per file:

| File | Class | Line to add |
|---|---|---|
| `species/photons.h` | `PhotonsSpecies` | `static constexpr std::string_view kTypeName = "photons";` |
| `species/baryons.h` | `BaryonsSpecies` | `static constexpr std::string_view kTypeName = "baryons";` |
| `species/cdm.h` | `CDMSpecies` | `static constexpr std::string_view kTypeName = "cdm";` |
| `species/ultra_relativistic.h` | `UltraRelativisticSpecies` | `static constexpr std::string_view kTypeName = "ur";` |
| `species/dcdm_dr_species.h` | `DCDM_DR_Species` | `static constexpr std::string_view kTypeName = "dcdm_dr";` |
| `species/ncdm_species.h` | `NCDMSpecies` | `static constexpr std::string_view kTypeName = "ncdm_standard";` |
| `species/greybody_ncdm_species.h` | `GreyBodyNCDMSpecies` | `static constexpr std::string_view kTypeName = "ncdm_greybody";` |
| `species/dncdm_species.h` | `DNCDMSpecies` | `static constexpr std::string_view kTypeName = "ncdm_decay_dr";` |
| `species/ncdm_interacting_species.h` | `NCDMInteractingSpecies` | `static constexpr std::string_view kTypeName = "ncdm_self_interacting";` |
| `species/idm_dr_idr_species.h` | `IDM_DR_IDR_Species` | `static constexpr std::string_view kTypeName = "idm_dr_idr";` |
| `species/idm_drmd_idr_drmd_species.h` | `IDM_DRMD_IDR_DRMD_Species` | `static constexpr std::string_view kTypeName = "idm_drmd_idr_drmd";` |
| `species/lambda.h` | `LambdaSpecies` | `static constexpr std::string_view kTypeName = "lambda";` |
| `species/fluid.h` | `FluidSpecies` | `static constexpr std::string_view kTypeName = "fluid";` |
| `species/scalar_field.h` | `ScalarFieldSpecies` | `static constexpr std::string_view kTypeName = "scalar_field";` |

- [ ] **Step 5: Route `kAllSpeciesFactories` and `ClosureSpeciesName` through `kTypeName`**

In `species/all_species.h`, replace the factory array (lines 40-55) so each row's name is the species constant (note the `ncdm_decay_dr` row keeps the `DNCDM_DR_Species::CreateAll` factory):
```cpp
inline constexpr std::array kAllSpeciesFactories = {
    SpeciesFactoryEntry{PhotonsSpecies::kTypeName, &PhotonsSpecies::CreateAll},
    SpeciesFactoryEntry{BaryonsSpecies::kTypeName, &BaryonsSpecies::CreateAll},
    SpeciesFactoryEntry{CDMSpecies::kTypeName, &CDMSpecies::CreateAll},
    SpeciesFactoryEntry{UltraRelativisticSpecies::kTypeName, &UltraRelativisticSpecies::CreateAll},
    SpeciesFactoryEntry{DCDM_DR_Species::kTypeName, &DCDM_DR_Species::CreateAll},
    SpeciesFactoryEntry{NCDMSpecies::kTypeName, &NCDMSpecies::CreateAll},
    SpeciesFactoryEntry{GreyBodyNCDMSpecies::kTypeName, &GreyBodyNCDMSpecies::CreateAll},
    SpeciesFactoryEntry{DNCDMSpecies::kTypeName, &DNCDM_DR_Species::CreateAll},
    SpeciesFactoryEntry{NCDMInteractingSpecies::kTypeName, &NCDMInteractingSpecies::CreateAll},
    SpeciesFactoryEntry{IDM_DR_IDR_Species::kTypeName, &IDM_DR_IDR_Species::CreateAll},
    SpeciesFactoryEntry{IDM_DRMD_IDR_DRMD_Species::kTypeName, &IDM_DRMD_IDR_DRMD_Species::CreateAll},
    SpeciesFactoryEntry{LambdaSpecies::kTypeName, &LambdaSpecies::CreateAll},
    SpeciesFactoryEntry{FluidSpecies::kTypeName, &FluidSpecies::CreateAll},
    SpeciesFactoryEntry{ScalarFieldSpecies::kTypeName, &ScalarFieldSpecies::CreateAll},
};
```
Then replace the `ClosureSpeciesName` body (lines 58-70) literals:
```cpp
inline constexpr std::string_view ClosureSpeciesName(ClosureSpecies cs) {
  switch (cs) {
    case ClosureSpecies::Lambda:
      return LambdaSpecies::kTypeName;
    case ClosureSpecies::Fluid:
      return FluidSpecies::kTypeName;
    case ClosureSpecies::ScalarField:
      return ScalarFieldSpecies::kTypeName;
    case ClosureSpecies::None:
      return "";
  }
  return "";
}
```

- [ ] **Step 6: Route the four `instances_with` literals + the input_module synthesis**

Replace each hardcoded type literal with `std::string(kTypeName)` (the call site is inside the owning class's `CreateAll`, so unqualified `kTypeName` resolves to that class — including the GreyBody/NCDMInt overrides):

- `species/ncdm_species.cpp:205`:
  ```cpp
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  ```
- `species/greybody_ncdm_species.cpp:298`:
  ```cpp
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  ```
- `species/dncdm_species.cpp:226`:
  ```cpp
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  ```
- `species/ncdm_interacting_species.cpp:85`:
  ```cpp
  const auto instances = ctx.pfc->instances_with("type", std::string(kTypeName));
  ```
- `source/input_module.cpp:594` (outside any species class — qualify it):
  ```cpp
  const auto standard_ncdm_instances =
      file_content_.instances_with("type", std::string(NCDMSpecies::kTypeName));
  ```

- [ ] **Step 7: Route `is_ncdm_family` through `kTypeName`**

In `source/input_module.cpp:259-260`, replace the CamelCase literals:
```cpp
    const bool is_ncdm_family = (entry.name == NCDMSpecies::kTypeName ||
                                 entry.name == DNCDMSpecies::kTypeName ||
                                 entry.name == NCDMInteractingSpecies::kTypeName);
```
(`ncdm_greybody` stays out of the family counter, preserving current behavior.)

- [ ] **Step 8: Build and pass the static test**

Run:
```bash
cmake --build build -j
ctest --test-dir build -R test-species-types --output-on-failure
```
Expected: clean build; `test-species-types` passes.

- [ ] **Step 9: Verify output-neutrality (`cmp`-identical) + full ctest**

Run:
```bash
for s in gauge_lcdm gauge_fluid gauge_ncdm gauge_dcdm dncdm_dr; do
  rm -f output/*; ./class test/scenarios/$s.ini || echo "RUN FAILED: $s"
  mkdir -p /tmp/dotsyn/partA/$s && mv output/* /tmp/dotsyn/partA/$s/ 2>/dev/null
  for f in /tmp/dotsyn/base/$s/*.dat; do
    cmp -s "$f" /tmp/dotsyn/partA/$s/$(basename "$f") || echo "DIFFERS: $s/$(basename "$f")"
  done
done
ctest --test-dir build --output-on-failure
```
Expected: no `DIFFERS` lines (pure rename); all ctest green. If anything differs, a `name_` was touched or a type string changed a discovery path — investigate before continuing.

- [ ] **Step 10: Commit**

```bash
git add CMakeLists.txt species/species_type_name_test.cpp species/all_species.h \
  species/photons.h species/baryons.h species/cdm.h species/ultra_relativistic.h \
  species/dcdm_dr_species.h species/ncdm_species.h species/greybody_ncdm_species.h \
  species/dncdm_species.h species/ncdm_interacting_species.h species/idm_dr_idr_species.h \
  species/idm_drmd_idr_drmd_species.h species/lambda.h species/fluid.h species/scalar_field.h \
  species/ncdm_species.cpp species/greybody_ncdm_species.cpp species/dncdm_species.cpp \
  species/ncdm_interacting_species.cpp source/input_module.cpp
git commit -m "refactor: unify factory name with canonical species kTypeName"
```

---

### Task 3: Part B core — `TranslateSingleInstanceDotSyntax` + unit tests

The translation logic, unit-tested in isolation against a hand-built `FileContent`. Depends on Task 2's `kTypeName` constants.

**Files:**
- Modify: `species/species_input.h` (declare the function)
- Modify: `species/species_input.cpp` (table + implementation; add `#include "species/all_species.h"`)
- Test: `tools/parser_test.cpp` (add test functions + main calls — builds into the existing `test-parser`)

**Interfaces:**
- Consumes: `<Species>::kTypeName` (Task 2); `FileContent::{instances_with,get,set,was_read}`; `SpeciesInput`.
- Produces: `void TranslateSingleInstanceDotSyntax(FileContent* pfc);`

- [ ] **Step 1: Write the failing unit tests**

In `tools/parser_test.cpp`, add these functions (after the existing `SpeciesInput` tests):
```cpp
static void test_dot_baryons_translation() {
  FileContent fc;
  fc.set("atom.type", "baryons");
  fc.set("atom.Omega", "0.05");
  TranslateSingleInstanceDotSyntax(&fc);
  auto ob = fc.get<std::string>("Omega_b");
  assert(ob && *ob == "0.05");
  assert(fc.was_read("atom.type"));
  assert(fc.was_read("atom.Omega"));
}

static void test_dot_name_is_ignored() {
  // Any instance name maps to the same legacy key.
  FileContent fc;
  fc.set("whatever.type", "cdm");
  fc.set("whatever.Omega", "0.25");
  TranslateSingleInstanceDotSyntax(&fc);
  auto oc = fc.get<std::string>("Omega_cdm");
  assert(oc && *oc == "0.25");
}

static void test_dot_fluid_full_set() {
  FileContent fc;
  fc.set("de.type", "fluid");
  fc.set("de.Omega", "0.7");
  fc.set("de.w0", "-0.9");
  fc.set("de.wa", "0.1");
  fc.set("de.cs2", "1.0");
  fc.set("de.use_ppf", "no");
  fc.set("de.equation_of_state", "CLP");
  TranslateSingleInstanceDotSyntax(&fc);
  assert(*fc.get<std::string>("Omega_fld") == "0.7");
  assert(*fc.get<std::string>("w0_fld") == "-0.9");
  assert(*fc.get<std::string>("wa_fld") == "0.1");
  assert(*fc.get<std::string>("cs2_fld") == "1.0");
  assert(*fc.get<std::string>("use_ppf") == "no");
  assert(*fc.get<std::string>("fluid_equation_of_state") == "CLP");
}

static void test_dot_ur_forms() {
  FileContent fc;
  fc.set("nu.type", "ur");
  fc.set("nu.N", "3.044");
  TranslateSingleInstanceDotSyntax(&fc);
  assert(*fc.get<std::string>("N_ur") == "3.044");
}

static void test_dot_duplicate_instance_throws() {
  FileContent fc;
  fc.set("a.type", "baryons");
  fc.set("a.Omega", "0.05");
  fc.set("b.type", "baryons");
  fc.set("b.Omega", "0.04");
  bool threw = false;
  try { TranslateSingleInstanceDotSyntax(&fc); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

static void test_dot_conflict_with_legacy_throws() {
  FileContent fc;
  fc.set("Omega_b", "0.04");
  fc.set("x.type", "baryons");
  fc.set("x.Omega", "0.05");
  bool threw = false;
  try { TranslateSingleInstanceDotSyntax(&fc); }
  catch (const std::invalid_argument&) { threw = true; }
  assert(threw);
}

static void test_dot_identical_legacy_ok() {
  FileContent fc;
  fc.set("Omega_b", "0.05");
  fc.set("x.type", "baryons");
  fc.set("x.Omega", "0.05");
  TranslateSingleInstanceDotSyntax(&fc);  // must not throw
  assert(*fc.get<std::string>("Omega_b") == "0.05");
}

static void test_dot_unknown_field_left_unread() {
  FileContent fc;
  fc.set("x.type", "baryons");
  fc.set("x.Omeega", "0.05");  // typo: not in the table
  TranslateSingleInstanceDotSyntax(&fc);
  assert(!fc.was_read("x.Omeega"));  // stays unread -> warned later
}
```
Then add the calls to `main()` (before `return 0;`):
```cpp
  test_dot_baryons_translation();
  test_dot_name_is_ignored();
  test_dot_fluid_full_set();
  test_dot_ur_forms();
  test_dot_duplicate_instance_throws();
  test_dot_conflict_with_legacy_throws();
  test_dot_identical_legacy_ok();
  test_dot_unknown_field_left_unread();
```

- [ ] **Step 2: Verify it fails to compile (function undeclared)**

Run:
```bash
cmake --build build --target test-parser -j 2>&1 | tail -15
```
Expected: error — `TranslateSingleInstanceDotSyntax was not declared`.

- [ ] **Step 3: Declare the function**

In `species/species_input.h`, after the existing free-function declarations (after `SynthesiseIdenticalScalarField`):
```cpp
/**
 * Translate dot-syntax for SINGLE-INSTANCE legacy species into their legacy
 * keys, in place. For each recognised single-instance type T (photons, baryons,
 * cdm, ur, lambda, fluid), finds the at-most-one instance N with "N.type == T",
 * rewrites each known "N.<dot_field>" to its legacy key via FileContent::set,
 * and marks the dot entries read. The instance name N is discarded, so output
 * is identical to the legacy-key form. Throws std::invalid_argument if a type
 * appears more than once, or if a translated legacy key is already present with
 * a different value. Unknown "N.<field>" entries are left untouched (and unread,
 * so the usual unrecognised-parameter warning still fires).
 */
void TranslateSingleInstanceDotSyntax(FileContent* pfc);
```

- [ ] **Step 4: Implement the table + function**

In `species/species_input.cpp`, add at the top (after the existing includes):
```cpp
#include <string_view>
#include <vector>

#include "species/all_species.h"
```
Then add, inside the existing anonymous `namespace { … }` block:
```cpp
struct DotFieldMap {
  std::string_view dot;
  std::string_view legacy;
};

struct SingleInstanceSpec {
  std::string_view type;
  std::vector<DotFieldMap> fields;
};

// Module-owned: which named species are single-instance, and how each clean
// dot-field maps to its legacy key. Keyed by the species' canonical kTypeName.
// ScalarField is intentionally absent (reserved for a multi-instance future).
const std::vector<SingleInstanceSpec>& SingleInstanceTable() {
  static const std::vector<SingleInstanceSpec> table = {
      {PhotonsSpecies::kTypeName, {{"Omega", "Omega_g"}, {"T_cmb", "T_cmb"}}},
      {BaryonsSpecies::kTypeName, {{"Omega", "Omega_b"}}},
      {CDMSpecies::kTypeName, {{"Omega", "Omega_cdm"}}},
      {UltraRelativisticSpecies::kTypeName,
       {{"N", "N_ur"}, {"Omega", "Omega_ur"}, {"omega", "omega_ur"}}},
      {LambdaSpecies::kTypeName, {{"Omega", "Omega_Lambda"}}},
      {FluidSpecies::kTypeName,
       {{"Omega", "Omega_fld"},
        {"w0", "w0_fld"},
        {"wa", "wa_fld"},
        {"cs2", "cs2_fld"},
        {"Omega_EDE", "Omega_EDE"},
        {"c_gamma_over_c", "c_gamma_over_c_fld"},
        {"use_ppf", "use_ppf"},
        {"equation_of_state", "fluid_equation_of_state"}}},
  };
  return table;
}
```
Then add the function definition (after the anonymous namespace, with the other free functions):
```cpp
void TranslateSingleInstanceDotSyntax(FileContent* pfc) {
  if (!pfc) {
    throw std::invalid_argument("TranslateSingleInstanceDotSyntax: null FileContent*");
  }
  for (const auto& spec : SingleInstanceTable()) {
    const std::string type(spec.type);
    const auto instances = pfc->instances_with("type", type);
    if (instances.empty()) {
      continue;
    }
    if (instances.size() > 1) {
      std::string joined;
      for (size_t i = 0; i < instances.size(); ++i) {
        joined += (i ? ", " : "") + instances[i];
      }
      throw std::invalid_argument("species type '" + type +
                                  "' is single-instance but was given " +
                                  std::to_string(instances.size()) +
                                  " times (" + joined + ")");
    }
    SpeciesInput input(pfc, instances.front());
    (void) input.get<std::string>("type");  // consume "<name>.type"
    for (const auto& fm : spec.fields) {
      const std::string dot(fm.dot);
      auto value = input.get<std::string>(dot);  // consumes "<name>.<dot>" if present
      if (!value) {
        continue;
      }
      const std::string legacy(fm.legacy);
      auto existing = pfc->get<std::string>(legacy);
      if (existing && *existing != *value) {
        throw std::invalid_argument("input sets both legacy key '" + legacy +
                                    "' and dot-syntax '" + instances.front() + "." + dot +
                                    "' with different values");
      }
      pfc->set(legacy, *value);  // set() marks the key unread for the real consumer
    }
  }
}
```

- [ ] **Step 5: Build and pass the unit tests**

Run:
```bash
cmake --build build -j
ctest --test-dir build -R test-parser --output-on-failure
```
Expected: clean build; `test-parser` passes (all 8 new cases plus the originals).

- [ ] **Step 6: Commit**

```bash
git add species/species_input.h species/species_input.cpp tools/parser_test.cpp
git commit -m "feat: TranslateSingleInstanceDotSyntax dot->legacy preprocessing"
```

---

### Task 4: Part B wiring + byte-identical integration

Call the translation once before all input reads, then prove a dot-form `.ini` produces `cmp`-identical output to the legacy form.

**Files:**
- Modify: `source/input_module.cpp` (invoke the pass before `input_read_precisions()`)
- Create: `test/scenarios/dotsyn_legacy.ini`, `test/scenarios/dotsyn_dot.ini`, `test/scenarios/dotsyn_fluid_legacy.ini`, `test/scenarios/dotsyn_fluid_dot.ini`

**Interfaces:**
- Consumes: `TranslateSingleInstanceDotSyntax(FileContent*)` (Task 3).

- [ ] **Step 1: Find the invocation point**

Run:
```bash
grep -n 'input_read_precisions();' source/input_module.cpp
```
Expected: the orchestration call (around line 208), immediately before `ReadContext();` / `ConstructSpecies();`. The pass must run before this so translated keys feed precision, closure, and species reads.

- [ ] **Step 2: Invoke the pass before the first read**

In `source/input_module.cpp`, immediately before the `input_read_precisions();` call (line ~208), add:
```cpp
    // Translate dot-syntax for single-instance legacy species (bary.Omega ->
    // Omega_b) before any consumer reads the file.
    TranslateSingleInstanceDotSyntax(&file_content_);
    input_read_precisions();
```
(`species/all_species.h` is already included at `source/input_module.cpp:25`, which transitively provides `species/species_input.h`; if the symbol is not found, add `#include "../species/species_input.h"`.)

- [ ] **Step 3: Create the paired LCDM .ini files (baryons/cdm/photons)**

Create `test/scenarios/dotsyn_legacy.ini`:
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_b = 0.05
Omega_cdm = 0.25
T_cmb = 2.7255
```
Create `test/scenarios/dotsyn_dot.ini` (same physics, dot-form):
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
b.type = baryons
b.Omega = 0.05
c.type = cdm
c.Omega = 0.25
g.type = photons
g.T_cmb = 2.7255
```

- [ ] **Step 4: Create the paired fluid .ini files**

Create `test/scenarios/dotsyn_fluid_legacy.ini`:
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_Lambda = 0
fluid_equation_of_state = CLP
w0_fld = -0.9
wa_fld = 0.1
cs2_fld = 1.0
use_ppf = no
```
Create `test/scenarios/dotsyn_fluid_dot.ini`:
```ini
output = tCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_Lambda = 0
de.type = fluid
de.equation_of_state = CLP
de.w0 = -0.9
de.wa = 0.1
de.cs2 = 1.0
de.use_ppf = no
```

- [ ] **Step 5: Build, run both forms, and `cmp`**

Run:
```bash
cmake --build build -j
for pair in "dotsyn_legacy dotsyn_dot" "dotsyn_fluid_legacy dotsyn_fluid_dot"; do
  set -- $pair
  rm -f output/*; ./class test/scenarios/$1.ini
  mkdir -p /tmp/dotsyn/$1 && mv output/* /tmp/dotsyn/$1/
  rm -f output/*; ./class test/scenarios/$2.ini
  mkdir -p /tmp/dotsyn/$2 && mv output/* /tmp/dotsyn/$2/
  for f in /tmp/dotsyn/$1/*.dat; do
    cmp -s "$f" /tmp/dotsyn/$2/$(basename "$f") \
      && echo "OK: $(basename "$f")" \
      || echo "DIFFERS: $1 vs $2 / $(basename "$f")"
  done
done
```
Expected: only `OK:` lines, no `DIFFERS`. Dot-form output is byte-identical to legacy-form.

- [ ] **Step 6: Spot-check a negative path end-to-end**

Run:
```bash
printf 'output = mPk\nb.type = baryons\nb.Omega = 0.05\natom.type = baryons\natom.Omega = 0.04\n' > /tmp/dotsyn/dup.ini
./class /tmp/dotsyn/dup.ini; echo "exit=$?"
```
Expected: non-zero exit with an error mentioning single-instance type `baryons` given twice (not a silent run).

- [ ] **Step 7: Commit**

```bash
git add source/input_module.cpp test/scenarios/dotsyn_legacy.ini test/scenarios/dotsyn_dot.ini \
  test/scenarios/dotsyn_fluid_legacy.ini test/scenarios/dotsyn_fluid_dot.ini
git commit -m "feat: wire single-instance dot-syntax translation into input_module"
```

---

### Task 5: Final verification + finish the branch

**Files:** none modified (verification + PR).

- [ ] **Step 1: Full ctest + Part A regression re-check**

Run:
```bash
ctest --test-dir build --output-on-failure
for s in gauge_lcdm gauge_fluid gauge_ncdm gauge_dcdm dncdm_dr; do
  rm -f output/*; ./class test/scenarios/$s.ini
  mkdir -p /tmp/dotsyn/final/$s && mv output/* /tmp/dotsyn/final/$s/
  for f in /tmp/dotsyn/base/$s/*.dat; do
    cmp -s "$f" /tmp/dotsyn/final/$s/$(basename "$f") || echo "DIFFERS: $s/$(basename "$f")"
  done
done
```
Expected: all ctest green; no `DIFFERS` (the whole PR is output-neutral for legacy inputs).

- [ ] **Step 2: Grep audit — no stale CamelCase factory names or hardcoded type literals remain**

Run:
```bash
grep -rn '"NCDM"\|"NCDMGreyBody"\|"DNCDM_DR"\|"NCDMInt"\|"Photons"\|"Baryons"\|"Lambda"\|"Fluid"\|"ScalarField"' species/ source/ | grep -v kTypeName
grep -rn 'instances_with("type", "' species/ source/
```
Expected: first grep returns nothing meaningful (no old factory-name string literals); second returns nothing (all type literals now route through `kTypeName`).

- [ ] **Step 3: Python wrapper smoke test**

Run:
```bash
pip install . 2>&1 | tail -3
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py 2>&1 | tail -15; cd ..
```
Expected: install succeeds (regenerates `cclassy.pxd` from the modified headers); pytest passes. If wrapper generation chokes on a header change, fix `python/generate_wrapper.py` — do NOT hand-edit `cclassy.pxd`.

- [ ] **Step 4: Finish the branch**

Use superpowers:verification-before-completion, then superpowers:finishing-a-development-branch. PR title: `Dot-notation for single-instance legacy species (+ factory name/type unification)`. PR body must include: the `cmp`-identical evidence (Part A scenarios + dot-vs-legacy pairs), the negative-path check, ctest + pytest results, and a link to the spec.

---

## Self-Review

- **Spec coverage:** Part A (factory/type unification, `kTypeName`, ClosureSpeciesName, is_ncdm_family) → Task 2. Part B (translation function, table, algorithm, conflict + single-instance guards, unknown-field passthrough) → Task 3. Wiring + byte-identical verification → Task 4. UR inclusion (`N`/`Omega`/`omega`) → Task 3 table + `test_dot_ur_forms`. ScalarField exclusion → table comment + absence. Negative tests → Task 3 (unit) + Task 4 Step 6 (e2e). Verification model (cmp, ctest, pytest) → Tasks 1/4/5. All spec sections mapped.
- **Placeholders:** none — every code/edit step shows complete content and exact commands with expected output.
- **Type consistency:** `kTypeName` is `static constexpr std::string_view` everywhere; `TranslateSingleInstanceDotSyntax(FileContent*)` signature identical in header, definition, and all call sites; `DotFieldMap`/`SingleInstanceSpec`/`SingleInstanceTable()` used consistently; `instances_with` always passed `std::string(...)`; table keys (`Omega`,`w0`,`wa`,`cs2`,`Omega_EDE`,`c_gamma_over_c`,`use_ppf`,`equation_of_state`,`N`,`omega`,`T_cmb`) and legacy targets match the spec table.
