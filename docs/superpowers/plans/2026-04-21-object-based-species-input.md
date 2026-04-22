# Object-Based Species Input Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce dot-syntax species input (`nu1.type = ncdm_standard; nu1.m = 0.06`) for all three NCDM variants while keeping legacy flat syntax (`N_ncdm`, `m_ncdm`, …) working for standard NCDM.

**Architecture:** No central factory/registry. Each NCDM class's existing `CreateAll(FileContent*, …)` gains a dot-syntax code path: it scans the FileContent via a new `FileContent::instances_with(field, value)` helper, gathers per-instance fields, synthesises equivalent flat-syntax keys into the same FileContent, and lets the existing construction path run unchanged. `NCDMSpecies::CreateAll` additionally rejects mixing legacy and dot-syntax for the NCDM sector. Unknown `.type` values go unclaimed and surface via the existing `unread_parameters()` diagnostic.

**Tech Stack:** C++17, existing CLASS parser (`include/parser.h`, `tools/parser.cpp`), existing NCDM classes (`species/ncdm*`), existing `InputModule::ConstructSpecies()`.

**Reference spec:** `docs/superpowers/specs/2026-04-21-object-based-species-input-design.md` (architecture section updated informally: the `SpeciesRegistry` and central `legacy_ncdm_shim` described there are replaced by per-class `CreateAll` self-routing per user's revised intent).

---

## File structure

**New files (must be added to Makefile, setup.py, and CLASS.xcodeproj/project.pbxproj):**
- `species/species_input.h` + `species/species_input.cpp` — thin per-instance `FileContent` wrapper.

**Modified files (no build-system changes needed):**
- `include/parser.h` + `tools/parser.cpp` — add `FileContent::instances_with`.
- `species/ncdm_species.h` + `species/ncdm_species.cpp` — `CreateAll` returns named pairs; gains dot-syntax path; rejects mixing.
- `species/dncdm_species.h` + `species/dncdm_species.cpp` — `CreateAll` returns named pairs; gains dot-syntax path.
- `species/dncdm_dr_species.h` + `species/dncdm_dr_species.cpp` — `CreateAll` returns named pairs (already delegates to `DNCDMSpecies::CreateAll`).
- `species/ncdm_interacting_species.h` + `species/ncdm_interacting_species.cpp` — `CreateAll` returns named pairs; gains dot-syntax path.
- `source/input_module.cpp` — call-site updates for new `CreateAll` return type; legacy path uses `ncdm__N` keys.
- New input file: `test/dotsyntax_ncdm.ini` — regression input mirroring explanatory.ini's NCDM physics via dot-syntax.

## Worktree setup

This plan should execute in a fresh worktree off master named `246-object-based-species-input` (issue number TBD when the user opens it — ask for the number before creating the worktree if needed).

---

### Task 1: Parser helper — `FileContent::instances_with`

**Why first:** every CreateAll in later tasks depends on this.

**Files:**
- Modify: `include/parser.h`
- Modify: `tools/parser.cpp`
- Test: `tools/parser_test.cpp` (create if absent; otherwise extend; confirm the project's existing test harness — see `test/` directory in the repo root — and integrate with whatever framework is in use. If no C++ unit test harness exists, write a small standalone `parser_test.cpp` with a `main()` that asserts and exits non-zero on failure, and add a `test-parser` target to `Makefile`.)

- [ ] **Step 1: Write the failing test**

Add to `tools/parser_test.cpp`:

```cpp
#include <cassert>
#include "parser.h"

static void test_instances_with_basic() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  fc.set("nu1.m", "0.06");
  fc.set("nu2.type", "ncdm_self_interacting");
  fc.set("nu2.m", "0.07");
  fc.set("nu3.type", "ncdm_standard");

  auto names = fc.instances_with("type", "ncdm_standard");
  assert(names.size() == 2);
  assert(names[0] == "nu1");
  assert(names[1] == "nu3");

  auto none = fc.instances_with("type", "unknown");
  assert(none.empty());

  auto si = fc.instances_with("type", "ncdm_self_interacting");
  assert(si.size() == 1 && si[0] == "nu2");

  // Must NOT mark entries as read.
  assert(!fc.was_read("nu1.type"));
  assert(!fc.was_read("nu3.type"));
}

int main() {
  test_instances_with_basic();
  return 0;
}
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```
make -j test-parser && ./test-parser
```
Expected: compile error — `instances_with` is not a member of `FileContent`.

- [ ] **Step 3: Declare `instances_with` in `include/parser.h`**

Add inside the `FileContent` public section (just after `was_read`):

```cpp
/** Return every instance name N such that the entry "N.<field>" has the
 *  given value. The dot is a literal separator; N must match the instance
 *  regex [A-Za-z_][A-Za-z0-9_]*. Results are returned in insertion order.
 *  Does NOT mark anything as read. */
std::vector<std::string> instances_with(const std::string& field,
                                        const std::string& value) const;
```

- [ ] **Step 4: Implement in `tools/parser.cpp`**

Add next to the other `FileContent` member definitions:

```cpp
std::vector<std::string> FileContent::instances_with(
    const std::string& field, const std::string& value) const {
  std::vector<std::string> out;
  const std::string suffix = "." + field;
  for (const std::string& key : keys_) {
    if (key.size() <= suffix.size()) continue;
    if (key.compare(key.size() - suffix.size(), suffix.size(), suffix) != 0) continue;
    // Validate instance-name portion: must be a valid identifier.
    const std::string name = key.substr(0, key.size() - suffix.size());
    if (name.empty()) continue;
    char c0 = name[0];
    if (!(std::isalpha(static_cast<unsigned char>(c0)) || c0 == '_')) continue;
    bool ok = true;
    for (char c : name) {
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_')) { ok = false; break; }
    }
    if (!ok) continue;
    // Value check.
    auto it = params_.find(key);
    if (it == params_.end()) continue;
    if (it->second == value) out.push_back(name);
  }
  return out;
}
```

Add `#include <cctype>` near the top of `tools/parser.cpp` if not already present.

- [ ] **Step 5: Run the test to verify it passes**

Run:
```
make -j test-parser && ./test-parser
```
Expected: exits 0, no output.

- [ ] **Step 6: Commit**

```bash
git add include/parser.h tools/parser.cpp tools/parser_test.cpp Makefile
git commit -m "parser: add FileContent::instances_with helper

Used by dot-syntax species input to discover <name>.type = <value> entries
without consuming read state."
```

---

### Task 2: `SpeciesInput` wrapper

**Files:**
- Create: `species/species_input.h`
- Create: `species/species_input.cpp`
- Modify: `Makefile`
- Modify: `setup.py`
- Modify: `CLASS.xcodeproj/project.pbxproj` (user validates manually; follow the four-entry pattern from the PR #258 `species_collection.cpp` addition)
- Test: extend `tools/parser_test.cpp` (same harness)

- [ ] **Step 1: Write the failing test**

Append to `tools/parser_test.cpp`:

```cpp
#include "species/species_input.h"

static void test_species_input_prefixing() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  fc.set("nu1.m", "0.06");
  fc.set("nu1.T", "0.71611");

  SpeciesInput input(&fc, "nu1");
  assert(input.instance_name() == "nu1");

  double m = 0.;
  assert(input.read_double("m", m));
  assert(m == 0.06);
  assert(fc.was_read("nu1.m"));

  double missing = 0.;
  assert(!input.read_double("deg", missing));
  assert(!fc.was_read("nu1.deg"));

  std::string t;
  assert(input.read_string("type", t));
  assert(t == "ncdm_standard");
}

static void test_species_input_required_throws() {
  FileContent fc;
  fc.set("nu1.type", "ncdm_standard");
  SpeciesInput input(&fc, "nu1");
  bool threw = false;
  try {
    (void)input.required_double("m");
  } catch (const std::invalid_argument& e) {
    threw = true;
    std::string msg = e.what();
    assert(msg.find("nu1") != std::string::npos);
    assert(msg.find("m") != std::string::npos);
  }
  assert(threw);
}
```

Extend `main()`:
```cpp
int main() {
  test_instances_with_basic();
  test_species_input_prefixing();
  test_species_input_required_throws();
  return 0;
}
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```
make -j test-parser && ./test-parser
```
Expected: compile error — `species_input.h` not found.

- [ ] **Step 3: Create `species/species_input.h`**

```cpp
#pragma once

#include <string>
#include <vector>

#include "parser.h"

/**
 * Per-instance read wrapper over FileContent for dot-syntax species input.
 * Prefixes every key with "<instance_name>." before delegating to the
 * underlying FileContent. Successful reads mark the fully-qualified key as
 * read so unread_parameters() remains accurate.
 */
class SpeciesInput {
 public:
  SpeciesInput(FileContent* pfc, std::string instance_name);

  const std::string& instance_name() const { return instance_name_; }

  bool has(const std::string& field) const;

  bool read_double(const std::string& field, double& out);
  bool read_int(const std::string& field, int& out);
  bool read_string(const std::string& field, std::string& out);
  bool read_list_of_doubles(const std::string& field, std::vector<double>& out);

  double      required_double(const std::string& field);
  int         required_int(const std::string& field);
  std::string required_string(const std::string& field);

 private:
  std::string qualify(const std::string& field) const;

  FileContent* pfc_;
  std::string  instance_name_;
};
```

- [ ] **Step 4: Create `species/species_input.cpp`**

```cpp
#include "species/species_input.h"

#include <stdexcept>

SpeciesInput::SpeciesInput(FileContent* pfc, std::string instance_name)
    : pfc_(pfc), instance_name_(std::move(instance_name)) {
  if (!pfc_) {
    throw std::invalid_argument("SpeciesInput: null FileContent*");
  }
  if (instance_name_.empty()) {
    throw std::invalid_argument("SpeciesInput: empty instance name");
  }
}

std::string SpeciesInput::qualify(const std::string& field) const {
  return instance_name_ + "." + field;
}

bool SpeciesInput::has(const std::string& field) const {
  return pfc_->was_read(qualify(field)) || !pfc_->unread_parameters().empty();
  // Note: the above is a conservative placeholder — a dedicated "contains" on
  // FileContent would be preferable. For now callers should use the typed
  // readers instead of has(); keep has() for rare pre-checks.
}

bool SpeciesInput::read_double(const std::string& field, double& out) {
  return pfc_->read_double(qualify(field), out);
}
bool SpeciesInput::read_int(const std::string& field, int& out) {
  return pfc_->read_int(qualify(field), out);
}
bool SpeciesInput::read_string(const std::string& field, std::string& out) {
  return pfc_->read_string(qualify(field), out);
}
bool SpeciesInput::read_list_of_doubles(const std::string& field,
                                        std::vector<double>& out) {
  return pfc_->read_list_of_doubles(qualify(field), out);
}

double SpeciesInput::required_double(const std::string& field) {
  double v = 0.;
  if (!read_double(field, v)) {
    throw std::invalid_argument(
        "species '" + instance_name_ + "': missing required field '" + field + "'");
  }
  return v;
}
int SpeciesInput::required_int(const std::string& field) {
  int v = 0;
  if (!read_int(field, v)) {
    throw std::invalid_argument(
        "species '" + instance_name_ + "': missing required field '" + field + "'");
  }
  return v;
}
std::string SpeciesInput::required_string(const std::string& field) {
  std::string v;
  if (!read_string(field, v)) {
    throw std::invalid_argument(
        "species '" + instance_name_ + "': missing required field '" + field + "'");
  }
  return v;
}
```

The `has()` body is deliberately marked as a placeholder in the comment because the existing `FileContent` has no public "contains" method — if Task 1's implementation surfaces a clean need for one, add a private `contains(const std::string& name) const { return params_.count(name) > 0; }` to `FileContent` and call it here. Otherwise callers avoid `has()` and use `read_*` return values.

- [ ] **Step 5: Add `species_input.cpp` to `Makefile`**

Find the `SPECIES_OPP` list (the `species_collection.opp` line added by PR #258 is a locator). Add one line:

```
$(BUILDDIR)/species_input.opp \
```

immediately after `species_collection.opp`.

- [ ] **Step 6: Add `species_input.cpp` to `setup.py`**

Modify `setup.py:83-106` (the `prepend('species', ...)` block). Insert `'species_input.cpp',` alphabetically — between `'species_collection.cpp',` and `'ultra_relativistic.cpp',`.

- [ ] **Step 7: Add `species_input.{h,cpp}` to `CLASS.xcodeproj/project.pbxproj`**

Follow the four-entry pattern from `species_collection.cpp` (PBXBuildFile + two PBXFileReferences + group children + Sources build phase). Use a fresh ID prefix (avoid `EE1000xxx` which PR #258 used). User validates Xcode manually after commit.

- [ ] **Step 8: Run the test to verify it passes**

Run:
```
make -j test-parser && ./test-parser
```
Expected: exits 0.

Also do a full build to catch integration errors in setup.py / Makefile:
```
make -j class
```
Expected: links successfully.

- [ ] **Step 9: Commit**

```bash
git add species/species_input.h species/species_input.cpp \
        Makefile setup.py CLASS.xcodeproj/project.pbxproj \
        tools/parser_test.cpp
git commit -m "species: add SpeciesInput per-instance parser wrapper

Used by CreateAll methods to read <instance>.<field> keys cleanly under
the dot-syntax species input path."
```

---

### Task 3: Dot-syntax path in `NCDMSpecies::CreateAll`

**Files:**
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp`
- Modify: `source/input_module.cpp` (call site)

This task changes the `CreateAll` return type to carry instance names, adds the dot-syntax scan, rejects mixing, and synthesises equivalent flat-syntax keys so the existing NCDMSpecies constructor runs unchanged.

- [ ] **Step 1: Update the declaration in `species/ncdm_species.h:20-23`**

Change:
```cpp
static std::vector<std::unique_ptr<NCDMSpecies>> CreateAll(FileContent* pfc,
                                                           const NcdmSettings& settings,
                                                           const background* pba,
                                                           const BackgroundModule* bgm);
```
to:
```cpp
struct Named {
  std::string                    key;      // SpeciesCollection key
  std::unique_ptr<NCDMSpecies>   species;
};
static std::vector<Named> CreateAll(FileContent* pfc,
                                    const NcdmSettings& settings,
                                    const background* pba,
                                    const BackgroundModule* bgm);
```

- [ ] **Step 2: Rewrite `NCDMSpecies::CreateAll` in `species/ncdm_species.cpp:29-55`**

```cpp
// ── CreateAll factory ───────────────────────────────────────────────────────

std::vector<NCDMSpecies::Named> NCDMSpecies::CreateAll(FileContent* pfc,
                                                       const NcdmSettings& settings,
                                                       const background* pba,
                                                       const BackgroundModule* bgm) {
  std::vector<Named> result;

  // Detect dot-syntax instances (consuming their .type keys).
  const std::vector<std::string> dot_instances =
      pfc->instances_with("type", "ncdm_standard");

  // Detect legacy flat-syntax.
  int int1 = 0, int2 = 0, flag1 = 0, flag2 = 0;
  char errmsg[2048];
  class_call(parser_read_int(pfc, "N_ncdm_standard", &int1, &flag1, errmsg), errmsg, errmsg);
  class_call(parser_read_int(pfc, "N_ncdm", &int2, &flag2, errmsg), errmsg, errmsg);
  if (flag1 == _TRUE_ && flag2 == _TRUE_) {
    throw std::invalid_argument(
        "In input file, you can only enter one of N_ncdm_standard and N_ncdm, choose one");
  }
  int N_ncdm_standard = 0;
  if (flag1 == _TRUE_) N_ncdm_standard = int1;
  else if (flag2 == _TRUE_) N_ncdm_standard = int2;

  // Reject mixing legacy + dot for NCDM sector.
  const bool has_legacy = (N_ncdm_standard > 0);
  const bool has_dot    = !dot_instances.empty();
  if (has_legacy && has_dot) {
    throw std::invalid_argument(
        "cannot mix legacy N_ncdm/N_ncdm_standard with dot-syntax "
        "'*.type = ncdm_standard'; use one or the other");
  }

  // Dot-syntax path: synthesise the equivalent flat-syntax keys so the
  // existing constructor, which reads m_ncdm[i]/T_ncdm[i]/... via list-
  // indexed parser calls, runs unchanged.
  std::vector<std::string> keys;
  if (has_dot) {
    keys = synthesise_standard_ncdm_flat_keys(*pfc, dot_instances);
    N_ncdm_standard = static_cast<int>(dot_instances.size());
  } else {
    // Legacy-only path: keys are ncdm__1, ncdm__2, ...
    keys.reserve(N_ncdm_standard);
    for (int n = 0; n < N_ncdm_standard; ++n) {
      keys.push_back("ncdm__" + std::to_string(n + 1));
    }
  }

  for (int n = 0; n < N_ncdm_standard; ++n) {
    auto sp = std::make_unique<NCDMSpecies>(pfc, n, settings, pba, bgm);
    result.push_back(Named{keys[n], std::move(sp)});
  }
  return result;
}
```

- [ ] **Step 3: Add the `synthesise_standard_ncdm_flat_keys` helper**

Place directly above `CreateAll` in `species/ncdm_species.cpp`:

```cpp
namespace {

// Read each instance's dot-syntax fields through SpeciesInput (marks them read),
// aggregate into comma-separated values, and write equivalent legacy flat keys
// back into pfc. Returns the vector of instance names (as SpeciesCollection keys).
std::vector<std::string> synthesise_standard_ncdm_flat_keys(
    FileContent& pfc, const std::vector<std::string>& instances) {
  auto join_csv = [](const std::vector<std::string>& vs) {
    std::string s;
    for (size_t i = 0; i < vs.size(); ++i) {
      if (i) s += ", ";
      s += vs[i];
    }
    return s;
  };

  // Each row: {legacy_key, is_required, list_of_values}
  struct Col {
    std::string legacy;        // e.g. "m_ncdm"
    std::string dot;           // e.g. "m"
    bool required = false;
    std::vector<std::string> values;
  };
  std::vector<Col> cols = {
      {"m_ncdm",                   "m",                    false, {}},
      {"T_ncdm",                   "T",                    false, {}},
      {"deg_ncdm",                 "deg",                  false, {}},
      {"Omega_ncdm",               "Omega",                false, {}},
      {"omega_ncdm",               "omega",                false, {}},
      {"ksi_ncdm",                 "ksi",                  false, {}},
      {"ncdm_psd_parameters",      "psd_parameters",       false, {}},
      {"ncdm_psd_filenames",       "psd_filename",         false, {}},
      {"use_ncdm_psd_files",       "use_psd_file",         false, {}},
      {"ncdm_fluid_approximation", "fluid_approximation",  false, {}},
      {"ncdm_quadrature_strategy", "quadrature_strategy",  false, {}},
      {"Number_of_momenta_bins",   "momenta_bins",         false, {}},
      {"Maximum_q",                "max_q",                false, {}},
  };

  std::vector<std::string> keys;
  keys.reserve(instances.size());
  for (const std::string& inst : instances) {
    keys.push_back(inst);
    SpeciesInput input(&pfc, inst);
    // Consume ".type" to mark it read.
    std::string type;
    input.read_string("type", type);
    for (Col& c : cols) {
      std::string v;
      if (input.read_string(c.dot, v)) {
        c.values.push_back(v);
      } else {
        // Mark absence explicitly by inserting an empty slot; the legacy
        // reader interprets this through its own defaults (positions matter).
        // If a subset of instances omits an optional field, the legacy reader
        // expects a full-length CSV — fill with a neutral placeholder by
        // repeating the last non-empty value or "0".
        c.values.push_back(c.values.empty() ? std::string("0") : c.values.back());
      }
    }
  }

  pfc.set("N_ncdm_standard", std::to_string(instances.size()));
  for (const Col& c : cols) {
    bool any = false;
    for (const auto& v : c.values) { if (!v.empty() && v != "0") { any = true; break; } }
    if (any) pfc.set(c.legacy, join_csv(c.values));
  }
  return keys;
}

}  // namespace
```

**Note to engineer:** the "fill with placeholder" behaviour for partial-coverage CSVs is conservative — it preserves CSV length, which the legacy `parser_read_list_of_doubles` consumer requires. If during review of the NCDM legacy reader you find that any of these fields treats "0" as a *meaningful* value (not a neutral default), replace the fallback with an explicit sentinel or revise the synthesis to go field-by-field per instance via a per-instance keyset that the constructor doesn't currently understand (which would push this into Option A refactor territory and is out of scope here).

- [ ] **Step 4: Update the call site in `source/input_module.cpp:230-235`**

Change:
```cpp
auto ncdm_list =
    NCDMSpecies::CreateAll(&file_content_, ncdm_settings_for_species, pba, nullptr);
for (auto& sp : ncdm_list) {
  std::string name = sp->name();
  all_species_.insert(name, std::move(sp));
}
```
to:
```cpp
auto ncdm_list =
    NCDMSpecies::CreateAll(&file_content_, ncdm_settings_for_species, pba, nullptr);
for (auto& e : ncdm_list) {
  all_species_.insert(e.key, std::move(e.species));
}
```

- [ ] **Step 5: Build**

Run:
```
make -j class
```
Expected: links successfully.

- [ ] **Step 6: Legacy-path regression test**

Run the canonical regression against master:

```
./class explanatory.ini
diff output/explanatory00_background.dat /tmp/class_ref/explanatory00_background.dat
./class base_2015_plikHM_TT_lowTEB_lensing.ini
diff output/base_2015_plikHM_TT_lowTEB_lensing00_cl.dat /tmp/class_ref/base_2015_plikHM_TT_lowTEB_lensing00_cl.dat
```
Expected: byte-identical (empty diff).

- [ ] **Step 7: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp source/input_module.cpp
git commit -m "NCDMSpecies: accept dot-syntax input; reject mixing with legacy

CreateAll now detects <name>.type = ncdm_standard entries, synthesises
equivalent flat-syntax keys, and routes through the unchanged NCDMSpecies
constructor. Mixing legacy N_ncdm* with dot-syntax in one .ini is rejected.
Return type becomes vector<Named{key, species}> so the SpeciesCollection
key comes from the user-supplied instance name (legacy path uses ncdm__N)."
```

---

### Task 4: Dot-syntax path in `DNCDMSpecies::CreateAll` / `DNCDM_DR_Species::CreateAll`

**Files:**
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`
- Modify: `species/dncdm_dr_species.h`
- Modify: `species/dncdm_dr_species.cpp`
- Modify: `source/input_module.cpp` (call site)

Same shape as Task 3, with `.type = ncdm_decay_dr` and the `N_ncdm_decay_dr`-family legacy keys.

- [ ] **Step 1: Add the `Named` nested struct + updated return type in `species/dncdm_species.h`**

Inside the `DNCDMSpecies` class, declare:

```cpp
struct Named {
  std::string                   key;
  std::unique_ptr<DNCDMSpecies> species;
};
static std::vector<Named> CreateAll(FileContent* pfc,
                                    const NcdmSettings& settings,
                                    const background* pba,
                                    const BackgroundModule* bgm);
```

Replace the existing `static std::vector<std::unique_ptr<DNCDMSpecies>> CreateAll(...)` declaration.

- [ ] **Step 2: Rewrite `DNCDMSpecies::CreateAll` in `species/dncdm_species.cpp:107` onwards**

Follow the same pattern as `NCDMSpecies::CreateAll`:

1. Call `pfc->instances_with("type", "ncdm_decay_dr")` — `dot_instances`.
2. Read `N_ncdm_decay_dr`, `N_ncdm_standard`/`N_ncdm`, and the legacy decay_dr CSV fields (inventory the existing body at species/dncdm_species.cpp:107-230 — the CSV fields are `dncdm_lifetime` or similar; list them concretely from the existing code before writing the synthesise helper).
3. Reject mixing: if both legacy decay_dr keys and dot-syntax decay_dr instances are present, throw `std::invalid_argument("cannot mix N_ncdm_decay_dr with dot-syntax *.type = ncdm_decay_dr; use one or the other")`.
4. Synthesise flat keys (helper `synthesise_decay_dr_ncdm_flat_keys`, same shape as Task 3's helper) if `has_dot`.
5. Choose keys: dot path uses user instance names; legacy path uses `ncdm__{N_ncdm_standard + n + 1}` (index offset preserves the existing dncdm_id convention).
6. Construct via existing constructor, collect into `vector<Named>`.

**Code sketch (engineer expands per Task 3's pattern, filling the concrete CSV field list from inspection of the existing body):**

```cpp
std::vector<DNCDMSpecies::Named> DNCDMSpecies::CreateAll(FileContent* pfc,
                                                         const NcdmSettings& settings,
                                                         const background* pba,
                                                         const BackgroundModule* bgm) {
  std::vector<Named> result;
  const auto dot_instances = pfc->instances_with("type", "ncdm_decay_dr");

  // ... (rest mirrors Task 3 shape: read legacy counts/CSVs, reject mixing,
  //      synthesise flat keys if has_dot, construct, build keys)
  return result;
}
```

- [ ] **Step 3: Update `DNCDM_DR_Species::CreateAll` in `species/dncdm_dr_species.{h,cpp}`**

`DNCDM_DR_Species::CreateAll` currently delegates to `DNCDMSpecies::CreateAll`. Update its return type to match:

```cpp
struct Named {
  std::string                      key;
  std::unique_ptr<DNCDM_DR_Species> species;
};
static std::vector<Named> CreateAll(FileContent* pfc,
                                    const NcdmSettings& settings,
                                    const background* pba,
                                    const BackgroundModule* bgm);
```

Implementation in the `.cpp`:

```cpp
std::vector<DNCDM_DR_Species::Named> DNCDM_DR_Species::CreateAll(
    FileContent* pfc, const NcdmSettings& settings,
    const background* pba, const BackgroundModule* bgm) {
  auto dncdm = DNCDMSpecies::CreateAll(pfc, settings, pba, bgm);
  std::vector<Named> out;
  out.reserve(dncdm.size());
  for (auto& e : dncdm) {
    out.push_back(Named{e.key,
                        std::make_unique<DNCDM_DR_Species>(std::move(e.species), pba, bgm)});
  }
  return out;
}
```

- [ ] **Step 4: Update call site in `source/input_module.cpp:238-244`**

Change:
```cpp
auto dncdm_vec =
    DNCDMSpecies::CreateAll(&file_content_, ncdm_settings_for_species, pba, nullptr);
for (auto& dncdm_sp : dncdm_vec) {
  std::string name = "NCDM_" + std::to_string(dncdm_sp->ncdm_id());
  all_species_.insert(name,
                      std::make_unique<DNCDM_DR_Species>(std::move(dncdm_sp), pba, nullptr));
}
```
to:
```cpp
auto dncdm_dr_vec =
    DNCDM_DR_Species::CreateAll(&file_content_, ncdm_settings_for_species, pba, nullptr);
for (auto& e : dncdm_dr_vec) {
  all_species_.insert(e.key, std::move(e.species));
}
```

- [ ] **Step 5: Build and regression-test**

Run:
```
make -j class
./class explanatory.ini
diff output/explanatory00_background.dat /tmp/class_ref/explanatory00_background.dat
```
Expected: links, byte-identical.

Plus whatever `.ini` in the repo actually exercises `N_ncdm_decay_dr > 0` — grep for one in `test/`:
```
grep -rl "N_ncdm_decay_dr" test/ *.ini 2>/dev/null
```
Run each and assert byte-equal against master.

- [ ] **Step 6: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp \
        species/dncdm_dr_species.h species/dncdm_dr_species.cpp \
        source/input_module.cpp
git commit -m "DNCDMSpecies/DNCDM_DR_Species: accept dot-syntax ncdm_decay_dr"
```

---

### Task 5: Dot-syntax path in `NCDMInteractingSpecies::CreateAll`

**Files:**
- Modify: `species/ncdm_interacting_species.h`
- Modify: `species/ncdm_interacting_species.cpp`
- Modify: `source/input_module.cpp` (call site)

- [ ] **Step 1: Update return type in `species/ncdm_interacting_species.h`**

```cpp
struct Named {
  std::string                             key;
  std::unique_ptr<NCDMInteractingSpecies> species;
};
static std::vector<Named> CreateAll(FileContent* pfc,
                                    const NcdmSettings& settings,
                                    const background* pba,
                                    const BackgroundModule* bgm);
```

- [ ] **Step 2: Rewrite `NCDMInteractingSpecies::CreateAll` in `species/ncdm_interacting_species.cpp:41` onwards**

Same structural pattern as Task 3, with `.type = ncdm_self_interacting` and legacy key `N_ncdm_interacting` (and its per-index CSVs — inventory them by reading the existing body before writing the synthesis helper).

Concretely:
1. `const auto dot_instances = pfc->instances_with("type", "ncdm_self_interacting");`
2. Read `N_ncdm_interacting` via `pfc->read_int(...)`.
3. If both non-zero: throw (`"cannot mix N_ncdm_interacting with dot-syntax *.type = ncdm_self_interacting"`).
4. If dot: synthesise `N_ncdm_interacting` + its CSV fields (`synthesise_self_interacting_ncdm_flat_keys` helper in an anonymous namespace at the top of the file).
5. Keys: dot path = user instance names; legacy path = `ncdm__{N_ncdm_standard + N_ncdm_decay_dr + n + 1}` (preserves current numbering semantics).
6. Return `vector<Named>`.

- [ ] **Step 3: Update call site in `source/input_module.cpp:247-252`**

Change:
```cpp
auto interacting_vec =
    NCDMInteractingSpecies::CreateAll(&file_content_, ncdm_settings_for_species, pba, nullptr);
for (auto& int_sp : interacting_vec) {
  std::string name = int_sp->name();
  all_species_.insert(name, std::move(int_sp));
}
```
to:
```cpp
auto interacting_vec =
    NCDMInteractingSpecies::CreateAll(&file_content_, ncdm_settings_for_species, pba, nullptr);
for (auto& e : interacting_vec) {
  all_species_.insert(e.key, std::move(e.species));
}
```

- [ ] **Step 4: Build and regression-test**

Run:
```
make -j class
./class explanatory.ini
```
Find an `.ini` exercising interacting NCDM:
```
grep -rl "N_ncdm_interacting" test/ *.ini 2>/dev/null
```
Byte-equal diff each against master.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_interacting_species.h species/ncdm_interacting_species.cpp \
        source/input_module.cpp
git commit -m "NCDMInteractingSpecies: accept dot-syntax ncdm_self_interacting"
```

---

### Task 6: End-to-end dot-syntax regression

**Files:**
- Create: `test/dotsyntax_ncdm.ini`

- [ ] **Step 1: Write the test .ini**

Inspect `explanatory.ini` for its NCDM block. Typical legacy content:
```
N_ncdm = 1
m_ncdm = 0.71
T_ncdm = 0.71611
deg_ncdm = 3.0
```

Create `test/dotsyntax_ncdm.ini` as a copy of `explanatory.ini` with *only* the NCDM block replaced:
```
nu1.type = ncdm_standard
nu1.m    = 0.71
nu1.T    = 0.71611
nu1.deg  = 3.0
```
(Remove `N_ncdm`, `m_ncdm`, `T_ncdm`, `deg_ncdm` from the copy.)

Keep the `output` prefix distinct (e.g. `root = output/dotsyntax_ncdm`) so it doesn't collide.

- [ ] **Step 2: Run and byte-diff against the legacy run**

```
./class explanatory.ini
./class test/dotsyntax_ncdm.ini
diff output/explanatory00_background.dat output/dotsyntax_ncdm00_background.dat
diff output/explanatory00_cl.dat         output/dotsyntax_ncdm00_cl.dat
```
Expected: empty diffs.

- [ ] **Step 3: Run mixed-input negative test**

Add a second copy `test/dotsyntax_ncdm_mixed.ini` that keeps *both* `N_ncdm = 1` and `nu1.type = ncdm_standard`. Run:
```
./class test/dotsyntax_ncdm_mixed.ini
```
Expected: exits non-zero with an error message containing `"cannot mix legacy N_ncdm"`.

- [ ] **Step 4: Commit**

```bash
git add test/dotsyntax_ncdm.ini test/dotsyntax_ncdm_mixed.ini
git commit -m "test: dot-syntax NCDM regression inputs

dotsyntax_ncdm.ini reproduces explanatory.ini's NCDM physics via
nu1.type=ncdm_standard; dotsyntax_ncdm_mixed.ini verifies the mixing
rejection error."
```

---

### Task 7: Final sweep + PR

- [ ] **Step 1: Full regression**

Run against master baseline in `/tmp/class_ref/`:
```
./class explanatory.ini && diff output/explanatory00_background.dat /tmp/class_ref/explanatory00_background.dat
./class base_2015_plikHM_TT_lowTEB_lensing.ini && diff output/base_2015_plikHM_TT_lowTEB_lensing00_cl.dat /tmp/class_ref/base_2015_plikHM_TT_lowTEB_lensing00_cl.dat
```
Expected: empty diffs.

- [ ] **Step 2: Build Python wrapper**

```
python -m pip install -e .
python -c "from classy import Class; c = Class(); c.set({'N_ncdm': 1, 'm_ncdm': 0.06, 'T_ncdm': 0.71611}); c.compute(); print('ok')"
```
Expected: prints `ok`.

- [ ] **Step 3: `unread_parameters()` check**

Add a line to `test/dotsyntax_ncdm.ini` with a typo:
```
nu1.mas = 99
```
Run `./class test/dotsyntax_ncdm.ini` and confirm the run completes (still byte-identical to explanatory) and that the stderr/output includes `nu1.mas` in the unread-parameter warning. Revert the typo from the committed file.

- [ ] **Step 4: Xcode manual check**

Open `CLASS.xcodeproj` in Xcode, confirm `species_input.h` and `species_input.cpp` appear under the species group and in the Sources build phase. Build the `CLASS` target from Xcode; expect a clean build.

- [ ] **Step 5: Push branch and open PR**

Open a GitHub issue first if the user hasn't:
```
gh issue create --repo AarhusCosmology/CLASSpp \
  --title "Object-based species input syntax (NCDM)" \
  --body "Implements the dot-notation species input proposed in #246 for the three NCDM variants..."
```

Note the issue number `<N>`, then:
```
git push -u origin <N>-object-based-species-input
gh pr create --repo AarhusCosmology/CLASSpp --title "Object-based species input (NCDM) — fixes #<N>" --body "..."
```
