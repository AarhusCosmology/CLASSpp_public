# FileContent `get<T>` Accessor API (#314) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the C-era `parser_read_*` out-param functions and the `class_read_*` macros with a single typed `std::optional<T> get<T>` accessor on `FileContent` (and a matching `get<T>`/`require<T>` on `SpeciesInput`), then delete the old surfaces.

**Architecture:** `FileContent::get<T>` becomes the single source of truth for typed reads (map lookup + parse + read-tracking). All other accessors are migrated to it, then deleted. The migration is staged so every intermediate build is green: add `get<T>` first (old `read_*` reimplemented on top of it), migrate all callers, reshape `SpeciesInput`, and only then delete `read_*` and `parser_read_int/double/string`. The refactor changes *how* values are read, not *which* — so a correct result is **byte-identical in output**, which is the regression witness.

**Tech Stack:** C++17 (`std::optional`, variadic templates, `if constexpr`/explicit specialization), CMake (in-source build at repo root; targets `class`, `test-parser`), `<cassert>`-based unit tests in `tools/parser_test.cpp`, pytest regression in `python/test_class.py`.

**Spec:** `docs/superpowers/specs/2026-06-18-issue-314-filecontent-accessors-design.md`

---

## File map

| File | Change |
|------|--------|
| `include/parser.h` | Add `#include <optional>`; add `template get<T>` decl + inline `get_or` (+ `const char*` overload); later delete `read_*` decls and `parser_read_int/double/string` decls |
| `tools/parser.cpp` | Add explicit `get<T>` specializations (own the parse logic); reimplement `read_*` on `get<T>` (transitional); later delete `read_*` + `parser_read_int/double/string` defs |
| `tools/parser_test.cpp` | Add `get<T>` tests; migrate `SpeciesInput` tests to `get`/`require` |
| `species/species_input.h` / `.cpp` | Replace `read_*`/`required_*` with `template get<T>`/`require<T>` |
| `source/input_module.h` | Delete all `class_read_*` + mutual-exclusion macros |
| `source/input_module.cpp` | Add `n_present`/`read_one_of` helpers; migrate 256 macro/raw sites + 175 precision-overload bodies + `readDoubleList`; delete shared `param/int/string/flag` locals |
| 9 species `.cpp` files | Migrate ~55 direct `ctx.pfc->read_*` callers + 35 `input.read_*` callers to `get<T>` |

## Conventions used in every task

- **Build (C++):** `cmake --build . -j --target test-parser class` from repo root (in-source build; `CMakeCache.txt` is in the root).
- **Run C++ unit tests:** `./test-parser && echo OK` (exit 0 = pass; asserts abort on failure). `ctest --output-on-failure` runs all three (`test-parser`, `test-bisection`, `test-photons`).
- **Smoke:** `./class explanatory.ini` should run to completion (fast sanity that input parsing still works end-to-end).
- **Commit** after each task with the message shown.

---

### Task 1: Add `FileContent::get<T>` + `get_or` as the new single source of truth

**Files:**
- Modify: `include/parser.h`
- Modify: `tools/parser.cpp`
- Test: `tools/parser_test.cpp`

- [ ] **Step 1: Write the failing tests** — append to `tools/parser_test.cpp` and call them from `main()`:

```cpp
#include <optional>
#include <vector>

static void test_get_typed() {
  FileContent fc;
  fc.set("an_int", "42");
  fc.set("a_double", "0.067");
  fc.set("a_string", "hello");
  fc.set("a_list", "1.0,2.5,3");

  assert(fc.get<int>("an_int") == 42);
  assert(fc.get<double>("a_double") == 0.067);
  assert(fc.get<std::string>("a_string").value() == "hello");
  auto l = fc.get<std::vector<double>>("a_list");
  assert(l && l->size() == 3 && (*l)[1] == 2.5);

  // Absent -> nullopt, and absence is not marked read.
  assert(!fc.get<int>("missing").has_value());
  assert(!fc.was_read("missing"));

  // A successful get marks the key read (parity with the old read_*).
  assert(fc.was_read("an_int"));

  // get_or: present returns value, absent returns fallback.
  assert(fc.get_or("a_double", 1.0) == 0.067);
  assert(fc.get_or("missing", 1.0) == 1.0);
  assert(fc.get_or("a_string", "fallback") == "hello");   // const char* overload
  assert(fc.get_or("missing_str", "fallback") == "fallback");
}

static void test_get_parse_error_throws() {
  FileContent fc;
  fc.set("bad_int", "not_a_number");
  bool threw = false;
  try { (void) fc.get<int>("bad_int"); }
  catch (const std::exception&) { threw = true; }
  assert(threw);
}
```
Add `test_get_typed();` and `test_get_parse_error_throws();` to `main()`.

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build . -j --target test-parser`
Expected: FAIL — compile error, `get` is not a member of `FileContent`.

- [ ] **Step 3: Declare the API in `include/parser.h`**

Add `#include <optional>` near the top. In the `public:` section of `FileContent` (next to the current `read_*` declarations), add:

```cpp
  /** Typed accessor: the value if @p name is present (marking it read), or
   *  std::nullopt if absent. Throws on a present-but-unparseable value.
   *  Supported T: int, double, std::string, std::vector<{int,double,std::string}>. */
  template <class T> std::optional<T> get(const std::string& name) const;

  /** Read-with-default: get<T>(name) if present, else @p fallback. */
  template <class T> T get_or(const std::string& name, T fallback) const {
    if (auto v = get<T>(name)) return *v;
    return fallback;
  }
  /** Overload so a string-literal fallback does not deduce T = const char*. */
  std::string get_or(const std::string& name, const char* fallback) const {
    if (auto v = get<std::string>(name)) return *v;
    return std::string(fallback);
  }
```

- [ ] **Step 4: Implement the specializations in `tools/parser.cpp`** (move the parse + read-tracking logic out of the current `read_*` bodies into `get<T>`, then make `read_*` thin wrappers).

Add these six explicit specializations (the parse/validation/`read_params_.insert` bodies are lifted verbatim from the current `read_*` definitions at `tools/parser.cpp:120-195`):

```cpp
template <>
std::optional<int> FileContent::get<int>(const std::string& name) const {
  auto it = params_.find(name);
  if (it == params_.end()) return std::nullopt;
  int value;
  if (std::sscanf(it->second.c_str(), "%d", &value) != 1)
    throw std::invalid_argument("Cannot read integer value of parameter '" + name + "' in file '" + filename_ + "'");
  read_params_.insert(name);
  return value;
}

template <>
std::optional<double> FileContent::get<double>(const std::string& name) const {
  auto it = params_.find(name);
  if (it == params_.end()) return std::nullopt;
  double value;
  if (std::sscanf(it->second.c_str(), "%lg", &value) != 1)
    throw std::invalid_argument("Cannot read double value of parameter '" + name + "' in file '" + filename_ + "'");
  read_params_.insert(name);
  return value;
}

template <>
std::optional<std::string> FileContent::get<std::string>(const std::string& name) const {
  auto it = params_.find(name);
  if (it == params_.end()) return std::nullopt;
  read_params_.insert(name);
  return it->second;
}

template <>
std::optional<std::vector<double>> FileContent::get<std::vector<double>>(const std::string& name) const {
  auto it = params_.find(name);
  if (it == params_.end()) return std::nullopt;
  const auto parts = split_csv(it->second);
  std::vector<double> values(parts.size());
  for (std::size_t i = 0; i < parts.size(); ++i)
    if (std::sscanf(parts[i].c_str(), "%lg", &values[i]) != 1)
      throw std::invalid_argument("Cannot read double entry " + std::to_string(i + 1) + " of parameter '" + name + "' in file '" + filename_ + "'");
  read_params_.insert(name);
  return values;
}

template <>
std::optional<std::vector<int>> FileContent::get<std::vector<int>>(const std::string& name) const {
  auto it = params_.find(name);
  if (it == params_.end()) return std::nullopt;
  const auto parts = split_csv(it->second);
  std::vector<int> values(parts.size());
  for (std::size_t i = 0; i < parts.size(); ++i)
    if (std::sscanf(parts[i].c_str(), "%d", &values[i]) != 1)
      throw std::invalid_argument("Cannot read integer entry " + std::to_string(i + 1) + " of parameter '" + name + "' in file '" + filename_ + "'");
  read_params_.insert(name);
  return values;
}

template <>
std::optional<std::vector<std::string>> FileContent::get<std::vector<std::string>>(const std::string& name) const {
  auto it = params_.find(name);
  if (it == params_.end()) return std::nullopt;
  read_params_.insert(name);
  return split_csv(it->second);
}
```

Then reimplement the existing members as one-line wrappers (kept until Task 5 so current callers still build):

```cpp
bool FileContent::read_int(const std::string& name, int& value) const {
  if (auto v = get<int>(name)) { value = *v; return true; }
  return false;
}
// same shape for read_double / read_string / read_list_of_doubles / read_list_of_integers / read_list_of_strings
```

Leave the primary `template <class T> std::optional<T> FileContent::get(...)` **undefined** (only the specializations are defined) so an unsupported `T` is a link error.

- [ ] **Step 5: Build + run tests**

Run: `cmake --build . -j --target test-parser && ./test-parser && echo OK`
Expected: builds, prints `OK` (all asserts pass, including the pre-existing `SpeciesInput` tests which still use `read_*`).

- [ ] **Step 6: Commit**

```bash
git add include/parser.h tools/parser.cpp tools/parser_test.cpp
git commit -m "#314: add FileContent::get<T>/get_or; read_* now wrap it"
```

---

### Task 2: Migrate `input_module.cpp` macro/raw reads to `get<T>`; delete the macros

**Files:**
- Modify: `source/input_module.h` (delete macros)
- Modify: `source/input_module.cpp` (helpers + call-site migration)

- [ ] **Step 1: Add the two policy helpers** to the anonymous namespace at the top of `source/input_module.cpp` (near `readDoubleList`):

```cpp
// Number of the given optionals that hold a value.
template <class... Opt>
int n_present(const Opt&... opts) { return (0 + ... + (opts.has_value() ? 1 : 0)); }

// Value of whichever of two names is present; error if BOTH present; nullopt if neither.
template <class T>
std::optional<T> read_one_of(const FileContent& fc, const char* n1, const char* n2) {
  auto a = fc.get<T>(n1);
  auto b = fc.get<T>(n2);
  class_test(a && b, "In input file, you can only enter one of %s, %s, choose one", n1, n2);
  return a ? a : b;
}
```

- [ ] **Step 2: Migrate the simple default-reads** (transformation, applied to every site). These are mechanical:

| Pattern | Replacement |
|---------|-------------|
| `class_read_double("X", DEST);` | `DEST = pfc->get_or("X", DEST);` |
| `class_read_int("X", DEST);` (DEST is `int`) | `DEST = pfc->get_or("X", DEST);` |
| `class_read_string("X", DEST);` | `DEST = pfc->get_or("X", DEST);` |
| `class_read_double_one_of_two("A", "B", DEST);` | `if (auto v = read_one_of<double>(*pfc, "A", "B")) DEST = *v;` |

The **one** enum destination is `pnl->extrapolation_method` (type `enum source_extrapolation`):
```cpp
pnl->extrapolation_method =
    (enum source_extrapolation) pfc->get_or<int>("extrapolation_method", pnl->extrapolation_method);
```

- [ ] **Step 3: Migrate the raw `parser_read_*` + flag blocks.** Each `parser_read_TYPE(pfc, "X", &paramN, &flagN); ... if (flagN == _TRUE_)` becomes a local optional. Canonical example — the IDR block (around lines 329–376):

```cpp
auto N_idr  = pfc->get<double>("N_idr");
auto N_dg   = pfc->get<double>("N_dg");
auto xi_idr = pfc->get<double>("xi_idr");
class_test(n_present(N_idr, N_dg, xi_idr) > 1,
           "In input file, you can only enter one of N_idr, N_dg or xi_idr, choose one");

double T_idr_local = 0.;
if (N_idr) {
  T_idr_local = pow(*N_idr / stat_f_idr * (7. / 8.) / pow(11. / 4., (4. / 3.)), (1. / 4.)) * pba->T_cmb;
  if (input_verbose > 1) printf("You passed N_idr = N_dg = %e, ...\n", *N_idr, T_idr_local / pba->T_cmb);
}
else if (N_dg) {
  T_idr_local = pow(*N_dg  / stat_f_idr * (7. / 8.) / pow(11. / 4., (4. / 3.)), (1. / 4.)) * pba->T_cmb;
  if (input_verbose > 2) printf("You passed N_dg = N_idr = %e, ...\n", *N_dg, T_idr_local / pba->T_cmb);
}
else if (xi_idr) {
  T_idr_local = *xi_idr * pba->T_cmb;
  if (input_verbose > 1) printf("You passed xi_idr = %e, ...\n", *xi_idr, stat_f_idr * pow(*xi_idr, 4.) / (7./8.) * pow(11./4., 4./3.));
}
if (N_idr || N_dg || xi_idr) coupled_inputs_.T_idr = T_idr_local;
// ... and the omega_budget_.idr line guarded by the same `if (N_idr || N_dg || xi_idr)`
```

Apply the same shape to every other raw block. Replace each remaining mutual-exclusion macro:
- `class_at_least_two_of_three(f1,f2,f3)` → `n_present(o1,o2,o3) > 1`
- `class_at_least_two_of_four(f1,f2,f3,f4)` → `n_present(o1,o2,o3,o4) > 1`
- `class_none_of_three(f1,f2,f3)` → `n_present(o1,o2,o3) == 0`

where `oN` are the corresponding `std::optional` locals.

- [ ] **Step 4: Rewire the local precision overloads** (bodies only; their 175 call sites are unchanged) near line 2490:

```cpp
void read(const FileContent& fc, const char* name, double& v)      { v = fc.get_or(name, v); }
void read(const FileContent& fc, const char* name, int& v)         { v = fc.get_or(name, v); }
void read(const FileContent& fc, const char* name, std::string& v) { v = fc.get_or(name, v); }
template <typename E>
void read_enum(const FileContent& fc, const char* name, E& v)      { v = static_cast<E>(fc.get_or<int>(name, static_cast<int>(v))); }
```
And `readDoubleList` (line 31): `if (auto l = pfc->get<std::vector<double>>(name)) { values = *l; *found = _TRUE_; } else *found = _FALSE_;` — or migrate its callers to `get<std::vector<double>>` directly and delete the helper.

- [ ] **Step 5: Delete the now-dead locals and macros.**
  - In `source/input_module.cpp`, remove the shared `param1/param2/param3`, `int1`, `string1`, `flag1/flag2/flag3` declarations from every Read function (they are now unreferenced).
  - In `source/input_module.h`, delete `class_read_double`, `class_read_int`, `class_read_string`, `class_read_double_one_of_two`, `class_at_least_two_of_three`, `class_at_least_two_of_four`, `class_none_of_three`, `class_any_nonzero_four`, `class_all_nonzero_four`.

- [ ] **Step 6: Verify completeness + build + smoke**

```bash
grep -nE "class_read_|class_at_least_|class_none_of|parser_read_(int|double|string)|\bflag[1-9]\b|\bparam[1-9]\b|\bint1\b|\bstring1\b" source/input_module.cpp source/input_module.h
```
Expected: **no matches** (all migrated/deleted).
```bash
cmake --build . -j --target class && ./class explanatory.ini && echo SMOKE_OK
```
Expected: builds clean, `explanatory.ini` runs to completion.

- [ ] **Step 7: Commit**

```bash
git add source/input_module.cpp source/input_module.h
git commit -m "#314: migrate input_module reads to get<T>/get_or; delete class_read_* macros"
```

---

### Task 3: Migrate direct `FileContent::read_*` callers in species files

**Files (modify):** `species/dcdm_dr_species.cpp`, `species/fluid.cpp`, `species/idm_dr_idr_species.cpp`, `species/ncdm_species.cpp`, `species/ncdm_interacting_species.cpp`, `species/dncdm_species.cpp`, `species/scalar_field.cpp`, `species/ultra_relativistic.cpp`, `species/greybody_ncdm_species.cpp`, `species/species_input.cpp` (the internal `pfc->read_string` at line ~118).

- [ ] **Step 1: Apply the transformation** to every direct `FileContent` read (`ctx.pfc->`, `fc.`, `pfc.`, `pfc->`):

| Pattern | Replacement |
|---------|-------------|
| `bool has_X = P->read_double("X", v);` then uses `has_X` + `v` | `auto X = P->get<double>("X");` then use `X.has_value()` / `*X` |
| `if (P->read_string("X", s)) { ... s ... }` | `if (auto s_ = P->get<std::string>("X")) { ... *s_ ... }` |
| `P->read_string(name + ".type", unused);  // mark consumed` | `(void) P->get<std::string>(name + ".type");  // mark consumed` |
| `P->read_double("X", v);` (value used after, default kept if absent) | `v = P->get_or("X", v);` |
| `P->read_list_of_doubles("X", vec)` (as bool) | `P->get<std::vector<double>>("X")` (as optional) |

Concrete example — `species/ultra_relativistic.cpp:438`:
```cpp
auto n_ur_opt  = ctx.pfc->get<double>("N_ur");
auto n_eff_opt = ctx.pfc->get<double>("N_eff");
bool flag_nur  = n_ur_opt.has_value();
bool flag_neff = n_eff_opt.has_value();
// then use *n_ur_opt / *n_eff_opt where n_ur / n_eff were used
```
Concrete example — `species/idm_dr_idr_species.cpp:598` (list):
```cpp
auto alpha = ctx.pfc->get<std::vector<double>>("alpha_idm_dr");
if (!alpha) alpha = ctx.pfc->get<std::vector<double>>("alpha_dark");
bool found = alpha.has_value();
if (found) alpha_idm_dr = *alpha;
```

- [ ] **Step 2: Verify no direct FileContent `read_*` remain, then build**

```bash
grep -rnE "(pfc|fc|pfc_|ctx\.pfc)(\.|->)read_(int|double|string|list)" source/ species/ | grep -v "input.read_"
```
Expected: **no matches**.
```bash
cmake --build . -j --target class test-parser && ./test-parser && echo OK
```
Expected: builds, `OK`.

- [ ] **Step 3: Commit**

```bash
git add species/*.cpp
git commit -m "#314: migrate direct FileContent::read_* callers in species to get<T>"
```

---

### Task 4: Reshape `SpeciesInput` to `get<T>`/`require<T>` and migrate its callers

**Files:**
- Modify: `species/species_input.h`, `species/species_input.cpp`
- Modify: `species/ncdm_interacting_species.cpp`, `species/dncdm_species.cpp`, `species/ncdm_base_species.cpp`, `species/greybody_ncdm_species.cpp` (the `input.read_*` sites)
- Test: `tools/parser_test.cpp`

- [ ] **Step 1: Update the `SpeciesInput` tests first** in `tools/parser_test.cpp` (`test_species_input_prefixing`, `test_species_input_required_throws`):

```cpp
// prefixing test:
auto m = input.get<double>("m");
assert(m && *m == 0.06);
assert(fc.was_read("nu1.m"));
assert(!input.get<double>("deg").has_value());
assert(!fc.was_read("nu1.deg"));
auto t = input.get<std::string>("type");
assert(t && *t == "ncdm_standard");

// required -> require test:
try { (void) input.require<double>("m"); }
catch (const std::invalid_argument& e) {
  std::string msg = e.what();
  assert(msg.find("nu1") != std::string::npos && msg.find("m") != std::string::npos);
  threw = true;
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cmake --build . -j --target test-parser`
Expected: FAIL — `get`/`require` not members of `SpeciesInput`.

- [ ] **Step 3: Reshape the header** `species/species_input.h` — replace the four `read_*` and three `required_*` declarations with templates (defined inline so they compile against `FileContent::get<T>` declared in `parser.h`):

```cpp
#include <optional>

  template <class T> std::optional<T> get(const std::string& field) const {
    return pfc_->get<T>(qualify(field));
  }
  template <class T> T require(const std::string& field) const {
    if (auto v = get<T>(field)) return *v;
    throw std::invalid_argument(
        "Required parameter '" + qualify(field) + "' not found in " + pfc_->get_filename());
  }
```
(Match the exact message text of the current `required_*` implementations in `species_input.cpp` so behavior — and the test's substring checks — are preserved. Add `#include <stdexcept>` if not already present.) Remove the now-obsolete `read_*`/`required_*` declarations.

- [ ] **Step 4: Update `species/species_input.cpp`** — delete the `read_double/read_int/read_string/read_list_of_doubles` and `required_double/required_int/required_string` definitions (now header templates). Update the internal users `CollectInstanceFieldValues` (`input.read_string(field, values[i])` → `if (auto v = input.get<std::string>(field)) values[i] = *v;`) and `SynthesiseIdenticalScalarField` (`pfc->read_string(legacy_key, existing_value)` → `auto existing = pfc->get<std::string>(legacy_key); if (existing && *existing != first) ...`).

- [ ] **Step 5: Migrate the 35 `input.read_*` call sites** in the four species files using the Task-3 patterns (`input.read_double("X", v)` → `auto x = input.get<double>("X")`).

- [ ] **Step 6: Verify, build, run tests**

```bash
grep -rnE "\binput\.(read_|required_)" source/ species/
grep -rnE "::read_(int|double|string|list)|::required_(double|int|string)" species/species_input.cpp
```
Expected: **no matches**.
```bash
cmake --build . -j --target test-parser class && ./test-parser && ./class explanatory.ini && echo OK
```
Expected: builds, `OK`, smoke completes.

- [ ] **Step 7: Commit**

```bash
git add species/species_input.* species/*.cpp tools/parser_test.cpp
git commit -m "#314: reshape SpeciesInput to get<T>/require<T> and migrate callers"
```

---

### Task 5: Delete the dead surfaces (`FileContent::read_*`, `parser_read_int/double/string`)

**Files:** `include/parser.h`, `tools/parser.cpp`

- [ ] **Step 1: Confirm there are zero remaining callers** (defence before deletion):

```bash
grep -rnE "\.read_(int|double|string|list)|->read_(int|double|string|list)|parser_read_(int|double|string)" source/ species/ tools/ main/ | grep -v "FileContent::read_\|::get<"
grep -rnE "read_int|read_double|read_string|read_list|parser_read_(int|double|string)" python/classy.pyx python/cclassy.pxd python/generate_wrapper.py
```
Expected: first command shows only the definitions in `tools/parser.cpp` (about to be deleted); second shows **no matches** (Cython does not bind them).

- [ ] **Step 2: Delete the declarations** in `include/parser.h`: the six `read_*` member decls; the `parser_read_int`, `parser_read_double` lines inside the `extern "C"` block; and the standalone `int parser_read_string(...)` decl below the block. **Keep** `parser_read_file` and `parser_cat`.

- [ ] **Step 3: Delete the definitions** in `tools/parser.cpp`: the six `FileContent::read_*` wrapper bodies and the `parser_read_int/double/string` free-function bodies. Keep `parser_read_file` / `parser_cat`. (`get<T>` already owns all the parse logic from Task 1.)

- [ ] **Step 4: Build everything + run all unit tests**

```bash
cmake --build . -j && ctest --output-on-failure
```
Expected: full build clean (no undefined-reference / missing-decl errors anywhere), all of `test-parser`/`test-bisection`/`test-photons` pass.

- [ ] **Step 5: Commit**

```bash
git add include/parser.h tools/parser.cpp
git commit -m "#314: delete FileContent::read_* and parser_read_int/double/string"
```

---

### Task 6: Byte-identical regression + targeted spot-checks

**Files:** none (verification only).

- [ ] **Step 1: Rebuild the Python extension** (the regression imports `classy`):

Run: `pip install . --no-build-isolation` (from repo root; matches the project's build recipe). Confirm `classyref` is installed for the comparison (per the lvl2 setup); if not, the COMPARE_OUTPUT_REF gate self-disables and the run is not a valid witness.

- [ ] **Step 2: Run the PR-level regression** (the exact CI command from `.github/workflows/test_on_pull_request.yml`):

```bash
cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py
```
Expected: all scenarios pass, 0 failures. Because parsing produces identical values, output is byte-identical to `classyref`.

- [ ] **Step 3: Targeted spot-check** on the migrated branches — run `./class` on inputs exercising (a) the IDR `N_idr/N_dg/xi_idr` mutual-exclusion, (b) a `class_read_double_one_of_two` parameter (e.g. `N_ur`/`N_eff`), (c) a list parameter (e.g. `z_pk`), and confirm the produced outputs match a pre-change build (e.g. compare against the same run on `master`). Also run a dot-syntax species `.ini` from `test/` (e.g. `test/dotsyntax_ncdm.ini`) to exercise the reshaped `SpeciesInput`.

- [ ] **Step 4: No commit** (verification task). If any diff appears, STOP and debug — a non-identical result means a migration changed semantics (likely a default-kept-vs-overwritten slip or a list-parse difference).

---

### Task 7: Format + final housekeeping

- [ ] **Step 1: clang-format the touched files** (master is clang-format clean):

```bash
clang-format -i include/parser.h tools/parser.cpp tools/parser_test.cpp \
  source/input_module.h source/input_module.cpp species/species_input.h species/species_input.cpp \
  species/dcdm_dr_species.cpp species/fluid.cpp species/idm_dr_idr_species.cpp species/ncdm_species.cpp \
  species/ncdm_interacting_species.cpp species/dncdm_species.cpp species/scalar_field.cpp \
  species/ultra_relativistic.cpp species/greybody_ncdm_species.cpp species/ncdm_base_species.cpp
```

- [ ] **Step 2: Rebuild once more to confirm formatting did not break anything**

Run: `cmake --build . -j --target class test-parser && ./test-parser && echo OK`
Expected: builds, `OK`.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "#314: clang-format touched files"
```

- [ ] **Step 4: Open the PR** referencing #314 and the spec; note "clean over bit-identical, but this refactor IS byte-identical (parsing only); full TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 green."

---

## Self-review notes (author)

- **Spec coverage:** API (Task 1), macro deletion + input_module migration (Task 2), direct-caller sweep (Task 3), SpeciesInput reshape (Task 4), deletions incl. Cython check (Task 5), byte-identical verification (Task 6), formatting (Task 7). All spec sections mapped.
- **Type consistency:** `get<T>` returns `std::optional<T>` everywhere; `get_or<T>` returns `T`; `require<T>` returns `T` and throws `std::invalid_argument`; `n_present(...)` returns `int`; `read_one_of<T>` returns `std::optional<T>`. Names identical across tasks.
- **Build-green invariant:** `read_*` survive as wrappers from Task 1 until all callers are gone (Tasks 2–4), deleted only in Task 5 — every intermediate task builds.
- **Known risk:** the "keep-default-if-absent" semantics. `get_or(name, dest)` reproduces `class_read_*`/precision `read()` exactly (no write on absence); the raw-block migration must preserve each original default. Task 6's spot-check is the guard.
