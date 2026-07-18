# Error-Severity Conventions Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `class_test_severe`/`class_stop_severe` (throwing `std::invalid_argument`), remove `InputModule`'s blanket `runtime_error`→`invalid_argument` rewrap, and reclassify every error site by the structure-vs-value rule so MCMC samplers abort only on structurally bad input and reject points on value/numerical failures.

**Architecture:** Spec: `docs/superpowers/specs/2026-07-17-error-severity-conventions-design.md` (read it first). Exception type is the severity channel: `invalid_argument` → classy `CosmoSevereError` (abort), `runtime_error` → `CosmoComputationError` (reject point), `logic_error` → programmer error (maps to severe via classy fallback). The governing rule: **a severe check must depend only on input *structure* (key presence, mutual exclusivity, unparseable strings, unknown enum names, list-length/count consistency, API-argument validation) — never on a parsed numeric value, because samplers vary numeric values.** All range checks on numbers and all numerical failures are computation errors, wherever they occur (input phase and shooting trial builds included).

**Tech Stack:** C++17 (no C-compat guards — this repo is C++-only except hyrec), CMake via `make` shim, assert-based unit tests, Cython wrapper (`classy.pyx` is hand-edited; `cclassy.pxd` is auto-generated — NEVER edit it).

## Global Constraints

- Never `git stash` (shared stack) and never `git checkout`/`git reset` — you share the working tree. Branch: `refactor/error-conventions` (already checked out).
- Never `git add -A` or `git add .` — stage explicit paths only.
- New test executables must be registered in BOTH `CMakeLists.txt` (both the `add_executable` line AND the `foreach(_t IN ITEMS ...)` list around line 235) and the `Makefile` `TEST_TARGETS` list (line 4). CI runs `make test` and only builds targets named there.
- Build: `make <target>` (CMake shim, build dir `build/cmake`). Full test suite: `make test`. Single test after build: `./build/cmake/test-errors` (run from repo root — tests assume the repo root as working directory).
- Message format is load-bearing and must be preserved exactly: `"<func>(L:<line>) :<prefix><body>"` with prefix `"condition (<stringified condition>) is true; "` for `class_test*` and `"error; "` for `class_stop*`.
- Error-path changes only: numerical output must not change. No golden/classyref regeneration.
- Line numbers below were taken at commit `df0e3824` and will drift as tasks land. Always locate sites by the quoted content anchor, not the line number.
- Commit messages end with:

```
Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01LMSkqYYndXZdTB6sy8D7co
```

---

### Task 1: Severe macros, `ThrowFormattedSevere`, and `errors_test`

**Files:**
- Modify: `include/errors.h`
- Modify: `tools/common.cpp` (contains `ThrowFormatted`, near the top)
- Create: `tools/errors_test.cpp`
- Modify: `CMakeLists.txt` (test registration ~lines 217–235)
- Modify: `Makefile` (TEST_TARGETS list, line 4)

**Interfaces:**
- Produces: `[[noreturn]] void ThrowFormattedSevere(const char* func, int line, const char* prefix, const char* fmt, ...)`; macros `class_test_severe(condition, args, ...)` and `class_stop_severe(args, ...)`. All later tasks use these macros.

- [ ] **Step 1: Write the failing test** — create `tools/errors_test.cpp`:

```cpp
#include "errors.h"

#include <cassert>
#include <stdexcept>
#include <string>

// class_test: runtime_error, message carries func name, condition text, %g-formatted arg.
static void test_class_test_type_and_message() {
  bool threw    = false;
  const double q = 1.5e-12;
  try {
    class_test(q > 0., "q = %g must not be positive", q);
  }
  catch (const std::invalid_argument&) {
    assert(false && "class_test must throw runtime_error, not invalid_argument");
  }
  catch (const std::runtime_error& e) {
    threw = true;
    const std::string msg = e.what();
    assert(msg.find("test_class_test_type_and_message") != std::string::npos);
    assert(msg.find("condition (q > 0.) is true") != std::string::npos);
    assert(msg.find("1.5e-12") != std::string::npos);
  }
  assert(threw);
}

// class_test_severe: identical format, but invalid_argument.
static void test_class_test_severe_type_and_message() {
  bool threw = false;
  try {
    class_test_severe(1 == 1, "keys '%s' and '%s' are mutually exclusive", "a", "b");
  }
  catch (const std::invalid_argument& e) {
    threw = true;
    const std::string msg = e.what();
    assert(msg.find("condition (1 == 1) is true") != std::string::npos);
    assert(msg.find("keys 'a' and 'b' are mutually exclusive") != std::string::npos);
  }
  assert(threw);
}

// invalid_argument IS-A logic_error, NOT runtime_error: severe must not be
// catchable as runtime_error (that would re-blur the severity channel).
static void test_severe_is_not_runtime_error() {
  bool caught_as_runtime = false;
  try {
    class_stop_severe("structural failure");
  }
  catch (const std::runtime_error&) {
    caught_as_runtime = true;
  }
  catch (const std::exception&) {
  }
  assert(!caught_as_runtime);
}

static void test_class_stop_and_stop_severe() {
  bool threw = false;
  try {
    class_stop("gave up after %d iterations", 20);
  }
  catch (const std::runtime_error& e) {
    threw = true;
    assert(std::string(e.what()).find("error; gave up after 20 iterations") != std::string::npos);
  }
  assert(threw);

  threw = false;
  try {
    class_stop_severe("incomprehensible input '%s'", "spam");
  }
  catch (const std::invalid_argument& e) {
    threw = true;
    assert(std::string(e.what()).find("error; incomprehensible input 'spam'") !=
           std::string::npos);
  }
  assert(threw);
}

static void test_false_condition_does_not_throw() {
  class_test(false, "unreachable");
  class_test_severe(false, "unreachable");
}

int main() {
  test_class_test_type_and_message();
  test_class_test_severe_type_and_message();
  test_severe_is_not_runtime_error();
  test_class_stop_and_stop_severe();
  test_false_condition_does_not_throw();
  return 0;
}
```

- [ ] **Step 2: Register the test.** In `CMakeLists.txt`, after `add_executable(test-parser tools/parser_test.cpp)` add:

```cmake
  add_executable(test-errors tools/errors_test.cpp)
```

and add `test-errors` to the `foreach(_t IN ITEMS ...)` list on the same line as the other test names. In `Makefile`, add `	test-errors \` as a new line in the `TEST_TARGETS :=` block (keep the trailing-backslash style; the last entry has no backslash).

- [ ] **Step 3: Run to verify it fails.** Run: `make test-errors 2>&1 | tail -20`
Expected: compile FAILURE — `class_test_severe` / `class_stop_severe` not defined.

- [ ] **Step 4: Implement.** In `tools/common.cpp`, replace the existing `ThrowFormatted` definition with a shared formatter plus two entry points (the format string logic is copied verbatim from the current implementation):

```cpp
namespace {

std::string FormatThrowMessage(
    const char* func, int line, const char* prefix, const char* fmt, va_list args) {
  va_list args_copy;
  va_copy(args_copy, args);
  const int n = std::vsnprintf(nullptr, 0, fmt, args_copy);
  va_end(args_copy);

  std::string body;
  if (n > 0) {
    body.resize(static_cast<std::size_t>(n));
    std::vsnprintf(body.data(), static_cast<std::size_t>(n) + 1, fmt, args);
  }

  // Preserves the historical "<func>(L:<line>) :<message>" format.
  return std::string(func) + "(L:" + std::to_string(line) + ") :" + prefix + body;
}

}  // namespace

[[noreturn]] void ThrowFormatted(
    const char* func, int line, const char* prefix, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  const std::string message = FormatThrowMessage(func, line, prefix, fmt, args);
  va_end(args);
  throw std::runtime_error(message);
}

[[noreturn]] void ThrowFormattedSevere(
    const char* func, int line, const char* prefix, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  const std::string message = FormatThrowMessage(func, line, prefix, fmt, args);
  va_end(args);
  throw std::invalid_argument(message);
}
```

In `include/errors.h`, add after the `ThrowFormatted` declaration:

```cpp
[[noreturn]] void ThrowFormattedSevere(
    const char* func, int line, const char* prefix, const char* fmt, ...);
```

and after the `class_test`/`class_stop` macros:

```cpp
/* Severe variants: throw std::invalid_argument (classy: CosmoSevereError, the
   sampler ABORTS). Use ONLY for checks that depend exclusively on the
   STRUCTURE of the input: key presence, mutual exclusivity, unparseable
   strings, unknown enum names, list-length/count consistency, API-argument
   validation. A severe check must never depend on possibly-varying (numeric)
   parameters: range checks on parsed values and all numerical failures use
   class_test/class_stop (std::runtime_error -> CosmoComputationError, the
   sampler rejects the point and the chain survives). */

#define class_test_severe(condition, args, ...)                    \
  {                                                                \
    if (condition) {                                               \
      ThrowFormattedSevere(__func__,                               \
                           __LINE__,                               \
                           "condition (" #condition ") is true; ", \
                           args,                                   \
                           ##__VA_ARGS__);                         \
    }                                                              \
  }

#define class_stop_severe(args, ...)                                          \
  {                                                                           \
    ThrowFormattedSevere(__func__, __LINE__, "error; ", args, ##__VA_ARGS__); \
  }
```

- [ ] **Step 5: Run to verify it passes.** Run: `make test-errors && ./build/cmake/test-errors && echo PASS`
Expected: `PASS`

- [ ] **Step 6: Full suite still green.** Run: `make test 2>&1 | tail -5`
Expected: `100% tests passed`

- [ ] **Step 7: Commit**

```bash
git add include/errors.h tools/common.cpp tools/errors_test.cpp CMakeLists.txt Makefile
git commit -m "errors: add class_test_severe/class_stop_severe (invalid_argument channel)"
```

---

### Task 2: Retire `ThrowRuntimeErrorIf`/`ThrowRuntimeError`

**Files:**
- Modify: `tools/exceptions.h`, `tools/exceptions.cpp`
- Modify: `source/lensing_module.cpp` (~lines 35, 58), `source/spectra_module.cpp` (~lines 100, 104, 156, 157)
- Modify: `species/ncdm_base_species.cpp` (~line 491)

**Interfaces:**
- Consumes: Task 1 macros.
- Produces: `tools/exceptions.h` keeps ONLY `get_my_py_error_message()` (classy needs it).

- [ ] **Step 1: Convert the 7 call sites.** These are mechanical: `ThrowRuntimeErrorIf(cond, msg, ...)` has the same polarity as `class_test(cond, msg, ...)`. Keep every condition and message string byte-identical. The 6 lensing/spectra sites are API-argument validation (fixed per-analysis configuration a sampler cannot vary; classy raises `CosmoSevereError` for the same conditions in its own Python-level checks) → `class_test_severe`. Examples (repeat the same transformation for all six):

```cpp
// lensing_module.cpp:35  (anchor: "lmax > l_lensed_max_")
class_test_severe((lmax > l_lensed_max_) || (lmax < 0), /* original message args unchanged */ ...);
// lensing_module.cpp:58  (anchor: "No lensed Cls was computed")
class_test_severe(!ple->has_lensed_cls, "No lensed Cls was computed, adjust your inputs.\n");
```

The NCDM Newton-mass failure is a numerical failure → computation error:

```cpp
// ncdm_base_species.cpp:491 (anchor: "Newton iteration could not converge on a mass")
class_stop("Newton iteration could not converge on a mass for some reason.");
```

Ensure each of the three .cpp files has `#include "errors.h"` (they may have relied on a transitive include of `exceptions.h`).

- [ ] **Step 2: Delete the two functions** from `tools/exceptions.h` (keep `get_my_py_error_message`) and their definitions from `tools/exceptions.cpp` (keep the include of `<typeinfo>` etc. needed by the surviving function; drop `<cstdarg>` if now unused).

- [ ] **Step 3: Verify.** Run: `make test 2>&1 | tail -5` — Expected: `100% tests passed`. Then `grep -rn "ThrowRuntimeError" --include="*.cpp" --include="*.h" source species tools main include` — Expected: no matches.

- [ ] **Step 4: Commit**

```bash
git add tools/exceptions.h tools/exceptions.cpp source/lensing_module.cpp source/spectra_module.cpp species/ncdm_base_species.cpp
git commit -m "errors: retire ThrowRuntimeErrorIf/ThrowRuntimeError third mechanism"
```

---

### Task 3: classy fallback maps unknown exception types to `CosmoSevereError`

**Files:**
- Modify: `classy.pyx` (`raise_my_py_error`, ~line 90–101)

- [ ] **Step 1: Edit.** In `classy.pyx`, replace the `else` branch of `raise_my_py_error`:

```python
    else:
        # logic_error (programmer-error invariants) and anything unrecognized:
        # abort loudly rather than surfacing a misleading NotImplementedError.
        raise CosmoSevereError(cpp_exception.second)
```

(The line being replaced is `raise NotImplementedError(cpp_exception.second)`.)

- [ ] **Step 2: Commit.** (Compile/e2e verification of the wrapper is deferred to Task 10, which builds classy.)

```bash
git add classy.pyx
git commit -m "classy: map unknown C++ exception types to CosmoSevereError"
```

---

### Task 4: InputModule — remove the rewrap, fix `readDoubleList`, hook guard → `logic_error`

**Files:**
- Modify: `source/input_module.cpp`

- [ ] **Step 1: Remove the constructor rewrap** (anchor: `catch (const std::runtime_error& e)` inside `InputModule::InputModule`). Replace the whole try/catch with the plain sequence:

```cpp
InputModule::InputModule(FileContent& fc) : file_content_(fc) {
  file_content_.mark_all_unread();
  // Translate dot-syntax for single-instance legacy species (bary.Omega ->
  // Omega_b) before any consumer reads the file.
  TranslateSingleInstanceDotSyntax(&file_content_);
  input_read_precisions();
  ReadContext();
  ConstructSpecies();
  ReadDerived();
  WriteParameterFiles();
}
```

- [ ] **Step 2: Fix `readDoubleList`** (anchor: `class_stop("%s", e.what())`). It currently catches the parser's correctly-typed `invalid_argument` and demotes it to `runtime_error`. Delete the try/catch, keep the logic:

```cpp
void readDoubleList(FileContent* pfc, const char* name, std::vector<double>& values, int* found) {
  if (auto l = pfc->get<std::vector<double>>(name)) {
    values = *l;
    *found = true;
  }
  else
    *found = false;
}
```

- [ ] **Step 3: Shooting hook-contract guard → `logic_error`** (anchor: `" reported "` + `" shooting target(s) but "`, in `DoShooting`). A mis-implemented species hook is a programmer error:

```cpp
    if (g.size() != tgts.size() || d.size() != tgts.size()) {
      throw std::logic_error("species '" + key + "' reported " + std::to_string(tgts.size()) +
                             " shooting target(s) but " + std::to_string(g.size()) +
                             " guess(es) / " + std::to_string(d.size()) + " Jacobian seed(s)");
    }
```

(`<stdexcept>` is already available transitively; add the include if the build says otherwise.)

- [ ] **Step 4: Verify.** Run: `make class && ./class explanatory.ini 2>&1 | tail -3` — Expected: normal completion (writes output files, no exception). If `./class` does not exist at the repo root, the binary is `./build/cmake/class`. Then `make test 2>&1 | tail -5` — Expected: `100% tests passed`.

- [ ] **Step 5: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: drop blanket runtime_error->invalid_argument rewrap"
```

---

### Task 5: `input_module.cpp` severity audit

**Files:**
- Modify: `source/input_module.cpp` only.

The rule (from the spec — read its "The convention" section first): structural → `_severe`; anything testing a parsed numeric value → stays plain `class_test` (unchanged). Locate every site by content anchor; line numbers are from `df0e3824` and have shifted slightly after Task 4.

- [ ] **Step 1: Convert every `class_stop(` in this file to `class_stop_severe(`.** After Task 4 there are 11, all "incomprehensible/unclear input" structural cases (anchors: `'on the spot'`, `'compute damping scale'`, ~1687, `'inflation behavior'`, ~1797, `In selection function input`, ~1923, ~1950, ~1977, `non_diagonal`, ~2087). Before converting each, confirm its message describes a parse/structure failure, not a numeric-value failure — if you find a value case, leave it and note it in the commit message.

- [ ] **Step 2: Convert these structural `class_test(` sites to `class_test_severe(`** (49 sites; anchor = message fragment):

| ~line | anchor |
|---|---|
| 56 | `you can only enter one of %s, %s` |
| 113, 118 | `!input_file.empty()` / `!precision_file.empty()` |
| 334 | `n_present(N_idr, N_dg, xi_idr) > 1` |
| 391 | `(Omega_cdm && omega_cdm)` |
| 417 | `n_present(Omega_idm_dr, omega_idm_dr, f_idm_dr) > 1` |
| 452 | `(Omega_dcdmdr && omega_dcdmdr)` |
| 474 | `any_drmd && !all_drmd` |
| 665 | `cannot enter both h and H0` |
| 680 | `100*theta_s and h (or H0)` |
| 687 | `n_present(T_cmb, Omega_g, omega_g) > 1` |
| 716 | `(Omega_b && omega_b)` |
| 738 | `G_eff_ur.has_value()` |
| 781, 816 | `Omega_Lambda && Omega_fld` / `Omega_scf` (the sign of `Omega_scf` is a mode-selection flag, not a sampled value → structural) |
| 918 | `!recognized` |
| 928 | `(z_reio && tau_reio)` |
| 951, 955, 967, 971, 982, 986 | binned_reio/many_tanh/reio_inter size mismatches |
| 1131 | `switch_sw == 0) && (ppt->switch_eisw == 0` |
| 1152 | `!ppt->has_nc_density` |
| 1187 | `!ppt->has_scalars && !ppt->has_vectors` |
| 1214 | `!ppt->has_ad && !ppt->has_bi` |
| 1222, 1225 | `has_cl_cmb_lensing_potential` / `has_pk_matter` |
| 1230, 1237 | `!ppt->has_cl_cmb_temperature` |
| 1272 | `!recognized` |
| 1396 | `n_present(A_s, ln_A_s, sigma8, S8) > 1` |
| 1703 | `ln_aH_ratio_str && N_star_str` |
| 1738 | `omitted to write a command` |
| 1756, 1759, 1761, 1764 | inflationary-module mode checks |
| 1812 | `full_limber` yes/no parse |
| 1822 | `P_k_max_h_Mpc && P_k_max_1_Mpc` |
| 1837 | `n_list > _Z_PK_NUM_MAX_` |
| 1881 | `n_list > _SELECTION_NUM_MAX_` |
| 2016, 2022 | `"bias"` / `"s_bias"` forbidden keys |
| 2097 | `!ppt->has_perturbations` |
| 2144 | `(feedback_model && (eta_0 \|\| c_min))` |
| 2283 | `l_switch_limber_for_cl_density_over_z` retired key |
| 2318 | `n_list > _MAX_NUMBER_OF_K_FILES_` |

- [ ] **Step 3: Leave these unchanged** (value checks and internal invariants — listed so you know they were classified, not missed): ~427, 432, 488, 497, 502, 593, 1282–1354 (k1/k2/prr/pii/pri positivity + the four `should NEVER happen` internal invariants), 1347, 1349, 1405 (`sigma8`), 1413 (`S8`), 1596–1674 (inflation potential values), 1888, 1896 (selection values/ordering), 2221–2240 (trigger tolerances), 2375–2401 (`l_max_*` floors).

- [ ] **Step 4: Verify.** `make class && ./class explanatory.ini 2>&1 | tail -3` (normal completion), `make test 2>&1 | tail -5` (`100% tests passed`). Sanity-check one severe path: `printf 'h = 0.67\nH0 = 67.0\n' > /tmp/bad_both.ini && ./class /tmp/bad_both.ini 2>&1 | tail -2` — Expected: error message containing `cannot enter both h and H0`.

- [ ] **Step 5: Commit**

```bash
git add source/input_module.cpp
git commit -m "input: classify class_test sites into severe (structural) vs computation (value)"
```

---

### Task 6: Species port — input plumbing and NCDM core

**Files:**
- Modify: `species/species_input.h`, `species/species_input.cpp`, `species/ncdm_family.cpp`, `species/ncdm_species.cpp`, `species/ncdm_base_species.cpp`

Conversion pattern for raw throws → macros (applies to this task and Tasks 7–8). Structural example, before:

```cpp
    throw std::invalid_argument("input sets both '" + legacy_key + "' and dot-syntax field '" +
                                dot_field + "'");
```

after (printf-style, `.c_str()` for std::string args, `%g` for doubles, `%d` for ints, `%zu` for size_t), keeping the existing `if` and replacing only the throw:

```cpp
    class_stop_severe("input sets both '%s' and dot-syntax field '%s'", legacy_key.c_str(), dot_field.c_str());
```

Never write `class_test_severe(true, ...)` — it stamps the useless prefix `condition (true) is true;` into the message. Where the `if` guards nothing but the throw, you may instead fold the condition into the macro:

```cpp
  class_test_severe(pfc->get<T>(legacy_key) && dot_present, "input sets both '%s' and dot-syntax field '%s'", legacy_key.c_str(), dot_field.c_str());
```

only when that is a pure mechanical inversion (no side effects in the branch); otherwise keep the `if` and use `class_stop_severe(...)` in its body. Keep message text byte-equivalent up to formatting.

- [ ] **Step 1: Classify and convert:**
  - `species_input.h` (anchor `missing required field`) → `class_stop_severe` inside the existing `if` (structural). `errors.h` include needed.
  - `species_input.cpp`: null-`FileContent` checks (anchors `SpeciesInput: null FileContent*`, `TranslateSingleInstanceDotSyntax: null FileContent*`) → `throw std::logic_error(...)` (programmer error, keep as raw throws). All others (empty instance name, dot-field integer parse, `must be identical for all`, `sets both`, `single-instance but was given`) → severe macros (structural).
  - `ncdm_family.cpp`: null-pointer check (anchor `SynthesiseNcdmFluidApproximation: null FileContent*`) → `logic_error`; type-not-a-consumer check (anchor `is not one of the types that reads`) → severe (structural).
  - `ncdm_species.cpp` (anchors `ncdm_psd_filenames has`, `has " + std::to_string(values.size())` legacy-list length) → severe (count consistency).
  - `ncdm_base_species.cpp` ~82, ~90: read the surrounding code and classify by the rule (key/flag presence → severe; numeric range → plain `class_test`).

- [ ] **Step 2: Verify.** `make test 2>&1 | tail -5` — Expected: `100% tests passed` (these files' tests — test-parser, test-species-types, test-ncdm-family — all provoke structural errors, which stay `invalid_argument`).

- [ ] **Step 3: Commit**

```bash
git add species/species_input.h species/species_input.cpp species/ncdm_family.cpp species/ncdm_species.cpp species/ncdm_base_species.cpp
git commit -m "species: port input plumbing + NCDM core to severity macros"
```

---

### Task 7: Species port — NCDM-family species files and their tests

**Files:**
- Modify: `species/axion_ncdm_species.cpp`, `species/greybody_ncdm_species.cpp`, `species/dncdm_species.cpp`, `species/dcdm_dr_species.cpp`, `species/dcdm_wdm_species.cpp`, `species/wdm_decay_product.cpp`, `species/ncdm_interacting_species.cpp`, `species/ultra_relativistic.cpp`, `species/idm_dr_idr_species.cpp`, `species/type3_species.cpp`
- Modify tests: `species/axion_ncdm_test.cpp`, `species/dcdm_wdm_test.cpp`, `species/ncdm_family_test.cpp`

Use the conversion pattern from Task 6. Known classifications (verify each in context; anything not listed: apply the rule):

- `axion_ncdm_species.cpp`: `give either T or gstar_dec` and `one of T ... is required` → severe. `T ... must be positive`, `gstar_dec must exceed` → **plain `class_test`** (values; type flips from invalid_argument to runtime_error — intended behaviour change). `ksi is not supported` → plain `class_test` (ksi is numeric and sampleable; strict rule wins). `use_psd_file is incompatible` → severe (0/1 capability flag, not sampler-varyable).
- `greybody_ncdm_species.cpp`: the two `runtime_error` bracketing/moment-ratio failures (anchors `moment ratio r =`, `failed to bracket alpha`) → `class_stop` (stay computation; M2/M3 are sampleable). All positivity/range checks on `alpha`, `x`, `q0`, `M2`, `M3` → plain `class_test` (type flips). Key-presence/exclusivity checks → severe. The `EvaluatePsdAnalytic called before gb_ready_` `logic_error` stays as-is.
- `wdm_decay_product.cpp`: unknown-parameter check (anchor `parameter '`) → severe if it tests key validity; all numeric ranges (`momenta_bins out of range`, `q_min_ratio must be in`, etc.) → plain `class_test` (type flips).
- Remaining files (`dncdm`, `dcdm_dr`, `dcdm_wdm`, `ncdm_interacting`, `ultra_relativistic`, `idm_dr_idr`, `type3`): classify each throw by the rule — exclusivity/missing/unknown-string → severe; numeric comparisons → plain `class_test`.

- [ ] **Step 1: Update tests FIRST (TDD).** In `species/axion_ncdm_test.cpp`, the `Throws()` helper catches `invalid_argument`. Add a second helper and reassign the value cases:

```cpp
static bool ThrowsComputation(FileContent& fc) {
  background pba{};
  pba.H0 = 2.2e-4;
  try {
    AxionNCDMSpecies sp(&fc, "ax", TestSettings(), &pba, nullptr);
    return false;
  }
  catch (const std::runtime_error&) {
    return true;
  }
}
```

Keep `Throws(fc)` (invalid_argument) for the both-T-and-gstar and neither-T-nor-gstar cases; switch the `gstar_dec` below 43/11 case (anchor `"3.0"`), and any T≤0 / ksi≠0 cases later in the file, to `ThrowsComputation(fc)`. Apply the same split in `dcdm_wdm_test.cpp` (its catch at ~line 46 — check what the case provokes: structural stays, value moves to a runtime_error catch). `ncdm_family_test.cpp` (~line 35) provokes a structural type-consumer error → unchanged.

- [ ] **Step 2: Run tests to verify the moved cases FAIL.** Run: `make test-axion-ncdm && ./build/cmake/test-axion-ncdm; echo "exit=$?"`
Expected: assertion failure (`exit != 0`) — the value cases still throw `invalid_argument` because the species code is not yet ported.

- [ ] **Step 3: Port the ten species files** per the classifications above.

- [ ] **Step 4: Run the full suite.** `make test 2>&1 | tail -5` — Expected: `100% tests passed`.

- [ ] **Step 5: Commit**

```bash
git add species/axion_ncdm_species.cpp species/greybody_ncdm_species.cpp species/dncdm_species.cpp species/dcdm_dr_species.cpp species/dcdm_wdm_species.cpp species/wdm_decay_product.cpp species/ncdm_interacting_species.cpp species/ultra_relativistic.cpp species/idm_dr_idr_species.cpp species/type3_species.cpp species/axion_ncdm_test.cpp species/dcdm_wdm_test.cpp species/ncdm_family_test.cpp
git commit -m "species: port NCDM-family species to severity macros (value checks -> computation)"
```

---

### Task 8: Species port — fluid, scalar field, and their tests

**Files:**
- Modify: `species/fluid.cpp`, `species/scalar_field.cpp`
- Modify tests: `species/axion_ede_fluid_test.cpp`, `species/axion_scf_factory_test.cpp`

Known classifications (verify in context):

- `fluid.cpp`: `incomprehensible input ... for the field 'fluid equation of state'` → severe (unparseable string). `give only one of 'n_pheno_axion' and 'w_fld_f'`, `pheno_axion needs 'n_pheno_axion'`, `give only one of 'a_c' and 'log10_axion_ac'`, `pheno_axion needs 'a_c'`, two-density-keys/no-density-key, `use_ppf`, `cs2_fld` checks → severe (structural). `n_pheno_axion must be >= 1`, `nu_fld must be > 0`, `w_fld_i must be >= -1`, `pheno_axion needs w_fld_f > w_fld_i`, `a_c must lie in (0, 1)`, `Theta_initial_fld in (0, pi)`, `fraction_fld_ac must lie in (0, 1)` → **plain `class_test`** (sampled EDE parameters; type flips).
- `scalar_field.cpp`: `unknown scf_potential` → severe. Both `scf_parameters parse error` sites → severe (unparseable string). Wrong-parameter-count checks → severe (count consistency). `n_axion must be >= 1` and other numeric range checks → plain `class_test` (type flips). `attractor_ic_scf`-style yes/no or exclusivity checks → severe.

- [ ] **Step 1: Update tests FIRST.** In `species/axion_ede_fluid_test.cpp`, add next to `expect_create_throw`:

```cpp
  auto expect_create_throw_computation = [&](FileContent fc) {
    bool threw = false;
    try {
      run_create(fc);
    }
    catch (const std::runtime_error&) {
      threw = true;
    }
    assert(threw);
  };
```

Move the value cases to it: `Theta out of range` (anchor `"3.15"`) and `nu_fld must be > 0` (anchor `fc.set("nu_fld", "0")`), plus any other numeric-range rejection cases in the file. Structural cases (both n and w_fld_f, missing Theta [missing key = structural], both a_c and log10, two density keys, no density key, use_ppf, cs2_fld, closure-budget override) stay on `expect_create_throw`. In `species/axion_scf_factory_test.cpp` check the three catch sites (~111, 137, 153): factory/config structural cases stay; numeric-range cases (if any) move to a runtime_error-catching variant of `expect_throw`.

- [ ] **Step 2: Verify the moved cases fail.** `make test-axion-ede-fluid && ./build/cmake/test-axion-ede-fluid; echo "exit=$?"` — Expected: `exit != 0`.

- [ ] **Step 3: Port `fluid.cpp` and `scalar_field.cpp`** per the classifications.

- [ ] **Step 4: Full suite.** `make test 2>&1 | tail -5` — Expected: `100% tests passed`.

- [ ] **Step 5: Commit**

```bash
git add species/fluid.cpp species/scalar_field.cpp species/axion_ede_fluid_test.cpp species/axion_scf_factory_test.cpp
git commit -m "species: port fluid + scalar field to severity macros (value checks -> computation)"
```

---

### Task 9: Documentation — STYLE.md §8 and spec cross-link

**Files:**
- Modify: `STYLE.md` (§8, anchor `# 8. Error Handling` — the current text still describes the deleted `error_message_` member and is stale)

- [ ] **Step 1: Replace §8 wholesale with:**

```markdown
# 8. Error Handling

Errors are C++ exceptions thrown at the error origin via the macros in
`include/errors.h`. The exception type encodes severity; the Python wrapper
maps it to sampler behaviour (`classy.pyx: raise_my_py_error`):

* `class_test_severe` / `class_stop_severe` → `std::invalid_argument` →
  `CosmoSevereError`: the sampler **aborts**. Use ONLY for checks that depend
  exclusively on the *structure* of the input: key presence, mutual
  exclusivity, unparseable strings, unknown enum names, list-length/count
  consistency, API-argument validation.
* `class_test` / `class_stop` → `std::runtime_error` →
  `CosmoComputationError`: the sampler **rejects the point** and the chain
  survives. Use for everything that depends on a parsed numeric value (range
  checks included) and for every numerical failure, wherever it occurs — the
  input phase and shooting trial builds included.
* The governing rule: a severe check must never depend on possibly-varying
  parameters. Aborting a long chain on one unlucky proposal destroys work;
  rejecting is always the safe direction.
* Programmer-error invariants throw `std::logic_error` (surfaces as
  `CosmoSevereError` via the classy fallback).
* Do not introduce alternative error mechanisms.

Design rationale: `docs/superpowers/specs/2026-07-17-error-severity-conventions-design.md`.
```

- [ ] **Step 2: Commit**

```bash
git add STYLE.md
git commit -m "docs: STYLE.md error-handling section reflects severity conventions"
```

---

### Task 10: End-to-end verification through classy

**Files:** none modified (verification only; fix regressions if found).

- [ ] **Step 1: Full C++ suite.** `make test 2>&1 | tail -5` — Expected: `100% tests passed`.

- [ ] **Step 2: CLI regression.** `./class explanatory.ini 2>&1 | tail -3` and `./class test_axion_scf.ini 2>&1 | tail -3` — Expected: both complete normally.

- [ ] **Step 3: Build classy.** `make classy 2>&1 | tail -5` — Expected: successful pip install (this also regenerates `cclassy.pxd`; takes several minutes).

- [ ] **Step 4: Severity round-trip.** Save as `scratchpad/severity_check.py` (use the session scratchpad directory, not the repo) and run with the venv/python that received the pip install:

```python
from classy import Class, CosmoSevereError, CosmoComputationError

# 1. Structural error -> CosmoSevereError (sampler would abort).
c = Class()
c.set({'h': 0.67, 'H0': 67.0})
try:
    c.compute()
    raise SystemExit('FAIL: structural error did not raise')
except CosmoSevereError:
    print('structural -> CosmoSevereError OK')
c.struct_cleanup(); c.empty()

# 2. Value error -> CosmoComputationError (sampler would reject the point).
c = Class()
c.set({'N_ncdm': 1, 'ax.type': 'ncdm_axion', 'ax.m': 1.0, 'ax.gstar_dec': 2.0})
try:
    c.compute()
    raise SystemExit('FAIL: value error did not raise')
except CosmoComputationError:
    print('value -> CosmoComputationError OK')
c.struct_cleanup(); c.empty()
print('severity round-trip PASS')
```

Expected output ends with `severity round-trip PASS`. (Case 2 exercises the axion `gstar_dec must exceed 43/11` check, a value error that previously surfaced as severe — this is the intended behaviour change.)

- [ ] **Step 5: No stray mechanisms.** `grep -rn "NotImplementedError" classy.pyx` → only unrelated hits (if any); `grep -rn "throw std::invalid_argument" species/ | grep -v _test` → only sites deliberately left raw (none expected; report leftovers).

- [ ] **Step 6: Commit anything fixed during verification** (explicit paths), otherwise nothing to commit.

---

## Self-Review Notes

- Spec coverage: convention/macros (T1, T9), rewrap+readDoubleList+hook guard (T4), input audit (T5), species port incl. type flips + tests (T6–T8), ThrowRuntimeError retirement (T2), classy fallback (T3), e2e severe/computation round-trip (T10). Out-of-scope items (fzero damped Newton) deliberately absent.
- The `errors_test` message assertions pin the exact prefixes (`condition (...) is true; `, `error; `) so the shared-formatter refactor in T1 cannot silently change the format.
- Type names used across tasks: `ThrowFormattedSevere`, `class_test_severe`, `class_stop_severe` — consistent everywhere.
