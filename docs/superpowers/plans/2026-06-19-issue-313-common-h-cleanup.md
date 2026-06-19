# Issue #313 — common.h cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the 1064-line `include/common.h` grab-bag into `constants.h` / `errors.h` / `precision.h` (with `common.h` a thin umbrella), replace the `sprintf`/`ErrorMsg`-buffer error macros with a typed `ThrowFormatted`, delete the now-dead C error path, swap `MIN/MAX/SIGN/NRSIGN` for `<algorithm>`/`<cmath>`, drop `typeof`, and make the C++-only `enum`s `enum class`.

**Architecture:** One PR, but sequenced into commits that each build green and pass tests — pure header *move* first, then the thrower + dead-code deletion, then the macro→std sweep, then the `enum class` conversions, then full verification. `common.h` stays an umbrella throughout so no downstream include breaks.

**Tech Stack:** C++17, CMake/Makefile build, pytest scenario suite (`python/test_class.py`), `clang-format-22` + `.clang-tidy`.

---

## Background facts (verified 2026-06-19, branch `issue-313-common-h-cleanup`)

- No `.c` translation unit includes `common.h`; `hyrec/*.c` is self-contained. The C (`return _FAILURE_`) macro branch and `*_parallel` macros (0 call sites) are **dead**.
- `ErrorMsg` (8 occ, all in `common.h`), `class_protect_sprintf` (only caller = macros), `class_protect_fprintf`/`class_protect_memcpy` (0 users) are all removable once the macros switch to `ThrowFormatted`.
- `typeof`: only non-`#define` occurrence is a *comment* in `input_module.cpp`.
- `MAX(` ~116, `MIN(` ~93 across the files listed in Task 3; `SIGN(`/`NRSIGN(` are the 6 sites listed in Task 3.
- `class_test` (~429), `class_stop` (~66), `class_open` (~16) call sites are **untouched** — only macro *bodies* change. Their signatures are preserved: `class_test(condition, args, ...)`, `class_stop(args, ...)`, `class_open(pointer, filename, mode)`.
- `class_call`: 5 sites in `transfer_module.cpp`, all `class_call(get_xmin_generic(...), "", "")`.
- Enums to convert: `evolver_type`, `pk_def`, `file_format`, `tca_method`, `rsa_method`, `rsa_idr_method`, `ufa_method`, `ncdmfa_method` (in `common.h` → moving to `precision.h`), and `possible_gauges` (in `source/perturbations.h`). `quadrature_method` is **out of scope**.

## Reference commands

```bash
# Build CLI + smoke test (fast inner loop):
make -j4 class && ./class explanatory.ini

# Build/install the Python wrapper the tests import:
pip install --no-build-isolation .

# Fast scenario suite (no reference compare):
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py

# Full PR reference-comparison suite (needs a `classyref` wheel from master — see Task 5):
cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py

# Format a changed file before committing:
clang-format-22 -i <file>
```

`git add` only the files you changed — **never `git add -A`** in this repo (in-source CMake/Xcode build artifacts).

---

## Task 0: Baseline green

**Files:** none (verification only).

- [ ] **Step 1: Confirm branch and clean build**

Run: `git branch --show-current` → expect `issue-313-common-h-cleanup`.
Run: `make -j4 class && ./class explanatory.ini`
Expected: builds, runs to "Writing output files in output/..." with no error.

- [ ] **Step 2: Confirm fast suite passes on the untouched tree**

Run: `pip install --no-build-isolation . && cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py`
Expected: all scenarios pass. This is the "before" reference.

---

## Task 1: Split `common.h` (pure move, no behavior change)

**Files:**
- Create: `include/constants.h`
- Create: `include/errors.h`
- Create: `include/precision.h`
- Modify: `include/common.h`

Move content verbatim. The only goal of this task is relocation; the thrower, dead-code deletion, macro sweep, and `enum class` all happen later. After it, every `#include "common.h"` still resolves the same symbols.

- [ ] **Step 1: Create `include/constants.h`**

Move the physics constants block (`common.h` lines 63–128: `_PI_`, `_PIHALF_`, `_TWOPI_`, `_SQRT2_`, `_SQRT6_`, `_SQRT_PI_`, `_E_`, `_Mpc_over_m_`, `_Gyr_over_Mpc_`, `_c_`, `_G_`, `_eV_`, `_k_B_`, `_h_P_`) **verbatim**, plus the `Omega0gFromTcmb`/`TcmbFromOmega0g` inline helpers (lines 1034–1047), into:

```cpp
/** @file constants.h Physics constants and conversion factors used across CLASS. */
#ifndef CLASS_CONSTANTS_H
#define CLASS_CONSTANTS_H

#include <cmath>  // pow(), used by the Omega <-> T_cmb helpers below

/* ---- moved verbatim from common.h: _PI_ ... _h_P_ (the physics constants) ---- */
// (paste lines 63–128 here unchanged)

/* ---- moved verbatim from common.h: Omega0gFromTcmb / TcmbFromOmega0g ---- */
// (paste the two inline functions, lines 1034–1047, unchanged)

#endif  // CLASS_CONSTANTS_H
```

- [ ] **Step 2: Create `include/errors.h` (move error machinery verbatim — transform in Task 2)**

Move, **unchanged for now**: the status defines (`_TRUE_`/`_FALSE_`/`_SUCCESS_`/`_FAILURE_`, lines 53–57), `_ERRORMSGSIZE_` + `ErrorMsg` typedef (lines 59–61), every error macro (`class_build_error_string`, `class_call`, `class_call_parallel`, `class_call_message`, `class_test_message`, `class_test`, `class_test_parallel`, `class_stop`, `class_open`, lines 143–227), the entire `#ifdef __cplusplus` C++ override block (lines 229–291), and the `class_protect_*` decls (lines 1049–1056):

```cpp
/** @file errors.h Error-reporting + status codes. */
#ifndef CLASS_ERRORS_H
#define CLASS_ERRORS_H

#include <stdio.h>   // fopen used by class_open
#ifdef __cplusplus
#include <stdexcept>
#include <string>
#endif

// (paste _TRUE_/_FALSE_/_SUCCESS_/_FAILURE_, _ERRORMSGSIZE_, ErrorMsg,
//  the C error macros, the __cplusplus override block, and the
//  class_protect_* extern "C" decls here, verbatim from common.h)

#endif  // CLASS_ERRORS_H
```

- [ ] **Step 3: Create `include/precision.h` (move enums + struct verbatim)**

Move the conversion enums (`evolver_type`, `pk_def`, `file_format`, `tca_method`, `rsa_method`, `rsa_idr_method`, `ufa_method`, `ncdmfa_method`, lines 369–406) and the entire `struct precision` (lines 408–1026) **verbatim** (still plain `enum`, still `int` fields — `enum class` is Task 4):

```cpp
/** @file precision.h Precision parameters and the enums they use. */
#ifndef CLASS_PRECISION_H
#define CLASS_PRECISION_H

#include <float.h>   // DBL_EPSILON
#include <string>

#include "constants.h"  // _TRUE_/_FALSE_ come via errors.h (umbrella) but precision uses _TRUE_
#include "errors.h"     // _TRUE_, _FALSE_

class FileContent;  // forward decl for precision::parse

// (paste the enums lines 369–406 and struct precision lines 408–1026 here, verbatim)

#endif  // CLASS_PRECISION_H
```

- [ ] **Step 4: Reduce `common.h` to a thin umbrella**

In `include/common.h`, delete the blocks moved in Steps 1–3, and add the three includes. What **remains** in `common.h`: the system includes (lines 3–23), the `#define typeof` (removed in Task 3, keep for now), the forward decls + `…ModulePtr` typedefs (lines 25–45), `_VERSION_` (50), `_TRUE_`/`_FALSE_`/`_SUCCESS_`/`_FAILURE_` now come from `errors.h` (delete the duplicates here), the generic defines `_MAX_IT_`/`_QUADRATURE_MAX_`/`_QUADRATURE_MAX_BG_`/`_TOLVAR_`/`_HUGE_`/`_EPSILON_`/`_OUTPUTPRECISION_`/`_COLUMNWIDTH_`/`_DELIMITER_`/`__CLASSDIR__` (78–104), `MIN`/`MAX`/`SIGN`/`NRSIGN` (130–133, removed in Task 3, keep for now), `index_symmetric_matrix` (134–139), `class_define_index` (294–300), the output macros `class_fprintf_*`/`class_store_*` (302–362), and `get_number_of_titles` (1061). Add near the top of the `#ifdef __cplusplus` region:

```cpp
#include "constants.h"
#include "errors.h"
#include "precision.h"
```

Note: `class_fprintf_int`/`class_fprintf_columntitle` use `MAX`/`MIN`, which still live in `common.h` at this stage — fine.

- [ ] **Step 5: Build + smoke + fast suite**

Run: `make -j4 class && ./class explanatory.ini`
Expected: clean build, normal run.
Run: `pip install --no-build-isolation . && cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py`
Expected: all pass (pure move ⇒ identical output).

- [ ] **Step 6: Format + commit**

```bash
clang-format-22 -i include/common.h include/constants.h include/errors.h include/precision.h
git add include/common.h include/constants.h include/errors.h include/precision.h
git commit -m "v4 #313: split common.h into constants.h/errors.h/precision.h (pure move)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Typed thrower + delete dead error infrastructure

**Files:**
- Modify: `tools/common.cpp` (add `ThrowFormatted`, delete `class_protect_*` defs)
- Modify: `include/errors.h` (rewrite to throw-only)
- Modify: `source/transfer_module.cpp:3490,3497,3510,3541,3575` (`class_call` single-arg)

`ThrowFormatted(func, line, prefix, fmt, ...)` keeps `prefix` (the compile-time `#condition`/`#function` text) separate from the user `fmt`, so call sites that pass a runtime (non-literal) format string still work, and the thrown message stays byte-for-byte the same as today.

- [ ] **Step 1: Add `ThrowFormatted` to `tools/common.cpp`; remove the `class_protect_*` definitions**

Open `tools/common.cpp`. Delete the bodies of `class_protect_sprintf`, `class_protect_fprintf`, `class_protect_memcpy`. Add:

```cpp
#include <cstdarg>
#include <cstdio>
#include <stdexcept>
#include <string>

[[noreturn]] void ThrowFormatted(const char* func, int line,
                                 const char* prefix, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  va_list args_copy;
  va_copy(args_copy, args);
  const int n = std::vsnprintf(nullptr, 0, fmt, args_copy);
  va_end(args_copy);

  std::string body;
  if (n > 0) {
    body.resize(static_cast<std::size_t>(n));
    std::vsnprintf(body.data(), static_cast<std::size_t>(n) + 1, fmt, args);
  }
  va_end(args);

  // Preserves the historical "<func>(L:<line>) :<message>" format.
  const std::string message =
      std::string(func) + "(L:" + std::to_string(line) + ") :" + prefix + body;
  throw std::runtime_error(message);
}
```

- [ ] **Step 2: Rewrite `include/errors.h` to throw-only**

Replace everything between the include guard with:

```cpp
/** @file errors.h Error-reporting (C++ exceptions) + status codes. */
#ifndef CLASS_ERRORS_H
#define CLASS_ERRORS_H

#include <cstdio>   // fopen, nullptr-compare in class_open

#define _TRUE_ 1
#define _FALSE_ 0
#define _SUCCESS_ 0
#define _FAILURE_ 1

// CLASSpp is C++-only (hyrec is the sole C code and never includes this header),
// so errors.h is a plain C++ header: no `#ifdef __cplusplus`, no `extern "C"`.

[[noreturn]] void ThrowFormatted(const char* func, int line,
                                 const char* prefix, const char* fmt, ...);

/* All error macros throw std::runtime_error from the error origin. */

#define class_call(function)                                          \
  {                                                                   \
    if ((function) == _FAILURE_) {                                    \
      ThrowFormatted(__func__, __LINE__, "error in " #function, "");  \
    }                                                                 \
  }

#define class_test(condition, args, ...)                                          \
  {                                                                               \
    if (condition) {                                                              \
      ThrowFormatted(__func__, __LINE__,                                          \
                     "condition (" #condition ") is true; ", args, ##__VA_ARGS__);\
    }                                                                             \
  }

#define class_stop(args, ...)                                              \
  {                                                                        \
    ThrowFormatted(__func__, __LINE__, "error; ", args, ##__VA_ARGS__);    \
  }

#define class_open(pointer, filename, mode)                                \
  {                                                                        \
    pointer = fopen(filename, mode);                                       \
    if (pointer == nullptr) {                                              \
      ThrowFormatted(__func__, __LINE__, "",                              \
                     "could not open %s with name %s and mode %s",         \
                     #pointer, filename, #mode);                           \
    }                                                                      \
  }

#endif  // CLASS_ERRORS_H
```

Deleted here: `_ERRORMSGSIZE_`, `ErrorMsg`, `class_build_error_string`, `class_call_message`, `class_test_message`, `class_call_parallel`, `class_test_parallel`, the C `return _FAILURE_` branch, the `class_protect_*` decls, **and the whole `#ifdef __cplusplus` / `extern "C"` C-compat scaffolding** — CLASSpp is C++-only (only `hyrec` is C, and it never includes this header), so `errors.h` is a plain C++ header.

- [ ] **Step 3: Simplify the 5 `class_call` sites**

In `source/transfer_module.cpp`, change each `class_call(get_xmin_generic(...), "", "")` to drop the trailing `, "", ""`. Example (`:3490`):

```cpp
// before:
class_call(get_xmin_generic(sgnK, lvec[0], nu, xtol, phiminabs, &x_nonzero, &fevals), "", "");
// after:
class_call(get_xmin_generic(sgnK, lvec[0], nu, xtol, phiminabs, &x_nonzero, &fevals));
```

Apply at lines 3490, 3497, 3510, 3541, 3575 (the multi-line ones drop the `, "", ""` on their final argument line).

- [ ] **Step 4: Guard — confirm the dead symbols are gone**

Run:
```bash
grep -rn 'ErrorMsg\|class_protect_\|class_build_error_string\|_parallel\|_ERRORMSGSIZE_' \
  source include tools species --include='*.cpp' --include='*.h' | grep -v 'in parallel\|OpenMP'
```
Expected: no hits referencing the removed symbols (only unrelated prose, if any).

- [ ] **Step 5: Build + smoke + fast suite**

Run: `make -j4 class && ./class explanatory.ini`
Expected: clean build and run. (Trigger an error path if convenient, e.g. a bad `.ini`, and confirm the message still reads `…(L:NN) :condition (…) is true; …`.)
Run: `pip install --no-build-isolation . && cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py`
Expected: all pass.

- [ ] **Step 6: Format + commit**

```bash
clang-format-22 -i include/errors.h tools/common.cpp source/transfer_module.cpp
git add include/errors.h tools/common.cpp source/transfer_module.cpp
git commit -m "v4 #313: typed ThrowFormatted; delete dead C error path + ErrorMsg/protect_*

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Replace `MIN/MAX/SIGN/NRSIGN/typeof` with the standard library

**Files (call sites):**
- `include/common.h` (delete the 5 `#define`s)
- `source/`: `background_module.cpp`, `input_module.cpp`, `lensing_module.cpp`, `nonlinear_module.cpp`, `output_module.cpp`, `perturbations_module.cpp`, `primordial_module.cpp`, `spectra_module.cpp`, `thermodynamics_module.cpp`, `transfer_module.cpp`
- `species/ncdm_base_species.cpp`
- `tools/`: `arrays.cpp`, `dei_rkck.cpp`, `evolver_ndf15.cpp`, `evolver_rkdp45.cpp`, `hyperspherical.cpp`, `quadrature.cpp`, `sparse.cpp`
- `tools/hermite3_interpolation_csource.h`, `tools/hermite4_interpolation_csource.h`, `tools/hermite6_interpolation_csource.h` (the compiled copies, `#include`d by `hyperspherical.cpp`)
- `include/hermite3_interpolation_csource.h`, `include/hermite4_interpolation_csource.h`, `include/hermite6_interpolation_csource.h` (unused duplicates — convert for consistency; deleting them is a possible follow-up, out of scope here)

**Transformation rules:**
- `MAX(a, b)` → `std::max(a, b)`; `MIN(a, b)` → `std::min(a, b)`.
  - **Type gotcha:** `std::max`/`std::min` require both args the same type. Where the macro mixed types (e.g. `MAX(0, _COLUMNWIDTH_ - ...)` with a `double`, or `MAX(0, <double>)`), type the literal: `std::max(0.0, x)` / `std::max<double>(0, x)`. The compiler flags every mismatch — fix until it builds.
- `SIGN(a)` → `std::copysign(1.0, a)`.
- `NRSIGN(a, b)` → `std::copysign(a, b)`.
- `typeof` → delete the `#define`; no call sites.

- [ ] **Step 1: Delete the macros from `common.h`**

Remove lines 130–133 (`MIN`/`MAX`/`SIGN`/`NRSIGN`) and line 16 (`#define typeof(x) ...`). Keep `index_symmetric_matrix`.

- [ ] **Step 2: Convert `SIGN`/`NRSIGN` (the 6 sites)**

`source/nonlinear_module.cpp:1654-1655`:
```cpp
double cosine_correlation = primordial_pk[index_ic1_ic2] *
    std::copysign(1.0, source_ic1) * std::copysign(1.0, source_ic2);
```
`tools/hyperspherical.cpp:1351`:
```cpp
} while (std::copysign(1.0, Fnew) == std::copysign(1.0, Fold));
```
`tools/hyperspherical.cpp:1437,1443,1447`:
```cpp
if (std::copysign(fm, fnew) != fm) { ... }
else if (std::copysign(fl, fnew) != fl) { ... }
else if (std::copysign(fh, fnew) != fh) { ... }
```
Ensure `#include <cmath>` is present in both files (add if missing).

- [ ] **Step 3: Convert `MAX`/`MIN` across the listed files**

Mechanical. A safe starting sweep per file (then hand-fix the type-mismatch sites the compiler reports):
```bash
# Run per file, NOT repo-wide, and review the diff:
sed -i '' -E 's/\bMAX\(/std::max(/g; s/\bMIN\(/std::min(/g' <file>
```
Add `#include <algorithm>` to any file that gains `std::max`/`std::min` and doesn't already include it. Representative mixed-type fix in `common.h`'s output macros — `class_fprintf_int` / `class_fprintf_columntitle` reference `MAX(0, _COLUMNWIDTH_ - _OUTPUTPRECISION_ - 5)` etc.; these are all `int` arithmetic, so `std::max(0, ...)` is fine, but `common.h` must then `#include <algorithm>`.

- [ ] **Step 4: Reword the `typeof` comment**

`source/input_module.cpp:2193` — change `// ... matches the old class_read_int (typeof) cast.` to drop the dead `typeof` reference, e.g. `// ... matches the old class_read_int narrowing cast.`

- [ ] **Step 5: Build until clean**

Run: `make -j4 class`
Expected: iterate on `std::max`/`std::min` type-mismatch errors until it builds, then `./class explanatory.ini` runs.

- [ ] **Step 6: Numeric verification (this is the only output-risk task)**

Run: `pip install --no-build-isolation . && cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py`
Expected: all pass. The `copysign` swaps differ from the old `SIGN`/`NRSIGN` only at exact `0.`/`-0.` inputs; if any scenario shifts, it is expected and acceptable (see Task 5).

- [ ] **Step 7: Format + commit**

```bash
clang-format-22 -i <every changed file>
git add include/common.h source/*.cpp species/ncdm_base_species.cpp tools/*.cpp \
        tools/hermite3_interpolation_csource.h tools/hermite4_interpolation_csource.h \
        tools/hermite6_interpolation_csource.h include/hermite3_interpolation_csource.h \
        include/hermite4_interpolation_csource.h include/hermite6_interpolation_csource.h
git commit -m "v4 #313: replace MIN/MAX with std::min/max, SIGN/NRSIGN with std::copysign; drop typeof

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: `enum` → `enum class`

Each sub-task is **atomic**: changing `enum X` to `enum class X` forces every value reference to a scoped name in the same commit, or it won't compile. Order: smallest blast radius first. After each, `make -j4 class && ./class explanatory.ini`, then `TEST_LEVEL=1 COMPARE_OUTPUT_REF=0` fast suite, then format + commit.

**Two categories (decided 2026-06-19 with the user, to keep the fragile Cython wrapper almost entirely untouched):**
- **Typed-field enums — `possible_gauges`, `pk_def`, `evolver_type`:** fields are already the enum type and NONE are Cython-exposed. Convert fully: `enum class` + scoped values + the field keeps the enum type with a scoped default. Where an integer input is read into such a field, read into a temp `int` and `static_cast<evolver_type>(tmp)`. Zero wrapper impact.
- **Int-stored approximation enums — `tca_method`, `rsa_method`, `rsa_idr_method`, `ufa_method`, `ncdmfa_method`:** their `precision` fields are `int` (user-facing `.ini` integer knobs) AND declared `int` in `cclassy.pxd`. **Keep the fields `int`.** Make only the enum *definitions* `enum class` (to scope the value names); at C++ compare/assign sites use `static_cast<int>(tca_method::first_order_CLASS)`. `classy.pyx` never touches these fields → pxd stays valid → zero wrapper impact. The `.ini` read path is unchanged (reads `int` into the `int` field).
- **Cython-exposed enum — `file_format`:** the ONE that needs wrapper work (a real Cython type in `perturb_output_*` signatures). See Step 2.

**Use `static_cast`, never C-style `(int)` casts, in all the new enum↔int boundary code.**

- [ ] **Step 1: `evolver_type` + `pk_def`** (typed fields, no Cython exposure)

`precision.h`: `enum class evolver_type { rk, ndf15, rkdp45 };`, `enum class pk_def { delta_m_squared, delta_tot_squared, delta_bc_squared, delta_tot_from_poisson_squared };`. Field `evolver_type evolver = evolver_type::ndf15;`. Scope all `ndf15`/`rk`/`rkdp45` and `delta_*_squared` value sites (evolvers, `nonlinear`, `spectra`). If `input_module` reads an integer `evolver` key, read into a temp `int` and `ppr->evolver = static_cast<evolver_type>(tmp);`. Build, fast-test, commit.

- [ ] **Step 2: `file_format`** (Cython-exposed — needs the wrapper change)

`precision.h`: `enum class file_format { class_format, camb_format };`. Scope all `class_format`/`camb_format` sites (~40, `species/*` column writers + `output`/`perturbations`). Then the **wrapper trio** (the Cython layer lives at the REPO ROOT: `classy.pyx`, `cclassy.pxd` [auto-generated + gitignored], `generate_wrapper.py`):
  - `generate_wrapper.py`: its enum scanner matches `'enum ' + name + ' '` (won't match `enum class file_format`) and emits `ctypedef enum file_format : …`. Update it to also detect `enum class <name>` and emit the Cython scoped form (`cdef enum class file_format:` with `class_format` / `camb_format` indented under it). Cython supports C++ scoped enums.
  - `classy.pyx` (~lines 1448-1465): `outf = camb_format` / `outf = class_format` → `outf = file_format.camb_format` / `file_format.class_format`. The `file_format outf` declaration and the `perturb_output_*` calls stay.
  - Rebuild with `pip install --no-build-isolation .` (regenerates `cclassy.pxd`); confirm `get_transfer(output_format='camb')` still works (`python/test_transfer_columns.py` + the fast suite exercise it).
  Build, fast-test, commit.

- [ ] **Step 3: approximation methods (`tca_method`, `rsa_method`, `rsa_idr_method`, `ufa_method`, `ncdmfa_method`) — enum-class defs, fields stay `int`**

`precision.h`: make all five `enum class`. **Do NOT retype the `precision` fields** (`tight_coupling_approximation`, `radiation_streaming_approximation`, `idr_streaming_approximation`, `ur_fluid_approximation`, `ncdm_fluid_approximation` stay `int`). Update defaults + compare sites to scoped, `static_cast`-ed values:
- default e.g. `int tight_coupling_approximation = static_cast<int>(tca_method::compromise_CLASS);`
- comparisons e.g. `input_module.cpp:2227-2228` `== (int) first_order_CLASS` → `== static_cast<int>(tca_method::first_order_CLASS)`; same for the `perturbations_module.cpp` switch/compare sites.
The `.ini` read path (`ReadPrecision`, `input_module.cpp:2520,2545,2555,2557`) is UNCHANGED (reads `int` into the `int` field). No `cclassy.pxd`/`classy.pyx` change. Build, fast-test, commit.

- [ ] **Step 4: `possible_gauges` (biggest — typed field, no Cython exposure)**

`source/perturbations.h`: `enum class possible_gauges { newtonian, synchronous };`, field `possible_gauges gauge = possible_gauges::synchronous;`. Scope every code reference (`== newtonian` → `== possible_gauges::newtonian`), concentrated in `perturbations_module.cpp`/`.h` and `species/perturb_source_context.h` plus any species using gauge. (`newtonian`/`synchronous` in comments/strings stay.) `gauge` is not Cython-exposed. Build until clean, fast-test, commit.

---

## Task 5: Full verification + baseline reconciliation

**Files:** none, unless baselines need regenerating.

- [ ] **Step 1: Clean rebuild + smoke**

Run: `make -j4 class && ./class explanatory.ini`. Expected: clean.

- [ ] **Step 2: Build the `classyref` reference wheel from the branch base (master)**

The full suite compares against a second extension named `classyref` built from master. Per `.github/workflows/test_on_pull_request.yml`, build it from a clean master checkout/worktree:
```bash
git worktree add ../classref-master master
cd ../classref-master
sed -i.bak "s/Extension('classy'/Extension('classyref'/g" setup.py
pip install . --config-settings=cmake.define.CLASS_PYTHON_MODULE_NAME=classyref
cd -
git worktree remove ../classref-master
```
(If a fresh `classyref` from the current master is already installed, skip.)

- [ ] **Step 3: Install this branch's `classy` and run the full reference suite**

Run:
```bash
pip install --no-build-isolation .
cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py
```
Expected: all scenarios pass within tolerance. This refactor is numerically inert except the `copysign` swaps, so a clean pass is the expected outcome.

- [ ] **Step 4: Reconcile any diffs**

If a scenario fails only on a tiny `copysign`-induced shift (handle `Cl^TE` zero-crossings; never a blind max-rel-diff), confirm it's that and accept it. If the project stores goldens that need regenerating, regenerate and commit them in a clearly-labelled commit; otherwise note that `classyref`=master self-heals after merge.

- [ ] **Step 5: Final format sweep + verify `common.h` shrank**

Run: `clang-format-22 --dry-run --Werror <all changed files>` → expect no diff.
Run: `wc -l include/common.h` → expect it down from 1064 to roughly the umbrella size (~150–200 lines).

- [ ] **Step 6: Push + open PR (only when the user asks)**

Open a PR referencing #313 summarizing: header split, `ThrowFormatted` + dead-C-path/`ErrorMsg`/`protect_*` deletion, `std::min/max`/`copysign` sweep, `enum class` conversions, and the verification result.

---

## Self-review notes

- **Spec coverage:** header split (Task 1) ✓; typed thrower + dead-code deletion (Task 2) ✓; `MIN/MAX/SIGN/NRSIGN`/`typeof` (Task 3) ✓; `enum class` incl. gauges + approximation methods + `.ini` integer interface preserved (Task 4) ✓; verification within tolerance, not bit-identical (Task 5) ✓; `quadrature_method` and `index_symmetric_matrix` left alone ✓.
- **Thrower signature:** `ThrowFormatted(func, line, prefix, fmt, ...)` is used identically in the decl (errors.h Step 2) and def (common.cpp Step 1).
- **Build manifests:** no new `.cpp` translation unit (thrower lives in `tools/common.cpp`), so Makefile/CMake/setup.py/Xcode source lists are untouched; new headers resolve via the include path.
