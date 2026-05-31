# Code Cleanup Pass: Tier-1 Modernization + NULL→nullptr Sweep

## Motivation

The post-species-refactor codebase still carries small bits of C-era code and
historical dead comments that no longer pay rent. Targeted, low-risk cleanup
reduces visual noise, removes a couple of one-off C helpers, and aligns the
codebase with C++17 idioms it already otherwise uses.

This spec covers a single bundled PR with four logically-separate commits.
Heavier refactors (slab allocator in `thermodynamics_module`, `const_cast` wall
around the `array_*` API) are explicitly out of scope and tracked for later.

## Scope

### In scope (four commits in one PR)

1. **Replace `compare_doubles` + `qsort` with `std::sort`.**
   - Remove `InputModule::compare_doubles` (`source/input_module.cpp:3240`,
     decl `source/input_module.h:87`).
   - Rewrite the call site (`source/input_module.cpp:2663`) as
     `std::sort(ppt->k_output_values, ppt->k_output_values + ppt->k_output_values_num);`
     (natural double ordering; no comparator needed).
   - Add `#include <algorithm>` to `input_module.cpp` if not already pulled in.

2. **Replace `file_exists` wrapper with `std::filesystem::exists`.**
   - Remove `InputModule::file_exists` (`source/input_module.cpp:3231`,
     decl `source/input_module.h:88`).
   - Rewrite the four call sites in `source/input_module.cpp:132–141` to
     `std::filesystem::exists(tmp_file, ec)` using the non-throwing
     `(path, error_code&)` overload (returns `bool`; the `== _TRUE_`
     comparisons go away). The throwing single-arg overload is NOT used —
     see "Error Handling" below for why.
   - Add `#include <filesystem>` to `input_module.cpp`.
   - C++17 is the project standard (`Makefile:37`, `setup.py:148`), so this
     is available everywhere we compile.

3. **Remove dead commented-out code blocks.**
   - `source/input_module.cpp:1331–1334`: stale `//free(ppt->alpha_idm_dr)` /
     `//free(ppt->beta_idr)` comments. The pointer fields themselves still
     exist (`source/perturbations.h:200–201`), but they are now non-owning
     views into `alpha_idm_dr_storage` / `beta_idr_storage` vectors — so
     calling `free()` on them was already wrong and the comments are dead.
   - `source/nonlinear_module.cpp:2293–2308`: the entire `/* ... */`-commented
     `class_test_except(...)` block with its now-irrelevant `free(pvecback);
     free(integrand_array)` cleanup.
   - `source/thermodynamics_module.cpp:4279–4335`: roughly a dozen
     `//class_store_columntitle(...)` and `//class_store_double(...)` lines
     for thermodynamics derivatives that are not exposed. Keep only the active
     calls; drop the commented siblings.

4. **`NULL` → `nullptr` sweep in our C++ code.**
   - Scope: `*.cpp` under `source/`, `tools/` and `*.h` under `source/`,
     `species/`, `tools/`, `include/`, `main/`.
   - **Excluded:**
     - `hyrec/` (third-party).
     - `include/hermite{3,4,6}_interpolation_csource.h` and the matching
       `tools/hermite{3,4,6}_interpolation_csource.h` — these are C source
       fragments spliced via `#include`; leaving them as `NULL` keeps them
       valid as C.
   - Expected ~170 substitutions: 145 in `source/*.cpp`, 25 in `tools/*.cpp`,
     2 in `include/common.h`. Species is already `NULL`-free.
   - Pure token replacement; no semantic change.

### Out of scope (flagged for future passes)

- The slab allocator with manual byte-offset arithmetic in
  `source/thermodynamics_module.cpp:3220–3260` (rate-table setup). Real
  refactor candidate but touches interop with hyrec.
- The wall of `const_cast<double*>(vec.data())` in `background_module`,
  `lensing_module`, `output_module`. Symptom of `tools/arrays.c` interfaces
  taking non-`const` pointers; fix at the API layer first, then strip casts.
- TODO/FIXME triage (~41 occurrences).
- `header-only-once` includes audit.

## Architecture / Files Touched

| Commit | Files | LOC change |
|---|---|---|
| 1 | `source/input_module.cpp`, `source/input_module.h` | −10 / +2 |
| 2 | `source/input_module.cpp`, `source/input_module.h` | −10 / +5 |
| 3 | `source/input_module.cpp`, `source/nonlinear_module.cpp`, `source/thermodynamics_module.cpp` | −30 / 0 |
| 4 | many .cpp/.h under source/, tools/, include/ | ~170 / ~170 |

## Error Handling

None of the four commits change error paths.

- `std::sort` is `noexcept`-equivalent here (sorts a `double*` range).
- `std::filesystem::exists(path)` can throw `filesystem_error` on
  permission/IO failures. The current `file_exists` silently returns
  `_FALSE_` in that case. We use the **`(std::error_code&)` overload** so
  behaviour matches: `std::filesystem::exists(tmp_file, ec)` returns `false`
  on error and never throws.

## Testing / Verification

For each commit:

1. Build cleanly: `make class` and `python setup.py build_ext --inplace`.
2. Run any existing test suite (e.g. `test/scenarios/`).
3. For numerical confidence, pick one moderate scenario (`gauge_lcdm.ini`
   is the lightest), run before and after the cleanup commit, and diff with
   `test/scenarios/compare_tol.py` at its default `RTOL = 1e-3` with the
   zero-crossing-aware column-peak floor. Per project convention, **do not
   require bit-identical output** — verify within ~0.1% tolerance.
4. None of these changes should perturb numerics:
   - `std::sort` is order-equivalent to `qsort` on `double` with strict ordering.
   - `std::filesystem::exists` is semantically identical to the `fopen`-probe.
   - Removed comments are non-executing.
   - `NULL`/`nullptr` are interchangeable for pointer comparison and
     assignment.

   If `compare_tol.py` reports any FAIL, treat it as a real regression and
   investigate before merging.

## Commit Plan (for the PR)

```
1. Replace compare_doubles + qsort with std::sort in input_module
2. Replace file_exists with std::filesystem::exists in input_module
3. Remove dead commented-out blocks (input, nonlinear, thermo)
4. NULL → nullptr in our C++ code (excludes hyrec, *_csource.h)
```

Each commit is independently buildable and revertable.

## Rollback

Each commit is independently revertable; commits 1–3 touch a single file each
or three discrete blocks; commit 4 is mechanical and revertable in one
operation.

## Risks

- **Commit 2 (`<filesystem>`):** introduces a new standard-library dependency.
  Already supported by the C++17 toolchain the project requires, but adds the
  `-lstdc++fs` link flag on some older Linux toolchains. **Mitigation:**
  verify the existing `Makefile` link line covers libstdc++ filesystem on the
  CI / dev platforms; on modern GCC ≥ 9 and clang ≥ 9 with libc++ it is
  bundled and needs no extra flag. If a flag is needed, add it in commit 2.
- **Commit 4 (NULL sweep):** wide diff. Risk that a `NULL` appears inside a
  macro that the C preprocessor stringifies (`#NULL` becoming `"NULL"`),
  which would change behaviour. **Mitigation:** grep for `#NULL` and
  `"NULL"` patterns first; review the diff for any non-pointer uses (e.g.
  format strings, log messages).
