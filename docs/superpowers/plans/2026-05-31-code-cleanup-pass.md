# Code Cleanup Pass Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace a few residual C-isms with C++17 idioms, delete commented-out dead code, and sweep `NULL` → `nullptr` across our C++ sources — as one PR with four logical commits.

**Architecture:** Pure refactor. No public API, behaviour, or numerical output changes. Each commit is independently buildable and revertable. Numerical regression check uses the project's existing `test/scenarios/compare_tol.py` harness against `gauge_lcdm.ini` at the default 1e-3 tolerance.

**Tech Stack:** C++17 (`-std=c++17` in `Makefile:37` and `setup.py:148`), `std::sort`, `std::filesystem`, GNU `sed` / BSD `sed` for the mechanical sweep.

**Spec:** `docs/superpowers/specs/2026-05-31-code-cleanup-pass-design.md`

---

## Task 0: Baseline build + reference output

**Files:** none (read-only baseline capture)

- [ ] **Step 1: Confirm clean working tree**

```bash
git status --short
```

Expected: empty or only untracked files unrelated to the cleanup. If there are unrelated unstaged changes, stop and ask the user.

- [ ] **Step 2: Build the project from master state**

```bash
make class -j4
```

Expected: builds cleanly to `./class` (or whatever `Makefile` names it). If the build fails on master without any of our changes, stop — the environment is broken and the rest of the plan is meaningless until that's fixed.

- [ ] **Step 3: Generate baseline numerical reference**

```bash
mkdir -p /tmp/cleanup-baseline /tmp/cleanup-new
./class test/scenarios/gauge_lcdm.ini
# CLASS writes outputs to wherever the .ini's `root =` points;
# the gauge_lcdm.ini in test/scenarios uses output/ relative to CWD.
# Move outputs to /tmp/cleanup-baseline so they survive subsequent runs.
mv output/gauge_lcdm*.dat /tmp/cleanup-baseline/
ls /tmp/cleanup-baseline/
```

Expected: a handful of `.dat` files (e.g. `gauge_lcdm00_cl.dat`, `gauge_lcdm00_pk.dat`, etc.). If `output/` doesn't exist, `mkdir -p output` first. If `gauge_lcdm.ini` is missing a `root =` line, briefly inspect it and adjust to write to `/tmp/cleanup-baseline/gauge_lcdm`.

- [ ] **Step 4: No commit (baseline-only setup)**

Skip — this task only captures a reference.

---

## Task 1: Replace `compare_doubles` + `qsort` with `std::sort`

**Files:**
- Modify: `source/input_module.cpp:2663` (call site)
- Modify: `source/input_module.cpp:3240–3248` (helper definition)
- Modify: `source/input_module.h:87` (helper declaration)

- [ ] **Step 1: Replace the `qsort` call**

Find this block at `source/input_module.cpp:2662–2663`:

```cpp
    /* Sort the k_array using qsort */
    qsort(ppt->k_output_values, ppt->k_output_values_num, sizeof(double), compare_doubles);
```

Replace with:

```cpp
    /* Sort k_output_values ascending */
    std::sort(ppt->k_output_values, ppt->k_output_values + ppt->k_output_values_num);
```

(`#include <algorithm>` is already present at `source/input_module.cpp:8`.)

- [ ] **Step 2: Remove the helper definition**

Delete `source/input_module.cpp:3240–3248`:

```cpp
int InputModule::compare_doubles(const void* a, const void* b) {
  double* x = (double*) a;
  double* y = (double*) b;
  if (*x < *y)
    return -1;
  else if (*x > *y)
    return 1;
  return 0;
}
```

If there's a blank line above or below that goes orphan, tidy up so the file ends with a single trailing newline.

- [ ] **Step 3: Remove the helper declaration**

Delete this line from `source/input_module.h` (at line 87 of the pre-edit file):

```cpp
  static int compare_doubles(const void* a, const void* b);
```

- [ ] **Step 4: Build**

```bash
make class -j4
```

Expected: clean build, no warnings about `compare_doubles`.

- [ ] **Step 5: Numerical regression check**

```bash
rm -rf /tmp/cleanup-new && mkdir -p /tmp/cleanup-new output
./class test/scenarios/gauge_lcdm.ini
mv output/gauge_lcdm*.dat /tmp/cleanup-new/
python test/scenarios/compare_tol.py /tmp/cleanup-baseline /tmp/cleanup-new
```

Expected: every file reported `OK`. If anything reports `FAIL`, do **not** commit — investigate. `std::sort` with `<` on a `double` range produces the same ordering as `qsort` with the deleted comparator, so a FAIL here indicates an unrelated issue (e.g. you accidentally touched a different line).

- [ ] **Step 6: Commit**

```bash
git add source/input_module.cpp source/input_module.h
git commit -m "$(cat <<'EOF'
refactor: replace compare_doubles + qsort with std::sort

Drop the static int(*)(const void*, const void*) comparator and the
qsort callsite in InputModule. std::sort on a (double*, double*+n) range
is order-equivalent and avoids a C-style helper that existed only to
sort k_output_values once.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Replace `file_exists` with `std::filesystem::exists`

**Files:**
- Modify: `source/input_module.cpp` (add `#include <filesystem>`; rewrite 4 call sites; delete helper)
- Modify: `source/input_module.h` (delete decl)

- [ ] **Step 1: Add `<filesystem>` include**

In `source/input_module.cpp` near the top, after `#include <algorithm>` (line 8), add:

```cpp
#include <filesystem>
```

Keep includes alphabetically grouped if the file already follows that style; otherwise just place it next to the other standard-library include.

- [ ] **Step 2: Rewrite the four call sites**

Find this block at `source/input_module.cpp:130–144`:

```cpp
      for (filenum = 0; filenum < 100; filenum++) {
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_cl.dat", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_pk.dat", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_tk.dat", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_parameters.ini", inifilename, filenum);
        if (file_exists(tmp_file) == _TRUE_)
          continue;
        break;
      }
```

Replace with:

```cpp
      for (filenum = 0; filenum < 100; filenum++) {
        std::error_code ec;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_cl.dat", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_pk.dat", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_tk.dat", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        snprintf(tmp_file, tmp_file_size, "output/%s%02d_parameters.ini", inifilename, filenum);
        if (std::filesystem::exists(tmp_file, ec))
          continue;
        break;
      }
```

Rationale for the `std::error_code` overload: the original `file_exists` silently returned `_FALSE_` on any failure to open. The throwing `std::filesystem::exists(path)` overload would change that to throw on permission/IO errors; the `(path, ec)` overload preserves the silent-false-on-error behaviour.

- [ ] **Step 3: Delete the helper definition**

At `source/input_module.cpp:3231–3238` (post-Task-1 line numbers will have shifted; search for `int InputModule::file_exists`):

```cpp
int InputModule::file_exists(const char* fname) {
  FILE* file = fopen(fname, "r");
  if (file != NULL) {
    fclose(file);
    return _TRUE_;
  }
  return _FALSE_;
}
```

Delete the whole function.

- [ ] **Step 4: Delete the helper declaration**

From `source/input_module.h`, remove:

```cpp
  static int file_exists(const char* fname);
```

- [ ] **Step 5: Build**

```bash
make class -j4
```

Expected: clean build. If the linker complains about `std::filesystem`, add `-lstdc++fs` to `LDFLAGS` in `Makefile`. On modern macOS (Apple LLVM ≥ 11) and GCC ≥ 9 this is not needed.

- [ ] **Step 6: Numerical regression check**

```bash
rm -rf /tmp/cleanup-new && mkdir -p /tmp/cleanup-new output
./class test/scenarios/gauge_lcdm.ini
mv output/gauge_lcdm*.dat /tmp/cleanup-new/
python test/scenarios/compare_tol.py /tmp/cleanup-baseline /tmp/cleanup-new
```

Expected: every file `OK`. (`file_exists` is only used to pick an unused output filename prefix; with a clean `output/` directory at start of each run, the picked prefix is the same as before.)

- [ ] **Step 7: Commit**

```bash
git add source/input_module.cpp source/input_module.h
git commit -m "$(cat <<'EOF'
refactor: replace file_exists with std::filesystem::exists

Drop the 8-line fopen-probe wrapper in favour of the C++17 standard.
Use the (path, error_code&) overload so the silent-false-on-error
behaviour of the original is preserved (the throwing overload would
have been a subtle behaviour change).

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Remove dead commented-out code blocks

**Files:**
- Modify: `source/input_module.cpp` (5-line block near old line 1331)
- Modify: `source/nonlinear_module.cpp` (block-comment near old line 2293)
- Modify: `source/thermodynamics_module.cpp` (~12 stale `//class_store_*` lines around old lines 4279–4335)

- [ ] **Step 1: input_module.cpp — drop stale `//free()` comment block**

Find this in `source/input_module.cpp` (around line 1330, search for `// The following lines make sure that if perturbations are not computed`):

```cpp
  /* The following lines make sure that if perturbations are not computed, IDR parameters are still freed */
  if (ppt->has_perturbations == _FALSE_) {
    //free(ppt->alpha_idm_dr);
    //free(ppt->beta_idr);
  }
```

The `alpha_idm_dr` / `beta_idr` pointer fields still exist
(`source/perturbations.h:200–201`), but they are now non-owning views into
`alpha_idm_dr_storage` / `beta_idr_storage` vectors; the raw `free()`
ownership pattern is gone, so the commented-out `free()` calls are dead.
Delete all 5 lines.

- [ ] **Step 2: nonlinear_module.cpp — drop block-commented `class_test_except`**

Search `source/nonlinear_module.cpp` for `special error handling here (using class_test_except` (around line 2293–2308):

```cpp
  /* the following error should not stop the code: it will arrive
     inevitably at some large redshift, and then the code should not
     stop, but just give up computing P_NL(k,z). This is why we have a
     special error handling here (using class_test_except and free()
     commands to avoid memory leaks, and calling this whole function
     not through a class_call) */

  /*
    class_test_except(sigma < 1.,
    error_message_,
    free(pvecback);free(integrand_array),
    "Your k_max=%g 1/Mpc is too small for Halofit to find the non-linearity scale z_nl at z=%g. Increase input parameter P_k_max_h/Mpc or P_k_max_1/Mpc",
    k_[k_size_-1],
    pba->a_today/pvecback[background_module_->index_bg_a_]-1.);
  */
```

The descriptive comment above also references the dead approach. Replace **both** blocks with a single, terser comment explaining what the code below still does:

```cpp
  /* Don't abort if sigma < 1: at high z we may legitimately fail to
     find the non-linearity scale, and the caller (computed via
     class_test/class_call) must be allowed to give up on P_NL(k,z)
     without halting the program. */
```

(The `if (sigma < 1.) { ... }` block that follows the deleted comments is unchanged.)

- [ ] **Step 3: thermodynamics_module.cpp — drop stale `//class_store_*` lines**

In `source/thermodynamics_module.cpp` around lines 4279–4335, two clusters of commented-out store calls (titles and values for derivatives that are not produced). Search and delete each of these lines (they are all entire-line comments — easy to spot):

In the "thermodynamics_titles" function (search `class_store_columntitle(titles, "z"`):

```cpp
  //class_store_columntitle(titles,"kappa''",_TRUE_);
  //class_store_columntitle(titles,"kappa'''",_TRUE_);
  //class_store_columntitle(titles,"g'",_TRUE_);
  //class_store_columntitle(titles,"g''",_TRUE_);
  //class_store_columntitle(titles,"max. rate",_TRUE_);
  //class_store_columntitle(titles,"ddmu_idm_dr",_TRUE_);
  //class_store_columntitle(titles,"dddmu_idm_dr",_TRUE_);
```

In `thermodynamics_output_data` (search `class_store_double(dataptr, pvecthermo[index_th_xe_]`):

```cpp
  //pth->number_of_thermodynamics_titles = get_number_of_titles(pth->thermodynamics_titles);
  //pth->size_thermodynamics_data = pth->number_of_thermodynamics_titles*tt_size_;
  //class_store_double(dataptr, pvecthermo[index_th_ddkappa_],_TRUE_, storeidx);
  //class_store_double(dataptr, pvecthermo[index_th_dddkappa_],_TRUE_, storeidx);
  //class_store_double(dataptr, pvecthermo[index_th_dg_],_TRUE_, storeidx);
  //class_store_double(dataptr, pvecthermo[index_th_ddg_],_TRUE_, storeidx);
  //class_store_double(dataptr, pvecthermo[index_th_rate_],_TRUE_, storeidx);
```

Delete every entire-line comment matching `//\s*class_store_*` and the two `//pth->...` lines listed above. Do **not** delete the active (non-commented) `class_store_*` calls.

- [ ] **Step 4: Build**

```bash
make class -j4
```

Expected: clean build with no warnings.

- [ ] **Step 5: Numerical regression check**

```bash
rm -rf /tmp/cleanup-new && mkdir -p /tmp/cleanup-new output
./class test/scenarios/gauge_lcdm.ini
mv output/gauge_lcdm*.dat /tmp/cleanup-new/
python test/scenarios/compare_tol.py /tmp/cleanup-baseline /tmp/cleanup-new
```

Expected: every file `OK`. Deleted lines are non-executing comments; numerics must not move.

- [ ] **Step 6: Commit**

```bash
git add source/input_module.cpp source/nonlinear_module.cpp source/thermodynamics_module.cpp
git commit -m "$(cat <<'EOF'
chore: remove stale commented-out code blocks

- input_module: drop dead //free(alpha_idm_dr)/(beta_idr) — the raw
  free() ownership pattern is gone; the pointers are now non-owning
  views into vector storage, so calling free() on them was already
  wrong and the comments are dead.
- nonlinear_module: replace a long block-commented class_test_except
  + matching prose with a one-paragraph comment explaining the
  current behaviour.
- thermodynamics_module: drop ~12 //class_store_columntitle/double
  lines for derivatives that are not produced.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: `NULL` → `nullptr` sweep

**Files:**
- Modify: every `.cpp` under `source/`, `tools/` (≈170 occurrences)
- Modify: `include/common.h` (2 occurrences)
- Excluded: `hyrec/` (third-party); `include/hermite{3,4,6}_interpolation_csource.h` and `tools/hermite{3,4,6}_interpolation_csource.h` (C-source splice fragments — leaving them as `NULL` keeps them valid as C)

- [ ] **Step 1: Sanity-check the scope**

```bash
# Count NULL occurrences in each in-scope directory:
for d in source tools; do
  cnt=$(find "$d" -maxdepth 1 -type f -name "*.cpp" -print0 | xargs -0 grep -h -E "\bNULL\b" 2>/dev/null | wc -l | tr -d ' ')
  echo "$d/*.cpp: $cnt"
done
echo "include/common.h: $(grep -cE "\bNULL\b" include/common.h)"
```

Expected (approximately):
- `source/*.cpp: 145`
- `tools/*.cpp: 25`
- `include/common.h: 2`

If these are wildly different, stop and investigate before sweeping.

- [ ] **Step 2: Verify no stringify-NULL or "NULL" string would be touched**

```bash
grep -rnE '#NULL\b|"NULL"' source/ tools/ include/ species/ main/ 2>/dev/null | grep -v "\.claude/"
```

Expected: empty output. If anything matches, inspect and decide line-by-line — `\bNULL\b` does not match inside a string literal anyway, but `#NULL` (a preprocessor stringify) would change meaning. Last sanity confirms no such surprise.

- [ ] **Step 3: Run the sweep**

Use Perl for portability between BSD and GNU sed:

```bash
# In-place replace whole-word NULL with nullptr in the in-scope files.
find source -maxdepth 1 -type f -name "*.cpp" -print0 | \
  xargs -0 perl -i -pe 's/\bNULL\b/nullptr/g'
find tools  -maxdepth 1 -type f -name "*.cpp" -print0 | \
  xargs -0 perl -i -pe 's/\bNULL\b/nullptr/g'
perl -i -pe 's/\bNULL\b/nullptr/g' include/common.h
```

- [ ] **Step 4: Verify no `NULL` remains in scope and no over-reach**

```bash
echo "=== should be 0 each ==="
find source -maxdepth 1 -name "*.cpp" -print0 | xargs -0 grep -c -E "\bNULL\b" | grep -v ":0$" || echo "source/*.cpp clean"
find tools  -maxdepth 1 -name "*.cpp" -print0 | xargs -0 grep -c -E "\bNULL\b" | grep -v ":0$" || echo "tools/*.cpp clean"
grep -cE "\bNULL\b" include/common.h
echo "=== unchanged dirs (NULL still present is fine here) ==="
grep -rE "\bNULL\b" hyrec/ 2>/dev/null | wc -l
grep -E "\bNULL\b" include/hermite3_interpolation_csource.h tools/hermite3_interpolation_csource.h 2>/dev/null | wc -l
```

Expected: `source/*.cpp clean`, `tools/*.cpp clean`, `include/common.h` shows `0`, and the unchanged dirs still report non-zero `NULL` counts (proof we didn't sweep them).

- [ ] **Step 5: Spot-check the diff**

```bash
git diff --stat
git diff include/common.h
```

`git diff --stat` should show ~170 lines changed across ~15 files. The `common.h` diff should be only two `nullptr` substitutions — no whitespace churn, no nearby line edits.

- [ ] **Step 6: Build**

```bash
make class -j4
```

Expected: clean build. (`nullptr` is in `<cstddef>` which is transitively pulled in everywhere, but if any file uses `NULL` without including `<cstddef>` or `<cstdlib>`, it might also lack `nullptr`. In practice `nullptr` is a keyword, not a macro, so no header is needed.)

- [ ] **Step 7: Numerical regression check**

```bash
rm -rf /tmp/cleanup-new && mkdir -p /tmp/cleanup-new output
./class test/scenarios/gauge_lcdm.ini
mv output/gauge_lcdm*.dat /tmp/cleanup-new/
python test/scenarios/compare_tol.py /tmp/cleanup-baseline /tmp/cleanup-new
```

Expected: every file `OK`. `NULL` is `(void*)0` (or `0L` depending on platform) and `nullptr` is a `std::nullptr_t`; both compare equal to any null pointer and convert identically in pointer contexts.

- [ ] **Step 8: Commit**

```bash
git add source/ tools/ include/common.h
git commit -m "$(cat <<'EOF'
style: NULL -> nullptr in our C++ sources

Pure mechanical sweep with `perl -i -pe 's/\bNULL\b/nullptr/g'` over
source/*.cpp, tools/*.cpp, and include/common.h. hyrec/ (third-party)
and the *_csource.h C-splice fragments are left untouched so they
remain valid as C.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Final verification + PR

**Files:** none (verification + PR)

- [ ] **Step 1: Run the full existing test discovery**

```bash
ls test/scenarios/*.ini
```

Pick at least one additional scenario beyond `gauge_lcdm.ini` for a final cross-check. `gauge_idmdr.ini` is a good choice because it exercises code paths in modules we touched (input, thermodynamics) and the `idm_dr` interaction.

- [ ] **Step 2: Cross-check on a second scenario**

```bash
rm -rf /tmp/cleanup-baseline-2 /tmp/cleanup-new-2
mkdir -p /tmp/cleanup-baseline-2 /tmp/cleanup-new-2 output
# Baseline from master
git stash
make class -j4 >/dev/null
./class test/scenarios/gauge_idmdr.ini
mv output/gauge_idmdr*.dat /tmp/cleanup-baseline-2/
git stash pop
# New from our branch
make class -j4 >/dev/null
./class test/scenarios/gauge_idmdr.ini
mv output/gauge_idmdr*.dat /tmp/cleanup-new-2/
python test/scenarios/compare_tol.py /tmp/cleanup-baseline-2 /tmp/cleanup-new-2
```

Expected: every file `OK`. If anything fails on this scenario, investigate by reverting commits one at a time (`git revert HEAD`, rerun) to identify which commit caused the regression.

- [ ] **Step 3: Confirm commit log shape**

```bash
git log --oneline -5
```

Expected: four new commits on top of the previous master HEAD, in the order:

```
<hash> style: NULL -> nullptr in our C++ sources
<hash> chore: remove stale commented-out code blocks
<hash> refactor: replace file_exists with std::filesystem::exists
<hash> refactor: replace compare_doubles + qsort with std::sort
<hash> [previous master HEAD]
```

- [ ] **Step 4: Push and open PR**

```bash
git push -u origin HEAD
gh pr create --title "Cleanup: drop two C helpers, dead comments, NULL -> nullptr" --body "$(cat <<'EOF'
## Summary
- Drop `InputModule::compare_doubles` + `qsort`; use `std::sort`.
- Drop `InputModule::file_exists`; use `std::filesystem::exists` with the
  `error_code&` overload to preserve silent-false-on-error semantics.
- Remove three stale comment blocks (`input_module`, `nonlinear_module`,
  `thermodynamics_module`).
- `NULL` → `nullptr` across our `.cpp` files and `include/common.h`.
  `hyrec/` and the `*_csource.h` C-splice fragments are intentionally
  excluded.

Spec: `docs/superpowers/specs/2026-05-31-code-cleanup-pass-design.md`

## Test plan
- [x] `make class -j4` after each commit
- [x] `compare_tol.py` against `test/scenarios/gauge_lcdm.ini` (1e-3) after each commit
- [x] `compare_tol.py` against `test/scenarios/gauge_idmdr.ini` (1e-3) after all four commits

Out of scope (flagged for follow-up): the slab allocator in
`thermodynamics_module:3220-3260`, the `const_cast` wall around the
`array_*` API.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

Return the PR URL.

---

## Self-Review

**Spec coverage:**
- Spec §Scope/in-scope item 1 (`compare_doubles` → `std::sort`) → Task 1 ✓
- Spec §Scope/in-scope item 2 (`file_exists` → `std::filesystem::exists`) → Task 2 ✓
- Spec §Scope/in-scope item 3 (dead commented blocks at three sites) → Task 3 ✓
- Spec §Scope/in-scope item 4 (`NULL` → `nullptr` sweep with named exclusions) → Task 4 ✓
- Spec §Verification (build + `compare_tol.py` on `gauge_lcdm.ini`) → present in every implementation task; plus a second scenario in Task 5 ✓
- Spec §Risks: `-lstdc++fs` flag and `#NULL` stringify → addressed in Task 2 Step 5 note and Task 4 Step 2 sanity check ✓

No gaps.

**Placeholders:** none — every step has an exact command or exact code block.

**Type consistency:** N/A (pure refactor, no new types). Function names referenced (`InputModule::compare_doubles`, `InputModule::file_exists`) match the codebase verbatim.
