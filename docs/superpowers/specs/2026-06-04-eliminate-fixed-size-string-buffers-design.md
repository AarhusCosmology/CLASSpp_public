# Eliminate fixed-size string buffers (issue #291, expanded)

**Date:** 2026-06-04
**Branch target:** new branch off `master`
**Issue:** #291 — originally "replace `FileName`/`FileArg` with `std::string`", expanded during
brainstorming to the full removal of fixed-size string buffers across the C++ codebase.

## Motivation

The codebase still uses pre-C++ fixed-size character buffers for paths, parameter values, lines,
and built strings:

- `typedef char FileName[_FILENAMESIZE_];`
- `typedef char FileArg[_ARGUMENT_LENGTH_MAX_];`
- bare `char buf[_LINE_LENGTH_MAX_]` / `char buf[N]` locals
- `snprintf`-based string construction (snprintf is everywhere only because `sprintf` is unsafe)
- `strcpy`/`strncpy`/`strcat`/`strstr`/`strcmp` on those buffers

These carry the classic hazards: silent truncation and buffer overflow on concatenation. The
fixed sizes are arbitrary and cause real, user-visible bugs — e.g. a `k_output_values` line longer
than `_LINE_LENGTH_MAX_` (1024) is split mid-line by `fgets` in `parser_read_file` and the tail is
mis-parsed as a separate line. Migrating to `std::string` removes the overflow/truncation class
entirely and modernizes the call sites.

Context: this came up in #290 (in-struct defaults), where `precision::ResolveDataPaths()` was
hardened with bounded `snprintf` as a localized stopgap and the broader migration was deferred to
this issue.

## Scope

**In scope** — eliminate every fixed-size `char[]` buffer used for a path, parameter value, input
line, or built string, and the constants that size them:

- `FileName` / `FileArg` typedefs and all their users (struct fields, function params, locals).
- The constants `_FILENAMESIZE_`, `_ARGUMENT_LENGTH_MAX_`, `_LINE_LENGTH_MAX_` — deleted once
  unused.
- The line-ingestion path in `tools/parser.cpp` (`fgets` into `char line_buf[_LINE_LENGTH_MAX_]`).
- `parser_read_string` and its 47 call sites.
- `snprintf` string-building (~73 sites, mostly `output_module.cpp` file-name construction).
- `strcpy`/`strncpy`/`strcat` on our buffers (~56 sites).
- `strstr`/`strcmp`/`strncmp` on our strings (~101 sites) — converted to idiomatic `std::string`
  (`.find(...) != std::string::npos`, `==`).
- The Cython wrapper for any wrapper-exposed buffer that changes type
  (`generate_wrapper.py`, regenerated `cclassy.pxd`, `classy.pyx`).
- The column-title mechanism: `char titles[_MAXTITLESTRINGLENGTH_]` and the constant
  `_MAXTITLESTRINGLENGTH_` (decided in scope during planning). The `*_output_titles()` methods
  across background/thermo/primordial/perturbations become `std::string&`, the 3 `perturbations`
  title fields and the `output_module` title locals become `std::string`, `class_store_columntitle`
  appends to a `std::string`, and the 4 `classy.pyx` title-reading sites drop the
  `resize`/`<char*> c_str()` contract. Delete the `ClassConstants::sMAXTITLESTRINGLENGTH` /
  `sARGUMENT_LENGTH_MAX` exports once unused.

**Numeric-formatting rule (avoids breaking bit-identicality):** when a number must be turned into
text, integers use `std::to_string` (exact, no format change); any value whose textual form must
match a `printf` conversion (`%g`, `%e`, `%.16e`, width/precision) is formatted with `snprintf` into
a small *scalar* numeric scratch and concatenated into the `std::string`. Such a single-scalar
scratch carries no truncation hazard and is not a "string/path/line buffer" — it is the one
permitted `char` local. Never introduce a fixed-size buffer to hold a path, parameter value, list,
or title.

**Out of scope** — explicitly untouched:

- `ErrorMsg` / `_ERRORMSGSIZE_` / `class_build_error_string` and the
  `class_call`/`class_test`/`class_stop` error macros (1761 occurrences, 1212 macro sites, C-varargs
  formatting, exposed through every module's `error_message_`). A separate future PR.
- `char** argv` program-entry interfaces (OS boundary, not a fixed-size buffer).

The work concentrates in `source/input_module.cpp` (~135 hits) and `source/output_module.cpp`
(~118 hits); the rest is spread thin across ~17 other files.

## Constraints / findings

1. **`FileContent` already stores values as `std::string` / `std::vector<double>` internally**
   (`read_string`, `read_list_of_doubles`). The truncation lives only in the C-API copy layer
   (`parser_read_string` → fixed buffer) and at line ingestion (`fgets` → `line_buf`). Fixing those
   two layers removes truncation everywhere downstream.

2. **None of `precision` / `transfers` / `output` are ever `memset`/`memcpy`'d**, so `std::string`
   members are safe in those structs.

3. **`class_open`/`fopen`/`%s` consumers** of the converted path strings just need `.c_str()`.

4. **Wrapper-exposed buffers** are the reason this needs wrapper work:
   - `precision::class_dir` (`FileArg`) appears in `cclassy.pxd` and is read in `classy.pyx`.
   - `perturb_output_firstline_and_ic_suffix(..., FileName ic_suffix)` is called from `classy.pyx`
     (`<bytes> ic_suffix`).
   - `generate_wrapper.py` emits method params verbatim (return type must be in `allowed_types`;
     that is why `FileName ic_suffix` currently leaks into the pxd), and filters struct fields by a
     separate `allowed_types` list that lacks `std::string`. Both lists need `std::string` added;
     the `FileArg`/`FileName` ctypedef emission is deleted (keep `ErrorMsg`).
   - `cclassy.pxd` is auto-generated — regenerate it, never hand-edit.

## The transformations

1. **Line ingestion (the `k_output_values` fix).** `tools/parser.cpp`: replace
   `char line_buf[_LINE_LENGTH_MAX_]` + `std::fgets` with `std::ifstream` + `std::getline(stream,
   std::string)`. Arbitrarily long lines parse correctly. No regression test is added — correctness
   here follows from principle (unbounded read), and the existing suites exercise the parser.

2. **`parser_read_string`** → `int parser_read_string(FileContent*, const char* name,
   std::string& value, int* found, ErrorMsg errmsg)`; fill `value` directly; drop the now-unneeded
   length `class_test`. (`ErrorMsg errmsg` param stays — error path is out of scope.)
   - 47 call sites: `&string1` → `string1`; `char string1[_ARGUMENT_LENGTH_MAX_]` →
     `std::string string1`.
   - Downstream `strstr(string1, "x")` → `string1.find("x") != std::string::npos`;
     `strcmp(string1, "x") == 0` → `string1 == "x"`.
   - `class_read_string` macro: `strcpy(destination, string1)` → `destination = string1;`.

3. **Stored path/param buffers → `std::string`:**
   `precision::{sBBN_file, hyrec_Alpha_inf_file, hyrec_R_inf_file, hyrec_two_photon_tables_file,
   class_dir}`, `transfers::{nz_file_name, nz_evo_file_name}`, `output::root`, and the `FileArg`/
   `char[]` locals (`stringoutput`, `inifilename`, …). Keep literal defaults.
   - `read(const FileContent&, const char*, std::string&)` overload replaces the `FileName&` one.
   - `precision::ResolveDataPaths()` collapses to `path = std::string(class_dir) + path;` per field
     (delete the `snprintf` lambda).
   - `class_dir` read: `parser_read_string(..., ppr->class_dir, ...)`; `__CLASSDIR__` fallback →
     `ppr->class_dir = __CLASSDIR__;`.
   - Thermo consumers (~10 in `thermodynamics_module.cpp`) and transfer consumers (`fopen`, `%s`)
     get `.c_str()`.

4. **`snprintf` string-building → C++ concatenation.** The ~9 `file_name` locals and similar in
   `output_module.cpp` become `std::string`, built by concatenation (`pop->root + "cl.dat"`,
   `pop->root + redshift_suffix + "tk.dat"`, etc.) instead of `snprintf(buf, SIZE-32, "%s%s", …)`.
   `output_open_cl_file`/`output_open_pk_file` params `FileName filename` → `const std::string&`;
   internal `class_open(*file, filename.c_str(), …)`. The `param_output_name`/`param_unused_name`
   builders in `input_module.cpp` likewise become concatenation.

5. **`strcpy`/`strncpy`/`strcat`** on our buffers → assignment / `+=`.

6. **Small/local `char[]`** → `std::string`: `extension`, `tmp`, `redshift_suffix`, `first_line`,
   `ic_info`, `ic_suffix`, the param-file-name buffers, and any other non-`ErrorMsg`, non-`argv`
   fixed buffer surfaced during implementation.

7. **`perturb_output_firstline_and_ic_suffix`** → `std::string& ic_suffix`;
   `strcpy(ic_suffix, "ad")` → `ic_suffix = "ad"`; the `output_module.cpp` caller's
   `char ic_suffix[4]` → `std::string` (snprintf using it → concatenation). `first_line` similarly.

8. **Delete typedefs + constants** from `include/common.h`: `FileName`, `FileArg`,
   `_FILENAMESIZE_`, `_ARGUMENT_LENGTH_MAX_`, `_LINE_LENGTH_MAX_`, once no references remain. Keep
   `ErrorMsg` / `_ERRORMSGSIZE_`.

## Wrapper tier

- `generate_wrapper.py`: add `std::string` to both `allowed_types` lists (struct fields ~line 193,
  method signatures ~line 266); delete the `FileArg` and `FileName` ctypedef `enums.append(...)`
  lines (keep `ErrorMsg`). The verbatim method-param emission then renders `std::string&` →
  `string&` via the existing `replace('std::','')`.
- Regenerate `cclassy.pxd`. Expect: `class_dir` and the now-`std::string` precision path fields to
  appear as `string` struct members (harmless read-only declarations); `ic_suffix` param as
  `string&`; `FileArg`/`FileName` ctypedefs gone.
- `classy.pyx`: `FileName ic_suffix` local → libcpp `string` (`from libcpp.string cimport string`
  already present); adjust the `<bytes> ic_suffix` extraction to the std::string form. Verify no
  struct-field access of `ppr.class_dir` exists (line 169 is a Python-dict entry, not a struct
  read).

## Commits (one PR)

1. Parser: line ingestion (`getline`) + `parser_read_string` → `std::string&` + call sites.
2. Precision / thermodynamics data-file paths + `class_dir` + `ResolveDataPaths`.
3. `transfers` nz fields + consumers.
4. `output_module` snprintf → string concatenation + `output_open_*` params + `output::root`.
5. Remaining locals / `perturb_output_firstline_and_ic_suffix` / misc `char[]`.
6. Wrapper: `generate_wrapper.py` + regenerated `cclassy.pxd` + `classy.pyx`.
7. Delete `FileName`/`FileArg` typedefs and the three constants from `common.h`; confirm zero
   references remain.

## Verification

- **After each non-wrapper commit:** build the `class` binary and confirm output is **bit-identical**
  to `master` for a reference run in **both gauges** (synchronous and newtonian). The only permitted
  difference is the `k_output_values` long-line case, which changes from broken to correct (not
  exercised by the bit-identical reference run).
- **After the wrapper commit:** `pip install .` succeeds and the `COMPARE_OUTPUT_REF` suite passes.
- **Final:** `grep -rn "FileName\|FileArg\|_FILENAMESIZE_\|_ARGUMENT_LENGTH_MAX_\|_LINE_LENGTH_MAX_"`
  over `source/`, `tools/`, `include/`, `main/`, `classy.pyx`, `cclassy.pxd`, `generate_wrapper.py`
  returns nothing (except inside this spec / docs).

## Risks

- **Hidden POD assumptions.** Audited: `precision`/`transfers`/`output` are not memset/memcpy'd.
  Re-check any newly touched struct before converting its fields.
- **Wrapper regeneration surface.** Adding `std::string` to the struct `allowed_types` exposes a few
  extra read-only fields in `cclassy.pxd`; harmless but widens the generated diff. Confirm the pxd
  still compiles and `classy` imports.
- **Semantic drift in `strstr`→`.find`.** `strstr` matches substrings; `.find != npos` is the exact
  equivalent. `strcmp == 0` ↔ `==`. Keep substring-vs-equality matching identical to the original at
  each site.
- **Scale.** ~318 sites; mitigated by small per-area commits each gated on bit-identical output.
