# Eliminate fixed-size string buffers (#291) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove every fixed-size `char[]` buffer used for paths, parameter values, input lines, built strings, and column-title strings — replacing them with `std::string` — and delete the constants/typedefs that sized them (`_FILENAMESIZE_`, `_ARGUMENT_LENGTH_MAX_`, `_LINE_LENGTH_MAX_`, `_MAXTITLESTRINGLENGTH_`, `FileName`, `FileArg`).

**Architecture:** A behavior-preserving mechanical sweep. The underlying `FileContent` store is already `std::string`-based; the buffers live in the C-API copy layer, in builders (`snprintf`/`strcat`), in line ingestion (`fgets`), and in the title/Cython contract. We convert area-by-area, gating each non-wrapper commit on **bit-identical `class` binary output in both gauges**. The Cython wrapper (`generate_wrapper.py` → regenerated `cclassy.pxd`, plus `classy.pyx`) is updated last, gated on `pip install .` + the `COMPARE_OUTPUT_REF` suite.

**Tech Stack:** C++17, the `class` Makefile build, the `tools/parser_test.cpp` unit test (`make test-parser`), Cython wrapper (`generate_wrapper.py`/`cclassy.pxd`/`classy.pyx`), `python/test_class.py` (`COMPARE_OUTPUT_REF`).

**Out of scope (explicit):** `ErrorMsg` / `_ERRORMSGSIZE_` and the C-varargs error/printf formatting buffers — `class_build_error_string`, `tools/exceptions.cpp` (`char error_message[...]`), `tools/common.cpp` (`class_protect_printf`/`class_protect_fprintf` `char dest[...]`). These stay (separate future PR). `char** argv` program-entry interfaces stay.

---

## File map

| File | Responsibility / change |
|---|---|
| `include/common.h` | Delete `FileName`, `FileArg`, `_FILENAMESIZE_`, `_ARGUMENT_LENGTH_MAX_`, `_LINE_LENGTH_MAX_`, `_MAXTITLESTRINGLENGTH_`; rewrite `class_store_columntitle` for `std::string`. |
| `include/parser.h` / `tools/parser.cpp` | `getline` ingestion; `parser_read_string(... std::string&)`; delete dead `parser_init`/`parser_read_line`. |
| `source/input_module.{h,cpp}` | 48 `parser_read_string` sites + downstream `strstr`/`strcmp` → `.find`/`==`; `class_read_string` macro; `class_dir`; all `char[]` locals; `param_*_name` builders. |
| `source/output_module.{h,cpp}` | `file_name`/`first_line`/`ic_suffix`/`titles`/`thetitle` locals → `std::string`; `snprintf`→concat; `output_open_*`/`print_*` params. |
| `source/thermodynamics_module.{h,cpp}` | path `.c_str()` consumers; `line` buffer; `output_titles` → `std::string&`. |
| `source/transfer.h` / `source/transfer_module.cpp` | `nz_*_file_name` fields → `std::string`; `fopen`/`%s` `.c_str()`. |
| `source/perturbations_module.{h,cpp}` | `ic_suffix`/`first_line`/`titles` params → `std::string&`; 3 title fields → `std::string`. |
| `source/primordial_module.{h,cpp}`, `source/background_module.{h,cpp}` | `output_titles` → `std::string&`; primordial command-builder buffers → `std::string`. |
| `source/output.h` | `root` field → `std::string`. |
| `tools/common.cpp` | `get_number_of_titles(const std::string&)`. |
| `tools/quadrature.cpp` | `method_chosen` buffer → `std::string`. |
| `main/class.h` | Drop `ClassConstants::sMAXTITLESTRINGLENGTH`, `sARGUMENT_LENGTH_MAX` (unused once buffers gone). |
| `generate_wrapper.py` | Add `std::string` to both `allowed_types`; delete `FileArg`/`FileName` ctypedef emission. |
| `cclassy.pxd` | Regenerated (never hand-edited). |
| `classy.pyx` | `ic_suffix`/`titles` locals → `string` (no `resize`/`<char*>`); drop `constvals.sMAXTITLESTRINGLENGTH`. |
| `scripts/verify_bitidentical.sh` (new, temp) | Build + run both gauges + diff vs golden. Deleted in final task. |

---

## Task 0: Verification harness (golden reference on master)

**Files:**
- Create: `scripts/verify_bitidentical.sh`
- Create (gitignored, local): `/tmp/scf_gauge/ref_sync/`, `/tmp/scf_gauge/ref_newt/`, `/tmp/scf_gauge/ref.ini`

- [ ] **Step 1: Write a reference ini that exercises every converted writer**

Create `/tmp/scf_gauge/ref.ini`:

```ini
output = tCl,pCl,lCl,mPk,mTk,vTk
root = /tmp/scf_gauge/out_
write background = yes
write thermodynamics = yes
write primordial = yes
write parameters = yeap
k_output_values = 0.001, 0.01, 0.1
P_k_max_h/Mpc = 1.
z_pk = 0, 1
lensing = yes
gauge = synchronous
```

- [ ] **Step 2: Write the verify script**

```bash
#!/usr/bin/env bash
# scripts/verify_bitidentical.sh — rebuild class and diff output vs golden refs (both gauges).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK=/tmp/scf_gauge
INI="$WORK/ref.ini"
cd "$ROOT"
make class -j >/dev/null
run_gauge () {  # $1=gauge $2=destdir
  rm -f "$WORK"/out_*
  sed "s/^gauge = .*/gauge = $1/" "$INI" > "$WORK/run_$1.ini"
  ./class "$WORK/run_$1.ini" >/dev/null
  mkdir -p "$2"; rm -f "$2"/*; cp "$WORK"/out_* "$2"/ 2>/dev/null || true
}
if [ "${1:-}" = "--save" ]; then
  run_gauge synchronous "$WORK/ref_sync"
  run_gauge newtonian   "$WORK/ref_newt"
  echo "Saved golden refs."
else
  run_gauge synchronous "$WORK/cur_sync"
  run_gauge newtonian   "$WORK/cur_newt"
  # Compare ignoring only the header line that prints the (changing) absolute ini path.
  d=0
  for g in sync newt; do
    for f in "$WORK/ref_$g"/*; do
      b="$(basename "$f")"
      if ! diff -q <(grep -av "$WORK/run_" "$f") <(grep -av "$WORK/run_" "$WORK/cur_$g/$b") >/dev/null; then
        echo "DIFF in $g/$b"; d=1
      fi
    done
  done
  [ "$d" = 0 ] && echo "BIT-IDENTICAL (both gauges)" || { echo "FAILED"; exit 1; }
fi
```

- [ ] **Step 3: Make it executable and capture the golden reference on master**

Run:
```bash
chmod +x scripts/verify_bitidentical.sh
git stash list   # ensure clean tree on the branch tip (no uncommitted source changes)
./scripts/verify_bitidentical.sh --save
```
Expected: `Saved golden refs.` and files present in `/tmp/scf_gauge/ref_sync` and `ref_newt`.

- [ ] **Step 4: Sanity-check the harness detects nothing yet**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 5: Commit**

```bash
git add scripts/verify_bitidentical.sh
git commit -m "Add bit-identical verification harness for #291 sweep"
```

**Note for all later non-wrapper tasks:** the verification step is always `./scripts/verify_bitidentical.sh` expecting `BIT-IDENTICAL (both gauges)`. The golden refs were captured on master, so they encode the pre-refactor behavior.

---

## Task 1: Parser layer — `getline` ingestion + `parser_read_string` → `std::string&`

**Files:**
- Modify: `tools/parser.cpp:14-45` (from_file), `:267-304` (dead wrappers), `:328-346` (parser_read_string)
- Modify: `include/parser.h:125-137`
- Modify: `source/input_module.h:107-112` (`class_read_string` macro)
- Modify: `source/input_module.cpp` — 48 `parser_read_string` call sites + their `char string1/string2[...]` locals + downstream `strstr`/`strcmp`
- Test: `make test-parser`, then `./scripts/verify_bitidentical.sh`

- [ ] **Step 1: Replace `fgets` ingestion with `std::getline`**

`tools/parser.cpp` `FileContent::from_file` (currently opens `FILE*` and loops `char line_buf[_LINE_LENGTH_MAX_]` + `std::fgets`). Replace the buffer+fgets loop with:

```cpp
/* static */ FileContent FileContent::from_file(const std::string& filename) {
  std::ifstream in(filename);
  if (!in) {
    throw std::invalid_argument("Cannot open file '" + filename + "': " + std::strerror(errno));
  }
  FileContent fc;
  fc.filename_ = filename;

  std::string line;
  while (std::getline(in, line)) {
    std::string name, value;
    if (!parse_line(line, name, value))
      continue;
    if (fc.params_.count(name)) {
      throw std::invalid_argument("Multiple entries of parameter '" + name + "' in file '" +
                                  filename + "'");
    }
    fc.keys_.push_back(name);
    fc.params_[name] = value;
  }
  if (fc.params_.empty()) {
    throw std::invalid_argument("No readable input in file '" + filename + "'");
  }
  return fc;
}
```

Add `#include <fstream>` near the top of `tools/parser.cpp` (keep `<cstring>` for `strerror`). Remove the now-unused `<cstdio>` only if nothing else uses it (the dead wrappers below also go).

- [ ] **Step 2: Delete the dead C wrappers `parser_init` and `parser_read_line`**

Both have **0 call sites** (verified). Delete their definitions in `tools/parser.cpp` (`:277-283` and `:285-304`) and declarations in `include/parser.h` (`:127`, `:129`). This removes their `_ARGUMENT_LENGTH_MAX_` references. Leave `parser_read_file` (2 callers) and `parser_cat`.

- [ ] **Step 3: Convert `parser_read_string` to fill `std::string&`**

`include/parser.h:136-137`:
```cpp
int parser_read_string(
    FileContent* pfc, const char* name, std::string& value, int* found, ErrorMsg errmsg);
```
`tools/parser.cpp:328-346`:
```cpp
int parser_read_string(
    FileContent* pfc, const char* name, std::string& value, int* found, ErrorMsg errmsg) {
  try {
    *found = pfc->read_string(name, value) ? _TRUE_ : _FALSE_;
  }
  catch (const std::exception& e) {
    class_stop(errmsg, "%s", e.what());
  }
  return _SUCCESS_;
}
```
(The length `class_test` is gone — `std::string` is unbounded.)

- [ ] **Step 4: Update the `class_read_string` macro**

`source/input_module.h:107-112` — change the copy line:
```cpp
#define class_read_string(name, destination)                                             \
  do {                                                                                   \
    class_call(parser_read_string(pfc, name, string1, &flag1, errmsg), errmsg, errmsg);  \
    if (flag1 == _TRUE_)                                                                 \
      destination = string1;                                                             \
  } while (0);
```
(`string1` is now `std::string`; passed by reference, no `&`. `destination` assignment works for the `std::string` nz fields — see Task 3.)

- [ ] **Step 5: Convert the `char stringN[_ARGUMENT_LENGTH_MAX_]` locals and all 48 call sites**

In `source/input_module.cpp`:
- Every `char string1[_ARGUMENT_LENGTH_MAX_];` (`:533`, `:646`, `:903`) and `char string2[_ARGUMENT_LENGTH_MAX_];` (`:904`) → `std::string string1;` / `std::string string2;`.
- Every `parser_read_string(pfc, "X", &string1, &flag1, errmsg)` → `parser_read_string(pfc, "X", string1, &flag1, errmsg)` (drop the `&`). Same for `string2`, `stringoutput`, `inifilename` (those locals are converted in Task 2/6 but the `&`-drop happens wherever the call appears). Sites: the 48 enumerated by `grep -n "parser_read_string(" source/input_module.cpp`.
- Downstream uses of these now-`std::string` locals (transform each, keeping substring-vs-equality identical):
  - `strstr(string1, "needle") != nullptr` → `string1.find("needle") != std::string::npos`
  - `(strstr(a,"x")||strstr(a,"y"))` → `(a.find("x")!=std::string::npos || a.find("y")!=std::string::npos)`
  - `strcmp(string1, "lit") == 0` → `string1 == "lit"`
  - `strcmp(string1, "lit") != 0` → `string1 != "lit"`
  - `strncmp(string1, "lit", n) == 0` → `string1.compare(0, n, "lit") == 0`
  - any `strcpy(dst, string1)` where `dst` is char[] staying char[] → use `string1.c_str()`; where `dst` becomes `std::string` → `dst = string1;`
  - Walk every `string1`/`string2` reference with `grep -n "string1\|string2" source/input_module.cpp` and convert. Do not leave any `char*`-API call taking a bare `std::string`.

- [ ] **Step 6: Build the parser unit test and verify**

Run: `make test-parser && ./test-parser`
Expected: parser tests PASS (no truncation regressions).

- [ ] **Step 7: Build class and verify bit-identical**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

(`_ARGUMENT_LENGTH_MAX_`/`_LINE_LENGTH_MAX_` may still be referenced elsewhere — they are deleted in Task 9. Do not delete them yet; the build must stay green.)

- [ ] **Step 8: Commit**

```bash
git add tools/parser.cpp include/parser.h source/input_module.h source/input_module.cpp
git commit -m "Parser: getline ingestion + parser_read_string(std::string&) + call sites (#291)"
```

---

## Task 2: precision data-file paths + `class_dir` + `ResolveDataPaths`

**Files:**
- Modify: `include/common.h:487` (`class_dir`), `:563`, `:663-667` (path fields)
- Modify: `source/input_module.cpp:595-604` (class_dir read+fallback), `:2638-2664` (read overload + ResolveDataPaths), `:2684,2722-2724` (reads)
- Modify: `source/thermodynamics_module.cpp:1688,1712,1719,1738,3276-3307` (`.c_str()` consumers)
- Test: `./scripts/verify_bitidentical.sh`

- [ ] **Step 1: Convert the fields to `std::string`**

`include/common.h`:
- `:487` `FileArg class_dir;` → `std::string class_dir;`
- `:563` `FileName sBBN_file = "/bbn/sBBN_2017.dat";` → `std::string sBBN_file = "/bbn/sBBN_2017.dat";`
- `:663-667` the three `FileName hyrec_*` fields → `std::string` (keep their literal default initializers verbatim).

- [ ] **Step 2: Replace the `FileName&` read overload with a `std::string&` one**

`source/input_module.cpp:2638-2642`:
```cpp
void read(const FileContent& fc, const char* name, std::string& v) {
  std::string s;
  if (fc.read_string(name, s))
    v = s;
}
```

- [ ] **Step 3: Simplify `ResolveDataPaths`**

`source/input_module.cpp:2651-2665` → replace the `snprintf` lambda body with direct concatenation:
```cpp
void precision::ResolveDataPaths() {
  // Prepend the runtime class_dir to each field's relative-path default.
  sBBN_file                   = class_dir + sBBN_file;
  hyrec_Alpha_inf_file        = class_dir + hyrec_Alpha_inf_file;
  hyrec_R_inf_file            = class_dir + hyrec_R_inf_file;
  hyrec_two_photon_tables_file = class_dir + hyrec_two_photon_tables_file;
}
```

- [ ] **Step 4: Convert the `class_dir` read + fallback**

`source/input_module.cpp:595-604`:
```cpp
  class_call(parser_read_string(&file_content_, "class_dir", ppr->class_dir, &flag1, error_message_),
             error_message_, error_message_);
  if (flag1 == _FALSE_) {
    ppr->class_dir = __CLASSDIR__;
  }
```

- [ ] **Step 5: Add `.c_str()` at the thermodynamics consumers**

`source/thermodynamics_module.cpp`:
- `:1688` `class_open(fA, ppr->sBBN_file, "r", error_message_);` → `... ppr->sBBN_file.c_str() ...`
- `:1712,1719,1738` the `%s` args `ppr->sBBN_file` → `ppr->sBBN_file.c_str()`
- `:3276` `class_open(fA, ppr->hyrec_Alpha_inf_file, ...)` → `.c_str()`; `:3277` `hyrec_R_inf_file` → `.c_str()`
- `:3283,3289,3307` `%s` args → `.c_str()` for the matching field
- `:3297` `class_open(fA, ppr->hyrec_two_photon_tables_file, ...)` → `.c_str()`

- [ ] **Step 6: Verify bit-identical**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 7: Commit**

```bash
git add include/common.h source/input_module.cpp source/thermodynamics_module.cpp
git commit -m "precision: data-file paths + class_dir to std::string; simplify ResolveDataPaths (#291)"
```

---

## Task 3: `transfers` nz file-name fields

**Files:**
- Modify: `source/transfer.h:102,107`
- Modify: `source/transfer_module.cpp:3557-3596`
- Modify: `source/input_module.cpp:2176,2190` (already use `class_read_string` macro from Task 1)
- Test: `./scripts/verify_bitidentical.sh`

- [ ] **Step 1: Convert the fields**

`source/transfer.h`: `:102` `FileName nz_file_name;` → `std::string nz_file_name;`; `:107` `FileName nz_evo_file_name;` → `std::string nz_evo_file_name;` (keep doc comments).

- [ ] **Step 2: Add `.c_str()` at the consumers**

`source/transfer_module.cpp`:
- `:3557` `input_file = fopen(ptr->nz_file_name, "r");` → `fopen(ptr->nz_file_name.c_str(), "r");`
- `:3558` `%s` arg `ptr->nz_file_name` → `.c_str()`
- `:3592` `fopen(ptr->nz_evo_file_name, "r")` → `.c_str()`
- `:3596` `%s` arg → `.c_str()`

(The two `class_read_string("dNdz_*", ptr->nz_*_file_name)` sites at `input_module.cpp:2176,2190` now work via the Task-1 macro `destination = string1;`.)

- [ ] **Step 3: Verify bit-identical**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 4: Commit**

```bash
git add source/transfer.h source/transfer_module.cpp
git commit -m "transfers: nz_file_name/nz_evo_file_name to std::string (#291)"
```

---

## Task 4: `output_module` file-name building (`snprintf` → concatenation) + `output::root`

**Files:**
- Modify: `source/output.h:25` (`root`)
- Modify: `source/output_module.h:29,32,34` (params)
- Modify: `source/output_module.cpp` — all `FileName file_name` locals + `snprintf(file_name, …)` + `output_open_*`/`print_*` bodies
- Modify: `source/input_module.cpp:542-545` (`param_*_name` builders)
- Test: `./scripts/verify_bitidentical.sh`

- [ ] **Step 1: Convert `output::root`**

`source/output.h:25` `char root[_FILENAMESIZE_ - 32] = "output/";` → `std::string root = "output/";` (keep doc comment). Then every `pop->root` use that goes into a C API (`printf("%s", pop->root)` at `output_module.cpp:101`, and the concatenations below) reads correctly as `std::string`. For bare `%s` printf of root, use `pop->root.c_str()`.

- [ ] **Step 2: Convert `output_open_*` signatures**

`source/output_module.h`:
```cpp
int output_open_cl_file(FILE** clfile, const std::string& filename, const char* first_line, int lmax);
int output_open_pk_file(FILE** pkfile, const std::string& filename, const char* first_line, double z);
```
`source/output_module.cpp:1226` and `:1457` parameter `FileName filename` → `const std::string& filename`; inside each body, `class_open(*clfile, filename, "w", …)` → `class_open(*clfile, filename.c_str(), "w", …)`.

- [ ] **Step 3: Convert the `file_name` locals and replace `snprintf` builders with concatenation**

Every `FileName file_name;` local (`output_module.cpp:183,465,890,974,1017,1073,1152`) → `std::string file_name;`. Replace each builder. Distinct shapes (apply to all matching lines found by `grep -n "snprintf(file_name" source/output_module.cpp`):

```cpp
// snprintf(file_name, _FILENAMESIZE_ - 32, "%s%s", pop->root, "cl.dat");
file_name = pop->root + "cl.dat";

// snprintf(file_name, _FILENAMESIZE_ - 32, "%s%s%s", pop->root, redshift_suffix, "tk.dat");
file_name = pop->root + redshift_suffix + "tk.dat";

// snprintf(file_name, ..., "%s%s%s%s%s", pop->root, redshift_suffix, "tk_", ic_suffix, ".dat");
file_name = pop->root + redshift_suffix + "tk_" + ic_suffix + ".dat";

// pk filename built from type_suffix + redshift_suffix (see :530-538):
//   snprintf(file_name, ..., "%s%s%s.dat", pop->root, redshift_suffix, type_suffix);
file_name = pop->root + redshift_suffix + type_suffix + ".dat";
```
For pk ic-decomposition builders (`:549..`), follow the same pattern: concatenate `pop->root + redshift_suffix + "<infix>_" + ic_suffix + ".dat"` matching the original literal pieces exactly. Because `pop->root` (the left operand) is now `std::string`, the chain works even while `type_suffix`/`ic_suffix` are still `char[]` at this point (they decay to `const char*`, and `std::string + const char*` is valid). Each `class_call(output_open_*_file(&out, file_name, …))` now passes a `std::string` to the `const std::string&` param — no change at the call.

Also convert `redshift_suffix` itself (`output_module.cpp:866-872`, currently `char redshift_suffix[redshift_suffix_size]` built by `snprintf(..., "z%d_", index_z + 1)`) — it is an integer-formatted string, so use exact `std::to_string`:
```cpp
std::string redshift_suffix;
if (pop->z_pk_num != 1)
  redshift_suffix = "z" + std::to_string(index_z + 1) + "_";
```
(Delete the `redshift_suffix_size` const and the `char` decl.) Likewise convert `type_suffix` (`:496-505`, set from string literals `"pk"`/`"pk_nl"`/`"pk_cb"`/`"pk_cb_nl"`) → `std::string type_suffix;` with `type_suffix = "pk";` etc.

- [ ] **Step 4: Convert the `param_*_name` builders in input_module**

`source/input_module.cpp:542-545`:
```cpp
    std::string param_output_name = std::string(pop->root) + "parameters.ini";
    std::string param_unused_name = std::string(pop->root) + "unused_parameters";
```
(Replace the two `char ...[_LINE_LENGTH_MAX_]` decls and the two `snprintf` lines.) Update the consumers of these two locals (the `fopen`/`%s` uses just below) with `.c_str()`.

- [ ] **Step 5: Verify bit-identical**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 6: Commit**

```bash
git add source/output.h source/output_module.h source/output_module.cpp source/input_module.cpp
git commit -m "output: build file names with std::string concatenation; root to std::string (#291)"
```

---

## Task 5: Column-title mechanism (`_MAXTITLESTRINGLENGTH_`)

**Files:**
- Modify: `include/common.h:414-420` (`class_store_columntitle` macro)
- Modify: `source/{background,thermodynamics,primordial}_module.{h,cpp}`, `source/perturbations_module.{h,cpp}` (`*_output_titles` + 3 title fields)
- Modify: `source/output_module.{h,cpp}` (`titles`/`thetitle` locals + `print_*` titles param)
- Modify: `tools/common.cpp` (`get_number_of_titles`)
- Test: `./scripts/verify_bitidentical.sh`

- [ ] **Step 1: Rewrite the title-append macro for `std::string`**

`include/common.h:414-420`:
```cpp
#define class_store_columntitle(titlestring, title, condition) \
  {                                                            \
    if (condition == _TRUE_) {                                 \
      titlestring += title;                                    \
      titlestring += _DELIMITER_;                              \
    }                                                          \
  }
```

- [ ] **Step 2: Convert the `*_output_titles` signatures to `std::string&`**

For each declaration/definition pair, `char titles[_MAXTITLESTRINGLENGTH_]` → `std::string& titles`:
- `source/background_module.h:16` / `source/background_module.cpp:1242`
- `source/thermodynamics_module.h:13` / `source/thermodynamics_module.cpp:4285`
- `source/primordial_module.h:17` / `source/primordial_module.cpp:3006`
- `source/perturbations_module.h:19-20` / `source/perturbations_module.cpp:377` (`perturb_output_titles(enum file_format, std::string& titles)`)

Inside each body the only writes are `class_store_columntitle(titles, ...)` (now string-appends) — no other change. Remove any initial `titles[0] = '\0';` (an empty `std::string` is already empty).

- [ ] **Step 3: Convert the 3 perturbations title fields**

`source/perturbations_module.h:193,195,197` `char scalar_titles_[_MAXTITLESTRINGLENGTH_]` (and vector/tensor) → `std::string scalar_titles_;` etc. Find their fill/read sites with `grep -n "scalar_titles_\|vector_titles_\|tensor_titles_" source/perturbations_module.cpp` and ensure: fills use `class_store_columntitle` (fine) or assignment; any `strlen`/`%s`/`strcpy` consumer gets `.c_str()` or `.size()`.

- [ ] **Step 4: Convert `output_module` title locals and `print_*` param**

`source/output_module.cpp:818,976,1018,1153` `char titles[_MAXTITLESTRINGLENGTH_] = {0};` → `std::string titles;`. `:1191` `char thetitle[_MAXTITLESTRINGLENGTH_];` → `std::string thetitle;`. `source/output_module.h:29` and `output_module.cpp:1181` `const char titles[_MAXTITLESTRINGLENGTH_]` param → `const std::string& titles`. Inside, any `strlen(titles)`/`%s titles`/`strcpy`/`strtok` on these → operate on `std::string` (`.size()`, `.c_str()`); if the code uses `strtok` to split on `_DELIMITER_`, replace with `std::string` find/substr or pass `titles.c_str()` to the existing `get_number_of_titles`. Match original splitting exactly.

- [ ] **Step 5: Convert `get_number_of_titles`**

`tools/common.cpp` `int get_number_of_titles(char* titlestring)` → `int get_number_of_titles(const std::string& titlestring)`; update its body to iterate over the `std::string` (count occurrences of `_DELIMITER_`); update its declaration in the corresponding header and all callers (`grep -rn "get_number_of_titles" source/ tools/ include/`). Callers passing a now-`std::string` pass it directly; callers passing `char*` pass as-is (implicit conversion) — verify each compiles.

- [ ] **Step 6: Verify bit-identical**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 7: Commit**

```bash
git add include/common.h source/background_module.* source/thermodynamics_module.* source/primordial_module.* source/perturbations_module.* source/output_module.* tools/common.cpp
git commit -m "titles: column-title strings to std::string across output_titles + fields (#291)"
```

---

## Task 6: Remaining locals + `perturb_output_firstline_and_ic_suffix`

**Files:**
- Modify: `source/input_module.cpp:65-72` (input_file/precision_file/tmp_file/extension/stringoutput/inifilename), `:2887,3012` (`buf[64]`)
- Modify: `source/perturbations_module.{h,cpp}` (`ic_suffix`/`first_line` params)
- Modify: `source/output_module.cpp:208,544,820-821` (first_line/ic_suffix locals + the snprintf using ic_suffix)
- Modify: `source/primordial_module.cpp:2876-2920` (command builders)
- Modify: `source/thermodynamics_module.cpp:1620,1691` (`line` + fgets)
- Modify: `tools/quadrature.cpp:94` (`method_chosen`)
- Test: `./scripts/verify_bitidentical.sh`

- [ ] **Step 1: `input_module.cpp` arg/temp buffers → `std::string`**

`:65` `char input_file[_ARGUMENT_LENGTH_MAX_];` → `std::string input_file;`; `:66` `precision_file` likewise; `:70` `char tmp_file[tmp_file_size];` → `std::string tmp_file;` (delete the `tmp_file_size` calc at `:67-69`); `:71` `char extension[5];` → `std::string extension;`; `:72` `FileArg stringoutput, inifilename;` → `std::string stringoutput, inifilename;`. Then walk `grep -n "input_file\|precision_file\|tmp_file\|extension\|stringoutput\|inifilename" source/input_module.cpp` and convert each use:
- `snprintf(tmp_file, ..., "%s%02d%s", a, n, b)` → build with `std::string`. For zero-padded integers use `std::ostringstream` (exact, no float-format risk since `n` is an `int`): `std::ostringstream os; os << a << std::setw(2) << std::setfill('0') << n << b; tmp_file = os.str();` (add `#include <sstream>`, `#include <iomanip>`). Do **not** use `ostringstream` for floating-point values — its default formatting differs from `printf` and can break bit-identicality.
- `strcpy`/`strcat` on these → `=`/`+=`.
- `fopen(tmp_file, ...)`/`%s` → `.c_str()`.
- `parser_read_string(&fc_input, "root", &stringoutput, ...)` (`:123`) → `parser_read_string(&fc_input, "root", stringoutput, ...)`.
- `fc_root.set("root", tmp_file)` (`:150`) — `set` takes `const std::string&`, so pass `tmp_file` directly.

- [ ] **Step 2: `buf[64]` locals**

`source/input_module.cpp:2887,3012` `char buf[64];` — inspect each (`grep -n "buf" source/input_module.cpp` near those lines). If used as `snprintf(buf,...,"%d"/"%g", x)` then `strcat`/`set`, replace with `std::string` built via a small numeric scratch + concatenation, matching the original format exactly.

- [ ] **Step 3: `perturb_output_firstline_and_ic_suffix` → `std::string& ic_suffix`**

`source/perturbations_module.h:22` and `source/perturbations_module.cpp:437-440`:
```cpp
int PerturbationsModule::perturb_output_firstline_and_ic_suffix(int index_ic,
                                                                char first_line[_LINE_LENGTH_MAX_],
                                                                std::string& ic_suffix) const {
  first_line[0] = '\0';
  ic_suffix.clear();
  ...
    ic_suffix = "ad";  // each strcpy(ic_suffix, "xx") becomes assignment
```
Keep `first_line` as a parameter for now, but convert it too in this task for consistency: change to `std::string& first_line`, replace `first_line[0]='\0'` → `first_line.clear()`, each `strcpy(first_line, "...")` → `first_line = "...";`. (This also touches the `output_module.cpp:884` caller — see next step — and the wrapper `ic_info` in Task 8.)

- [ ] **Step 4: `output_module.cpp` `first_line`/`ic_suffix` locals**

`:208,544,820` `char first_line[_LINE_LENGTH_MAX_];` → `std::string first_line;`; `:821` `char ic_suffix[4];` → `std::string ic_suffix;`. Every `strcpy(first_line, "literal")` → `first_line = "literal";`. The caller at `:884` now passes two `std::string&` (matches Step 3). The `snprintf(file_name, ..., ic_suffix, ...)` built in Task 4 already concatenates `ic_suffix` as a `std::string`. The `fprintf(..., "%s", first_line)` uses → `first_line.c_str()`.

- [ ] **Step 5: `primordial_module.cpp` command builders**

`:2876-2878` `char arguments[_ARGUMENT_LENGTH_MAX_]`, `char line[_LINE_LENGTH_MAX_]`, `char command_with_arguments[2*_ARGUMENT_LENGTH_MAX_]` → `std::string`. The builder (`:2895` `snprintf(arguments," ")`, `:2900` appends each param, `:2915` `snprintf(command_with_arguments, ..., "%s %s", command, arguments)`) → build with `std::string`: `arguments = " ";` then `arguments += ...;` per param (matching original numeric formatting via a numeric scratch), and `command_with_arguments = std::string(command) + " " + arguments;`. The `line` buffer is read from `popen`/`fgets` — convert that read loop to `std::getline` over the stream if it uses `fgets`, else `.c_str()` at consumers. Inspect `:2876-2940` and match exactly.

- [ ] **Step 6: `thermodynamics_module.cpp` `line` + fgets**

`:1620` `char line[_LINE_LENGTH_MAX_];` and `:1691` `while (fgets(line, _LINE_LENGTH_MAX_ - 1, fA) != nullptr)` — this reads the sBBN file. Convert to `std::string line;` + read via the existing `FILE* fA` using `std::fgets` is not available for `std::string`; simplest behavior-preserving option: keep a local numeric/line scratch ONLY if the parsing uses `sscanf(line, ...)`. **Preferred:** since this loop `sscanf`s fixed columns, replace the `FILE*`+`fgets` with `std::ifstream`+`std::getline(in, line)` and `std::sscanf(line.c_str(), ...)`. Inspect `:1680-1740` and convert the open (`class_open`→`std::ifstream`) and loop together; keep the numeric parsing identical.

- [ ] **Step 7: `quadrature.cpp` `method_chosen`**

`tools/quadrature.cpp:94` `char method_chosen[_METHOD_CHOSEN_SIZE_];` → `std::string method_chosen;`. Convert its `strcpy`/`%s` uses; if `_METHOD_CHOSEN_SIZE_` becomes unused, delete its `#define` (find with `grep -rn "_METHOD_CHOSEN_SIZE_"`).

- [ ] **Step 8: Verify bit-identical**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 9: Commit**

```bash
git add source/input_module.cpp source/perturbations_module.* source/output_module.cpp source/primordial_module.cpp source/thermodynamics_module.cpp tools/quadrature.cpp
git commit -m "Convert remaining char[] locals + ic_suffix/first_line to std::string (#291)"
```

---

## Task 7: Audit — confirm no non-wrapper `char[]`/constant references remain

**Files:** none (verification only)

- [ ] **Step 1: Grep for stragglers (excluding ErrorMsg/varargs/argv families)**

Run:
```bash
grep -rIn "FileName\|FileArg\|_FILENAMESIZE_\|_ARGUMENT_LENGTH_MAX_\|_LINE_LENGTH_MAX_\|_MAXTITLESTRINGLENGTH_" \
  source/ tools/ include/ main/ | grep -v "define _ARGUMENT_LENGTH_MAX_\|define _LINE_LENGTH_MAX_\|define _FILENAMESIZE_\|define _MAXTITLESTRINGLENGTH_\|typedef char FileName\|typedef char FileArg\|sMAXTITLESTRINGLENGTH\|sARGUMENT_LENGTH_MAX"
```
Expected: **no output** (the only remaining references are the definitions themselves in `common.h` and the `ClassConstants` exports, handled in Tasks 8–9).

- [ ] **Step 2: If any straggler appears, convert it** using the patterns from Tasks 1–6, rebuild, `./scripts/verify_bitidentical.sh`, and fold into the appropriate prior commit area (or a small follow-up commit). Do not proceed until the grep is clean.

---

## Task 8: Cython wrapper — `generate_wrapper.py`, regenerate `cclassy.pxd`, `classy.pyx`

**Files:**
- Modify: `generate_wrapper.py:127-129,193,266-272`
- Modify: `classy.pyx:1105-1114,1150-1160,1195-1205,1317-1346,169` (titles + ic_suffix)
- Modify: `main/class.h:33,35` (`ClassConstants`)
- Regenerate: `cclassy.pxd`
- Test: `pip install .` + `python/test_class.py` (`COMPARE_OUTPUT_REF`)

- [ ] **Step 1: Teach the generator about `std::string`; drop dead ctypedefs**

`generate_wrapper.py`:
- Delete `:127` (`ctypedef char FileArg[...]`) and `:129` (`ctypedef char FileName[...]`) `enums.append(...)` lines. Keep `:128` `ErrorMsg`.
- `:193` struct `allowed_types` → add `'std::string'`: `allowed_types = ['double', 'int', 'short', 'std::string']`.
- `:266` method `allowed_types` → add `'std::string'` (and remove `'FileArg'` since the typedef is gone): keep `'ErrorMsg'`, add `'std::string'`.

- [ ] **Step 2: Drop the `ClassConstants` exports that no longer exist**

`main/class.h`: delete `:33` `sMAXTITLESTRINGLENGTH` and `:35` `sARGUMENT_LENGTH_MAX` lines (both constants are deleted in Task 9; `sMAXTITLESTRINGLENGTH` was only used by the titles `resize`, removed below; `sARGUMENT_LENGTH_MAX` is unused in `classy.pyx`, verified). Keep `sFALSE`, `sFAILURE`.

- [ ] **Step 3: Update `classy.pyx` title-reading sites (4)**

At each of the title-reading blocks (background `:1108-1114`, thermo `~:1150`, perturbations-density/`get_perturbations` `~:1195`, and `get_transfer` `:1344`), change:
```python
            string titles
            ...
            titles.resize(constvals.sMAXTITLESTRINGLENGTH)
            status = deref(MOD).MOD_output_titles(<char*> titles.c_str())
            ...
            tmp = <bytes> titles.c_str()
```
to:
```python
            string titles
            ...
            status = deref(MOD).MOD_output_titles(titles)
            ...
            tmp = <bytes> titles
```
(Drop the `resize` line and the `<char*> ... .c_str()` cast — `output_titles` now takes `std::string&` and Cython passes `titles` by reference; `<bytes> titles` converts the filled `std::string` directly.)

- [ ] **Step 4: Update the `ic_suffix` site in `get_transfer`**

`classy.pyx:1320-1372`: change `FileName ic_suffix` (`:1321`) → `string ic_suffix`; the call `perturb_output_firstline_and_ic_suffix(index_ic, ic_info, ic_suffix)` now passes `string& first_line` and `string& ic_suffix`. Change `char ic_info[1024]` (`:1320`) → `string ic_info` (Task 6 made `first_line` a `std::string&`). Replace `ic_key = <bytes> ic_suffix` (`:1372`) → `ic_key = <bytes> ic_suffix` works on `string` (keep), and the `first_line` use (`ic_info`) accordingly. Verify the surrounding `<char*>`/`.c_str()` casts are removed.

- [ ] **Step 5: Regenerate `cclassy.pxd` (do not hand-edit)**

Run: `python generate_wrapper.py`
Expected: `cclassy.pxd` rewritten; `git diff cclassy.pxd` shows `FileArg`/`FileName` ctypedefs gone, `class_dir` and the precision path fields now `string`, `ic_suffix`/`titles`-related method params now `string&`/`string`, `ErrorMsg` unchanged. Confirm it still contains the `precision`/`transfers` structs and all module classes.

- [ ] **Step 6: Build the extension and run the comparison suite**

Run:
```bash
pip install . 2>&1 | tail -5
cd python && COMPARE_OUTPUT_REF=1 python -m pytest test_class.py -x -q 2>&1 | tail -20
```
Expected: `pip install .` succeeds; `COMPARE_OUTPUT_REF` suite passes (no P(k)/Cl regressions). If the suite flags scenarios, investigate before proceeding — wrapper output must match.

- [ ] **Step 7: Commit**

```bash
git add generate_wrapper.py cclassy.pxd classy.pyx main/class.h
git commit -m "Wrapper: std::string for class_dir/ic_suffix/titles; regenerate cclassy.pxd (#291)"
```

---

## Task 9: Delete the typedefs and constants; final verification

**Files:**
- Modify: `include/common.h` (delete `:65` FileName, `:103` `_MAXTITLESTRINGLENGTH_`, `:112` `_LINE_LENGTH_MAX_`, `:114` `_ARGUMENT_LENGTH_MAX_`, `:117` FileArg, and the `_FILENAMESIZE_` define)
- Delete: `scripts/verify_bitidentical.sh`
- Test: full `make class`, `make test-parser`, `pip install .`, `COMPARE_OUTPUT_REF`

- [ ] **Step 1: Confirm zero references before deleting the definitions**

Run:
```bash
grep -rIn "FileName\|FileArg\|_FILENAMESIZE_\|_ARGUMENT_LENGTH_MAX_\|_LINE_LENGTH_MAX_\|_MAXTITLESTRINGLENGTH_" \
  source/ tools/ include/ main/ classy.pyx cclassy.pxd generate_wrapper.py \
  | grep -v "include/common.h:"
```
Expected: **no output** (every remaining hit must be the definition lines in `common.h`).

- [ ] **Step 2: Delete the definitions in `include/common.h`**

Remove the `#define _FILENAMESIZE_ ...`, `#define _LINE_LENGTH_MAX_ ...`, `#define _ARGUMENT_LENGTH_MAX_ ...`, `#define _MAXTITLESTRINGLENGTH_ 8000`, `typedef char FileName[...]` (`:65`), and `typedef char FileArg[...]` (`:117`). Keep `_ERRORMSGSIZE_`, `ErrorMsg`, `_DELIMITER_`.

- [ ] **Step 3: Full rebuild + bit-identical (last run of the harness)**

Run: `./scripts/verify_bitidentical.sh`
Expected: `BIT-IDENTICAL (both gauges)`

- [ ] **Step 4: Parser unit test + wrapper suite**

Run:
```bash
make test-parser && ./test-parser
pip install . 2>&1 | tail -3
cd python && COMPARE_OUTPUT_REF=1 python -m pytest test_class.py -x -q 2>&1 | tail -10; cd ..
```
Expected: parser tests PASS; `pip install .` succeeds; suite passes.

- [ ] **Step 5: Remove the temporary harness and commit**

```bash
git rm scripts/verify_bitidentical.sh
git add include/common.h
git commit -m "Delete FileName/FileArg typedefs and fixed-size string constants (#291)"
```

- [ ] **Step 6: Final grep sanity**

Run:
```bash
grep -rIn "FileName\|FileArg\|_FILENAMESIZE_\|_ARGUMENT_LENGTH_MAX_\|_LINE_LENGTH_MAX_\|_MAXTITLESTRINGLENGTH_" \
  source/ tools/ include/ main/ classy.pyx cclassy.pxd generate_wrapper.py
```
Expected: **no output anywhere.**

---

## Self-review notes (for the executor)

- **The "test" for this refactor is bit-identical output**, not new unit tests — the spec deliberately adds no `k_output_values` regression test (`getline` correctness is structural). Run `./scripts/verify_bitidentical.sh` after every non-wrapper task; never batch two tasks before verifying.
- **`strstr`→`.find`, `strcmp`→`==` must preserve substring-vs-equality semantics** at each site. Re-read the original before converting.
- **Numeric formatting:** integers → `std::to_string` (exact) or `std::ostringstream` with `setw`/`setfill` when zero-padding is needed. Floating-point that must match a `printf` conversion (`%g`/`%e`/`%.16e`) → format with `snprintf` into a small *scalar* numeric scratch (e.g. `char nb[32]`) and concatenate; `ostringstream` default float formatting differs from `printf` and **will** break bit-identicality. The scalar numeric scratch is the one permitted `char` local — never a path/param/list/title.
- **`cclassy.pxd` is generated** — only edit `generate_wrapper.py` and rerun it; never hand-edit the pxd.
- **Constants are deleted only in Task 9**, after Task 7's audit proves no references remain. The build must stay green at every commit.
