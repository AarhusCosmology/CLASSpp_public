# Remove `error_message_` Member Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the shared mutable `error_message_` module member (and the CLASSpp C-struct `error_message` scratch fields) by formatting and throwing error strings from stack locals, closing the issue #38 data race.

**Architecture:** In the C++ build the error macros (`include/common.h`, `#ifdef __cplusplus` block) already throw `std::runtime_error`. We change them to format into a *macro-local* `ErrorMsg` instead of a caller-owned buffer, then mechanically strip the now-dead buffer arguments from ~1180 call sites, convert throwing-member `class_call` wrappers to bare calls, and route the ~28 `tools/`-leaf calls through a new `class_call_failure(function, buf)` macro with a local buffer. The `tools/` numerical layer keeps its C return-code contract. HyRec `.c` and the C macro definitions are untouched.

**Tech Stack:** C++17, CLASS build via `make` (binary) and `pip install .` (classy), pytest (`python/test_class.py` with `COMPARE_OUTPUT_REF`), clang-format.

**Spec:** [`docs/superpowers/specs/2026-06-08-remove-error-message-member-design.md`](../specs/2026-06-08-remove-error-message-member-design.md)

---

## Conventions used in every task

- **Fast build check:** `make class -j` (compiles the C++ binary; ~minutes). A green compile is the primary guard — after the member is removed, *any* stray `error_message_` reference is a hard compile error that names the file and line.
- **Bit-identical check (fast):** this refactor touches no numerics, so output must be byte-identical to the pre-change baseline created in Task 0. Run:
  ```bash
  ./class python/baseline/ref.ini
  diff -rq output/ python/baseline/out_master/ && echo "BIT-IDENTICAL OK"
  ```
- **Authoritative test (slower):** `pip install . && cd python && python -m pytest -q test_class.py test_greybody.py test_hyperspherical.py`. Full reference comparison (`COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -m test_scenario test_class.py`) runs once in Task 9.
- **Formatting:** after any scripted sweep, run `clang-format -i` on changed files (a CI check enforces it). The repo's `.clang-format` governs style.
- **Commit** after each task once its checks pass.
- The helper script `scripts/drop_macro_buffer_arg.py` (created in Task 4) is a token-aware rewriter; it is deleted in Task 9.

---

## File Structure

Files created/modified across the plan:

- Modify: `include/common.h` — the `#ifdef __cplusplus` macro block (Tasks 1, 4, 5).
- Modify: `source/base_module.h` — remove the member + ctor init (Task 6).
- Modify: `source/input_module.h` — remove the member (Task 6).
- Modify: `source/spectra.h`, `include/common.h`, `include/background.h` — remove C-struct `error_message` fields (Task 7).
- Modify (scripted, many): `source/**/*.cpp`, `species/**/*.cpp`, `tools/**/*.cpp` — call-site sweeps (Tasks 2–5, 7).
- Create then delete: `scripts/drop_macro_buffer_arg.py` (Tasks 4, 9).
- Create: `python/baseline/` — pre-change reference output (Task 0; git-ignored, not committed).

---

## Task 0: Establish the bit-identical baseline

**Files:**
- Create: `python/baseline/ref.ini` (a small, fast, representative run)
- Create: `python/baseline/out_master/` (reference output — not committed)

- [ ] **Step 1: Build the current (master) binary**

Run:
```bash
git status   # confirm on branch fix-issue-38-remove-error-message-member, clean
make class -j
```
Expected: builds `class` with no errors.

- [ ] **Step 2: Write a small reference ini that exercises the threaded paths**

Create `python/baseline/ref.ini`:
```ini
output = tCl,pCl,lCl,mPk
lensing = yes
l_max_scalars = 1000
P_k_max_h/Mpc = 10
root = python/baseline/out_master/ref_
write_warnings = yes
```

- [ ] **Step 3: Generate the reference output**

Run:
```bash
mkdir -p python/baseline/out_master
./class python/baseline/ref.ini
ls python/baseline/out_master/
```
Expected: several `ref_*.dat` files written.

- [ ] **Step 4: Record the baseline commit and add baseline dir to .gitignore**

Run:
```bash
echo "python/baseline/" >> .gitignore
git rev-parse HEAD   # note this as the bit-identical reference commit
```

- [ ] **Step 5: Commit the .gitignore change**

```bash
git add .gitignore
git commit -m "test: ignore local bit-identical baseline dir (#38)"
```

---

## Task 1: Throw from a macro-local buffer; add `class_call_failure`

This decouples the throw buffer from any caller-owned buffer **without changing
any macro signature or call site**. After this task the member is no longer the
throw target, so the most acute part of the #38 race is already mitigated, and
the build remains green and bit-identical.

**Files:**
- Modify: `include/common.h` (the `#ifdef __cplusplus` block, ~lines 279–348)

- [ ] **Step 1: Rewrite the C++ macro bodies to use a local buffer**

In `include/common.h`, replace the entire `#ifdef __cplusplus … #endif` override
block with the following. Signatures are unchanged; the `error_message_output`
parameter is simply no longer referenced (a macro arg that is not pasted is not
evaluated, so existing call sites still compile and incur no extra cost). A new
2-arg `class_call_failure` is added for `tools/`-leaf calls.

```cpp
#ifdef __cplusplus

#undef class_call_except
#define class_call_except(function,                                                     \
                          error_message_from_function,                                  \
                          error_message_output,                                         \
                          list_of_commands)                                             \
  {                                                                                     \
    try {                                                                               \
      if ((function) == _FAILURE_) {                                                    \
        ErrorMsg _class_err_;                                                           \
        class_call_message(_class_err_, #function, error_message_from_function);        \
        throw std::runtime_error(_class_err_);                                          \
      }                                                                                 \
    }                                                                                   \
    catch (...) {                                                                       \
      list_of_commands;                                                                 \
      throw;                                                                            \
    }                                                                                   \
  }

#undef class_call
#define class_call(function, error_message_from_function, error_message_output) \
  {                                                                             \
    if ((function) == _FAILURE_) {                                              \
      ErrorMsg _class_err_;                                                     \
      class_call_message(_class_err_, #function, error_message_from_function);  \
      throw std::runtime_error(_class_err_);                                    \
    }                                                                           \
  }

/* New: for tools/ leaf functions that return _FAILURE_ and wrote into buf. */
#define class_call_failure(function, error_message_from_function)              \
  {                                                                            \
    if ((function) == _FAILURE_) {                                             \
      ErrorMsg _class_err_;                                                    \
      class_call_message(_class_err_, #function, error_message_from_function); \
      throw std::runtime_error(_class_err_);                                   \
    }                                                                          \
  }

#undef class_test_except
#define class_test_except(condition, error_message_output, list_of_commands, args, ...) \
  {                                                                                     \
    if (condition) {                                                                    \
      ErrorMsg _class_err_;                                                             \
      class_test_message(_class_err_, #condition, args, ##__VA_ARGS__);                 \
      list_of_commands;                                                                 \
      throw std::runtime_error(_class_err_);                                            \
    }                                                                                   \
  }

#undef class_test
#define class_test(condition, error_message_output, args, ...)         \
  {                                                                    \
    if (condition) {                                                   \
      ErrorMsg _class_err_;                                            \
      class_test_message(_class_err_, #condition, args, ##__VA_ARGS__); \
      throw std::runtime_error(_class_err_);                           \
    }                                                                  \
  }

#undef class_stop
#define class_stop(error_message_output, args, ...)                    \
  {                                                                    \
    ErrorMsg _class_err_, _class_opt_args_;                            \
    class_protect_sprintf(_class_opt_args_, args, ##__VA_ARGS__);      \
    class_build_error_string(_class_err_, "error; %s", _class_opt_args_); \
    throw std::runtime_error(_class_err_);                             \
  }

#undef class_open
#define class_open(pointer, filename, mode, error_output)                    \
  {                                                                          \
    pointer = fopen(filename, mode);                                         \
    if (pointer == nullptr) {                                                \
      ErrorMsg _class_err_;                                                  \
      class_build_error_string(_class_err_,                                  \
                               "could not open %s with name %s and mode %s", \
                               #pointer,                                     \
                               filename,                                     \
                               #mode);                                       \
      throw std::runtime_error(_class_err_);                                 \
    }                                                                        \
  }

#endif
```

- [ ] **Step 2: Build**

Run: `make class -j`
Expected: compiles with no errors (signatures unchanged, all call sites still valid).

- [ ] **Step 3: Bit-identical check**

Run:
```bash
./class python/baseline/ref.ini
diff -rq output/ python/baseline/out_master/ 2>/dev/null; \
  diff -q <(ls output) <(ls python/baseline/out_master) ; \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo "BIT-IDENTICAL OK"
```
Expected: `BIT-IDENTICAL OK` (no diffs). *(Note: the run writes to `root` from the ini; adjust `root` to `output/ref_` for this comparison or compare against the same `root`.)*

- [ ] **Step 4: Run the fast test subset**

Run: `pip install . -q && cd python && python -m pytest -q test_class.py -k "test_parameters or test_input" ; cd ..`
Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add include/common.h
git commit -m "Throw error macros from a local buffer; add class_call_failure (#38)"
```

---

## Task 2: Route `tools/`-leaf `class_call` sites through `class_call_failure`

The ~28 sites where `error_message_` (or a struct field) is passed *into* a
`tools/` leaf function that returns `_FAILURE_`. These must keep a real buffer
for the leaf to write into — a stack local.

**Files:**
- Modify: `source/*.cpp`, `species/*.cpp` (leaf-call sites only)

- [ ] **Step 1: Enumerate the leaf function names (functions that take an ErrorMsg out-param)**

Run:
```bash
grep -rhoE "\b[a-z_]+[a-z0-9_]*\s*\([^;]*ErrorMsg[^;]*\)" tools/*.h \
  | grep -oE "\b[a-z_]+[a-z0-9_]*\s*\(" | tr -d ' (' | sort -u > /tmp/leaf_names.txt
cat /tmp/leaf_names.txt
```
Expected: a list like `array_spline_table_lines`, `array_interpolate_spline`,
`hyperspherical_HIS_free`, `array_integrate_all_spline`, etc.

- [ ] **Step 2: Find every module call site that passes a buffer INTO one of those leaves**

Run:
```bash
grep -rnE "class_call\(" source/ species/ --include="*.cpp" \
  | grep -Ef <(sed 's/.*/\\b&\\s*(/' /tmp/leaf_names.txt) > /tmp/leaf_sites.txt
wc -l /tmp/leaf_sites.txt ; cat /tmp/leaf_sites.txt
```
Expected: ~28 sites. Inspect the list; these are the only ones converted in this task.

- [ ] **Step 3: Convert each leaf site to a local buffer + `class_call_failure`**

For each site, transform (example from `perturbations_module.cpp:816`):

```cpp
// BEFORE
class_call(array_spline_table_lines(ln_tau_.data(),
                                    ln_tau_size_,
                                    late_sources_[index_md][...],
                                    k_size_[index_md],
                                    ddlate_sources_[index_md][...].data(),
                                    _SPLINE_EST_DERIV_,
                                    error_message_),
           error_message_,
           error_message_);

// AFTER
ErrorMsg buf;
class_call_failure(array_spline_table_lines(ln_tau_.data(),
                                            ln_tau_size_,
                                            late_sources_[index_md][...],
                                            k_size_[index_md],
                                            ddlate_sources_[index_md][...].data(),
                                            _SPLINE_EST_DERIV_,
                                            buf),
                   buf);
```

Rules:
- Replace the leaf's trailing `error_message_` argument **and** the macro's two
  trailing buffer args with a single local `buf`.
- Declare `ErrorMsg buf;` in the smallest enclosing scope (often just above the
  call; inside a loop body is fine — it is cheap stack scratch).
- If a site passes a struct field (e.g. `psp->error_message`) into the leaf,
  treat it identically (local `buf`).

- [ ] **Step 4: Build**

Run: `make class -j`
Expected: compiles cleanly.

- [ ] **Step 5: Bit-identical + format**

Run:
```bash
clang-format -i $(git diff --name-only)
make class -j && ./class python/baseline/ref.ini && \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo OK
```
Expected: `OK`.

- [ ] **Step 6: Commit**

```bash
git add source/ species/
git commit -m "Route tools/-leaf class_call sites through class_call_failure with local buffers (#38)"
```

---

## Task 3: Convert remaining `class_call` member sites to bare calls

After Task 2, every remaining `class_call(...)` in `source/`/`species/` wraps a
**throwing** member (or cross-module member) call. Those propagate exceptions on
their own, so the wrapper is a no-op and is removed.

**Files:**
- Modify: `source/*.cpp`, `species/*.cpp`

- [ ] **Step 1: Sanity-check that no leaf sites remain**

Run:
```bash
grep -rnE "class_call\(" source/ species/ --include="*.cpp" \
  | grep -Ef <(sed 's/.*/\\b&\\s*(/' /tmp/leaf_names.txt)
```
Expected: **no output** (all leaf sites handled in Task 2). If any appear, convert them per Task 2 first.

- [ ] **Step 2: Rewrite `class_call(EXPR, error_message_, error_message_);` → `EXPR;`**

Because `EXPR` can span multiple lines and contain nested parens/commas, use a
token-aware unwrap. Run this Python one-shot over the module sources:

```bash
python3 - <<'PY'
import re, glob
# Unwrap class_call(EXPR, A, B);  ->  EXPR;
# where the two trailing args are error-buffer tokens we are removing.
BUF = re.compile(r'(error_message_|[A-Za-z_][A-Za-z0-9_]*->error_message_?|[A-Za-z_][A-Za-z0-9_]*\.error_message_?)$')
def split_top(s):
    args, depth, cur, i = [], 0, [], 0
    instr=False; inchr=False; esc=False
    for ch in s:
        if esc: cur.append(ch); esc=False; continue
        if ch=='\\': cur.append(ch); esc=True; continue
        if instr:
            cur.append(ch); instr = ch!='"'; continue
        if inchr:
            cur.append(ch); inchr = ch!="'"; continue
        if ch=='"': instr=True; cur.append(ch); continue
        if ch=="'": inchr=True; cur.append(ch); continue
        if ch in '([{': depth+=1
        elif ch in ')]}': depth-=1
        if ch==',' and depth==0:
            args.append(''.join(cur)); cur=[]; continue
        cur.append(ch)
    args.append(''.join(cur))
    return args
for path in glob.glob('source/**/*.cpp', recursive=True)+glob.glob('species/**/*.cpp', recursive=True):
    src=open(path).read(); out=[]; i=0; n=len(src); changed=False
    while i<n:
        m=re.compile(r'\bclass_call\s*\(').match(src, i)
        if not m:
            out.append(src[i]); i+=1; continue
        # find matching close paren
        depth=0; j=m.end()-1; instr=inchr=esc=False
        while j<n:
            ch=src[j]
            if esc: esc=False
            elif ch=='\\': esc=True
            elif instr: instr= ch!='"'
            elif inchr: inchr= ch!="'"
            elif ch=='"': instr=True
            elif ch=="'": inchr=True
            elif ch=='(': depth+=1
            elif ch==')':
                depth-=1
                if depth==0: break
            j+=1
        inner=src[m.end():j]
        args=split_top(inner)
        if len(args)==3 and BUF.search(args[1].strip()) and BUF.search(args[2].strip()):
            out.append(args[0].strip()); i=j+1; changed=True
        else:
            out.append(src[i:m.end()]); i=m.end()  # leave untouched (e.g. except/try forms)
    if changed: open(path,'w').write(''.join(out))
PY
```

- [ ] **Step 3: Verify zero `class_call(` left in modules (except the cleanup variant)**

Run:
```bash
grep -rnE "\bclass_call\(" source/ species/ --include="*.cpp"
```
Expected: **no output**. (`class_call_except`/`class_call_failure` are different tokens and are fine.)

- [ ] **Step 4: Build + format**

Run:
```bash
clang-format -i $(git diff --name-only)
make class -j
```
Expected: compiles cleanly. If a residual `class_call(` failed to unwrap (e.g. unusual layout), the compiler/grep pinpoints it — fix by hand.

- [ ] **Step 5: Bit-identical**

Run:
```bash
./class python/baseline/ref.ini && \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo OK
```
Expected: `OK`.

- [ ] **Step 6: Commit**

```bash
git add source/ species/
git commit -m "Convert throwing-member class_call wrappers to bare calls (#38)"
```

---

## Task 4: Drop the buffer argument from `class_test`/`class_stop`/`class_open`

These macros are shared by `source/`, `species/`, **and** `tools/`. Changing
their signatures forces a sweep of all `.cpp` call sites (372 + 124 + a few). In
`tools/` this drops only the *throw-buffer* argument from `class_test`/`class_stop`;
the manual `return _FAILURE_` paths (which sprintf into the `ErrorMsg` param) are
untouched, so the `tools/` return-code contract is preserved.

**Files:**
- Create: `scripts/drop_macro_buffer_arg.py`
- Modify: `include/common.h`
- Modify: `source/**/*.cpp`, `species/**/*.cpp`, `tools/**/*.cpp`

- [ ] **Step 1: Create the token-aware arg-dropper script**

Create `scripts/drop_macro_buffer_arg.py`:
```python
#!/usr/bin/env python3
"""Remove the positional argument at INDEX from every MACRO(...) invocation.
Token-aware: respects nested parens, string and char literals (so commas inside
format strings are not mistaken for argument separators). Handles multi-line
invocations. Usage: drop_macro_buffer_arg.py MACRO INDEX file1.cpp [file2.cpp ...]
"""
import re, sys

def split_top(s):
    args, depth, cur = [], 0, []
    instr = inchr = esc = False
    for ch in s:
        if esc: cur.append(ch); esc = False; continue
        if ch == '\\': cur.append(ch); esc = True; continue
        if instr: cur.append(ch); instr = ch != '"'; continue
        if inchr: cur.append(ch); inchr = ch != "'"; continue
        if ch == '"': instr = True; cur.append(ch); continue
        if ch == "'": inchr = True; cur.append(ch); continue
        if ch in '([{': depth += 1
        elif ch in ')]}': depth -= 1
        if ch == ',' and depth == 0:
            args.append(''.join(cur)); cur = []; continue
        cur.append(ch)
    args.append(''.join(cur))
    return args

def process(macro, index, text):
    pat = re.compile(r'\b' + re.escape(macro) + r'\s*\(')
    out, i, n, changed = [], 0, len(text), False
    while i < n:
        m = pat.search(text, i)
        if not m:
            out.append(text[i:]); break
        out.append(text[i:m.end()])
        depth, j = 0, m.end() - 1
        instr = inchr = esc = False
        while j < n:
            ch = text[j]
            if esc: esc = False
            elif ch == '\\': esc = True
            elif instr: instr = ch != '"'
            elif inchr: inchr = ch != "'"
            elif ch == '"': instr = True
            elif ch == "'": inchr = True
            elif ch == '(': depth += 1
            elif ch == ')':
                depth -= 1
                if depth == 0: break
            j += 1
        inner = text[m.end():j]
        args = split_top(inner)
        if len(args) > index:
            del args[index]
            out.append(','.join(args)); changed = True
        else:
            out.append(inner)
        out.append(')')
        i = j + 1
    return ''.join(out), changed

def main():
    macro, index, files = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    for path in files:
        text = open(path).read()
        new, changed = process(macro, index, text)
        if changed:
            open(path, 'w').write(new)
            print("rewrote", path)

if __name__ == '__main__':
    main()
```

- [ ] **Step 2: Drop `class_test` arg 1 (the `error_message_output`)**

Run:
```bash
python3 scripts/drop_macro_buffer_arg.py class_test 1 \
  $(grep -rlE "\bclass_test\(" source/ species/ tools/ --include="*.cpp")
```
Expected: prints `rewrote <file>` for each affected file.

- [ ] **Step 3: Drop `class_test_except` arg 1**

Run:
```bash
python3 scripts/drop_macro_buffer_arg.py class_test_except 1 \
  $(grep -rlE "\bclass_test_except\(" source/ species/ tools/ --include="*.cpp")
```
Expected: rewrites any `class_test_except` sites (signature `(condition, error_message_output, list_of_commands, args, ...)` → arg index 1 removed).

- [ ] **Step 4: Drop `class_stop` arg 0**

Run:
```bash
python3 scripts/drop_macro_buffer_arg.py class_stop 0 \
  $(grep -rlE "\bclass_stop\(" source/ species/ tools/ --include="*.cpp")
```

- [ ] **Step 5: Drop `class_open` arg 3 (the `error_output`)**

Run:
```bash
python3 scripts/drop_macro_buffer_arg.py class_open 3 \
  $(grep -rlE "\bclass_open\(" source/ species/ tools/ --include="*.cpp")
```

- [ ] **Step 6: Update the macro signatures in `include/common.h`**

In the `#ifdef __cplusplus` block, change the signatures to drop the buffer
parameter (bodies already use `_class_err_` from Task 1):

```cpp
#undef class_test_except
#define class_test_except(condition, list_of_commands, args, ...)      \
  {                                                                    \
    if (condition) {                                                   \
      ErrorMsg _class_err_;                                            \
      class_test_message(_class_err_, #condition, args, ##__VA_ARGS__); \
      list_of_commands;                                                \
      throw std::runtime_error(_class_err_);                           \
    }                                                                  \
  }

#undef class_test
#define class_test(condition, args, ...)                              \
  {                                                                   \
    if (condition) {                                                  \
      ErrorMsg _class_err_;                                           \
      class_test_message(_class_err_, #condition, args, ##__VA_ARGS__); \
      throw std::runtime_error(_class_err_);                          \
    }                                                                 \
  }

#undef class_stop
#define class_stop(args, ...)                                         \
  {                                                                   \
    ErrorMsg _class_err_, _class_opt_args_;                           \
    class_protect_sprintf(_class_opt_args_, args, ##__VA_ARGS__);     \
    class_build_error_string(_class_err_, "error; %s", _class_opt_args_); \
    throw std::runtime_error(_class_err_);                            \
  }

#undef class_open
#define class_open(pointer, filename, mode)                          \
  {                                                                  \
    pointer = fopen(filename, mode);                                 \
    if (pointer == nullptr) {                                        \
      ErrorMsg _class_err_;                                          \
      class_build_error_string(_class_err_,                          \
                               "could not open %s with name %s and mode %s", \
                               #pointer, filename, #mode);           \
      throw std::runtime_error(_class_err_);                         \
    }                                                                \
  }
```

- [ ] **Step 7: Build + format**

Run:
```bash
clang-format -i $(git diff --name-only --diff-filter=M | grep -E '\.(cpp|h)$')
make class -j
```
Expected: compiles cleanly. Arity mismatches (a missed site, or a `class_stop`/`class_test` whose first/second arg was *not* a buffer) become compile errors that name the line — fix by hand and rebuild.

- [ ] **Step 8: Bit-identical**

Run:
```bash
./class python/baseline/ref.ini && \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo OK
```
Expected: `OK`.

- [ ] **Step 9: Commit**

```bash
git add include/common.h source/ species/ tools/ scripts/drop_macro_buffer_arg.py
git commit -m "Drop throw-buffer arg from class_test/class_stop/class_open macros and call sites (#38)"
```

---

## Task 5: Convert the 6 `class_call_except` sites to the buffer-less form

`class_call_except` has only 6 uses (thermodynamics + primordial), each wrapping
a throwing member with cleanup commands.

**Files:**
- Modify: `include/common.h`
- Modify: `source/thermodynamics_module.cpp` (2 sites), `source/primordial_module.cpp` (4 sites)

- [ ] **Step 1: Update the macro signature to drop both buffer args**

In `include/common.h`:
```cpp
#undef class_call_except
#define class_call_except(function, list_of_commands)                 \
  {                                                                   \
    try {                                                             \
      (function);                                                     \
    }                                                                 \
    catch (...) {                                                     \
      list_of_commands;                                              \
      throw;                                                          \
    }                                                                 \
  }
```
(The member throws on its own, so the explicit `== _FAILURE_` check is no longer
needed; the `try`/cleanup/`throw` is what matters.)

- [ ] **Step 2: Rewrite the 6 call sites**

Find them: `grep -rn "class_call_except" source/ --include="*.cpp"`.
For each, drop the two error-buffer arguments, keeping `function` and the
`list_of_commands` cleanup. Example (`primordial_module.cpp:288`):

```cpp
// BEFORE
class_call_except(primordial_analytic_spectrum_init(),
                  error_message_,
                  error_message_,
                  primordial_free(); /* cleanup */);
// AFTER
class_call_except(primordial_analytic_spectrum_init(),
                  primordial_free(); /* cleanup */);
```
Inspect each site to preserve its exact cleanup statement(s).

- [ ] **Step 3: Build + format**

Run:
```bash
clang-format -i source/thermodynamics_module.cpp source/primordial_module.cpp include/common.h
make class -j
```
Expected: compiles cleanly.

- [ ] **Step 4: Bit-identical**

Run:
```bash
./class python/baseline/ref.ini && \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo OK
```
Expected: `OK`.

- [ ] **Step 5: Commit**

```bash
git add include/common.h source/thermodynamics_module.cpp source/primordial_module.cpp
git commit -m "Convert class_call_except to buffer-less form (#38)"
```

---

## Task 6: Remove the `error_message_` member

At this point no macro uses a caller buffer, so the member should be referenced
only by its own declaration/initialization. The compiler will flag any straggler.

**Files:**
- Modify: `source/base_module.h`
- Modify: `source/input_module.h`

- [ ] **Step 1: Confirm remaining references**

Run: `grep -rn "error_message_\b" source/ species/ --include="*.cpp" --include="*.h" | grep -v "error_message_from_function"`
Expected: only the declarations in `base_module.h` and `input_module.h` (and possibly the ctor init line). If any `.cpp` reference remains, it is a missed site — convert it (member call → bare; leaf → `class_call_failure`; test/stop → drop arg) before continuing.

- [ ] **Step 2: Remove the member from `base_module.h`**

In `source/base_module.h`:
- Delete `mutable ErrorMsg error_message_;` (line ~38).
- Delete the `error_message_[0] = '\n';` line from the constructor body (line ~34).

- [ ] **Step 3: Remove the member from `input_module.h`**

In `source/input_module.h`: delete `ErrorMsg error_message_;` (line ~57). If
`input_module.cpp`'s constructor initializes it, remove that too (grep to check).

- [ ] **Step 4: Build**

Run: `make class -j`
Expected: compiles cleanly. Any remaining `error_message_` use is now a hard
error naming the exact file:line — fix each and rebuild until clean.

- [ ] **Step 5: Bit-identical**

Run:
```bash
./class python/baseline/ref.ini && \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo OK
```
Expected: `OK`.

- [ ] **Step 6: Commit**

```bash
git add source/base_module.h source/input_module.h source/input_module.cpp
git commit -m "Remove the mutable error_message_ module member (#38)"
```

---

## Task 7: Remove the CLASSpp C-struct `error_message` fields

Three CLASSpp struct fields (`spectra`, `precision`, `background`) used as scratch
at 6 sites. `generic_integrator_workspace::error_message` (dei_rkck) is **kept**.

**Files:**
- Modify: `source/spectra.h`, `include/common.h`, `include/background.h`
- Modify: the ~6 `.cpp` call sites (`psp->error_message`, `ppr->error_message`, `pba->error_message`)

- [ ] **Step 1: List the field uses to convert**

Run:
```bash
grep -rnE "(psp|ppr|pba)->error_message\b|\bgi\.error_message\b" source/ tools/ --include="*.cpp" | grep -v "gi.error_message"
```
Expected: ~6 lines (the `gi.error_message` dei_rkck ones are excluded and stay).
Note the exact files/lines.

- [ ] **Step 2: Convert each site to a local buffer**

For each of the ~6 sites, the field is passed into a `tools/` leaf or used as a
throw buffer. Convert exactly as in Task 2: introduce `ErrorMsg buf;` in scope
and pass `buf` (wrapping with `class_call_failure` if it feeds a leaf that
returns `_FAILURE_`). Example (`output_module.cpp:73`, `psp->error_message`
passed into a spectra leaf):
```cpp
// BEFORE
class_call_failure(spectra_pk_at_z(..., psp->error_message), psp->error_message);
// AFTER
ErrorMsg buf;
class_call_failure(spectra_pk_at_z(..., buf), buf);
```

- [ ] **Step 3: Remove the field declarations**

- `source/spectra.h:50` — delete `ErrorMsg error_message;` (and its doc comment).
- `include/common.h:1080` — delete `ErrorMsg error_message;` from the `precision`
  struct (and the surrounding `@name - zone for writing error messages` doxygen
  block).
- `include/background.h` — locate and delete the `background` struct's
  `ErrorMsg error_message;` field. (Confirm its exact declaration; the spec flags
  that grep matched it indirectly via the single `pba->error_message` use.)

- [ ] **Step 4: Build**

Run: `make class -j`
Expected: compiles cleanly. A straggler `->error_message` use becomes a compile
error naming the line — fix and rebuild.

- [ ] **Step 5: Bit-identical**

Run:
```bash
./class python/baseline/ref.ini && \
  for f in python/baseline/out_master/ref_*; do diff -q "$f" "output/$(basename $f)"; done && echo OK
```
Expected: `OK`.

- [ ] **Step 6: Commit**

```bash
git add source/ include/ tools/
git commit -m "Remove spectra/precision/background C-struct error_message fields (#38)"
```

---

## Task 8: Verify the error path itself still produces good messages

A pure compile + bit-identical pass does not exercise the *error* path (errors
are exceptional). Confirm a deliberately triggered failure throws a well-formed,
complete `std::runtime_error` (not a garbled/empty buffer).

**Files:** none modified.

- [ ] **Step 1: Trigger a known input error via classy and inspect the message**

Run:
```bash
pip install . -q
python3 - <<'PY'
from classy import Class, CosmoComputationError, CosmoSevereError
c = Class()
c.set({'h': -1.0})  # invalid Hubble parameter -> should raise with a clear message
try:
    c.compute()
    print("NO ERROR RAISED — unexpected")
except (CosmoComputationError, CosmoSevereError) as e:
    msg = str(e)
    print("RAISED:", repr(msg[:300]))
    assert msg.strip(), "empty error message"
    assert "\x00" not in msg, "embedded NUL in message"
    print("ERROR-PATH OK")
PY
```
Expected: `ERROR-PATH OK` with a readable, non-empty message naming the failing
condition/function.

- [ ] **Step 2 (optional): Threaded error sanity**

If feasible, run a scenario known to fail inside a threaded module (e.g. a
pathological `P_k_max`/`l_max` combination) with `number_of_threads > 1` and
confirm a single clean exception surfaces from `future.get()` without crashing.
Document the result; this directly validates the #38 fix.

- [ ] **Step 3: Commit (if any test scaffolding was added; otherwise skip)**

---

## Task 9: Full verification across all build manifests + finish

**Files:**
- Modify: `setup.py`, `CLASS.xcodeproj/...` only if the build lists changed (no
  files were added/removed in this refactor, so expect no manifest edits).
- Delete: `scripts/drop_macro_buffer_arg.py`

- [ ] **Step 1: Confirm no manifest changes needed**

This refactor adds/removes no source files, so `setup.py`, the Makefile object
lists, and `CLASS.xcodeproj` need no edits. Verify:
```bash
git diff --name-only origin/master | grep -E "setup.py|Makefile|xcodeproj" || echo "no manifest changes (expected)"
```

- [ ] **Step 2: Clean build of binary + classy**

Run:
```bash
make clean && make -j        # builds class + classy
```
Expected: both build with no errors.

- [ ] **Step 3: Full nose_tests suite**

Run:
```bash
cd python && python -m pytest -q test_class.py test_greybody.py test_hyperspherical.py ; cd ..
```
Expected: all pass.

- [ ] **Step 4: Reference-comparison test (the CI gate)**

Build a `master` reference as `classyref` (per `.github/workflows/test_on_pull_request.yml`)
and run:
```bash
cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -q -m test_scenario test_class.py ; cd ..
```
Expected: pass with no output diffs (bit-identical refactor; no `COMPARE_OUTPUT_REF`
ref regeneration needed).

- [ ] **Step 5: Remove the helper script**

Run:
```bash
git rm scripts/drop_macro_buffer_arg.py
git commit -m "Remove one-shot macro-arg rewrite helper (#38)"
```

- [ ] **Step 6: Push and open the PR**

Run:
```bash
git push -u origin fix-issue-38-remove-error-message-member
gh pr create --base master \
  --title "Remove error_message_ member; throw from locals (closes #38)" \
  --body "$(cat <<'EOF'
Closes #38.

Removes the shared `mutable ErrorMsg error_message_` module member (and the
CLASSpp C-struct `error_message` scratch fields), the source of the documented
thread-safety data race: concurrent thread-pool tasks shared a single 2048-byte
format buffer when throwing. The C++ error macros now format into a macro-local
`ErrorMsg` and throw, mirroring the existing `tools/exceptions.cpp` pattern, so
the buffer is per-call and race-free.

- `class_call` around throwing members → bare calls.
- New `class_call_failure(function, buf)` for the ~28 `tools/`-leaf calls that
  return `_FAILURE_`; callers pass a stack-local buffer.
- `class_test`/`class_stop`/`class_open`/`class_call_except` drop their buffer
  argument.
- `tools/` numerical layer (incl. `generic_integrator`/dei_rkck) keeps its C
  return-code contract; HyRec `.c` and the C macros are untouched.

Pure error-path refactor: output is bit-identical to master, verified with the
reference-comparison suite.
EOF
)"
```

- [ ] **Step 7: Update memory**

Add/update a project memory note recording: #38 root cause (shared throw buffer
in thread-pool tasks; symptom = garbled message, not crash), the fix (throw from
locals + `class_call_failure`), what was intentionally left on the C contract
(`tools/` leaves + dei_rkck `generic_integrator`), and the PR number.

---

## Self-Review

**Spec coverage:**
- Macro layer changes (class_test/stop/open/call_except/call_failure) → Tasks 1, 4, 5. ✓
- ~1455 member sites → bare calls → Task 3. ✓
- ~28 leaf-boundary sites → `class_call_failure` + local → Task 2. ✓
- class_test/class_stop arg drop → Task 4. ✓
- 6 CLASSpp struct-field sites + field removal → Task 7. ✓
- Member removal (base_module.h, input_module.h) → Task 6. ✓
- tools/ + dei_rkck + HyRec left untouched → stated in Tasks 2/4/7; dei_rkck explicitly excluded in Task 7. ✓
- Verification: 3 manifests, classy, nose_tests, COMPARE_OUTPUT_REF, error-path → Tasks 8, 9. ✓

**Placeholder scan:** No TBD/TODO; every code/script step shows concrete content. The one acknowledged unknown (exact `background` struct field declaration) is handled procedurally in Task 7 Step 3 with a compile-driven fallback, matching the spec's flagged risk.

**Type/name consistency:** `_class_err_` local name is used consistently across all macro bodies (Tasks 1, 4, 5). `class_call_failure(function, error_message_from_function)` signature is defined in Task 1 and used with a local `buf` in Tasks 2 and 7. `scripts/drop_macro_buffer_arg.py` is created in Task 4 and removed in Task 9. The leaf-name inventory `/tmp/leaf_names.txt` is built in Task 2 and reused in Task 3.
