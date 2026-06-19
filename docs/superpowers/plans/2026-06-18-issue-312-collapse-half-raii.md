# Issue #312 — Collapse the half-RAII pattern — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the owning `std::vector` the single member across `perturb_vector`, `perturb_workspace`, and the ndf15 `jacobian` / `numjac_workspace`, deleting the raw-pointer twins — with byte-identical output.

**Architecture:** For each twin, the backing vector is renamed to the bare member name and the raw pointer (plus its `raw = vec.data()` aliasing line in the init function) is deleted. Indexing (`m[i]`, `m[i][j]`) keeps working through `std::vector::operator[]`; sites that pass a member where a `double*`/`int*`/`double**` is expected get `.data()`. For the three 2D row-pointer arrays we keep the contiguous `*_data_vec` flat backing and rename only the `*_rows_vec` (`std::vector<double*>`), so `m[i][j]` is untouched. **No arithmetic changes ⇒ bit-identical results; the compiler enumerates the `.data()` sites and a captured reference run is the oracle.**

**Tech Stack:** C++17; CMake (in-source build, generated `Makefile` + `./class` binary); spec at `docs/superpowers/specs/2026-06-18-issue-312-collapse-half-raii-design.md`.

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `source/perturbations.h` | `perturb_vector`, `perturb_workspace` struct defs | replace raw+storage twins with single vectors |
| `source/perturbations_module.cpp` | sizing/wiring + consumers of `pvec*`, `pv->y/dy/used_in_sources` | drop aliasing lines; `.data()` at pointer-passing sites |
| `include/evolver_ndf15.h` | `jacobian`, `numjac_workspace` struct defs | replace raw+backing twins with single vectors |
| `tools/evolver_ndf15.cpp` | init functions + all evolver call sites | drop aliasing lines; `.data()` at pointer-passing sites |

The `jacobian` / `numjac_workspace` structs are evolver-local (declared in `include/evolver_ndf15.h`, used only in `tools/evolver_ndf15.cpp`). The perturbation structs are used across `source/perturbations_module.cpp`.

**Verification model used by every task:** the byte-identical reference captured in Task 0 (`/tmp/ref312/`) is the oracle. After each task, rebuild `./class`, re-run the same ini, and `diff` the `.dat` outputs against the reference — they must be identical. A non-empty diff means the refactor changed behavior and must be fixed before committing.

---

### Task 0: Capture the byte-identical reference (oracle)

The branch currently equals master `f0939a2b` plus the spec commit (which does not affect the build), so building now yields the true pre-refactor reference.

**Files:** none (creates `/tmp/ref312/`).

- [ ] **Step 1: Write the reference ini**

Create `/tmp/ref312.ini` with exactly:

```ini
output = tCl,pCl,lCl,mPk
lensing = yes
root = /tmp/ref312/run_
```

(Plain ΛCDM with lensed Cls + matter P(k). The default scalar perturbation vector has >15 equations, so this exercises the dense **and** sparse `jacobian` members, `numjac_workspace`, `perturb_vector`, and `perturb_workspace`.)

- [ ] **Step 2: Build the `class` binary at the current (pre-refactor) tree**

Run:
```bash
mkdir -p /tmp/ref312
cmake --build . --target class -j 2>&1 | tail -5
ls -la ./class
```
Expected: build succeeds; `./class` exists. (If CMake is not yet configured in-tree, run `cmake .` first.)

- [ ] **Step 3: Produce the reference output**

Run:
```bash
./class /tmp/ref312.ini
ls /tmp/ref312/
```
Expected: files `run_cl.dat`, `run_cl_lensed.dat`, `run_pk.dat` (plus `run_parameters.ini`) are written.

- [ ] **Step 4: Snapshot the reference data files**

Run:
```bash
cp /tmp/ref312/run_cl.dat        /tmp/ref312/REF_cl.dat
cp /tmp/ref312/run_cl_lensed.dat /tmp/ref312/REF_cl_lensed.dat
cp /tmp/ref312/run_pk.dat        /tmp/ref312/REF_pk.dat
echo "reference captured"
```
Expected: `reference captured`. These `REF_*` files are the oracle for Tasks 1–5.

- [ ] **Step 5: Define the compare helper (used verbatim in later tasks)**

The byte-compare step in every task is:
```bash
./class /tmp/ref312.ini >/dev/null \
 && diff /tmp/ref312/REF_cl.dat        /tmp/ref312/run_cl.dat \
 && diff /tmp/ref312/REF_cl_lensed.dat /tmp/ref312/run_cl_lensed.dat \
 && diff /tmp/ref312/REF_pk.dat        /tmp/ref312/run_pk.dat \
 && echo "BYTE-IDENTICAL"
```
Expected (after each refactor task): `BYTE-IDENTICAL` with no `diff` output.

---

### Task 1: Collapse `perturb_vector` (y, dy, used_in_sources)

**Files:**
- Modify: `source/perturbations.h:298-306`
- Modify: `source/perturbations_module.cpp:3047-3053` (wiring) + `.data()` sites flagged by the compiler

- [ ] **Step 1: Replace the struct members**

In `source/perturbations.h`, replace the current block:
```cpp
  double* y;  /**< vector of perturbations to be integrated */
  double* dy; /**< time-derivative of the same vector */
  std::vector<double> y_storage;
  std::vector<double> dy_storage;

  int* used_in_sources; /**< boolean array specifying which
                           perturbations enter in the calculation of
                           source functions */
  std::vector<int> used_in_sources_storage;
```
with:
```cpp
  std::vector<double> y;  /**< vector of perturbations to be integrated */
  std::vector<double> dy; /**< time-derivative of the same vector */

  std::vector<int> used_in_sources; /**< boolean array specifying which
                           perturbations enter in the calculation of
                           source functions */
```

- [ ] **Step 2: Drop the aliasing lines in the wiring**

In `source/perturbations_module.cpp`, replace:
```cpp
  ppv->y_storage.assign(ppv->pt_size, 0.0);
  ppv->dy_storage.resize(ppv->pt_size);
  ppv->used_in_sources_storage.resize(ppv->pt_size);
  ppv->y               = ppv->y_storage.data();
  ppv->dy              = ppv->dy_storage.data();
  ppv->used_in_sources = ppv->used_in_sources_storage.data();
```
with:
```cpp
  ppv->y.assign(ppv->pt_size, 0.0);
  ppv->dy.resize(ppv->pt_size);
  ppv->used_in_sources.resize(ppv->pt_size);
```

- [ ] **Step 3: Build — the compiler lists the pointer-passing sites**

Run:
```bash
cmake --build . --target class -j 2>&1 | grep -E "error:|perturbations" | head -40
```
Expected: errors of the form "no viable conversion from `std::vector<double>` to `double*`" at each site that passes `ppv->y` / `pv->y` / `pv->dy` / `pv->used_in_sources` as a pointer argument (e.g. the `generic_evolver` call, `memcpy`/copy helpers, and any `double* … = pv->y;` locals). Pure `pv->y[i]` indexing produces **no** error.

- [ ] **Step 4: Add `.data()` at each flagged site**

For every compiler-flagged site, append `.data()` to the offending member access (e.g. `pv->y` → `pv->y.data()`, `ppv->used_in_sources` → `ppv->used_in_sources.data()`). Do **not** touch `[i]` indexing sites. Re-run the build until clean:
```bash
cmake --build . --target class -j 2>&1 | tail -5
```
Expected: build succeeds with no errors.

- [ ] **Step 5: Byte-compare against the reference**

Run the Task 0 Step 5 compare helper.
Expected: `BYTE-IDENTICAL`, no `diff` output.

- [ ] **Step 6: Commit**

```bash
git add source/perturbations.h source/perturbations_module.cpp
git commit -m "$(printf '%s\n' 'refactor(#312): collapse perturb_vector half-RAII twins' '' 'y/dy/used_in_sources become single std::vector members; raw' 'pointers and their .data() aliasing dropped. Byte-identical.' '' 'Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 2: Collapse `perturb_workspace` (pvecback, pvecthermo, pvecmetric)

**Files:**
- Modify: `source/perturbations.h:343-345` (member decls) and `:452-454` (storage block — delete)
- Modify: `source/perturbations_module.cpp:2099-2104` (wiring) + `.data()` sites flagged by the compiler

- [ ] **Step 1: Replace the member declarations**

In `source/perturbations.h`, replace:
```cpp
  double* pvecback;   /**< background quantities */
  double* pvecthermo; /**< thermodynamics quantities */
  double* pvecmetric; /**< metric quantities */
```
with:
```cpp
  std::vector<double> pvecback;   /**< background quantities */
  std::vector<double> pvecthermo; /**< thermodynamics quantities */
  std::vector<double> pvecmetric; /**< metric quantities */
```

- [ ] **Step 2: Delete the now-redundant storage block**

In `source/perturbations.h`, delete the three lines:
```cpp
  std::vector<double> pvecback_storage;
  std::vector<double> pvecthermo_storage;
  std::vector<double> pvecmetric_storage;
```

- [ ] **Step 3: Drop the aliasing lines in the wiring**

In `source/perturbations_module.cpp`, replace:
```cpp
  ppw->pvecback_storage.resize(background_module_->bg_size_normal_);
  ppw->pvecthermo_storage.resize(thermodynamics_module_->th_size_);
  ppw->pvecmetric_storage.resize(ppw->mt_size);
  ppw->pvecback   = ppw->pvecback_storage.data();
  ppw->pvecthermo = ppw->pvecthermo_storage.data();
  ppw->pvecmetric = ppw->pvecmetric_storage.data();
```
with:
```cpp
  ppw->pvecback.resize(background_module_->bg_size_normal_);
  ppw->pvecthermo.resize(thermodynamics_module_->th_size_);
  ppw->pvecmetric.resize(ppw->mt_size);
```

- [ ] **Step 4: Build — fix the pointer-passing sites**

Run:
```bash
cmake --build . --target class -j 2>&1 | grep -E "error:" | head -40
```
Expected: errors at the `double* pvecback = ppw->pvecback;` style locals (around `perturbations_module.cpp:4435, 5138-5140, 5549-5551, 5800-5802, 6188-6190`) and `switch_ctx.pvecback = ppw->pvecback;` (~`:3141`). Append `.data()` to each (`ppw->pvecback` → `ppw->pvecback.data()`, etc.). `pvecback[i]` indexing is unaffected. Rebuild until clean:
```bash
cmake --build . --target class -j 2>&1 | tail -5
```
Expected: build succeeds.

- [ ] **Step 5: Byte-compare against the reference**

Run the Task 0 Step 5 compare helper.
Expected: `BYTE-IDENTICAL`.

- [ ] **Step 6: Commit**

```bash
git add source/perturbations.h source/perturbations_module.cpp
git commit -m "$(printf '%s\n' 'refactor(#312): collapse perturb_workspace pvec half-RAII twins' '' 'pvecback/pvecthermo/pvecmetric become single std::vector members;' 'raw pointers, *_storage vectors and the .data() aliasing dropped.' 'Byte-identical.' '' 'Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 3: Collapse `jacobian` (8 1D twins + dfdy, LU)

**Files:**
- Modify: `include/evolver_ndf15.h:12-51` (struct def)
- Modify: `tools/evolver_ndf15.cpp` (`initialize_jacobian` + call sites)

- [ ] **Step 1: Rewrite the `jacobian` struct definition**

In `include/evolver_ndf15.h`, replace the whole `struct jacobian { … };` (lines 12-51) with:
```cpp
struct jacobian {
  /*Stuff for normal method: */
  std::vector<double*> dfdy;   /* row pointers into dfdy_data_vec */
  std::vector<double> jacvec;  /*Stores experience gained from subsequent calls */
  std::vector<double*> LU;     /* row pointers into LU_data_vec */
  std::vector<double> LUw;
  std::vector<int> luidx;
  /*Sparse stuff:*/
  int use_sparse;
  int sparse_stuff_initialized;
  int max_nonzero; /*Maximal number of non-zero entries to be considered sparse */
  int repeated_pattern;
  int trust_sparse; /* Number of times a pattern is repeated (actually included) before we trust it. */
  int has_grouping;
  int has_pattern;
  int new_jacobian; /* True if sp_ludcmp has not been run on the current jacobian. */
  int cnzmax;
  std::vector<int> col_group;  /* Column grouping. Groups go from 0 to max_group*/
  std::vector<int> col_wi;     /* Workarray for column grouping*/
  int max_group;               /*Number of columngroups -1 */
  std::unique_ptr<sp_mat> spJ; /* Stores the matrix we want to decompose */
  std::vector<double> xjac;    /*Stores the values of the sparse jacobian. (Same pattern as spJ) */
  std::unique_ptr<sp_num> Numerical; /*Stores the LU decomposition.*/
  std::vector<int> Cp; /* Stores the column pointers of the spJ+spJ' sparsity pattern. */
  std::vector<int> Ci; /* Stores the row indices of the  spJ+spJ' sparsity pattern. */

  /* Contiguous flat backing for the 2D row-pointer arrays above */
  std::vector<double> dfdy_data_vec;
  std::vector<double> LU_data_vec;
};
```
(The raw `double**`/`double*`/`int*` members and the redundant `*_rows_vec` / 1D `*_vec` backing vectors are gone; `dfdy`/`LU` are now the row-pointer vectors themselves, still backed by the flat `*_data_vec`.)

- [ ] **Step 2: Update `initialize_jacobian`**

In `tools/evolver_ndf15.cpp`, in `initialize_jacobian`, replace the dense-setup block:
```cpp
  jac->dfdy_data_vec.assign(neq * neq + 1, 0.0);
  jac->dfdy_rows_vec.resize(neq + 1);
  jac->dfdy    = jac->dfdy_rows_vec.data();
  jac->dfdy[0] = nullptr;
  jac->dfdy[1] = jac->dfdy_data_vec.data();
  for (i = 2; i <= neq; i++)
    jac->dfdy[i] = jac->dfdy[i - 1] + neq; /* Set row pointers... */

  jac->LU_data_vec.assign(neq * neq + 1, 0.0);
  jac->LU_rows_vec.resize(neq + 1);
  jac->LU    = jac->LU_rows_vec.data();
  jac->LU[0] = nullptr;
  jac->LU[1] = jac->LU_data_vec.data();
  for (i = 2; i <= neq; i++)
    jac->LU[i] = jac->LU[i - 1] + neq; /* Set row pointers... */

  jac->LUw_vec.resize(neq + 1);
  jac->LUw = jac->LUw_vec.data();
  jac->jacvec_vec.resize(neq + 1);
  jac->jacvec = jac->jacvec_vec.data();
  jac->luidx_vec.resize(neq + 1);
  jac->luidx = jac->luidx_vec.data();
```
with:
```cpp
  jac->dfdy_data_vec.assign(neq * neq + 1, 0.0);
  jac->dfdy.resize(neq + 1);
  jac->dfdy[0] = nullptr;
  jac->dfdy[1] = jac->dfdy_data_vec.data();
  for (i = 2; i <= neq; i++)
    jac->dfdy[i] = jac->dfdy[i - 1] + neq; /* Set row pointers... */

  jac->LU_data_vec.assign(neq * neq + 1, 0.0);
  jac->LU.resize(neq + 1);
  jac->LU[0] = nullptr;
  jac->LU[1] = jac->LU_data_vec.data();
  for (i = 2; i <= neq; i++)
    jac->LU[i] = jac->LU[i - 1] + neq; /* Set row pointers... */

  jac->LUw.resize(neq + 1);
  jac->jacvec.resize(neq + 1);
  jac->luidx.resize(neq + 1);
```
and in the sparse block replace:
```cpp
    jac->xjac_vec.resize(jac->max_nonzero);
    jac->xjac = jac->xjac_vec.data();
    jac->col_group_vec.resize(neq);
    jac->col_group = jac->col_group_vec.data();
    jac->col_wi_vec.resize(neq);
    jac->col_wi = jac->col_wi_vec.data();
    jac->Cp_vec.resize(neq + 1);
    jac->Cp = jac->Cp_vec.data();
    jac->Ci_vec.resize(jac->cnzmax);
    jac->Ci = jac->Ci_vec.data();
```
with:
```cpp
    jac->xjac.resize(jac->max_nonzero);
    jac->col_group.resize(neq);
    jac->col_wi.resize(neq);
    jac->Cp.resize(neq + 1);
    jac->Ci.resize(jac->cnzmax);
```
The trailing `for (i = 1; i <= neq; i++) jac->jacvec[i] = …;` loop is unchanged (`operator[]` still works). `uninitialize_jacobian` is already pure RAII and needs **no** change.

- [ ] **Step 3: Build — fix the pointer-passing sites**

Run:
```bash
cmake --build . --target class -j 2>&1 | grep -E "error:" | head -40
```
Expected: errors only where a member is passed as a pointer, e.g.:
- `ludcmp(jac->LU, neq, jac->luidx, &luparity, jac->LUw)` (~`:987`) → `ludcmp(jac->LU.data(), neq, jac->luidx.data(), &luparity, jac->LUw.data())`
- `dFdy = jac->dfdy;` (local `double** dFdy`, ~`:1239`) → `dFdy = jac->dfdy.data();`
- any `calc_C` / sparse-path call passing `jac->xjac` / `jac->col_group` / `jac->col_wi` / `jac->Cp` / `jac->Ci` / `jac->jacvec` as a pointer → append `.data()`.

`jac->dfdy[i][j]`, `jac->LU[i][j]`, `jac->jacvec[i]` indexing produce **no** error. Append `.data()` only at flagged sites; rebuild until clean:
```bash
cmake --build . --target class -j 2>&1 | tail -5
```
Expected: build succeeds.

- [ ] **Step 4: Byte-compare against the reference**

Run the Task 0 Step 5 compare helper.
Expected: `BYTE-IDENTICAL`. (This task touches the dense + sparse linear-algebra path, so the byte-compare is the key safety check.)

- [ ] **Step 5: Commit**

```bash
git add include/evolver_ndf15.h tools/evolver_ndf15.cpp
git commit -m "$(printf '%s\n' 'refactor(#312): collapse jacobian half-RAII twins' '' '8 1D members become single std::vector members; dfdy/LU keep [i][j]' 'by renaming the row-pointer vectors and dropping the redundant' 'double** twins (flat *_data_vec backing retained). Byte-identical.' '' 'Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 4: Collapse `numjac_workspace` (12 1D twins + ydel_Fdel)

**Files:**
- Modify: `include/evolver_ndf15.h:53-86` (struct def)
- Modify: `tools/evolver_ndf15.cpp` (`initialize_numjac_workspace` + call sites)

- [ ] **Step 1: Rewrite the `numjac_workspace` struct definition**

In `include/evolver_ndf15.h`, replace the whole `struct numjac_workspace { … };` (lines 53-86) with:
```cpp
struct numjac_workspace {
  std::vector<double> yscale;
  std::vector<double> del;
  std::vector<double> Difmax;
  std::vector<double> absFdelRm;
  std::vector<double> absFvalue;
  std::vector<double> absFvalueRm;
  std::vector<double> Fscale;
  std::vector<double> ffdel;
  std::vector<double> yydel;
  std::vector<double> tmp;

  std::vector<double*> ydel_Fdel;  /* row pointers into ydel_Fdel_data_vec */

  std::vector<int> logj;
  std::vector<int> Rowmax;

  /* Contiguous flat backing for ydel_Fdel */
  std::vector<double> ydel_Fdel_data_vec;
};
```

- [ ] **Step 2: Update `initialize_numjac_workspace`**

In `tools/evolver_ndf15.cpp`, replace the body's allocation block:
```cpp
  nj_ws->yscale_vec.resize(neqp);
  nj_ws->yscale = nj_ws->yscale_vec.data();
  nj_ws->del_vec.resize(neqp);
  nj_ws->del = nj_ws->del_vec.data();
  nj_ws->Difmax_vec.resize(neqp);
  nj_ws->Difmax = nj_ws->Difmax_vec.data();
  nj_ws->absFdelRm_vec.resize(neqp);
  nj_ws->absFdelRm = nj_ws->absFdelRm_vec.data();
  nj_ws->absFvalue_vec.resize(neqp);
  nj_ws->absFvalue = nj_ws->absFvalue_vec.data();
  nj_ws->absFvalueRm_vec.resize(neqp);
  nj_ws->absFvalueRm = nj_ws->absFvalueRm_vec.data();
  nj_ws->Fscale_vec.resize(neqp);
  nj_ws->Fscale = nj_ws->Fscale_vec.data();
  nj_ws->ffdel_vec.resize(neqp);
  nj_ws->ffdel = nj_ws->ffdel_vec.data();
  nj_ws->yydel_vec.resize(neqp);
  nj_ws->yydel = nj_ws->yydel_vec.data();
  nj_ws->tmp_vec.resize(neqp);
  nj_ws->tmp = nj_ws->tmp_vec.data();

  nj_ws->ydel_Fdel_data_vec.assign(neq * neq + 1, 0.0);
  nj_ws->ydel_Fdel_rows_vec.resize(neq + 1);
  nj_ws->ydel_Fdel    = nj_ws->ydel_Fdel_rows_vec.data();
  nj_ws->ydel_Fdel[0] = nullptr;
  nj_ws->ydel_Fdel[1] = nj_ws->ydel_Fdel_data_vec.data();
  for (int i = 2; i <= neq; i++)
    nj_ws->ydel_Fdel[i] = nj_ws->ydel_Fdel[i - 1] + neq; /* Set row pointers... */

  nj_ws->logj_vec.resize(neqp);
  nj_ws->logj = nj_ws->logj_vec.data();
  nj_ws->Rowmax_vec.resize(neqp);
  nj_ws->Rowmax = nj_ws->Rowmax_vec.data();
```
with:
```cpp
  nj_ws->yscale.resize(neqp);
  nj_ws->del.resize(neqp);
  nj_ws->Difmax.resize(neqp);
  nj_ws->absFdelRm.resize(neqp);
  nj_ws->absFvalue.resize(neqp);
  nj_ws->absFvalueRm.resize(neqp);
  nj_ws->Fscale.resize(neqp);
  nj_ws->ffdel.resize(neqp);
  nj_ws->yydel.resize(neqp);
  nj_ws->tmp.resize(neqp);

  nj_ws->ydel_Fdel_data_vec.assign(neq * neq + 1, 0.0);
  nj_ws->ydel_Fdel.resize(neq + 1);
  nj_ws->ydel_Fdel[0] = nullptr;
  nj_ws->ydel_Fdel[1] = nj_ws->ydel_Fdel_data_vec.data();
  for (int i = 2; i <= neq; i++)
    nj_ws->ydel_Fdel[i] = nj_ws->ydel_Fdel[i - 1] + neq; /* Set row pointers... */

  nj_ws->logj.resize(neqp);
  nj_ws->Rowmax.resize(neqp);
```
`uninitialize_numjac_workspace` is already pure RAII and needs **no** change.

- [ ] **Step 3: Build — fix the pointer-passing sites**

Run:
```bash
cmake --build . --target class -j 2>&1 | grep -E "error:" | head -40
```
Expected: errors only where a `numjac_workspace` member is passed as a pointer argument or assigned to a `double*`/`int*`/`double**` local (e.g. inside `numjac`, where `nj_ws->yydel`, `nj_ws->ffdel`, `nj_ws->del`, `nj_ws->ydel_Fdel`, etc. are handed to helpers). Append `.data()` at each. Indexing (`nj_ws->ydel_Fdel[i][j]`, `nj_ws->del[i]`, `nj_ws->Rowmax[j]`) produces **no** error. Rebuild until clean:
```bash
cmake --build . --target class -j 2>&1 | tail -5
```
Expected: build succeeds.

- [ ] **Step 4: Byte-compare against the reference**

Run the Task 0 Step 5 compare helper.
Expected: `BYTE-IDENTICAL`. (Numerical-Jacobian path — byte-compare is the key safety check.)

- [ ] **Step 5: Confirm no twins remain**

Run:
```bash
grep -nE "_storage|_rows_vec|_vec\b" include/evolver_ndf15.h
grep -nE "pvec.*_storage|y_storage|dy_storage|used_in_sources_storage" source/perturbations.h
```
Expected: the only matches in `evolver_ndf15.h` are the two intended flat backings `dfdy_data_vec`, `LU_data_vec`, `ydel_Fdel_data_vec`; **no** matches in `perturbations.h`. (If any other `*_vec`/`*_storage` twin remains, collapse it the same way before committing.)

- [ ] **Step 6: Commit**

```bash
git add include/evolver_ndf15.h tools/evolver_ndf15.cpp
git commit -m "$(printf '%s\n' 'refactor(#312): collapse numjac_workspace half-RAII twins' '' '12 1D members become single std::vector members; ydel_Fdel keeps' '[i][j] via the renamed row-pointer vector (flat backing retained).' 'Byte-identical. Closes #312.' '' 'Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>')"
```

---

### Task 5: Full-suite verification + PR

**Files:** none (verification + PR).

- [ ] **Step 1: Build the Python module**

Run:
```bash
pip install . --no-build-isolation 2>&1 | tail -5
```
Expected: `classy` builds and installs cleanly.

- [ ] **Step 2: Run the full test suite at TEST_LEVEL=2**

Run:
```bash
cd python && TEST_LEVEL=2 python -m pytest test_class.py -q 2>&1 | tail -20; cd ..
```
Expected: all tests pass against the existing `classyref` (no regressions; output unchanged, so no golden/classyref regeneration).

- [ ] **Step 3: Final byte-compare**

Run the Task 0 Step 5 compare helper one more time after the full build.
Expected: `BYTE-IDENTICAL`.

- [ ] **Step 4: Push and open the PR**

```bash
git push -u origin issue-312-collapse-half-raii
gh pr create --title "v4 #312: collapse the raw-pointer + backing-vector half-RAII pattern" --body "$(printf '%s\n' 'Closes #312.' '' 'Makes the owning std::vector the single member across perturb_vector,' 'perturb_workspace, and the ndf15 jacobian / numjac_workspace; deletes the' 'raw-pointer twins and their .data() aliasing wiring. 2D row-pointer arrays' '(dfdy/LU/ydel_Fdel) keep [i][j] by renaming the row-pointer vectors and' 'dropping the redundant double** twins (flat *_data_vec backing retained).' '' 'Pure mechanical/storage refactor — byte-identical output, verified by' 'bit-compare against pre-refactor HEAD and full TEST_LEVEL=2 suite.' '' 'Second-to-last v4-prep item; only #313 remains.' '' '🤖 Generated with [Claude Code](https://claude.com/claude-code)')"
```
Expected: PR created.

---

## Self-Review

**Spec coverage:** perturb_vector (Task 1), perturb_workspace (Task 2), jacobian incl. dfdy/LU (Task 3), numjac_workspace incl. ydel_Fdel (Task 4) — all four struct families in the spec's scope table are covered. The "keep `[i][j]`, drop only the redundant raw member" decision is realized in Tasks 3 & 4 (rename `*_rows_vec` → bare name, retain `*_data_vec`). Out-of-scope items (2D flattening, `unique_ptr` members, #313) are left untouched. Verification model (byte-identical vs HEAD, no golden regeneration, full TEST_LEVEL=2) is in Task 0 + Task 5.

**Placeholder scan:** none — every edit shows the exact before/after; the only "find the rest yourself" steps are the `.data()` fixes, which are intentionally compiler-driven (the build error list is the exhaustive worklist) and bounded by the named representative sites.

**Type consistency:** `dfdy`/`LU`/`ydel_Fdel` are `std::vector<double*>` (so `m[i][j]` and `m.data()` as `double**` both hold); 1D members are `std::vector<double>`/`std::vector<int>` (so `m[i]` and `m.data()` hold); the retained flat backings keep the `_data_vec` suffix. `initialize_*` writes through the renamed members; `uninitialize_*` unchanged (RAII). Consistent across tasks.

---

## ADDENDUM 2026-06-19 — full sweep (Tasks 5–8 added; Task 5→9)

Per the owner decision (see spec addendum), #312 became a full sweep of the half-RAII pattern. Task 2 expands to also cover `perturb_workspace::approx` and `s_l`. Four new file-tasks are inserted; the verify+PR task is renumbered to Task 9. The expanded byte-compare oracle is captured at HEAD in `/tmp/ref312/` (REF_*.dat) — see the spec's verification matrix.

**Generalized recipe (identical for every file/struct below):**
1. In the header, for each twin: rename the backing vector (`_storage`, or `_rows` for 2D) to the bare raw member name; delete the raw pointer member (`double*`/`int*`/`double**`). For 2D in nonlinear, retain the jagged `*_storage` (`vector<vector<double>>`) and rename `*_rows` (`vector<double*>`) to the bare name. Keep doc comments.
2. In the init/alloc function, delete every `X = X_*.data()` (and `rows[i] = storage[i].data()` collapses to writing into the renamed rows member) aliasing line; size the renamed vectors directly.
3. Build (`cmake --build . --target class -j 2>&1 | grep -E "error:"`) — the error list is the exhaustive `.data()` worklist. Append `.data()` ONLY at flagged pointer-passing sites; never touch `[i]`/`[i][j]` indexing. Rebuild until clean.
4. Byte-compare the relevant reference run(s) — must be `BYTE-IDENTICAL`.
5. Commit.

### Task 2 (expanded): `perturb_workspace` — pvecback/pvecthermo/pvecmetric **+ approx + s_l**
Files: `source/perturbations.h` (struct + delete the `*_storage` block incl. `approx_storage`/`s_l_storage`), `source/perturbations_module.cpp` (wiring + `.data()` sites). Verify: ref_default (`BYTE-IDENTICAL` on REF_cl/cl_lensed/pk).

### Task 5 (NEW): `include/dei_rkck.h` — 11 1D twins
`yscal`, `y`, `dydx`, `yerr`, `ytempo`, `ak2`..`ak6`, `ytemp` (`double*` + `_storage`). Init in `tools/dei_rkck.cpp` (`initialize_generic_integrator` / the Cash-Karp workspace alloc). Call sites in `tools/dei_rkck.cpp`. Verify: ref_default (exercised via thermo RECFAST + primordial).

### Task 6 (NEW): `source/transfer.h` — 14 1D twins
`interpolated_sources`, `sources`, `tau0_minus_tau`, `w_trapz`, `chi`, `cscKgen`, `cotKgen`, `chireverse`, `rescale_function`, `radial_function`, `chi_full_reverse`, **`Phi`** (storage `phi_storage`), **`dPhi`** (`dphi_storage`), **`d2Phi`** (`d2phi_storage`). GOTCHA: rename the lowercase storage to the CamelCase raw name for the three Phi members. Init/alloc + call sites in `source/transfer_module.cpp`. Verify: ref_default.

### Task 7 (NEW): `source/thermodynamics.h` — 6 1D input-config twins
`binned_reio_z`, `binned_reio_xe`, `many_tanh_z`, `many_tanh_xe`, `reio_inter_z`, `reio_inter_xe` (`double* X = nullptr` + `X_storage`). These are parsed from input; collapsing drops the `= nullptr` default (empty vector). Find readers/writers in `source/input_module.cpp` and `source/thermodynamics_module.cpp`. Verify: ref_reiobins + ref_reiomanytanh + ref_reiointer (each `BYTE-IDENTICAL` on its REF_*_cl/pk).

### Task 8 (NEW): `source/nonlinear.h` — 6 1D + 4 jagged-2D twins
1D: `rtab`, `stab`, `ddstab`, `growtable`, `ztable`, `tautable`. 2D (jagged): `sigma_8`, `sigma_disp`, `sigma_disp_100`, `sigma_prime` — each `double**` + `vector<vector<double>> *_storage` + `vector<double*> *_rows`. Treatment: drop raw `double**`, rename `*_rows`→bare member (`vector<double*>`), RETAIN the jagged `*_storage`; init populates rows as `member[i] = storage[i].data()`. `[i][j]` preserved; `double**`-passing sites get `.data()`. Init/alloc + call sites in `source/nonlinear_module.cpp`. Verify: ref_hmcode (`BYTE-IDENTICAL` on REF_hmcode_cl/cl_lensed/pk/pk_nl).

### Task 9 (was Task 5): full-suite verify + PR
`pip install . --no-build-isolation`; `TEST_LEVEL=2` pytest green; re-run ALL reference scenarios (default + hmcode + 3 reio) and confirm `BYTE-IDENTICAL`; push; open PR closing #312.
