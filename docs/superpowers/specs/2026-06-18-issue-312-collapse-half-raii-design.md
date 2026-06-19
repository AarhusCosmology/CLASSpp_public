# Issue #312 — collapse the raw-pointer + backing-vector half-RAII pattern

**Date:** 2026-06-18
**Issue:** [#312](https://github.com/AarhusCosmology/CLASSpp/issues/312) — collapse the raw-pointer + backing-vector half-RAII pattern
**Tech Stack:** C++17 (`include/evolver_ndf15.h`, `tools/evolver_ndf15.cpp`, `source/perturbations.h`, `source/perturbations_module.cpp`); build via CMake (`CMakeLists.txt`); Python `classy` Cython wrapper; reference comparison via `classyref` + the `TEST_LEVEL`/`COMPARE_OUTPUT_REF` env knobs in `python/test_class.py`.

## Goal

Several structs carry a raw-pointer API *and* the `std::vector` that owns the storage, with the pointer aliasing into the vector (`raw = vec.data()` set in an init function). It is leak-safe but transitional: it doubles the member count and leaves a re-allocation hazard (any resize silently invalidates the raw twin).

This issue makes the `std::vector` the single owning member and deletes the raw twin. Call sites that index (`m[i]`, `m[i][j]`) keep working through `vector::operator[]`; the few sites that pass a member as a `double*`/`int*`/`double**` argument use `.data()`. **Output is byte-identical** — this is a pure mechanical/storage refactor with no change to floating-point arithmetic or evaluation order.

This is the second-to-last item of the v4 design-review programme (#307–#314); only #313 remains after it.

## Scope — the affected members

The half-RAII twins live in exactly two places (both struct families are evolver-local: declared in `include/evolver_ndf15.h`, used only in `tools/evolver_ndf15.cpp`, plus the perturbation structs in `source/perturbations.{h,cpp}`).

### A. 1D twins (~26 members) — `raw + backing vector` → single vector

| Struct | Members |
|---|---|
| `perturb_vector` (`perturbations.h:298–306`) | `y`, `dy`, `used_in_sources` |
| `perturb_workspace` (`perturbations.h:343–345`, storage `:452–454`) | `pvecback`, `pvecthermo`, `pvecmetric` |
| `jacobian` (`evolver_ndf15.h:12–51`) | `jacvec`, `LUw`, `luidx`, `xjac`, `col_group`, `col_wi`, `Cp`, `Ci` |
| `numjac_workspace` (`evolver_ndf15.h:53–86`) | `yscale`, `del`, `Difmax`, `absFdelRm`, `absFvalue`, `absFvalueRm`, `Fscale`, `ffdel`, `yydel`, `tmp`, `logj`, `Rowmax` |

**Treatment:** rename the backing vector to the bare member name (drop the `_storage` / `_vec` suffix), delete the raw pointer, and delete the `raw = vec.data()` aliasing line from the init function. `m[i]` indexing is unchanged via `operator[]`; sites passing the member as a pointer argument get `.data()`.

### B. 2D row-pointer arrays (3 members) — keep `[i][j]`, drop only the redundant raw twin

| Struct | Member | Backing |
|---|---|---|
| `jacobian` | `dfdy` | `dfdy_data_vec` (flat) + `dfdy_rows_vec` (`vector<double*>`) |
| `jacobian` | `LU` | `LU_data_vec` + `LU_rows_vec` |
| `numjac_workspace` | `ydel_Fdel` | `ydel_Fdel_data_vec` + `ydel_Fdel_rows_vec` |

**Treatment (per the #312 decision — "keep `[i][j]`, drop only the redundant raw member"):** keep the contiguous `*_data_vec` (flat storage) and rename the `*_rows_vec` (`std::vector<double*>`) to the bare member name; delete the raw `double**` twin. `m[i][j]` then resolves as `vector<double*>::operator[](i)` → `double*`, then `[j]` — **byte-identical, same memory, no edit** at the 19/22 indexing sites. The init still populates the row pointers (`rows[i] = data.data() + i*neq`). Only the ~3 sites that pass the whole array as a `double**` (e.g. `ludcmp(jac->LU, …)` at `evolver_ndf15.cpp:987`, the local `dFdy = jac->dfdy` at `:1239`) get `.data()`.

We deliberately do **not** flatten to `data[i*stride+j]` — it rewrites every indexing site for the same byte output and adds risk in numerically-sensitive code.

## Init-function changes

- `initialize_jacobian` / `initialize_numjac_workspace` (`tools/evolver_ndf15.cpp`): drop every `raw = vec.data()` aliasing line for 1D members; for the 2D arrays, keep the `*_data_vec.assign(...)` sizing and the row-pointer population loop, now writing into the renamed rows vector.
- `perturb_vector` setup (`perturbations_module.cpp:3047–3053`): the three `ppv->… = ppv->…_storage.data()` aliasing lines are deleted; the `.assign()/.resize()` sizing of the (renamed) vectors stays.
- `uninitialize_*`: any explicit raw-pointer nulling becomes unnecessary and is removed.

## Out of scope

- Flattening the 2D arrays to stride arithmetic (rejected above).
- The `std::unique_ptr<sp_mat>` / `std::unique_ptr<sp_num>` members of `jacobian` — already proper RAII, untouched.
- Everything in #313 (common.h split, typed thrower, enum class) — that is the next and final v4-prep item.

## Verification

The bar is **byte-identical output vs current master HEAD (`f0939a2b`)**. There is no behavior change, so nothing in `classyref` or the goldens is regenerated; a single differing bit is a bug in the refactor, not a new baseline.

1. From clean master HEAD, build and capture a representative `TEST_LEVEL=2` reference run.
2. Apply the refactor on the branch; rebuild.
3. Bit-compare the same run's output against the captured reference — must match exactly.
4. Run the full test suite (`python/test_class.py`, `TEST_LEVEL=2`) — must stay green with the existing `classyref`.

## Files touched

- `include/evolver_ndf15.h` — struct definitions (delete raw twins, rename backing vectors).
- `tools/evolver_ndf15.cpp` — init functions + ~158 call sites (mostly unchanged indexing; `.data()` at the handful of pointer-passing sites).
- `source/perturbations.h` — `perturb_vector` / `perturb_workspace` struct members.
- `source/perturbations_module.cpp` — drop aliasing lines at `:3047`; `.data()` at the `pv->y`-as-`double*` sites.

## Branch / PR

Branch `issue-312-collapse-half-raii` off master; one PR closing #312.

---

## ADDENDUM 2026-06-19 — scope expanded to a FULL sweep (user decision)

During execution, a codebase sweep (`grep -rE "std::vector<…>\s+\w+_storage"`) revealed the half-RAII pattern is **far more widespread** than issue #312's bullet list (which named only `perturb_vector`/`jacobian`/`numjac_workspace`/`perturb_workspace` as *examples*). The owner chose to make this PR a **full sweep** so the issue's stated Value — "makes the RAII story consistent before the v4 presentation" — actually lands. Full inventory (~45 twins, 6 files):

| File | Struct | Twins | Notes |
|---|---|---|---|
| `source/perturbations.h` | `perturb_vector` | `y`, `dy`, `used_in_sources` | ✅ done (commit 4fd4e91b) |
| `source/perturbations.h` | `perturb_workspace` | `pvecback`, `pvecthermo`, `pvecmetric`, **`approx`**, **`s_l`** | originally under-scoped; `approx`/`s_l` added |
| `include/evolver_ndf15.h` | `jacobian` | 8 1D + `dfdy`/`LU` (2D) | 2D: keep `[i][j]` (see §B) |
| `include/evolver_ndf15.h` | `numjac_workspace` | 12 1D + `ydel_Fdel` (2D) | 2D: keep `[i][j]` |
| `include/dei_rkck.h` | `dei_rkck` | 11 1D (`yscal`,`y`,`dydx`,`yerr`,`ytempo`,`ak2..ak6`,`ytemp`) | the Cash-Karp `generic_integrator` workspace (used by thermo RECFAST + primordial **by default**) |
| `source/transfer.h` | `transfer_workspace` | 14 1D | **name mismatch:** raw `Phi`/`dPhi`/`d2Phi` ↔ storage `phi_storage`/`dphi_storage`/`d2phi_storage`; rename storage to the CamelCase raw name |
| `source/thermodynamics.h` | thermo input struct | 6 1D (`binned_reio_z/xe`, `many_tanh_z/xe`, `reio_inter_z/xe`) | input-parsed config arrays (`double* X = nullptr` + `X_storage`); reachable only with exotic `reio_parametrization` |
| `source/nonlinear.h` | `nonlinear_workspace` | 6 1D + **4 2D** (`sigma_8`,`sigma_disp`,`sigma_disp_100`,`sigma_prime`) | 2D is **jagged**: `double**` + `vector<vector<double>> *_storage` + `vector<double*> *_rows`. Keep `[i][j]`: rename `*_rows`→bare member, retain `*_storage` jagged backing, drop raw `double**` |

The treatment is identical to §A/§B for every file: rename the backing vector to the bare member name (drop `_storage`/`_rows`), delete the raw twin and its `raw = vec.data()` aliasing in the init function; `[i]`/`[i][j]` indexing is preserved; pointer-passing sites get `.data()`. For the jagged 2D in nonlinear, the `vector<vector<double>>` storage is retained (jagged, not flat) and the renamed `vector<double*>` rows member makes `[i][j]` work.

### Verification matrix (expanded oracle, all captured byte-deterministic at HEAD)

| Scenario (ini) | Covers |
|---|---|
| `/tmp/ref312.ini` (default LCDM, lensed Cls + mPk) | perturb_vector, perturb_workspace, jacobian, numjac_workspace, transfer, **dei_rkck** (thermo RECFAST + primordial), thermo (default reio) |
| `/tmp/ref312_hmcode.ini` (`non linear = hmcode`) | nonlinear.h (1D `rtab`/`stab`/`ddstab`/`growtable` + 2D `sigma_*`) |
| `/tmp/ref312_reiobins.ini` (`reio_bins_tanh`) | thermo `binned_reio_z/xe` |
| `/tmp/ref312_reiomanytanh.ini` (`reio_many_tanh`) | thermo `many_tanh_z/xe` |
| `/tmp/ref312_reiointer.ini` (`reio_inter`) | thermo `reio_inter_z/xe` |

Each task byte-compares the run(s) that exercise its file against the captured `REF_*.dat`; the final task runs all scenarios + the full `TEST_LEVEL=2` python suite. (The opt-in `evolver = rk` *perturbation* path segfaults on master — a pre-existing, unrelated bug; `dei_rkck` is still fully covered via the default thermo/primordial usage.)
