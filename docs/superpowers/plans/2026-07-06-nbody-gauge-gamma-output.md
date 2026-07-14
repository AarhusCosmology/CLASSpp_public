# N-body gauge `k2gamma_Nb` output — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the N-body-gauge `γ` transfer output by supplying the missing `−H_T_Nb''` term via a per-`k` cubic-spline derivative of the stored `H_T_Nb'(τ)` source, and remove the `class_stop` that currently aborts all N-body-gauge output.

**Architecture:** Two source-only edits to `source/perturbations_module.cpp`. (1) Remove the `class_stop` in `perturb_sources` so N-body-gauge output runs (density/velocity gauge conversion needs only `H_T_Nb'`, already available). (2) In `perturb_solve`, after the evolver finishes each `k`, spline the stored `H_T_Nb'(τ)` time series over `tau_sampling_`, read the analytic derivative `H_T_Nb''` off that spline, and subtract it from the stored (partial) `k2gamma_Nb`. No `dτ` step to tune; reuses `array_spline_table_line_to_line` + `array_derive_spline_table_line_to_line` (same idiom as `dκ/dτ` in thermodynamics).

**Tech Stack:** C++17 (CLASSpp), CMake / scikit-build-core, `classy` Python wrapper, NumPy for verification scripts.

## Global Constraints

- C++-only; no C-compatibility guards (`extern "C"`, `#ifdef __cplusplus`). No new headers needed — `array_spline_table_line_to_line`, `array_derive_spline_table_line_to_line`, and `_SPLINE_EST_DERIV_` are declared in `include/arrays.h`, already included and used in this file (`source/perturbations_module.cpp:695`, `:703`).
- **No committed golden or regression test** (user decision). All verification scripts are throwaway, written under the scratchpad dir, never `git add`-ed.
- Never `git add -A` in this repo (in-source build artifacts get swept in). Stage the single explicit source file only.
- Both commits land on branch `nbody-gauge-gamma-output` (already checked out; the design spec is committed there as `fbdf76ac`).
- Fast rebuild of `classy` after a C++ edit: `pip install --no-build-isolation .` from the repo root.
- Verification tolerance is loose (≈10% on the smooth low-`k` plateau): the reference reconstruction uses a 2-point finite difference while the code uses a spline, so they agree only where `H_T_Nb''` is smooth. Large-`k` oscillatory disagreement is expected and out of scope for the automated check (validated visually instead).
- Scratchpad dir (verification scripts + baseline): `/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/e892d42f-6641-4ba9-afcc-04d82d0fab80/scratchpad`

---

### Task 1: Remove the `class_stop` — unblock N-body-gauge output

Removing the stop lets `compute()` succeed with N-body-gauge flags. `k2gamma_Nb` is emitted but still **partial** (missing `−H_T_Nb''`) — Task 2 completes it. This is an independently reviewable step: it unblocks the density/velocity gauge conversion (which never needed `γ`) and can be accepted even if a reviewer wants changes to the derivative approach in Task 2.

**Files:**
- Modify: `source/perturbations_module.cpp:4947-4959` (remove `class_stop`, update comment)
- Verify (throwaway, NOT committed): `<scratchpad>/check_unblock.py`, `<scratchpad>/save_baseline.py`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: for Task 2 — the stored `k2gamma_Nb` column now holds the **partial** value `−ℋ·H_T_Nb' + (9/2)a²(ρ+p)σ`; and a saved reference baseline `<scratchpad>/gamma_baseline.npz` containing, per test redshift, the arrays `kh` (k in h/Mpc), `partial` (partial `k2gamma_Nb`), `phi`, and `gamma_ref` (the notebook's full-γ reconstruction `partial − H_T_Nb''_fd`).

- [ ] **Step 1: Write the throwaway smoke test**

Create `<scratchpad>/check_unblock.py`:

```python
from classy import Class
import numpy as np

default = {'output': 'dTk, vTk, mPk', 'P_k_max_1/Mpc': 10, 'z_pk': '0, 1000',
           'extra metric transfer functions': 'yes', 'k_per_decade_for_pk': 100,
           'radiation_streaming_approximation': 3, 'k_output_values': 1e-6,
           'gauge': 'Newtonian'}

cosmo = Class()
cosmo.set(default)
cosmo.compute()                      # currently aborts at the class_stop
tk = cosmo.get_transfer(z=49)
assert 'k2gamma_Nb' in tk, "k2gamma_Nb column missing"
assert np.all(np.isfinite(tk['k2gamma_Nb'])), "k2gamma_Nb has non-finite entries"
print("UNBLOCK OK: k2gamma_Nb present and finite")
```

- [ ] **Step 2: Build current code and run the test to verify it FAILS**

Run:
```bash
pip install --no-build-isolation . && \
python "<scratchpad>/check_unblock.py"
```
Expected: FAIL — a `CosmoComputationError`/`CosmoSevereError` whose message contains
"We need to compute the derivative of H_T_Nb_prime numerically". (This is the `class_stop` firing.)

- [ ] **Step 3: Remove the `class_stop` and update the comment**

In `source/perturbations_module.cpp`, replace this block (`:4947-4959`):

```cpp
      /** gamma in Nbody gauge, see Eq. A.2 in 1811.00904. */
      if (has_source_k2gamma_Nb_) {
        _set_source_(index_tp_k2gamma_Nb_) = -a_prime_over_a * H_T_Nb_prime +
                                             9. / 2. * a2 * ppw->rho_plus_p_shear;
      }
    }

    if (has_source_k2gamma_Nb_) {
      class_stop(
          "We need to compute the derivative of H_T_Nb_prime numerically. Written by T. "
          "Tram but not yet propagated here. See devel branch prior to merging with hmcode "
          "branch");
    }
```

with (drop the `class_stop`; keep the partial computation; document where the rest comes from):

```cpp
      /** gamma in Nbody gauge, see Eq. A.2 in 1811.00904. This stores only the
          partial value; the remaining -H_T_Nb'' term (conformal-time derivative
          of H_T_Nb') is added as a post-processing spline derivative in
          perturb_solve, because it cannot be formed analytically here (it needs
          the total shear derivative). */
      if (has_source_k2gamma_Nb_) {
        _set_source_(index_tp_k2gamma_Nb_) = -a_prime_over_a * H_T_Nb_prime +
                                             9. / 2. * a2 * ppw->rho_plus_p_shear;
      }
    }
```

- [ ] **Step 4: Rebuild and run the smoke test to verify it PASSES**

Run:
```bash
pip install --no-build-isolation . && \
python "<scratchpad>/check_unblock.py"
```
Expected: PASS — prints `UNBLOCK OK: k2gamma_Nb present and finite`.

- [ ] **Step 5: Save the reference baseline for Task 2**

Create `<scratchpad>/save_baseline.py` (captures the partial column and the notebook-method full-γ reconstruction, from THIS partial build):

```python
from classy import Class
import numpy as np

default = {'output': 'dTk, vTk, mPk', 'P_k_max_1/Mpc': 10, 'z_pk': '0, 1000',
           'extra metric transfer functions': 'yes', 'k_per_decade_for_pk': 100,
           'radiation_streaming_approximation': 3, 'k_output_values': 1e-6,
           'gauge': 'Newtonian'}

cosmo = Class()
cosmo.set(default)
cosmo.compute()

def H_T_Nb_prime_prime(z_val):
    # notebook get_H_T_Nb_prime_prime: 2-point z-derivative of H_T_Nb_prime -> d/dtau
    dz = 0.1
    zl, zr = z_val + dz, max(0., z_val - dz)
    delta_z = zr - zl
    HTpl = cosmo.get_transfer(z=zl)['H_T_Nb_prime']
    HTpr = cosmo.get_transfer(z=zr)['H_T_Nb_prime']
    return -(HTpr - HTpl) / delta_z * cosmo.Hubble(z_val)

out = {}
for z in (49, 100):
    tk = cosmo.get_transfer(z=z)
    out[f'kh_{z}'] = tk['k (h/Mpc)']
    out[f'partial_{z}'] = tk['k2gamma_Nb']
    out[f'phi_{z}'] = tk['phi']
    out[f'gamma_ref_{z}'] = tk['k2gamma_Nb'] - H_T_Nb_prime_prime(z)  # full gamma, notebook method

np.savez("<scratchpad>/gamma_baseline.npz", **out)
print("baseline saved")
```

Run:
```bash
python "<scratchpad>/save_baseline.py"
```
Expected: prints `baseline saved`; file `<scratchpad>/gamma_baseline.npz` exists.

- [ ] **Step 6: Commit (source only)**

```bash
git add source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
Unblock N-body gauge output: remove k2gamma_Nb class_stop

Requesting 'extra metric transfer functions' or 'Nbody gauge transfer
functions' previously aborted compute() at a class_stop that guarded the
not-yet-computed H_T_Nb'' term. Remove it; the stored k2gamma_Nb is still
the partial value (-H H_T' + shear). The -H_T'' term is added as a spline
post-process in perturb_solve in the following commit. This also unblocks the
N-body-gauge density/velocity transfers, which only ever needed H_T'.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Uzbwf8RtwEDNs9vJNE4dDM
EOF
)"
```

---

### Task 2: Add the spline `H_T_Nb''` post-process — complete `γ`

**Files:**
- Modify: `source/perturbations_module.cpp` — insert a block in `perturb_solve` after the evolver interval loop closes (`:2297`) and before the zero-fill tail comment (`:2304`).
- Verify (throwaway, NOT committed): `<scratchpad>/check_gamma.py`; visual check via `ComputeGamma.ipynb`.

**Interfaces:**
- Consumes from Task 1: stored `k2gamma_Nb` holds the partial value; `<scratchpad>/gamma_baseline.npz` (`kh_{z}`, `gamma_ref_{z}` for `z ∈ {49, 100}`).
- Produces: stored `k2gamma_Nb` now equals the full `k²γ = −ℋ·H_T_Nb' + (9/2)a²(ρ+p)σ − H_T_Nb''`; every mTk output path (`get_transfer`, `*_tk.dat`) returns the full `γ`.

- [ ] **Step 1: Write the throwaway consistency test**

Create `<scratchpad>/check_gamma.py` (asserts the in-code column matches the notebook full-γ reconstruction on the smooth low-`k` plateau):

```python
from classy import Class
import numpy as np

default = {'output': 'dTk, vTk, mPk', 'P_k_max_1/Mpc': 10, 'z_pk': '0, 1000',
           'extra metric transfer functions': 'yes', 'k_per_decade_for_pk': 100,
           'radiation_streaming_approximation': 3, 'k_output_values': 1e-6,
           'gauge': 'Newtonian'}

ref = np.load("<scratchpad>/gamma_baseline.npz")

cosmo = Class()
cosmo.set(default)
cosmo.compute()

ok = True
for z in (49, 100):
    kh = cosmo.get_transfer(z=z)['k (h/Mpc)']
    code = cosmo.get_transfer(z=z)['k2gamma_Nb']
    kh_ref = ref[f'kh_{z}']
    gamma_ref = ref[f'gamma_ref_{z}']
    assert np.allclose(kh, kh_ref), f"k grid mismatch at z={z}"

    mask = kh < 1e-3                      # smooth plateau: spline and FD agree here
    rel = np.abs(code[mask] - gamma_ref[mask]) / (np.abs(gamma_ref[mask]) + 1e-300)
    worst = float(np.max(rel))
    print(f"z={z}: max rel diff on k<1e-3 = {worst:.3e}")
    ok = ok and (worst < 0.1)            # 10% (FD-vs-spline on the plateau)

assert ok, "in-code k2gamma_Nb does not match notebook full-gamma reconstruction on the plateau"
print("GAMMA OK: full k2gamma_Nb reproduces the notebook reconstruction on the low-k plateau")
```

- [ ] **Step 2: Run the test on the Task-1 build to verify it FAILS**

Run (no rebuild — the current build still has the partial column from Task 1):
```bash
python "<scratchpad>/check_gamma.py"
```
Expected: FAIL — the stored column is the partial value, so it differs from `gamma_ref` by `H_T_Nb''` on the plateau (the `z≈50` large-scale cancellation makes this term non-negligible). The assertion trips with a large `max rel diff`.

- [ ] **Step 3: Insert the spline post-processing block in `perturb_solve`**

In `source/perturbations_module.cpp`, immediately after the evolver interval loop closes at `:2297` (the line `  }` following `perhaps_print_variables);`) and before the comment `/** - fill the source terms array with zeros ... */` at `:2304`, insert:

```cpp
  /** - complete the N-body-gauge gamma source. The stored k2gamma_Nb still
      lacks the -H_T_Nb'' term (the conformal-time derivative of H_T_Nb'),
      which cannot be formed analytically in the RHS because it needs the
      total shear derivative. Compute it here, output-only, as the analytic
      derivative of a cubic spline through the H_T_Nb'(tau) time series for
      this k over the integrated region [0, tau_actual_size), and subtract it
      from the partial k2gamma_Nb. */
  if (has_source_k2gamma_Nb_ && index_md == index_md_scalars_ && tau_actual_size >= 2) {
    const int n = tau_actual_size;
    const int stride = k_size_[index_md];
    const double* HTp =
        sources_[index_md][index_ic * tp_size_[index_md] + index_tp_H_T_Nb_prime_].data();
    double* k2g =
        sources_[index_md][index_ic * tp_size_[index_md] + index_tp_k2gamma_Nb_].data();

    /* scratch table with columns (y=0, ddy=1, dy=2) */
    std::vector<double> work(static_cast<size_t>(n) * 3);
    for (int i = 0; i < n; i++)
      work[i * 3 + 0] = HTp[i * stride + index_k];

    array_spline_table_line_to_line(tau_sampling_.data(), n, work.data(), 3,
                                    0, 1, _SPLINE_EST_DERIV_);
    array_derive_spline_table_line_to_line(tau_sampling_.data(), n, work.data(), 3,
                                           0, 1, 2);

    for (int i = 0; i < n; i++)
      k2g[i * stride + index_k] -= work[i * 3 + 2];  /* subtract H_T_Nb'' */
  }
```

- [ ] **Step 4: Rebuild and run the consistency test to verify it PASSES**

Run:
```bash
pip install --no-build-isolation . && \
python "<scratchpad>/check_gamma.py"
```
Expected: PASS — prints per-`z` `max rel diff` values below 0.1 and
`GAMMA OK: full k2gamma_Nb reproduces the notebook reconstruction on the low-k plateau`.

- [ ] **Step 5: Visual validation against the published plot**

The notebook's plotting cell currently subtracts the derivative externally
(`k2gamma = tkz['k2gamma_Nb'] - get_H_T_Nb_prime_prime(zval)`). Post-fix the column IS the full
`γ`, so that line would double-subtract. Edit the notebook cell (`ComputeGamma.ipynb`, cell
`b5dca360-...`) to use the column directly:

```python
    k2gamma = tkz['k2gamma_Nb']
```

Run the notebook and compare its `|γ|/Φ` figure to `arXiv-1505.04756v2/gammaPhiRatio.pdf`. Confirm:
- low-`k` plateaus per redshift (`z=1000 ≈ 0.1–0.2`; the `z=49` curve crosses the 1% line near
  `k ≈ 5×10⁻³ h/Mpc`; lower-`z` curves progressively lower),
- the falloff toward large `k`, with oscillatory noise no worse than the published figure.

`ComputeGamma.ipynb` is untracked; do not `git add` it (leave committing/tidying the notebook to the user).

- [ ] **Step 6: Commit (source only)**

```bash
git add source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
Complete N-body gauge gamma: add H_T_Nb'' spline post-process

After each k integration, perturb_solve now splines the stored H_T_Nb'(tau)
time series over tau_sampling_ and subtracts its analytic derivative H_T_Nb''
from the partial k2gamma_Nb, so get_transfer / *_tk.dat return the full
N-body-gauge k^2 gamma. Reproduces arXiv:1505.04756 gammaPhiRatio; the spline
derivative avoids any finite-difference step size. Reuses the array spline
helpers already used for dkappa/dtau.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01Uzbwf8RtwEDNs9vJNE4dDM
EOF
)"
```

---

## Notes for the implementer

- Replace `<scratchpad>` with the literal path from Global Constraints in every command.
- If Step 2 of Task 2 unexpectedly PASSES (partial already matches within 10%), the plateau chosen is not discriminating — tighten the mask toward `kh < 3e-4` and/or lower the tolerance to 3%, where the `H_T_Nb''` contribution to `γ` is unambiguous. Do not weaken the test to make it pass.
- The build recompiles only `perturbations_module.cpp` plus a relink, so each `pip install --no-build-isolation .` is quick.
- Do not touch `python/cclassy.pxd` — it is auto-generated. No wrapper changes are needed (the `k2gamma_Nb` column already exists in the output).
