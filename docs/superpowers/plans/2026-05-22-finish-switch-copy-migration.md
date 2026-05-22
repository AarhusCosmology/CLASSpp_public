# Finish Switch-Copy Migration (N_decay_dr-free species) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the remaining inline "always-reconduct" perturbation slot copies in `perturb_vector_init` into per-species `CopyPerturbationsAcrossSwitch` overrides (Baryons, CDM, Fluid, ScalarField), and wire the `IDM_DR_IDR` / `IDM_DRMD_IDR_DRMD` composites to delegate to their existing child overrides — which restores the `idr`/`idr_drmd`/`idm_drmd` reconduct dropped during #268.

**Architecture:** PR #268 added `BaseSpecies::CopyPerturbationsAcrossSwitch` and a per-species dispatch loop at the top of `perturb_vector_init`'s switch (`else`) branch. Migrated species do their copy through that loop; the rest are still copied inline. This PR finishes the migration for the `N_decay_dr`-free species and composites. The simple-species moves are behavior-neutral (verified against master); the composite delegation restores dropped copies (verified against the pre-#268 reference via two new full-Cl scenarios).

**Tech Stack:** C++17, `make class -j`, `pytest` in `python/`.

**Spec:** `docs/superpowers/specs/2026-05-22-finish-switch-copy-migration-design.md`.

**Worktree:** `/Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration` (branch `worktree-finish-switch-copy-migration`, off merged #268 `d204e3a2`).

---

## Regression cycle (`<regression-check>`)

Run after every behaviour-touching task. Default thread count (do NOT set `OMP_NUM_THREADS=1`).

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration
rm -f output/scen_*.dat                       # 1. clear stale outputs (CRITICAL)
make class -j 2>&1 | tail -10                  # 2. clean build
set -e                                         # 3. run each scenario, check exit code
for f in test/scenarios/ncdm_single.ini \
         test/scenarios/ncdm_multi_unsorted.ini \
         test/scenarios/ncdm_self_interacting.ini \
         test/scenarios/dncdm_dr.ini \
         test/scenarios/ncdm_dncdm_idmdr_combined.ini \
         test/scenarios/idm_dr_full.ini \
         test/scenarios/idm_drmd_full.ini; do
  echo "===== $f ====="; ./class "$f" || { echo "FAIL: $f"; exit 1; }
done
set +e
./class explanatory.ini > /tmp/exp.log 2>&1 || { tail -30 /tmp/exp.log; exit 1; }
```

Then diff with the **tolerance-based, zero-crossing-aware** comparator (NOT bit-identical — this is a
numerical ODE integration; rtol=1e-6 amplifies ULP noise, and Cℓ^TE crosses zero so blind max-relative
is meaningless). Save this script once as `test/scenarios/compare_tol.py`:

```python
import pathlib, numpy as np, sys
# usage: python compare_tol.py <baseline_dir> <new_dir> [glob]
# Tolerance-based, zero-crossing-aware comparator for CLASS .dat outputs.
# A column entry passes if |a-b| <= atol + RTOL*|a|, with atol = RTOL*column_peak.
# Near a zero-crossing (|a|~0) the absolute floor (0.1% of the column's peak
# magnitude) dominates, so meaningless blown-up relative errors at TE/Ephi
# crossings are ignored; peak-region points are judged at ~0.1% relative.
if len(sys.argv) < 3:
    sys.exit("usage: python compare_tol.py <baseline_dir> <new_dir> [glob]")
base, new = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2])
glob = sys.argv[3] if len(sys.argv) > 3 else "*.dat"
RTOL = 1e-3
ok = True
for f in sorted(base.glob(glob)):
    n = new / f.name
    if not n.exists(): print(f"MISSING {f.name}"); ok = False; continue
    # ndmin=2 keeps single-row files 2D so colpeak stays per-column.
    a, b = np.loadtxt(f, comments="#", ndmin=2), np.loadtxt(n, comments="#", ndmin=2)
    if a.shape != b.shape: print(f"SHAPE {f.name}: {a.shape} vs {b.shape}"); ok = False; continue
    if a.size == 0: print(f"OK   {f.name}: empty"); continue
    scale = np.maximum(np.abs(a), np.abs(b))
    colpeak = scale.max(axis=0, keepdims=True)          # per-column peak (always 2D via ndmin=2)
    atol = RTOL * np.where(colpeak > 0, colpeak, 1.0)   # absolute floor = 0.1% of column peak
    diff = np.abs(a - b)
    bad = diff > (atol + RTOL * np.abs(a))
    nbad = int(bad.sum())
    worst = float((diff / np.where(colpeak > 0, colpeak, 1.0)).max())  # worst |a-b|/colpeak
    flag = "OK  " if nbad == 0 else "FAIL"
    if nbad: ok = False
    print(f"{flag} {f.name}: n_exceed={nbad} worst_vs_colpeak={worst:.2e}")
sys.exit(0 if ok else 1)
```

Acceptance per task:
- Five #268 scenarios + `explanatory` → compare `output/scen_*.dat` / relevant outputs against the
  **master baseline** (`/tmp/swcopy_master/`). Expect `relmax <= 1e-3` (simple-species moves are
  exact; existing scenarios don't exercise idm_dr/idm_drmd perturbations).
- `idm_dr_full` / `idm_drmd_full` → compare against the **pre-#268 baseline** (`/tmp/swcopy_pre268/`)
  once the composite delegation lands (Tasks 6–7). Expect `relmax <= 1e-3`.

---

## Task 1: New full-Cl idm_dr / idm_drmd scenarios + capture baselines

**Files:**
- Create: `test/scenarios/idm_dr_full.ini`
- Create: `test/scenarios/idm_drmd_full.ini`
- Create: `test/scenarios/compare_tol.py`

- [ ] **Step 1: Write `test/scenarios/idm_dr_full.ini`**

```ini
# Full-Cl IDM_DR + IDR (free-streaming) scenario. Integrates idm_dr/idr
# perturbations through approximation switches -- coverage absent from the five
# #268 regression scenarios. Used to verify the IDM_DR_IDR composite delegation
# restores the idr reconduct dropped in #268.
h = 0.67556
T_cmb = 2.7255
omega_b = 0.022032
omega_cdm = 0.10
N_ur = 2.0
YHe = BBN
recombination = RECFAST
output = tCl,pCl,lCl,mPk
l_max_scalars = 1000
lensing = yes
# IDM_DR + IDR, free-streaming dark radiation (exercises idr shear/l3 hierarchy)
Omega_idm_dr = 0.01
xi_idr = 0.5
stat_f_idr = 0.875
a_idm_dr = 1e3
idr_nature = free_streaming
input_verbose = 1
perturbations_verbose = 2
root = output/scen_idm_dr_full_
```

- [ ] **Step 2: Write `test/scenarios/idm_drmd_full.ini`**

`f_idm_drmd > 0` brings in IDM_DRMD; `delta_Neff_drmd > 0` brings in IDR_DRMD — together they form the
`IDM_DRMD_IDR_DRMD` composite. (The pytest `test_drmd_without_idr_drmd_computes` uses these keys with
`delta_Neff_drmd = 0`; we set it positive to include the radiation child.)

```ini
# Full-Cl IDM_DRMD + IDR_DRMD scenario. Integrates idm_drmd/idr_drmd
# perturbations through approximation switches. Used to verify the
# IDM_DRMD_IDR_DRMD composite delegation restores the idm_drmd/idr_drmd
# reconduct dropped in #268.
h = 0.67556
T_cmb = 2.7255
omega_b = 0.022032
omega_cdm = 0.10
N_ur = 2.0
YHe = BBN
recombination = RECFAST
output = tCl,pCl,lCl,mPk
l_max_scalars = 1000
lensing = yes
# IDM_DRMD + IDR_DRMD
z_stop = 1.0e4
G_over_aH_drmd_ini = 1.0
f_idm_drmd = 0.1
delta_Neff_drmd = 0.5
input_verbose = 1
perturbations_verbose = 2
root = output/scen_idm_drmd_full_
```

- [ ] **Step 3: Write `test/scenarios/compare_tol.py`**

Copy the comparator script from the `<regression-check>` section above into this file.

- [ ] **Step 4: Build current branch (== master, no code changes yet) and capture master baseline**

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration
make class -j 2>&1 | tail -3
rm -f output/scen_*.dat
mkdir -p /tmp/swcopy_master
set -e
for f in test/scenarios/*.ini; do echo "== $f =="; ./class "$f" || echo "  (exit non-zero on $f)"; done
set +e
cp output/scen_*.dat /tmp/swcopy_master/ 2>/dev/null
ls /tmp/swcopy_master/
```

Expected: the five #268 scenarios + `explanatory` succeed. The `idm_dr_full` / `idm_drmd_full` runs
complete and write `scen_idm_dr_full_*` / `scen_idm_drmd_full_*` — but with the **dropped** reconduct
(idr/idr_drmd/idm_drmd reset at each switch). Record whether they error; their outputs here are the
"buggy master" reference (they should later DIFFER from the pre-#268 reference, proving coverage).

- [ ] **Step 5: Capture pre-#268 baseline for the two new scenarios (throwaway worktree)**

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration
git worktree add /tmp/class-pre268 d204e3a2^
cp test/scenarios/idm_dr_full.ini test/scenarios/idm_drmd_full.ini /tmp/class-pre268/test/scenarios/ 2>/dev/null || \
  { mkdir -p /tmp/class-pre268/test/scenarios && cp test/scenarios/idm_dr_full.ini test/scenarios/idm_drmd_full.ini /tmp/class-pre268/test/scenarios/; }
cd /tmp/class-pre268
make class -j 2>&1 | tail -3
rm -f output/scen_*.dat
./class test/scenarios/idm_dr_full.ini   || echo "pre268 idm_dr_full exit non-zero"
./class test/scenarios/idm_drmd_full.ini || echo "pre268 idm_drmd_full exit non-zero"
mkdir -p /tmp/swcopy_pre268
cp output/scen_idm_dr_full_*.dat output/scen_idm_drmd_full_*.dat /tmp/swcopy_pre268/ 2>/dev/null
ls /tmp/swcopy_pre268/
cd /Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration
git worktree remove /tmp/class-pre268 --force
```

Expected: pre-#268 builds and runs both new scenarios (it performed the idr/idr_drmd reconduct). These
outputs are the **correct reference**. Sanity check the regression is real:

```bash
python test/scenarios/compare_tol.py /tmp/swcopy_pre268 /tmp/swcopy_master "scen_idm_dr_full_*.dat"
python test/scenarios/compare_tol.py /tmp/swcopy_pre268 /tmp/swcopy_master "scen_idm_drmd_full_*.dat"
```

Expected: at least one **FAIL** (master differs from pre-#268), confirming the new scenarios actually
exercise the dropped reconduct. If everything is `OK`, the scenarios don't hit the switch path — raise
`l_max_scalars`, confirm `perturbations_verbose=2` prints "switch ... approximation" lines, and adjust
before proceeding. (If they're OK because the transitions genuinely never co-occur, note it and flag to
the user — the "regression" may be unreachable.)

- [ ] **Step 6: Commit**

```bash
git add test/scenarios/idm_dr_full.ini test/scenarios/idm_drmd_full.ini test/scenarios/compare_tol.py
git commit -m "test: full-Cl idm_dr & idm_drmd scenarios + tolerance comparator"
```

---

## Task 2: BaryonsSpecies CopyPerturbationsAcrossSwitch override

**Files:**
- Modify: `species/baryons.h`
- Modify: `species/baryons.cpp`

Replaces inline block (perturb_vector_init section a.2): `{ const size_t b_sw_i = all_species_.index_of("Baryons"); ... ppv->y[b_new.idx_delta] = ...; ppv->y[b_new.idx_theta] = ...; }`. (Inline deleted in Task 8.)

- [ ] **Step 1: Declare the override in `species/baryons.h`**

In the `public:` section of `class BaryonsSpecies`, immediately after the
`void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout, ...)` declaration (~line 96):

```cpp
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 2: Define it in `species/baryons.cpp`** (append at file scope)

```cpp
void BaryonsSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                   const BaseSpecies::PerturbLayout& new_base,
                                                   const double* old_y,
                                                   double* new_y,
                                                   const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  // Scalar-mode only: idx_delta is registered (>=0) in scalar pv and -1 in
  // vector/tensor pv, so this no-ops outside scalar mode (the vector-mode
  // baryon theta copy stays inline -- out of scope for this PR).
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}
```

- [ ] **Step 3: Build** — `make class -j 2>&1 | tail -5`. Expected: clean.

- [ ] **Step 4: Run `<regression-check>`; diff existing scenarios + explanatory vs master**

```bash
python test/scenarios/compare_tol.py /tmp/swcopy_master output "scen_ncdm_*.dat"
python test/scenarios/compare_tol.py /tmp/swcopy_master output "scen_combined_*.dat"
```

Expected: all `OK` (the override duplicates the still-present inline copy → identical values).

- [ ] **Step 5: Commit**

```bash
git add species/baryons.h species/baryons.cpp
git commit -m "baryons: CopyPerturbationsAcrossSwitch override (delta/theta)"
```

---

## Task 3: CDMSpecies CopyPerturbationsAcrossSwitch override

**Files:**
- Modify: `species/cdm.h`
- Modify: `species/cdm.cpp`

Replaces inline block: `if (all_species_.count("CDM")) { ... ppv->y[cdm_new_lay.idx_delta] = ...; if (ppt->gauge == newtonian) ppv->y[cdm_new_lay.idx_theta] = ...; }`. CDM's layout registers `idx_theta` only in newtonian gauge, so the gauge condition is expressed as `idx_theta >= 0`.

- [ ] **Step 1: Declare in `species/cdm.h`** (after `ApplyInitialConditions(layout,...)` ~line 120)

```cpp
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 2: Define in `species/cdm.cpp`** (append at file scope)

```cpp
void CDMSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                               const BaseSpecies::PerturbLayout& new_base,
                                               const double* old_y,
                                               double* new_y,
                                               const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0)
    return;
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  // theta is registered only in newtonian gauge (idx_theta >= 0); synchronous
  // gauge derives cdm velocity from the metric, so idx_theta == -1 there.
  if (old_l.idx_theta >= 0 && new_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
}
```

- [ ] **Step 3: Build** — `make class -j 2>&1 | tail -5`. Expected: clean.

- [ ] **Step 4: Run `<regression-check>`; diff vs master.** Expected: all `OK`.

- [ ] **Step 5: Commit**

```bash
git add species/cdm.h species/cdm.cpp
git commit -m "cdm: CopyPerturbationsAcrossSwitch override (delta; theta in newtonian)"
```

---

## Task 4: FluidSpecies CopyPerturbationsAcrossSwitch override

**Files:**
- Modify: `species/fluid.h`
- Modify: `species/fluid.cpp`

Replaces inline block: `if (all_species_.count("Fluid")) { if (pba->use_ppf == _FALSE_) { copy idx_delta, idx_theta } else { copy idx_Gamma } }`. The layout registers `idx_delta`/`idx_theta` in standard mode and `idx_Gamma` in PPF mode (the unused ones are -1), so per-slot `idx >= 0` guards reproduce the `use_ppf` branch without `ctx`. **Note:** no early `return` on `idx_delta` here — PPF mode has `idx_delta == -1` but a valid `idx_Gamma`.

- [ ] **Step 1: Declare in `species/fluid.h`** (after `ApplyInitialConditions(layout,...)` ~line 110)

```cpp
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 2: Define in `species/fluid.cpp`** (append at file scope)

```cpp
void FluidSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                 const BaseSpecies::PerturbLayout& new_base,
                                                 const double* old_y,
                                                 double* new_y,
                                                 const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  // Standard fluid registers idx_delta/idx_theta; PPF registers idx_Gamma. The
  // unused slots are -1. Vector/tensor pv: all -1 -> no-op.
  if (old_l.idx_delta >= 0 && new_l.idx_delta >= 0)
    new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  if (old_l.idx_theta >= 0 && new_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
  if (old_l.idx_Gamma >= 0 && new_l.idx_Gamma >= 0)
    new_y[new_l.idx_Gamma] = old_y[old_l.idx_Gamma];
}
```

- [ ] **Step 3: Build** — `make class -j 2>&1 | tail -5`. Expected: clean.

- [ ] **Step 4: Run `<regression-check>`; diff vs master.** Expected: all `OK`. (The five scenarios use
  no Fluid, but `explanatory.ini` does — confirm it still completes; the override no-ops where Fluid is
  absent.)

- [ ] **Step 5: Commit**

```bash
git add species/fluid.h species/fluid.cpp
git commit -m "fluid: CopyPerturbationsAcrossSwitch override (delta/theta or PPF Gamma)"
```

---

## Task 5: ScalarFieldSpecies CopyPerturbationsAcrossSwitch override

**Files:**
- Modify: `species/scalar_field.h`
- Modify: `species/scalar_field.cpp`

Replaces inline block: `if (all_species_.count("ScalarField")) { ppv->y[scf_new_lay.idx_phi] = ...; ppv->y[scf_new_lay.idx_phi_prime] = ...; }`.

- [ ] **Step 1: Declare in `species/scalar_field.h`** (after `ApplyInitialConditions(layout,...)` ~line 100)

```cpp
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 2: Define in `species/scalar_field.cpp`** (append at file scope)

```cpp
void ScalarFieldSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                       const BaseSpecies::PerturbLayout& new_base,
                                                       const double* old_y,
                                                       double* new_y,
                                                       const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_phi < 0 || new_l.idx_phi < 0)
    return;
  new_y[new_l.idx_phi]       = old_y[old_l.idx_phi];
  new_y[new_l.idx_phi_prime] = old_y[old_l.idx_phi_prime];
}
```

- [ ] **Step 3: Build** — `make class -j 2>&1 | tail -5`. Expected: clean.

- [ ] **Step 4: Run `<regression-check>`; diff vs master.** Expected: all `OK`. (`explanatory.ini` does
  not use a scalar field by default; the override no-ops where absent. If you want positive coverage,
  confirm a scf run still completes — out of the required scope.)

- [ ] **Step 5: Commit**

```bash
git add species/scalar_field.h species/scalar_field.cpp
git commit -m "scalar_field: CopyPerturbationsAcrossSwitch override (phi/phi_prime)"
```

---

## Task 6: Wire IDM_DR_IDR_Species composite delegation (restores idr reconduct)

**Files:**
- Modify: `species/idm_dr_idr_species.h`
- Modify: `species/idm_dr_idr_species.cpp`

The composite's nested layout is `PerturbLayout { IDM_DRSpecies::PerturbLayout idm_dr; IDRSpecies::PerturbLayout idr; }`; children are `idm_dr_` / `idr_`. The child overrides already exist (`IDM_DRSpecies` copies delta/theta; `IDRSpecies` copies delta/theta/shear/l3). Delegating activates them — restoring the idr reconduct dropped in #268. (The inline idm_dr delta/theta copy is removed in Task 8.)

- [ ] **Step 1: Declare the override in `species/idm_dr_idr_species.h`** (public section, after `CreatePerturbLayout`)

```cpp
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 2: Define in `species/idm_dr_idr_species.cpp`** (append at file scope)

```cpp
void IDM_DR_IDR_Species::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                       const BaseSpecies::PerturbLayout& new_base,
                                                       const double* old_y,
                                                       double* new_y,
                                                       const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  idm_dr_->CopyPerturbationsAcrossSwitch(old_l.idm_dr, new_l.idm_dr, old_y, new_y, ctx);
  idr_->CopyPerturbationsAcrossSwitch(old_l.idr, new_l.idr, old_y, new_y, ctx);
}
```

- [ ] **Step 3: Build** — `make class -j 2>&1 | tail -5`. Expected: clean.

- [ ] **Step 4: Run `<regression-check>`; verify the new idm_dr scenario against the pre-#268 reference**

```bash
python test/scenarios/compare_tol.py /tmp/swcopy_pre268 output "scen_idm_dr_full_*.dat"
# And confirm idm_dr now DIFFERS from buggy master (the fix is active):
python test/scenarios/compare_tol.py /tmp/swcopy_master output "scen_idm_dr_full_*.dat"
```

Expected: vs **pre-#268** → all `OK` (restored idr reconduct matches the correct reference, within
0.1%). vs **master** → at least one `FAIL` (we changed behavior — that is the regression fix). The five
#268 scenarios + explanatory remain `OK` vs master.

- [ ] **Step 5: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp
git commit -m "idm_dr_idr: delegate CopyPerturbationsAcrossSwitch to children (restores idr reconduct dropped in #268)"
```

---

## Task 7: Wire IDM_DRMD_IDR_DRMD_Species composite delegation

**Files:**
- Modify: `species/idm_drmd_idr_drmd_species.h`
- Modify: `species/idm_drmd_idr_drmd_species.cpp`

Nested layout `PerturbLayout { IDM_DRMDSpecies::PerturbLayout idm_drmd; IDR_DRMDSpecies::PerturbLayout idr_drmd; }`; children `idm_drmd_` / `idr_drmd_`. Both child overrides copy delta/theta. The current code copies idm_drmd/idr_drmd delta/theta **only** inside the tca_idm_drmd-off block; delegating restores the general always-reconduct.

- [ ] **Step 1: Declare in `species/idm_drmd_idr_drmd_species.h`** (public section, after `CreatePerturbLayout`)

```cpp
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
```

- [ ] **Step 2: Define in `species/idm_drmd_idr_drmd_species.cpp`** (append at file scope)

```cpp
void IDM_DRMD_IDR_DRMD_Species::CopyPerturbationsAcrossSwitch(
    const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base,
    const double* old_y,
    double* new_y,
    const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  idm_drmd_->CopyPerturbationsAcrossSwitch(old_l.idm_drmd, new_l.idm_drmd, old_y, new_y, ctx);
  idr_drmd_->CopyPerturbationsAcrossSwitch(old_l.idr_drmd, new_l.idr_drmd, old_y, new_y, ctx);
}
```

- [ ] **Step 3: Build** — `make class -j 2>&1 | tail -5`. Expected: clean.

- [ ] **Step 4: Run `<regression-check>`; verify the new idm_drmd scenario against pre-#268**

```bash
python test/scenarios/compare_tol.py /tmp/swcopy_pre268 output "scen_idm_drmd_full_*.dat"
python test/scenarios/compare_tol.py /tmp/swcopy_master output "scen_idm_drmd_full_*.dat"
```

Expected: vs **pre-#268** → all `OK`. vs **master** → expect a `FAIL` if master dropped the reconduct
for the transitions this scenario hits (regression fix active). If idm_drmd matches master *and*
pre-#268 (no transition exercised), note it — the DRMD reconduct may be unreachable in this config; in
that case the delegation is still correct (matches pre-#268) but unverified for behavior change.

- [ ] **Step 5: Commit**

```bash
git add species/idm_drmd_idr_drmd_species.h species/idm_drmd_idr_drmd_species.cpp
git commit -m "idm_drmd_idr_drmd: delegate CopyPerturbationsAcrossSwitch to children (restores reconduct dropped in #268)"
```

---

## Task 8: Delete the now-redundant inline section-a.2 copy blocks

**Files:**
- Modify: `source/perturbations_module.cpp`

The dispatch loop (top of the switch `else` branch) now performs all of these via species/composite
overrides. Delete the corresponding inline blocks in section a.2 (the "some variables ... are not
affected by any approximation" block, ~lines 3856–3935). **Do NOT delete:** the DCDM_DR block
(`if (all_species_.count("DCDM_DR")) { ... }` incl. the `N_decay_dr` DR-multipole loop — PR B), the
metric copies (`index_pt_eta` / `index_pt_phi`), perturbed-recombination, or the vector-mode baryon
(`b_vsw`) copy.

- [ ] **Step 1: Delete the Baryons inline block**

```cpp
      {
        const size_t b_sw_i = all_species_.index_of("Baryons");
        const auto& b_old   = static_cast<const BaryonsSpecies::PerturbLayout&>(
            *ppw->pv->species_layouts[b_sw_i]);
        auto& b_new = static_cast<BaryonsSpecies::PerturbLayout&>(*ppv->species_layouts[b_sw_i]);
        ppv->y[b_new.idx_delta] = ppw->pv->y[b_old.idx_delta];
        ppv->y[b_new.idx_theta] = ppw->pv->y[b_old.idx_theta];
      }
```

- [ ] **Step 2: Delete the CDM inline block**

```cpp
      if (all_species_.count("CDM")) {
        const size_t cdm_i      = all_species_.index_of("CDM");
        const auto& cdm_old_lay = static_cast<const CDMSpecies::PerturbLayout&>(
            *ppw->pv->species_layouts[cdm_i]);
        auto& cdm_new_lay = static_cast<CDMSpecies::PerturbLayout&>(*ppv->species_layouts[cdm_i]);
        ppv->y[cdm_new_lay.idx_delta] = ppw->pv->y[cdm_old_lay.idx_delta];
        if (ppt->gauge == newtonian) {
          ppv->y[cdm_new_lay.idx_theta] = ppw->pv->y[cdm_old_lay.idx_theta];
        }
      }
```

- [ ] **Step 3: Delete the IDM_DR inline block**

```cpp
      if (all_species_.count("IDM_DR_IDR")) {
        const size_t idm_dr_sw_i         = all_species_.index_of("IDM_DR_IDR");
        const auto& old_idm_dr_lay       = static_cast<const IDM_DR_IDR_Species::PerturbLayout&>(
                                               *ppw->pv->species_layouts[idm_dr_sw_i])
                                               .idm_dr;
        const auto& new_idm_dr_lay       = static_cast<const IDM_DR_IDR_Species::PerturbLayout&>(
                                               *ppv->species_layouts[idm_dr_sw_i])
                                               .idm_dr;
        ppv->y[new_idm_dr_lay.idx_delta] = ppw->pv->y[old_idm_dr_lay.idx_delta];
        ppv->y[new_idm_dr_lay.idx_theta] = ppw->pv->y[old_idm_dr_lay.idx_theta];
      }
```

- [ ] **Step 4: Delete the Fluid inline block**

```cpp
      if (all_species_.count("Fluid")) {
        const size_t fld_i      = all_species_.index_of("Fluid");
        const auto& fld_old_lay = static_cast<const FluidSpecies::PerturbLayout&>(
            *ppw->pv->species_layouts[fld_i]);
        auto& fld_new_lay = static_cast<FluidSpecies::PerturbLayout&>(*ppv->species_layouts[fld_i]);
        if (pba->use_ppf == _FALSE_) {
          ppv->y[fld_new_lay.idx_delta] = ppw->pv->y[fld_old_lay.idx_delta];

          ppv->y[fld_new_lay.idx_theta] = ppw->pv->y[fld_old_lay.idx_theta];
        }
        else {
          ppv->y[fld_new_lay.idx_Gamma] = ppw->pv->y[fld_old_lay.idx_Gamma];
        }
      }
```

- [ ] **Step 5: Delete the ScalarField inline block**

```cpp
      if (all_species_.count("ScalarField")) {
        const size_t scf_i      = all_species_.index_of("ScalarField");
        const auto& scf_old_lay = static_cast<const ScalarFieldSpecies::PerturbLayout&>(
            *ppw->pv->species_layouts[scf_i]);
        auto& scf_new_lay = static_cast<ScalarFieldSpecies::PerturbLayout&>(
            *ppv->species_layouts[scf_i]);
        ppv->y[scf_new_lay.idx_phi]       = ppw->pv->y[scf_old_lay.idx_phi];
        ppv->y[scf_new_lay.idx_phi_prime] = ppw->pv->y[scf_old_lay.idx_phi_prime];
      }
```

- [ ] **Step 6: Confirm the surviving inline content is intact**

```bash
grep -n "index_pt_eta\|index_pt_phi\b\|all_species_.count(\"DCDM_DR\")\|b_vsw" source/perturbations_module.cpp | head
```

Expected: the metric (`index_pt_eta`/`index_pt_phi`), DCDM_DR, and vector-mode `b_vsw` blocks are still
present. Also verify the deleted-species layout types are no longer referenced in the a.2 region:

```bash
awk 'NR>=3840 && NR<=3960' source/perturbations_module.cpp | grep -nE "BaryonsSpecies::PerturbLayout|CDMSpecies::PerturbLayout|FluidSpecies::PerturbLayout|ScalarFieldSpecies::PerturbLayout|\.idm_dr" || echo "clean: no leftover simple-species/idm_dr inline copies"
```

- [ ] **Step 7: Build** — `make class -j 2>&1 | tail -5`. Expected: clean (no unused-variable warnings
  from half-deleted blocks).

- [ ] **Step 8: Run full `<regression-check>` + both comparators**

```bash
python test/scenarios/compare_tol.py /tmp/swcopy_master output "scen_ncdm_*.dat"
python test/scenarios/compare_tol.py /tmp/swcopy_master output "scen_combined_*.dat"
python test/scenarios/compare_tol.py /tmp/swcopy_pre268 output "scen_idm_dr_full_*.dat"
python test/scenarios/compare_tol.py /tmp/swcopy_pre268 output "scen_idm_drmd_full_*.dat"
```

Expected: existing scenarios `OK` vs master; new scenarios `OK` vs pre-#268.

- [ ] **Step 9: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: delete inline switch-copy for Baryons/CDM/Fluid/ScalarField/IDM_DR (now species hooks)"
```

---

## Task 9: Final verification + finish branch

**Files:** none (verification + PR).

- [ ] **Step 1: Clean rebuild + full `<regression-check>`**

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration
make clean && make class -j 2>&1 | tail -10
# then the full <regression-check> loop + both comparators (existing vs master, new vs pre-#268)
```

- [ ] **Step 2: pytest scenario suite**

```bash
cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
```

Expected: all pass.

- [ ] **Step 3: Confirm the inline switch branch is reduced as intended**

```bash
cd /Users/au192734/Projects/class_claude/.claude/worktrees/finish-switch-copy-migration
# Simple-species + idm_dr inline copies gone from the a.2 region; DCDM_DR / metric / b_vsw remain:
grep -n "CopyPerturbationsAcrossSwitch" species/baryons.cpp species/cdm.cpp species/fluid.cpp \
   species/scalar_field.cpp species/idm_dr_idr_species.cpp species/idm_drmd_idr_drmd_species.cpp
```

Expected: each file has the new override.

- [ ] **Step 4: Finish the branch**

Use the superpowers:finishing-a-development-branch skill to open the PR. PR body must call out: (a) the
four simple-species CopyPerturbationsAcrossSwitch overrides (behavior-neutral relocations, verified vs
master); (b) the IDM_DR_IDR / IDM_DRMD_IDR_DRMD composite delegation **restoring idr/idr_drmd/idm_drmd
reconduct dropped in #268** (behavior change, verified vs the pre-#268 reference via the two new
full-Cl scenarios); (c) DR composites + `N_decay_dr` remain for PR B; (d) verification is
tolerance-based (~0.1%, zero-crossing-aware), not bit-identical.

---

## Plan summary

9 tasks: (1) new idm_dr/idm_drmd full-Cl scenarios + tolerance comparator + master & pre-#268
baselines; (2–5) Baryons / CDM / Fluid / ScalarField `CopyPerturbationsAcrossSwitch` overrides
(behavior-neutral); (6–7) IDM_DR_IDR / IDM_DRMD_IDR_DRMD composite delegation (restores #268-dropped
reconduct); (8) delete the redundant inline a.2 blocks; (9) final verification + PR.

**No build-system changes:** all tasks modify existing `.cpp`/`.h` files; no new source files, so the
Makefile / setup.py / Xcode project lists are untouched. The new files are `.ini` scenarios and a
`.py` comparator (not compiled).

## Spec coverage check

| Spec section | Task(s) |
|---|---|
| Goal 1 (Baryons/CDM/Fluid/ScalarField overrides) | 2, 3, 4, 5 |
| Goal 2 (new full-Cl idm_dr/idm_drmd scenarios) | 1 |
| Goal 3 (composite delegation; restores dropped reconduct) | 6, 7 |
| Goal 4 (delete inline a.2 blocks) | 8 |
| Behavior-change boundary (master vs pre-#268 baselines) | 1 (baselines), 2–5 (vs master), 6–7 (vs pre-#268) |
| Verification (tolerance-based, zero-crossing-aware) | `<regression-check>` + `compare_tol.py`, every task |
| Non-goal: DR composites / N_decay_dr untouched | 8 Step 6 (explicit "do not delete") |
