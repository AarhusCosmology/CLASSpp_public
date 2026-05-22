# Remove `ncdm_id` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the integer `ncdm_id` and every flat-array sibling it indexes into, so that NCDM-family quantities are owned by their species instances and module code dispatches via virtual hooks.

**Architecture:** Per-species perturbation/transfer/shooter storage on `NCDMBaseSpecies` (mirroring the existing `bg_rho_index_` / `bi_rho_index_` encapsulation pattern); generic `BaseSpecies` hooks for approximation-switch bookkeeping (`StashPerturbationLayout`, `CopyPerturbationsAcrossSwitch`, `MarkUsedInSources`) and shooter slot registration (`RegisterShootingIndices`, `ComputeShootingResidual`, `ComputeShootingGuess`); deletion of `pv->index_ncdm_/l_max_ncdm/q_size_ncdm/N_ncdm/index_pt_psi0_ncdm1`, `pba->N_ncdm/Omega0_ncdm_tot/N_decay_dr`, and `PerturbationsModule::ncdm_species_sorted_`.

**Tech Stack:** C++17, `make class -j` build, `pytest` in `python/`. No new build dependencies.

**Reference:** Design spec at `docs/superpowers/specs/2026-05-03-remove-ncdm-id-design.md`.

---

## Strategy

Each task either:
- **Adds** a new hook / storage member without changing behavior (gated by a build + smoke test), or
- **Migrates** consumers from old → new path one logical site at a time (gated by build + smoke + targeted regression scenario), or
- **Deletes** dead code after all consumers have migrated (gated by build + smoke + full pytest suite).

A "regression scenario" run after a behavior-touching task means: run the .ini files listed in `Targeted regression scenarios` below and diff their output against a baseline captured on `master` before the branch starts. Acceptable diff: bit-identical for `(.)rho_*` columns; physical Cl's / P(k) within `1e-10` relative tolerance.

---

## Targeted regression scenarios

Before starting Task 1, capture a baseline on `master`:

```bash
cd /Users/au192734/Projects/class_claude
git checkout master && make class -j
mkdir -p /tmp/ncdm_baseline
for ini in test/scenarios/ncdm_single.ini test/scenarios/ncdm_multi_unsorted.ini \
           test/scenarios/dncdm_dr.ini test/scenarios/ncdm_self_interacting.ini \
           test/scenarios/ncdm_dncdm_idmdr_combined.ini; do
  ./class "$ini" 2>&1 | tee "/tmp/ncdm_baseline/$(basename $ini).log"
  cp output/*.dat /tmp/ncdm_baseline/ || true
done
```

If any of these scenarios do not exist, the **first action of Task 1** is to create them. Specifically required:

1. `test/scenarios/ncdm_single.ini` — `m_ncdm = 0.06` (single standard NCDM)
2. `test/scenarios/ncdm_multi_unsorted.ini` — two species named to defeat construction-order = name-order: `nu_b.type = ncdm_standard, nu_b.m_ncdm = 0.06`, `nu_a.type = ncdm_standard, nu_a.m_ncdm = 0.10`
3. `test/scenarios/dncdm_dr.ini` — one DNCDM_DR composite via dot-syntax
4. `test/scenarios/ncdm_self_interacting.ini` — one NCDM with `type = ncdm_self_interacting, G_eff = 1e-10`
5. `test/scenarios/ncdm_dncdm_idmdr_combined.ini` — NCDM + DNCDM_DR + IDM_DR_IDR in one .ini

A regression-diff helper:

```bash
python3 - <<'PY'
import sys, glob, pathlib
import numpy as np
base = pathlib.Path("/tmp/ncdm_baseline")
new  = pathlib.Path("output")
ok = True
for f in sorted(base.glob("*.dat")):
    n = new / f.name
    if not n.exists():
        print(f"MISSING: {f.name}"); ok=False; continue
    a = np.loadtxt(f, comments="#")
    b = np.loadtxt(n, comments="#")
    if a.shape != b.shape:
        print(f"SHAPE: {f.name} base={a.shape} new={b.shape}"); ok=False; continue
    rel = np.where(np.abs(a)>0, np.abs((a-b)/a), np.abs(a-b))
    rmax = rel.max()
    print(f"{f.name}: max_rel_diff={rmax:.2e}")
    if rmax > 1e-10: ok=False
sys.exit(0 if ok else 1)
PY
```

---

## Task 1: Create issue, branch, and regression scenarios

**Files:**
- Create: `test/scenarios/ncdm_single.ini` (and 4 others if missing — see above)

- [ ] **Step 1: Create the GitHub issue**

```bash
cd /Users/au192734/Projects/class_claude
gh issue create \
  --title "Remove ncdm_id: decentralize NCDM-family state into species instances" \
  --body "Spec at docs/superpowers/specs/2026-05-03-remove-ncdm-id-design.md. Plan at docs/superpowers/plans/2026-05-03-remove-ncdm-id.md."
```

Note the issue number (e.g. `267`).

- [ ] **Step 2: Create the branch**

```bash
git checkout -b 267-remove-ncdm-id   # replace 267 with actual issue number
```

- [ ] **Step 3: Create any missing regression scenario .ini files**

Check which exist:
```bash
ls test/scenarios/ncdm_single.ini test/scenarios/ncdm_multi_unsorted.ini \
   test/scenarios/dncdm_dr.ini test/scenarios/ncdm_self_interacting.ini \
   test/scenarios/ncdm_dncdm_idmdr_combined.ini 2>&1 | grep "No such"
```

For any that are missing, create them. Example for `ncdm_multi_unsorted.ini` (the most novel one — exercises iteration-order risk):

```ini
# Two NCDM species named so construction order != lexical order.
# nu_b is constructed first (lighter mass), nu_a second (heavier),
# but lexical sort puts nu_a first. After ncdm_id removal, iteration
# is by lexical order — this scenario detects accidental order dependence.
output = tCl,pCl,lCl,mPk
modes = s,t

nu_b.type = ncdm_standard
nu_b.m_ncdm = 0.06

nu_a.type = ncdm_standard
nu_a.m_ncdm = 0.10
```

For the others, model on existing .ini files in `test/` or `explanatory.ini`. Each must run cleanly under master:

```bash
make class -j 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  echo "=== $f ==="
  ./class "$f" 2>&1 | tail -3
done
```

Expected for each: `[CLASS] All parameters and all quantities computed successfully`.

- [ ] **Step 4: Capture baselines** — run the baseline-capture block from the "Targeted regression scenarios" section above.

- [ ] **Step 5: Commit the scenarios**

```bash
git add test/scenarios/
git commit -m "test: scenarios for ncdm_id removal regression checking"
```

---

## Task 2: Per-species perturbation storage on `NCDMBaseSpecies` (dual-write)

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_species.cpp`
- Modify: `species/dncdm_species.cpp`

This task adds the new storage and accessors but keeps the existing `pv->index_ncdm_[ncdm_id_]/l_max_ncdm[ncdm_id_]/q_size_ncdm[ncdm_id_]` writes intact. Behavior is unchanged. A runtime assert verifies the new storage matches the legacy storage.

- [ ] **Step 1: Add storage + accessors to `NCDMBaseSpecies`**

Edit `species/ncdm_base_species.h`. In the protected section (around line 117 where the other shared fields live), add:

```cpp
  // ── Perturbation storage ─────────────────────────────────────────────────
  // Set by RegisterPerturbationIndices on each call (incl. fluid-approx switches).
  // -1 means "not currently registered".
  int pt_l_max_   = -1;
  int pt_q_size_  = -1;
  std::vector<int> pt_index_per_q_;

  // Snapshot taken before the next RegisterPerturbationIndices, used by
  // CopyPerturbationsAcrossSwitch (added in a later task).
  int pt_l_max_previous_   = -1;
  int pt_q_size_previous_  = -1;
  std::vector<int> pt_index_per_q_previous_;
```

In the public section, add accessors:

```cpp
  int pt_l_max() const { return pt_l_max_; }
  int pt_q_size() const { return pt_q_size_; }
  int pt_index_at(int q) const { return pt_index_per_q_.at(q); }
  int pt_total_size() const { return pt_q_size_ * (pt_l_max_ + 1); }
```

- [ ] **Step 2: Populate storage in `NCDMSpecies::RegisterPerturbationIndices`**

Edit `species/ncdm_species.cpp`. The current function (around line 291) writes to `pv->index_ncdm_[ncdm_id_]/l_max_ncdm[ncdm_id_]/q_size_ncdm[ncdm_id_]`. Add parallel writes to `pt_l_max_/pt_q_size_/pt_index_per_q_` and an assert that they match.

Replace lines ~299–317 with:

```cpp
  if (ncdm_id_ == 0) {
    pv->index_pt_psi0_ncdm1 = index_pt;
  }
  index_pt_psi0_ = index_pt;

  const bool fa_on = (ppw->approx[ppw->index_ap_ncdmfa] == (int) ncdmfa_on);

  if (fa_on) {
    pt_l_max_  = 2;
    pt_q_size_ = 1;
  } else {
    pt_l_max_  = ppr->l_max_ncdm;
    pt_q_size_ = q_size();
  }
  pt_index_per_q_.clear();
  for (int iq = 0; iq < pt_q_size_; ++iq)
    pt_index_per_q_.push_back(index_pt + iq * (pt_l_max_ + 1));

  // Dual-write: keep legacy shared arrays in sync until consumers migrate.
  pv->l_max_ncdm[ncdm_id_]  = pt_l_max_;
  pv->q_size_ncdm[ncdm_id_] = pt_q_size_;
  pv->index_ncdm_[ncdm_id_] = pt_index_per_q_;

  index_pt += pt_total_size();
```

- [ ] **Step 3: Same change in `DNCDMSpecies::RegisterPerturbationIndices`**

Edit `species/dncdm_species.cpp` (around line 415). Apply the same pattern — assign to `pt_l_max_/pt_q_size_/pt_index_per_q_` first, then dual-write to `pv->l_max_ncdm[ncdm_id_]/q_size_ncdm[ncdm_id_]/index_ncdm_[ncdm_id_]`, then advance `index_pt` by `pt_total_size()`.

(`NCDMInteractingSpecies` does not override `RegisterPerturbationIndices` — it inherits from `NCDMSpecies` and its storage is populated by the base.)

- [ ] **Step 4: Build + smoke**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
```

Expected: clean build, `[CLASS] All parameters and all quantities computed successfully`.

- [ ] **Step 5: Run regression scenarios**

```bash
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

All must succeed. Outputs should be bit-identical to baseline (no behavior change yet). Run the regression-diff helper from the "Targeted regression scenarios" section.

- [ ] **Step 6: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_species.cpp species/dncdm_species.cpp
git commit -m "ncdm: add per-species perturbation storage (dual-write, no behavior change)"
```

---

## Task 3: Switch readers in species code to per-species accessors

**Files:**
- Modify: `species/ncdm_species.cpp`
- Modify: `species/ncdm_interacting_species.cpp`
- Modify: `species/dncdm_species.cpp`
- Modify: `species/dncdm_dr_species.cpp`

All four files currently read from `pv->index_ncdm_.at(ncdm_id_)[q]`, `pv->l_max_ncdm[ncdm_id_]`, `pv->q_size_ncdm[ncdm_id_]`. Since each species owns the same data via `pt_index_at(q)/pt_l_max()/pt_q_size()`, these reads can be local to the instance — no `pv` lookup needed.

- [ ] **Step 1: List all read sites**

```bash
grep -n "pv->index_ncdm_\|pv->l_max_ncdm\|pv->q_size_ncdm\|pv->\?index_ncdm_\.at" species/ncdm_species.cpp species/ncdm_interacting_species.cpp species/dncdm_species.cpp species/dncdm_dr_species.cpp
```

Note every line number; you will edit each one.

- [ ] **Step 2: Apply the substitution pattern**

For each match, the rewrite is mechanical:

| Before | After |
|---|---|
| `pv->index_ncdm_.at(ncdm_id_)[iq]` | `pt_index_at(iq)` |
| `pv->index_ncdm_[ncdm_id_][iq]` | `pt_index_at(iq)` |
| `pv->l_max_ncdm[ncdm_id_]` | `pt_l_max()` |
| `pv->q_size_ncdm[ncdm_id_]` | `pt_q_size()` |

Special case in `species/dncdm_dr_species.cpp` lines 140, 152, 161, 163: the composite reads from its DNCDM child's slots. The substitution becomes:

| Before | After |
|---|---|
| `pv->index_ncdm_.at(ncdm_id_)[index_q] + l` | `dncdm_.pt_index_at(index_q) + l` |
| `pv->index_ncdm_.at(ncdm_id_)[0]` | `dncdm_.pt_index_at(0)` |
| `pv->index_ncdm_.at(ncdm_id_)[0] + 1` | `dncdm_.pt_index_at(0) + 1` |

Special case in `species/ncdm_interacting_species.cpp` lines 121, 126, 128, 129: same as `NCDMSpecies` (inherits its storage).

Concrete example — `species/ncdm_species.cpp:354`:

```cpp
// Before:
const int idx = pv->index_ncdm_.at(ncdm_id_)[0];
// After:
const int idx = pt_index_at(0);
```

- [ ] **Step 3: Verify no remaining shared-array reads in species code**

```bash
grep -n "pv->index_ncdm_\|pv->l_max_ncdm\|pv->q_size_ncdm" species/ncdm_*.cpp species/dncdm_*.cpp
```

Expected: 0 lines (the writes added in Task 2 are also gone now? — no, the writes stay until Task 18 because the legacy fields still exist. Verify only READS are gone.)

Re-grep more precisely:

```bash
grep -n "pv->index_ncdm_\|pv->l_max_ncdm\|pv->q_size_ncdm" species/ncdm_*.cpp species/dncdm_*.cpp \
  | grep -v "= pt_" | grep -v "= pv->index_ncdm_\[ncdm_id_\] = pt_"
```

The remaining hits should only be the dual-write lines from Task 2.

- [ ] **Step 4: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run the regression-diff helper. Expect bit-identical output.

- [ ] **Step 5: Commit**

```bash
git add species/
git commit -m "ncdm: species read perturbation slots from own storage instead of shared map"
```

---

## Task 4: Switch readers in `perturbations_module.cpp` to per-species accessors

**Files:**
- Modify: `source/perturbations_module.cpp`

This is the biggest mechanical task: 25+ sites read `pv->index_ncdm_[n][q]` / `pv->l_max_ncdm[n]` / `pv->q_size_ncdm[n]` where `n = ncdm_sp->ncdm_id()`. Each site already has a species pointer in scope (from a `for sp ... if NCDMSpecies` filter); replace the array access with the species accessor.

- [ ] **Step 1: List all read sites**

```bash
grep -n "ppv->index_ncdm_\|ppw->pv->index_ncdm_\|ppv->l_max_ncdm\|ppw->pv->l_max_ncdm\|ppv->q_size_ncdm\|ppw->pv->q_size_ncdm" source/perturbations_module.cpp
```

Expect ~30 hits. Each is in a loop body that already has a species pointer named `ncdm_sp` or `sp` or `nsp` (or similar).

- [ ] **Step 2: Apply substitution at each site**

For each site, identify the species pointer in scope and replace:

| Before | After |
|---|---|
| `ppv->index_ncdm_[n][index_q]` (where `n = ncdm_sp->ncdm_id()`) | `ncdm_base_sp->pt_index_at(index_q)` (cast pointer to `NCDMBaseSpecies*` once at top of loop body if not already) |
| `ppw->pv->index_ncdm_[n][index_q]` | same |
| `ppv->l_max_ncdm[n]` | `ncdm_base_sp->pt_l_max()` |
| `ppv->q_size_ncdm[n]` | `ncdm_base_sp->pt_q_size()` |

For sites that work with a "previous pv" during approximation switches (lines ~3870–4620), the species's `pt_*_previous_` fields are not yet populated (Task 5 adds the `StashPerturbationLayout` hook). For now, leave those switch-copy blocks reading `ppw->pv->index_ncdm_[n][q]` (the old pv) and `ppv->index_ncdm_[n][q]` (the new pv) untouched — they will be replaced wholesale in Task 6.

Concrete example for the source-accumulation block (line ~6403):

```cpp
// Before:
for (auto& [name, sp] : all_species_) {
  auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
  if (!ncdm_sp) continue;
  const int n = ncdm_sp->ncdm_id();
  for (int index_q = 0; index_q < ppw->pv->q_size_ncdm[n]; index_q++) {
    const int idx = ppw->pv->index_ncdm_[n][index_q];
    /* ... use idx ... */
  }
}

// After:
for (auto& [name, sp] : all_species_) {
  auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
  if (!ncdm_sp) continue;
  for (int index_q = 0; index_q < ncdm_sp->pt_q_size(); index_q++) {
    const int idx = ncdm_sp->pt_index_at(index_q);
    /* ... use idx ... */
  }
}
```

- [ ] **Step 3: Verify post-conditions**

After this task, the only remaining reads of `pv->index_ncdm_/l_max_ncdm/q_size_ncdm` in `perturbations_module.cpp` should be inside the approximation-switch copy blocks (lines ~3870–4620). Verify:

```bash
grep -n "index_ncdm_\|l_max_ncdm\|q_size_ncdm" source/perturbations_module.cpp \
  | grep -v "ncdm_id\|register\|storage" | wc -l
```

You should see only the switch-copy blocks (~12 hits, all between lines 3850 and 4630).

- [ ] **Step 4: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run regression-diff. Bit-identical expected.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: read NCDM perturbation slots via per-species accessors"
```

---

## Task 5: Add bookkeeping hooks on `BaseSpecies`

**Files:**
- Modify: `species/base_species.h`
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp`
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`

Add three hooks to `BaseSpecies` with no-op defaults; override on NCDM-family species. No module-side change yet — Task 6 wires them in.

- [ ] **Step 1: Add hooks to `BaseSpecies`**

Edit `species/base_species.h`. After `RegisterTensorPerturbationIndices` (~line 208), add:

```cpp
  // ── Perturbation switch bookkeeping ───────────────────────────────────────
  // Called during perturb_vector_init at every approximation switch.
  // Order is: StashPerturbationLayout() on every species, then a new
  // perturb_vector is built and RegisterPerturbationIndices() runs on every
  // species, then CopyPerturbationsAcrossSwitch() reads the old y at the
  // stashed (previous) layout and writes the new y at the current layout,
  // then MarkUsedInSources() flags slots whose values must NOT propagate to
  // source functions (e.g., NCDM multipoles l > 2).
  // Defaults: no-op (most species don't need migration logic).

  virtual void StashPerturbationLayout() {}
  virtual void CopyPerturbationsAcrossSwitch(const double* /*old_y*/,
                                              double* /*new_y*/) const {}
  virtual void MarkUsedInSources(int* /*used_in_sources*/) const {}
```

- [ ] **Step 2: Implement `StashPerturbationLayout` on `NCDMBaseSpecies`**

Edit `species/ncdm_base_species.h`. In the public section, add the override:

```cpp
  void StashPerturbationLayout() override {
    pt_l_max_previous_       = pt_l_max_;
    pt_q_size_previous_      = pt_q_size_;
    pt_index_per_q_previous_ = pt_index_per_q_;
  }
```

Since this only manipulates members on `NCDMBaseSpecies`, both `NCDMSpecies` and `DNCDMSpecies` inherit it for free.

- [ ] **Step 3: Implement `CopyPerturbationsAcrossSwitch` on `NCDMSpecies`**

The migration logic for ncdmfa-on / ncdmfa-off is currently spread across lines ~4355–4470 of `perturbations_module.cpp`. We are about to lift it wholesale into the species. For this task, write the species-side version so it produces identical results.

Edit `species/ncdm_species.h` to declare:

```cpp
  void CopyPerturbationsAcrossSwitch(const double* old_y, double* new_y) const override;
  void MarkUsedInSources(int* used_in_sources) const override;
```

Edit `species/ncdm_species.cpp`. After `RegisterPerturbationIndices`, add:

```cpp
void NCDMSpecies::CopyPerturbationsAcrossSwitch(const double* old_y, double* new_y) const {
  if (pt_q_size_previous_ < 0 || pt_q_size_ < 0) return;

  const bool fa_was_on = (pt_q_size_previous_ == 1 && pt_l_max_previous_ == 2);
  const bool fa_is_on  = (pt_q_size_ == 1 && pt_l_max_ == 2);

  if (fa_was_on == fa_is_on) {
    // No transition (e.g., this species's slots were just re-registered at
    // the same layout). Copy slot-by-slot.
    for (int iq = 0; iq < pt_q_size_; ++iq) {
      const int new_base = pt_index_per_q_[iq];
      const int old_base = pt_index_per_q_previous_[iq];
      const int lmax     = pt_l_max_;
      for (int l = 0; l <= lmax; ++l)
        new_y[new_base + l] = old_y[old_base + l];
    }
    return;
  }

  if (!fa_was_on && fa_is_on) {
    // Moving INTO fluid approximation: collapse all momentum bins to a single
    // (delta, theta, shear) triple in the new vector. Use the same formulas
    // perturbations_module.cpp:4355 used.
    // PORT: lines 4355–4467 of perturbations_module.cpp verbatim, with the
    // substitutions listed below this code block.
    return;
  }

  if (fa_was_on && !fa_is_on) {
    // Moving OUT of fluid approximation: rebuild the full hierarchy from the
    // collapsed (delta, theta, shear) triple via RescaledPerturbations.
    // PORT: lines 4469–4620 of perturbations_module.cpp verbatim, with the
    // substitutions listed below this code block.
    return;
  }
}

void NCDMSpecies::MarkUsedInSources(int* used_in_sources) const {
  if (pt_q_size_ < 0) return;
  // Multipoles l > 2 do NOT enter source functions.
  for (int iq = 0; iq < pt_q_size_; ++iq) {
    const int base = pt_index_per_q_[iq];
    for (int l = 0; l <= pt_l_max_; ++l)
      if (l > 2) used_in_sources[base + l] = _FALSE_;
  }
}
```

The `PORT:` markers indicate where you need to **read the corresponding module-side block in `perturbations_module.cpp`** (lines 4355–4467 for the ncdmfa-on transition, 4469–4620 for ncdmfa-off) and **copy the arithmetic into the species method verbatim**, with these substitutions:

- `ppv->q_size_ncdm[n]` → `pt_q_size_`
- `ppv->index_ncdm_[n][iq]` → `pt_index_per_q_[iq]`
- `ppw->pv->q_size_ncdm[n]` → `pt_q_size_previous_`
- `ppw->pv->index_ncdm_[n][iq]` → `pt_index_per_q_previous_[iq]`
- `ppv->y[...]` → `new_y[...]`
- `ppw->pv->y[...]` → `old_y[...]`
- `ncdm_sp->q()`, `factor()`, `M()`, `w()`, `dlnf0_dlnq()` → same names without the `ncdm_sp->` prefix (they are members of `NCDMBaseSpecies`)
- For DNCDM-specific branches in the module code, do NOT include them in `NCDMSpecies::CopyPerturbationsAcrossSwitch` — they are handled by the override on `DNCDMSpecies` (next step).

- [ ] **Step 4: Implement `CopyPerturbationsAcrossSwitch` and `MarkUsedInSources` on `DNCDMSpecies`**

Same pattern as `NCDMSpecies` but uses `dq()` (DNCDM-specific quadrature differential) and reads `dlnf0_dlnq` from `pvecback[bg_dlnfdlnq_index_ + iq]` instead of `dlnf0_dlnq_[iq]`. Port the corresponding DNCDM branches from `perturbations_module.cpp` lines ~4355–4620.

Edit `species/dncdm_species.h` and `species/dncdm_species.cpp` accordingly.

`MarkUsedInSources` for `DNCDMSpecies` is identical to `NCDMSpecies` — consider declaring it once on `NCDMBaseSpecies`. To do that, move the implementation from `NCDMSpecies` to `NCDMBaseSpecies` (add to header as inline, since `NCDMBaseSpecies` is currently header-only for inline methods). Then `DNCDMSpecies` inherits it.

- [ ] **Step 5: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

No behavior change yet (hooks added but module hasn't called them). Bit-identical regression expected.

- [ ] **Step 6: Commit**

```bash
git add species/base_species.h species/ncdm_base_species.h species/ncdm_species.{h,cpp} species/dncdm_species.{h,cpp}
git commit -m "species: add Stash/Copy/MarkUsedInSources hooks for approximation-switch bookkeeping"
```

---

## Task 6: Replace NCDM switch-copy blocks in `perturb_vector_init` with hook calls

**Files:**
- Modify: `source/perturbations_module.cpp`

The function `perturb_vector_init` (around line 3220–4640) constructs a new `perturb_vector` at every approximation switch and copies y values from the old one. The NCDM-shaped blocks (lines ~3870, 3960, 4056, 4141, 4246, 4355, 4470, 4619) all do the same thing: "for each NCDM-family species, copy its slots". After this task they go away in favour of `for sp ... sp->CopyPerturbationsAcrossSwitch(...)`.

- [ ] **Step 1: Find the construction site**

In `perturb_vector_init`, find the block that allocates the new `perturb_vector* ppv`, calls all `RegisterPerturbationIndices`, allocates `ppv->y_storage`, and currently does the copy loops. The structure is:

```cpp
struct perturb_vector* ppv = new perturb_vector();
// ... index_pt = 0 ...
// ... call RegisterPerturbationIndices on each species ...
// ... allocate ppv->y_storage ...
// ... NCDM/UR/IDR/IDM_DRMD copy blocks at lines 3870–4620 ...
```

- [ ] **Step 2: Insert `StashPerturbationLayout` BEFORE the new `perturb_vector` construction**

Find the line that allocates `ppv = new perturb_vector()` and immediately above it add:

```cpp
// Snapshot every species's current pt layout so that
// CopyPerturbationsAcrossSwitch can read old y at the previous indices.
for (auto& [_, sp] : all_species_) sp->StashPerturbationLayout();
```

- [ ] **Step 3: Insert `CopyPerturbationsAcrossSwitch` and `MarkUsedInSources` AFTER the y_storage allocation**

After `ppv->y_storage.resize(...); ppv->y = ppv->y_storage.data();`, add:

```cpp
for (auto& [_, sp] : all_species_) sp->CopyPerturbationsAcrossSwitch(ppw->pv->y, ppv->y);
for (auto& [_, sp] : all_species_) sp->MarkUsedInSources(ppv->used_in_sources);
```

- [ ] **Step 4: Delete the NCDM-shaped switch-copy blocks**

For each of the 8 sites listed (the ones containing `for (int index_q = 0; index_q < ppv->q_size_ncdm[n]; ...)` followed by `ppv->y[ppv->index_ncdm_[n][...]]`), delete the entire NCDM block. Leave the UR / IDR / IDM_DRMD / DR blocks untouched — those are out-of-scope per the spec.

For the `used_in_sources` masking block (~line 3580), delete the NCDM-specific block — it is now handled by `MarkUsedInSources`.

After deletion, search for any remaining NCDM-specific code in `perturb_vector_init`:

```bash
grep -n "NCDM\|ncdm_id\|ncdm_sp\|index_ncdm_\[\|q_size_ncdm\[\|l_max_ncdm\[" source/perturbations_module.cpp \
  | awk -F: '$2 >= 3220 && $2 <= 4640'
```

Expected: 0 NCDM-specific lines in `perturb_vector_init`.

- [ ] **Step 5: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run regression-diff. **This is the highest-risk task in the plan** — if the species-side `CopyPerturbationsAcrossSwitch` does not exactly reproduce the module-side block from Task 5, regression will fail. Iterate until bit-identical.

- [ ] **Step 6: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturbations: NCDM approximation-switch copy via species hooks"
```

---

## Task 7: Move tensor-mode NCDM emission into species hooks

**Files:**
- Modify: `species/base_species.h` (only if a new column section enum value is needed)
- Modify: `species/ncdm_species.{h,cpp}`
- Modify: `species/dncdm_species.{h,cpp}`
- Modify: `source/perturbations_module.cpp`

Tensor-mode handling for NCDM today lives in three places in `perturbations_module.cpp`:

1. Tensor source title emission (line ~2540) — emits `delta_ncdm[n]/theta_ncdm[n]/shear_ncdm[n]` titles.
2. Tensor `PerturbDerivs` (line ~7455) — per-momentum-bin derivative loop.
3. Tensor source accumulation (line ~7055) — sums per-q contributions to `delta/theta/shear_ncdm`.

After this task, all three live in NCDM-family overrides of existing `BaseSpecies` virtuals.

- [ ] **Step 1: Override `RegisterTensorPerturbationIndices` on `NCDMSpecies`**

`BaseSpecies::RegisterTensorPerturbationIndices` is already a no-op virtual (base_species.h:205). Override on `NCDMSpecies` to allocate the tensor slots the species needs. Look at how the module currently allocates tensor NCDM slots and lift the logic into the species. Store the resulting indices as `pt_tensor_*_` members.

- [ ] **Step 2: Override `PerturbTensorDerivs` on `NCDMSpecies`**

`BaseSpecies::PerturbTensorDerivs` is already a no-op virtual (base_species.h:227). Override on `NCDMSpecies` and port the loop from `perturbations_module.cpp:~7455` using the species's own tensor slot indices.

- [ ] **Step 3: Override tensor source emission**

For tensor titles + tensor source accumulation, add the emission to the species's `WriteOutputColumns` / `FillSources` overrides (or, if the existing virtuals do not differentiate scalar vs tensor, add a small `is_tensor` mode flag — same pattern as `TransferColumnSection`). The naming uses the instance name (per Task 17), but for now keep `delta_ncdm[<id>]` style — Task 17 will rewrite the title strings.

- [ ] **Step 4: Delete the corresponding module-side blocks**

Delete the tensor NCDM blocks at lines ~2540, ~7055, ~7455 from `perturbations_module.cpp`. Replace each with a generic dispatch:

```cpp
for (auto& [_, sp] : all_species_) sp->PerturbTensorDerivs(tau, y, dy, ppaw);
// (etc. for the other sites)
```

- [ ] **Step 5: Same for `DNCDMSpecies`** (if DNCDM has tensor-mode equations — check the module-side block; if it iterates `dynamic_cast<DNCDMSpecies>`, port it; if it only iterates `NCDMSpecies`, skip).

- [ ] **Step 6: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run regression-diff. The tensor scenarios in `ncdm_multi_unsorted.ini` (which has `modes = s,t`) will exercise the tensor paths.

- [ ] **Step 7: Commit**

```bash
git add species/base_species.h species/ncdm_species.{h,cpp} species/dncdm_species.{h,cpp} source/perturbations_module.cpp
git commit -m "perturbations: NCDM tensor-mode emission via species hooks"
```

---

## Task 8: Move NCDM source accumulation block into `FillSources`

**Files:**
- Modify: `species/ncdm_species.{h,cpp}`
- Modify: `species/dncdm_species.{h,cpp}`
- Modify: `source/perturbations_module.cpp`

`BaseSpecies::FillSources` is already the dispatch hook for source-table writes. Today the NCDM block in `perturb_sources_member` (around line 6390) reads NCDM contributions outside this hook because the legacy code pre-dates the unified dispatch. Move it in.

- [ ] **Step 1: Inspect the current NCDM source block**

Find the block in `perturbations_module.cpp` around line 6390 that iterates NCDM-family species and writes to `p_mod->index_tp_delta_ncdm1_ + n` / `index_tp_theta_ncdm1_ + n` / `index_tp_shear_ncdm1_ + n` source slots. Capture the formulas verbatim.

- [ ] **Step 2: Port into `NCDMSpecies::FillSources`**

The species already has its own `FillSources` (look at how `CDMSpecies::FillSources` writes to `index_tp_delta_cdm_`). Add NCDM equivalents that write to the per-species transfer indices. Note: `index_tp_delta_ncdm1_ + n` is currently a flat-block addressing — Task 9 introduces per-species `tp_delta_index_` members. For this task, use the legacy flat-block addressing (`p_mod->index_tp_delta_ncdm1_ + n` where `n` is now obtained somehow — see next step).

Actually: since Task 9 (per-species transfer indices) cannot land before this task without losing the addressability altogether, this task and Task 9 merge. **Update plan**: combine into a single task. See **Task 9 (revised)** below for the merged version. SKIP this task and proceed to Task 9.

- [ ] **Step 3: Skip — merged into Task 9.**

---

## Task 9: Per-species transfer indices + source accumulation move

**Files:**
- Modify: `species/base_species.h`
- Modify: `species/ncdm_species.{h,cpp}`
- Modify: `species/dncdm_species.{h,cpp}`
- Modify: `source/perturbations_module.cpp`
- Modify: `source/perturbations_module.h`

This task delivers what spec §5 calls "Transfer functions": each NCDM-family species owns its own `tp_delta_index_` / `tp_theta_index_` / `tp_shear_index_` slots, the flat-block `class_define_index(index_tp_delta_ncdm1_, ..., pba->N_ncdm)` is deleted, and `FillSources` writes to the per-species slot instead of `index_tp_delta_ncdm1_ + n`.

- [ ] **Step 1: Add the hook to `BaseSpecies`**

Edit `species/base_species.h`. After `WriteOutputColumns` (line ~272), add:

```cpp
  /**
   * Claim slots in the transfer-function index space (index_tp_*) used by
   * source/transfer accumulation. Default: no-op (most species have a fixed
   * number of transfer slots already declared in the module).
   * @param index_type  in/out counter; advance by however many slots claimed.
   */
  virtual void RegisterTransferIndices(int& /*index_type*/,
                                       const perturbs* /*ppt*/) {}
```

- [ ] **Step 2: Override on `NCDMSpecies` and `DNCDMSpecies`**

Add to each species class:

```cpp
  void RegisterTransferIndices(int& index_type, const perturbs* ppt) override;

 private:
  int tp_delta_index_ = -1;
  int tp_theta_index_ = -1;
  int tp_shear_index_ = -1;  // tensor mode if applicable
```

Implementation in `species/ncdm_species.cpp`:

```cpp
void NCDMSpecies::RegisterTransferIndices(int& index_type, const perturbs* ppt) {
  if (ppt->has_source_delta_ncdm == _TRUE_) {
    tp_delta_index_ = index_type++;
  }
  if (ppt->has_source_theta_ncdm == _TRUE_) {
    tp_theta_index_ = index_type++;
  }
  // (tensor shear, if any)
}
```

Same for `DNCDMSpecies`.

- [ ] **Step 3: Replace flat-block allocation in `perturbations_module.cpp`**

Find the lines (1089, 1105) that do `class_define_index(index_tp_delta_ncdm1_, has_source_delta_ncdm_, index_type, pba->N_ncdm)` and delete them. Replace with a per-species iteration:

```cpp
for (auto& [_, sp] : all_species_) sp->RegisterTransferIndices(index_type, ppt);
```

The shared `index_tp_delta_ncdm1_` / `index_tp_theta_ncdm1_` member fields on `PerturbationsModule` can be deleted from `source/perturbations_module.h` once nothing reads them; verify with a grep before deleting.

- [ ] **Step 4: Move source accumulation into `FillSources`**

In `species/ncdm_species.cpp::FillSources` (and `species/dncdm_species.cpp::FillSources`), implement what was at `perturbations_module.cpp:~6390`. Use the species's own `tp_delta_index_` / `tp_theta_index_` instead of `index_tp_delta_ncdm1_ + n`. The arithmetic for `delta_ncdm` / `theta_ncdm` / `shear_ncdm` from the per-q sums is the same; just the slot name changes.

- [ ] **Step 5: Delete the module-side NCDM source block**

Delete the block at `perturbations_module.cpp:~6390`. The generic `for (auto& [_, sp] : all_species_) sp->FillSources(...)` already runs each species's `FillSources` — no new dispatch needed.

- [ ] **Step 6: Update consumers of the flat transfer slots**

Search for any code that does `index_tp_delta_ncdm1_ + n`:

```bash
grep -n "index_tp_delta_ncdm\|index_tp_theta_ncdm\|index_tp_shear_ncdm" source/ include/
```

Each hit is a consumer that previously addressed by `+ n`. Replace with iteration over species and use of the per-species `tp_*_index_` accessor. Add public accessors on `NCDMSpecies` / `DNCDMSpecies`:

```cpp
int tp_delta_index() const { return tp_delta_index_; }
int tp_theta_index() const { return tp_theta_index_; }
```

- [ ] **Step 7: Drop the scratch vectors `delta_ncdm.resize(pba->N_ncdm)`** at lines ~6903–6905. Replace with locally-scoped accumulation inside the per-species loop.

- [ ] **Step 8: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run regression-diff. Output column ORDER may shift if species iteration order differs from old construction order — for the unsorted scenario this is expected and the diff helper will flag it. Verify by reading the diff: column values should match column-by-column even if column indices differ.

- [ ] **Step 9: Commit**

```bash
git add species/ source/perturbations_module.{h,cpp}
git commit -m "perturbations: per-species transfer indices and FillSources for NCDM"
```

---

## Task 10: Remove `pba->N_ncdm`

**Files:**
- Modify: `source/background.h`
- Modify: `source/input_module.cpp`
- Modify: `source/output_module.cpp`
- Modify: `source/background_module.cpp`
- Modify: `source/perturbations_module.cpp`

- [ ] **Step 1: List every read site**

```bash
grep -n "pba->N_ncdm\|ba\.N_ncdm\|ba->N_ncdm" source/*.cpp source/*.h species/*.cpp species/*.h
```

Expect ~12 hits.

- [ ] **Step 2: Add a count helper on `SpeciesCollection`**

Edit `species/species_collection.h` (locate the file and inspect existing helpers like `count(name)`). Add:

```cpp
int count_ncdm_family() const {
  int n = 0;
  for (const auto& [_, sp] : *this)
    if (dynamic_cast<NCDMBaseSpecies*>(sp.get())) ++n;
  return n;
}
```

(This is a temporary helper. The clean answer is "no module-side species-type filtering" — but the alternative requires a `BaseSpecies::IsNCDMFamily()` virtual or similar, which carves out NCDM-specific category-membership at the architectural level. We'll keep the typed helper here behind a single function, and the spec's TODO marker covers redesign.)

- [ ] **Step 3: Apply rewrites at each read site**

| Site (today) | Rewrite |
|---|---|
| `if (pba->N_ncdm > 0) { ... }` as a presence guard around code that does `for sp ... NCDMSpecies` | Delete the guard. Iteration is naturally empty. |
| `if (pba->N_ncdm > 0) { ... }` around code that uses pba->* legacy fields | Replace with `if (all_species_.count_ncdm_family() > 0)` until the body is migrated, then delete. |
| `pba->N_ncdm` as array allocation (lines 1089/1105/6903–6905 of perturbations_module.cpp) | Already removed in Task 9. |
| `class_test(pba->N_ncdm > 1, ...)` in output_module.cpp:823 | Replace with `class_test(all_species_.count_ncdm_family() > 1, ...)` |

- [ ] **Step 4: Remove the write site in `input_module.cpp:346`**

```cpp
// Before:
pba->N_ncdm = n_ncdm;
// After: delete the line.
```

- [ ] **Step 5: Delete the field from `background.h:133`**

```cpp
// Delete:
int N_ncdm = 0;  /**< Number of distinguishable ncdm species */
```

- [ ] **Step 6: Build**

```bash
make class -j 2>&1 | tail -10
```

Expected: clean. If any compile errors remain, finish migrating the missed read sites.

- [ ] **Step 7: Smoke + regression**

```bash
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run regression-diff. Bit-identical expected.

- [ ] **Step 8: Commit**

```bash
git add source/background.h source/input_module.cpp source/output_module.cpp \
        source/background_module.cpp source/perturbations_module.cpp species/species_collection.h
git commit -m "background: remove pba->N_ncdm; consumers iterate species directly"
```

---

## Task 11: Remove `pba->Omega0_ncdm_tot` and `pba->N_decay_dr`

**Files:**
- Modify: `source/background.h`
- Modify: `source/input_module.cpp`
- Modify: `source/nonlinear_module.cpp`
- Modify: `source/thermodynamics_module.cpp`
- Modify: `source/perturbations_module.cpp`
- Modify: `source/background_module.cpp`

- [ ] **Step 1: List read sites**

```bash
grep -n "Omega0_ncdm_tot\|N_decay_dr" source/*.cpp source/*.h species/*.cpp species/*.h
```

Expect ~8 hits between the two fields.

- [ ] **Step 2: Replace `pba->Omega0_ncdm_tot` consumers**

Each consumer becomes a sum-over-NCDM-family. Concrete pattern:

```cpp
// Before:
fnu = pba->Omega0_ncdm_tot / Omega0_m;

// After:
double Omega0_ncdm_tot = 0.;
for (auto& [_, sp] : all_species_) {
  // TODO(architecture): redesign once Omega0_b/Omega0_cdm get the same treatment.
  if (auto* nsp = dynamic_cast<NCDMBaseSpecies*>(sp.get()))
    Omega0_ncdm_tot += nsp->GetOmega0();
}
fnu = Omega0_ncdm_tot / Omega0_m;
```

The 5 sites are: `nonlinear_module.cpp:2133`, `nonlinear_module.cpp:2661`, `nonlinear_module.cpp:2666`, `thermodynamics_module.cpp:3150`, `input_module.cpp:3096`. Wherever the consumer is in a tight loop, hoist the `Omega0_ncdm_tot` computation outside the loop.

The `TODO(architecture)` comment is required at every site. This is the only place in the PR that breaks the "no species-type picking" rule, and it must be visible to whoever picks up the follow-up.

- [ ] **Step 3: Replace `pba->N_decay_dr` consumers**

```bash
grep -n "N_decay_dr" source/*.cpp species/*.cpp
```

Each becomes `count_dncdm_dr_family()` or equivalent — likely a similar helper to `count_ncdm_family()` but typed on `DNCDM_DR_Species`. Add to `SpeciesCollection`:

```cpp
int count_dncdm_dr_family() const {
  int n = 0;
  for (const auto& [_, sp] : *this)
    if (dynamic_cast<DNCDM_DR_Species*>(sp.get())) ++n;
  return n;
}
```

The classic site is `background_module.cpp` where `class_store_columntitle` iterates `n < pba->N_decay_dr` — replace `pba->N_decay_dr` with the helper or, better, make each `DNCDM_DR_Species` write its own column titles directly.

- [ ] **Step 4: Remove writes**

In `input_module.cpp:348`:

```cpp
// Before:
pba->Omega0_ncdm_tot = omega0_ncdm_tot;
// After: delete the line.
```

Find and delete the analogous write for `pba->N_decay_dr`.

- [ ] **Step 5: Delete fields from `background.h:134, 135`**

```cpp
// Delete:
double Omega0_ncdm_tot = 0.;
int N_decay_dr = 0;
```

- [ ] **Step 6: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/ncdm_*.ini test/scenarios/dncdm_dr.ini; do
  ./class "$f" 2>&1 | tail -3
done
```

Run regression-diff. Bit-identical expected (the TOTAL Omega is the same, just computed from species instead of cached on pba).

- [ ] **Step 7: Commit**

```bash
git add source/background.h source/input_module.cpp source/nonlinear_module.cpp \
        source/thermodynamics_module.cpp source/perturbations_module.cpp \
        source/background_module.cpp species/species_collection.h
git commit -m "background: remove Omega0_ncdm_tot and N_decay_dr; sum from species"
```

---

## Task 12: Add shooter hooks and `ShootingContext` header

**Files:**
- Create: `source/shooting_context.h`
- Modify: `species/base_species.h`
- Modify: `Makefile`
- Modify: `setup.py`
- Modify: `CLASS.xcodeproj/project.pbxproj` (the user verifies Xcode manually per project memory)

- [ ] **Step 1: Inspect today's shooter API**

Open `source/input_module.cpp` and read `input_try_unknown_parameters` (around line 3493) and `input_get_guess` (around line 3675). Note what each `case` needs from the surrounding context: cosmology pointer (`ba`), pre-built `BackgroundModulePtr`, verbose flag, `pfzw->target_values[idx]`. Bundle these into `ShootingContext`.

- [ ] **Step 2: Create `source/shooting_context.h`**

```cpp
#pragma once

class BackgroundModule;
struct background;
struct precision;

/**
 * Bundle of read-only context passed to species during the shooter iteration.
 * Mirrors what input_try_unknown_parameters / input_get_guess need today.
 */
struct ShootingContext {
  const background* pba;
  const precision*  ppr;
  const BackgroundModule* bgm;  // may be nullptr during initial guess
  int verbose;
  // Per-call data filled by InputModule before dispatching:
  const double* xvalues;   // current trial parameter vector
  const double* target_values;  // user-requested target values
};
```

- [ ] **Step 3: Add shooter hooks to `BaseSpecies`**

Edit `species/base_species.h`. Add `#include "../source/shooting_context.h"` (or put `ShootingContext` in `species/` if that's cleaner — pick whichever matches the existing header layout). After the matter-tally section, add:

```cpp
  // ── Shooter ───────────────────────────────────────────────────────────────
  // Mirrors RegisterIntegrationIndices: each species claims slots and stores
  // its own indices internally.

  /**
   * Claim slots in the global xvalues/target_values arrays. Read pfc/ba to
   * decide which target(s) (if any) this species needs to shoot for, and
   * store the captured target value internally.
   */
  virtual void RegisterShootingIndices(int& /*index_sh*/,
                                       FileContent* /*pfc*/,
                                       const background& /*ba*/) {}

  /**
   * Compute this species' residual contribution: output[sh_*_index_] = ...
   */
  virtual void ComputeShootingResidual(double* /*output*/,
                                       const ShootingContext& /*ctx*/) const {}

  /**
   * Compute this species' initial-guess + Jacobian-diagonal contribution.
   */
  virtual void ComputeShootingGuess(double* /*xguess*/,
                                    double* /*dxdy*/,
                                    const ShootingContext& /*ctx*/) const {}
```

Add a forward declaration of `FileContent` near the top of `base_species.h` if not already present.

- [ ] **Step 4: Register the new header in build systems**

`Makefile`:

```bash
grep -n "shooting_context\|species_build_context" Makefile
```

Find the existing pattern that lists header dependencies. If header-only files are auto-discovered, no Makefile change is needed — verify by building. If they are explicitly listed, add `source/shooting_context.h` to the right list.

`setup.py`:

```bash
grep -n "shooting_context\|species_build_context" setup.py
```

Same logic — add if there is an explicit list.

`CLASS.xcodeproj/project.pbxproj`: per project memory, the user verifies Xcode manually. Note in the commit message that the Xcode project may need `shooting_context.h` added.

- [ ] **Step 5: Build + smoke**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
```

No behavior change expected (hooks added but not yet wired in).

- [ ] **Step 6: Commit**

```bash
git add source/shooting_context.h species/base_species.h Makefile setup.py
git commit -m "shooter: add ShootingContext + Register/Residual/Guess hooks on BaseSpecies"
```

---

## Task 13: Move DNCDM_DR shooting cases into `DNCDM_DR_Species`

**Files:**
- Modify: `species/dncdm_dr_species.{h,cpp}`
- Modify: `source/input_module.cpp`

The current `omega_dncdmdr` / `Omega_dncdmdr` / `Neff_ini_dncdm` / `deg_ncdm_decay_dr` / `omega_ini_dncdm` / `Omega_ini_dncdm` cases live in `input_try_unknown_parameters` (lines ~3616–3666) and `input_get_guess` (lines ~3771–3845). Each iterates `dncdm_id` from 0 to `dncdm_dr_species.size() - 1` and writes to `output[idx + dncdm_id]`. Each instance now writes to its own slot.

- [ ] **Step 1: Add storage to `DNCDM_DR_Species`**

Edit `species/dncdm_dr_species.h`. In the private section:

```cpp
  // Shooter slots (-1 if this target is not active for this instance).
  int sh_omega_dncdmdr_index_ = -1;
  int sh_Omega_dncdmdr_index_ = -1;
  int sh_Neff_ini_dncdm_index_ = -1;
  int sh_deg_ncdm_decay_dr_index_ = -1;
  int sh_omega_ini_dncdm_index_ = -1;
  int sh_Omega_ini_dncdm_index_ = -1;

  // Captured target value for whichever target is active.
  double sh_target_value_ = 0.;
  enum class ShootingTarget { None, omega_dncdmdr, Omega_dncdmdr, Neff_ini_dncdm,
                              deg_ncdm_decay_dr, omega_ini_dncdm, Omega_ini_dncdm };
  ShootingTarget sh_active_target_ = ShootingTarget::None;
```

- [ ] **Step 2: Implement `RegisterShootingIndices`**

```cpp
void DNCDM_DR_Species::RegisterShootingIndices(int& index_sh, FileContent* pfc, const background& ba) {
  // Decide which target this DNCDM instance shoots for, by inspecting input.
  // The same input rules that input_module.cpp uses today must be applied.
  // Concrete: if user supplied "omega_dncdmdr" via dot-syntax for THIS instance,
  // then shoot for omega_dncdmdr; if "Omega_dncdmdr" then that; etc.
  // SpeciesInput already provides per-instance reads.

  SpeciesInput inp(pfc, name());  // uses the composite's instance name
  double v;
  if (inp.read_double("omega_dncdmdr", v)) {
    sh_active_target_ = ShootingTarget::omega_dncdmdr;
    sh_omega_dncdmdr_index_ = index_sh++;
    sh_target_value_ = v;
  } else if (inp.read_double("Omega_dncdmdr", v)) {
    /* ... */
  } else if (...) {
    /* ... etc. */
  }
  // If none, sh_active_target_ stays None and no slot is registered.
}
```

The exact list of target keys + their interpretation rules must match what `input_module.cpp:input_read_parameters_general_targets` does today. Consult the existing target-name → input-string mapping table near line ~570 of `input_module.cpp` (look for `target_namestrings`).

- [ ] **Step 3: Implement `ComputeShootingResidual`**

Lift the body of `case omega_dncdmdr: case Omega_dncdmdr: ...` from `input_try_unknown_parameters` into the species method. The current body's per-instance loop becomes a single computation (the species IS the instance).

```cpp
void DNCDM_DR_Species::ComputeShootingResidual(double* output, const ShootingContext& ctx) const {
  if (sh_active_target_ == ShootingTarget::None) return;

  const double* bg_today = ctx.bgm->background_table_.data() +
                           (ctx.bgm->bt_size_ - 1) * ctx.bgm->bg_size_;
  double rho_dr_today    = dr_.Rho(bg_today);
  double rho_dncdm_today = dncdm_.Rho(bg_today);

  const double H0_2 = ctx.pba->H0 * ctx.pba->H0;

  switch (sh_active_target_) {
    case ShootingTarget::omega_dncdmdr: {
      double Omega_target = sh_target_value_ / ctx.pba->h / ctx.pba->h;
      output[sh_omega_dncdmdr_index_] = (rho_dr_today + rho_dncdm_today) / H0_2 - Omega_target;
      if (ctx.verbose > 0)
        printf(" -> Shooting iteration: Omega_dncdmdr (input) = %g, computed = %g\n",
               Omega_target, (rho_dr_today + rho_dncdm_today) / H0_2);
      break;
    }
    case ShootingTarget::Omega_dncdmdr: { /* ... */ break; }
    /* ... etc. */
    default: break;
  }
}
```

- [ ] **Step 4: Implement `ComputeShootingGuess`**

Same pattern — port the `case omega_dncdmdr: case Omega_dncdmdr:` body from `input_get_guess` (line ~3771). The `for (int dncdm_id = 0; dncdm_id < ...; ++dncdm_id)` collapses to a single computation per species.

- [ ] **Step 5: Delete the legacy cases from `input_module.cpp`**

Delete the `case omega_dncdmdr: case Omega_dncdmdr: case Neff_ini_dncdm: case deg_ncdm_decay_dr: case omega_ini_dncdm: case Omega_ini_dncdm: { ... }` blocks from both `input_try_unknown_parameters` (~lines 3616–3666) and `input_get_guess` (~lines 3771–3845).

But do NOT yet delete the `target_name` / `target_size` infrastructure — Task 15 will refactor the dispatch to call species methods.

- [ ] **Step 6: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
./class test/scenarios/dncdm_dr.ini 2>&1 | tail -3
```

Run regression-diff specifically on `dncdm_dr.dat` outputs. Bit-identical expected.

- [ ] **Step 7: Commit**

```bash
git add species/dncdm_dr_species.{h,cpp} source/input_module.cpp
git commit -m "shooter: DNCDM_DR shooting moved into DNCDM_DR_Species"
```

---

## Task 14: Move DCDM_DR and ScalarField shooting cases

**Files:**
- Modify: `species/dcdm_dr_species.{h,cpp}`
- Modify: `species/scalar_field.{h,cpp}`
- Modify: `source/input_module.cpp`

Same pattern as Task 13 but simpler: each species owns one shooter slot.

- [ ] **Step 1: DCDM_DR — add `sh_*_index_` storage + Register/Residual/Guess overrides for `omega_dcdmdr`, `Omega_dcdmdr`, `omega_ini_dcdm`, `Omega_ini_dcdm`**

Lift the corresponding bodies from `input_module.cpp:~3565–3615` (residual) and `~3704–3769` (guess).

- [ ] **Step 2: ScalarField — same for `Omega_scf`**

Lift from `input_module.cpp:~3728–3746` (guess) and the corresponding residual case.

- [ ] **Step 3: Delete the legacy cases from `input_module.cpp`**

- [ ] **Step 4: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
# Verify SCF and DCDM scenarios:
./class test/dcdm_dr.ini 2>&1 | tail -3 || true   # if exists
./class test/scf.ini 2>&1 | tail -3 || true       # if exists
```

If standalone DCDM/SCF scenarios do not exist, `explanatory.ini` exercises both. Regression-diff against baseline.

- [ ] **Step 5: Commit**

```bash
git add species/dcdm_dr_species.{h,cpp} species/scalar_field.{h,cpp} source/input_module.cpp
git commit -m "shooter: DCDM_DR + ScalarField shooting moved into species"
```

---

## Task 15: Refactor `InputModule` shooter dispatch

**Files:**
- Modify: `source/input_module.cpp`
- Modify: `source/input_module.h`

After Tasks 13 and 14, the only `case` left in the giant `switch(target_name[counter])` blocks is `theta_s`. Refactor the dispatch to iterate species instead of switching on target enum.

- [ ] **Step 1: Replace `input_try_unknown_parameters` body**

The function currently:
1. Calls `input_set_pfzw_target_values_into_ba(...)` to copy target values into background.
2. Computes background.
3. Loops `counter` from 0 to `unknown_parameters_size`, switches on `target_name[counter]`, writes to `output[idx ... idx + target_size[counter]]`.

After:
1. Same.
2. Same.
3. `output[idx_theta_s] = ... theta_s residual ...; for (auto& [_, sp] : all_species_) sp->ComputeShootingResidual(output, ctx);`

```cpp
ShootingContext ctx{
  /*pba=*/&ba, /*ppr=*/&pr, /*bgm=*/cosmology.GetBackgroundModule().get(),
  /*verbose=*/input_verbose,
  /*xvalues=*/unknown_values,
  /*target_values=*/pfzw->target_values.data(),
};

// theta_s: handled inline (only non-species cosmology target).
if (theta_s_index_ >= 0) {
  output[theta_s_index_] = 100. * thm->rs_rec_ / thm->ra_rec_ - pfzw->target_values[theta_s_index_];
}

for (auto& [_, sp] : input_module->all_species_)
  sp->ComputeShootingResidual(output, ctx);
```

- [ ] **Step 2: Replace `input_get_guess` body**

Same structure: keep theta_s inline, dispatch the rest to species.

- [ ] **Step 3: Delete `target_name` / `target_size` plumbing**

Once no `case` in `input_try_unknown_parameters` reads `pfzw->target_name[counter]` other than `theta_s`, the per-counter array of enum target names becomes redundant. Replace `pfzw->target_name`, `pfzw->target_sizes`, `pfzw->target_values`, `pfzw->unknown_parameter_names` with whatever minimal state is still needed for theta_s + species iteration. Likely:

- `pfzw->theta_s_target_value` — single double, only set if theta_s is active
- `pfzw->theta_s_index` — int slot
- (the species own their own captured target values via `RegisterShootingIndices`)

- [ ] **Step 4: Update `RegisterShootingIndices` driver**

The function that today builds `pfzw->target_name[counter]` etc. is `input_read_parameters_general_targets` (or similar) around line ~557. Replace with:

```cpp
int index_sh = 0;
RegisterCosmologicalShootingIndices(index_sh);  // theta_s
for (auto& [_, sp] : all_species_) sp->RegisterShootingIndices(index_sh, &file_content_, background_);
unknown_parameters_size = index_sh;
```

- [ ] **Step 5: Build + smoke + regression — full pytest**

```bash
make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/*.ini; do ./class "$f" 2>&1 | tail -3; done

cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

This is the first task that touches the shooter dispatch — the full pytest suite is the safety net. If any scenario fails, find the missing or mis-ported target and fix it.

- [ ] **Step 6: Commit**

```bash
git add source/input_module.{h,cpp}
git commit -m "shooter: dispatch via per-species hooks; theta_s remains inline"
```

---

## Task 16: Rename DNCDM_DR composite + DR child

**Files:**
- Modify: `species/dncdm_dr_species.cpp`
- Modify: `species/dncdm_decay_radiation_species.h`
- Modify: any code that does `all_species_.at("DNCDM_DR_<id>")` or `all_species_.at("DNCDM_DecayRadiation_<id>")`

- [ ] **Step 1: Find all literal references to the legacy composite name**

```bash
grep -rn "DNCDM_DR_\|DNCDM_DecayRadiation_" source/ species/ include/
```

Expect ~5 hits: the `CompositeSpecies(...)` constructor in `dncdm_dr_species.cpp:24`, the `BaseSpecies("DNCDM_DecayRadiation_" + ...)` in `dncdm_decay_radiation_species.h:16`, and any module-side lookups.

- [ ] **Step 2: Rewrite the constructors**

`species/dncdm_dr_species.cpp:24`:

```cpp
// Before:
: CompositeSpecies("DNCDM_DR_" + std::to_string(dncdm_arg->ncdm_id()), ...)
// After:
: CompositeSpecies(dncdm_arg->name(), ...)   // composite takes the DNCDM child's instance name
```

`species/dncdm_decay_radiation_species.h:15-17`:

```cpp
// Before:
DNCDM_DecayRadiationSpecies(int ncdm_id, const background* pba, const BackgroundModule* bgm)
    : BaseSpecies("DNCDM_DecayRadiation_" + std::to_string(ncdm_id), EnergyType::Radiation),
      ncdm_id_(ncdm_id), pba_(pba), bgm_(bgm) {}
// After:
DNCDM_DecayRadiationSpecies(const std::string& parent_name, const background* pba, const BackgroundModule* bgm)
    : BaseSpecies(parent_name + "_DR", EnergyType::Radiation),
      pba_(pba), bgm_(bgm) {}
```

The `ncdm_id_` member of the DR child is no longer needed — drop it. Update the call site in `dncdm_dr_species.cpp:27`:

```cpp
// Before:
auto dr_sp = std::make_unique<DNCDM_DecayRadiationSpecies>(ncdm_id_, pba, bgm);
// After:
auto dr_sp = std::make_unique<DNCDM_DecayRadiationSpecies>(dncdm_arg->name(), pba, bgm);
```

- [ ] **Step 3: Update SpeciesCollection insertion key**

In whatever `CreateAll` factory inserts the composite into `all_species_`, the key string used today is likely `"DNCDM_DR_" + std::to_string(...)` or similar. Find it:

```bash
grep -n "DNCDM_DR" species/dncdm_dr_species.cpp source/input_module.cpp
```

Change the insertion key to the composite's `BaseSpecies::name()` (which is now the DNCDM instance name).

- [ ] **Step 4: Update any module-side lookups**

```bash
grep -rn "all_species_.at(\"DNCDM_DR\|all_species_.count(\"DNCDM_DR\|all_species_.find(\"DNCDM_DR" source/ species/
```

Each hit needs to either iterate species by type (since "all DNCDM_DR composites" is no longer one literal key) or be replaced by some other mechanism. Most likely candidates: a `for sp ... if dynamic_cast<DNCDM_DR_Species>` loop.

- [ ] **Step 5: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
./class test/scenarios/dncdm_dr.ini 2>&1 | tail -3
./class test/scenarios/ncdm_dncdm_idmdr_combined.ini 2>&1 | tail -3
```

Output column ORDER may shift again (the composite is now keyed by instance name). The regression-diff helper will flag new column positions; verify column-by-column values still match.

- [ ] **Step 6: Commit**

```bash
git add species/dncdm_dr_species.cpp species/dncdm_decay_radiation_species.h source/
git commit -m "dncdm: composite + DR child use the DNCDM instance name (no ncdm_id)"
```

---

## Task 17: Rewrite NCDM-family output column titles to instance-name based

**Files:**
- Modify: `species/ncdm_species.cpp`
- Modify: `species/dncdm_species.cpp`
- Modify: `species/dncdm_decay_radiation_species.{h,cpp}` (may need an .cpp file)
- Modify: `source/background_module.cpp` (verbose budget prints)
- Modify: `source/perturbations_module.cpp` (tensor titles, if not already moved by Task 7)

- [ ] **Step 1: Background-output titles in `species/ncdm_species.cpp:232–246`**

```cpp
// Before:
snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
// After:
snprintf(tmp, 80, "(.)number_%s", name().c_str());
```

Apply the same pattern to `(.)rho_ncdm[%d]`, `(.)p_ncdm[%d]`, `(.)pseudo_p_ncdm[%d]`, etc. Bump buffer size to 80 to accommodate longer instance names.

Repeat for the `WriteBackgroundData` block (lines 242–246).

- [ ] **Step 2: Same in `species/dncdm_species.cpp:376–408`**

```cpp
// Before:
snprintf(tmp, 40, "(.)number_ncdm[%d]", ncdm_id_);
snprintf(tmp, 40, "lnf_dncdm[%d][%d]", ncdm_id_, i);
// After:
snprintf(tmp, 80, "(.)number_%s", name().c_str());
snprintf(tmp, 80, "lnf_%s[%d]", name().c_str(), i);
```

Note: the inner `[%d]` for momentum-bin index `i` stays — only the species-id `[%d]` is removed.

- [ ] **Step 3: DR child column titles**

In `species/dncdm_decay_radiation_species.h` (or its .cpp if one exists), find any `(.)rho_dr_species` or similar emission. Per the spec, rename to `(.)rho_<instance>_DR` — but `name()` is already `<instance>_DR` after Task 16, so:

```cpp
// Before:
snprintf(tmp, 40, "(.)rho_dr_species");
// After:
snprintf(tmp, 80, "(.)rho_%s", name().c_str());
```

- [ ] **Step 4: Perturbation-output titles (scalar `delta_ncdm[%d]`, etc.)**

These live in each species's `WriteOutputColumns` and `PrintVariables`. Apply the same `[%d]` → `_<instance>` substitution. After this task, output strings should look like `delta_nu1`, `theta_nu1`, `shear_nu1` for each instance.

- [ ] **Step 5: Tensor titles in `perturbations_module.cpp:~2540`**

If Task 7 already moved tensor title emission into species, the rename happens in the species; verify. If not, do the rename in module code.

- [ ] **Step 6: Verbose budget prints in `background_module.cpp:~559, ~2115ff`**

Find prints of the form `printf("-> N_eff = %g for ncdm species %d\n", ..., ncdm_id_);` and change `%d`→`%s`, passing `name().c_str()`.

- [ ] **Step 7: Final grep audit**

```bash
grep -rn "ncdm\[%d\]\|ncdm species %d\|ncdm_id_\|ncdm_id()" species/ source/
```

Remaining hits should only be in code that explicitly stores/reads the integer (Task 18 will delete those). No printf/snprintf format string should still embed `ncdm_id`.

- [ ] **Step 8: Build + smoke + regression**

```bash
make class -j 2>&1 | tail -5
./class explanatory.ini 2>&1 | tail -3
for f in test/scenarios/*.ini; do ./class "$f" 2>&1 | tail -3; done
```

**Output COLUMN NAMES change in this task.** The numeric values of physical quantities are unchanged but column headers are different. The regression-diff helper will report shape/column mismatches — that is expected. Verify by manually inspecting one output file:

```bash
head -3 output/ncdm_single_background.dat   # check headers
diff <(head -3 /tmp/ncdm_baseline/ncdm_single_background.dat) <(head -3 output/ncdm_single_background.dat)
# Expect: header lines differ (rho_ncdm[0] → rho_nu1); data lines unchanged
```

Run pytest:

```bash
cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

If pytest scenarios assert specific column names, update them as part of this task. Search:

```bash
grep -rn "rho_ncdm\[\|delta_ncdm\[\|theta_ncdm\[\|shear_ncdm\[" python/
```

- [ ] **Step 9: Commit**

```bash
git add species/ source/background_module.cpp source/perturbations_module.cpp python/
git commit -m "output: NCDM-family columns use instance name instead of ncdm_id integer"
```

---

## Task 18: Final removal — delete `ncdm_id` field and dead pv/Module fields

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_species.{h,cpp}`
- Modify: `species/ncdm_interacting_species.cpp`
- Modify: `species/dncdm_species.{h,cpp}`
- Modify: `species/dncdm_dr_species.{h,cpp}`
- Modify: `species/species_build_context.h`
- Modify: `source/perturbations.h`
- Modify: `source/perturbations_module.{h,cpp}`
- Modify: `source/input_module.cpp`

After all consumers have migrated to per-species storage, the legacy fields can be deleted.

- [ ] **Step 1: Delete `ncdm_id_` and accessors from species classes**

`species/ncdm_species.h`:

```cpp
// Delete:
int ncdm_id() const { return ncdm_id_; }
void SetNcdmId(int id) override { ncdm_id_ = id; }
int ncdm_id_ = -1;
```

`species/dncdm_species.h`: same.

`species/ncdm_base_species.h`: delete the pure virtual `virtual void SetNcdmId(int id) = 0;` declaration (line 115).

- [ ] **Step 2: Delete `SetNcdmId` calls in CreateAll factories**

`species/ncdm_species.cpp:193`: `sp->SetNcdmId((*ctx.ncdm_id_next)++);` — delete.
`species/ncdm_interacting_species.cpp:80`: same — delete.
`species/dncdm_species.cpp:191`: same — delete.

- [ ] **Step 3: Delete `ncdm_id_next` from `SpeciesBuildContext`**

`species/species_build_context.h:25`:

```cpp
// Delete:
int* ncdm_id_next;  // non-null; mutable counter, advanced by NCDM-family CreateAll
```

`source/input_module.cpp:215, 222`: delete the local counter and the field initializer.

- [ ] **Step 4: Delete legacy fields from `perturb_vector` struct**

`source/perturbations.h:307–314`:

```cpp
// Delete all of:
int index_pt_psi0_ncdm1;
int N_ncdm;
int* l_max_ncdm;
int* q_size_ncdm;
std::map<int, std::vector<int>> index_ncdm_;
std::vector<int> l_max_ncdm_storage;
std::vector<int> q_size_ncdm_storage;
```

Remove the corresponding `ppv->index_pt_psi0_ncdm1 = index_pt;` write at `perturbations_module.cpp:3450` and `species/ncdm_species.cpp:300`.

Remove the `ppv->N_ncdm = pba->N_ncdm; ppv->l_max_ncdm_storage.resize(...); ...` block at `perturbations_module.cpp:3357–3360, 3452–3455`.

Also remove the dual-write block introduced in Task 2 from `species/ncdm_species.cpp:RegisterPerturbationIndices` (the lines that did `pv->l_max_ncdm[ncdm_id_] = pt_l_max_;` etc.) and the analogous one in `species/dncdm_species.cpp`.

- [ ] **Step 5: Delete `ncdm_species_sorted_` from `PerturbationsModule`**

`source/perturbations_module.h`: find the `ncdm_species_sorted_` member declaration and delete it.

`source/perturbations_module.cpp`: find where it is populated (around lines 640–660 — the loop that pushes NCDMSpecies/DNCDMSpecies pointers and sorts by `ncdm_id`) and delete the entire block.

Replace any remaining iteration over `ncdm_species_sorted_` (e.g. lines 5220–5235, the IDR-related block) with iteration over `all_species_` filtered by `dynamic_cast<NCDMBaseSpecies>` — or, ideally, by virtual dispatch with no filter. If the body of the loop is a physics computation, prefer moving it into the species (smaller scope but adds another hook).

- [ ] **Step 6: Change `RescaledNCDMPerturbations` API**

`source/perturbations_module.h`: change signature.

```cpp
// Before:
std::tuple<double, double, double> RescaledNCDMPerturbations(int n_ncdm, double a, double k, perturb_workspace* ppw);
// After:
std::tuple<double, double, double> RescaledNCDMPerturbations(BaseSpecies* sp, double a, double k, perturb_workspace* ppw);
```

`source/perturbations_module.cpp:7965`: implementation.

```cpp
std::tuple<double, double, double> PerturbationsModule::RescaledNCDMPerturbations(
    BaseSpecies* sp, double a, double k, perturb_workspace* ppw) {
  if (auto* composite = dynamic_cast<DNCDM_DR_Species*>(sp))
    return composite->dncdm().RescaledPerturbations(a, k, ppw);
  throw std::runtime_error("RescaledNCDMPerturbations: species is not a DNCDM_DR composite");
}
```

Find the single call site of `RescaledNCDMPerturbations` and update it:

```bash
grep -n "RescaledNCDMPerturbations" source/ species/
```

- [ ] **Step 7: Delete `dncdm_decay_radiation_species.h:84` (`int ncdm_id_;`)**

`species/dncdm_dr_species.h:53`: same.

Delete every remaining read of `ncdm_id_` in `species/`.

- [ ] **Step 8: Final grep audit — must be EMPTY**

```bash
grep -rn "ncdm_id" source/ species/ include/ 2>/dev/null
```

Expected output: zero hits. If any remain, fix them in this task.

```bash
grep -rn "index_ncdm_\|l_max_ncdm\|q_size_ncdm\|N_ncdm\|index_pt_psi0_ncdm1" source/ species/ include/ 2>/dev/null
```

Expected: zero hits except for `pba->has_ncdm` / `pba->has_ncdm_decay_dr` (out of scope; has-guards) and per-species local fields like `pt_q_size_` (different identifier).

```bash
grep -rn "ncdm_species_sorted_\|SetNcdmId" source/ species/ 2>/dev/null
```

Expected: zero hits.

- [ ] **Step 9: Build + smoke**

```bash
make clean && make class -j 2>&1 | tail -10
./class explanatory.ini 2>&1 | tail -3
```

- [ ] **Step 10: Full regression + pytest**

```bash
for f in test/scenarios/*.ini; do
  echo "=== $f ==="
  ./class "$f" 2>&1 | tail -3
done

cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

All scenarios must pass cleanly.

- [ ] **Step 11: Commit**

```bash
git add source/ species/
git commit -m "ncdm: remove ncdm_id, ncdm_species_sorted_, and pv legacy fields"
```

---

## Task 19: Final audit + create PR

- [ ] **Step 1: Audit grep — final**

```bash
grep -rn "ncdm_id\|index_ncdm_\|l_max_ncdm\|q_size_ncdm\|index_pt_psi0_ncdm1\|ncdm_species_sorted_\|SetNcdmId\|Omega0_ncdm_tot\|N_decay_dr" source/ species/ include/ 2>/dev/null
```

Expected: zero hits (modulo `pba->has_ncdm` / `has_ncdm_decay_dr` which are out of scope, and any per-species member named `pt_q_size_` etc. which has a different identifier).

```bash
grep -rn "pba->N_ncdm" source/ species/ include/ 2>/dev/null
```

Expected: zero hits.

- [ ] **Step 2: Verify TODO markers landed where the spec called for them**

```bash
grep -rn "TODO(architecture)" source/ species/
```

Expected: ~5 hits in the consumers of the removed `Omega0_ncdm_tot` (per the spec, this is the only acceptable remaining type-pick).

- [ ] **Step 3: Clean build + smoke + full pytest**

```bash
cd /Users/au192734/Projects/class_claude
make clean && make class -j 2>&1 | tail -20
./class explanatory.ini 2>&1 | tail -3

for f in test/scenarios/*.ini; do
  echo "=== $f ==="
  ./class "$f" 2>&1 | tail -3
done

cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -30
cd ..
```

Expected: all green.

- [ ] **Step 4: Update CHANGELOG**

Open `CHANGELOG.md` (or wherever release notes live in this repo). Add an entry under the upcoming release:

```
- Removed `ncdm_id` integer handle. NCDM-family species (NCDM, DNCDM, DNCDM_DR
  composite, NCDM-self-interacting) now own their perturbation, transfer, and
  shooter indices internally. Background struct fields `N_ncdm`,
  `Omega0_ncdm_tot`, `N_decay_dr` are gone — sum over species directly. Output
  column titles changed: `rho_ncdm[0]` → `rho_nu1` (instance-name based).
  Public API change: `PerturbationsModule::RescaledNCDMPerturbations(int, ...)`
  → `(BaseSpecies*, ...)`.
```

- [ ] **Step 5: Push branch + create PR**

```bash
cd /Users/au192734/Projects/class_claude
git push -u origin 267-remove-ncdm-id

gh pr create \
  --title "Remove ncdm_id: decentralize NCDM-family state into species instances" \
  --body "$(cat <<'EOF'
## Summary

- Removes the integer `ncdm_id` and every flat-array sibling it indexes into.
- Each NCDM-family species (NCDM, DNCDM, DNCDM_DR composite, NCDM-self-interacting) now owns its own perturbation, transfer, and shooter indices internally.
- Module code iterates `all_species_` and dispatches via virtual hooks instead of reaching into shared per-id arrays.
- New `BaseSpecies` hooks: `StashPerturbationLayout`, `CopyPerturbationsAcrossSwitch`, `MarkUsedInSources` (perturbation switch bookkeeping); `RegisterShootingIndices`, `ComputeShootingResidual`, `ComputeShootingGuess` (shooter slot registration); `RegisterTransferIndices` (transfer-function slots).
- Background struct fields `N_ncdm`, `Omega0_ncdm_tot`, `N_decay_dr` deleted — sum from species directly.
- DNCDM_DR composite + DR child renamed to `<dncdm_instance_name>` and `<dncdm_instance_name>_DR`.
- Output columns renamed: `(.)rho_ncdm[0]` → `(.)rho_<instance>`, `delta_ncdm[0]` → `delta_<instance>`, etc. (Breaking change for downstream parsers.)
- Public API change: `PerturbationsModule::RescaledNCDMPerturbations(int, ...)` → `(BaseSpecies*, ...)`.
- Spec: `docs/superpowers/specs/2026-05-03-remove-ncdm-id-design.md`.
- Plan: `docs/superpowers/plans/2026-05-03-remove-ncdm-id.md`.

Out of scope (follow-ups):
- Migrating non-NCDM approximation-switch copy blocks (UR ufa, IDR rsa/tca, IDM_DRMD tca) into the new `CopyPerturbationsAcrossSwitch` pattern. Hook is general; follow-up is mechanical.
- Removing the `dynamic_cast<NCDMBaseSpecies*>` filter in `Omega0_ncdm_tot` consumers — flagged with `TODO(architecture)` at each site, to be redesigned alongside `Omega0_b/Omega0_cdm`.
- `pba->has_ncdm` / `pba->has_ncdm_decay_dr` has-guards — opportunistic cleanup track.

## Test plan

- [x] `make clean && make class -j` builds cleanly with no new warnings
- [x] `./class explanatory.ini` runs successfully
- [x] All five regression scenarios in `test/scenarios/ncdm_*.ini` + `dncdm_dr.ini` pass
- [x] `TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py` passes
- [x] Numerical regression diff against `master`-baseline within `1e-10` relative tolerance for physical quantities
- [x] `grep -rn "ncdm_id" source/ species/ include/` returns zero hits
EOF
)"
```

- [ ] **Step 6: Note the PR URL** — return it to the user.

---

## Spec coverage check (run after writing all tasks)

| Spec section | Implementing task(s) |
|---|---|
| §1 Per-species perturbation storage | Task 2, Task 3, Task 4 |
| §2(a) Physics loops — extend BaseSpecies virtuals | Task 7, Task 9 |
| §2(b) Bookkeeping loops — new generic hooks | Task 5, Task 6 |
| §2(c) Identity / API changes (RescaledNCDMPerturbations, ncdm_species_sorted_) | Task 18 |
| §3 Shooter index registration | Task 12, Task 13, Task 14, Task 15 |
| §4 Background-struct cleanup | Task 10, Task 11 |
| §5 Output column titles + DNCDM_DR composite naming | Task 16, Task 17 |
| §5 Transfer functions per-species | Task 9 |
| Risks: iteration-order change | Task 1 (multi-unsorted scenario), regression-diff at every behavior-touching task |
| Risks: column titles renamed | Task 17 (Step 8 explicitly handles this) |
| Risks: background.h ABI change | Task 10, Task 11 (rebuild Python wrapper validated by Task 15 + Task 19 pytest) |
| Risks: switch y-copy correctness | Task 5 (TODO markers for executor to port arithmetic), Task 6 (regression check) |
| Risks: DNCDM_DR composite-name rename ABI | Task 16 (Step 4 finds all literal lookups) |
| Risks: index_pt_psi0_ncdm1 sentinel | Confirmed dead — no consumers; Task 18 deletes it |
