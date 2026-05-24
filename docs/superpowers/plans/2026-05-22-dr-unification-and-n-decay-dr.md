# DR species unification + remove `N_decay_dr` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Collapse `DarkRadiationSpecies` to a single self-contained DR channel and make it the sole DR
class (deleting `DNCDM_DecayRadiationSpecies`), removing `pba->N_decay_dr`, the shared
`index_pt_F0_dr_sum`/`index_pt_F0_dr_species` perturb_vector fields, and the per-channel layout — then
migrate the deferred DCDM_DR switch-copy/IC into the species hook, and migrate the 3 inherited
idm_dr/drmd guards via composite presence accessors.

**Architecture:** Each decay channel owns one `DarkRadiationSpecies` (a composite needing N products
holds N instances); no shared DR. The DR species evolves a single free-streaming hierarchy + its own
decay source (injected by the parent composite's `AddCouplingDerivs`), so the total DR is the
stress-energy loop's per-species sum (no materialized "sum"). The synchronous-gauge output correction
`(ρ̇_dr/ρ_dr)·α` is computed generically from a decay-source functor the composite supplies, so the DR
species holds no decay-channel-specific code. Design + the latent shared-sum bug are in
`docs/superpowers/specs/2026-05-22-dr-unification-and-n-decay-dr-design.md`.

**Tech Stack:** C++17, CLASSpp Makefile + classy, pytest scenario suite, `test/scenarios/compare_tol.py`.

**Verification model (per spec):** **dcdm-only is the strict ~0.1% behavior baseline** (the old "sum"
was a redundant copy of the one channel). **dncdm-only and combined are correct-by-construction
(option a)** — the migration fixes the legacy double-count, so we do NOT match the old reference there;
sanity-check shape/magnitude and regenerate those references (author eyeballs).

---

## File structure

| File | Responsibility after this PR |
|---|---|
| `species/dark_radiation_species.{h,cpp}` | The sole, generic, single-channel DR species: dilution + free-streaming hierarchy + ICs + `Delta`/`Theta`/`DeltaP`/`RhoPlusPShear`, layout `{idx_F0, l_max}`, generic gauge-correction via a decay-source functor, per-instance output names, `CopyPerturbationsAcrossSwitch` |
| `species/dncdm_decay_radiation_species.{h,cpp}` | **deleted** |
| `species/dcdm_dr_species.{h,cpp}` | Owns one `DarkRadiationSpecies`; supplies its decay-source functor; `AddCouplingDerivs` writes only its own `idx_F0`; drops `N_decay_dr` output loops |
| `species/dncdm_dr_species.{h,cpp}` | Owns one `DarkRadiationSpecies` (instead of `DNCDM_DecayRadiationSpecies`); `AddCouplingDerivs` writes only its own `idx_F0` |
| `source/background.h` | `N_decay_dr` field removed |
| `source/input_module.cpp` | `N_decay_dr` setter removed |
| `source/background_module.{h,cpp}` | `index_bg_rho_dr_species_`/`index_bg_rho_dr_` collapse to a single total slot owned by each DR species; guard `:532` → composite accessors |
| `source/perturbations.h` | `index_pt_F0_dr_sum`, `index_pt_F0_dr_species`, `l_max_dr`, `l_max_dr_col` perturb_vector fields removed |
| `source/perturbations_module.cpp` | DR registration via the single layout; switch-copy/IC blocks removed (now species hooks) |
| `species/idm_dr_idr_species.h`, `species/idm_drmd_idr_drmd_species.h` | presence accessors |
| `source/nonlinear_module.cpp`, `source/output_module.cpp` | 2 `has_idm_dr` guards → composite accessors |
| `Makefile`, `setup.py`, `CLASS.xcodeproj/project.pbxproj` | drop `dncdm_decay_radiation_species` |

---

# Chunk A — Single-channel generic DarkRadiationSpecies + remove `N_decay_dr`

**This chunk lands as one commit** — intermediate states don't compile (the shared sum is referenced
across the species, both composites, and the module). Do the edits in order, then build/verify/commit
at the end (Task A9).

### Task A1: Rewrite `DarkRadiationSpecies` header to single-channel + generic gauge hook

**Files:** Modify `species/dark_radiation_species.h`

- [ ] **Step 1: Replace the `PerturbLayout`, constructor, members, and accessors.**

Use `DNCDM_DecayRadiationSpecies` (the already-modern single-channel class) as the structural template.
Key changes vs the current `DarkRadiationSpecies`:
- Layout becomes `struct PerturbLayout : BaseSpecies::PerturbLayout { int idx_F0 = -1; int l_max = -1; };`
  (drop `idx_F0_sum`/`idx_F0_species`).
- Drop the `const DCDMSpecies* dcdm_` member and the `(dcdm, pba, bgm)` constructor. New constructor
  takes an instance name + `pba`/`bgm`:
  ```cpp
  DarkRadiationSpecies(const std::string& name, const background* pba, const BackgroundModule* bgm)
      : BaseSpecies(name, EnergyType::Radiation), pba_(pba), bgm_(bgm) {}
  ```
- Add the generic decay-source hook (used only for the output gauge correction), defaulting to none:
  ```cpp
  // Returns the decay-source rate added to dy[rho_dr] by the parent composite
  // (e.g. a*Gamma_dcdm*rho_dcdm). Used to form the synchronous-gauge output
  // correction (rho_dot_dr/rho_dr)*alpha. Default (unset) => pure dilution.
  void SetDecaySourceRate(std::function<double(const double* pvecback)> f) {
    decay_source_rate_ = std::move(f);
  }
  ```
  (add `#include <functional>` and `#include <string>`).
- Replace the per-channel index members/accessors with single-slot ones, reusing the base-class
  `index_bg_rho_`:
  ```cpp
  int bg_rho_index() const { return index_bg_rho_; }
  int bi_rho_index() const { return index_bi_rho_; }
  int pt_F0_index()  const { return index_pt_F0_; }   // for any external reference
  ...
 private:
  const background* pba_;
  const BackgroundModule* bgm_;
  std::function<double(const double* pvecback)> decay_source_rate_;
  int index_bi_rho_ = -1;
  int index_pt_F0_  = -1;
  ```
  Remove `index_bg_rho_dr_species_`, `index_bi_rho_dr_species_`, `index_pt_F0_dr_species_`,
  `bg_rho_dr_species_index()`, `bi_rho_dr_species_index()`.
- Keep the same set of overridden virtuals as `DNCDM_DecayRadiationSpecies` exposes (the layout-based
  perturbation interface + legacy no-op scalar overloads), plus `WriteOutputColumns`/`PrintVariables`
  and `CopyPerturbationsAcrossSwitch` (added in A4). Keep `Rho`/`P`/`DpDloga` reading `index_bg_rho_`.

### Task A2: Rewrite `DarkRadiationSpecies` background + perturbation bodies (single channel)

**Files:** Modify `species/dark_radiation_species.cpp`

- [ ] **Step 1: Replace background index/compute/derivs with the single-slot forms** (identical to
`dncdm_decay_radiation_species.cpp:10-33`):
```cpp
void DarkRadiationSpecies::RegisterBackgroundIndices(int& index_bg)  { index_bg_rho_ = index_bg++; }
void DarkRadiationSpecies::RegisterIntegrationIndices(int& index_bi) { index_bi_rho_ = index_bi++; }
void DarkRadiationSpecies::ComputeBackground(double, const double* pvecback_B, double* pvecback) {
  pvecback[index_bg_rho_] = pvecback_B[index_bi_rho_];
}
void DarkRadiationSpecies::BackgroundDerivs(double, const double* y, double* dy, const double* pvecback) {
  const double a = pvecback[bgm_->index_bg_a_];
  const double H = pvecback[bgm_->index_bg_H_];
  dy[index_bi_rho_] = -4. * a * H * y[index_bi_rho_];   // dilution; source added by parent composite
}
```

- [ ] **Step 2: Replace `RegisterPerturbationIndices` with a single hierarchy** (mirrors
`dncdm_decay_radiation_species.cpp:37-47`), writing only the layout (no `pv->index_pt_*`):
```cpp
void DarkRadiationSpecies::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base, perturb_vector* pv,
    const precision* ppr, int& index_pt, const perturb_workspace*, int) {
  auto& layout   = static_cast<PerturbLayout&>(base);
  layout.idx_F0  = index_pt;
  layout.l_max   = ppr->l_max_dr;
  index_pt_F0_   = index_pt;
  index_pt      += ppr->l_max_dr + 1;
}
```
(Read `l_max_dr` from `ppr` — the precision config — since the `pv->l_max_dr` copy is removed in Task A6.)

- [ ] **Step 3: Replace `PerturbDerivs` with the single-hierarchy free-streaming body** — identical to
`dncdm_decay_radiation_species.cpp:49-95` (uses `layout.idx_F0`, `layout.l_max`, `r_dr` from
`pvecback[index_bg_rho_]`). Remove the old sum/accumulate and the `if (pba_->N_decay_dr > 0)` guard.

- [ ] **Step 4: Replace `Delta`/`Theta`/`DeltaP`/`RhoPlusPShear` (layout overloads)** with the
`idx_F0` forms — identical to `dncdm_decay_radiation_species.cpp:97-152`. Keep the legacy scalar
overloads as no-ops/`0.` returns (as in the dncdm header).

- [ ] **Step 5: Replace `ApplyInitialConditions` with the single-hierarchy IC** (the `idx_F0` branch of
the current `dark_radiation_species.cpp:222-256`, dropping the per-channel loop and the sum branch):
```cpp
void DarkRadiationSpecies::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base, double* y,
                                                  const PerturbIcContext& ctx) {
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_) return;
  const auto& layout     = static_cast<const PerturbLayout&>(base);
  if (layout.idx_F0 < 0) return;
  const double* pvecback = ctx.ppw->pvecback;
  const double r_prefactor = std::pow(std::pow(ctx.a / pba_->a_today, 2) / pba_->H0, 2);
  const double r_dr = r_prefactor * pvecback[index_bg_rho_];
  y[layout.idx_F0 + 0] = ctx.delta_dr * r_dr;
  if (layout.l_max >= 1) y[layout.idx_F0 + 1] = 4. / (3. * ctx.k) * ctx.theta_ur * r_dr;
  if (layout.l_max >= 2) y[layout.idx_F0 + 2] = 2. * ctx.shear_ur * r_dr;
  if (layout.l_max >= 3) y[layout.idx_F0 + 3] = ctx.l3_ur * r_dr;
}
```

### Task A3: Generic gauge correction + per-instance output in `DarkRadiationSpecies::PrintVariables`

**Files:** Modify `species/dark_radiation_species.cpp`, `species/dark_radiation_species.h`

- [ ] **Step 1: Rewrite `WriteOutputColumns`/`PrintVariables`** for a single channel, deriving column
names from `name()` and computing the gauge correction generically:
```cpp
void DarkRadiationSpecies::PrintVariables(PerturbColumnWriter& w, double, const double* y,
                                          const PerturbationsModule& mod, const perturb_workspace* ppw) const {
  const int l_max_dr = mod.GetPrecision()->l_max_dr;
  double k = 0., a = 0., H = 0., alpha_corr = 0., rho_dot_over_rho = 0.;
  int base = -1;
  if (!w.IsTitleMode()) {
    const perturb_vector* pv = ppw->pv;
    base = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]).idx_F0;
    k = ppw->scalar_ctx.k;
    a = ppw->pvecback[bgm_->index_bg_a_];
    H = ppw->pvecback[bgm_->index_bg_H_];
    const double rho_dr = ppw->pvecback[index_bg_rho_];
    if (mod.GetPerturbs()->gauge == synchronous && rho_dr > 0.) {
      alpha_corr = ppw->pvecmetric[ppw->index_mt_alpha];
      const double src = decay_source_rate_ ? decay_source_rate_(ppw->pvecback) : 0.;
      rho_dot_over_rho = -4. * a * H + src / rho_dr;   // = rho_dot_dr / rho_dr
    }
  }
  // delta/theta/shear + raw F multipoles, named from this instance
  ... (single channel: drop the `for n_dr` loops; delta = y[base]/r_dr + rho_dot_over_rho*alpha_corr, etc.)
}
```
Use `name()` (e.g. `"DR"`, `"<flavor>_DR"`) to build column titles instead of the hardcoded
`delta_ncdm[%d]`/`F_dr[%d][%d]`. (`collection_index_` is the base-class member already used by other
species to find their layout in `pv->species_layouts`.)

- [ ] **Step 2:** Confirm no remaining `pba_->N_decay_dr`, `index_bg_rho_dr_species_`,
`index_pt_F0_dr_species`, or `idx_F0_sum` references in the file:
Run: `grep -nE "N_decay_dr|rho_dr_species_|F0_dr_sum|idx_F0_species|idx_F0_sum" species/dark_radiation_species.cpp`
Expected: no matches.

### Task A4: Add `CopyPerturbationsAcrossSwitch` to `DarkRadiationSpecies`

**Files:** Modify `species/dark_radiation_species.{h,cpp}`

- [ ] **Step 1:** Declare the override (mirror `species/idm_dr.h:127-131`) and implement it as a
straight copy of the single hierarchy (this absorbs `perturbations_module.cpp:3866-3877`):
```cpp
void DarkRadiationSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
    const BaseSpecies::PerturbLayout& new_base, const double* old_y, double* new_y,
    const PerturbSwitchContext&) const {
  const auto& o = static_cast<const PerturbLayout&>(old_base);
  const auto& n = static_cast<const PerturbLayout&>(new_base);
  if (o.idx_F0 < 0 || n.idx_F0 < 0) return;
  for (int l = 0; l <= n.l_max; ++l) new_y[n.idx_F0 + l] = old_y[o.idx_F0 + l];
}
```
(Inspect `PerturbSwitchContext` for the exact reconduct convention used by the sibling species; match it.)

### Task A5: Rewire `DCDM_DR_Species` and `DNCDM_DR_Species` to the single DR class

**Files:** Modify `species/dcdm_dr_species.{h,cpp}`, `species/dncdm_dr_species.{h,cpp}`

- [ ] **Step 1 (DCDM_DR construct + decay-source functor):** in `DCDM_DR_Species` ctor, construct
`DarkRadiationSpecies("DR", pba, bgm)` and set its decay-source rate to the same expression the
composite adds in `BackgroundDerivs`:
```cpp
auto dr_sp = std::make_unique<DarkRadiationSpecies>("DR", pba, bgm);
dr_sp->SetDecaySourceRate([this](const double* pv) {
  return pv[bgm_->index_bg_a_] * pba_->Gamma_dcdm * dcdm_->Rho(pv);
});
```
(Capture is valid for the composite's lifetime; `dcdm_`/`bgm_`/`pba_` are members.)

- [ ] **Step 2 (DCDM_DR derivs/IC/output):** update `DCDM_DR_Species::BackgroundDerivs` to add the
source to `dr_sp_->bi_rho_index()` (was `bi_rho_dr_species_index()`); in `SetBackgroundInitialConditions`
use `bi_rho_index()`; in `AddCouplingDerivs` write only `dy[base + 0/1]` where `base = my.dr.idx_F0`
(remove the `pv->index_pt_F0_dr_sum` writes); delete the `for (j < N_decay_dr)` loops in
`WriteBackgroundColumnTitles`/`WriteBackgroundData` (single `(.)rho_dr` line) and the
`pvecback[bgm_->index_bg_rho_dr_]` reads (use `dr().Rho(pvecback)` / `bg_rho_index()`).

- [ ] **Step 3 (DNCDM_DR):** replace `DNCDM_DecayRadiationSpecies` with `DarkRadiationSpecies` in
`dncdm_dr_species.{h,cpp}` (member type, include, ctor `make_unique<DarkRadiationSpecies>(name+"_DR", pba, bgm)`),
set its decay-source functor to `a*Gamma*M_ncdm*number` (matching `DNCDM_DR::BackgroundDerivs`), and in
`AddCouplingDerivs` write the collision term only to `dy[base + l]` (`base = my.dr.idx_F0`) — **remove the
`dy[pv->index_pt_F0_dr_sum + l]` write** (`dncdm_dr_species.cpp:262`). The `PerturbLayout` nested member
type changes from `DNCDM_DecayRadiationSpecies::PerturbLayout dr;` to `DarkRadiationSpecies::PerturbLayout dr;`.

### Task A6: Remove `N_decay_dr` + parallel perturb_vector fields

**Files:** Modify `source/background.h`, `source/input_module.cpp`, `source/perturbations.h`,
`source/perturbations_module.cpp`, `source/background_module.{h,cpp}`

- [ ] **Step 1:** Delete `int N_decay_dr = 0;` (`background.h:133`) and its setter
(`input_module.cpp:341-343`, incl. the stale comment).
- [ ] **Step 2:** Delete the perturb_vector fields `index_pt_F0_dr_sum`, `index_pt_F0_dr_species`,
`l_max_dr`, `l_max_dr_col` (`perturbations.h:284-287`). Keep the **precision** `l_max_dr`/`l_max_dr_col`
(`common.h:736-738`, `input_module.cpp:3257-3258`) and the input validation
(`perturbations_module.cpp:619-623`).
- [ ] **Step 3:** In `perturbations_module.cpp`, the DCDM_DR/DNCDM_DR DR registration (`:3350-3401`)
collapses: each composite registers its DR child's single-hierarchy layout via the normal
species-dispatch (no `index_pt_F0_dr_sum` repointing, no `ppv->l_max_dr` assignment). Remove the
`need_sum_index` logic (`:3380-3398`) and the legacy second-overload call (`:3357`) where it only served
the sum.
- [ ] **Step 4:** In `background_module.{h,cpp}`, remove `index_bg_rho_dr_species_`/`index_bg_rho_dr_`
(`background_module.h:52-53`, `.cpp:696-701`) — each DR species owns its single `index_bg_rho_`; the
`dcdm_dr` reads at `dcdm_dr_species.cpp:219/232` use `dr().bg_rho_index()`/`dr().Rho()` instead.

### Task A7: Migrate the DCDM_DR switch-copy + IC blocks out of `perturbations_module`

**Files:** Modify `source/perturbations_module.cpp`

- [ ] **Step 1:** Delete the `if (all_species_.count("DCDM_DR"))` switch-copy block
(`perturbations_module.cpp:3866-~3878`, the `index_pt_F0_dr_sum`/`index_pt_F0_dr_species` copy loops) —
now handled by `DarkRadiationSpecies::CopyPerturbationsAcrossSwitch` (Task A4) via the generic
per-species reconduct dispatch. Confirm the generic dispatch invokes child species hooks for composites
(PR A wired `IDM_DR_IDR`/`IDM_DRMD_IDR_DRMD` to delegate; verify `DCDM_DR`/`DNCDM_DR` delegate to their
DR child — add delegation if missing, mirroring the idm composites).
- [ ] **Step 2:** Delete the `if (all_species_.count("DCDM_DR"))` IC block
(`perturbations_module.cpp:4931-~4952`, the `for n_dr < N_decay_dr` + `index_pt_F0_dr_sum` writes) — now
handled by `DarkRadiationSpecies::ApplyInitialConditions` (Task A2 Step 5) via the composite's
`ApplyInitialConditions` delegation.
- [ ] **Step 3:** `grep -n "N_decay_dr\|index_pt_F0_dr" source/ -r` → expect **no matches** anywhere.

### Task A8: Delete `DNCDM_DecayRadiationSpecies` + update build systems

**Files:** Delete `species/dncdm_decay_radiation_species.{h,cpp}`; modify `Makefile`, `setup.py`,
`CLASS.xcodeproj/project.pbxproj`

- [ ] **Step 1:** `git rm species/dncdm_decay_radiation_species.h species/dncdm_decay_radiation_species.cpp`
- [ ] **Step 2:** Remove `dncdm_decay_radiation_species` from `Makefile` (SPECIES_OPP), `setup.py`
(source_files), and `CLASS.xcodeproj/project.pbxproj` (user verifies Xcode). Confirm no remaining
references: `grep -rn "dncdm_decay_radiation" . --include=Makefile --include=setup.py --include=*.pbxproj`
→ none.

### Task A9: Build + verify + commit (Chunk A)

- [ ] **Step 1:** `make class -j 2>&1 | tail -20` → builds, no new warnings. Fix compile errors before
proceeding (the cross-file edits must all be consistent).
- [ ] **Step 2:** `make classy 2>&1 | tail -5` → wrapper rebuilt.
- [ ] **Step 3 (dcdm-only strict baseline):** run the dcdm_dr reference regression —
`cd python && python -m pytest -v test_class.py -k "dcdm_dr_matches_reference"` → **PASS** (physics
identical; if it fails, the dcdm path is not behavior-preserving — debug before continuing).
- [ ] **Step 4 (no crashes across the grid):**
`cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py` →
all pass (no crashes). dncdm/combined scenarios run; their *values* may change (expected).
- [ ] **Step 5 (dncdm/combined — option a):** run `./class test/scenarios/dncdm_dr.ini` and a combined
dcdm+dncdm ini; sanity-check `d_dr`/`t_dr` shape+magnitude (standard free-streaming + decay source per
channel). Surface the dncdm_dr output diff for the author to eyeball; regenerate the `dncdm_dr` (and any
combined) reference if accepted.
- [ ] **Step 6: Commit.**
```bash
git add -A
git commit -m "DR: single generic DarkRadiationSpecies; remove N_decay_dr + shared sum"
```

---

# Chunk B — IDM composite presence accessors + 3 inherited guards

Independent of Chunk A. Behaviour-neutral (the accessors return exactly the `pba` flags).

### Task B1: Add presence accessors to the IDM composites

**Files:** Modify `species/idm_dr_idr_species.h`, `species/idm_drmd_idr_drmd_species.h`

- [ ] **Step 1 (IDM_DR_IDR):** capture the flags at construction and expose accessors:
```cpp
// in ctor body (idm_dr_idr_species.cpp): has_idm_dr_ = pba.has_idm_dr; has_idr_ = pba.has_idr;
bool has_idm_dr() const { return has_idm_dr_; }   // == pba->has_idm_dr at build time
bool has_idr()    const { return has_idr_; }
// private: bool has_idm_dr_ = false, has_idr_ = false;
```
- [ ] **Step 2 (IDM_DRMD_IDR_DRMD):** likewise `has_idm_drmd()` / `has_idr_drmd()` capturing
`pba.has_idm_drmd` / `pba.has_idr_drmd`.

### Task B2: Migrate `background_module.cpp:532`

**Files:** Modify `source/background_module.cpp`

- [ ] **Step 1:** Replace
```cpp
  if ((pba->has_idr_drmd) && (pba->has_idm_drmd)) {
```
with (guarded by composite presence, then the per-sub-species flags):
```cpp
  if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
    auto& c = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
    if (c.has_idr_drmd() && c.has_idm_drmd()) {
```
(close the extra brace; include the composite header if not already). Equivalent to the old
`has_idr_drmd && has_idm_drmd`.

### Task B3: Migrate the two `has_idm_dr` guards

**Files:** Modify `source/nonlinear_module.cpp` (~:1019), `source/output_module.cpp` (~:1047)

- [ ] **Step 1:** Each `if (pba->has_idm_dr ...)` → guard on the composite presence + accessor:
```cpp
  if (all_species_.count("IDM_DR_IDR") &&
      static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).has_idm_dr()) {
```
(include `species/idm_dr_idr_species.h` where needed). This restores the exact `has_idm_dr` semantics
(true only when idm_dr is present, not for idr-only).

### Task B4: Build + verify + commit (Chunk B)

- [ ] **Step 1:** `make class -j 2>&1 | tail -20` → clean.
- [ ] **Step 2:** `grep -rn "pba->has_idm_dr\|pba->has_idr_drmd\|pba->has_idm_drmd" source/background_module.cpp source/nonlinear_module.cpp source/output_module.cpp` → no *dispatch* matches remain (only `class_test` validation, if any).
- [ ] **Step 3:** `cd python && python -m pytest -v test_class.py -k "idm or idr or drmd or computes"` →
PASS (idr-only / drmd-partial confirm the guards weren't widened).
- [ ] **Step 4: Commit.**
```bash
git add source/
git commit -m "idm composites: presence accessors; migrate 3 idm_dr/drmd dispatch guards off pba->has_*"
```

---

## Execution findings (discovered while implementing — supersede earlier assumptions)

Investigation of the live dispatch corrected three plan assumptions:

1. **`DarkRadiationSpecies::PrintVariables` is dead code.** `PrintVariables` is dispatched on
   top-level `all_species_` only (`perturbations_module.cpp:6672`), and `DCDM_DR`/`DNCDM_DR` do **not**
   override it (composites don't delegate it). So the DR perturbation-dump columns aren't produced
   today. ⇒ stub `DarkRadiationSpecies::PrintVariables` as a no-op (behavior-preserving); **no
   layout-aware PrintVariables dispatch rework is needed.**
2. **The synchronous-gauge `decay_corr` was in that dead path**, so the **decay-source functor is
   unnecessary** — drop `SetDecaySourceRate`. The *live* transfer output already carries its own
   N-body gauge correction (`+ 4·a'/a·θ/k²`) in the composite's `FillSources`, which is channel-agnostic.
3. **`DarkRadiationSpecies::WriteOutputColumns` is also dead** (same top-level-only dispatch at
   `:301-305`, base no-op, composites don't override). The DR transfer *value* is written by
   `DCDM_DR::FillSources` (live) but its title isn't written at all today (half-dead output).

**Corrected output approach (replaces Task A3 + the output parts of A5):** the **composites own DR
output**. `DCDM_DR` and `DNCDM_DR` each:
- override `WriteOutputColumns` to write `"d_" + dr().name()` / `"t_" + dr().name()` at
  `mod.index_tp_delta_dr_ + dr().source_slot()` and `index_tp_theta_dr_ + dr().source_slot()`
  (gated on `has_density_transfers`/`has_velocity_transfers` + section), and
- write the value in `FillSources` at the same per-slot index (from their own DR `idx_F0`, keeping the
  `+ 4·a'/a·θ/k²` N-body correction).
- `DarkRadiationSpecies::WriteOutputColumns`/`PrintVariables` become no-op stubs.

**`source_slot_` assignment (mirror NCDM `:680-688`):** in `perturb_indices_of_perturbs`, after the
ncdm slot loop, assign DR slots in lex order — `for (auto& [name, sp] : all_species_) { if DCDM_DR →
d->dr().SetSourceSlot(slot++); else if DNCDM_DR → d->dr().SetSourceSlot(slot++); }` — and size
`class_define_index(index_tp_delta_dr_, has_source_delta_dr_, index_type, n_dr_species)` (and
`index_tp_theta_dr_`) where `n_dr_species` = that DR count (replaces the old single-slot define at
`:1114/:1133`). This count is the same enumeration, **not** `pba->N_decay_dr`.

## Self-Review

**Spec coverage:**
- Single-channel collapse + sole class + delete `DNCDM_DecayRadiationSpecies` → A1–A5, A8. ✓
- Generic `ρ̇_dr/ρ_dr` gauge correction (functor) → A1 Step (hook), A3, A5 Step 1. ✓
- Remove `N_decay_dr` → A6 Step 1. ✓
- Remove parallel perturb_vector fields + sum → A6 Step 2-3, A7. ✓
- Switch-copy/IC into species hooks → A4, A2 Step 5, A7. ✓
- No-sharing (each composite owns its DR; collision into own hierarchy) → A5 Steps 2-3. ✓
- Total `d_dr`/`t_dr` via stress-energy loop → consequence of removing the sum; verify in A9 Step 5. ✓
- IDM accessors + 3 guards → B1–B3. ✓
- Build systems → A8. ✓
- Verification model (dcdm strict / dncdm-combined option a) → A9 Steps 3-5. ✓

**Placeholder scan:** The DR evolution bodies are specified by exact equivalence to the existing
`dncdm_decay_radiation_species.cpp` methods (a real, in-repo template — not a vague "similar to"); the
novel logic (gauge functor, single-channel ICs, switch-copy hook, accessors) is given as concrete code.
`PrintVariables` Step shows the gauge-correction core and points at the per-channel loops to drop.

**Type consistency:** `DarkRadiationSpecies::PerturbLayout {idx_F0, l_max}`, `bg_rho_index()`/
`bi_rho_index()`/`pt_F0_index()`, `SetDecaySourceRate(std::function<double(const double*)>)`,
`has_idm_dr()/has_idr()/has_idm_drmd()/has_idr_drmd()` are used consistently across tasks.

**Open items to resolve during A (flagged inline):** the exact `PerturbSwitchContext` reconduct
convention (A4); whether `DCDM_DR`/`DNCDM_DR` already delegate `CopyPerturbationsAcrossSwitch`/
`ApplyInitialConditions` to children or need wiring (A7 Step 1); per-instance output column names (A3).
