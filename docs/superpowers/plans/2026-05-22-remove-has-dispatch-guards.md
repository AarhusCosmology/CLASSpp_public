# Remove `pba->has_*` dispatch guards Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 9 behavior-neutral `pba->has_*` species-presence *dispatch* guards in
`nonlinear_module.cpp` and `output_module.cpp` with container/accessor queries, so module code stops
type-picking on background flags. Keep all `class_test` validation guards and `input_module` setters.

**Architecture:** Pure control-flow-equivalent swaps in two existing files — no new files, helpers, or
headers. `has_ncdm` consumption sites use the existing `has_pk_cb_` member; the `has_pk_cb_`
derivation point uses the public `BackgroundModule::GetNcdmCount()`; `has_fld`/`has_scf` use
`all_species_.count("Fluid"|"ScalarField")`. The 3 idm_dr/idr/drmd-family guards are out of scope
(deferred to PR B — see the design doc). Equivalences are proved in
`docs/superpowers/specs/2026-05-22-remove-has-dispatch-guards-design.md`.

**Tech Stack:** C++17, CLASSpp Makefile build, classy Python wrapper, pytest scenario suite.

---

### Task 1: Capture the baseline

**Files:** none (read-only)

- [ ] **Step 1: Confirm clean tree on the PR-C branch**

Run: `git status && git branch --show-current`
Expected: clean tree, branch `remove-has-dispatch-guards`.

- [ ] **Step 2: Confirm the 9 target lines are present and unchanged**

Run:
```bash
grep -n "pba->has_ncdm\|pba->has_fld" source/nonlinear_module.cpp
grep -n "pba->has_scf" source/output_module.cpp
```
Expected (nonlinear): `has_ncdm` at 1007, 1383, 2695, 3583, 3647, 3711, 3775; `has_fld` at 3173.
Expected (output): `has_scf` at 997.
(If line numbers drift, locate by content — the edits below match on text, not line number.)

---

### Task 2: nonlinear_module — `has_ncdm` consumption sites → `has_pk_cb_`

**Files:**
- Modify: `source/nonlinear_module.cpp` (the 5 pk_cb-consumption sites: 2695, 3583, 3647, 3711, 3775)

- [ ] **Step 1: Replace the `index_pk_cb` site (~:2695)**

Old:
```cpp
  /** Test whether pk_cb has to be taken into account (only if we have massive neutrinos)*/
  if (pba->has_ncdm == _TRUE_) {
    index_pk_cb = index_pk_cb_;
  }
```
New:
```cpp
  /** Test whether pk_cb has to be taken into account (only if we have massive neutrinos)*/
  if (has_pk_cb_ == _TRUE_) {
    index_pk_cb = index_pk_cb_;
  }
```

- [ ] **Step 2: Replace the four `sigma_*_cb` sites (~:3583, 3647, 3711, 3775)**

Each of these four reads `pnw->sigma_*[index_pk_cb_]` under `if (pba->has_ncdm) {`. Replace the
guard only. Because the four are textually identical (`  if (pba->has_ncdm) {` followed by
`    if (tau_size_ == 1) {`), use a scoped replace that disambiguates by the following line, or edit
each occurrence individually. Target result for all four:
```cpp
  if (has_pk_cb_) {
    if (tau_size_ == 1) {
```
Verify exactly four `if (pba->has_ncdm) {` → `if (has_pk_cb_) {` conversions were made here (sites
3583/3647/3711/3775), and that 1007/1383 are NOT among them (handled in Task 3).

- [ ] **Step 3: Build**

Run: `make class -j 2>&1 | tail -20`
Expected: builds to completion, no new warnings referencing `nonlinear_module.cpp`.

- [ ] **Step 4: Commit**

```bash
git add source/nonlinear_module.cpp
git commit -m "nonlinear: gate pk_cb consumption on has_pk_cb_ instead of pba->has_ncdm"
```

---

### Task 3: nonlinear_module — `has_ncdm` derivation + warning loop

**Files:**
- Modify: `source/nonlinear_module.cpp` (1383 derivation; 1007 warning loop)

- [ ] **Step 1: Replace the `has_pk_cb_` derivation (~:1383)**

Old:
```cpp
  has_pk_m_ = _TRUE_;
  if (pba->has_ncdm == _TRUE_) {
    has_pk_cb_ = _TRUE_;
  }
  else {
    has_pk_cb_ = _FALSE_;
  }
```
New:
```cpp
  has_pk_m_ = _TRUE_;
  if (background_module_->GetNcdmCount() > 0) {
    has_pk_cb_ = _TRUE_;
  }
  else {
    has_pk_cb_ = _FALSE_;
  }
```
(`GetNcdmCount()` is declared public at `source/background_module.h:25`; `background_module_` is the
already-injected dependency used throughout this file.)

- [ ] **Step 2: Drop the redundant `has_ncdm` wrapper around the Halofit warning loop (~:1007)**

Old:
```cpp
    if (pba->has_ncdm) {
      for (auto& [name, sp] : all_species_) {
        auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
        if (!ncdm_sp)
          continue;
        double m_ncdm_in_electronvolt = ncdm_sp->GetMassInElectronvolt();
        if (m_ncdm_in_electronvolt > _M_EV_TOO_BIG_FOR_HALOFIT_)
          fprintf(stdout,
                  "Warning: Halofit and HMcode are proved to work for CDM, and also with a small "
                  "HDM component. But it sounds like you are running with a WDM component of mass "
                  "%f eV, which makes the use of Halofit suspicious.\n",
                  m_ncdm_in_electronvolt);
      }
    }
```
New (remove the outer `if`; the loop self-filters via `dynamic_cast` and is a no-op with no NCDM):
```cpp
    for (auto& [name, sp] : all_species_) {
      auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
      if (!ncdm_sp)
        continue;
      double m_ncdm_in_electronvolt = ncdm_sp->GetMassInElectronvolt();
      if (m_ncdm_in_electronvolt > _M_EV_TOO_BIG_FOR_HALOFIT_)
        fprintf(stdout,
                "Warning: Halofit and HMcode are proved to work for CDM, and also with a small "
                "HDM component. But it sounds like you are running with a WDM component of mass "
                "%f eV, which makes the use of Halofit suspicious.\n",
                m_ncdm_in_electronvolt);
    }
```
Note: this loop sits inside `if (pnl->method > nl_none) { ... }`, immediately before the
`if (pba->has_idm_dr)` warning at ~:1021 which is **deferred to PR B** — leave that one untouched.

- [ ] **Step 3: Confirm only the intended `has_ncdm` sites remain changed**

Run: `grep -n "pba->has_ncdm" source/nonlinear_module.cpp`
Expected: **no matches** (all 7 nonlinear `has_ncdm` sites are now migrated).
Run: `grep -n "pba->has_idm_dr\|pba->has_fld" source/nonlinear_module.cpp`
Expected: `has_idm_dr` at ~:1021 (deferred — still present); `has_fld` at ~:3173 (Task 4, still present).

- [ ] **Step 4: Build**

Run: `make class -j 2>&1 | tail -20`
Expected: builds to completion, no new warnings.

- [ ] **Step 5: Commit**

```bash
git add source/nonlinear_module.cpp
git commit -m "nonlinear: derive has_pk_cb_ from GetNcdmCount(); drop redundant has_ncdm guard on Halofit warning"
```

---

### Task 4: nonlinear_module `has_fld` + output_module `has_scf`

**Files:**
- Modify: `source/nonlinear_module.cpp` (~:3173)
- Modify: `source/output_module.cpp` (~:997)

- [ ] **Step 1: Replace `has_fld` (nonlinear ~:3173)**

Old:
```cpp
  if (pba->has_fld == _TRUE_) {
```
New:
```cpp
  if (all_species_.count("Fluid")) {
```

- [ ] **Step 2: Replace `has_scf` (output ~:997)**

Old:
```cpp
    if (pba->has_scf == _TRUE_) {
```
New:
```cpp
    if (all_species_.count("ScalarField")) {
```

- [ ] **Step 3: Confirm the deferred guards are still present and nothing else changed**

Run:
```bash
grep -n "pba->has_" source/nonlinear_module.cpp source/output_module.cpp
```
Expected: only the deferred-to-B guards remain — `nonlinear:~1021 pba->has_idm_dr`,
`output:~1047 pba->has_idm_dr`. No `has_ncdm`, `has_fld`, or `has_scf` matches.

- [ ] **Step 4: Build**

Run: `make class -j 2>&1 | tail -20`
Expected: builds to completion, no new warnings.

- [ ] **Step 5: Commit**

```bash
git add source/nonlinear_module.cpp source/output_module.cpp
git commit -m "nonlinear/output: gate fld/scf on all_species_.count() instead of pba->has_*"
```

---

### Task 5: Verify behavior is unchanged

**Files:** none (verification)

- [ ] **Step 1: Rebuild the wrapper**

Run: `make classy 2>&1 | tail -15`
Expected: builds `python/build/lib.*` and installs the classy module successfully.

- [ ] **Step 2: Run the per-PR scenario gate**

Run: `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py 2>&1 | tail -40`
Expected: all scenarios PASS (same count as on `origin/master`). The grid exercises ncdm, fld, scf.

- [ ] **Step 3: Run the dedicated reference tests covering the deferred-to-B species**

Run: `cd python && python -m pytest -v test_class.py -k "matches_reference or computes" 2>&1 | tail -40`
Expected: PASS — confirms idm_dr/idr/drmd/dcdm_dr paths are untouched by this PR.

- [ ] **Step 4: Sanity-run the binary on a representative ini**

Run: `./class explanatory.ini 2>&1 | tail -5`
Expected: completes without error.

---

### Task 6: Update the roadmap + cross-PR record

**Files:**
- Modify: `docs/superpowers/plans/2026-05-22-pr268-followups-roadmap.md` (section C; note handoff to B)

- [ ] **Step 1: Rewrite the roadmap C section** to reflect: the original plan file never existed;
the real scope was 9 clean guards (done here) + 3 idm_dr/idr/drmd-family guards handed to B; point at
this plan + design doc. Add to the B section that B now also owns `background:532`, `nonlinear:~1021`,
`output:~1047` and must add per-sub-species presence accessors to the composites.

- [ ] **Step 2: Commit**

```bash
git add docs/superpowers/
git commit -m "docs: PR C done (9 has_* dispatch guards); hand 3 idm_dr/drmd guards to B"
```

---

## Self-Review notes

- **Spec coverage:** all 9 in-scope guards (nonlinear 1007/1383/2695/3583/3647/3711/3775/3173,
  output 997) have a task; the 3 deferred guards are explicitly left in place and documented.
- **Type consistency:** `has_pk_cb_` (member, `nonlinear_module.h:72`), `GetNcdmCount()` (public,
  `background_module.h:25`), `all_species_.count(std::string)` (existing API, used 187× in source/).
- **No new files** → build-system lists untouched.
