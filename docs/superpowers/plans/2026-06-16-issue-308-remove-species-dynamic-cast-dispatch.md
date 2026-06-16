# Remove `dynamic_cast` species dispatch from modules (#308) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the per-species `dynamic_cast` *contribution* loops in the modules with neutral-default polymorphic methods on `BaseSpecies` that composites forward/sum to their children, and delete the dead `get_ncdm()` introspection helper — so module code loops `for sp : all_species_` and never downcasts.

**Architecture:** Each offending site is a loop that downcasts to an NCDM/DR family type to read a per-species contribution. For each, add a `BaseSpecies` virtual with a neutral default (0 / no-op / true), override it on the concrete species (`NCDMBaseSpecies`, `DarkRadiationSpecies`), and override it on `CompositeSpecies` to forward/sum over `children_` (mirroring the existing `GetRadiationOmega0` pattern at `species/composite_species.h:49-54`). The module loop then calls the method on every species directly. The NCDM *detection/enumeration* layer (`HasNcdm`, `NcdmFamily`, `GetNcdmSpecies` counts, `DrSpeciesCount`, the background index-reg loop, the `SetSourceSlot` loops) is deliberately left untouched for #309/#310.

**Tech Stack:** C++17 (CLASSpp species/modules), Cython wrapper (`classy.pyx`), build via scikit-build-core (`pip install .`), tests via pytest (`python/test_class.py`, `-m test_scenario`).

**Spec:** `docs/superpowers/specs/2026-06-16-issue-308-remove-species-dynamic-cast-dispatch-design.md`

---

## Test strategy (read first)

This is a **behaviour-preserving refactor** in all but one spot, so the existing
integration suite is the safety net (characterization testing), not new unit tests:

- The background-side dissolves are exactly behaviour-preserving because
  `GetNcdmSpecies` (background_module.cpp:95-104) *already* reaches the DNCDM child;
  composite-sum-children forwarding hits the identical set.
- The **only** behaviour change is the perturbations-side relativistic-IC check
  (Task 5): the old narrow `dynamic_cast<NCDMSpecies*>` excluded the wrapped DNCDM
  child; forwarding now includes it. Task 5 adds a new DNCDM_DR scenario test as the
  guard/decision instrument: if it fails (the child genuinely isn't relativistic at
  IC), fall back to top-level-only for that site.
- Tensor sites cannot change behaviour: tensors + DNCDM is forbidden by a
  `class_test` (perturbations_module.cpp:5916).

**Build (run from repo root after each task's edits):**
```bash
pip install . --no-build-isolation
```
**Targeted test (from repo root):**
```bash
cd python && python -m pytest test_class.py -k "<pattern>" -v
```
Reference-comparison tests (`*_matches_reference`) **skip** if `classyref` is not
installed; that is expected for local runs. The final task runs the full scenario
suite.

No C++ unit tests are added: this codebase verifies species behaviour through the
Python integration suite (the prior greybody C++ unit test was deliberately removed
in favour of `test_class.py`).

---

## File structure

**Modified — headers (add neutral-default virtuals / overrides):**
- `species/base_species.h` — new `BaseSpecies` virtual hooks (one block).
- `species/composite_species.h` — `CompositeSpecies` forward/sum overrides.
- `species/dark_radiation_species.h` — `DarkRadiationRhoToday` override.
- `species/ncdm_base_species.h` — NCDM overrides (decls).

**Modified — implementations:**
- `species/ncdm_base_species.cpp` — NCDM override bodies.
- `source/background_module.cpp` / `source/background_module.h` — dissolve
  `Omega0_dr`, `GetOmega0NcdmTot`, N_eff / mass-info loops; delete 4 indexed
  accessors.
- `source/perturbations_module.cpp` — dissolve relativistic-IC checks,
  tensor-massless `3·P`, `WriteTensorOutputColumnTitles` loop.
- `source/nonlinear_module.cpp` — dissolve Halofit WDM warning.
- `classy.pyx` — delete `get_ncdm()`.

**Modified — tests:**
- `python/test_class.py` — new DNCDM_DR compute test (Task 5).

**Left untouched (deferred to #309/#310):** `HasNcdm`, `NcdmFamily`,
`GetNcdmSpecies`, `DrSpeciesCount`, the background NCDM index-reg loop
(background_module.cpp:633-651), the `SetSourceSlot` loops
(perturbations_module.cpp:642-664), the nested-layout tensor-GW sites
(perturbations_module.cpp:5285-5293, 5918-5934), and the `Omega0_dcdm_` by-name
`DCDM_DR` read.

---

## Task 1: Delete the dead `get_ncdm()` introspection helper (Bucket 2)

`get_ncdm()` has zero callers in the repo; deleting it makes its four backing
accessors dead. `GetNcdmCount` stays (used by nonlinear_module.cpp:1249).

**Files:**
- Modify: `classy.pyx` (remove `get_ncdm`, ~lines 1068-1097)
- Modify: `source/background_module.h:27-30` (remove 4 decls)
- Modify: `source/background_module.cpp:1515-1529` (remove 4 defs)

- [ ] **Step 1: Delete `get_ncdm()` from `classy.pyx`**

Remove the whole method (from the `cpdef get_ncdm(self):` line through its final
`return ncdm_dict`):

```cython
    cpdef get_ncdm(self):
        """
        Return an array of data related to the NonColdDarkMatter module.
        ...
        """
        cdef:
            dict ncdm_dict
            int q_size

        ibg = deref(self._thisptr).GetBackgroundModule()

        ncdm_dict = {}
        for ncdm_id in range(deref(ibg).GetNcdmCount()):
            q_size = deref(ibg).GetNcdmQSize(ncdm_id)
            ncdm_dict[f"deg[{ncdm_id}]"]      = deref(ibg).GetNcdmDeg(ncdm_id)
            ncdm_dict[f"m_ncdm[{ncdm_id}]"]   = deref(ibg).GetNcdmMassInEV(ncdm_id)
            ncdm_dict[f"q_size[{ncdm_id}]"]   = q_size
            for q_id in range(q_size):
                ncdm_dict[f"q[{ncdm_id}][{q_id}]"] = deref(ibg).GetNcdmQ(ncdm_id, q_id)

        return ncdm_dict
```

- [ ] **Step 2: Remove the 4 indexed accessor declarations from `background_module.h`**

Delete these lines (keep `int GetNcdmCount() const;`):

```cpp
  double GetNcdmDeg(int n) const;
  double GetNcdmMassInEV(int n) const;
  int GetNcdmQSize(int n) const;
  double GetNcdmQ(int n, int q_id) const;
```

- [ ] **Step 3: Remove the 4 accessor definitions from `background_module.cpp`**

Delete (keep `GetNcdmCount` at 1511-1513 and `GetOmega0NcdmTot` at 1531-1536):

```cpp
double BackgroundModule::GetNcdmDeg(int n) const {
  return GetNcdmSpecies(all_species_).at(n)->GetDeg();
}

double BackgroundModule::GetNcdmMassInEV(int n) const {
  return GetNcdmSpecies(all_species_).at(n)->GetMassInElectronvolt();
}

int BackgroundModule::GetNcdmQSize(int n) const {
  return GetNcdmSpecies(all_species_).at(n)->q_size();
}

double BackgroundModule::GetNcdmQ(int n, int q_id) const {
  return GetNcdmSpecies(all_species_).at(n)->q()[q_id];
}
```

- [ ] **Step 4: Build**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly; `cclassy.pxd` regenerates without the deleted methods.

- [ ] **Step 5: Verify the wrapper imports and computes**

Run: `cd python && python -c "from classy import Class; c=Class(); c.set({'output':'tCl'}); c.compute(); print('ok')"`
Expected: prints `ok` (and `get_ncdm` is gone: `python -c "from classy import Class; print(hasattr(Class(),'get_ncdm'))"` prints `False`).

- [ ] **Step 6: Commit**

```bash
git add classy.pyx source/background_module.h source/background_module.cpp
git commit -m "$(cat <<'EOF'
species: delete dead get_ncdm() introspection helper and indexed accessors (#308)

get_ncdm() had no callers in-tree; removing it makes GetNcdmQ/GetNcdmDeg/
GetNcdmMassInEV/GetNcdmQSize dead. GetNcdmCount stays (nonlinear_module use).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: `DarkRadiationRhoToday` — dissolve the `Omega0_dr` cast cluster

**Files:**
- Modify: `species/base_species.h` (after `GetRadiationOmega0`, ~line 510)
- Modify: `species/composite_species.h` (after `GetRadiationOmega0`, ~line 54)
- Modify: `species/dark_radiation_species.h` (near `bi_rho_index`, ~line 122)
- Modify: `source/background_module.cpp:858-867`

- [ ] **Step 1: Add the `BaseSpecies` hook block**

In `species/base_species.h`, immediately after the `GetRadiationOmega0()` method
(the block ending at line 510), insert:

```cpp
  // ── NCDM/DR-family hooks (#308) ───────────────────────────────────────────
  // Neutral defaults let modules loop over all_species_ and call these directly
  // instead of downcasting; composites forward/sum over their children.

  /** Dark-radiation energy density today, read from the integrated background
   *  vector. Default 0; DarkRadiationSpecies returns its own; composites sum. */
  virtual double DarkRadiationRhoToday(const double* /*pvecback_integration*/) const {
    return 0.;
  }
```

- [ ] **Step 2: Add the `CompositeSpecies` sum-over-children override**

In `species/composite_species.h`, after the `GetRadiationOmega0()` override
(ends line 54), insert:

```cpp
  /** Sums DarkRadiationRhoToday() over all children. */
  double DarkRadiationRhoToday(const double* pvecback_integration) const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->DarkRadiationRhoToday(pvecback_integration);
    return sum;
  }
```

- [ ] **Step 3: Add the `DarkRadiationSpecies` override**

In `species/dark_radiation_species.h`, after `bi_rho_index()` (line 120-122),
insert:

```cpp
  double DarkRadiationRhoToday(const double* pvecback_integration) const override {
    return pvecback_integration[index_bi_rho_];
  }
```

- [ ] **Step 4: Dissolve the `Omega0_dr` loop in `background_module.cpp`**

Replace the current block (lines 858-867):

```cpp
  Omega0_dr_ = 0.;
  if (all_species_.count("DCDM_DR")) {
    auto& dcdm_dr  = dynamic_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
    Omega0_dcdm_   = pvecback_integration[dcdm_dr.dcdm().bi_rho_index()] / pba->H0 / pba->H0;
    Omega0_dr_    += pvecback_integration[dcdm_dr.dr().bi_rho_index()] / pba->H0 / pba->H0;
  }
  for (auto& [key, sp] : all_species_) {
    if (auto* dncdm_dr = dynamic_cast<DNCDM_DR_Species*>(sp.get()))
      Omega0_dr_ += pvecback_integration[dncdm_dr->dr().bi_rho_index()] / pba->H0 / pba->H0;
  }
```

with (DR is now polymorphic; the DCDM *matter* read stays by-name, out of scope):

```cpp
  Omega0_dr_ = 0.;
  for (auto& sp : all_species_)
    Omega0_dr_ += sp->DarkRadiationRhoToday(pvecback_integration) / pba->H0 / pba->H0;
  if (all_species_.count("DCDM_DR")) {
    auto& dcdm_dr = dynamic_cast<DCDM_DR_Species&>(*all_species_.at("DCDM_DR"));
    Omega0_dcdm_  = pvecback_integration[dcdm_dr.dcdm().bi_rho_index()] / pba->H0 / pba->H0;
  }
```

- [ ] **Step 5: Build**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly.

- [ ] **Step 6: Verify the DCDM_DR path is unchanged**

Run: `cd python && python -m pytest test_class.py -k "dcdm_dr" -v`
Expected: PASS (or SKIP if `classyref` absent). If skipped, additionally run the
compute check:
`cd python && python -c "from classy import Class; c=Class(); c.set({'Omega_dcdmdr':0.12,'Gamma_dcdm':10.0,'output':'tCl','l_max_scalars':200}); c.compute(); print('Omega_dr ok')"`
Expected: prints `Omega_dr ok`.

- [ ] **Step 7: Commit**

```bash
git add species/base_species.h species/composite_species.h species/dark_radiation_species.h source/background_module.cpp
git commit -m "$(cat <<'EOF'
species: DarkRadiationRhoToday hook replaces Omega0_dr dynamic_cast loop (#308)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: Background-side NCDM aggregates — `NeutrinoOmega0`, `NeffContribution`, `PrintNeffInfo`/`PrintMassInfo`

Dissolves `GetOmega0NcdmTot` and the N_eff / mass-info verbose loops. The
`GetNcdmSpecies(...).empty()` *gates* stay (NCDM detection, deferred); only the
inner contribution loops change. Behaviour-preserving (`GetNcdmSpecies` already
includes the DNCDM child; composite-sum-children matches).

**Files:**
- Modify: `species/base_species.h` (extend the hook block from Task 2)
- Modify: `species/composite_species.h` (forwarders)
- Modify: `species/ncdm_base_species.h` (decls / `override`)
- Modify: `source/background_module.cpp` (1531-1536, 503-508, 551-553)

- [ ] **Step 1: Add the `BaseSpecies` hooks**

In `species/base_species.h`, inside the `#308` hook block, add:

```cpp
  /** This species' contribution to Omega0 of the neutrino/NCDM sector (fnu).
   *  Default 0; NCDM returns GetOmega0(); composites sum children. */
  virtual double NeutrinoOmega0() const {
    return 0.;
  }

  /** Contribution to N_eff at redshift z. Default 0; NCDM returns GetNeff(z). */
  virtual double NeffContribution(double /*z*/) const {
    return 0.;
  }

  /** Verbose N_eff / mass info (background_verbose). Default no-op; NCDM prints. */
  virtual void PrintNeffInfo() const {}
  virtual void PrintMassInfo() const {}
```

- [ ] **Step 2: Add the `CompositeSpecies` forwarders**

In `species/composite_species.h`, after the Task-2 `DarkRadiationRhoToday`
override, add:

```cpp
  /** Sums NeutrinoOmega0() over children (matter NCDM children contribute). */
  double NeutrinoOmega0() const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->NeutrinoOmega0();
    return sum;
  }

  /** Sums NeffContribution() over children. */
  double NeffContribution(double z) const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->NeffContribution(z);
    return sum;
  }

  void PrintNeffInfo() const override {
    for (const auto& c : children_)
      c->PrintNeffInfo();
  }
  void PrintMassInfo() const override {
    for (const auto& c : children_)
      c->PrintMassInfo();
  }
```

- [ ] **Step 3: Override on `NCDMBaseSpecies`**

In `species/ncdm_base_species.h`, the lines currently reading (127-128):

```cpp
  void PrintNeffInfo() const;
  void PrintMassInfo() const;
```

become (add `override`, and the two new inline overrides):

```cpp
  void PrintNeffInfo() const override;
  void PrintMassInfo() const override;

  double NeutrinoOmega0() const override {
    return GetOmega0();
  }
  double NeffContribution(double z) const override {
    return GetNeff(z);
  }
```

- [ ] **Step 4: Dissolve `GetOmega0NcdmTot` (background_module.cpp:1531-1536)**

Replace:

```cpp
double BackgroundModule::GetOmega0NcdmTot() const {
  double total = 0.;
  for (auto* sp : GetNcdmSpecies(all_species_))
    total += sp->GetOmega0();
  return total;
}
```

with:

```cpp
double BackgroundModule::GetOmega0NcdmTot() const {
  double total = 0.;
  for (auto& sp : all_species_)
    total += sp->NeutrinoOmega0();
  return total;
}
```

- [ ] **Step 5: Dissolve the N_eff loop (background_module.cpp:503-508)**

Replace:

```cpp
      if (!GetNcdmSpecies(all_species_).empty()) {
        for (auto* sp : GetNcdmSpecies(all_species_)) {
          Neff += sp->GetNeff(0.);
          sp->PrintNeffInfo();
        }
      }
```

with (gate kept; inner loop polymorphic):

```cpp
      if (!GetNcdmSpecies(all_species_).empty()) {
        for (auto& sp : all_species_) {
          Neff += sp->NeffContribution(0.);
          sp->PrintNeffInfo();
        }
      }
```

- [ ] **Step 6: Dissolve the mass-info loop (background_module.cpp:551-553)**

Replace:

```cpp
  if ((pba->background_verbose > 0) && (!GetNcdmSpecies(all_species_).empty())) {
    for (auto* sp : GetNcdmSpecies(all_species_))
      sp->PrintMassInfo();
  }
```

with:

```cpp
  if ((pba->background_verbose > 0) && (!GetNcdmSpecies(all_species_).empty())) {
    for (auto& sp : all_species_)
      sp->PrintMassInfo();
  }
```

- [ ] **Step 7: Build**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly.

- [ ] **Step 8: Verify an NCDM scenario (N_eff, masses, fnu unchanged)**

Run: `cd python && python -m pytest test_class.py -k "ncdm or massive_ncdm" -v`
Expected: PASS / SKIP as applicable. Also a direct check:
`cd python && python -c "from classy import Class; c=Class(); c.set({'N_ur':0.0,'N_ncdm':1,'m_ncdm':0.06,'deg_ncdm':3.0,'output':'mPk','background_verbose':1}); c.compute(); print('Neff', c.Neff())"`
Expected: prints an `Neff` value (~3.04) with no errors; the verbose `-> total
N_eff` and mass-info lines print as before.

- [ ] **Step 9: Commit**

```bash
git add species/base_species.h species/composite_species.h species/ncdm_base_species.h source/background_module.cpp
git commit -m "$(cat <<'EOF'
species: polymorphic NeutrinoOmega0/NeffContribution/Print* replace NCDM casts (#308)

Dissolves GetOmega0NcdmTot and the N_eff / mass-info verbose loops; the
GetNcdmSpecies().empty() detection gates stay (deferred to #309/#310).

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Task 4: Tensor sites — `WriteTensorOutputColumnTitles` loop + `TensorMasslessRelativisticRho`

No behaviour change possible (tensors + DNCDM forbidden). The
`WriteTensorOutputColumnTitles` loop just swaps `NcdmFamily` for `all_species_`
(the base no-op covers non-NCDM and the composite); the `3·P` loop gets a hook.

**Files:**
- Modify: `species/base_species.h` (extend hook block)
- Modify: `species/composite_species.h` (forwarder)
- Modify: `species/ncdm_base_species.h` (override)
- Modify: `source/perturbations_module.cpp` (2734-2737, 5263-5271)

- [ ] **Step 1: Add the `BaseSpecies` hook**

In `species/base_species.h`, inside the `#308` block, add:

```cpp
  /** Contribution to the tensor-mode "relativistic" density in the massless
   *  approximation. Default 0; NCDM returns 3*P; composites sum. */
  virtual double TensorMasslessRelativisticRho(const double* /*pvecback*/) const {
    return 0.;
  }
```

- [ ] **Step 2: Add the `CompositeSpecies` forwarder**

In `species/composite_species.h`, add:

```cpp
  double TensorMasslessRelativisticRho(const double* pvecback) const override {
    double sum = 0.0;
    for (const auto& c : children_)
      sum += c->TensorMasslessRelativisticRho(pvecback);
    return sum;
  }
```

- [ ] **Step 3: Add the `NCDMBaseSpecies` override**

In `species/ncdm_base_species.h`, alongside the other Task-3 overrides, add:

```cpp
  double TensorMasslessRelativisticRho(const double* pvecback) const override {
    return 3. * P(pvecback);
  }
```

- [ ] **Step 4: Dissolve `WriteTensorOutputColumnTitles` (perturbations_module.cpp:2734-2737)**

Replace:

```cpp
      if (evolve_tensor_ncdm_ == _TRUE_) {
        for (auto* sp : NcdmFamily(all_species_))
          sp->WriteTensorOutputColumnTitles(tensor_titles_);
      }
```

with:

```cpp
      if (evolve_tensor_ncdm_ == _TRUE_) {
        for (auto& sp : all_species_)
          sp->WriteTensorOutputColumnTitles(tensor_titles_);
      }
```

- [ ] **Step 5: Dissolve the tensor-massless `3·P` loop (perturbations_module.cpp:5263-5271)**

Replace:

```cpp
        if (HasNcdm(all_species_)) {
          for (auto& [name, sp] : all_species_) {
            auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
            if (!ncdm_sp)
              continue;
            /* (3 p_ncdm1) is the "relativistic" contribution to rho_ncdm1 */
            rho_relativistic += 3. * ncdm_sp->P(ppw->pvecback);
          }
        }
```

with (HasNcdm gate kept; inner loop polymorphic):

```cpp
        if (HasNcdm(all_species_)) {
          for (auto& sp : all_species_)
            rho_relativistic += sp->TensorMasslessRelativisticRho(ppw->pvecback);
        }
```

- [ ] **Step 6: Build**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly.

- [ ] **Step 7: Verify the tensor massive-NCDM regression**

Run: `cd python && python -m pytest test_class.py -k "tensor_massive_ncdm" -v`
Expected: PASS (or SKIP without `classyref`). Direct compute check:
`cd python && python -c "from classy import Class; c=Class(); c.set({'N_ur':0.0,'N_ncdm':1,'m_ncdm':0.06,'deg_ncdm':3.0,'output':'tCl','modes':'t','tensor method':'massless'}); c.compute(); print('tensor ncdm ok')"`
Expected: prints `tensor ncdm ok`.

- [ ] **Step 8: Commit**

```bash
git add species/base_species.h species/composite_species.h species/ncdm_base_species.h source/perturbations_module.cpp
git commit -m "$(cat <<'EOF'
perturbations: TensorMasslessRelativisticRho + all_species_ tensor-titles loop (#308)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Task 5: Relativistic-IC check — `CheckUltraRelativisticAtIc` / `IsUltraRelativisticAtIc` (the behaviour change)

The old narrow `dynamic_cast<NCDMSpecies*>` excluded the wrapped DNCDM child;
forwarding now includes it (the composite forwards to its `dncdm` child). This is
the one intentional behaviour change. A new DNCDM_DR compute test guards it.

**Files:**
- Modify: `species/base_species.h` (extend hook block)
- Modify: `species/composite_species.h` (forwarders)
- Modify: `species/ncdm_base_species.h` (decls)
- Modify: `species/ncdm_base_species.cpp` (bodies)
- Modify: `source/perturbations_module.cpp` (2466-2481, 2501-2510)
- Modify: `python/test_class.py` (new DNCDM test)

- [ ] **Step 1: Write the failing/guard test first**

In `python/test_class.py`, add this method to the reference-style test class that
contains `test_dcdm_dr_matches_reference` (the class around line 755):

```python
    def test_dncdm_dr_computes(self):
        """#308 guard: a decaying-NCDM (DNCDM_DR) run still computes after the
        relativistic-IC check is forwarded to the wrapped DNCDM child."""
        cosmo = Class()
        try:
            cosmo.set({
                'output': 'tCl',
                'l_max_scalars': 200,
                'N_ur': 3.046,
                'omega_b': 0.022032,
                'omega_cdm': 0.12038,
                # Pin YHe to skip BBN; this is a compute guard, not a precision test.
                'YHe': 0.25,
                'dncdm1.type': 'ncdm_decay_dr',
                'dncdm1.m': 1.0,
                'dncdm1.T': 0.71611,
                'dncdm1.Gamma': 1e3,
                'dncdm1.Omega_ini': 0.001,  # ncdm_decay_dr rejects .Omega; use .Omega_ini
            })
            cosmo.compute()
            self.assertTrue(cosmo.raw_cl(100)['tt'].size > 0)
        finally:
            cosmo.struct_cleanup()
            cosmo.empty()
```

- [ ] **Step 2: Run the test on the current (pre-change) code to confirm it passes**

Run: `cd python && python -m pytest test_class.py -k "dncdm_dr_computes" -v`
Expected: PASS (the DNCDM run computes today; this establishes the baseline the
behaviour change must preserve).

- [ ] **Step 3: Add the `BaseSpecies` hooks**

In `species/base_species.h`, inside the `#308` block, add:

```cpp
  /** Initial-time ultra-relativistic check (throws via class_test on failure).
   *  Default no-op; NCDM checks |w-1/3|; composites forward to children. */
  virtual void CheckUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const {}

  /** Initial-time ultra-relativistic predicate. Default true; NCDM checks
   *  |w-1/3| <= tol; composites AND over children. */
  virtual bool IsUltraRelativisticAtIc(const double* /*pvecback*/, double /*tol*/) const {
    return true;
  }
```

- [ ] **Step 4: Add the `CompositeSpecies` forwarders**

In `species/composite_species.h`, add:

```cpp
  void CheckUltraRelativisticAtIc(const double* pvecback, double tol) const override {
    for (const auto& c : children_)
      c->CheckUltraRelativisticAtIc(pvecback, tol);
  }
  bool IsUltraRelativisticAtIc(const double* pvecback, double tol) const override {
    for (const auto& c : children_)
      if (!c->IsUltraRelativisticAtIc(pvecback, tol))
        return false;
    return true;
  }
```

- [ ] **Step 5: Declare the overrides on `NCDMBaseSpecies`**

In `species/ncdm_base_species.h`, alongside the other overrides, add:

```cpp
  void CheckUltraRelativisticAtIc(const double* pvecback, double tol) const override;
  bool IsUltraRelativisticAtIc(const double* pvecback, double tol) const override;
```

- [ ] **Step 6: Implement the bodies in `ncdm_base_species.cpp`**

Add (near the other small accessors; `class_test`, `fabs`, `<cmath>` are already
available in this file):

```cpp
void NCDMBaseSpecies::CheckUltraRelativisticAtIc(const double* pvecback, double tol) const {
  const double p_ncdm   = P(pvecback);
  const double rho_ncdm = Rho(pvecback);
  class_test(fabs(p_ncdm / rho_ncdm - 1. / 3.) > tol,
             "your choice of initial time for integrating wavenumbers is inappropriate: it "
             "corresponds to a time at which the ncdm species '%s' is not "
             "ultra-relativistic anymore, with w=%g, p=%g and rho=%g\n",
             name().c_str(),
             p_ncdm / rho_ncdm,
             p_ncdm,
             rho_ncdm);
}

bool NCDMBaseSpecies::IsUltraRelativisticAtIc(const double* pvecback, double tol) const {
  return fabs(P(pvecback) / Rho(pvecback) - 1. / 3.) <= tol;
}
```

- [ ] **Step 7: Dissolve the hard check (perturbations_module.cpp:2466-2481)**

Replace:

```cpp
  if (HasNcdm(all_species_)) {
    for (auto& [name, sp] : all_species_) {
      auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
      if (!ncdm_sp)
        continue;
      const double p_ncdm   = ncdm_sp->P(ppw->pvecback);
      const double rho_ncdm = ncdm_sp->Rho(ppw->pvecback);
      class_test(fabs(p_ncdm / rho_ncdm - 1. / 3.) > ppr->tol_ncdm_initial_w,
                 "your choice of initial time for integrating wavenumbers is inappropriate: it "
                 "corresponds to a time at which the ncdm species '%s' is not "
                 "ultra-relativistic anymore, with w=%g, p=%g and rho=%g\n",
                 ncdm_sp->name().c_str(),
                 p_ncdm / rho_ncdm,
                 p_ncdm,
                 rho_ncdm);
    }
  }
```

with:

```cpp
  if (HasNcdm(all_species_)) {
    for (auto& sp : all_species_)
      sp->CheckUltraRelativisticAtIc(ppw->pvecback, ppr->tol_ncdm_initial_w);
  }
```

- [ ] **Step 8: Dissolve the soft check (perturbations_module.cpp:2501-2510)**

Replace:

```cpp
    if (HasNcdm(all_species_)) {
      for (auto& [name, sp] : all_species_) {
        auto* ncdm_sp = dynamic_cast<NCDMSpecies*>(sp.get());
        if (!ncdm_sp)
          continue;
        if (fabs(ncdm_sp->P(ppw->pvecback) / ncdm_sp->Rho(ppw->pvecback) - 1. / 3.) >
            ppr->tol_ncdm_initial_w)
          is_early_enough = _FALSE_;
      }
    }
```

with:

```cpp
    if (HasNcdm(all_species_)) {
      for (auto& sp : all_species_) {
        if (!sp->IsUltraRelativisticAtIc(ppw->pvecback, ppr->tol_ncdm_initial_w))
          is_early_enough = _FALSE_;
      }
    }
```

- [ ] **Step 9: Build**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly.

- [ ] **Step 10: Run the DNCDM guard test (the decision instrument)**

Run: `cd python && python -m pytest test_class.py -k "dncdm_dr_computes" -v`
Expected: PASS — the DNCDM run still computes; the decaying child is now included
in the check but is ultra-relativistic at `tau_ini`, so the assertion holds.

> **If this test FAILS** with the "not ultra-relativistic anymore" message, the
> behaviour change is harmful for DNCDM. Fall back: in Steps 7-8 keep the loop but
> skip the composite — e.g. guard with `if (sp->AsNcdm-equivalent)`... instead,
> simplest fallback that preserves old behaviour: leave these two sites calling the
> old `dynamic_cast<NCDMSpecies*>` loop (revert Steps 7-8 only), keep the new
> methods unused here, and note in the PR that the relativistic-IC check stays
> top-level-only. Re-run the test to confirm green, then commit.

- [ ] **Step 11: Verify ordinary NCDM still passes**

Run: `cd python && python -m pytest test_class.py -k "ncdm or massive_ncdm" -v`
Expected: PASS / SKIP as applicable; no new failures.

- [ ] **Step 12: Commit**

```bash
git add species/base_species.h species/composite_species.h species/ncdm_base_species.h species/ncdm_base_species.cpp source/perturbations_module.cpp python/test_class.py
git commit -m "$(cat <<'EOF'
perturbations: forward relativistic-IC check to species; cover DNCDM child (#308)

Replaces the narrow dynamic_cast<NCDMSpecies*> IC checks with forwarded
Check/IsUltraRelativisticAtIc hooks. The wrapped DNCDM child is now included
(via composite forwarding); a new test_dncdm_dr_computes guards this.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Task 6: Halofit WDM warning — `WarnIfTooHeavyForHalofit`

**Files:**
- Modify: `species/base_species.h` (extend hook block)
- Modify: `species/composite_species.h` (forwarder)
- Modify: `species/ncdm_base_species.h` (decl)
- Modify: `species/ncdm_base_species.cpp` (body)
- Modify: `source/nonlinear_module.cpp` (917-928)

- [ ] **Step 1: Add the `BaseSpecies` hook**

In `species/base_species.h`, inside the `#308` block, add:

```cpp
  /** Warn (stdout) if this species is too heavy for Halofit/HMcode. The policy
   *  threshold is passed by the nonlinear module. Default no-op; NCDM compares
   *  its mass; composites forward. */
  virtual void WarnIfTooHeavyForHalofit(double /*m_ev_threshold*/) const {}
```

- [ ] **Step 2: Add the `CompositeSpecies` forwarder**

In `species/composite_species.h`, add:

```cpp
  void WarnIfTooHeavyForHalofit(double m_ev_threshold) const override {
    for (const auto& c : children_)
      c->WarnIfTooHeavyForHalofit(m_ev_threshold);
  }
```

- [ ] **Step 3: Declare + implement on `NCDMBaseSpecies`**

In `species/ncdm_base_species.h`, add the declaration:

```cpp
  void WarnIfTooHeavyForHalofit(double m_ev_threshold) const override;
```

In `species/ncdm_base_species.cpp`, add the body:

```cpp
void NCDMBaseSpecies::WarnIfTooHeavyForHalofit(double m_ev_threshold) const {
  const double m_ncdm_in_electronvolt = GetMassInElectronvolt();
  if (m_ncdm_in_electronvolt > m_ev_threshold)
    fprintf(stdout,
            "Warning: Halofit and HMcode are proved to work for CDM, and also with a small "
            "HDM component. But it sounds like you are running with a WDM component of mass "
            "%f eV, which makes the use of Halofit suspicious.\n",
            m_ncdm_in_electronvolt);
}
```

- [ ] **Step 4: Dissolve the warning loop (nonlinear_module.cpp:917-928)**

Replace:

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

with:

```cpp
    for (auto& sp : all_species_)
      sp->WarnIfTooHeavyForHalofit(_M_EV_TOO_BIG_FOR_HALOFIT_);
```

(Leave the `IDM_DR_IDR` `has_idm_dr()` block immediately below unchanged — IDM is
out of scope.)

- [ ] **Step 5: Build**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly.

- [ ] **Step 6: Verify a heavy-NCDM + nonlinear run computes (and warns)**

Run: `cd python && python -c "from classy import Class; c=Class(); c.set({'N_ur':2.0308,'N_ncdm':1,'m_ncdm':2.0,'output':'mPk','non linear':'halofit'}); c.compute(); print('halofit ncdm ok')"`
Expected: prints the WDM "Warning:" line on stdout, then `halofit ncdm ok`.

- [ ] **Step 7: Commit**

```bash
git add species/base_species.h species/composite_species.h species/ncdm_base_species.h species/ncdm_base_species.cpp source/nonlinear_module.cpp
git commit -m "$(cat <<'EOF'
nonlinear: WarnIfTooHeavyForHalofit hook replaces NCDM dynamic_cast loop (#308)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Task 7: Full-suite verification

**Files:** none (verification only).

- [ ] **Step 1: Confirm no stray `dynamic_cast`/helper regressions remain in scope**

Run: `grep -n "dynamic_cast" source/perturbations_module.cpp source/background_module.cpp source/nonlinear_module.cpp source/output_module.cpp`
Expected: the remaining casts are only the deferred ones — `HasNcdm`/`NcdmFamily`
(perturbations 53/65/76/78), `DrSpeciesCount`, the `SetSourceSlot` loops (642-664),
the nested-layout tensor sites (5285-5293, 5918-5934), `GetNcdmSpecies`
(background 98/100), the background index-reg loop (636-651), `Omega0_dcdm_` by-name
(860), and `output_module.cpp:668`. No casts remain at the sites dissolved in Tasks
2-6.

- [ ] **Step 2: Build clean from scratch**

Run: `pip install . --no-build-isolation`
Expected: builds cleanly, no warnings about the new methods.

- [ ] **Step 3: Run the full scenario suite**

Run: `cd python && TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py`
Expected: all PASS (same set as on `master`). If `classyref` is installed, also run
the reference comparison:
`cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest -v -m test_scenario test_class.py`
Expected: all PASS.

- [ ] **Step 4: Run the named regression/guard tests together**

Run: `cd python && python -m pytest test_class.py -k "dcdm_dr or tensor_massive_ncdm or dncdm_dr_computes" -v`
Expected: PASS (reference ones SKIP without `classyref`; `dncdm_dr_computes` PASS).

- [ ] **Step 5: Update the spec status**

In `docs/superpowers/specs/2026-06-16-issue-308-remove-species-dynamic-cast-dispatch-design.md`,
change the Status line to `Implemented`.

- [ ] **Step 6: Commit**

```bash
git add docs/superpowers/specs/2026-06-16-issue-308-remove-species-dynamic-cast-dispatch-design.md
git commit -m "$(cat <<'EOF'
docs: mark #308 spec implemented

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
EOF
)"
```

---

## Self-review notes

- **Spec coverage:** Bucket 2 deletion → Task 1. Bucket 1 dissolves: Omega0_dr →
  Task 2; GetOmega0NcdmTot + N_eff/mass-info → Task 3; tensor titles + 3·P → Task 4;
  relativistic-IC check → Task 5; Halofit warning → Task 6. Deferred items
  (HasNcdm, NcdmFamily/GetNcdmSpecies counts/gates, DrSpeciesCount, SetSourceSlot,
  nested-layout tensor, Omega0_dcdm by-name, IDM) explicitly left, verified in
  Task 7 Step 1. Behaviour change (relativistic-IC) carries its own guard test
  (Task 5) and fallback.
- **Type consistency:** every new method name is used identically across base decl,
  composite forwarder, NCDM/DR override, and call site: `DarkRadiationRhoToday`,
  `NeutrinoOmega0`, `NeffContribution`, `PrintNeffInfo`, `PrintMassInfo`,
  `TensorMasslessRelativisticRho`, `CheckUltraRelativisticAtIc`,
  `IsUltraRelativisticAtIc`, `WarnIfTooHeavyForHalofit`.
- **No placeholders:** every code/edit step shows the exact before/after text and
  every test step an exact command + expected result.
