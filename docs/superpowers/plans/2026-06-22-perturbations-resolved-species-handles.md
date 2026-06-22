# Perturbations Resolved-Species-Handles Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate every repeated (per-step / per-sample / per-k-bisection) string-keyed `SpeciesCollection` lookup from the perturbations hot path by resolving each once at construction, capturing the prototype's ~1–3% Perturbations-stage win while keeping output byte-identical.

**Architecture:** Introduce one private `ResolvedSpecies resolved_` struct on `PerturbationsModule`, populated once by `ResolveSpecies()` in the constructor (the `all_species_` collection is frozen after construction, so every lookup result is final). Replace per-step `all_species_.count/at/find/index_of(...)` and `HasNcdm(...)` with reads of `resolved_`; typed handles also eliminate the per-step `static_cast` downcasts. Add a cached `perturb_vector::cdm_layout` (mirroring the existing `photon_layout`/`baryon_layout`). The existing scattered `idr_nature_`/`ppf_fluid_` caches are consolidated into `resolved_`. No new virtual hooks; the approximation framework stays module-owned.

**Tech Stack:** C++17, CMake (build tree `build_prof/`, binaries land in repo root), the existing `class`/`class_profiled` executables, `python/test_class.py` regression suite, `.superpowers/sdd/tol_check.py` physics gate.

## Global Constraints

- **Byte-identical output is expected** (a cached value equals the scan result for a frozen collection; no floating-point reassociation is introduced). Every code task must verify `cl_lensed.dat` + `pk.dat` are **byte-identical** to the master baseline. If any task is not byte-identical, STOP and investigate — it indicates a logic error, not benign FP drift.
- **Physics gate** (final): `.superpowers/sdd/tol_check.py CL_OUT CL_REF PK_OUT PK_REF` at default tol `1e-3` (0.1% on Cl^TT + P(k), no TE/EE/BB zero-crossing max-rel), per `feedback_no_bit_identical_requirement`.
- **Final correctness gate**: `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` full suite, 0 failures vs master classyref.
- **Never `git add -A`** in this repo (build artifacts leak in); always `git add` explicit paths.
- **clang-format 22.1.3** every touched file before the final commit.
- **No new species-type picking in the generic loops** and **no new typed accessors on `SpeciesCollection`** (`feedback_no_species_picking_in_modules`); typed handles only hoist downcasts that already exist at these module-owned sites to construction time.
- Benchmark rigor (Task 6): `OMP_NUM_THREADS=1`, fresh isolated worktrees, identical flags, order-alternating paired A/B, min/median + bootstrap CI, `-falign-functions=64` layout probe.
- Line numbers below are as of master `eba6f7d7`/branch base and **will shift** as edits land — match on the quoted code, not the line number.

---

### Task 0: Capture the master golden baseline

**Files:** none modified (baseline capture only).

**Interfaces:**
- Produces: golden reference files in `/tmp/perturb_handles_golden/` used by every later task's byte-identity check. The current branch HEAD (spec + plan commits only) is code-identical to master, so building it produces the master output.

- [ ] **Step 1: Build the current tree (sanity + baseline binary)**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: builds `class` and `class_profiled` into the repo root with no errors (last line a link step / `Built target class_profiled`).

- [ ] **Step 2: Produce baseline output and stash it as golden**

Run:
```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
mkdir -p /tmp/perturb_handles_golden
cp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat
cp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat        /tmp/perturb_handles_golden/pk.dat
cp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl.dat        /tmp/perturb_handles_golden/cl.dat
ls -l /tmp/perturb_handles_golden/
```
Expected: three non-empty `.dat` files copied.

- [ ] **Step 3: Record the baseline Perturbations timing (for Task 2's checkpoint)**

Run:
```bash
OMP_NUM_THREADS=1 ./class_profiled base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini 8 | grep -i Perturbations
```
Expected: a line with min/median/mean/stddev ms for the Perturbations stage. Note the median for later comparison (informal — the rigorous A/B is Task 6).

- [ ] **Step 4: No commit** (nothing changed).

---

### Task 1: Introduce `ResolvedSpecies` + `ResolveSpecies()`, consolidate `idr_nature_`/`ppf_fluid_`

**Files:**
- Modify: `source/perturbations_module.h` (forward decls, struct, member, method decl, `ppf_fluid()` body; remove `idr_nature_` and `ppf_fluid_` members)
- Modify: `source/perturbations_module.cpp` (ctor blocks → `ResolveSpecies()` call + definition; repoint `idr_nature_` usage)

**Interfaces:**
- Produces, for all later tasks, the member `ResolvedSpecies resolved_` with fields:
  - `const IDM_DR_IDR_Species* idm_dr_idr` (nullptr if absent)
  - `std::size_t idm_dr_idr_index` (`all_species_.size()` if absent)
  - `const IDM_DRMD_IDR_DRMD_Species* idm_drmd` (nullptr if absent)
  - `FluidSpecies* ppf_fluid` (nullptr unless a PPF fluid is present)
  - `const BaseSpecies* ur` (nullptr if absent)
  - `const BaseSpecies* lambda` (nullptr if absent)
  - `bool has_ncdm`
  - `int idr_nature` (default `idr_free_streaming`)
  - `std::size_t cdm_index` (`all_species_.size()` if absent)
- Accessor `ppf_fluid()` now returns `resolved_.ppf_fluid` (call sites unchanged).

- [ ] **Step 1: Add forward declarations and the struct/member to the header**

In `source/perturbations_module.h`, replace the lone forward declaration:
```cpp
class FluidSpecies;
```
with:
```cpp
class FluidSpecies;
class IDM_DR_IDR_Species;
class IDM_DRMD_IDR_DRMD_Species;
```

- [ ] **Step 2: Point the `ppf_fluid()` accessor at the struct**

In `source/perturbations_module.h`, change the accessor body:
```cpp
  FluidSpecies* ppf_fluid() const {
    return ppf_fluid_;
  }
```
to:
```cpp
  FluidSpecies* ppf_fluid() const {
    return resolved_.ppf_fluid;
  }
```

- [ ] **Step 3: Replace the standalone `idr_nature_` member with the struct + declare `ResolveSpecies()`**

In `source/perturbations_module.h`, replace:
```cpp
  int idr_nature_ =
      idr_free_streaming; /**< cached IDR nature (free-streaming or fluid); set once in constructor */
```
with:
```cpp
  /* Resolved-once views into the frozen all_species_ collection. The module owns
     several numerical approximation schemes (TCA/RSA/UFA/ncdmfa), the
     Einstein/metric equations, and output conventions that legitimately reference
     specific named species; these handles give that module-owned logic O(1)
     access to the relevant species' data, replacing per-step/per-sample/per-k
     O(n) string-keyed scans of all_species_ (frozen after construction) and the
     per-step downcasts that accompanied them. Resolved by ResolveSpecies() in the
     ctor. */
  struct ResolvedSpecies {
    const IDM_DR_IDR_Species* idm_dr_idr = nullptr;            // TCA-idm_dr, RSA-idr
    std::size_t idm_dr_idr_index = static_cast<std::size_t>(-1);// its species_layouts[] index
    const IDM_DRMD_IDR_DRMD_Species* idm_drmd = nullptr;       // TCA-idm_drmd
    FluidSpecies* ppf_fluid = nullptr;                         // PPF fluid (absorbs ppf_fluid_)
    const BaseSpecies* ur = nullptr;                           // UFA, RSA-ur, tensor rho
    const BaseSpecies* lambda = nullptr;                       // excluded from rho_tot (delta_tot)
    bool has_ncdm = false;                                     // ncdmfa scheme present
    int idr_nature = idr_free_streaming;                       // absorbs idr_nature_
    std::size_t cdm_index = static_cast<std::size_t>(-1);      // h source (sync gauge)
  };
  ResolvedSpecies resolved_;

  void ResolveSpecies(); /**< resolve `resolved_` once; called from the ctor */
```

- [ ] **Step 4: Remove the old `ppf_fluid_` member from the header**

In `source/perturbations_module.h`, delete the line:
```cpp
  FluidSpecies* ppf_fluid_ = nullptr;
```

- [ ] **Step 5: Replace the ctor resolution blocks with a `ResolveSpecies()` call**

In `source/perturbations_module.cpp`, in the constructor, replace:
```cpp
  /* Cache idr_nature once (avoids per-RHS string-map lookup). */
  if (all_species_.count("IDM_DR_IDR"))
    idr_nature_ =
        static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).idr().idr_nature();

  /* Cache the single PPF fluid (the one downcast, localized here). use_ppf is a
     construction-fixed input; at most one "Fluid" can exist. */
  if (auto* p = all_species_.find("Fluid")) {
    auto* f = static_cast<FluidSpecies*>(p->get());
    if (f->use_ppf())
      ppf_fluid_ = f;
  }
  perturb_init();
```
with:
```cpp
  ResolveSpecies();
  perturb_init();
```

- [ ] **Step 6: Add the `ResolveSpecies()` definition**

In `source/perturbations_module.cpp`, immediately after the constructor's closing brace (before `PerturbationsModule::~PerturbationsModule()`), add:
```cpp
void PerturbationsModule::ResolveSpecies() {
  if (auto* p = all_species_.find("IDM_DR_IDR")) {
    resolved_.idm_dr_idr       = static_cast<const IDM_DR_IDR_Species*>(p->get());
    resolved_.idm_dr_idr_index = all_species_.index_of("IDM_DR_IDR");
    resolved_.idr_nature       = resolved_.idm_dr_idr->idr().idr_nature();
  }
  if (auto* p = all_species_.find("IDM_DRMD_IDR_DRMD"))
    resolved_.idm_drmd = static_cast<const IDM_DRMD_IDR_DRMD_Species*>(p->get());

  if (auto* p = all_species_.find("UR"))
    resolved_.ur = p->get();
  if (auto* p = all_species_.find("Lambda"))
    resolved_.lambda = p->get();

  resolved_.has_ncdm  = HasNcdm(all_species_);
  resolved_.cdm_index = all_species_.index_of("CDM");

  if (auto* p = all_species_.find("Fluid")) {
    auto* f = static_cast<FluidSpecies*>(p->get());
    if (f->use_ppf())
      resolved_.ppf_fluid = f;
  }
}
```

- [ ] **Step 7: Repoint the remaining `idr_nature_` usage**

In `source/perturbations_module.cpp` (in `perturb_total_stress_energy`, ~line 4817), replace:
```cpp
    ppw->scalar_ctx.idr_nature = idr_nature_;
```
with:
```cpp
    ppw->scalar_ctx.idr_nature = resolved_.idr_nature;
```

- [ ] **Step 8: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: builds cleanly, no errors. (`resolved_`'s handle fields are populated but not yet read — that is expected; struct members do not warn when unused.)

- [ ] **Step 9: Run base_2018 and verify byte-identical**

Run:
```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  && cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat /tmp/perturb_handles_golden/pk.dat \
  && echo "BYTE-IDENTICAL OK"
```
Expected: `BYTE-IDENTICAL OK` (no `cmp` differ output).

- [ ] **Step 10: Commit**

```bash
git add source/perturbations_module.h source/perturbations_module.cpp
git commit -m "perturb: add ResolvedSpecies + ResolveSpecies(); consolidate idr_nature_/ppf_fluid_"
```

---

### Task 2: Site 1 + 1b — `perturb_einstein` IDR corrections (perf-critical)

**Files:**
- Modify: `source/perturbations_module.cpp` (`perturb_einstein`, the three IDM_DR_IDR blocks)

**Interfaces:**
- Consumes: `resolved_.idm_dr_idr` (`const IDM_DR_IDR_Species*`), `resolved_.idm_dr_idr_index` (`std::size_t`) from Task 1.

- [ ] **Step 1: Replace the Newtonian-gauge RSA-IDR guard**

In `perturb_einstein`, replace:
```cpp
      if ((all_species_.count("IDM_DR_IDR")) &&
          (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)) {
        perturb_rsa_idr_delta_and_theta(k, y, a_prime_over_a, ppw->pvecthermo.data(), ppw);
      }
```
with (this is the first occurrence, in the `newtonian` branch):
```cpp
      if (resolved_.idm_dr_idr &&
          (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)) {
        perturb_rsa_idr_delta_and_theta(k, y, a_prime_over_a, ppw->pvecthermo.data(), ppw);
      }
```

- [ ] **Step 2: Replace the synchronous-gauge RSA-IDR guard + theta correction**

In `perturb_einstein` (the `synchronous` branch), replace:
```cpp
      if ((all_species_.count("IDM_DR_IDR")) &&
          (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)) {
        perturb_rsa_idr_delta_and_theta(k, y, a_prime_over_a, ppw->pvecthermo.data(), ppw);

        ppw->rho_plus_p_theta += 4. / 3. *
                                 static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"))
                                     .idr()
                                     .Rho(ppw->pvecback.data()) *
                                 ppw->rsa_theta_idr;
      }
```
with:
```cpp
      if (resolved_.idm_dr_idr &&
          (ppw->approx[ppw->index_ap_rsa_idr] == (int) rsa_idr_on)) {
        perturb_rsa_idr_delta_and_theta(k, y, a_prime_over_a, ppw->pvecthermo.data(), ppw);

        ppw->rho_plus_p_theta += 4. / 3. *
                                 resolved_.idm_dr_idr->idr().Rho(ppw->pvecback.data()) *
                                 ppw->rsa_theta_idr;
      }
```

- [ ] **Step 3: Replace the TCA-IDR shear guard + layout index (Site 1b)**

In `perturb_einstein`, replace:
```cpp
      if ((all_species_.count("IDM_DR_IDR")) &&
          (ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_on)) {
        const size_t idr_th_i  = all_species_.index_of("IDM_DR_IDR");
        const auto& idr_th_lay = static_cast<const IDM_DR_IDR_Species::PerturbLayout&>(
                                     *ppw->pv->species_layouts[idr_th_i])
                                     .idr;
        shear_idr = 0.5 * 8. / 15. / ppw->pvecthermo[thermodynamics_module_->index_th_dmu_idm_dr_] /
                    static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"))
                        .idr()
                        .alpha_idm_dr()[0] *
                    (y[idr_th_lay.idx_theta] + k2 * ppw->pvecmetric[ppw->index_mt_alpha]);

        ppw->rho_plus_p_shear += 4. / 3. *
                                 static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"))
                                     .idr()
                                     .Rho(ppw->pvecback.data()) *
                                 shear_idr;
      }
```
with:
```cpp
      if (resolved_.idm_dr_idr &&
          (ppw->approx[ppw->index_ap_tca_idm_dr] == (int) tca_idm_dr_on)) {
        const auto& idr_th_lay = static_cast<const IDM_DR_IDR_Species::PerturbLayout&>(
                                     *ppw->pv->species_layouts[resolved_.idm_dr_idr_index])
                                     .idr;
        shear_idr = 0.5 * 8. / 15. / ppw->pvecthermo[thermodynamics_module_->index_th_dmu_idm_dr_] /
                    resolved_.idm_dr_idr->idr().alpha_idm_dr()[0] *
                    (y[idr_th_lay.idx_theta] + k2 * ppw->pvecmetric[ppw->index_mt_alpha]);

        ppw->rho_plus_p_shear += 4. / 3. *
                                 resolved_.idm_dr_idr->idr().Rho(ppw->pvecback.data()) *
                                 shear_idr;
      }
```

- [ ] **Step 4: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: clean build.

- [ ] **Step 5: Verify byte-identical on base_2018 (the IDM_DR_IDR-absent common path)**

Run:
```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  && cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat /tmp/perturb_handles_golden/pk.dat \
  && echo "BYTE-IDENTICAL OK"
```
Expected: `BYTE-IDENTICAL OK`. (base_2018 has no IDM_DR_IDR, so this confirms the perf-critical absent path; the IDM_DR_IDR-present path is gated by the full TEST_LEVEL=2 suite in Task 6.)

- [ ] **Step 6: Perf checkpoint (informal) — confirm the win lands**

Run:
```bash
OMP_NUM_THREADS=1 ./class_profiled base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini 8 | grep -i Perturbations
```
Expected: Perturbations median at or below the Task 0 baseline (the +0.85% einstein-IDR no-op-scan removal). This is informal; the rigorous paired A/B is Task 6. Do not block on a small number here — byte-identity is the gate.

- [ ] **Step 7: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturb: resolve IDM_DR_IDR handle in perturb_einstein (Site 1/1b)"
```

---

### Task 3: Site 2 + Site 3 — `perturb_rsa_delta_and_theta` UR and tensor relativistic-rho

**Files:**
- Modify: `source/perturbations_module.cpp` (`perturb_rsa_delta_and_theta` UR blocks; `perturb_total_stress_energy` tensor branch)

**Interfaces:**
- Consumes: `resolved_.ur` (`const BaseSpecies*`), `resolved_.has_ncdm` (`bool`) from Task 1.

- [ ] **Step 1: Replace the Newtonian-gauge UR streaming guard**

In `perturb_rsa_delta_and_theta`, replace the first:
```cpp
    if (all_species_.count("UR")) {
      if (ppr->radiation_streaming_approximation == static_cast<int>(rsa_method::rsa_null)) {
        ppw->rsa_delta_ur = 0.;
        ppw->rsa_theta_ur = 0.;
      }
      else {
        ppw->rsa_delta_ur = -4. * y[ppw->pv->index_pt_phi];
        ppw->rsa_theta_ur = 6. * ppw->pvecmetric[ppw->index_mt_phi_prime];
      }
    }
```
with (only the guard line changes):
```cpp
    if (resolved_.ur) {
      if (ppr->radiation_streaming_approximation == static_cast<int>(rsa_method::rsa_null)) {
        ppw->rsa_delta_ur = 0.;
        ppw->rsa_theta_ur = 0.;
      }
      else {
        ppw->rsa_delta_ur = -4. * y[ppw->pv->index_pt_phi];
        ppw->rsa_theta_ur = 6. * ppw->pvecmetric[ppw->index_mt_phi_prime];
      }
    }
```

- [ ] **Step 2: Replace the synchronous-gauge UR streaming guard**

In `perturb_rsa_delta_and_theta`, replace:
```cpp
    if (all_species_.count("UR")) {
      if (ppr->radiation_streaming_approximation == static_cast<int>(rsa_method::rsa_null)) {
        ppw->rsa_delta_ur = 0.;
        ppw->rsa_theta_ur = 0.;
      }
      else {
        ppw->rsa_delta_ur = 4. / k2 *
                            (a_prime_over_a * ppw->pvecmetric[ppw->index_mt_h_prime] -
                             k2 * y[ppw->pv->index_pt_eta]);
        ppw->rsa_theta_ur = -0.5 * ppw->pvecmetric[ppw->index_mt_h_prime];
      }
    }
```
with (only the guard line changes):
```cpp
    if (resolved_.ur) {
      if (ppr->radiation_streaming_approximation == static_cast<int>(rsa_method::rsa_null)) {
        ppw->rsa_delta_ur = 0.;
        ppw->rsa_theta_ur = 0.;
      }
      else {
        ppw->rsa_delta_ur = 4. / k2 *
                            (a_prime_over_a * ppw->pvecmetric[ppw->index_mt_h_prime] -
                             k2 * y[ppw->pv->index_pt_eta]);
        ppw->rsa_theta_ur = -0.5 * ppw->pvecmetric[ppw->index_mt_h_prime];
      }
    }
```

- [ ] **Step 3: Replace the UR total-contribution block (collapses 3 scans + 3 Rho calls)**

In `perturb_rsa_delta_and_theta`, replace:
```cpp
  if (all_species_.count("UR")) {
    ppw->delta_rho += all_species_.at("UR")->Rho(ppw->pvecback.data()) * ppw->rsa_delta_ur;
    ppw->delta_p += 1. / 3. * all_species_.at("UR")->Rho(ppw->pvecback.data()) * ppw->rsa_delta_ur;
    ppw->rho_plus_p_theta += 4. / 3. * all_species_.at("UR")->Rho(ppw->pvecback.data()) *
                             ppw->rsa_theta_ur;
  }
```
with:
```cpp
  if (resolved_.ur) {
    const double rho_ur = resolved_.ur->Rho(ppw->pvecback.data());
    ppw->delta_rho += rho_ur * ppw->rsa_delta_ur;
    ppw->delta_p += 1. / 3. * rho_ur * ppw->rsa_delta_ur;
    ppw->rho_plus_p_theta += 4. / 3. * rho_ur * ppw->rsa_theta_ur;
  }
```
Note: hoisting `Rho()` into `rho_ur` is the same value used three times; multiplication order is unchanged, so the result is byte-identical.

- [ ] **Step 4: Replace the tensor-branch UR/NCDM lookups (Site 3)**

In `perturb_total_stress_energy` (tensor branch), replace:
```cpp
      if (ppt->tensor_method == tm_exact) {
        if (auto* ur = all_species_.find("UR"))
          rho_relativistic += (*ur)->Rho(ppw->pvecback.data());
      }

      if (ppt->tensor_method == tm_massless_approximation) {
        if (auto* ur = all_species_.find("UR"))
          rho_relativistic += (*ur)->Rho(ppw->pvecback.data());

        if (HasNcdm(all_species_)) {
          for (auto& sp : all_species_)
            rho_relativistic += sp->TensorMasslessRelativisticRho(ppw->pvecback.data());
        }
      }
```
with:
```cpp
      if (ppt->tensor_method == tm_exact) {
        if (resolved_.ur)
          rho_relativistic += resolved_.ur->Rho(ppw->pvecback.data());
      }

      if (ppt->tensor_method == tm_massless_approximation) {
        if (resolved_.ur)
          rho_relativistic += resolved_.ur->Rho(ppw->pvecback.data());

        if (resolved_.has_ncdm) {
          for (auto& sp : all_species_)
            rho_relativistic += sp->TensorMasslessRelativisticRho(ppw->pvecback.data());
        }
      }
```
(The ncdm `TensorMasslessRelativisticRho` delegation loop is intentionally left unchanged.)

- [ ] **Step 5: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: clean build.

- [ ] **Step 6: Verify byte-identical on base_2018 (exercises Site 2 — base has N_ur)**

Run:
```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  && cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat /tmp/perturb_handles_golden/pk.dat \
  && echo "BYTE-IDENTICAL OK"
```
Expected: `BYTE-IDENTICAL OK`. (Site 3 is tensor-only; it is gated by the tensor scenarios in TEST_LEVEL=2, Task 6.)

- [ ] **Step 7: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturb: resolve UR handle in rsa + tensor stress-energy (Site 2/3)"
```

---

### Task 4: Site 4a + 4b — CDM h-source layout cache and Lambda rho_tot exclusion

**Files:**
- Modify: `source/perturbations.h` (add `cdm_layout` to `perturb_vector`)
- Modify: `source/perturbations_module.cpp` (`perturb_vector_init` wiring; `perturb_sources_member` CDM + Lambda sites)

**Interfaces:**
- Consumes: `resolved_.cdm_index` (`std::size_t`), `resolved_.lambda` (`const BaseSpecies*`) from Task 1.
- Produces: `perturb_vector::cdm_layout` (`BaseSpecies::PerturbLayout*`, nullptr if CDM absent).

- [ ] **Step 1: Add the `cdm_layout` field to `perturb_vector`**

In `source/perturbations.h`, replace:
```cpp
  BaseSpecies::PerturbLayout* photon_layout = nullptr;
  BaseSpecies::PerturbLayout* baryon_layout = nullptr;
```
with:
```cpp
  BaseSpecies::PerturbLayout* photon_layout = nullptr;
  BaseSpecies::PerturbLayout* baryon_layout = nullptr;
  BaseSpecies::PerturbLayout* cdm_layout    = nullptr;
```

- [ ] **Step 2: Wire `cdm_layout` in `perturb_vector_init`**

In `source/perturbations_module.cpp`, after:
```cpp
  ppv->photon_layout = ppv->species_layouts[all_species_.photons_index()].get();
  ppv->baryon_layout = ppv->species_layouts[all_species_.baryons_index()].get();
```
add:
```cpp
  if (resolved_.cdm_index < all_species_.size())
    ppv->cdm_layout = ppv->species_layouts[resolved_.cdm_index].get();
```

- [ ] **Step 3: Replace the CDM h-source lookup (Site 4a)**

In `perturb_sources_member`, replace:
```cpp
      if (has_source_h_) {
        const size_t cdm_i  = all_species_.index_of("CDM");
        const auto& cdm_lay = static_cast<const CDMSpecies::PerturbLayout&>(
            *ppw->pv->species_layouts[cdm_i]);
        _set_source_(index_tp_h_) = -2 * y[cdm_lay.idx_delta];
      }
```
with:
```cpp
      if (has_source_h_) {
        const auto& cdm_lay = static_cast<const CDMSpecies::PerturbLayout&>(*ppw->pv->cdm_layout);
        _set_source_(index_tp_h_) = -2 * y[cdm_lay.idx_delta];
      }
```

- [ ] **Step 4: Replace the Lambda rho_tot exclusion (Site 4b)**

In `perturb_sources_member`, replace:
```cpp
      double rho_tot;
      if (all_species_.count("Lambda")) {
        rho_tot = pvecback[background_module_->index_bg_rho_tot_] -
                  all_species_.at("Lambda")->Rho(pvecback);
      }
      else {
        rho_tot = pvecback[background_module_->index_bg_rho_tot_];
      }
```
with:
```cpp
      double rho_tot;
      if (resolved_.lambda) {
        rho_tot = pvecback[background_module_->index_bg_rho_tot_] -
                  resolved_.lambda->Rho(pvecback);
      }
      else {
        rho_tot = pvecback[background_module_->index_bg_rho_tot_];
      }
```

- [ ] **Step 5: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: clean build.

- [ ] **Step 6: Verify byte-identical on base_2018**

Run:
```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  && cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat /tmp/perturb_handles_golden/pk.dat \
  && echo "BYTE-IDENTICAL OK"
```
Expected: `BYTE-IDENTICAL OK`. (The `has_source_h_` / `delta_tot` source paths are exercised by the `nCl sCl` / density scenarios in TEST_LEVEL=2, Task 6; base_2018 confirms the common path is untouched.)

- [ ] **Step 7: Commit**

```bash
git add source/perturbations.h source/perturbations_module.cpp
git commit -m "perturb: cache cdm_layout + resolve Lambda handle in sources (Site 4a/4b)"
```

---

### Task 5: Site 5 — `perturb_approximations` triggers

**Files:**
- Modify: `source/perturbations_module.cpp` (`perturb_approximations`, the five lookup sites)

**Interfaces:**
- Consumes: `resolved_.idm_dr_idr`, `resolved_.idm_drmd`, `resolved_.ur`, `resolved_.has_ncdm` from Task 1.

- [ ] **Step 1: Replace the TCA-idm_dr trigger block**

In `perturb_approximations`, replace:
```cpp
    if (all_species_.count("IDM_DR_IDR")) {
      const auto& idm_dr_idr = static_cast<const IDM_DR_IDR_Species&>(
          *all_species_.at("IDM_DR_IDR"));
```
with:
```cpp
    if (resolved_.idm_dr_idr) {
      const auto& idm_dr_idr = *resolved_.idm_dr_idr;
```
(The block body that follows — using `idm_dr_idr.idm_dr().nindex_idm_dr()` and `idm_dr_idr.idr().idr_nature()` — is unchanged; both accessors are const.)

- [ ] **Step 2: Replace the TCA-idm_drmd trigger block**

In `perturb_approximations`, replace:
```cpp
    if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
      double Rint, csp2, Gint;
      double conformalH = ppw->pvecback[background_module_->index_bg_H_] *
                          ppw->pvecback[background_module_->index_bg_a_];

      background_module_->background_idm_drmd(ppw->pvecback[background_module_->index_bg_a_],
                                              static_cast<IDM_DRMD_IDR_DRMD_Species&>(
                                                  *all_species_.at("IDM_DRMD_IDR_DRMD"))
                                                      .idm_drmd()
                                                      .Rho(ppw->pvecback.data()) /
                                                  static_cast<IDM_DRMD_IDR_DRMD_Species&>(
                                                      *all_species_.at("IDM_DRMD_IDR_DRMD"))
                                                      .idr_drmd()
                                                      .Rho(ppw->pvecback.data()),
                                              &Rint,
                                              &csp2,
                                              &Gint);
```
with:
```cpp
    if (resolved_.idm_drmd) {
      double Rint, csp2, Gint;
      double conformalH = ppw->pvecback[background_module_->index_bg_H_] *
                          ppw->pvecback[background_module_->index_bg_a_];

      background_module_->background_idm_drmd(
          ppw->pvecback[background_module_->index_bg_a_],
          resolved_.idm_drmd->idm_drmd().Rho(ppw->pvecback.data()) /
              resolved_.idm_drmd->idr_drmd().Rho(ppw->pvecback.data()),
          &Rint,
          &csp2,
          &Gint);
```

- [ ] **Step 3: Replace the RSA-IDR streaming trigger guard**

In `perturb_approximations`, replace:
```cpp
    if (all_species_.count("IDM_DR_IDR")) {
      auto& idm_idr = static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
      if (idm_idr.idm_dr().IsPresent()) {
```
with:
```cpp
    if (resolved_.idm_dr_idr) {
      const auto& idm_idr = *resolved_.idm_dr_idr;
      if (idm_idr.idm_dr().IsPresent()) {
```
(The remainder of the block — reading `idm_idr.idm_dr().nindex_idm_dr()` — is unchanged; `IsPresent()` and `nindex_idm_dr()` are const.)

- [ ] **Step 4: Replace the UFA trigger guard**

In `perturb_approximations`, replace:
```cpp
    if (all_species_.count("UR")) {
      if ((tau / tau_k > ppr->ur_fluid_trigger_tau_over_tau_k) &&
```
with:
```cpp
    if (resolved_.ur) {
      if ((tau / tau_k > ppr->ur_fluid_trigger_tau_over_tau_k) &&
```

- [ ] **Step 5: Replace the ncdmfa trigger guard**

In `perturb_approximations`, replace:
```cpp
    if (HasNcdm(all_species_)) {
      if ((tau / tau_k > ppr->ncdm_fluid_trigger_tau_over_tau_k) &&
```
with:
```cpp
    if (resolved_.has_ncdm) {
      if ((tau / tau_k > ppr->ncdm_fluid_trigger_tau_over_tau_k) &&
```

- [ ] **Step 6: Build**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
```
Expected: clean build.

- [ ] **Step 7: Verify byte-identical on base_2018**

Run:
```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  && cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat /tmp/perturb_handles_golden/pk.dat \
  && echo "BYTE-IDENTICAL OK"
```
Expected: `BYTE-IDENTICAL OK`. (base_2018 exercises the UFA/ncdmfa presence guards via N_ur + ncdm; idm_dr/drmd paths are gated by TEST_LEVEL=2, Task 6.)

- [ ] **Step 8: Confirm no per-step/per-sample/per-k lookups remain in the hot functions**

Run:
```bash
grep -nE 'all_species_\.(count|at|find|index_of)|HasNcdm\(all_species_\)' source/perturbations_module.cpp
```
Then confirm by eye that **none** of the listed matches fall inside these five hot functions: `perturb_approximations`, `perturb_einstein`, `perturb_total_stress_energy`, `perturb_sources_member`, `perturb_rsa_delta_and_theta`. (Find each function's current line range with `grep -n 'PerturbationsModule::perturb_' source/perturbations_module.cpp`.) Matches in setup/registration functions — `perturb_init`, `perturb_indices_of_perturbs`, `perturb_get_k_list`, `perturb_workspace_init`, `perturb_solve`, `perturb_find_approximation_switches`, `perturb_vector_init`, `perturb_initial_conditions`, `perturb_print_variables_member` — are expected and out of scope (spec §6).

- [ ] **Step 9: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturb: resolve handles in perturb_approximations (Site 5)"
```

---

### Task 6: Format, full regression gate, and rigorous benchmark

**Files:**
- Modify: any touched file (formatting only)

- [ ] **Step 1: clang-format every touched file**

Run:
```bash
clang-format --version   # expect 22.1.3
clang-format -i source/perturbations_module.h source/perturbations_module.cpp source/perturbations.h
git diff --stat
```
Expected: only whitespace/formatting changes, if any.

- [ ] **Step 2: Rebuild and re-confirm byte-identical**

Run:
```bash
cmake --build build_prof -j 2>&1 | tail -3
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  && cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat /tmp/perturb_handles_golden/pk.dat \
  && echo "BYTE-IDENTICAL OK"
```
Expected: clean build, `BYTE-IDENTICAL OK`.

- [ ] **Step 3: Physics-tolerance gate (sanity; should pass trivially given byte-identity)**

Run:
```bash
python3 .superpowers/sdd/tol_check.py \
  output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/perturb_handles_golden/cl_lensed.dat \
  output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat       /tmp/perturb_handles_golden/pk.dat
echo "exit=$?"
```
Expected: `exit=0`.

- [ ] **Step 4: Full regression suite vs master classyref (the comprehensive gate)**

This covers tensor modes (Site 3), `nCl sCl` density sources (Site 4a/b), both gauges, and ncdm — the paths base_2018 does not exercise. Run the project's standard:
```bash
TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python3 -m pytest python/test_class.py -q 2>&1 | tail -20
```
Expected: all pass, 0 failures (the established `2538 passed / 0 failed` shape against the current master classyref). If classyref is stale relative to the current master, regenerate it first per the project's classyref procedure, then re-run.

- [ ] **Step 5: Rigorous paired A/B benchmark (the real perf number)**

Per spec §8 and `feedback_no_bit_identical_requirement` / the memory's benchmark protocol: build BOTH master and this branch fresh from **isolated git worktrees** with **identical flags** (`-O3 -DNDEBUG`, same compiler) into separate dirs; run order-alternating paired A/B on `base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini`, `OMP_NUM_THREADS=1`, ≥40+40 rounds; report min/median + bootstrap CI and paired mean/se on the Perturbations stage. Confirm algorithmic-vs-layout with `-falign-functions=64` variants of both binaries.
```bash
# Reuse the /tmp/bench_*.sh A/B harnesses + analyzers from the memory if present.
git worktree add /tmp/rh_master "$(git merge-base HEAD master)"
git worktree add /tmp/rh_branch HEAD
# In each worktree: configure+build class_profiled with identical flags
#   (cmake -DCMAKE_BUILD_TYPE=Release; -O3 -DNDEBUG; same compiler) into its own
#   build dir, then run the order-alternating paired A/B harness on base_2018.
```
Expected: branch Perturbations stage reproducibly faster (target ~1–3%), surviving the layout probe. Record the number in the PR description. Clean up worktrees: `git worktree remove /tmp/rh_master /tmp/rh_branch`.

- [ ] **Step 6: Commit formatting (if any) and finalize**

```bash
git add source/perturbations_module.h source/perturbations_module.cpp source/perturbations.h
git commit -m "perturb: clang-format resolved-species-handles changes" || echo "nothing to format-commit"
```

- [ ] **Step 7: Open the PR**

Summarize: the resolved-handle design, byte-identical across the suite, the measured Perturbations-stage speedup + layout-probe confirmation, and a pointer to the spec `docs/superpowers/specs/2026-06-22-perturbations-resolved-species-handles-design.md`.

---

## Notes for the implementer

- If a const accessor unexpectedly fails to compile off a `const` handle, do **not** drop to a non-const handle — the spec relies on const correctness, and all required accessors (`idr()`, `idm_dr()`, `idr_drmd()`, `idm_drmd()`, `alpha_idm_dr()`, `nindex_idm_dr()`, `idr_nature()`, `IsPresent()`, `Rho()`) were verified const-callable. A compile failure means a genuine non-const accessor; add a trivial `const` overload to that getter rather than weakening the handle.
- Each task is independently byte-identical against the same Task 0 golden; a task that fails `cmp` has a logic error — debug it before moving on (see `superpowers:systematic-debugging`).
- The hot functions are the only legitimate edit sites; lookups in setup/registration functions are explicitly out of scope (spec §6).
