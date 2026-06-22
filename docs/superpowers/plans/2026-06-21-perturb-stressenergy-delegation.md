# StressEnergy struct-delegation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the four per-species stress-energy virtuals (`DeltaRho`, `RhoPlusPTheta`, `DeltaP`, `RhoPlusPShear`) on `BaseSpecies` with a single struct-returning `StressEnergy` virtual, cutting the perturbations RHS from ~7 virtual dispatches per species per ODE step to one.

**Architecture:** Introduce a value struct `StressEnergyContribution` and a `StressEnergy(...)` virtual. Land it behind a temporary non-pure base default that delegates to the four existing virtuals (so every intermediate commit compiles and is byte-identical), route all call sites through it, fold the anisotropic-stress pre-pass into the single loop, give every concrete species a direct override (NCDM/DNCDM a single fused q-pass), then delete the base default and the four old virtuals to make `StressEnergy` pure.

**Tech Stack:** C++17, CMake (Unix Makefiles, `-O3 -DNDEBUG`), pytest test harness (`python/test_class.py`), the `class_profiled` per-stage timer.

## Global Constraints

- **Physics agreement to ~0.1% is the bar — NOT bit-identity.** This is a high-performance research code; write clean, idiomatic, fast code and do NOT contort it (e.g. splitting a fused loop into four, or forcing a particular FP operand order) to preserve byte-for-byte output. Benign vectorization/reduction ULP drift is expected and fine. Each task must keep `Cl^TT` and `P(k)` within `1e-3` relative of the clean-master reference (`/tmp/ref_master`, captured in Task 1) via the tolerance check below. (`Cl^TT`/`P(k)` are checked because they have no zero-crossings; never apply a blind max-rel metric to `Cl^TE`/`EE`/`BB`.) Baselines/`classyref` are regenerated at the end (Task 9), not defended per-step.
- **Scope is the StressEnergy delegation only.** Do NOT touch the `all_species_.count("IDM_DR_IDR")` string scans in `perturb_einstein` — that count-cache harvest is a separate follow-up PR.
- **No by-reference accumulator.** The struct is returned by value (measured: RVO'd into registers at `-O3`, the by-ref variant gave no gain and broke bit-identity).
- **`Rho`/`P` stay as their own virtuals.** They have many background callers outside stress-energy; the struct only carries copies.
- Field order is fixed everywhere: `rho, p, delta_rho, rho_plus_p_theta, delta_p, rho_plus_p_shear`.
- Build: `cmake --build build -j` → `./class`. Benchmark ini: `base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini` (LCDM + 1 massive ncdm + lensing + mPk).
- Reference prototype: a verified byte-identical prototype of the hot-species overrides and the NCDM fused pass lives in `git stash@{0}` ("perturb perf experiment…"). Use `git stash show -p stash@{0}` to read it. **Ignore its `AddStressEnergy(...)` / `StressEnergyAccumulator` methods — that is the rejected by-ref design; port only the `StressEnergy(...)` methods and the `StressEnergyContribution` struct.**

### Standard verification block (referenced as "run VERIFY")

```bash
cmake --build build -j 2>&1 | tail -5            # must succeed
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/tmp/class_run.log 2>&1
python3 .superpowers/sdd/tol_check.py \
  output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/ref_master/cl_lensed.dat \
  output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat        /tmp/ref_master/pk.dat 1e-3
```
Expected: `-> PASS` (worst `Cl^TT`/`P(k)` max_rel <= 1e-3 vs master). A reported `max_rel` of a few `e-5` is normal ULP drift and PASSES; only a result above `1e-3` (or a build failure / shape mismatch) is a real failure. Wherever a task below says "cmp byte-identical" or "cmp exit 0", use THIS tolerance check instead — bit-identity is not required.

---

### Task 1: Introduce `StressEnergyContribution` + scaffold base-default `StressEnergy`; capture reference

**Files:**
- Modify: `species/base_species.h` (add struct + non-pure virtual, just before the existing `DeltaRho` declaration ~line 330)

**Interfaces:**
- Produces: `struct BaseSpecies::StressEnergyContribution { double rho, p, delta_rho, rho_plus_p_theta, delta_p, rho_plus_p_shear; StressEnergyContribution& operator+=(const StressEnergyContribution&); }` and `virtual StressEnergyContribution BaseSpecies::StressEnergy(const PerturbLayout&, const perturb_vector*, const double* y, const double* pvecback, const perturb_workspace*) const`.

- [ ] **Step 1: Capture the clean-master reference** (the branch HEAD is master + the spec doc, so its output == master)

```bash
cmake --build build -j 2>&1 | tail -3
mkdir -p /tmp/ref_master
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null 2>&1
cp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/ref_master/cl_lensed.dat
cp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat        /tmp/ref_master/pk.dat
```
Expected: build succeeds; two files copied.

- [ ] **Step 2: Add the struct + scaffold default** in `species/base_species.h`, immediately above the `/** Density contribution ... DeltaRho` declaration (~line 330):

```cpp
  /** The scalar stress-energy contribution of one species at one ODE step
   *  (Ma & Bertschinger). Returned by value: six doubles, RVO'd into registers
   *  at -O3 (the by-reference accumulator variant was measured to give nothing). */
  struct StressEnergyContribution {
    double rho              = 0.;  // background ρ
    double p                = 0.;  // background P
    double delta_rho        = 0.;  // δρ
    double rho_plus_p_theta = 0.;  // (ρ+P)θ
    double delta_p          = 0.;  // δp
    double rho_plus_p_shear = 0.;  // (ρ+P)σ

    StressEnergyContribution& operator+=(const StressEnergyContribution& o) {
      rho += o.rho; p += o.p; delta_rho += o.delta_rho;
      rho_plus_p_theta += o.rho_plus_p_theta; delta_p += o.delta_p;
      rho_plus_p_shear += o.rho_plus_p_shear;
      return *this;
    }
  };

  /** All scalar stress-energy perturbations in one call.
   *  TEMPORARY base default delegates to the four individual virtuals; this
   *  default is removed (and the virtual made pure) once every species
   *  overrides it (final task). Species override this with a direct, single
   *  computation. */
  virtual StressEnergyContribution StressEnergy(const PerturbLayout& layout,
                                                const perturb_vector* pv,
                                                const double* y,
                                                const double* pvecback,
                                                const perturb_workspace* ppw) const {
    return {Rho(pvecback), P(pvecback),
            DeltaRho(layout, pv, y, pvecback, ppw),
            RhoPlusPTheta(layout, pv, y, pvecback, ppw),
            DeltaP(layout, pv, y, pvecback, ppw),
            RhoPlusPShear(layout, pv, y, pvecback, ppw)};
  }
```

- [ ] **Step 3: Build** — `cmake --build build -j 2>&1 | tail -5`. Expected: succeeds (header-only change, nothing calls `StressEnergy` yet).

- [ ] **Step 4: Commit**

```bash
git add species/base_species.h
git commit -m "perturb: add StressEnergyContribution struct + scaffold base StressEnergy"
```

---

### Task 2: Route `perturb_total_stress_energy` through `StressEnergy`; fold shear; drop pre-pass

**Files:**
- Modify: `source/perturbations_module.cpp` (scalar block of `perturb_total_stress_energy`, ~lines 4822-4865)

**Interfaces:**
- Consumes: `BaseSpecies::StressEnergy(...)` (Task 1).

- [ ] **Step 1: Replace the pre-pass + per-species loop.** In `perturb_total_stress_energy`, delete the shear pre-pass (the comment block + `for (const auto& e : pv->active_species) ppw->rho_plus_p_shear += e.species->RhoPlusPShear(...)`, ~lines 4832-4839) and replace the main loop body so it makes one `StressEnergy` call and folds shear in. The block from `ppw->delta_rho = 0.;` through the end of the species `for` loop becomes:

```cpp
    ppw->delta_rho        = 0.;
    ppw->rho_plus_p_theta = 0.;
    ppw->rho_plus_p_shear = 0.;
    ppw->delta_p          = 0.;
    ppw->rho_plus_p_tot   = 0.;
    struct Tally {
      double drho = 0, rho = 0, rppt = 0, rpp = 0;
    } cold, warm;
    const bool tally = has_source_delta_m_ || has_source_theta_m_;

    /* Single pass: each species returns all six stress-energy quantities in one
       virtual call. rho_plus_p_shear is accumulated here too (the former
       separate pre-pass was needed only by ScalarFieldSpecies' Newtonian-gauge
       DeltaRho/DeltaP, which is a hard class_test error at init — see issue
       #285; any future Newtonian-scf fix must restore the complete shear before
       the scalar field reads it, e.g. by ordering it last or reinstating a
       pre-pass). Accumulation order over active_species is unchanged, so totals
       are bit-identical. */
    for (const auto& e : pv->active_species) {
      const BaseSpecies::StressEnergyContribution se =
          e.species->StressEnergy(*e.layout, pv, y, pb, ppw);
      const double rpp = se.rho + se.p;

      ppw->delta_rho        += se.delta_rho;
      ppw->rho_plus_p_theta += se.rho_plus_p_theta;
      ppw->rho_plus_p_shear += se.rho_plus_p_shear;
      ppw->delta_p          += se.delta_p;
      ppw->rho_plus_p_tot   += rpp;

      if (tally && e.clusters_as_matter) {
        const double m_rho   = se.rho - 3. * se.p;
        const double m_drho  = se.delta_rho - 3. * se.delta_p;
        const double m_rppt  = (rpp > 0.) ? m_rho * se.rho_plus_p_theta / rpp : 0.;
        Tally& t             = e.is_cold ? cold : warm;
        t.drho              += m_drho;
        t.rho               += m_rho;
        t.rppt              += m_rppt;
        t.rpp               += m_rho;
      }
    }
```
Leave everything after the loop (`delta_cb`/`theta_cb`, the `ppf_fluid()` block, `delta_m`/`theta_m`) exactly as-is.

- [ ] **Step 2: Run VERIFY** (see Global Constraints). Expected: build succeeds; both `cmp` exit 0.

If a `cmp` differs, the most likely cause is an accidental change to accumulation order or to the tally arithmetic — diff against the original loop and restore exact operation order.

- [ ] **Step 3: Commit**

```bash
git add source/perturbations_module.cpp
git commit -m "perturb: single StressEnergy call per species; fold shear into the main loop"
```

---

### Task 3: Switch the remaining (cold-path) callers to `StressEnergy`

**Files:**
- Modify: `source/perturbations_module.cpp:3984-3985` (IC analogue)
- Modify: `species/ncdm_species.cpp` (`FillSources` ~285/297, `PrintVariables` ~479-482)
- Modify: `species/fluid.cpp:256-258` (`PrintVariables`, non-PPF branch)

**Interfaces:**
- Consumes: `BaseSpecies::StressEnergy(...)` (still the base default here).

- [ ] **Step 1: IC analogue** (`perturbations_module.cpp` ~3982-3986). Replace the separate `DeltaRho` + `RhoPlusPTheta` calls with one `StressEnergy` call. Read the surrounding lines first; the edit is: compute `const auto se = sp->StressEnergy(layout, ppw->pv.get(), ppw->pv->y.data(), ppw->pvecback.data(), ppw);` once, then use `se.delta_rho` where `DeltaRho(...)` was and `se.rho_plus_p_theta` where `RhoPlusPTheta(...)` was.

- [ ] **Step 2: NCDM `PrintVariables`** (`ncdm_species.cpp` ~479-482). Replace the four separate calls:

```cpp
    const auto se = StressEnergy(layout, pv, y, pvecback, ppw);
    delta_ncdm = (rho_ncdm_bg > 0.) ? se.delta_rho / rho_ncdm_bg : 0.;
    theta_ncdm = (rho_plus_p > 0.) ? se.rho_plus_p_theta / rho_plus_p : 0.;
    shear_ncdm = (rho_plus_p > 0.) ? se.rho_plus_p_shear / rho_plus_p : 0.;
    const double delta_p_ncdm   = se.delta_p;
```
(Keep the exact divisor variables `rho_ncdm_bg`/`rho_plus_p` already computed above those lines.)

- [ ] **Step 3: NCDM `FillSources`** (`ncdm_species.cpp` ~285/297). Read the two lines; replace the `DeltaRho(...)` and `RhoPlusPTheta(...)` calls with one `const auto se = StressEnergy(...);` and use `se.delta_rho` / `se.rho_plus_p_theta` in the same expressions (preserve the `/ pvecback[index_bg_rho_]` and `/ (rho+p)` divisions exactly).

- [ ] **Step 4: Fluid `PrintVariables`** (`fluid.cpp:256-258`, inside the `else` non-PPF branch):

```cpp
      const auto se        = StressEnergy(layout, ppw->pv.get(), y, ppw->pvecback.data(), ppw);
      delta_rho_fld        = se.delta_rho;
      rho_plus_p_theta_fld = se.rho_plus_p_theta;
      delta_p_fld          = se.delta_p;
```

- [ ] **Step 5: Run VERIFY.** Expected: build succeeds; both `cmp` exit 0. (These are output/IC paths; the benchmark exercises NCDM `FillSources`. The output-column paths are covered fully by `TEST_LEVEL=2` in Task 9.)

- [ ] **Step 6: Commit**

```bash
git add source/perturbations_module.cpp species/ncdm_species.cpp species/fluid.cpp
git commit -m "perturb: route IC and output callers through StressEnergy"
```

---

### Task 4: Direct overrides for the four hot leaf species (Photons, Baryons, CDM, UR)

**Files:**
- Modify: `species/cdm.h` + `species/cdm.cpp`
- Modify: `species/photons.h` + `species/photons.cpp`
- Modify: `species/baryons.h` + `species/baryons.cpp`
- Modify: `species/ultra_relativistic.h` + `species/ultra_relativistic.cpp`

**Interfaces:**
- Produces: `StressEnergyContribution <Species>::StressEnergy(...) const override` for each of the four.

**Recipe (applies to each species):** add a `StressEnergy(...) const override` declaration to the header (next to the existing `RhoPlusPShear` declaration), and a definition that returns a `StressEnergyContribution` whose six fields reproduce — bit-for-bit, same expressions and operand order — `Rho`, `P`, and the four existing methods, computing shared quantities once. The verified bodies are in `git stash@{0}` (read with `git stash show -p stash@{0}`); port the `StressEnergy` methods only (skip the `AddStressEnergy` ones).

- [ ] **Step 1: CDM** — header declaration in `species/cdm.h` (after `RhoPlusPShear`):

```cpp
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;
```
Definition in `species/cdm.cpp` (verbatim from the verified prototype):

```cpp
// Direct fused override: one dispatch, no inner virtual calls. Each field
// reproduces its individual function bit-for-bit (P, delta_p, shear ≡ 0).
BaseSpecies::StressEnergyContribution CDMSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* /*pv*/,
    const double* y,
    const double* pvecback,
    const perturb_workspace* /*ppw*/) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho              = pvecback[index_bg_rho_cdm_];
  se.delta_rho        = pvecback[index_bg_rho_] * y[layout.idx_delta];
  se.rho_plus_p_theta = (layout.idx_theta >= 0) ? pvecback[index_bg_rho_] * y[layout.idx_theta] : 0.;
  return se;
}
```

- [ ] **Step 2: Photons, Baryons, UltraRelativistic** — same recipe. Add the override declaration to each header, and port each `StressEnergy` body from `stash@{0}`. Each body must reproduce that species' existing `DeltaRho`/`RhoPlusPTheta`/`DeltaP`/`RhoPlusPShear` and `Rho`/`P` with identical expressions (photons/baryons carry `delta_p`; photons carries shear; UR carries shear; check each against its current four methods).

- [ ] **Step 3: Run VERIFY.** Expected: build succeeds; both `cmp` exit 0. (These four species are all active in the benchmark, so any arithmetic drift shows immediately.)

- [ ] **Step 4: Commit**

```bash
git add species/cdm.* species/photons.* species/baryons.* species/ultra_relativistic.*
git commit -m "perturb: direct StressEnergy overrides for photons/baryons/cdm/ur"
```

---

### Task 5: Fused `StressEnergy` for the NCDM family (`NCDMSpecies` + `DNCDMSpecies`)

**Files:**
- Modify: `species/ncdm_species.h` + `species/ncdm_species.cpp`
- Modify: `species/dncdm_species.h` + `species/dncdm_species.cpp`

**Interfaces:**
- Produces: `StressEnergyContribution NCDMSpecies::StressEnergy(...) const override` (inherited by `GreyBodyNCDMSpecies`, `NCDMInteractingSpecies`) and `StressEnergyContribution DNCDMSpecies::StressEnergy(...) const override`.

- [ ] **Step 1: `NCDMSpecies::StressEnergy`** — add the override declaration to `species/ncdm_species.h` (after `RhoPlusPShear`, ~line 91) and the definition to `species/ncdm_species.cpp` (after `RhoPlusPShear`, ~line 628), as ONE clean fused pass over the q grid (from the `stash@{0}` prototype): compute `epsilon` and `factor = factor_ * pow(a0/a,4)` ONCE per bin, accumulate all four quantities in a single loop, handle the `ncdmfa_on` branch, set `se.rho`/`se.p` from `pvecback[index_bg_rho_]`/`[index_bg_p_]`. Do NOT split into four loops to preserve bit-identity — a single fused pass that drifts at the ULP level (~1e-5 on `Cl^TT`, validated by the tolerance check) is exactly what's wanted. (Already implemented in commit `fadd1f36`; this step documents intent.)

- [ ] **Step 2: `DNCDMSpecies::StressEnergy`** — `DNCDMSpecies` (decaying ncdm) has its **own** four methods (`dncdm_species.cpp:538+`) that differ from `NCDMSpecies` (decay factor), so it needs its own override. Read its current `DeltaRho`/`RhoPlusPTheta`/`DeltaP`/`RhoPlusPShear`, then write a `StressEnergy` that fuses them into one q-pass using the same per-term expressions and operand order, setting `se.rho`/`se.p` as its background `Rho`/`P` do. Add the override declaration to `species/dncdm_species.h`.

- [ ] **Step 3: Run VERIFY** plus a DNCDM check:

```bash
# 1-massive-nu (exercises NCDMSpecies::StressEnergy):
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null 2>&1
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/ref_master/cl_lensed.dat
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat        /tmp/ref_master/pk.dat
```
For DNCDM byte-identity, rely on `TEST_LEVEL=2`'s `test_dncdm_dr_computes` reference check in Task 9 (the standalone dncdm inis under `test/` write a different root; a dedicated cmp is optional — if desired, capture a master reference for `test/dotsyntax_dncdm.ini` the same way as Task 1 before editing).

Expected: both 1-massive-ν `cmp` exit 0.

- [ ] **Step 4: Commit**

```bash
git add species/ncdm_species.* species/dncdm_species.*
git commit -m "perturb: fused single q-pass StressEnergy for NCDM and DNCDM"
```

---

### Task 6: Overrides for the remaining simple leaves (Fluid, DCDM, DarkRadiation, ScalarField, Lambda)

**Files:**
- Modify: `species/fluid.{h,cpp}`, `species/dcdm.{h,cpp}`, `species/dark_radiation_species.{h,cpp}`, `species/scalar_field.{h,cpp}`, `species/lambda.h`

**Interfaces:**
- Produces: `StressEnergyContribution <Species>::StressEnergy(...) const override` for each.

**Recipe (self-contained, mandatory):** for each species, add the override declaration to its header and a definition that **inlines** the exact bodies of its current `Rho`, `P`, `DeltaRho`, `RhoPlusPTheta`, `DeltaP`, `RhoPlusPShear` into the six struct fields. The override **must not call the species' own four methods** — those four methods are deleted in Task 9, so a body that calls them would stop compiling. Read each species' four current method bodies first, then transcribe each expression (same operands, same order — byte-identity is the bar) into the matching field, computing any shared quantity once. `LambdaSpecies` declares its four methods inline in `lambda.h` (no `.cpp`); add `StressEnergy` inline there too — its perturbation fields are all `0.`, with `se.rho = Rho(pvecback)`, `se.p = P(pvecback)`.

Worked example — **DCDM** (read `species/dcdm.cpp` first; the shape mirrors CDM in Task 4, but transcribe DCDM's own `Rho`/`P`/`DeltaRho`/... bodies — fill each field with that method's actual returned expression, e.g. `se.delta_rho = <body of DCDMSpecies::DeltaRho>;`):
```cpp
BaseSpecies::StressEnergyContribution DCDMSpecies::StressEnergy(
    const BaseSpecies::PerturbLayout& base, const perturb_vector* pv,
    const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& layout = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se;
  se.rho              = /* body of DCDMSpecies::Rho(pvecback)            */;
  se.p                = /* body of DCDMSpecies::P(pvecback)              */;
  se.delta_rho        = /* body of DCDMSpecies::DeltaRho(...)            */;
  se.rho_plus_p_theta = /* body of DCDMSpecies::RhoPlusPTheta(...)       */;
  se.delta_p          = /* body of DCDMSpecies::DeltaP(...)              */;
  se.rho_plus_p_shear = /* body of DCDMSpecies::RhoPlusPShear(...)       */;
  return se;
}
```
This self-contained form survives the Task 9 deletion of the four methods. **Do not** route through the *base default* and **do not** call the four methods — inline their bodies so the override is concrete and standalone.

- [ ] **Step 1:** Add the override to Fluid, DCDM, DarkRadiation, ScalarField (each `cpp`+`h`) and Lambda (`h` inline), per the recipe.
- [ ] **Step 2: Run VERIFY.** Expected: both `cmp` exit 0. (DCDM/DR/Lambda inactive in the benchmark but must still compile and not perturb the common path; full coverage lands in Task 9's `TEST_LEVEL=2`.)
- [ ] **Step 3: Commit**

```bash
git add species/fluid.* species/dcdm.* species/dark_radiation_species.* species/scalar_field.* species/lambda.h
git commit -m "perturb: StressEnergy overrides for fluid/dcdm/dr/scalar-field/lambda"
```

---

### Task 7: Overrides for the interacting leaves (IDM_DR, IDR, IDM_DRMD, IDR_DRMD)

**Files:**
- Modify: `species/interacting_species.cpp` (defs at ~52, ~206, ~358, ~437) and the headers `species/idm_dr.h`, `species/idr.h`, `species/idm_drmd.h`, `species/idr_drmd.h`

**Interfaces:**
- Produces: `StressEnergyContribution <Species>::StressEnergy(...) const override` for `IDM_DRSpecies`, `IDRSpecies`, `IDM_DRMDSpecies`, `IDR_DRMDSpecies`.

- [ ] **Step 1:** For each of the four classes, add the override declaration to its header and a definition in `interacting_species.cpp` that **inlines** the bodies of its existing `Rho`/`P`/`DeltaRho`/`RhoPlusPTheta`/`DeltaP`/`RhoPlusPShear` into one struct (the self-contained form from the Task 6 recipe — **do not call** the four methods, since Task 9 deletes them; transcribe each expression byte-identically).
- [ ] **Step 2: Build only** — `cmake --build build -j 2>&1 | tail -5` (these species are not in the benchmark; correctness is gated by `TEST_LEVEL=2` in Task 9). Then run VERIFY to confirm the common path is still byte-identical. Expected: build succeeds; both `cmp` exit 0.
- [ ] **Step 3: Commit**

```bash
git add species/interacting_species.cpp species/idm_dr.h species/idr.h species/idm_drmd.h species/idr_drmd.h
git commit -m "perturb: StressEnergy overrides for the interacting (idm/idr) leaves"
```

---

### Task 8: Composite overrides via struct sum (DCDM_DR, DNCDM_DR, IDM_DR_IDR, IDM_DRMD_IDR_DRMD)

**Files:**
- Modify: `species/dcdm_dr_species.{h,cpp}`, `species/dncdm_dr_species.{h,cpp}`, `species/idm_dr_idr_species.{h,cpp}`, `species/idm_drmd_idr_drmd_species.{h,cpp}`

**Interfaces:**
- Consumes: child species' `StressEnergy` (Tasks 5-7) and `StressEnergyContribution::operator+=` (Task 1).
- Produces: `StressEnergyContribution <Composite>::StressEnergy(...) const override` for each composite.

**Pattern:** each composite currently delegates its four methods to its two children, summing where both contribute and taking only one child where the other is physically zero (e.g. `idm_dr_idr_species.cpp:305` "IDM_DR has DeltaP == 0"). Collapse all four into one `StressEnergy` that sums the children's structs, starting from the **first** child (to match today's left-operand order):

```cpp
BaseSpecies::StressEnergyContribution IDM_DR_IDR_Species::StressEnergy(
    const BaseSpecies::PerturbLayout& base, const perturb_vector* pv,
    const double* y, const double* pvecback, const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se = idm_dr_->StressEnergy(my.idm_dr, pv, y, pvecback, ppw);
  se += idr_->StressEnergy(my.idr, pv, y, pvecback, ppw);
  return se;
}
```
Use each composite's existing layout-member names (`my.dcdm`/`my.dr`, `my.dncdm`/`my.dr`, `my.idm_drmd`/`my.idr_drmd`) and child-pointer names (read the current four methods in each file). This is byte-identical iff each cold child's `StressEnergy` returns literal `0.0` for the fields the old code omitted (`delta_p`, `rho_plus_p_shear`) — which holds by physics and is asserted by the existing code comments. Verify those child overrides (Tasks 6-7) return `0.` for those fields.

- [ ] **Step 1:** Add the override to each of the four composites (declaration in `.h`, definition in `.cpp`), removing nothing yet.
- [ ] **Step 2: Capture composite references, then verify byte-identical:**

```bash
# Build current branch first so refs are from this in-progress (still byte-identical) state? NO —
# references must be clean master. Capture them ONCE here from /tmp/ref_master-style runs of the
# committed master baseline. Easiest: check out the ini outputs from a clean class build of master.
git stash list   # confirm experiment still safe
./class test/scenarios/gauge_idmdr.ini  >/dev/null 2>&1
./class test/scenarios/gauge_dcdm.ini   >/dev/null 2>&1
```
Then compare each scenario's `cl`/`pk` (whatever roots those inis set — inspect with `grep '^root' test/scenarios/gauge_idmdr.ini`) against a master baseline captured the same way as Task 1. If composite master references were not captured up front, defer the strict cmp to Task 9's `TEST_LEVEL=2` (`test_idm_dr_idr_perturbations_match_reference`, `test_dcdm_dr_matches_reference`, `test_dncdm_dr_computes`), which compare these scenarios against `classyref`. Always also run VERIFY for the 1-massive-ν path.

Expected: 1-massive-ν `cmp` exit 0; composite scenarios match their references (here or in Task 9).

- [ ] **Step 3: Commit**

```bash
git add species/dcdm_dr_species.* species/dncdm_dr_species.* species/idm_dr_idr_species.* species/idm_drmd_idr_drmd_species.*
git commit -m "perturb: composite StressEnergy via struct sum of children"
```

---

### Task 9: Make `StressEnergy` pure; delete the four old virtuals + base default; full verification

**Files:**
- Modify: `species/base_species.h` (make `StressEnergy` pure, delete the four virtual declarations + the base default)
- Modify: every species `.h`/`.cpp` listed in Tasks 4-8 (delete the four method declarations + definitions)
- Modify: comment-only references in `species/base_species.h:70`, `species/composite_species.h:21`, `species/idr.h:122`, `species/ultra_relativistic.cpp:73`, `species/photons.cpp:164`

**Interfaces:**
- Produces: `virtual StressEnergyContribution BaseSpecies::StressEnergy(...) const = 0;` (pure); the four old virtuals no longer exist.

- [ ] **Step 1: Make the base virtual pure and delete the old four** in `species/base_species.h`: replace the scaffold default body with `= 0;`, and delete the `DeltaRho`/`RhoPlusPTheta`/`DeltaP`/`RhoPlusPShear` pure-virtual declarations (~lines 330-356).

- [ ] **Step 2: Delete the four method declarations + definitions from every species** updated in Tasks 4-8 (cdm, photons, baryons, ultra_relativistic, ncdm_species, dncdm_species, fluid, dcdm, dark_radiation_species, scalar_field, lambda, interacting_species + its four headers, and the four composites). The override `StressEnergy` added in those tasks is now the only stress-energy method each carries.

- [ ] **Step 3: Update the stale comments** (the five comment-only sites above) to refer to `StressEnergy` instead of the removed names.

- [ ] **Step 4: Build** — `cmake --build build -j 2>&1 | tail -10`. Expected: succeeds. A compile error here means a species or caller still references a removed method — fix that site (it should have been migrated in Tasks 3-8).

- [ ] **Step 5: Byte-identical + full suite + perf:**

```bash
# (a) byte-identical common path
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini >/dev/null 2>&1
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_cl_lensed.dat /tmp/ref_master/cl_lensed.dat
cmp output/base_2018_plikHM_TTTEEE_lowl_lowE_lensing_pk.dat        /tmp/ref_master/pk.dat
# (b) full regression suite incl. all composite reference comparisons
TEST_LEVEL=2 python -m pytest python/test_class.py -q 2>&1 | tail -20
# (c) perf: confirm the Perturbations stage improved vs master (no regression)
./class_profiled base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini 5 2>&1 | grep -i perturb
```
Expected: both `cmp` exit 0; pytest reports all pass / 0 failures; the Perturbations stage time is at/below master (target ≈ −6%).

- [ ] **Step 6: Commit**

```bash
git add -u species/ source/
git commit -m "perturb: make StressEnergy pure; remove the four legacy stress-energy virtuals"
```

---

## Self-Review

**Spec coverage:**
- Remove four virtuals → add one `StressEnergy` (pure): Tasks 1 (add) + 9 (remove/pure). ✓
- Struct + `operator+=`: Task 1. ✓
- Hot leaves direct overrides: Task 4. ✓
- NCDM fused single q-pass (and DNCDM): Task 5. ✓
- Composites via sum: Task 8. ✓
- Fluid + all remaining leaves: Tasks 6-7. ✓
- Shear pre-pass folded + #285 note: Task 2. ✓
- Cold callers (IC, ncdm output, fluid output): Task 3. ✓
- `Rho`/`P` retained: never removed (only the four perturbation virtuals are). ✓
- Bit-identity + composite scenario + `TEST_LEVEL=2` + perf: Task 9 (and per-task VERIFY). ✓
- Count-cache excluded: stated in Global Constraints; not in any task. ✓

**Placeholder scan:** the cold-species override bodies in Tasks 6-7 are specified as a deterministic transformation with a complete worked example (call-own-four-methods form), not "TBD"; the hot/ncdm bodies are verbatim or sourced from the named stash. Composite layout/child member names are to be read from each file (the four current methods name them) — the pattern code is complete.

**Type consistency:** `StressEnergyContribution` field names (`rho, p, delta_rho, rho_plus_p_theta, delta_p, rho_plus_p_shear`) and the `StressEnergy(layout, pv, y, pvecback, ppw)` signature are identical across Tasks 1-9. ✓

**Known soft spot (call out at execution):** the byte-identical *gates* for the cold/niche and composite species (Tasks 6-8) lean on `TEST_LEVEL=2` (Task 9) rather than a per-task `cmp`, because those species aren't in the 1-massive-ν benchmark. If a reviewer wants per-task gating, capture master references for `test/scenarios/gauge_idmdr.ini`, `test/scenarios/gauge_dcdm.ini`, and `test/dotsyntax_dncdm.ini` up front (same method as Task 1, Step 1) and add `cmp`s to those tasks.
