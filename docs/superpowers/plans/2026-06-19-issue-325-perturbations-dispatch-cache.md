# Perturbations Hot-Path Dispatch Cache Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove per-RHS string lookups, the deferred second pass, name-string filtering, and no-op derivs calls from the perturbations scalar/vector/tensor RHS loops, recovering the regression documented in issue #325.

**Architecture:** Each `perturb_vector` (single-mode by construction) gains a non-owning `active_species` view (`{species*, layout*}`, filtered to species that registered ≥1 perturbation variable in that mode) plus base-typed `photon_layout`/`baryon_layout` pointers — all non-owning views into the still-sole-owner `species_layouts`, built once per pv. The three RHS loops iterate `active_species`; ~18 hot photon/baryon reads switch from `all_species_.index_of(...)` to the cached pointers. PPF's "deferred" pass is deleted (its `dy[idx_Gamma]=Gamma_prime_fld` is valid in-loop because `ComputePpf` already ran in `perturb_einstein`).

**Tech Stack:** C++ (CLASS++), CMake, `./class` CLI, `class_profiled`, Python pytest (`TEST_LEVEL`).

**Correctness bar:** byte-identical output (verified by golden diff after every task), with the project's ~0.1% tolerance + `TEST_LEVEL=2` as the final gate.

**Design spec:** `docs/superpowers/specs/2026-06-19-issue-325-perturbations-dispatch-cache-design.md`

> **Line numbers are pre-change anchors.** All line numbers reflect the branch
> tip before Task 1. Earlier tasks insert/delete lines, so by Tasks 5–8 the
> targets have shifted (e.g. Task 4 adds ~10 lines above line 6018). **Locate by
> the quoted surrounding content or by grepping a distinctive token** (e.g.
> `PerturbVectorDerivs`, `perturb_tca_slip_and_shear`), not by the bare number.

---

## File Structure

| File | Change |
| --- | --- |
| `species/species_collection.h` | Add `photons_index_`/`baryons_index_` members + `photons_index()`/`baryons_index()` accessors. |
| `species/species_collection.cpp` | Set the two indices in `freeze()`. |
| `source/perturbations.h` | Add `ActiveSpecies` struct + `active_species`, `photon_layout`, `baryon_layout` to `perturb_vector`. |
| `source/perturbations_module.cpp` | Build the cache in `perturb_vector_init`; unified scalar/vector/tensor RHS loops; move recombination block; convert ~18 hot photon/baryon reads. |
| `species/base_species.h` | Delete `RequiresDeferredPerturbDerivs`. |
| `species/fluid.h` | Delete `RequiresDeferredPerturbDerivs` override. |

---

## Task 1: Baseline golden output + verify harness

No code change. Establishes the characterization baseline every later task diffs against. **Run on the current branch tip before touching any source.**

**Files:**
- Create: `/tmp/issue325/s1_ppf.ini`, `s2_fluid.ini`, `s3_ncdm.ini`, `s4_tensors.ini`
- Create: `/tmp/issue325_verify.sh`

- [ ] **Step 1: Write the four scenario inis** (cover PPF fluid, non-PPF fluid, massive NCDM, tensors/B-modes)

`/tmp/issue325/s1_ppf.ini`:
```
root = /tmp/issue325_out/s1_
output = tCl,pCl,lCl,mPk
l_max_scalars = 1200
P_k_max_1/Mpc = 1.0
Omega_Lambda = 0
fluid_equation_of_state = CLP
w0_fld = -0.9
wa_fld = 0.1
cs2_fld = 1.0
use_ppf = yes
```

`/tmp/issue325/s2_fluid.ini`: identical to `s1_ppf.ini` except `root = /tmp/issue325_out/s2_` and `use_ppf = no`.

`/tmp/issue325/s3_ncdm.ini`:
```
root = /tmp/issue325_out/s3_
output = tCl,pCl,lCl,mPk
l_max_scalars = 1200
P_k_max_1/Mpc = 1.0
N_ncdm = 1
m_ncdm = 0.06
```

`/tmp/issue325/s4_tensors.ini`:
```
root = /tmp/issue325_out/s4_
output = tCl,pCl,lCl
modes = s,t
l_max_scalars = 1200
l_max_tensors = 1000
r = 0.1
```

- [ ] **Step 2: Write the verify script**

`/tmp/issue325_verify.sh`:
```bash
#!/usr/bin/env bash
set -euo pipefail
cd /Users/au192734/Projects/class_claude
make class >/tmp/issue325_build.log 2>&1 || { echo "BUILD FAILED"; tail -30 /tmp/issue325_build.log; exit 1; }
rm -rf /tmp/issue325_out && mkdir -p /tmp/issue325_out
for s in s1_ppf s2_fluid s3_ncdm s4_tensors; do
  ./class /tmp/issue325/$s.ini >/dev/null
done
if [ "${1:-}" = "--save" ]; then
  rm -rf /tmp/issue325_golden && cp -r /tmp/issue325_out /tmp/issue325_golden
  echo "GOLDEN SAVED"
else
  diff -r /tmp/issue325_golden /tmp/issue325_out && echo "IDENTICAL" || { echo "DIFF FOUND"; exit 1; }
fi
```
Make executable: `chmod +x /tmp/issue325_verify.sh`

- [ ] **Step 3: Build + capture golden**

Run: `bash /tmp/issue325_verify.sh --save`
Expected: ends with `GOLDEN SAVED`, and `/tmp/issue325_golden/` contains `s1_*.dat … s4_*.dat`.

- [ ] **Step 4: Sanity re-run (no change ⇒ identical)**

Run: `bash /tmp/issue325_verify.sh`
Expected: `IDENTICAL`.

- [ ] **Step 5: Commit** (only the spec/plan are tracked; the harness lives in /tmp — nothing to commit here). Skip.

---

## Task 2: Cache photon/baryon indices in `SpeciesCollection`

**Files:**
- Modify: `species/species_collection.h`
- Modify: `species/species_collection.cpp`

- [ ] **Step 1: Add members + accessors**

In `species/species_collection.h`, after the baryon accessor block (the `baryons() const` returning `*baryons_;`, ~line 122), add:
```cpp
  std::size_t photons_index() const {
    assert(frozen_);
    return photons_index_;
  }
  std::size_t baryons_index() const {
    assert(frozen_);
    return baryons_index_;
  }
```
In the private section, after `BaseSpecies* baryons_ = nullptr;` (~line 159), add:
```cpp
  std::size_t photons_index_ = 0;
  std::size_t baryons_index_ = 0;
```

- [ ] **Step 2: Set the indices in `freeze()`**

In `species/species_collection.cpp`, in `freeze()`, after `baryons_ = baryons_slot->get();` and before `frozen_ = true;`, add:
```cpp
  for (std::size_t i = 0; i < species_.size(); ++i) {
    if (species_[i].key == "Photons")
      photons_index_ = i;
    if (species_[i].key == "Baryons")
      baryons_index_ = i;
  }
```
(Photons/Baryons presence is already guaranteed by the `photons_slot`/`baryons_slot` checks above.)

- [ ] **Step 3: Build**

Run: `make class`
Expected: builds clean (accessors unused so far — no behavior change).

- [ ] **Step 4: Verify byte-identical**

Run: `bash /tmp/issue325_verify.sh`
Expected: `IDENTICAL`.

- [ ] **Step 5: Commit**
```bash
git add species/species_collection.h species/species_collection.cpp
git commit -m "v4 #325: cache photons_index_/baryons_index_ in SpeciesCollection

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Add `ActiveSpecies` view + cached pointers to `perturb_vector`

**Files:**
- Modify: `source/perturbations.h:254`

- [ ] **Step 1: Add the struct + members**

In `source/perturbations.h`, immediately after the `species_layouts` declaration (line 254, `std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> species_layouts;`), add:
```cpp

  // Non-owning view into species_layouts: the species that registered ≥1
  // perturbation variable in THIS pv's mode (no-ops like Lambda excluded), in
  // lex-key order. Consumed by the scalar/vector/tensor RHS dispatch loops.
  struct ActiveSpecies {
    BaseSpecies*                species;
    BaseSpecies::PerturbLayout* layout;
  };
  std::vector<ActiveSpecies> active_species;

  // Always-present species, resolved once per pv (non-owning, into
  // species_layouts). Base-typed because perturbations.h cannot see the nested
  // PhotonsSpecies/BaryonsSpecies::PerturbLayout types (photons.h/baryons.h
  // include perturbations.h). Read sites static_cast to the concrete type.
  BaseSpecies::PerturbLayout* photon_layout = nullptr;
  BaseSpecies::PerturbLayout* baryon_layout = nullptr;
```

- [ ] **Step 2: Build**

Run: `make class`
Expected: builds clean (members unused so far).

- [ ] **Step 3: Verify byte-identical**

Run: `bash /tmp/issue325_verify.sh`
Expected: `IDENTICAL`.

- [ ] **Step 4: Commit**
```bash
git add source/perturbations.h
git commit -m "v4 #325: add ActiveSpecies view + photon/baryon layout pointers to perturb_vector

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Populate the cache in `perturb_vector_init`

Build the pointers and `active_species` but do **not** consume them yet (loops unchanged) — isolates the build step. Output stays identical because nothing reads the new state.

**Files:**
- Modify: `source/perturbations_module.cpp:2939` (after the `species_layouts` fill), and the three registration loops at `2950–2961` (scalars), `2985–2996` (vectors), `3012–3022` (tensors).

- [ ] **Step 1: Resolve the photon/baryon pointers**

In `source/perturbations_module.cpp`, immediately after the `species_layouts` fill loop (line 2939, the `}` closing the `push_back(...CreatePerturbLayout())` loop) and before `int index_pt = 0;` (line 2943), add:
```cpp

  // Hot-path: resolve the always-present photon/baryon layouts once per pv
  // (non-owning, base-typed views into species_layouts).
  ppv->photon_layout = ppv->species_layouts[all_species_.photons_index()].get();
  ppv->baryon_layout = ppv->species_layouts[all_species_.baryons_index()].get();
  assert(ppv->photon_layout && ppv->baryon_layout);
```
(If `<cassert>` is not already included in this TU, add `#include <cassert>` near the top with the other standard includes.)

- [ ] **Step 2: Populate `active_species` in the scalar registration loop**

Replace the scalar loop body (lines 2950–2961) with:
```cpp
    {
      size_t i = 0;
      for (auto& entry : all_species_) {
        const int before = index_pt;
        entry->RegisterPerturbationIndices(*ppv->species_layouts[i],
                                           ppv.get(),
                                           ppr,
                                           index_pt,
                                           ppw,
                                           static_cast<int>(ppt->gauge));
        if (index_pt > before)
          ppv->active_species.push_back({entry.get(), ppv->species_layouts[i].get()});
        ++i;
      }
    }
```

- [ ] **Step 3: Populate `active_species` in the vector registration loop**

Replace the vector loop body (lines 2985–2996) with:
```cpp
    {
      size_t i = 0;
      for (auto& entry : all_species_) {
        const int before = index_pt;
        entry->RegisterVectorPerturbationIndices(*ppv->species_layouts[i],
                                                 ppv.get(),
                                                 ppr,
                                                 index_pt,
                                                 ppw,
                                                 static_cast<int>(ppt->gauge));
        if (index_pt > before)
          ppv->active_species.push_back({entry.get(), ppv->species_layouts[i].get()});
        ++i;
      }
    }
```

- [ ] **Step 4: Populate `active_species` in the tensor registration loop**

Replace the tensor loop body (lines 3012–3022) with:
```cpp
    {
      size_t i = 0;
      for (auto& entry : all_species_) {
        const int before = index_pt;
        entry->RegisterTensorPerturbationIndices(*ppv->species_layouts[i],
                                                 ppv.get(),
                                                 ppr,
                                                 index_pt,
                                                 ppw,
                                                 static_cast<int>(ppt->gauge));
        if (index_pt > before)
          ppv->active_species.push_back({entry.get(), ppv->species_layouts[i].get()});
        ++i;
      }
    }
```

- [ ] **Step 5: Build (Debug, asserts on) + verify**

Run: `make class && bash /tmp/issue325_verify.sh`
Expected: builds clean, asserts pass, `IDENTICAL` (new state is built but unused).

- [ ] **Step 6: Commit**
```bash
git add source/perturbations_module.cpp
git commit -m "v4 #325: build per-pv active_species view + photon/baryon pointers

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Unified scalar RHS loop + move recombination block

Replace the explicit baryon/photon calls, both passes, and reposition the perturbed-recombination derivative block (it reads `dy[b_dr_lay.idx_delta]`).

**Files:**
- Modify: `source/perturbations_module.cpp:6018–6085`

- [ ] **Step 1: Replace the region**

Replace lines 6018–6085 (from `if (ppw->approx[ppw->index_ap_tca] == (int) tca_on) {` through the closing of the newtonian-gauge metric `if`) with:
```cpp
    if (ppw->approx[ppw->index_ap_tca] == (int) tca_on) {
      perturb_tca_slip_and_shear(y, pppaw);
    }

    /* Species contributions to the scalar perturbation ODE — single pass over
       the active (non-no-op) species in lex order. Photons/baryons are ordinary
       entries; the PPF fluid too (its dy[idx_Gamma] = Gamma_prime_fld is valid
       here because ComputePpf already ran in perturb_einstein() above). Order is
       free: each species reads scalar_ctx/y and writes only its own dy. */
    for (const auto& [species, layout] : pv->active_species)
      species->PerturbDerivs(*layout, tau, y, dy, *pppaw);

    /* perturbed recombination — derivatives of delta x_e and delta T_b. Runs
       after the species loop because it reads dy[baryon.idx_delta]. */
    if ((ppt->has_perturbed_recombination) && (ppw->approx[ppw->index_ap_tca] == (int) tca_off)) {
      /* alpha * n_H is in inverse seconds, so we have to multiply it by Mpc_in_sec */
      dy[ppw->pv->index_pt_perturbed_recombination_delta_chi] =
          -alpha_rec * a * chi * n_H * (delta_alpha_rec + delta_chi + delta_b) * _Mpc_over_m_ / _c_;

      /* see the documentation for this formula */
      dy[ppw->pv->index_pt_perturbed_recombination_delta_temp] =
          2. / 3. * dy[b_dr_lay.idx_delta] -
          a * Compton_CR * pow(pba->T_cmb / a, 4) * chi / (1. + chi + fHe) *
              ((1. -
                pba->T_cmb * pba->a_today / a / pvecthermo[thermodynamics_module_->index_th_Tb_]) *
                   (delta_g + delta_chi * (1. + fHe) / (1. + chi + fHe)) +
               pba->T_cmb * pba->a_today / a / pvecthermo[thermodynamics_module_->index_th_Tb_] *
                   (delta_temp - 1. / 4. * delta_g));
    }

    /** - --> metric */

    if (ppt->gauge == possible_gauges::synchronous) {
      dy[pv->index_pt_eta] = pvecmetric[ppw->index_mt_eta_prime];
    }

    if (ppt->gauge == possible_gauges::newtonian) {
      dy[pv->index_pt_phi] = pvecmetric[ppw->index_mt_phi_prime];
    }
```

This removes: the two explicit `baryons()/photons()` `PerturbDerivs` blocks, the old `/* perturbed recombination */` comment + block at its original position, the first generic loop (with `sp->name()` filtering), and the deferred second loop. `b_dr_lay`, `delta_g`, `delta_b`, `delta_temp`, `delta_chi`, `chi`, `n_H`, `fHe`, `alpha_rec`, `delta_alpha_rec`, `a` are all defined earlier in the `_scalars_` block and remain in scope.

- [ ] **Step 2: Build + verify byte-identical**

Run: `make class && bash /tmp/issue325_verify.sh`
Expected: `IDENTICAL`. (Same calls, same lex order — baryons before photons; recombination still reads baryon's freshly-written dy; fluid's `dy[idx_Gamma]` is position-independent.)

- [ ] **Step 3: Commit**
```bash
git add source/perturbations_module.cpp
git commit -m "v4 #325: unify scalar RHS loop over active_species; drop deferred pass

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Unified vector & tensor RHS loops

**Files:**
- Modify: `source/perturbations_module.cpp` — vector loop `6118–6125`, tensor loop `6139–6146`

- [ ] **Step 1: Replace the vector loop**

Replace lines 6118–6125 with:
```cpp
    /** - --> species Boltzmann hierarchies (vector modes) */
    for (const auto& [species, layout] : pv->active_species)
      species->PerturbVectorDerivs(*layout, tau, y, dy, *pppaw);
```

- [ ] **Step 2: Replace the tensor loop**

Replace lines 6139–6146 (the `{ size_t i = 0; for (auto& sp : all_species_) { sp->PerturbTensorDerivs(...); ++i; } }` block) with:
```cpp
    /** - --> species Boltzmann hierarchies (tensor modes) */
    for (const auto& [species, layout] : pv->active_species)
      species->PerturbTensorDerivs(*layout, tau, y, dy, *pppaw);
```

- [ ] **Step 3: Build + verify byte-identical**

Run: `make class && bash /tmp/issue325_verify.sh`
Expected: `IDENTICAL` (tensor scenario `s4` exercises this; the previously-iterated no-op species write nothing).

- [ ] **Step 4: Commit**
```bash
git add source/perturbations_module.cpp
git commit -m "v4 #325: unify vector/tensor RHS loops over active_species

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: Delete `RequiresDeferredPerturbDerivs`

Now unused. Removing it is what makes the deferral deletion complete.

**Files:**
- Modify: `species/base_species.h:256–262`
- Modify: `species/fluid.h:83–85`

- [ ] **Step 1: Delete the base declaration**

In `species/base_species.h`, delete the block (lines 256–262):
```cpp
  /**
   * Returns true if this species' PerturbDerivs must run AFTER all other
   * species in a second pass. Used for PPF fluid (FluidSpecies).
   */
  virtual bool RequiresDeferredPerturbDerivs() const {
    return false;
  }
```

- [ ] **Step 2: Delete the fluid override**

In `species/fluid.h`, delete (lines 83–85):
```cpp
  bool RequiresDeferredPerturbDerivs() const override {
    return true;
  }
```

- [ ] **Step 3: Confirm no remaining references**

Run: `grep -rn RequiresDeferredPerturbDerivs species/ source/`
Expected: no matches.

- [ ] **Step 4: Build + verify byte-identical**

Run: `make class && bash /tmp/issue325_verify.sh`
Expected: builds clean, `IDENTICAL`.

- [ ] **Step 5: Commit**
```bash
git add species/base_species.h species/fluid.h
git commit -m "v4 #325: delete now-unused RequiresDeferredPerturbDerivs

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8a: Convert hot photon/baryon reads — einstein, stress-energy, sources

The mechanical sweep that delivers the measured win. **Conversion rule** for each hot site:

> Delete the line `const size_t <X>_i = all_species_.index_of("Photons");` (or `"Baryons"`), and in the following `static_cast<…PerturbLayout&>(*<H>->species_layouts[<X>_i])` replace `<H>->species_layouts[<X>_i]` with `<H>->photon_layout` (Photons) / `<H>->baryon_layout` (Baryons), keeping the existing pv handle `<H>` (`pv` or `ppw->pv`).
>
> Edge case: if `<X>_i` is used for anything other than that single `species_layouts[<X>_i]` deref, instead replace it with `all_species_.photons_index()` / `baryons_index()` (cached, no string search) and keep the variable.

**Worked example** (the scalar-derivs setup, applied in Task 8b):
```cpp
// before
const size_t g_dr_i  = all_species_.index_of("Photons");
const auto& g_dr_lay = static_cast<const PhotonsSpecies::PerturbLayout&>(
    *pv->species_layouts[g_dr_i]);
// after
const auto& g_dr_lay = static_cast<const PhotonsSpecies::PerturbLayout&>(
    *pv->photon_layout);
```

**Files:**
- Modify: `source/perturbations_module.cpp` — sites at lines `4663` (`g_ein_i`, einstein), `4795` (`g_tse_i`), `4834` (`b_tse_pre_i`), `4880` (`g_tse2_i`), `4881` (`b_tse_i`) (total_stress_energy), `5035` (`g_vec_i`), `5068` (`g_tens_i`) (einstein vec/tens), `5201` (`g_src_i`), `5204` (`b_src_i`), `5502` (`g_tsrc_i`) (sources).

- [ ] **Step 1: Locate the sites**

Run: `grep -n 'index_of("Photons")\|index_of("Baryons")' source/perturbations_module.cpp`
Expected: the lines above plus the cold construction sites (`~3209–3527`) and the Task-8b sites — convert only the Task-8a list in this task.

- [ ] **Step 2: Apply the conversion rule to each Task-8a site**

Photons sites (`g_ein_i`, `g_tse_i`, `g_tse2_i`, `g_vec_i`, `g_tens_i`, `g_src_i`, `g_tsrc_i`) → `*<H>->photon_layout`. Baryons sites (`b_tse_pre_i`, `b_tse_i`, `b_src_i`) → `*<H>->baryon_layout`. Use the exact pv handle already present on each line.

- [ ] **Step 3: Build + verify byte-identical**

Run: `make class && bash /tmp/issue325_verify.sh`
Expected: builds clean, `IDENTICAL`.

- [ ] **Step 4: Commit**
```bash
git add source/perturbations_module.cpp
git commit -m "v4 #325: hot photon/baryon reads via cached pointers (einstein/stress-energy/sources)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8b: Convert hot photon/baryon reads — derivs setup, tca, rsa, vector

Same conversion rule as Task 8a.

**Files:**
- Modify: `source/perturbations_module.cpp` — sites at lines `5669` (`g_pv_i`), `5889` (`g_dr_i`), `5897` (`b_dr_i`) (scalar derivs setup), `6095` (`g_vec_dr_i`), `6099` (`b_vec_dr_i`) (vector derivs), `6226` (`g_tca_i`), `6229` (`b_tca_i`) (tca slip/shear), `6513` (`b_rsa_i`) (rsa).

- [ ] **Step 1: Apply the conversion rule to each Task-8b site**

Photons sites (`g_pv_i`, `g_dr_i`, `g_vec_dr_i`, `g_tca_i`) → `*<H>->photon_layout`. Baryons sites (`b_dr_i`, `b_vec_dr_i`, `b_tca_i`, `b_rsa_i`) → `*<H>->baryon_layout`. Use the exact pv handle already present on each line.

- [ ] **Step 2: Confirm only cold construction sites remain on `index_of`**

Run: `grep -n 'index_of("Photons")\|index_of("Baryons")' source/perturbations_module.cpp`
Expected: only the `CopyPerturbationsAcrossSwitch`/construction sites in the `~3209–3527` region remain (cold, per-approximation-switch — intentionally left).

- [ ] **Step 3: Build + verify byte-identical**

Run: `make class && bash /tmp/issue325_verify.sh`
Expected: builds clean, `IDENTICAL`.

- [ ] **Step 4: Commit**
```bash
git add source/perturbations_module.cpp
git commit -m "v4 #325: hot photon/baryon reads via cached pointers (derivs/tca/rsa/vector)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 9: Full verification + re-profile

**Files:** none (validation only).

- [ ] **Step 1: C++ unit tests**

Run: `make test`
Expected: all CTest targets pass.

- [ ] **Step 2: Full `TEST_LEVEL=2` suite vs classyref**

Build the Python wrapper and run the level-2 suite (classyref must be the post-merge master build per project convention):
```bash
pip install --no-build-isolation . >/tmp/issue325_pipinstall.log 2>&1
cd python && TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest test_class.py -q
```
Expected: all pass (report the pass/fail count, e.g. `2538/0`).

- [ ] **Step 3: Re-profile vs the issue baseline**

Build `class_profiled` and run the issue's Planck-2018 single-NCDM benchmark single-threaded; take the perturbations-module median over 5 loops:
```bash
cmake --build build/cmake --target class_profiled --parallel
OMP_NUM_THREADS=1 ./class_profiled base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini
```
(Set `threads = 1` in the ini if not already.) Expected: perturbations median materially below the `5633 ms` baseline and at/under the `~5075 ms` cached-index experiment — record the number.

- [ ] **Step 4: Report**

Summarize for the PR: byte-identical golden diff across the four scenarios, `make test` green, `TEST_LEVEL=2` count, and the before/after perturbations-median profile. No commit (validation only).

---

## Notes / risks

- **Asserts:** Tasks 4–8 should be built at least once with asserts enabled (Debug / non-`NDEBUG`) so the `frozen_`, `photon_layout`/`baryon_layout` non-null, and `index_of` debug checks actually fire. Per [[project_cmake_build_system]], release builds can render NDEBUG asserts vacuous.
- **Never `git add -A`** in this repo (in-source CMake/Xcode artifacts) — each commit lists explicit paths. See [[project_v4_design_review_issues]].
- **If any golden diff fails:** stop and debug that task in isolation before proceeding; the change is designed to be byte-identical, so a diff is a real defect (use `diff /tmp/issue325_golden/sN_cl.dat /tmp/issue325_out/sN_cl.dat` to localize).
- **classyref staleness:** regenerated against post-merge master per [[project_v4_design_review_issues]]; since output is byte-identical, no regeneration is expected to be needed.
