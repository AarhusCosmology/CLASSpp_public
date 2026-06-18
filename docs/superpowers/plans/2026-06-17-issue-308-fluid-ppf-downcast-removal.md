# #308 Fluid PPF Downcast Removal — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the fluid's module-level `dynamic_cast`/`static_cast<FluidSpecies>` downcasts by making the fluid an ordinary species everywhere except the genuine PPF modified-gravity coupling, which is reached through one module-local `ppf_fluid()` accessor.

**Architecture:** (1) Fluid self-computes `w(a)` in `ComputeBackground` and rejoins the main background loop. (2) `PerturbationsModule` caches a `FluidSpecies* ppf_fluid_` (one downcast in its ctor; non-null iff a fluid exists and `use_ppf()`). (3) The perturbation stress-energy loops drop the fluid skip — a PPF fluid self-zeroes through the generic `Delta`/`Theta`/`DeltaP` interface (its `idx_delta`/`idx_theta` are unregistered), a non-PPF fluid contributes its real values; the only fluid-specific handling left is a small post-loop "PPF extra" calling `ppf_fluid()->ComputePpf(...)`, and `ComputePpf` self-subtracts its own ρ+p from the now-inclusive `rho_plus_p_tot`.

**Tech Stack:** C++17 (CLASSpp), Python/Cython `classy` bindings, pytest goldens, CMake build (`pip install . --no-build-isolation`).

**Policy:** Clean code over bit-identical (project policy). Output changes at the ≤0.1% level by design; per-task gate is a `classy`-vs-frozen-`classyref` spot-check ≤0.1%, and the exact goldens + `classyref` are regenerated at the end (Task 7).

**Spec:** `docs/superpowers/specs/2026-06-17-issue-308-fluid-ppf-downcast-removal-design.md`

**Branch:** `issue-308-fluid-ppf-downcast-removal` (already created).

---

### Task 0: Verification harness

**Files:**
- Create: `/tmp/spotcheck_308.py`

The frozen `classyref` (current pre-slice master build, already installed) is the
oracle. This script measures `classy` (this branch, rebuilt each task) against it.
It is the per-task gate. Do **not** reinstall `classyref` until Task 7.

- [ ] **Step 1: Create the spot-check script**

```python
# /tmp/spotcheck_308.py — classy (this branch) vs frozen classyref, expect <= 0.1%
import numpy as np
import classy
import classyref

SCENARIOS = {
    "lcdm": {},
    "ppf_fluid":    {"Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.1, "use_ppf": "yes"},
    "nonppf_fluid": {"Omega_Lambda": 0.0, "w0_fld": -0.9, "wa_fld": 0.0, "use_ppf": "no"},
}

def run(Mod, extra, gauge):
    p = {"output": "tCl,pCl,lCl,mPk", "lensing": "yes", "gauge": gauge,
         "l_max_scalars": 2000, "P_k_max_1/Mpc": 5.0, "z_pk": 0.0}
    p.update(extra)
    c = Mod.Class(); c.set(p); c.compute()
    tt = c.lensed_cl(2000)["tt"].copy()
    kk = np.logspace(-3, np.log10(5.0), 200)
    pk = np.array([c.pk(k, 0.0) for k in kk])
    c.struct_cleanup(); c.empty()
    return tt, pk

def relmax(a, b):
    a, b = np.asarray(a), np.asarray(b)
    denom = (np.abs(a) + np.abs(b)) / 2.0
    m = denom > 0
    return float(np.max(np.abs(a - b)[m] / denom[m])) if m.any() else 0.0

worst = 0.0
for name, extra in SCENARIOS.items():
    for gauge in ("synchronous", "newtonian"):
        tt_n, pk_n = run(classy, extra, gauge)
        tt_r, pk_r = run(classyref, extra, gauge)
        rtt, rpk = relmax(tt_n, tt_r), relmax(pk_n, pk_r)
        worst = max(worst, rtt, rpk)
        print(f"{name:14s} {gauge:11s}  relmax TT={rtt:.2e}  Pk={rpk:.2e}")
print(f"\nWORST = {worst:.2e}  ({'PASS <=1e-3' if worst <= 1e-3 else 'FAIL'})")
```

- [ ] **Step 2: Baseline run (before any code change)**

Run: `cd /Users/au192734/Projects/class_claude && python /tmp/spotcheck_308.py`
Expected: every line `TT=0.00e+00  Pk=0.00e+00`, `WORST = 0.00e+00 (PASS ...)`.
(Confirms `classy == classyref` at the start; the harness works.)

---

### Task 1: Fluid self-computes `w` in `ComputeBackground`; drop `WriteWFld` and the dead deferral hook

**Files:**
- Modify: `species/fluid.cpp` (`ComputeBackground`, remove `WriteWFld` def)
- Modify: `species/fluid.h` (remove `WriteWFld` decl, remove `RequiresDeferredBackground` override)
- Modify: `species/base_species.h` (remove dead `RequiresDeferredBackground` virtual)
- Modify: `source/background_module.cpp` (`background_functions()` — delete fluid special-case)

- [ ] **Step 1: Make `ComputeBackground` self-compute `w`**

In `species/fluid.cpp`, replace the body of `ComputeBackground` (currently fluid.cpp:63–67):

```cpp
void FluidSpecies::ComputeBackground(double a_rel, const double* pvecback_B, double* pvecback) {
  // a == a_rel (a_today == 1; see project note). Compute w(a) ourselves rather
  // than having the module pre-fill the slots.
  double w_fld, dw_over_da_fld, integral_fld;
  ComputeWFld(a_rel, &w_fld, &dw_over_da_fld, &integral_fld);
  pvecback[index_bg_w_fld_]          = w_fld;
  pvecback[index_bg_dw_over_da_fld_] = dw_over_da_fld;
  pvecback[index_bg_rho_fld_]        = pvecback_B[index_bi_rho_fld_];
}
```

- [ ] **Step 2: Remove `WriteWFld`**

In `species/fluid.cpp`, delete the whole `WriteWFld` definition (fluid.cpp:107–110):

```cpp
void FluidSpecies::WriteWFld(double w_fld, double dw_over_da_fld, double* pvecback) const {
  pvecback[index_bg_w_fld_]          = w_fld;
  pvecback[index_bg_dw_over_da_fld_] = dw_over_da_fld;
}
```

In `species/fluid.h`, delete its declaration (fluid.h:93–96, the Doxygen comment + signature):

```cpp
  /** Called by BackgroundModule::background_functions() before ComputeBackground().
   ... */
  void WriteWFld(double w_fld, double dw_over_da_fld, double* pvecback) const;
```

- [ ] **Step 3: Remove the dead `RequiresDeferredBackground`**

In `species/fluid.h`, delete the override (fluid.h:83–88 region):

```cpp
  bool RequiresDeferredBackground() const override {
    return true;
  }
```

In `species/base_species.h`, delete the base virtual + its doc comment (base_species.h:264–270 region):

```cpp
  /**
   * Returns true if this species' ComputeBackground must be deferred.
   * Used for FluidSpecies which needs w_fld evaluated before it can run.
   */
  virtual bool RequiresDeferredBackground() const {
    return false;
  }
```

- [ ] **Step 4: Delete the fluid special-case in `background_functions()`**

In `source/background_module.cpp`, replace the loop + fluid block (background_module.cpp:351–369) with a single unconditional loop:

```cpp
  /* Compute background for all species. Each species writes its own pvecback
     slots; the Fluid computes its own w(a) in ComputeBackground. */
  for (const auto& [name, sp] : all_species_) {
    sp->ComputeBackground(a_rel, pvecback_B, pvecback);
    accumulate(*sp);
  }
```

(Removes the `if (name == "Fluid") continue;`, the `background_w_fld(a, …)` call,
the `static_cast<FluidSpecies&>(…).WriteWFld(…)`, and the separate fluid block.
`background_w_fld` itself is untouched.)

- [ ] **Step 5: Build**

Run: `cd /Users/au192734/Projects/class_claude && pip install . --no-build-isolation 2>&1 | tail -5`
Expected: `Successfully installed classy-community-...` (no compile errors).

- [ ] **Step 6: Verify ≤0.1%**

Run: `python /tmp/spotcheck_308.py`
Expected: all scenarios `WORST <= 1e-3` → `PASS`. (Non-fluid scenarios stay
`0.00e+00`; the fluid scenarios shift at ULP because the fluid now accumulates
into `rho_tot` at its lex position instead of last.)

- [ ] **Step 7: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add species/fluid.cpp species/fluid.h species/base_species.h source/background_module.cpp
git commit -m "v4 #308: fluid self-computes w; rejoin main background loop

ComputeBackground now calls ComputeWFld itself; deletes the WriteWFld
downcast in background_functions() and the dead RequiresDeferredBackground
hook. Fluid is an ordinary species in the background loop.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 2: Add the `ppf_fluid()` accessor to `PerturbationsModule`

**Files:**
- Modify: `source/perturbations_module.h` (forward-decl, accessor, member)
- Modify: `source/perturbations_module.cpp` (ctor resolution)

- [ ] **Step 1: Forward-declare `FluidSpecies` and add accessor + member to the header**

In `source/perturbations_module.h`, add a forward declaration after the includes
(after line 7), before `class PerturbationsModule`:

```cpp
class FluidSpecies;
```

Add the public accessor in the `public:` section (e.g. just after the opening
`public:` at line 10):

```cpp
  // The single PPF fluid, or nullptr. PPF is a modified-gravity closure defined
  // relative to the whole universe (intrinsically singular), so the module owns
  // it; non-PPF fluids are ordinary species. Resolved once in the ctor.
  FluidSpecies* ppf_fluid() const { return ppf_fluid_; }
```

Add the member in the `private:` section (after line 183):

```cpp
  FluidSpecies* ppf_fluid_ = nullptr;
```

- [ ] **Step 2: Resolve `ppf_fluid_` in the constructor**

In `source/perturbations_module.cpp`, in the constructor (perturbations_module.cpp:61–71),
insert the resolution after the `SetThermodynamicsModule`/`SetPerturbs` loop and
**before** `perturb_init();`. `fluid.h` is already included (line 35).

```cpp
  /* Cache the single PPF fluid (the one downcast, localized here). use_ppf is a
     construction-fixed input; at most one "Fluid" can exist. */
  if (auto* p = all_species_.find("Fluid")) {
    auto* f = static_cast<FluidSpecies*>(p->get());
    if (f->use_ppf())
      ppf_fluid_ = f;
  }
```

- [ ] **Step 3: Build**

Run: `cd /Users/au192734/Projects/class_claude && pip install . --no-build-isolation 2>&1 | tail -5`
Expected: `Successfully installed ...` (the accessor is unused so far; this only
proves it compiles and links).

- [ ] **Step 4: Verify unchanged**

Run: `python /tmp/spotcheck_308.py`
Expected: identical to Task 1 Step 6 (pure addition, no behavior change).

- [ ] **Step 5: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add source/perturbations_module.h source/perturbations_module.cpp
git commit -m "v4 #308: add module-local ppf_fluid() accessor

Caches the single PPF fluid (one downcast in the ctor). Non-null iff a
fluid exists and use_ppf(); encodes both invariants for call sites.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 3: `FillSources`/`PrintVariables` compute non-PPF fluid quantities locally

This removes the non-PPF fluid's dependence on the `ppw->*_fld` workspace handoff,
so Task 5 can stop populating it for non-PPF. Behavior-preserving here (the
dedicated block still populates `*_fld`; these methods just stop reading it for
non-PPF).

**Files:**
- Modify: `species/fluid.cpp` (`FillSources`, `PrintVariables`)

- [ ] **Step 1: Use the layout in `FillSources` for non-PPF**

In `species/fluid.cpp`, change the `FillSources` signature to bind the layout
(currently `const BaseSpecies::PerturbLayout& /*layout*/` at fluid.cpp:179) and
source the delta/theta from `y` for non-PPF. Replace the two `SetSourceValue`
blocks (fluid.cpp:189–198) so each chooses its base value by `use_ppf_`:

```cpp
void FluidSpecies::FillSources(const BaseSpecies::PerturbLayout& base,
                               const double* y,
                               const double* /*dy*/,
                               PerturbSourceContext& ctx) {
  const auto& layout          = static_cast<const PerturbLayout&>(base);
  PerturbationsModule* p_mod  = ctx.p_mod;
  perturb_workspace* ppw      = ctx.ppw;
  const double* pvecback      = ppw->pvecback;

  // Fluid sources are scalar-only
  if (ctx.index_md != p_mod->index_md_scalars_)
    return;

  // ── delta_fld ─────────────────────────────────────────────────────────────
  if (index_tp_delta_ >= 0) {
    const double w_fld = W(pvecback);
    // PPF: delta_rho_fld is published by ComputePpf. Non-PPF: delta = y[idx_delta].
    const double delta_fld =
        use_ppf_ ? ppw->delta_rho_fld / Rho(pvecback) : y[layout.idx_delta];
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp_delta_,
                          ctx.index_tau, ctx.index_k,
                          delta_fld + 3. * ctx.a_prime_over_a * (1. + w_fld) *
                                          ctx.theta_over_k2);  // N-body gauge correction
  }

  // ── theta_fld ─────────────────────────────────────────────────────────────
  if (index_tp_theta_ >= 0) {
    const double w_fld = W(pvecback);
    // PPF: rho_plus_p_theta_fld is published by ComputePpf. Non-PPF: theta = y[idx_theta].
    const double theta_fld =
        use_ppf_ ? ppw->rho_plus_p_theta_fld / (1. + w_fld) / Rho(pvecback)
                 : y[layout.idx_theta];
    p_mod->SetSourceValue(ctx.index_md, ctx.index_ic, index_tp_theta_,
                          ctx.index_tau, ctx.index_k,
                          theta_fld + ctx.theta_shift);  // N-body gauge correction
  }
}
```

- [ ] **Step 2: Compute non-PPF quantities locally in `PrintVariables`**

In `species/fluid.cpp`, replace the value-capture block of `PrintVariables`
(fluid.cpp:245–249) so non-PPF derives the three diagnostics from the generic
interface instead of `ppw->*_fld`:

```cpp
  if (!w.IsTitleMode()) {
    if (use_ppf_) {
      delta_rho_fld        = ppw->delta_rho_fld;
      rho_plus_p_theta_fld = ppw->rho_plus_p_theta_fld;
      delta_p_fld          = ppw->delta_p_fld;
    }
    else {
      const auto& layout =
          static_cast<const PerturbLayout&>(*ppw->pv->species_layouts[collection_index_]);
      const double rho = Rho(ppw->pvecback);
      delta_rho_fld        = rho * Delta(layout, ppw->pv.get(), y, ppw->pvecback, ppw);
      rho_plus_p_theta_fld = (rho + P(ppw->pvecback)) *
                             Theta(layout, ppw->pv.get(), y, ppw->pvecback, ppw);
      delta_p_fld          = DeltaP(layout, ppw->pv.get(), y, ppw->pvecback, ppw);
    }
  }
```

- [ ] **Step 3: Build**

Run: `cd /Users/au192734/Projects/class_claude && pip install . --no-build-isolation 2>&1 | tail -5`
Expected: `Successfully installed ...`.

- [ ] **Step 4: Verify ≤0.1%**

Run: `python /tmp/spotcheck_308.py`
Expected: `WORST <= 1e-3`. The non-PPF fluid transfer source now uses
`y[idx_delta]` directly instead of `Rho·y[idx_delta]/Rho` (a round-trip), so the
`nonppf_fluid` numbers may move at ULP; PPF path is unchanged.

- [ ] **Step 5: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add species/fluid.cpp
git commit -m "v4 #308: FillSources/PrintVariables compute non-PPF fluid locally

Non-PPF fluid no longer reads the ppw->*_fld workspace handoff (which
becomes PPF-only). Sources/diagnostics derive from y / the generic
interface. Prep for dropping the dedicated stress-energy block.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 4: Site 5 — IC validity tests use `ppf_fluid()`

**Files:**
- Modify: `source/perturbations_module.cpp` (the fluid IC-validity block, ~514–535)

- [ ] **Step 1: Replace the downcast + `use_ppf()` with `ppf_fluid()`**

In `source/perturbations_module.cpp`, in the fluid IC-validity block
(perturbations_module.cpp:514–535), remove the `static_cast<const FluidSpecies&> fluid`
local and change the inner guard from `!fluid.use_ppf()` to `!ppf_fluid()`:

```cpp
  if (all_species_.count("Fluid")) {
    /* check values of w_fld at initial time and today */
    background_module_->background_w_fld(0., &w_fld_ini, &dw_over_da_fld, &integral_fld);
    background_module_->background_w_fld(pba->a_today, &w_fld_0, &dw_over_da_fld, &integral_fld);

    class_test(w_fld_ini >= 0.,
               "The fluid is meant to be negligible at early time, and unimportant for defining "
               "the initial conditions of other species. You are using parameters for which this "
               "assumption may break down, since at early times you have w_fld(a--->0) = %e >= 0",
               w_fld_ini);

    if (!ppf_fluid()) {  // a non-PPF fluid is present (single-fluid: count && !ppf)
      class_test((w_fld_ini + 1.0) * (w_fld_0 + 1.0) <= 0.0,
                 "w crosses -1 between the infinite past and today, and this would lead to "
                 "divergent perturbation equations for the fluid perturbations. Try to switch to "
                 "PPF scheme: use_ppf = yes");
```

(Leave the rest of that `if (!ppf_fluid())` body — the `(w0,wa)=(-1,0)` test —
unchanged; only the `static_cast` local and the `if (!fluid.use_ppf())` head change.)

- [ ] **Step 2: Build**

Run: `cd /Users/au192734/Projects/class_claude && pip install . --no-build-isolation 2>&1 | tail -5`
Expected: `Successfully installed ...`.

- [ ] **Step 3: Verify ≤0.1% (validation-only change)**

Run: `python /tmp/spotcheck_308.py`
Expected: identical to Task 3 (these are `class_test` guards; no numeric output
change). Confirm both `ppf_fluid` and `nonppf_fluid` scenarios still run (the
non-PPF guard still fires correctly).

- [ ] **Step 4: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add source/perturbations_module.cpp
git commit -m "v4 #308: fluid IC-validity tests gate on ppf_fluid() not a downcast

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 5: Unconditional stress-energy loops + PPF-extra + `ComputePpf` self-subtract

The core change. The fluid stops being skipped; the dedicated block collapses to
the PPF-extra; `ComputePpf` subtracts its own ρ+p from the now-inclusive
`rho_plus_p_tot`.

**Files:**
- Modify: `source/perturbations_module.cpp` (shear pass ~4888–4903; pass-2 skip ~4914–4920; dedicated block ~4960–5004)
- Modify: `species/fluid.cpp` (`ComputePpf` — rest-total)

- [ ] **Step 1: `ComputePpf` computes the rest-of-universe ρ+p total**

In `species/fluid.cpp`, in `ComputePpf`, after `w_fld` is computed (after
fluid.cpp:423, `const double w_prime_fld = ...`), add:

```cpp
  // The PPF closure is defined relative to the rest of the universe. The module's
  // runtime rho_plus_p_tot now INCLUDES this fluid (it flows through the generic
  // loop), so subtract our own background rho+p here.
  const double rho_plus_p_tot_rest =
      ppw->rho_plus_p_tot - Rho(ppw->pvecback) * (1. + w_fld);
```

Then replace **every** `ppw->rho_plus_p_tot` inside `ComputePpf` with
`rho_plus_p_tot_rest`. There are five occurrences:
- fluid.cpp:458 (`S_fld` denominator)
- fluid.cpp:470 (`rho_plus_p_theta_fld`, first ratio)
- fluid.cpp:473 (`rho_plus_p_theta_fld`, `1 + 4.5·a2/k2/s2sq·…`)
- fluid.cpp:500 (`theta_t`)
- fluid.cpp:503 (`theta_t_prime` divisor)

Leave `ppw->rho_plus_p_theta`, `ppw->delta_p`, `ppw->rho_plus_p_shear`, and the
background-table `rho_t`/`p_t` (fluid.cpp:485–486) **unchanged** — the PPF fluid
contributes 0 to those accumulators, and `rho_t` already subtracts `rho_fld` from
the background table.

- [ ] **Step 2: Shear pass skips only Photons/Baryons**

In `source/perturbations_module.cpp`, replace the shear pass
(perturbations_module.cpp:4888–4903):

```cpp
    /* Pass 1: accumulate rho_plus_p_shear over all sectors except photons and
       baryons (handled in the dedicated block above). Both PPF (always 0) and
       non-PPF fluid contribute 0 from RhoPlusPShear, so no fluid skip needed. */
    {
      size_t i = 0;
      for (const auto& sp : all_species_) {
        if (sp->name() != "Photons" && sp->name() != "Baryons")
          ppw->rho_plus_p_shear +=
              sp->RhoPlusPShear(*ppw->pv->species_layouts[i], ppw->pv.get(), y, ppw->pvecback, ppw);
        ++i;
      }
    }
```

- [ ] **Step 3: Pass 2 skips only Photons/Baryons**

In `source/perturbations_module.cpp`, in pass 2 (perturbations_module.cpp:4916–4920),
change the skip condition to drop `Fluid`:

```cpp
      for (const auto& sp : all_species_) {
        if (sp->name() == "Photons" || sp->name() == "Baryons") {
          ++i;
          continue;  // Photons/Baryons have dedicated paths above
        }
```

(The fluid now flows through: non-PPF contributes real `Delta`/`Theta`/`DeltaP`;
PPF contributes 0 to those but its real `rho + P()` to `rho_plus_p_tot`. The
matter-tally guard `if (sp->IsMatterSpecies())` already excludes the DarkEnergy
fluid.)

- [ ] **Step 4: Replace the dedicated fluid block with the PPF-extra**

In `source/perturbations_module.cpp`, replace the entire `/* fluid contribution */`
block (perturbations_module.cpp:4967–5002) and delete the trailing
"fluid must be last" comment (perturbations_module.cpp:5004) with:

```cpp
    /* PPF fluid extra: a modified-gravity closure defined relative to the rest of
       the universe, so it runs after the totals are assembled. A PPF fluid
       contributed 0 to the perturbation totals in the loop above (its idx_delta/
       idx_theta are unregistered); ComputePpf now fills delta_rho_fld etc. A
       non-PPF fluid is an ordinary species handled entirely in the loop above. */
    if (auto* f = ppf_fluid()) {
      f->ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
      ppw->delta_rho        += ppw->delta_rho_fld;
      ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
      ppw->delta_p          += ppw->delta_p_fld;
    }
```

(The old `rho_plus_p_tot += (1+w)·Rho()` is intentionally gone — the loop already
counted the fluid's ρ+p.)

- [ ] **Step 5: Build**

Run: `cd /Users/au192734/Projects/class_claude && pip install . --no-build-isolation 2>&1 | tail -5`
Expected: `Successfully installed ...`. If a compile error names a now-unused
local (e.g. `w_prime_fld`/`fluid_layout`/`fluid_tse` left over in the old block),
ensure the whole old block was removed in Step 4.

- [ ] **Step 6: Verify ≤0.1% across PPF / non-PPF / ΛCDM, both gauges**

Run: `python /tmp/spotcheck_308.py`
Expected: `WORST <= 1e-3` → `PASS`. The `ppf_fluid` scenario exercises the
`ComputePpf` rest-total change; `nonppf_fluid` exercises the fluid flowing
through the generic loop; `lcdm` should remain `0.00e+00`.

If any scenario exceeds 1e-3, debug with superpowers:systematic-debugging before
proceeding — likely a missed `ppw->rho_plus_p_tot` → `rho_plus_p_tot_rest`
replacement in Step 1, or a sign/skip error in the loops.

- [ ] **Step 7: Commit**

```bash
cd /Users/au192734/Projects/class_claude
git add source/perturbations_module.cpp species/fluid.cpp
git commit -m "v4 #308: unconditional stress-energy loop; PPF via ppf_fluid()

Fluid is no longer skipped in the perturbation stress-energy passes; a
PPF fluid self-zeroes through the generic interface and is handled by a
post-loop ComputePpf extra reached via ppf_fluid(). ComputePpf subtracts
its own rho+p from the now-inclusive rho_plus_p_tot. Removes the last
fluid downcasts from the stress-energy assembly.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

### Task 6: Confirm no fluid downcasts remain in scope; clang-format

**Files:**
- (verification only) all touched files

- [ ] **Step 1: Grep for residual in-scope fluid downcasts**

Run:
```bash
cd /Users/au192734/Projects/class_claude
grep -nE "static_cast<[^>]*FluidSpecies|dynamic_cast<[^>]*FluidSpecies" \
  source/background_module.cpp source/perturbations_module.cpp species/fluid.cpp
```
Expected: the **only** matches are the deferred concern-(B) sites — `background_w_fld`
(background_module.cpp ~437) and the Pk_equal reads (input_module.cpp, nonlinear_module.cpp,
which are not in this command) — plus the one localized ctor downcast in
perturbations_module.cpp (Task 2). No `FluidSpecies` downcast should remain in the
background loop, the stress-energy assembly, or the IC-validity block.

- [ ] **Step 2: clang-format the touched files (master is clang-format clean)**

Run:
```bash
cd /Users/au192734/Projects/class_claude
clang-format -i species/fluid.cpp species/fluid.h species/base_species.h \
  source/background_module.cpp source/perturbations_module.cpp source/perturbations_module.h
git diff --stat
```
Expected: only whitespace/format diffs, if any.

- [ ] **Step 3: Commit if clang-format changed anything**

```bash
cd /Users/au192734/Projects/class_claude
git add -A
git commit -m "v4 #308: clang-format touched files" || echo "nothing to format"
```

---

### Task 7: Regenerate goldens + `classyref`; full-suite gate

Output changed by design, so the exact characterization oracles are regenerated
**after** the ≤0.1% spot-checks confirmed physics is preserved.

**Files:**
- Modify: `python/background_golden/*.npz` (regenerated)
- (reinstall) `classyref`

- [ ] **Step 1: Regenerate the background goldens**

Run:
```bash
cd /Users/au192734/Projects/class_claude
python python/gen_background_golden.py
python -m pytest python/test_background_columns.py -q
```
Expected: generator writes `python/background_golden/*.npz`; the test then passes
(it now compares against the regenerated oracle).

- [ ] **Step 2: Regenerate + run the transfer-columns golden**

`gen_transfer_golden.py` has a `"fluid"` case whose values shift at ≤0.1% (the
column set/order — `d_fld`/`t_fld` — is unchanged; only values move). Regenerate
the stored `.npz` oracle, then run the test:
```bash
cd /Users/au192734/Projects/class_claude
python python/gen_transfer_golden.py
python -m pytest python/test_transfer_columns.py -q
```
Expected: generator writes `python/transfer_golden/*.npz`; test PASS. (Non-fluid
cases regenerate byte-identical.)

- [ ] **Step 3: Regenerate `classyref` (frozen reference is now stale)**

Run, in a throwaway clean worktree of this branch's HEAD (per the project recipe):
```bash
cd /Users/au192734/Projects/class_claude
git worktree add /tmp/classyref-build HEAD
cd /tmp/classyref-build
sed -i '' 's/classy-community/classy-community-ref/' pyproject.toml
pip install . --no-build-isolation \
  --config-settings=cmake.define.CLASS_PYTHON_MODULE_NAME=classyref 2>&1 | tail -5
cd /Users/au192734/Projects/class_claude
git worktree remove --force /tmp/classyref-build
```
Expected: `Successfully installed classy-community-ref-...`; `python -c "import classyref"`
works. (On macOS `sed -i ''`; on Linux use `sed -i`.)

- [ ] **Step 4: Full suite**

Run:
```bash
cd /Users/au192734/Projects/class_claude
TEST_LEVEL=2 COMPARE_OUTPUT_REF=1 python -m pytest python/test_class.py -q 2>&1 | tail -20
```
Expected: all green against the regenerated `classyref` (the COMPARE_OUTPUT_REF
gate now self-consistent). Investigate any non-fluid failure with
superpowers:systematic-debugging.

- [ ] **Step 5: Commit regenerated goldens**

```bash
cd /Users/au192734/Projects/class_claude
git add python/background_golden python/transfer_golden
git commit -m "v4 #308: regenerate background/transfer goldens after fluid refactor

Output shifts at <=0.1% by design (fluid summation reordered, ComputePpf
rest-total). Goldens regenerated after physics confirmed preserved.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

- [ ] **Step 6: Final completeness check**

Per superpowers:verification-before-completion, confirm before declaring done:
- `git log --oneline` shows Tasks 1–7 committed.
- `python /tmp/spotcheck_308.py` → `WORST <= 1e-3`.
- `python -m pytest python/test_background_columns.py python/test_transfer_columns.py -q` → PASS.
- The Task 6 grep shows no in-scope fluid downcasts remain.
