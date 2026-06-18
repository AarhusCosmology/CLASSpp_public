# Issue #308 (slice): remove the fluid's module-level downcasts via `ppf_fluid()`

Date: 2026-06-17
Issue: #308 — "remove dynamic_cast species-dispatch sites in modules → capability
queries" (v4-prep arc). This spec covers **one slice**: the residual **fluid**
casts. Family-detection (`HasNcdm`) and the IDM/DRMD/DCDM/SCF named-composite
downcasts are explicitly *out of scope* (see §6).

## 1. Background & principle

Part of the v4-prep design-review arc. The governing rule
([[feedback_no_species_picking_in_modules]]):

> Module code must loop + dispatch, never downcast. No typed accessors beyond
> photons/baryons.

Earlier #308 slices (#317, #318/#309) already converted the *dispatch* casts —
`BaseSpecies` now carries a block of neutral-default capability queries
(`DarkRadiationRhoToday`, `NeutrinoOmega0`, `NeffContribution`, …). What remains
falls in three categories:

1. **Family-detection** — `HasNcdm()` is still a `dynamic_cast` loop (+ siblings).
2. **Named-composite downcasts** — ~50 `static_cast<IDM_DR_IDR_Species&>(…)`-style
   sites across 7 modules. The bulk.
3. **Fluid PPF special-case** — the block #310 explicitly tagged "#308 territory".

This slice does **category 3 only**.

## 2. Site-by-site investigation (the seven fluid sites)

| # | Location | Reaches for | Verdict |
|---|----------|-------------|---------|
| 1 | `background_module.cpp:365` | `WriteWFld(w,dw,pvecback)` | **Dissolves** — see §3.A |
| 2 | `background_module.cpp:437` | `ComputeWFld(a,…)` (the `background_w_fld` service) | **Deferred (B)** — §6 |
| 3 | `input_module.cpp:2440` | `fluid_eos()`, `wa_fld()` (Pk_equal gate) | **Deferred (B)** — §6 |
| 4 | `nonlinear_module.cpp:3473` | `w0_fld()` (Pk_equal seed) | **Deferred (B)** — §6 |
| 5 | `perturbations_module.cpp:515,526` | `use_ppf()` (IC validity tests) | **In scope** — §3.C |
| 6 | `perturbations_module.cpp:4892–4898` | `use_ppf()` (shear-pass skip) | **In scope** — §3.C |
| 7 | `perturbations_module.cpp:4967–5002` | `use_ppf()`, `Rho()`, `cs2_fld()`, `ComputePpf(…)` | **In scope** — §3.C |

Key facts established during investigation (do not re-derive):

- **`FluidSpecies::ComputeBackground` (fluid.cpp:63–67) ignores `a_rel`** and only
  copies the ODE-integrated `rho_fld` out of `pvecback_B`. The `w`/`dw` slots are
  filled *separately* by the module (`background_w_fld` → `ComputeWFld` →
  `WriteWFld`) *before* the call. The comment's "checked error handling"
  justification is **vestigial**: `ComputeWFld`'s only failure paths are
  `class_stop` (EDE-unfinished :398, EDE+DCDM :350), which *throw* post-#38/#311.
  There is no `_FAILURE_` return to propagate.
- **`FluidSpecies` already implements the generic `Delta`/`Theta`/`DeltaP`/
  `RhoPlusPShear`** (fluid.cpp:256–305). `DeltaP` computes exactly the non-PPF
  `cs2·δρ + (cs2−ca2)·3ℋ(ρ+p)θ/k²` the module's site-7 block recomputes inline.
- **`ComputePpf` (fluid.cpp:412–) reads the accumulated global totals**
  `ppw->delta_rho`, `rho_plus_p_theta`, `rho_plus_p_shear`, `rho_plus_p_tot` —
  the non-fld stress-energy. This is the genuine modified-gravity coupling: PPF is
  defined relative to the whole universe, so the PPF fluid *must* run last and is
  *legitimately* module-orchestrated. (Same finding as #310 re the `*_fld`
  workspace fields.)
- **`use_ppf_` is construction-fixed** (fluid.cpp:26 ctor init only; never
  reassigned). **At most one fluid** can exist (single factory key `"Fluid"`,
  all_species.h:53). ⇒ at most one PPF fluid.
- `FluidSpecies::FillSources` (fluid.cpp:179–198) and `PrintVariables`
  (fluid.cpp:238–253) consume `ppw->{delta_rho,rho_plus_p_theta,delta_p}_fld`. In
  the final design these workspace fields become **PPF-only** (written by
  `ComputePpf`); for the non-PPF fluid both methods compute their quantities
  locally from `y`/`pvecback` (e.g. `delta_rho_fld = Rho()·y[idx_delta]`), which
  retires the vestigial general handoff #310 flagged.
- **The fluid's generic perturbation methods already return 0 for PPF** — a PPF
  fluid registers `idx_Gamma`, never `idx_delta`/`idx_theta`
  (fluid.cpp:112–126), so `Delta`/`Theta`/`DeltaP` (guarded on `idx_delta>=0`)
  and `RhoPlusPShear` (always 0) all yield 0. **No fluid-method edits needed**
  for the PPF fluid to self-zero in the generic loop.

## 3. Design

### A. Background dissolution (site 1) — fluid rejoins the main loop

Make the fluid an ordinary species in `background_functions()`:

- `FluidSpecies::ComputeBackground(a_rel, pvecback_B, pvecback)` calls its own
  `ComputeWFld(a_rel, …)` and writes `index_bg_w_fld_` /
  `index_bg_dw_over_da_fld_` itself (folding in `WriteWFld`), then copies
  `rho_fld` as today. **Assume `a == a_rel`** (i.e. `a_today == 1`): `a_today`
  defaults to 1 (background.h:103, "arbitrary and irrelevant") and `a/a_rel` is
  already used inconsistently across the codebase (longer-term direction: remove
  `a_today` entirely — out of scope here).
- **Delete the fluid special-case** (background_module.cpp:355 `if (name ==
  "Fluid") continue;` and the entire post-loop block 361–369). Fluid is computed
  and accumulated inside the single `for (sp : all_species_)` loop like every
  other species. This reorders the `rho_tot += rho` sum (background_module.cpp:333)
  → ULP shift in `rho_tot`/`H`; **acceptable** (clean over bit-identical — §5).
- `WriteWFld` becomes unused → remove it (decl fluid.h, def fluid.cpp).
- `RequiresDeferredBackground()` is **dead code** — declared on `BaseSpecies`
  (base_species.h:268) and overridden by `FluidSpecies` (fluid.h:86) but **read
  nowhere** (the old loop skipped by name, not by predicate). Remove both the
  override and the base virtual.
- `background_w_fld` (the service, §6) is *unchanged*; its caller at
  background_module.cpp:364 is removed (folded into ComputeBackground).

### B. `ppf_fluid()` — one privileged, localized accessor in `PerturbationsModule`

- Add a private `FluidSpecies* ppf_fluid_ = nullptr;` to `PerturbationsModule`,
  resolved **once** in the constructor (the single, localized downcast):
  ```cpp
  if (auto* p = all_species_.find("Fluid")) {
    auto* f = static_cast<FluidSpecies*>(p->get());
    if (f->use_ppf()) ppf_fluid_ = f;
  }
  ```
- Expose `FluidSpecies* ppf_fluid() const { return ppf_fluid_; }` (nullable).
- Rationale for module-local (not on `SpeciesCollection`): the PPF coupling is
  perturbations-only; keeping it here avoids giving the generic collection a
  concrete `fluid.h` dependency (`photons_`/`baryons_` are stored as
  `BaseSpecies*`). `if (auto* f = ppf_fluid())` encodes both invariants — "the
  one PPF fluid is present" and "`use_ppf()` is true" — so `use_ppf()` is never
  called at a module site again.

### C. The stress-energy loop becomes unconditional in the fluid

The fluid stops being skipped in the accumulating loops; only the (genuinely
special) Photons/Baryons keep dedicated paths. The PPF fluid self-zeroes through
the generic interface (see §2), so the only fluid-specific handling left is a
small post-loop "PPF extra" reached through `ppf_fluid()`.

- **Site 5** (perturbations IC validity, ~514–535): inside the existing
  `if (all_species_.count("Fluid"))`, replace the
  `static_cast<const FluidSpecies&> fluid` + `!fluid.use_ppf()` with `!ppf_fluid()`
  (single-fluid ⇒ `count("Fluid") && !ppf_fluid()` ⇔ a non-PPF fluid is present).
  The `background_w_fld` calls stay (module service).
- **Site 6** (shear pass, 4896–4901): skip only `Photons` / `Baryons` (drop the
  `fluid_tse`/`use_ppf()` downcast entirely). Both PPF (always 0) and non-PPF
  (0) fluid `RhoPlusPShear` add nothing.
- **Site 7** (pass 2, 4914–4952): skip only `Photons` / `Baryons` (drop the
  `"Fluid"` skip). The fluid contributes through the generic `Delta`/`Theta`/
  `DeltaP` + `rho_plus_p_tot += rho + P()`: a **non-PPF** fluid contributes its
  real values; a **PPF** fluid contributes 0 to `delta_rho`/`rho_plus_p_theta`/
  `delta_p` (self-zeroing) but its real `rho + p` to `rho_plus_p_tot`.
- **PPF extra** (replaces the old dedicated block 4967–5002):
  ```cpp
  if (auto* f = ppf_fluid()) {
    f->ComputePpf(k, a, a_prime_over_a, ppr, y, ppw);
    ppw->delta_rho        += ppw->delta_rho_fld;
    ppw->rho_plus_p_theta += ppw->rho_plus_p_theta_fld;
    ppw->delta_p          += ppw->delta_p_fld;
  }
  ```
  The old `rho_plus_p_tot += (1+w)·rho` line is **removed** (the loop already
  counted it).
- **`ComputePpf` self-subtracts its own background ρ+p.** `ComputePpf` needs the
  *rest-of-universe* ρ+p total (fluid.cpp:458,473 — "the whole point of the PPF
  scheme"), which `ppw->rho_plus_p_tot` no longer excludes. Inside `ComputePpf`,
  compute once and use at the 3 sites:
  ```cpp
  const double rho_plus_p_tot_rest =
      ppw->rho_plus_p_tot - Rho(ppw->pvecback) * (1. + w_fld);
  ```
  `delta_rho`/`rho_plus_p_theta`/`rho_plus_p_shear` need no adjustment — the PPF
  fluid contributed 0 to those.

## 4. Components & interfaces touched

- `species/fluid.h` / `species/fluid.cpp`: `ComputeBackground` self-computes w;
  remove `WriteWFld`; remove the dead `RequiresDeferredBackground` override;
  `ComputePpf` self-subtracts its own ρ+p from `rho_plus_p_tot`;
  `FillSources`/`PrintVariables` compute non-PPF fluid quantities locally.
- `species/base_species.h`: remove the dead `RequiresDeferredBackground` virtual.
- `source/background_module.cpp`: delete the fluid special-case in
  `background_functions()`; fluid rejoins the main loop.
- `source/perturbations_module.{h,cpp}`: add `ppf_fluid_` + accessor + ctor
  resolution; make the stress-energy loops unconditional in the fluid; replace
  the dedicated block with the PPF-extra; rewrite site 5.

No `BaseSpecies` virtuals are added (the PPF surface is fluid-specific and stays
behind the module-local accessor; we do not pollute the base class). One dead
`BaseSpecies` virtual (`RequiresDeferredBackground`) is *removed*.

## 5. Verification

**Clean code over bit-identical** (per [[feedback_no_bit_identical_requirement]]):
this slice intentionally changes floating-point summation order in two places
(fluid folded into the background loop; fluid folded into the perturbation
stress-energy loop; `ComputePpf` subtraction), so output will differ at the ULP
level and beyond (the ODE amplifies). The bar is **≤ 0.1%**, verified with `Cl^TE`
zero-crossing handling, never blind max-rel-diff.

Because output changes intentionally, the characterization goldens and the
`classyref` reference are **regenerated** after confirming physics is preserved
at tolerance (the workflow from [[feedback_vectorization_reduction_drift]]: verify
≤0.1%, *then* regenerate the baseline):
- `python/background_golden/*.npz` (via `gen_background_golden.py`) after Task A.
- `classyref` (rebuild current branch under module name `classyref`) before the
  final full-suite run, since the newtonian/sync output shifts.

Gauntlet (in order):
1. **During development** — `classy` (this branch) vs the **frozen pre-slice
   `classyref`** spot-checks across **both gauges** (lensed TT + P(k)), for: ΛCDM,
   **PPF CPL fluid** (`use_ppf=yes`, `wa_fld≠0`), **non-PPF CPL fluid**
   (`use_ppf=no`). Expect **≤ 0.1%** (not 0 — this measures the intended drift).
2. **After ≤0.1% confirmed** — regenerate `background_golden` + `classyref`.
3. **Final gate** — `python/test_background_columns.py`,
   `python/test_transfer_columns.py`, and full `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1`
   against the regenerated `classyref`.

## 6. Explicitly out of scope (tracked follow-ups)

- **Concern (B) — fluid EoS reads** (sites 2/3/4): `background_w_fld`'s internal
  cast, and the Pk_equal `w0/wa/cs2/fluid_eos` reads in input/nonlinear. These
  apply to *any* fluid (PPF or not) and bake in single-fluid assumptions →
  resolve alongside the eventual **multi-fluid** work, when the right "w(a) with N
  fluids" shape is knowable. YAGNI now.
- **Category 1 — family-detection** (`HasNcdm` + NCDM/DNCDM `dynamic_cast`
  loops): the immediate next #308 slice.
- **Category 2 — IDM/DRMD/DCDM/SCF named-composite downcasts**: later slice(s);
  some are genuine module-level physics coupling needing per-case judgement.

## 7. Multi-fluid forward-compatibility

This design is chosen to age well: when multiple fluids land, ordinary (non-PPF)
fluids are *already* generic-loop species (they implement the full generic
interface), and only the single PPF fluid stays privileged. `ppf_fluid()` remains
well-defined precisely because PPF is intrinsically singular; a generic `fluid()`
accessor would have become ambiguous.
