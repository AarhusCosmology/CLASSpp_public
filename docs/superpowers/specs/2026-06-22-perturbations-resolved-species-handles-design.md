# Perturbations hot-path: resolve frozen species lookups once (design)

Date: 2026-06-22
Status: design — approved scope, not yet implemented
Branch (proposed): `perturb-resolved-species-handles`
Related: PR #327 (active-species dispatch cache), #328 / #329 (StressEnergy
delegation), memory `project-perturbations-regression-vs-c`.

## 1. Background & motivation

The perturbations performance programme has progressively removed per-ODE-step
overhead from the right-hand side. After #328/#329 the generic stress-energy and
derivative paths are fully delegated (`StressEnergy()` / `PerturbDerivs()` virtuals
looped over `pv->active_species`, no name-picking, no downcasts).

What remains in the hot path is a scattering of **string-keyed lookups into the
`SpeciesCollection`** — `all_species_.count("UR")`, `.at("IDM_DR_IDR")`,
`.find("UR")`, `.index_of("CDM")`, and the free function `HasNcdm(all_species_)`.
Each of these is an **O(n) linear scan doing `std::string` comparisons**
(`species/species_collection.cpp:53-81`), executed per ODE step / per source
sample / per approximation-bisection probe. Crucially, **the collection is frozen
after construction** (`SpeciesCollection::freeze()`), so every one of these lookups
resolves to the *same answer for the module's entire lifetime*.

A prototype that cached these results measured roughly **+3% on the Perturbations
stage, bit-identical** (memory `project-perturbations-regression-vs-c`,
experiments #3/#4). The single largest chunk (`+0.85%`, mean/se 3.48) came from
`perturb_einstein`'s three `count("IDM_DR_IDR")` guards — and `IDM_DR_IDR` is
**absent** in the standard benchmark, so those are pure no-op scans of all species
on every step.

This work was deliberately scoped out of #329 as the follow-up. The prototype
cached results as scattered loose `bool`/pointer members ("fine for the
experiment"). This spec defines the **clean final design**: resolve each lookup
once, behind a principled choice of mechanism, with no per-step scan and no
per-step downcast.

## 2. Goal & non-goals

**Goal.** Eliminate every *repeated* (per-step / per-sample / per-k-bisection)
string-keyed `SpeciesCollection` lookup from the perturbations hot path, capturing
the ~1–3% Perturbations-stage win, while *improving* design cleanliness rather than
bolting on caches.

**Non-goals.**
- Touching the **one-shot-per-run / per-k setup** lookups (see §6) — not hot,
  explicitly out of scope.
- Adding a per-species `UpdateApproximationFlags` hook / making species own their
  approximation schemes. That is a larger, separate (#280-style) programme; this PR
  keeps the approximation *framework* in the module.
- Changing physics. Output must stay within the 0.1% bar and is expected to be
  **byte-identical** (a cached value equals the scan result for a frozen
  collection; no floating-point reassociation is introduced).

## 3. Guiding principle (the litmus)

Each residual lookup is replaced by one of two mechanisms, chosen by **ownership**:

> **Resolved handle** when the surrounding logic *belongs to the module* — an
> Einstein/metric equation, an approximation scheme, an output convention.
> **Delegation** (a per-species virtual summed in a loop) when it computes a
> *species' own contribution to a sum*.

Applying the litmus honestly: because #278/#328/#329 already delegated the
genuinely species-owned physics, **the residual lookups sit overwhelmingly inside
module-owned machinery** (Einstein equations; the TCA / RSA / UFA / ncdmfa
approximation framework, whose flags `ppw->approx[]`, indices `index_ap_*`, trigger
thresholds `ppr->*`, and switch-point bisection are *all* module state; metric
sources; the δ_tot output convention). The principled fix for these is therefore a
**resolved, typed handle**, not delegation — and not a bare `bool`, because a typed
handle also eliminates the *per-step downcast* that currently accompanies each
lookup.

**Refinement (approximation state).** The module *owns and decides* the
approximation state and **passes** it where physics needs it; the module must not
*compute a species' physical quantity by branching on a flag* (`if (approx)
delta_rho = A else delta_rho = B`) — that physics belongs in the species. The
handle fix in this spec **does not introduce** that anti-pattern: the IDR
shear/RSA additions in `perturb_einstein` mirror the *already-inline photon*
treatment (the photon TCA shear correction directly above them uses the
`photons()` handle the same way), and any genuinely flag-conditional species
physics already lives in the species (e.g. baryon δp under TCA in
`BaryonsSpecies::StressEnergy`). The single genuinely species-owned sum in the
residual (tensor relativistic ρ) is *already* delegated for ncdm; we only tidy its
stray `find("UR")`.

## 4. Site inventory & per-site fix

Line numbers are as of master `eba6f7d7` and **will shift**; functions are the
stable reference.

### Per-RHS-step (`perturb_derivs_member` → `perturb_einstein` → {stress-energy, rsa})

**Site 1 — `perturb_einstein`, IDR RSA/TCA shear corrections.**
`count("IDM_DR_IDR")` ×3 guards (≈4639, 4660, 4699), `at("IDM_DR_IDR")` ×3
(≈4665, 4706, 4712), `index_of("IDM_DR_IDR")` (≈4701).
*Owner: module* (Einstein eqs + TCA-idm_dr / RSA-idr scheme; parallels the inline
photon shear correction).
*Fix:* typed handle `resolved_.idm_dr_idr` (`const IDM_DR_IDR_Species*`). Guards
become `if (resolved_.idm_dr_idr && approx_on)`; bodies use the handle directly
(no `at`, no per-step `static_cast`). The `index_of` for `species_layouts[i]` is
replaced by a resolved index (see Site 1b). **This is the +0.85% win**: absent →
null-check, zero work.

**Site 1b — `perturb_einstein`, IDR layout index.** `index_of("IDM_DR_IDR")`
(≈4701) feeds `ppw->pv->species_layouts[idr_th_i]`. Resolve once as
`resolved_.idm_dr_idr_index` (`std::size_t`), or read it from the cached species'
`collection_index_` (every species is stamped with its sorted index in
`freeze()`, `species_collection.cpp:30-31`) — prefer the latter to avoid a second
resolved field: `resolved_.idm_dr_idr->collection_index_`.

**Site 2 — `perturb_rsa_delta_and_theta`, UR.** `count("UR")` ×3 (≈6331, 6372,
6393), `at("UR")->Rho()` ×3 (≈6394-6396).
*Owner: module* (RSA scheme; `rsa_delta_ur`/`rsa_theta_ur` are module workspace
fields the UR species reads back; the photon path is inline-with-`photons()`).
*Fix:* base handle `resolved_.ur` (`const BaseSpecies*`). `if (count("UR"))` →
`if (resolved_.ur)`; the three `at("UR")->Rho(...)` collapse to one
`resolved_.ur->Rho(...)`.

**Site 3 — `perturb_total_stress_energy`, tensor branch.** `find("UR")` ×2
(≈4934, 4939), `HasNcdm(all_species_)` (≈4942). Runs only for **tensor** modes.
*Owner: mixed* — the ncdm relativistic-ρ sum is a *species contribution* and is
**already delegated** (`TensorMasslessRelativisticRho` loop, ≈4943-4944); the
`tm_exact`/`tm_massless` selection is module-owned.
*Fix:* keep the ncdm delegation untouched. Replace the stray `find("UR")` with
`resolved_.ur` and `HasNcdm(...)` with `resolved_.has_ncdm`, inside the existing
module-owned method selection. (A deeper unification — UR opting into
`FreestreamingRho`/`TensorMasslessRelativisticRho` so the whole sum is one loop —
is **not** done here: the `tm_exact` case must exclude ncdm, so a blunt loop would
change behavior. Noted as a possible future cleanup.)

### Per-source-sample (`perturb_sources_member`, also calls `perturb_einstein`)

**Site 4a — CDM h-source (synchronous gauge).** `index_of("CDM")` (≈5238) feeds
`species_layouts[cdm_i]` for `source_h = -2·δ_cdm`.
*Owner: module* (it's a metric source).
*Fix:* a cached `cdm_layout` pointer on `perturb_vector`, wired in
`perturb_vector_init` exactly like `photon_layout`/`baryon_layout`
(`perturbations_module.cpp:2949-2950`). Site becomes
`static_cast<const CDMSpecies::PerturbLayout&>(*ppw->pv->cdm_layout)`. CDM is
guaranteed present in synchronous gauge (the existing `class_test`), so the wiring
sets `cdm_layout` whenever `resolved_.cdm_index` is valid.

**Site 4b — Lambda excluded from ρ_tot.** `count("Lambda")` (≈5297),
`at("Lambda")->Rho()` (≈5299): `ρ_tot = ρ_tot_bg − ρ_Λ` (CMBFAST δ_tot
convention).
*Owner: module* (output convention).
*Fix:* base handle `resolved_.lambda`. `if (count("Lambda"))` →
`if (resolved_.lambda)`; `at("Lambda")->Rho` → `resolved_.lambda->Rho`. Behavior
(Lambda-only exclusion) preserved exactly.

### Per-k approximation bisection (`perturb_approximations` via `perturb_find_approximation_switches`)

**Site 5 — approximation triggers.** `count("IDM_DR_IDR")` ×2 (≈4267, 4334) +
`at` (≈4269, 4335), `count("IDM_DRMD_IDR_DRMD")` (≈4296) + `at` (≈4303, 4307),
`count("UR")` (≈4364), `HasNcdm(...)` (≈4374).
*Owner: module* — the whole approximation framework (flags, `index_ap_*`,
`ppr->*_trigger_*`, switch bisection) is module state; species supply only data
values (`nindex_idm_dr`, `idr_nature`) and presence.
*Fix:* reuse the same resolved handles — `resolved_.idm_dr_idr`,
`resolved_.idm_drmd` (`const IDM_DRMD_IDR_DRMD_Species*`), `resolved_.ur`,
`resolved_.has_ncdm`. Presence guards become handle null-checks; the typed handles
also remove the per-probe `static_cast`s. No new hook; trigger logic stays in the
module.

## 5. Mechanism

### 5.1 One resolved-handle struct on the module

A private struct, one member, populated once by `ResolveSpecies()` called from the
constructor (after `all_species_` is frozen). It **consolidates** the two existing
scattered caches (`idr_nature_`, `ppf_fluid_`) into a single documented home.

```cpp
// perturbations_module.h, private section.
//
// Resolved-once views into the frozen all_species_ collection. The perturbations
// module owns several numerical approximation schemes (TCA/RSA/UFA/ncdmfa), the
// Einstein/metric equations, and output conventions that legitimately reference
// specific named species. These handles give that module-owned logic O(1) access
// to the relevant species' data, replacing per-step/per-sample/per-k O(n)
// string-keyed scans of all_species_ (frozen after construction) and the
// per-step downcasts that accompanied them. Resolve-by-handle vs delegate-by-hook
// follows the ownership litmus (see design doc): a handle when the surrounding
// logic belongs to the module; delegation when summing a species' own
// contribution (already done for stress-energy / tensor ncdm).
struct ResolvedSpecies {
  // Typed handles — module-owned schemes read these species' data directly,
  // hoisting the existing per-step static_casts to construction time.
  const IDM_DR_IDR_Species*        idm_dr_idr = nullptr;  // TCA-idm_dr, RSA-idr
  const IDM_DRMD_IDR_DRMD_Species* idm_drmd   = nullptr;  // TCA-idm_drmd
  FluidSpecies*                    ppf_fluid  = nullptr;  // PPF fluid (absorbs ppf_fluid_)

  // Base handles — only presence + Rho() are needed.
  const BaseSpecies* ur     = nullptr;  // UFA, RSA-ur, tensor relativistic rho
  const BaseSpecies* lambda = nullptr;  // excluded from rho_tot (delta_tot convention)

  // Presence / scalars.
  bool        has_ncdm   = false;                      // ncdmfa scheme present
  int         idr_nature = idr_free_streaming;          // (absorbs idr_nature_)
  std::size_t cdm_index  = static_cast<std::size_t>(-1);// h source (sync gauge); ResolveSpecies
                                                        // sets it (== all_species_.size() if absent)
};
ResolvedSpecies resolved_;

void ResolveSpecies();  // called once from the ctor
```

`ResolveSpecies()` (replaces the ad-hoc ctor blocks at `perturbations_module.cpp:72-83`):

```cpp
void PerturbationsModule::ResolveSpecies() {
  if (auto* p = all_species_.find("IDM_DR_IDR"))
    resolved_.idm_dr_idr = static_cast<const IDM_DR_IDR_Species*>(p->get());
  if (resolved_.idm_dr_idr)
    resolved_.idr_nature = resolved_.idm_dr_idr->idr().idr_nature();

  if (auto* p = all_species_.find("IDM_DRMD_IDR_DRMD"))
    resolved_.idm_drmd = static_cast<const IDM_DRMD_IDR_DRMD_Species*>(p->get());

  if (auto* p = all_species_.find("UR"))     resolved_.ur     = p->get();
  if (auto* p = all_species_.find("Lambda")) resolved_.lambda = p->get();

  resolved_.has_ncdm  = HasNcdm(all_species_);
  resolved_.cdm_index = all_species_.index_of("CDM");  // size() if absent

  if (auto* p = all_species_.find("Fluid")) {
    auto* f = static_cast<FluidSpecies*>(p->get());
    if (f->use_ppf()) resolved_.ppf_fluid = f;
  }
}
```

The `ppf_fluid()` accessor returns `resolved_.ppf_fluid`; usages of `idr_nature_`
are repointed to `resolved_.idr_nature`. The standalone members `idr_nature_` and
`ppf_fluid_` are removed.

### 5.2 `perturb_vector::cdm_layout`

```cpp
// source/perturbations.h, next to photon_layout / baryon_layout
BaseSpecies::PerturbLayout* cdm_layout = nullptr;
```

Wired in `perturb_vector_init` (mirrors lines 2949-2950):

```cpp
if (resolved_.cdm_index < all_species_.size())
  ppv->cdm_layout = ppv->species_layouts[resolved_.cdm_index].get();
```

### 5.3 Why not delegation here, and why not bare bools

- **Not delegation:** by the litmus these are module-owned schemes/conventions;
  delegating them would either (a) require new hooks the user explicitly excluded
  from this PR, or (b) for the loop-style sites, replace a cheap cached null-check
  with N usually-no-op virtual calls per step — a *regression* in the common
  (species-absent) case that the +0.85% win depends on. Delegation stays where it
  already is (generic stress-energy/derivs, tensor ncdm).
- **Not bare bools:** typed handles (`idm_dr_idr`, `idm_drmd`) remove the per-step
  `static_cast` downcasts as well as the scans, and follow the existing
  `ppf_fluid_` (typed `FluidSpecies*`) precedent. This does **not** add new
  species-type picking to the generic loops, and adds **no** typed accessors to
  `SpeciesCollection` (respecting `feedback_no_species_picking_in_modules`): it
  hoists casts that already exist at these module-owned sites to construction time.

## 6. Explicitly out of scope (one-shot lookups — not touched)

All lookups that run once per run or once per k at setup, not in any repeated hot
loop: `perturb_init`, `perturb_indices_of_perturbs`, `perturb_get_k_list`,
`perturb_workspace_init` (including the `class_define_index` approximation-index
*registration*), `perturb_solve`, the non-bisection parts of
`perturb_find_approximation_switches`, `perturb_vector_init`'s own
switch-time/index scans, `perturb_initial_conditions`, every `class_test` gauge
guard, and `perturb_print_variables_member`. The generic `StressEnergy()` /
`PerturbDerivs()` loops are unchanged.

## 7. Correctness & verification

**Physics bar (gating):** `.superpowers/sdd/tol_check.py` at ~0.1% on Cl^TT and
P(k) (no TE/EE/BB zero-crossing max-rel-diff), per
`feedback_no_bit_identical_requirement`. Run the benchmark
`base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini` plus a few composites covering the
touched species (a config with `IDM_DR_IDR`, one with `Fluid`/Lambda, one tensor
run exercising Site 3, a synchronous-gauge run exercising Site 4a).

**Byte-identity (free bonus, expected):** because each resolved value equals the
scan result for a frozen collection and no arithmetic is reassociated, `cl_*.dat`
and `pk*.dat` are expected byte-identical to master on every scenario. Verify and
report it; if any scenario is *not* byte-identical, stop and investigate (it would
indicate a logic error in resolution, not benign FP drift).

**Full suite:** `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` against master classyref.

## 8. Benchmark methodology (mandatory rigor)

This is a ~1–3% effect near the laptop noise floor; the #329 win was nearly lost to
a contaminated baseline. Therefore:

- Build **both** master and branch **fresh from isolated git worktrees** with
  **identical flags** (`-O3 -DNDEBUG`, same `/usr/bin/c++`) into **separate**
  build dirs.
- **Order-alternating, paired A/B** on
  `base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini`, `OMP_NUM_THREADS=1`, the
  per-stage Perturbations timing from `./class_profiled`.
- Many rounds (≥40+40); report **min and median** plus a **bootstrap CI** and
  paired mean/se on the per-round deltas.
- Since the number is small, confirm **algorithmic vs layout** with
  `-falign-functions=64` variants of both binaries (the gain should survive a
  perturbed layout).

## 9. Risks & mitigations

| Risk | Mitigation |
|---|---|
| A resolved handle is `nullptr` where code assumed presence | The only mandatory-presence site is CDM in synchronous gauge, already guarded by an existing `class_test`; `cdm_layout` is wired only when present and used only on that gauge path. |
| `IDM_DRMD_IDR_DRMD` typed handle type/const mismatch with existing `static_cast<...&>` sites | Match the existing cast's const-ness; the header is already included. |
| Win is layout noise, not algorithmic | The `-falign-functions=64` probe in §8. |
| Output drift from accidental logic change | Byte-identity check (§7) catches any non-trivial change immediately. |
| Over-reach into one-shot setup sites | §6 lists the exclusions explicitly; reviewer checks the diff touches only the §4 functions + the struct/wiring. |

## 10. Suggested implementation order (single PR)

1. Add `ResolvedSpecies` + `ResolveSpecies()`; fold in `idr_nature_`/`ppf_fluid_`;
   call from ctor. (No call-site changes yet — compiles, behavior identical.)
2. Repoint `ppf_fluid()` and `idr_nature_` usages to `resolved_`; delete the old
   members. Verify byte-identical.
3. Site 1 + 1b (`perturb_einstein`) — the perf-critical change. Benchmark here to
   confirm the win lands before proceeding.
4. Sites 2, 3 (`perturb_rsa_delta_and_theta`, tensor branch).
5. Site 4a (`cdm_layout` on `perturb_vector` + wiring) and 4b (Lambda).
6. Site 5 (`perturb_approximations`).
7. clang-format (22.1.3) every touched file. Full TEST_LEVEL=2; final paired A/B
   benchmark + layout probe.

`clang-format` every touched file. **Never `git add -A`** in this repo (build
artifacts leak in).
