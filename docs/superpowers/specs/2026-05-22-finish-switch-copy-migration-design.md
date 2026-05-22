# Finish Switch-Copy Migration (N_decay_dr-free species)

## Background

This is **PR A** of a four-PR follow-up sequence to PR #268 ("Encapsulate species perturbation
layouts; remove `ncdm_id`"). The agreed order is **A → C → B → D**:

- **A** (this spec): finish moving inline switch-copy logic into species `CopyPerturbationsAcrossSwitch` overrides, for the species not entangled with `N_decay_dr`.
- **C**: execute the (stale, needs-refresh) has-guards cleanup plan.
- **B**: remove `pba->N_decay_dr`; also finishes the DR-bearing composites' switch-copy migration deferred from A.
- **D**: shooter hooks (encapsulate-species-layouts plan Tasks 37–40).

PR #268 introduced `BaseSpecies::CopyPerturbationsAcrossSwitch(old_layout, new_layout, old_y, new_y, ctx)`
and a per-species dispatch loop at the top of `perturb_vector_init`'s "switching approximation"
(`else`) branch:

```cpp
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->CopyPerturbationsAcrossSwitch(*ppw->pv->species_layouts[i],
                                                 *ppv->species_layouts[i],
                                                 ppw->pv->y, ppv->y, switch_ctx);
}
```

#268 migrated the **approximation-driven** copies (UR ufa, IDR rsa/tca, IDM_DRMD tca, Photons,
NCDM fa-collapse) into top-level species overrides that fire through this loop. As a result the
old approximation-switch sections for those species are now empty except for verbose `fprintf`s.

### What still lives inline (the subject of this PR)

The "always-reconducted" copy block (the section commented *"some variables (b, cdm, fld, ...)
are not affected by any approximation. They need to be reconducted whatever the approximation
switching is"*) still does manual, layout-based slot copies for:

- **Baryons** — `idx_delta`, `idx_theta` (unconditional).
- **CDM** — `idx_delta` always; `idx_theta` only when `ppt->gauge == newtonian`.
- **Fluid** — non-PPF: `idx_delta`, `idx_theta`; PPF: `idx_Gamma`.
- **ScalarField** — `idx_phi`, `idx_phi_prime`.
- **IDM_DR** (child of the `IDM_DR_IDR` composite) — `idx_delta`, `idx_theta`.
- **DCDM / DR** (children of `DCDM_DR`) — DCDM `idx_delta`/`idx_theta`, plus the DR multipole
  copy that loops over `pba->N_decay_dr` channels. **Deferred to PR B.**

Additionally, the child species `IDM_DRSpecies`, `IDRSpecies`, `IDM_DRMDSpecies`,
`IDR_DRMDSpecies` already *have* `CopyPerturbationsAcrossSwitch` overrides, but their
composites (`IDM_DR_IDR`, `IDM_DRMD_IDR_DRMD`) do **not** override/delegate, so those child
methods are currently **dead code**.

### Latent #268 regression (idr / idr_drmd / idm_drmd reconduct dropped)

Investigation (git archaeology + scenario coverage) found that wiring those composites is **not**
a neutral relocation. Pre-#268, `idr` delta/theta/shear/l3 were reconducted across the general
approximation transitions (old `perturb_vector_init` lines ~3840/3922/4016/4162); `idr_drmd` and
`idm_drmd` likewise. #268 consolidated that logic into the child `CopyPerturbationsAcrossSwitch`
overrides but never wired the composites to call them. The current `master` therefore **drops**
those always-reconduct copies: `idr` gets only a tca_idm_dr-off shear/l3 *recompute*, and
`idm_drmd`/`idr_drmd` get copies *only* inside the tca_idm_drmd-off block. On a general switch
(e.g. photon tca-off) while those species are integrated, their slots are not reconducted.

This is invisible to the five regression scenarios because the only idm_dr scenario
(`ncdm_dncdm_idmdr_combined.ini`) is **background-only** (no `output=`), so it never integrates
idm_dr/idm_drmd perturbations through a switch. Wiring the composite delegation therefore
**restores the dropped copies** — a genuine behavior change (a regression fix) that the current
suite cannot verify. PR A therefore adds full-Cl idm_dr and idm_drmd scenarios so the restored
behavior is exercised and checked against the pre-#268 reference (which had the copies).

## Goals

1. Add `CopyPerturbationsAcrossSwitch` overrides to **Baryons, CDM, Fluid, ScalarField**, each
   reproducing exactly what its inline block does today, guarded by `idx >= 0` so gauge/mode/PPF
   variation falls out of the layout (the idiom the existing `UltraRelativisticSpecies` override
   uses).
2. Add full-Cl **idm_dr** and **idm_drmd** regression scenarios (`output=tCl,pCl,lCl,mPk`) that
   integrate those species' perturbations through approximation switches — coverage that does not
   exist today — so the composite-delegation behavior can be verified.
3. Wire **`IDM_DR_IDR_Species`** and **`IDM_DRMD_IDR_DRMD_Species`** to override
   `CopyPerturbationsAcrossSwitch` and delegate to their children's existing overrides via the
   nested layout members (`my.idm_dr`, `my.idr`; `my.idm_drmd`, `my.idr_drmd`). This **restores
   the idr/idr_drmd/idm_drmd reconduct dropped in #268** (a behavior change, verified by goal 2
   against the pre-#268 reference).
4. Delete the corresponding inline always-reconduct blocks from `perturb_vector_init` once the
   species/composite overrides cover them.

## Non-goals

- **DR-bearing composites** (`DCDM_DR`, `DNCDM_DR`) and the DR multipole / `N_decay_dr` copy →
  **PR B**. DR genuinely needs the decay-channel count to copy its per-channel sub-hierarchies;
  PR B reworks that storage while removing `N_decay_dr`. (`DNCDM_DR` perturbation-switch copying
  also can't be exercised by the background-only `dncdm_dr.ini` scenario, another reason to handle
  it where the storage is being reworked.)
- **Approximation-IC recompute** — tca-off → photon shear/l3/pol IC, tca_idm_dr-off → IDR
  shear/l3 (`TcaShearIdr`), tca_idm_drmd-off → DRMD. These *compute* new initial conditions from
  tight-coupling formulas; they are not old→new copies and do not fit the
  `CopyPerturbationsAcrossSwitch` shape. Out of scope (would need a separate hook).
- **Metric and perturbed-recombination** copies (`index_pt_eta`/`phi`/`gw`/`gwdot`/`hv_prime`/`V`,
  recomb temp/chi) — not species-owned; stay module-level.

## Architecture

Each species override has the shape established by `UltraRelativisticSpecies`:

```cpp
void XxxSpecies::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                               const BaseSpecies::PerturbLayout& new_base,
                                               const double* old_y, double* new_y,
                                               const PerturbSwitchContext& /*ctx*/) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  if (old_l.idx_delta < 0 || new_l.idx_delta < 0) return;   // not allocated in this mode
  new_y[new_l.idx_delta] = old_y[old_l.idx_delta];
  if (new_l.idx_theta >= 0 && old_l.idx_theta >= 0)
    new_y[new_l.idx_theta] = old_y[old_l.idx_theta];
  // ... remaining slots, each guarded by idx >= 0
}
```

Per-species slot sets:

| Species | Slots copied | Conditional |
|---|---|---|
| Baryons | `idx_delta`, `idx_theta` | both always present |
| CDM | `idx_delta`, `idx_theta` | `idx_theta` present only in newtonian gauge (`idx_theta >= 0`) |
| Fluid | `idx_delta`+`idx_theta` (non-PPF) **or** `idx_Gamma` (PPF) | copy whichever are `>= 0` |
| ScalarField | `idx_phi`, `idx_phi_prime` | both always present |

Composite delegation:

```cpp
void IDM_DR_IDR_Species::CopyPerturbationsAcrossSwitch(/* bases, y, ctx */) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  idm_dr_->CopyPerturbationsAcrossSwitch(old_l.idm_dr, new_l.idm_dr, old_y, new_y, ctx);
  idr_->CopyPerturbationsAcrossSwitch(old_l.idr,    new_l.idr,    old_y, new_y, ctx);
}
```

(Same shape for `IDM_DRMD_IDR_DRMD_Species` delegating to `idm_drmd_` / `idr_drmd_`.)

Because the dispatch loop already iterates `all_species_` and calls each entry's override,
no new module wiring is needed — only the deletion of the now-redundant inline blocks.

## Behavior-change boundary

Two classes of change, verified differently:

- **Neutral relocations** — Baryons, CDM, Fluid, ScalarField overrides, and the `IDM_DR` child
  (whose override copies exactly the `delta`/`theta` its inline block copies today). These must be
  behavior-preserving: the new override copies exactly the slot set the inline code copies, no more.
  Verified against the **master** baseline within tolerance.
- **Intentional regression fix** — the `idr` / `idr_drmd` / `idm_drmd` reconduct restored by
  composite delegation. These copies do *not* happen on `master` today; wiring delegation adds them
  back. This is expected to change idm_dr/idm_drmd full-Cl output relative to `master`. Verified
  against the **pre-#268** reference (commit `d204e3a2^`, which performed these copies) using the
  new full-Cl scenarios — confirming we restore the correct pre-regression behavior, not invent new
  behavior. (Running the new scenarios on `master` is expected to differ or fail; that *is* the
  regression.)

## Verification

No bit-identical requirement (see project memory: this is a numerical ODE integration; the
rtol=1e-6 solver amplifies any ULP difference, and forcing zero drift produces contorted
order-preserving code for no physical gain). Acceptance is a numerically-meaningful tolerance.

Two baselines, captured before editing:

- **master baseline** (`d204e3a2`): build and run the five existing scenarios + the two new
  idm_dr/idm_drmd scenarios + `explanatory.ini`; save `.dat` outputs. (The new scenarios may NaN or
  differ here — that is the regression; record what happens.)
- **pre-#268 baseline** (`d204e3a2^`): build and run the two new idm_dr/idm_drmd scenarios; save
  outputs. This is the *correct* reference for the restored-copy behavior.

Cycle on the branch, run at default thread count (multi-thread is the point of the per-thread `pv`
design):

1. `rm -f output/scen_*.dat` (clear stale outputs — a real failure mode in past attempts).
2. `make class -j`.
3. Run each `test/scenarios/*.ini` (the five from #268 plus the two new idm_dr/idm_drmd scenarios),
   checking exit code per run, plus `./class explanatory.ini`.
4. Diff each `.dat` with the zero-crossing-aware criterion (TT/EE/Pk ≤ ~0.1% relative;
   Cℓ^TE/TPhi/Ephi via `|a−b| < atol + rtol·|a|`, `atol` at the column's ~1e-16 floor — never blind
   max-relative):
   - Five existing scenarios + simple-species changes → diff against **master** baseline.
   - New idm_dr/idm_drmd scenarios → diff against **pre-#268** baseline (restored behavior must
     match the pre-regression reference).
5. `TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py`.

The simple-species edits are relocations of slot *assignments* (no reductions reordered), so drift
vs master is expected at the noise floor. No order-preserving contortions.

## Risks

| Risk | Mitigation |
|---|---|
| Restored idr/idr_drmd/idm_drmd copy is wrong (over/under-copies vs pre-#268) | New full-Cl idm_dr/idm_drmd scenarios diffed against the **pre-#268** reference; must match within tolerance. |
| Simple-species override accidentally changes behavior | Each copies exactly its inline slot set, guarded by `idx >= 0`; diffed against master baseline. |
| New idm_dr/idm_drmd scenarios mis-specified (don't actually integrate perturbations through a switch) | Use `output=tCl,pCl,lCl,mPk`; confirm the run hits the approximation-switch path (verbose log shows switch messages) and that master vs pre-#268 differ on them (proving coverage). |
| CDM `idx_theta` / Fluid PPF-vs-standard mode mishandled | Guard every slot by `idx >= 0`; the layout already encodes which slots exist in the current gauge/mode. |
| Stale `output/scen_*.dat` masking a crash | `rm -f output/scen_*.dat` is step 1 of every cycle; exit code checked per scenario. |
| Touching DR / `N_decay_dr` by accident | DR composites explicitly out of scope; the only `N_decay_dr` reference in the inline switch block (DR multipole loop) is left untouched for PR B. |
