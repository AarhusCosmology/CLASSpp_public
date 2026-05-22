# PR #268 Follow-up Roadmap (A → C → B → D)

Tracks the four follow-up PRs spun out of #268 ("Encapsulate species perturbation layouts; remove
`ncdm_id`"). Agreed sequence: **A → C → B → D**, one PR at a time, reviewed/merged between each.
This doc is the cross-PR record so context survives between sessions; each PR gets its own
spec + plan when started.

| PR | Title | Status | Spec / Plan |
|----|-------|--------|-------------|
| **A** | Finish switch-copy migration (N_decay_dr-free species) | **Specced + planned; implementing** | `specs/2026-05-22-finish-switch-copy-migration-design.md`, `plans/2026-05-22-finish-switch-copy-migration.md` |
| **C** | Remove `pba->has_*` dispatch guards | **DONE** — 9 clean guards migrated; 3 idm_dr/drmd guards handed to B | `specs/2026-05-22-remove-has-dispatch-guards-design.md`, `plans/2026-05-22-remove-has-dispatch-guards.md` |
| **B** | Remove `pba->N_decay_dr` | Not started; **design gap** (setter pattern) | (to write) |
| **D** | Shooter hooks (encapsulate-layouts plan Tasks 37–40) | Not started | `plans/2026-05-07-encapsulate-species-layouts.md` (Tasks 37–40) |

Task 34 (`RescaledNCDMPerturbations(int)→(BaseSpecies*)`) is **moot**: the function had no callers and
was deleted entirely in #268.

---

## A. Finish switch-copy migration  *(in progress)*

Move the remaining inline always-reconduct copies in `perturb_vector_init` into species
`CopyPerturbationsAcrossSwitch` overrides. Scope: Baryons/CDM/Fluid/ScalarField overrides
(behavior-neutral) + wire `IDM_DR_IDR` / `IDM_DRMD_IDR_DRMD` composites to delegate to their existing
child overrides.

**Key finding (carries into B):** #268 wrote child `CopyPerturbationsAcrossSwitch` overrides for
`IDM_DR`/`IDR`/`IDM_DRMD`/`IDR_DRMD` but never wired the composites to call them, so the
`idr`/`idr_drmd`/`idm_drmd` reconduct was **dropped** on `master` (a latent regression, invisible
because no scenario integrates those perturbations through a switch). PR A restores it and adds two
new full-Cl scenarios (`test/scenarios/idm_dr_full.ini`, `idm_drmd_full.ini`) verified against the
**pre-#268** reference (`d204e3a2^`).

**Deferred to B:** the DR-bearing composites `DCDM_DR` and `DNCDM_DR`, and the DR multipole copy
(the `pba->N_decay_dr` loop in section a.2). Note `DCDM_DR`'s inline copy is *intact* on master (dcdm
delta/theta + DR multipoles still copied), so unlike idr it is **not** part of the dropped-reconduct
regression — B just needs to migrate it cleanly while removing `N_decay_dr`.

---

## C. Remove `pba->has_*` dispatch guards  *(DONE)*

The premise was badly stale. The referenced `plans/2026-04-15-remove-has-guards.md` **never existed**
in the repo (only referenced by other docs; the "PR 247-remove-has-guards" it alludes to is not on
GitHub either), and the intervening merges (#258/#260/#264/#266/#268/#269) had already removed nearly
all the `pba->has_*` dispatch guards. Reality on master @ 5741e002: of the surviving `pba->has_*`,
the `perturbations_module` ones (14) are all `class_test` *validation* (kept), `thermodynamics_module`
has none, `input_module` (21) are the flag *setters* (kept), and only **12 dispatch guards** remained
— in `nonlinear_module` (9), `output_module` (2), `background_module` (1).

PR C migrated the **9 behavior-neutral** ones (verified equivalences in the design doc):
`has_ncdm`×7 in `nonlinear` → `has_pk_cb_` (5 consumption sites) / `background_module_->GetNcdmCount() > 0`
(derivation) / dropped guard (Halofit warning loop self-filters); `has_fld` → `count("Fluid")`;
`has_scf` (`output`) → `count("ScalarField")`. No new files/helpers. Verified by the scenario suite
(84 scenarios + 6 reference tests) — behavior unchanged.

The **3 idm_dr/idr/drmd-family guards** (`background:532` `has_idr_drmd && has_idm_drmd`,
`nonlinear` `has_idm_dr`, `output` `has_idm_dr`) were **handed to B** (see below): the composites are
created on *either* sub-species, so a `count()` swap is not behavior-neutral, and there is no clean
species-based per-sub-species presence signal today.

---

## B. Remove `pba->N_decay_dr`

`N_decay_dr` is the count of DR-bearing decay channels (`(Omega0_dcdmdr>0?1:0) + n_dncdm`, set in
`input_module.cpp`). It is used as **both a count and an indexing offset** across `background_module`
(`index_bg_rho_dr_ = index_bg_rho_dr_species_ + pba->N_decay_dr`), `perturbations_module`,
`dcdm_dr_species.cpp`, and especially `dark_radiation_species.cpp` (multipole/background slot
indexing via the species' private `pba_` pointer).

**Design gap (resolve via brainstorming when started):** the `dark_radiation_species` consumers read
`pba_->N_decay_dr` through the species' private background pointer; removing the field needs a
"setter pattern" (give each DR species its own channel count/offset, set at construction/registration
time). Also note: **DR is only half-migrated** — `DRPerturbLayout` has `idx_F0_sum/idx_F0_species/
l_max`, but `perturb_vector` *still* has the parallel `index_pt_F0_dr_sum/index_pt_F0_dr_species/
l_max_dr/l_max_dr_col` fields; B should finish that. B also absorbs the DR-composite switch-copy
migration deferred from A (DCDM_DR / DNCDM_DR + DR multipole copy).

When B reworks this section, also merge the two consecutive `if (all_species_.count("DCDM_DR"))`
blocks in `perturb_vector_init`'s always-reconduct region into one (flagged by the #269 Copilot
review — left as-is in PR A because that block is the deferred N_decay_dr code being rewritten here).

**Inherited from C — the 3 idm_dr/idr/drmd-family dispatch guards.** B now also owns
`background_module.cpp:532` (`if (pba->has_idr_drmd && pba->has_idm_drmd)`), the `pba->has_idm_dr`
guard in `nonlinear_module.cpp` (the Halofit/HMcode "you requested idm_dr" warning), and the
`pba->has_idm_dr` guard in `output_module.cpp` (the `dmu_idm_dr` thermo-file header). These can't be
swapped for `all_species_.count("IDM_DR_IDR")` / `count("IDM_DRMD_IDR_DRMD")` because the composites
are created on `has_idm_dr || has_idr` / `has_idm_drmd || has_idr_drmd` (so `count()` ⊋ the flag, and
divergence is observable for idr-only / drmd-partial models — `test_idr_without_idm_dr_computes`,
`test_drmd_without_idr_drmd_computes`). `IDM_DR_IDR_Species` even builds *both* children
unconditionally and distinguishes them by reading `pba->has_idm_dr`/`has_idr` internally. Migrating
them needs **per-sub-species presence accessors** on the composites (e.g. `bool has_idm_dr() const`
capturing the flag at construction) — the same setter-pattern work this PR already requires.

---

## D. Shooter hooks (Tasks 37–40)

Move species-specific shooting cases out of `InputModule` into per-species
`RegisterShootingIndices` / `ComputeShootingResidual` / `ComputeShootingGuess` hooks + a
`ShootingContext`; cosmological `theta_s` stays module-level. Per the encapsulate-layouts plan
Tasks 37–40, the species are `DNCDM_DR`, `DCDM_DR`, `ScalarField`.

**Reality vs the plan sketch:** the plan idealizes the shooter. The real `InputModule` uses an
enum-based dispatch (`enum target_names`: `theta_s`, `Omega_dcdmdr`, `Omega_scf`, `Omega_ini_dcdm`,
`Omega_dncdmdr`, `Omega_ini_dncdm`) over a `fzerofun_workspace`, spread across
`input_try_unknown_parameters` (residual), `input_get_guess` (guess), and a default-target-value
helper. D needs a real design pass against that machinery, not a verbatim application of the sketch.
