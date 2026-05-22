# Remove `pba->has_*` dispatch guards — design (PR C of #268 follow-ups)

## Goal

Replace species-presence **dispatch** guards written as `pba->has_*` in module code with
queries against the species container (`all_species_` / existing module accessors), so module
code stops type-picking on background flags. Behavior-neutral: the executed code paths must be
identical, so the existing scenario suite must pass unchanged.

## Why the roadmap is stale

The roadmap (`plans/2026-05-22-pr268-followups-roadmap.md`, section C) describes "~97 `pba->has_*`
dispatch guards in `background_module.cpp` / `perturbations_module.cpp` / `thermodynamics_module.cpp`"
and points at `plans/2026-04-15-remove-has-guards.md`. Both are stale:

- **That plan file does not exist** in the repo — it is only *referenced* by three other docs
  (and one mentions a "previous PR 247-remove-has-guards" that is also not on GitHub). There is
  nothing to "refresh"; this design replaces it.
- The intervening merges (#258, #260, #264, #266, #268, #269) already removed the bulk of the
  `pba->has_*` dispatch guards as a side effect of the species-container work.

## Full `pba->has_*` inventory (master @ 5741e002)

| File | count | kind |
|---|---|---|
| `input_module.cpp` | 21 | **setters** — `pba->has_cdm = (Omega0_cdm != 0) ? ...` etc. These *define* the flags; keep (flags stay on the struct). |
| `perturbations_module.cpp` | 14 | **validation** — all 14 are inside `class_test(...)` (unsupported-combination errors). Keep per roadmap. |
| `thermodynamics_module.cpp` | 0 | — |
| `nonlinear_module.cpp` | 9 | **dispatch** — `has_ncdm`×7, `has_idm_dr`×1, `has_fld`×1 |
| `output_module.cpp` | 2 | **dispatch** — `has_scf`×1, `has_idm_dr`×1 |
| `background_module.cpp` | 1 | **dispatch** — `has_idr_drmd && has_idm_drmd` (`:532`) |

So the only remaining *dispatch* guards are 12, in `nonlinear`/`output`/`background` — not the
three modules the roadmap named.

## In scope for PR C: the 9 behavior-neutral guards

Each replacement below is provably equivalent because the target species/flag is created/derived
under exactly the same condition.

| Site | Current | Replacement | Equivalence |
|---|---|---|---|
| `nonlinear:1007` | `if (pba->has_ncdm) { for (...) {dynamic_cast<NCDMSpecies*>...} }` | **drop the `if`**, keep the loop | The inner loop already self-filters via `dynamic_cast<NCDMSpecies*>`; it is a no-op when no NCDM species exist. Removes the guard entirely (loop + dispatch, no gate). |
| `nonlinear:1383` | `if (pba->has_ncdm == _TRUE_) has_pk_cb_=_TRUE_; else _FALSE_;` | `background_module_->GetNcdmCount() > 0` | `GetNcdmCount()` (public, `background_module.h:25`) counts `NCDMBaseSpecies` + wrapped `DNCDM_DR_Species`, i.e. exactly the set that makes `has_ncdm == _TRUE_` (`input_module.cpp:347-349`). |
| `nonlinear:2695` | `if (pba->has_ncdm == _TRUE_)` | `if (has_pk_cb_ == _TRUE_)` | `has_pk_cb_` is set `== has_ncdm` at 1383 and gates `index_pk_cb_` (defined only when `has_pk_cb_`, `:1391`); other sites already gate this data on `has_pk_cb_` (`:1732`, `:2135`). More correct semantically (the site consumes pk_cb). |
| `nonlinear:3583/3647/3711/3775` | `if (pba->has_ncdm)` | `if (has_pk_cb_)` | Same as 2695 — these read `pnw->sigma_*[index_pk_cb_]`, valid only when `has_pk_cb_`. |
| `nonlinear:3173` | `if (pba->has_fld == _TRUE_)` | `if (all_species_.count("Fluid"))` | `FluidSpecies::CreateAll` pushes `"Fluid"` iff `has_fld == _TRUE_` (`fluid.cpp:512`). |
| `output:997` | `if (pba->has_scf == _TRUE_)` | `if (all_species_.count("ScalarField"))` | `ScalarFieldSpecies::CreateAll` pushes `"ScalarField"` iff `has_scf == _TRUE_` (`scalar_field.cpp:371`). |

No new helpers, headers, or files are required.

## Deferred to PR B: the 3 idm_dr/idr/drmd-family guards

`background:532` (`has_idr_drmd && has_idm_drmd`), `nonlinear:1021` (`has_idm_dr`),
`output:1047` (`has_idm_dr`) are **not** behavior-neutral as a `count()` swap:

- The composites are created when **either** sub-species is present:
  `IDM_DR_IDR::CreateAll` fires on `has_idm_dr || has_idr` (`idm_dr_idr_species.cpp:414`);
  `IDM_DRMD_IDR_DRMD::CreateAll` fires on `has_idm_drmd || has_idr_drmd` (`idm_drmd_idr_drmd_species.cpp:368`).
  So `count("IDM_DR_IDR")` ⊋ `has_idm_dr` and `count("IDM_DRMD_IDR_DRMD")` ⊋ `has_idr_drmd && has_idm_drmd`.
- The divergence is observable for idr-only / drmd-partial models, which have dedicated tests
  (`test_idr_without_idm_dr_computes`, `test_drmd_without_idr_drmd_computes`).
- There is no clean species-based "is idm_dr individually present" signal: `IDM_DR_IDR_Species`'s
  constructor **unconditionally** builds both children and then distinguishes them by reading
  `pba->has_idm_dr`/`has_idr` *internally* (`idm_dr_idr_species.cpp:32,38,44`). `pba->has_*` is the
  actual source of truth.

Migrating these needs per-sub-species presence accessors on the composites — the same
"setter pattern" design work the roadmap assigns to **PR B** (the idr/drmd/`N_decay_dr` machinery).
Doing it here would pre-empt and duplicate B. PR C therefore leaves these three reads in place and
hands them to B.

## Verification

The changes are control-flow-equivalent, so output should be unchanged.

- `make class -j` — builds with no new warnings.
- `make classy` — rebuild the wrapper.
- `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py`
  — the per-PR CI gate (`build.yml`). The cosmology grid exercises ncdm, fld, and scf paths.
- The dedicated reference tests (`test_idm_dr_idr_perturbations_match_reference`,
  `test_dcdm_dr_matches_reference`, `test_idr_without_idm_dr_computes`,
  `test_drmd_without_idr_drmd_computes`) confirm the deferred-to-B species are untouched.

No new `.cpp`/`.h` files → no Makefile / setup.py / pbxproj changes.
