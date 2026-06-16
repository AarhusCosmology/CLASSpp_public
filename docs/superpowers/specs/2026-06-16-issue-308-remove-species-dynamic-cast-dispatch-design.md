# Remove `dynamic_cast` species dispatch from modules (#308)

- **Issue:** [#308](https://github.com/AarhusCosmology/CLASSpp/issues/308)
- **Date:** 2026-06-16
- **Status:** Design — reviewed (HasNcdm deferred per review); ready for implementation plan
- **Branch:** `308-remove-species-dynamic-cast-dispatch`

## Problem

The stated plugin rule is "module code loops + dispatches, never downcasts." Modules
still violate it via `dynamic_cast` in ~20 sites that type-pick NCDM- and
dark-radiation-family species (`NCDMBaseSpecies`, `NCDMSpecies`, `DCDM_DR_Species`,
`DNCDM_DR_Species`, `DarkRadiationSpecies`), plus the `NcdmFamily` /
`GetNcdmSpecies` helpers that exist only to reach *into* a composite for its child.

## Design principle (what this PR is really about)

The composite-species design is: take uncoupled species and add couplings; the
composite owns the coupling, each child owns its free-streaming. A composite is a
`BaseSpecies` like any other — the rest of the code should ask *it* a question and
get an answer, never reach through it to its children.

The codebase already follows this in the hot paths:

- `background_module.cpp:1030` — `Omega_rad += sp->GetRadiationOmega0();`
- `perturbations_module.cpp:3789` — `rho_nu += sp->FreestreamingRho(...);`
- `perturbations_module.cpp:5108` — `mrd = sp->MatterRhoDelta(...);` (ask the
  composite for δρ_m; it answers like any species)

Each of these is a **neutral-default polymorphic method**: plain species that don't
participate return 0 / no-op; composites forward to / sum over their children
internally. **#308 extends this existing pattern to the remaining type-querying
sites.** It is therefore mostly *deletion* of dispatch machinery, not new API.

Explicitly rejected during design:

- **`As<Type>()` capability queries** (`AsNcdm()`, `AsDarkRadiation()`, …) — these
  are `dynamic_cast` wearing a polite hat; they keep type-dispatch in the modules.
- **A module-level `children()` accessor / flatten helper** — this reaches through
  the composite's encapsulation, which is exactly the split we are protecting.

## The two buckets

Every offending site is one of two kinds:

### Bucket 1 — the composite simply fails to forward (the real work)

The module reaches around the composite (via `dynamic_cast` or `NcdmFamily`) only
because the composite doesn't forward a call to its child. Canonical example,
`perturbations_module.cpp:2734`:

```cpp
if (evolve_tensor_ncdm_ == _TRUE_)
  for (auto* sp : NcdmFamily(all_species_))      // reaches THROUGH DNCDM_DR to its child
    sp->WriteTensorOutputColumnTitles(tensor_titles_);
```

`WriteTensorOutputColumnTitles` is *already* a `BaseSpecies` virtual (base_species.h:306,
default no-op). The fix is to make `DNCDM_DR_Species` forward it to its `dncdm_`
child (it already forwards `WriteBackgroundColumnTitles`, `WriteBackgroundData`, …)
and then loop all species:

```cpp
if (evolve_tensor_ncdm_ == _TRUE_)
  for (auto& sp : all_species_)
    sp->WriteTensorOutputColumnTitles(tensor_titles_);   // non-NCDM no-op; DNCDM_DR forwards
```

No cast, no `NcdmFamily`, no module-level `children()`.

### Bucket 2 — genuine subset enumeration (DELETED)

The only true "give me the list of NCDM species, indexed" need was the public
introspection helper `Class.get_ncdm()` (classy.pyx:1068-1097), which dumps each
NCDM species' `deg` / `m_ncdm` / `q_size` / `q[i][j]` into a dict. It has **zero
callers anywhere in the repo** (no tests, notebooks, or examples) and is a
long-obsolete leftover.

**Decision:** delete `get_ncdm()` and its four backing accessors
(`GetNcdmQ`, `GetNcdmDeg`, `GetNcdmMassInEV`, `GetNcdmQSize` in
`background_module.{h,cpp}`). `cclassy.pxd` regenerates on install, so no manual
binding edit. This removes Bucket 2 entirely — after it, *every* remaining site is
Bucket 1.

## API changes

### Methods made to forward in composites

For each `BaseSpecies` virtual currently reached via `NcdmFamily`/`dynamic_cast`,
add the missing forward in the relevant composite (`DNCDM_DR_Species` to its
`dncdm_`/`dr_sp_` children; `DCDM_DR_Species` to its `dr_sp_`):

- `WriteTensorOutputColumnTitles` — forward to dncdm child (DNCDM_DR).

### New neutral-default `BaseSpecies` methods (default 0 / no-op)

Added where no existing virtual fits; each composite overrides to forward/sum:

- `double DarkRadiationRhoToday(const double* pvecback_integration) const` —
  default 0; `DarkRadiationSpecies` returns `pvecback_integration[bi_rho_index()]`;
  `DCDM_DR`/`DNCDM_DR` forward to their `dr` child. Replaces the `Omega0_dr`
  cast cluster (background_module.cpp:858-867). (Final name/shape TBD in plan —
  may reuse the integrated-vector convention already used nearby.)
- `double NeffContribution(double z) const` — default 0; NCDM returns `GetNeff(z)`;
  composites forward. Replaces the N_eff accumulation loop
  (background_module.cpp:503-508).
- `void PrintNeffInfo() const` / `void PrintMassInfo() const` — promote to
  `BaseSpecies` no-op defaults (already exist on `NCDMBaseSpecies`); composites
  forward. Replaces the verbose-print loops (background_module.cpp:506, 551-553).
- `double NeutrinoOmega0() const` — default 0; NCDM returns `GetOmega0()` of the
  NCDM sector; composites forward to the ncdm child only. Replaces
  `GetOmega0NcdmTot` (background_module.cpp:1531-1536), which must *exclude* a
  DNCDM_DR composite's DR child (so it cannot use the composite's own
  `GetOmega0()`).
- A relativistic-IC check hook and a Halofit-WDM-warning hook (small no-op-default
  predicates/methods) for perturbations_module.cpp:2466-2510 and
  nonlinear_module.cpp:917-928. Exact shape decided in the plan; both default to
  no-op and composites forward to the ncdm child.

### Deletions

- `get_ncdm()` (classy.pyx); `GetNcdmQ`/`GetNcdmDeg`/`GetNcdmMassInEV`/`GetNcdmQSize`
  (background_module.{h,cpp}). `GetNcdmCount` stays — still used by
  nonlinear_module.cpp:1249.
- `NcdmFamily` (perturbations_module.cpp:73-82) and `GetNcdmSpecies`
  (background_module.cpp:95-104) are **not** deleted in #308. Their
  contribution-forwarding consumers are converted (so each helper shrinks to its
  detection/count uses), and the helpers themselves are deferred to #309/#310 — see
  the decision below.

## Decision: NCDM detection stays (deferred to #309/#310)

`HasNcdm(all_species_)` (perturbations_module.cpp:51) is called in ~25 places, almost
all boolean gates around NCDM-specific perturbation machinery that is **not yet
species-encapsulated** (the fluid-approximation index `index_ap_ncdmfa` at 2266,
RSA handling, tensor-evolution flags `evolve_tensor_ncdm_`, etc.). Fully removing
the module's *need* to know "is there NCDM" means moving that machinery into the
species — that is #309/#310 work.

**Decision (reviewed): leave the NCDM detection/enumeration layer untouched in #308
and fold it into #309/#310.** These stay exactly as-is, internal `dynamic_cast`
included:

- `HasNcdm` and all ~25 boolean-gate call sites.
- `NcdmFamily(...).size()` (perturbations_module.cpp:1088, 1107) and `DrSpeciesCount`
  (812) — counts that feed `index_tp_ncdm` / `index_tp_dr_` sizing (#309 territory).
- `GetNcdmSpecies(...).empty()` gates (background_module.cpp:498, 503, 551, 1362, 1390).
- The NCDM background-index registration loop (background_module.cpp:633-646), which
  is entangled with index ordering and `bg_pseudo_p_index` (#309).

What #308 *does* remove is the type-dispatch **inside** these gated regions — the
per-species `dynamic_cast` contribution loops, which dissolve into
`for sp : all_species_` + a neutral-default forwarded method (see the API section).
The outer `HasNcdm` gate is preserved purely as an optimization (skip the loop when
no NCDM is present); the loop is correct without it. Example — the relativistic IC
check (perturbations_module.cpp:2466):

```cpp
if (HasNcdm(all_species_)) {                  // gate kept (deferred to #309/#310)
  for (auto& sp : all_species_)
    sp->AssertRelativisticAtIc(ppw->pvecback, ppr->tol_ncdm_initial_w);  // no-op default; composite forwards to child
}
```

## Scope boundaries (with reasons)

- **Source-slot assignment loops** (perturbations_module.cpp:642-664, `SetSourceSlot`)
  → **#309**. They *are* the source-type mechanism #309 replaces.
- **Tensor-GW-source sites needing the nested per-pv layout**
  (perturbations_module.cpp:5285-5293, 5918-5934) → **#310**. The `dynamic_cast`
  there is only to fetch the child's *nested* `PerturbLayout` (e.g.
  `comp_lay.dncdm`); de-typing it needs per-species layout slots, which is #310.
  `ContributeTensorGwSource` is already a forwarded virtual; only the layout access
  blocks it.
- **IDM `has_idm_dr` / `has_idr_drmd` by-name `static_cast`** (nonlinear_module.cpp:929,
  background_module.cpp:524-531) → **out of scope** (different sector; not NCDM/DR).
- **`pba->has_X` / `pba->Omega0_X`** → already removed from master; the stale
  `275-remove-has-flags` branch is superseded and irrelevant here.

## Behaviour changes

All stem from composites now forwarding to children, so a wrapped **DNCDM child**
participates where the old narrow `dynamic_cast<NCDMSpecies*>` excluded it. These
affect only runs containing a `DNCDM_DR` species:

1. **Relativistic IC check** (perturbations_module.cpp:2466-2510) now also asserts the
   decaying NCDM child is ultra-relativistic at `tau_ini`. Physically it should be;
   this is an *assertion* that could newly fail an existing DNCDM run — verify
   against the DNCDM-DR tests; if it breaks, keep this one site top-level-only.
2. **Tensor-GW massless-approximation `3·P` term** (perturbations_module.cpp:5263-5271)
   now includes the child — matching what the *exact* branch (5285) already does.
   Removes an existing inconsistency.
3. **Halofit WDM-mass warning** (nonlinear_module.cpp:917-928) now also warns for a
   heavy decaying child. Benign (warning only).

This is consistent with the design: the composite decides whether its child
participates, not the module.

## Testing

- Full build + the standard test suite, with emphasis on the scenarios that
  exercise these paths: NCDM (massive neutrinos), DNCDM_DR closure, DCDM_DR, and
  dark-radiation runs, in both gauges.
- Confirm **no output diffs** on non-DNCDM scenarios (the forwarding is a no-op
  there).
- Confirm the three flagged behaviour changes behave as intended on DNCDM_DR runs;
  in particular that the relativistic-IC assertion still passes for the existing
  DNCDM-DR test inputs.
- Confirm `classy` still imports and builds after `get_ncdm()` removal.

## Relationship to other v4 issues

- **#309** (species-registered source types) builds on the cleaned base; it owns the
  `index_tp_*` / `SetSourceSlot` / `has_source_*` restructuring and may absorb
  `HasNcdm` if we leave it for then.
- **#310** (per-species workspace/layout) unblocks the two deferred tensor-GW sites.
