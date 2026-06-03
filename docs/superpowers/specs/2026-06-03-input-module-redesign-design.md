# Input Module Redesign — Design

**Date:** 2026-06-03
**Status:** Approved (ready for implementation plan)

## Problem

`source/input_module.cpp` does not follow a clear flow. The intended mental model is:

1. Read non-species-specific quantities needed to build the species construction context.
2. Construct species — each species reads its own inputs.
3. Read additional inputs that may depend on the constructed species.

The current code approximates this but is turned inside-out: **construction happens last**, so
everything that should be phase (iii) is forced to run in phase (i), before any species exists.
That inversion is the root of several workarounds. The constructor today is:

```cpp
input_init();        // input_read_precisions(), input_read_parameters() [~2100 lines], write-params dump
ConstructSpecies();  // builds all species, bolted on at the very end
```

### Blockers identified

1. **Construction runs once, at the end.** Species-dependent reads cannot see the species, so they
   fake it: `S8→sigma8` reads `omega_budget_.cdm` instead of the CDM species
   (`input_module.cpp:1689`); the Halofit `pk_eq` gate *peeks* `wa_fld`/`fluid_equation_of_state`
   from `pfc` because the Fluid species is not built yet (`input_module.cpp:2693–2719`).
2. **The "write parameters" dump runs before construction.** `input_init()` emits the read/unread
   parameter report at `input_module.cpp:547–586` and returns; `ConstructSpecies()` only runs
   afterwards. Species inputs (`N_ur`, `Omega_scf`, `m_ncdm`, dot-syntax ncdm fields, …) are marked
   read inside `CreateAll`, which has not executed at dump time — so `write parameters = yes`
   mis-reports every species parameter as **unused**. Concrete correctness bug caused purely by the
   ordering.
3. **Photons/Baryons density is read in the monolith, not by their factory.** `PhotonsSpecies` and
   `BaryonsSpecies` `CreateAll` read nothing — they wrap `*ctx.pba`; the `T_cmb`/`Omega_g`/`omega_g`
   and `Omega_b` parsing lives at `input_module.cpp:758–810`. This is *not* a misplacement to fix by
   forcing self-parsing: `Omega0_g` is genuine shared background context (see "Photon density"
   below). The fix is to give phase (i) an explicit home and single-source the `T_cmb⇄Omega0_g`
   formula on the species class.
4. **Coupled-cluster params are parsed twice.** `ReadCoupledOmegaBudget()` parses `T_idr` from
   `stat_f_idr`/`N_idr`/`N_dg`/`xi_idr` to derive `omega_budget_.idr`, then *discards* `T_idr`
   (`input_module.cpp:305–353`). `IDM_DR_IDR_Species::CreateAll` re-parses the identical block
   (`idm_dr_idr_species.cpp:437–456`) because the species constructor needs `T_idr` as physics. The
   budget keeps only the Ω it computed; nothing carries the intermediate.
5. **IDM/IDR physics params live on `thermo`/`perturbs`, read in the monolith.** `a_idm_dr`,
   `nindex_idm_dr`, `m_idm`, `b_idr` (on `thermo`) and `idr_nature`, `alpha_idm_dr`, `beta_idr` (on
   `perturbs`) are parsed at `input_module.cpp:841–962`. They are consumed by `thermodynamics_module`
   and `perturbations_module`, and the `IDM_DR_IDR` species reaches back into the globals for its
   *own* data (`idm_dr_idr_species.cpp:372`, `:388`). The species owns the behavior but not the data.

## Decisions (locked)

- **Blocker 5 scope:** full ownership migration of the seven `IDM_DR_IDR` params onto the species.
- **Access pattern:** typed accessors on a **sanctioned `IDM_DR_IDR` downcast exception**. The
  thermo and perturbation modules already do `static_cast<IDM_DR_IDR_Species&>(...)`
  (`thermodynamics_module.cpp:203/506/830`, `perturbations_module.cpp:3659`); this formalizes that
  exception with typed getters (`comp.idm_dr().a_idm_dr()`, `comp.idr().idr_nature()`), mirroring the
  existing `comp.idr().T_idr()` precedent.
- **Photon density stays on `pba`.** `pba->Omega0_g` is a shared background scalar read by the
  background module itself (`Neff = Omega0_ur/… / pba->Omega0_g` at `background_module.cpp:537`,
  `N_dark` at `:549`, `Omega_rad` at `:1105`, budget print `:1436`) and by `FluidSpecies`
  (`fluid.cpp:333`) and `UltraRelativisticSpecies` (`ultra_relativistic.cpp:460`). `PhotonsSpecies`
  itself is a thin reader of `pba_.Omega0_g` (`photons.h:30/38`). The value cannot move onto the
  photon object without creating a second source of truth, and `Neff` is inherently a cross-species
  ratio with no single-species home. The asymmetry "photons read `pba`, UR reads `pfc`" reflects a
  real difference and is kept. We single-source only the **formula** (below).

## Design — Approach A: three named phases + one shared coupled-inputs carrier

### Control flow

```cpp
InputModule::InputModule(FileContent& fc) : file_content_(fc) {
  file_content_.mark_all_unread();
  try {
    input_read_precisions();   // unchanged
    input_default_params();    // all defaults up front (cheap)
    ReadContext();             // phase i  — only what building species needs
    ConstructSpecies();        // phase ii — registry loop (mostly unchanged)
    ReadDerived();             // phase iii — everything else + species-dependent reads
    WriteParameterFiles();     // moved here → sees a fully-read pfc (fixes Blocker 2)
  } catch (const std::runtime_error& e) {
    throw std::invalid_argument(e.what());
  }
}
```

- `input_read_parameters()` is **dissolved** into `ReadContext()` + `ReadDerived()`; nothing else
  calls it. `input_init()` is removed (its orchestration moves to the constructor; the
  `write parameters`/`write warnings` dump becomes `WriteParameterFiles()`).
- Validation tests move with their inputs: e.g. "cannot enter both `h` and `100*theta_s`" → phase i;
  "only one of `A_s`/`sigma8`/`S8`" → phase iii. The `runtime_error→invalid_argument` wrapper spans
  all three phases.

### Phase i — `ReadContext()`

Owns exactly the inputs that `SpeciesBuildContext` and the budget consume:

- `input_verbose`, `threads`
- `gauge`
- `a_today`
- `h` / `H0` / `100*theta_s` (theta_s consumed here; shooting itself handled by `DoShooting`)
- `T_cmb` / `Omega_g` / `omega_g` → `pba->Omega0_g`, via a new **static** helper on the photon
  species (see "Photon formula")
- `Omega_b` / `omega_b` → `pba->Omega0_b`
- `Omega_k` → `pba->Omega0_k`, `K`, `sgnK`
- Coupled-cluster parse → fills `omega_budget_` **and** `coupled_inputs_` (Blocker 4)
- Closure selection → `pba->closure_species`

Approximately the current `input_module.cpp:688–1044`, minus the parts that migrate to species
(Blocker 5) and minus the `fluid_present_pfc`/`wa_fld_peek` snapshots (no longer needed — phase iii
asks the Fluid species directly).

### Phase ii — `ConstructSpecies()`

Mostly unchanged. The build context gains `coupled_inputs` (Blocker 4). Coupled factories read
intermediates from it instead of re-parsing `pfc`. `IDM_DR_IDR_Species::CreateAll` additionally
parses the migrated interaction params and stores them on the species (Blocker 5). The closure
two-pass and the NCDM precision-consistency check are unchanged.

### Phase iii — `ReadDerived()`

Owns everything from the thermo block onward:

- thermo: `YHe`, recombination, reionization parametrizations + params, energy injection, compute
  damping scale
- perturbation/output configuration (the large species-independent block, current
  `input_module.cpp:~1255–2505`), relocated intact
- primordial (`A_s`/`ln10^{10}A_s`/`sigma8`/`S8`): the **S8 branch asks the CDM species**
  (`all_species_.count("CDM") ? all_species_.at("CDM").GetOmega0() : 0` + `Omega0_b`) instead of
  `omega_budget_.cdm`
- transfer, spectra, lensing configuration
- tensor method; trigger-consistency tests
- **Halofit `pk_eq` gate asks the Fluid species** for its EoS type and `wa_fld` instead of peeking
  `pfc`; the `fluid_equation_of_state`/`wa_fld` peek and `fluid_present_pfc` snapshot are deleted
- `l_max` consistency tests (depend on `ppr` only)

`WriteParameterFiles()` runs after this, so the read/unread report is now correct.

### Photon formula (Blocker 3, lightweight capture)

Add a `static` method pair on `PhotonsSpecies` that owns the `T_cmb ⇄ Omega0_g` conversion
(currently the `sigma_B` block at `input_module.cpp:766–798`), e.g.
`PhotonsSpecies::Omega0gFromTemp(T_cmb, h)` and its inverse. `ReadContext()` calls it and writes the
result to `pba->Omega0_g`. The species class owns the math (single source of truth, no duplicated
`sigma_B` expression); the value stays where the background module needs it; no redundant live
object and no registry special-case.

### Blocker 4 — `CoupledClusterInputs` carrier

New struct alongside `SpeciesOmegaBudget` in `species/species_build_context.h`:

```cpp
// Raw intermediates parsed once in phase i that the coupled factories need for
// physics construction (not just budget math). Lets a factory read e.g. T_idr
// instead of re-parsing N_idr/N_dg/xi_idr identically to the budget resolver.
struct CoupledClusterInputs {
  std::optional<double> T_idr;       // resolved from N_idr | N_dg | xi_idr (+ stat_f_idr)
  double stat_f_idr = 7. / 8.;
  // DRMD interaction physics: parsed once by the budget resolver, consumed by
  // IDM_DRMD_IDR_DRMD_Species::CreateAll (previously re-parsed there verbatim).
  double f_idm_drmd = 0.;
  double G_over_aH_drmd_ini = 0.;
  double delta_Neff_drmd = 0.;
  double z_stop = 0.;
};
```

- `ReadCoupledOmegaBudget()` is renamed to `ReadCoupledCluster()` and fills **both**
  `omega_budget_` and a new `coupled_inputs_` member; the `T_idr` and DRMD values it already parses
  (currently discarded after the budget math) are retained.
- `SpeciesBuildContext` gains `const CoupledClusterInputs* coupled_inputs = nullptr;`.
- `IDM_DR_IDR_Species::CreateAll` deletes its duplicate `stat_f_idr`/`N_idr`/`N_dg`/`xi_idr` block
  and reads `ctx.coupled_inputs->T_idr`.
- `IDM_DRMD_IDR_DRMD_Species::CreateAll` deletes its `read_double` block for `f_idm_drmd`/
  `G_over_aH_drmd_ini`/`delta_Neff_drmd`/`z_stop` (`idm_drmd_idr_drmd_species.cpp:392–395`) and reads
  them from `ctx.coupled_inputs`.

**Carrier covers all coupled double-parses.** Both `T_idr` and the four DRMD keys
(`z_stop`/`G_over_aH_drmd_ini`/`f_idm_drmd`/`delta_Neff_drmd`) are currently parsed in the budget
resolver *and* re-parsed in their factory (`input_module.cpp:444–449` vs.
`idm_drmd_idr_drmd_species.cpp:392–395`). The carrier eliminates both duplications and is the single
source of truth for the coupled cluster's raw intermediates. Presence-based validation
(`any_drmd && !all_drmd`) stays in the budget resolver, which holds the parse flags.

### Blocker 5 — ownership migration

`IDM_DR_IDR` is the **only** consumer of the seven global fields. The DRMD sector
(`IDM_DRMD_IDR_DRMD`) uses **none** of them — it is a physically different model (`G/aH` interaction
shutting off at `z_stop`) and already self-owns its params (`f_idm_drmd`, `G_over_aH_drmd_ini`,
`delta_Neff_drmd`, `z_stop`) by reading `pfc` in its own `CreateAll`. So the migration is total, not
partial: after rewiring, the seven fields are **deleted outright** from `thermo`/`perturbs`.

Fields and new homes (each param on its natural child, reached through the composite, mirroring
`comp.idr().T_idr()`):

| Field | From | To | Accessor |
|---|---|---|---|
| `a_idm_dr`, `nindex_idm_dr`, `m_idm` | `thermo` | `IDM_DRSpecies` | `comp.idm_dr().a_idm_dr()` … |
| `b_idr` | `thermo` | `IDRSpecies` | `comp.idr().b_idr()` |
| `idr_nature`, `alpha_idm_dr` (+storage), `beta_idr` (+storage) | `perturbs` | `IDRSpecies` | `comp.idr().idr_nature()` … |

**Parsing** moves out of the deleted monolith block (`input_module.cpp:841–962`) into
`IDM_DR_IDR_Species::CreateAll`, which already parses the sector. It reads `T_idr` from
`coupled_inputs`, parses the interaction params from `pfc` (including the `a_idm_dr`/`a_dark`/
`Gamma_0_nadm` aliasing and the ETHOS/NADM `input_verbose` printouts — the factory reads
`input_verbose` from `pfc`).

**Consumer rewrites** (all already hold, or already downcast to, the composite):

- `thermodynamics_module.cpp` (`:209`, `:221`, `:242`, and the `Gamma_heat_idm_dr` block
  `~:878–1041`): `pth->a_idm_dr` → `comp.idm_dr().a_idm_dr()`, etc. `comp` is already bound at
  `:203`/`:506`/`:830`.
- `perturbations_module.cpp`: the four `ppt->idr_nature` sites (`:2301`, `:3658`, `:4638`, `:5265`).
  The decisive one is `:5265` (`ppw->scalar_ctx.idr_nature = ppt->idr_nature;`) — re-source from the
  composite; the deep hot-path consumers reading `ppw->scalar_ctx.idr_nature`
  (`interacting_species.cpp:135/264/326`) then need **no change**.
- `idm_dr_idr_species.cpp`: its own reads (`:52`, `:166`, `:197`, `:372`, `:388`) become reads of
  its **own members** — the species stops reaching into globals for its own data.

**Highest-risk item: `idr_nature`** (eight read sites across two modules plus the species). The
`scalar_ctx` indirection contains the blast radius, but this sub-task lands last and is verified in
isolation.

## Error handling

Unchanged patterns: `class_test`/`class_call` and `throw std::invalid_argument`. Each validation
moves to the phase where its inputs are available. The single `runtime_error→invalid_argument`
wrapper now spans precisions + all three phases.

## Verification

Behavior-preserving refactor → characterization-based verification, **~0.1% agreement, not
bit-identical**, with care around Cl^TE zero-crossings (per standing project guidance).

**Characterization fixtures**, captured *before* any change, covering every branch that moves
between phases: ΛCDM; +massive ncdm; +idm_dr (ETHOS *and* NADM/`Gamma_0_nadm` paths); +idm_drmd;
+dcdm; +scf; +fld; `S8` input; `100*theta_s` shooting.

**Landing sequence** — three independently-shippable, independently-verifiable steps, low→high risk:

1. **Reordering only.** Split into `ReadContext`/`ConstructSpecies`/`ReadDerived`, move the dump,
   convert the S8/Halofit peeks to "ask the species." No data moves. *Verify:* full reference-output
   comparison + a `write parameters = yes` run confirming species params now land in
   `parameters.ini`, not `unused_parameters` (the Blocker 2 proof).
2. **Blocker 4.** `CoupledClusterInputs` carrier; delete the duplicate `T_idr` parse and the
   duplicate DRMD-key parse. *Verify:* IDM_DR/IDR *and* IDM_DRMD scenarios stable (same numbers from
   one parse instead of two).
3. **Blocker 5.** Field migration, sub-staged: thermo params first (`a_idm_dr`/`nindex`/`m_idm`/
   `b_idr`), then `alpha`/`beta`, then `idr_nature` last and alone. *Verify after each sub-stage:*
   the IDM-DR-interacting scenarios specifically (the only hot-path touch).

## Out of scope

- Eliminating `pba->Omega0_g` / routing the background module's radiation reads through the species
  collection (the "modules loop + dispatch, species own all densities" end-state). Noted as a future
  direction; `Neff` resists single-species ownership.
- Any change to the `IDM_DRMD_IDR_DRMD` sector's **physics/runtime** — it shares none of the seven
  global fields and is untouched by the Blocker 5 migration. (Its `CreateAll` *is* touched by
  Blocker 4, to read the four DRMD keys from `coupled_inputs` instead of re-parsing `pfc` — a
  parse-site change only.)
