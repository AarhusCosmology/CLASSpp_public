# NCDM-family key-name-driven input

**Status:** design  
**Date:** 2026-04-30  
**Issue:** simplification follow-up to #246 / PR #260 (object-based NCDM input)

## Motivation

PR #260 introduced object-based (dot-syntax) input for NCDM species, but kept full
backward compatibility with the legacy comma-separated CSV form. The compatibility
glue is heavy:

- `NCDMSpecies::CreateAll` runs **dot → CSV** synthesis (`synthesise_standard_ncdm_flat_keys`),
  then constructs each species with an integer offset `n`, then `NCDMBaseSpecies::ReadParameters`
  reads the CSV and indexes it with `n`. Net: dot → CSV → indexed-CSV-read.
- The constructor takes an `int species_index` whose only job is "where in the CSV to look".
- The same dual-mode dance is repeated in `DNCDMSpecies::CreateAll` and
  `NCDMInteractingSpecies::CreateAll`.
- Mixed legacy + dot input is rejected wholesale, even though the conceptual model would
  let them coexist.

This spec flips the canonical internal form to dot-syntax. Legacy CSV is normalised into
per-instance dot keys *before* any species is constructed; from that point on, only
dot-syntax exists. Constructors take an instance-name string, not an integer.

## Goals

- `NCDMSpecies`, `DNCDMSpecies`, `NCDMInteractingSpecies` constructors take an
  **instance-name string**, not an integer index.
- Parameters are read via `SpeciesInput` keyed on that name. No CSV indexing inside the
  constructor.
- `NCDMSpecies::CreateAll` is a flat: discover-instances → construct loop. Legacy CSV is
  normalised to per-instance dot keys up front.
- `DNCDMSpecies` and `NCDMInteractingSpecies` accept dot-syntax only — legacy CSV is
  rejected with a clear migration message.
- Mixed legacy + dot for standard NCDM works automatically (no "cannot mix" error).
- `ncdm_id_` (the perturbation-array slot index) becomes a separate concern from the
  constructor key — set by the orchestrator/factory after construction.

## Non-goals

- Removing `ncdm_id_` itself. Perturbations module still indexes
  `pv->index_ncdm_[n]`, `pv->q_size_ncdm[n]`, `pv->l_max_ncdm[n]`. That's a separate
  follow-up.
- Touching non-NCDM species.

## Design

### Constructor & read path

#### `NCDMBaseSpecies` — common parameters only

Single constructor:

```cpp
NCDMBaseSpecies(std::string name,
                EnergyType energy_type,
                FileContent* pfc,
                const std::string& instance_name,
                const NcdmSettings& settings);
```

Reads only the universal NCDM-family fields via `SpeciesInput`:
`m`, `T`, `deg`, `Omega`, `omega`, `ksi`, `psd_parameters`, `psd_filename`,
`use_psd_file`, `quadrature_strategy`, `momenta_bins`, `max_q`.

Each member declares its default at the class definition
(`double m_in_eV_ = 0.;`, `double T_ = 0.71611;`, `double deg_ = 1.;`,
`int quadrature_strategy_ = 3;` etc.). On a missed read, the member keeps its in-class
default — no `kStandardFields` table, no `default_value` strings. Snippet pattern:

```cpp
SpeciesInput input(pfc, instance_name);
input.read_double("m", m_in_eV_);  // unchanged on miss
input.read_double("T", T_);
// ...
```

The two existing constructor overloads (standard + decay_dr) collapse into one. No
`Variant` enum is introduced on the base; variant-specific reads happen in the subclass
constructor.

Defaults are unified across variants: the existing `m_ncdm_decay_dr` default of `1.0` is
dropped; everything uses the in-class `m_in_eV_ = 0.` default.

#### `NCDMSpecies`

```cpp
NCDMSpecies(FileContent* pfc,
            const std::string& instance_name,
            const NcdmSettings& settings,
            const background* pba,
            const BackgroundModule* bgm);
```

No `species_index`, no `suffix`. `name_` is the instance name itself (e.g. `"ncdm_1"`)
instead of `"NCDM_<index>"` — matches what the user wrote in the input file and makes
log/error messages self-explanatory. Reads no extra fields beyond what the base reads.

#### `DNCDMSpecies`

```cpp
DNCDMSpecies(FileContent* pfc,
             const std::string& instance_name,
             const NcdmSettings& settings,
             const background* pba,
             const BackgroundModule* bgm);
```

Reads its variant-specific fields directly via `SpeciesInput`:

- `Gamma` / `log10Gamma` / `lifetime` / `log10lifetime` (exactly one required;
  conversion to Mpc⁻¹ unchanged).
- `Omega` / `omega` / `Omega_ini` / `omega_ini` / `Neff_ini` — the per-species closure
  inputs that today live in `DNCDMSpecies::CreateAll`'s bulk loop at lines 236–340.

That bulk loop in `CreateAll` is deleted; each species reads its own.

#### `NCDMInteractingSpecies`

```cpp
NCDMInteractingSpecies(FileContent* pfc,
                       const std::string& instance_name,
                       const NcdmSettings& settings,
                       const background* pba,
                       const BackgroundModule* bgm);
```

Reads `G_eff` / `log10G_eff` directly via `SpeciesInput`.

### `ncdm_id_` assignment

`ncdm_id_` stays an internal field on `NCDMSpecies` / `DNCDMSpecies`. Declared with an
in-class default `int ncdm_id_ = -1;` (consistent with the rest of the field-default
strategy). New setter:

```cpp
void NCDMBaseSpecies::SetNcdmId(int id) { ncdm_id_ = id; }
```

`SpeciesBuildContext` gains a single field:

```cpp
struct SpeciesBuildContext {
  // ... existing fields ...
  int* ncdm_id_next;  // mutable counter, advanced by each NCDM-family CreateAll
};
```

Initialised to `0` by the orchestrator before the factory loop. Each NCDM-family
`CreateAll` calls `SetNcdmId((*ctx.ncdm_id_next)++)` on each species it constructs. The
order entries appear in `kAllSpeciesFactories` defines the global ncdm slot ordering,
which is already the case today (NCDM before DNCDM_DR before NCDMInt) — just made
explicit.

### `NCDMSpecies::CreateAll`

```cpp
std::vector<Named> NCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  TransmuteLegacyStandardNcdmToDot(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_standard");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    auto sp = std::make_unique<NCDMSpecies>(
        ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
```

That's the whole factory. No dot-vs-legacy branching, no `synthesise_*_flat_keys`, no
"cannot mix" error case, no `make_unique<NCDMSpecies>(pfc, n, ...)`.

### `TransmuteLegacyStandardNcdmToDot(FileContent& pfc)`

Free function in `ncdm_species.cpp`.

1. Read `N_ncdm_standard` (or fall back to `N_ncdm`). If neither set or value is 0,
   return — nothing to transmute.
2. If both `N_ncdm` and `N_ncdm_standard` are set, throw the existing "choose one"
   error.
3. Before writing anything, check whether any of the synthesised instance prefixes
   (`ncdm__1` … `ncdm__N`) already appear as a dot-instance prefix in PFC (i.e. any
   key of the form `ncdm__i.*` exists). If so, throw a clear error pointing the user
   to either rename their dot instances or drop the legacy keys.
4. For each legacy CSV key — `m_ncdm`, `T_ncdm`, `deg_ncdm`, `Omega_ncdm`, `omega_ncdm`,
   `ksi_ncdm`, `ncdm_psd_parameters`, `ncdm_psd_filenames`, `use_ncdm_psd_files`,
   `quadrature_strategy_ncdm_standard`, `N_momentum_bins_ncdm_standard`,
   `maximum_q_ncdm_standard` — if present:
   - Parse as a CSV / list. Exactly one value broadcasts to all `N` instances; `N`
     values assigns one per instance; anything else is a clear error.
   - Write per-instance dot keys: `ncdm__1.<dot-name>`, ..., `ncdm__N.<dot-name>`
     where `<dot-name>` is the `kStandardFields[*].dot` name from today's table.
5. Set `ncdm__i.type = "ncdm_standard"` for `i = 1..N`.
6. Mark legacy CSV keys as read so `unread_parameters()` stays accurate.

The `ncdm__N` (double-underscore) naming matches the convention already used at
`ncdm_species.cpp:134` for the legacy path, so no new convention is introduced.

After this helper runs, the `instances_with("type", "ncdm_standard")` query naturally
returns the **union** of legacy-synthesised and user-supplied dot instances — making
mixed legacy + dot input work without any special-case logic.

### `DNCDMSpecies::CreateAll` and `NCDMInteractingSpecies::CreateAll`

Same shape as NCDM's, **without** the transmutation step. Instead, a strict rejector:

```cpp
std::vector<Named> DNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyDecayDrKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_decay_dr");
  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    auto sp = std::make_unique<DNCDMSpecies>(
        ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
```

`RejectLegacyDecayDrKeys` detects any of the old CSV keys: `N_ncdm_decay_dr`,
`m_ncdm_decay_dr`, `T_ncdm_decay_dr`, `ksi_ncdm_decay_dr`, `deg_ncdm_decay_dr`,
`quadrature_strategy_ncdm_decay_dr`, `N_momentum_bins_ncdm_decay_dr`,
`maximum_q_ncdm_decay_dr`, `Gamma_ncdm_decay_dr`, `log10Gamma_ncdm_decay_dr`,
`lifetime_ncdm_decay_dr`, `log10lifetime_ncdm_decay_dr`, `Omega_dncdmdr`,
`omega_dncdmdr`, `Omega_ini_dncdm`, `omega_ini_dncdm`, `Neff_ini_dncdm`. If any is
present, throws with a migration message: "`<key>` is no longer supported; use
dot-syntax `<instance>.<dot-name> = ...` with `<instance>.type = ncdm_decay_dr`."

`RejectLegacyInteractingKeys` does the equivalent for `N_ncdm_interacting`,
`G_eff_ncdm_interacting`, `log10G_eff_ncdm_interacting`,
`quadrature_strategy_ncdm_interacting`, `N_momentum_bins_ncdm_interacting`,
`maximum_q_ncdm_interacting`.

## Files touched

| File | Change |
|---|---|
| `species/ncdm_base_species.h` | Single constructor taking `instance_name`. Members get in-class default initialisers. Add `SetNcdmId(int)`. |
| `species/ncdm_base_species.cpp` | Replace `ReadParameters` + `ReadDecayDrParameters` with one `SpeciesInput`-based read of the universal fields. |
| `species/ncdm_species.h/cpp` | New constructor signature. `CreateAll` becomes a flat loop. Add `TransmuteLegacyStandardNcdmToDot`. Drop `synthesise_standard_ncdm_flat_keys`. `name_` becomes the instance name. |
| `species/dncdm_species.h/cpp` | New constructor signature reads decay_dr-specific fields directly via `SpeciesInput`. `CreateAll` becomes flat loop + `RejectLegacyDecayDrKeys`. Drop the lines 236–340 bulk Omega/deg/Omega_ini loop. |
| `species/dncdm_dr_species.cpp` | Adjust to delegate to the new DNCDM CreateAll signature. |
| `species/ncdm_interacting_species.h/cpp` | New constructor signature reads `G_eff` / `log10G_eff` directly. `CreateAll` becomes flat loop + `RejectLegacyInteractingKeys`. Drop `synthesise_self_interacting_ncdm_flat_keys`. |
| `species/species_build_context.h` | Add `int* ncdm_id_next`. |
| `source/input_module.cpp` (orchestrator) | Initialise `ncdm_id_next = 0` before the factory loop. |
| `species/species_input.h/cpp` | The dot → CSV helpers (`CsvWithDefaults`, `CsvForPsdFilenames`, `SynthesiseIdenticalScalarField`, `CollectInstanceFieldValues`, `AnyInstanceFieldValue`) become unused — delete after confirming no external callers. |
| Tests / fixtures | Standard-NCDM legacy fixtures keep working via transmutation. DNCDM / NCDMInt legacy fixtures (if any) migrate to dot-syntax — error message guides them. |

## Behavioural changes visible to users

1. `m` for decay_dr species defaults to `0` rather than `1.0`. In practice always
   set explicitly; this just keeps the unified-defaults invariant.
2. DNCDM legacy syntax (`N_ncdm_decay_dr`, `Gamma_ncdm_decay_dr`, ...) is no longer
   accepted. Clear error directs users to dot-syntax.
3. `NCDMInteracting` legacy syntax (`N_ncdm_interacting`, `G_eff_ncdm_interacting`) is
   no longer accepted. Same.
4. Standard NCDM legacy + dot is no longer an error — they coexist via transmutation.
5. Species' `name_` strings change from `"NCDM_<index>"` / `"DNCDM_<index>"` to the
   user's instance name (e.g. `"ncdm_1"`). Logs / error messages reflect this.

## Documented fallback

If the implementer discovers that mutating PFC during
`TransmuteLegacyStandardNcdmToDot` is unexpectedly invasive — for example, if PFC
consumers downstream are broken by the new dot keys, or if the parser cannot mark them
as consumed cleanly — fall back to **(b)**: legacy NCDM is restricted to a single
species. `TransmuteLegacyStandardNcdmToDot` then handles only the `N_ncdm = 1` case and
throws on `N_ncdm >= 2` with a migration message. The mix scenario still works because
legacy synthesises one instance and dot adds others. Net code savings: only the
CSV-splitting branch in the helper. The user has signed off on this fallback under the
condition that a real implementation problem is discovered — not as a default.

## Out of scope / future work

- **Removing `ncdm_id_` entirely.** The user has confirmed this is the longer-term
  direction: no perturbation array should be indexed by `ncdm_id_`. That requires
  changes in `perturbations_module.cpp` — replacing `pv->index_ncdm_[n]`,
  `q_size_ncdm[n]`, `l_max_ncdm[n]` with per-species owned data — and is not part of
  this refactor. After that work, `SetNcdmId`, `ncdm_id_`, and the `ncdm_id_next`
  counter on `SpeciesBuildContext` can all be deleted.
- **Generalising the transmutation pattern to other species families** (DCDM, IDM_*,
  ...) if they grow legacy/dot dual-mode. This refactor sets the template
  (`Transmute*Legacy*ToDot`).

## Acceptance criteria

- All existing standard-NCDM legacy-input fixtures pass without changes.
- A fixture mixing legacy `N_ncdm = 1, m_ncdm = 0.06` with a dot-syntax
  `something.type = ncdm_standard` passes (was an error before).
- DNCDM and NCDMInt legacy fixtures, if any exist, are migrated to dot-syntax. New
  fixtures using legacy DNCDM/NCDMInt keys produce the migration error.
- `NCDMSpecies::CreateAll`, `DNCDMSpecies::CreateAll`,
  `NCDMInteractingSpecies::CreateAll` are each ≤ ~15 lines (the discover-and-loop
  shape).
- No constructor in the NCDM family takes an integer species index.
- `unread_parameters()` reports no false positives or negatives for either input form.
