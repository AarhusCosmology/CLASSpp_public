# Dot-notation for single-instance legacy species

**Status:** design approved, ready for implementation plan
**Date:** 2026-06-26
**Scope:** one PR, two parts (A: factory-name/type-string unification; B: the dot-notation feature)

## Goal

Let the single-instance legacy species accept the same dot-notation that
multi-instance species use, so a user can write

```ini
bary.type  = baryons
bary.Omega = 0.05
```

and have it behave exactly as `Omega_b = 0.05`. The instance name (`bary`) is
discarded — only the `type` selects the translation table, and the result is
the legacy keys, so output is byte-identical to the legacy form.

Single-instance species (one allowed instance, ever, by design decision):
**photons, baryons, cdm, ur, lambda, fluid**. ScalarField is deliberately
excluded — it is destined to become multi-instance.

This is base-code logic: *we* decided these named species are single-instance,
so the translation table belongs in the input module, not in a per-species
registry.

## Motivation / two parts

The feature needs a canonical `type` string per species. Today the user-facing
type strings (`"ncdm_standard"`, `"ncdm_greybody"`, …) are hardcoded inside each
factory's `instances_with("type", …)` call, while `kAllSpeciesFactories` carries
a *separate* CamelCase identifier (`"NCDM"`, `"NCDMGreyBody"`, …). The two are
maintained in different places, and there is no single list of "what types
exist." Adding a species means inventing two names.

So the PR has two parts:

- **Part A (output-neutral refactor):** make the `kAllSpeciesFactories` entry
  name *be* the canonical type string, sourced from a `kTypeName` constant owned
  by each species and reused inside its `CreateAll`. `kAllSpeciesFactories`
  becomes the registry of all type strings.
- **Part B (the feature):** the dot-notation → legacy-key preprocessing pass,
  keyed by the unified `kTypeName`s.

## Background facts (verified 2026-06-26)

- **Factory table:** `species/all_species.h:40-55`. Entries: `{name, create_all}`,
  `create_all` is `std::vector<Named> (*)(const SpeciesBuildContext&)`.
  Current names are CamelCase: `Photons, Baryons, CDM, UR, DCDM_DR, NCDM,
  NCDMGreyBody, DNCDM_DR, NCDMInt, IDM_DR_IDR, IDM_DRMD_IDR_DRMD, Lambda, Fluid,
  ScalarField`.
- **`entry.name` is used only twice in `ConstructSpecies`** (`source/input_module.cpp`):
  closure match `entry.name == closure_name` (lines 256, 288, where `closure_name`
  comes from `ClosureSpeciesName`, returning `"Lambda"/"Fluid"/"ScalarField"`),
  and ncdm-family detection `entry.name == "NCDM" || "DNCDM_DR" || "NCDMInt"`
  (lines 259-260). `ClosureSpeciesName` is in `species/all_species.h:58-70`.
- **`CreateAll` is invoked only through the table** (`entry.create_all(ctx)` at
  `input_module.cpp:258` and `290`); `DNCDM_DR_Species::CreateAll` additionally
  calls `DNCDMSpecies::CreateAll` internally (`species/dncdm_dr_species.cpp:12`).
  No other external callers — Part A's changes are contained.
- **Hardcoded type strings** (the multi-instance discovery literals):
  - `species/ncdm_species.cpp:205` → `"ncdm_standard"`
  - `species/greybody_ncdm_species.cpp:298` → `"ncdm_greybody"`
  - `species/dncdm_species.cpp:226` → `"ncdm_decay_dr"`
  - `species/ncdm_interacting_species.cpp:85` → `"ncdm_self_interacting"`
  - `source/input_module.cpp:594` → `"ncdm_standard"` (the `fluid_approximation`
    synthesis, a non-factory call site that also needs the unified constant)
- **Composites** (`dcdm_dr`, `idm_dr_idr`, `idm_drmd_idr_drmd`) are built from
  coupled inputs / `has_*` checks, **not** from a `type` field. They have no
  user-facing type string; their `kTypeName` is a synthetic internal id that is
  never a valid `.type` value.
- **Input phase order** (`source/input_module.cpp:208-210`):
  `input_read_precisions()` → `ReadContext()` (gauge/h/closure) →
  `ConstructSpecies()`. The translation pass must run **before** all three so
  translated keys feed precision, closure detection, and species reads.
- **Existing preprocessing precedent:** `SynthesiseIdenticalScalarField`
  (`species/species_input.{h,cpp}`, called at `input_module.cpp:594-599`). Despite
  the name it has nothing to do with the ScalarField species — it collapses a
  per-instance dot-field (which must be identical across instances) into one
  legacy key, with conflict detection. Part B reuses this exact validation shape.
- **`FileContent` semantics** (`include/parser.h`):
  - `set(name, value)` *silently overwrites* and marks the key unread (line 36-37).
  - `instances_with(field, value)` returns instance names `N` where `N.field ==
    value`, in insertion order, marking nothing read (lines 84-86).
  - `was_read` / `get<T>` mark reads; unread keys trip the
    "input line not recognized" warning at `input_module.cpp:575`.
- **Legacy keys per single-instance species** (translation targets):
  - cdm: `Omega_cdm` (`input_module.cpp:395`)
  - photons: `T_cmb`, `Omega_g` (`input_module.cpp:695-696`)
  - baryons: `Omega_b` (`input_module.cpp:725`)
  - lambda: `Omega_Lambda` (`input_module.cpp:788`)
  - fluid: `Omega_fld` (`input_module.cpp:789`) plus `w0_fld`, `wa_fld`,
    `cs2_fld`, `Omega_EDE`, `c_gamma_over_c_fld`, `use_ppf`,
    `fluid_equation_of_state` (`species/fluid.cpp:435-482`)
  - ur: `N_ur` (deprecated alias `N_eff`), `Omega_ur`, `omega_ur`
    (`species/ultra_relativistic.cpp:422-475`); only one of the density forms may
    be given (existing UR validation).

## Part A — unify factory name ↔ canonical type string

**Output-neutral.** Pure internal rename; outputs must be `cmp`-identical.

1. Add `static constexpr std::string_view kTypeName` to every species class.
   Values:
   - `photons`, `baryons`, `cdm`, `ur`, `lambda`, `fluid`, `scalar_field`
   - `ncdm_standard`, `ncdm_greybody`, `ncdm_decay_dr`, `ncdm_self_interacting`
   - composites (synthetic): `dcdm_dr`, `idm_dr_idr`, `idm_drmd_idr_drmd`
2. `kAllSpeciesFactories` rows reference `Species::kTypeName` instead of string
   literals.
3. Replace the hardcoded `instances_with("type", "…")` literals (4 factory call
   sites + `input_module.cpp:594`) with the corresponding `kTypeName`.
4. Update `ClosureSpeciesName` to return `lambda` / `fluid` / `scalar_field`, and
   the `is_ncdm_family` comparison to `ncdm_standard` / `ncdm_decay_dr` /
   `ncdm_self_interacting`. (Note: `ncdm_greybody` is intentionally *not* in the
   ncdm-family counter today — preserve that.)

**Decision:** source of truth is the per-species `kTypeName` constant (the species
owns its identity — consistent with "all species are equal"), not a
signature change to `CreateAll`. This avoids churning 14 function signatures.

**Critical:** `kTypeName` is the factory/registry/`type`-input identifier. It is
**distinct from** the species' runtime `name_` (the display/output name set in
constructors, e.g. `BaseSpecies("UR", …)` at `species/ultra_relativistic.cpp:14`).
Part A must **not** touch `name_` — output column names derive from `name_`, and
changing it would break output-neutrality. Only the factory-table entry name and
the `instances_with` literals move to `kTypeName`.

## Part B — dot-notation → legacy translation

A new preprocessing pass (e.g. `TranslateSingleInstanceDotSyntax()`) invoked
once at the very start of input processing, before `input_read_precisions()`.

### Translation table (module-owned, keyed by `kTypeName`)

| type | dot-field → legacy key |
|---|---|
| `photons` | `Omega`→`Omega_g`, `T_cmb`→`T_cmb` |
| `baryons` | `Omega`→`Omega_b` |
| `cdm` | `Omega`→`Omega_cdm` |
| `ur` | `N`→`N_ur`, `Omega`→`Omega_ur`, `omega`→`omega_ur` |
| `lambda` | `Omega`→`Omega_Lambda` |
| `fluid` | `Omega`→`Omega_fld`, `w0`→`w0_fld`, `wa`→`wa_fld`, `cs2`→`cs2_fld`, `Omega_EDE`→`Omega_EDE`, `c_gamma_over_c`→`c_gamma_over_c_fld`, `use_ppf`→`use_ppf`, `equation_of_state`→`fluid_equation_of_state` |

ScalarField is excluded (not in the table) despite having a `kTypeName`.

### Naming convention (one rule)

The dot-field is the legacy key with the **species-identifying prefix/suffix
stripped**:
- `Omega_b`→`Omega`, `Omega_cdm`→`Omega`, `Omega_g`→`Omega`,
  `Omega_Lambda`→`Omega`, `Omega_fld`→`Omega`, `Omega_ur`→`Omega`
- `w0_fld`→`w0`, `wa_fld`→`wa`, `cs2_fld`→`cs2`,
  `c_gamma_over_c_fld`→`c_gamma_over_c`
- `N_ur`→`N`, `omega_ur`→`omega`
- `fluid_equation_of_state`→`equation_of_state` (prefix stripped)

Keys with no species marker stay verbatim: `T_cmb`, `use_ppf`, `Omega_EDE`.

The deprecated UR alias `N_eff` is **not** exposed in dot-form (use `ur.N`); the
legacy `N_eff` key still works on its own.

### Algorithm (per single-instance table entry)

1. `instances = file_content_.instances_with("type", type)`.
2. `instances.empty()` → skip.
3. `instances.size() > 1` → **throw** `std::invalid_argument`: single-instance
   type given more than once, listing the offending instance names.
   *Required for correctness:* `set()` overwrites silently, so without this guard
   two `baryons` instances would clobber `Omega_b` and run with one value, no
   error.
4. Otherwise take the single instance name (the name itself is irrelevant). For
   each `dot_field → legacy_key` in the table where `<name>.<dot_field>` is
   present:
   - **Conflict check:** if `legacy_key` is already present with a *different*
     value → throw (same message style as `SynthesiseIdenticalScalarField`);
     identical value is fine.
   - `pfc->set(legacy_key, value)`.
   - Mark `<name>.<dot_field>` read (via the `get` that fetched it).
   - Mark `<name>.type` read.
5. Unknown `<name>.<field>` not in the table is left unread → trips the existing
   "input line not recognized" warning (`input_module.cpp:575`), so typos still
   surface. No extra code.

### Why output is unchanged

The translated species is built by its normal factory from legacy keys, named by
its factory; the dot-instance name never reaches the species. A dot-form `.ini`
and the equivalent legacy-form `.ini` produce identical `FileContent` going into
`ReadContext`/`ConstructSpecies`.

## Verification model

**Bar: byte-identical output** (this is a translation + rename, not an
arithmetic change).

- **Part A:** run the existing scenarios on master vs branch; `cmp`-identical.
  Confirm every type-discovered species (ncdm family) is still found.
- **Part B (per species):** author a dot-form `.ini` for each of the six
  single-instance species and assert `cmp`-identical output to the legacy-form
  `.ini`. At minimum cover: baryons, cdm, photons (`T_cmb` and `Omega`), ur
  (`N`, `Omega`, `omega`), lambda, fluid (full field set).
- **Negative tests:** two instances of the same single-instance type → error;
  `Omega_b` and `bary.Omega` set to different values → error; identical values →
  OK.
- `ctest --test-dir build --output-on-failure` green.
- Python smoke: `pip install .` (regenerates `cclassy.pxd`) + `TEST_LEVEL=1`
  scenario pytest.

## Out of scope

- ScalarField dot-notation (reserved for its multi-instance future).
- Making composite `type` strings user-facing (they remain synthetic ids).
- Reducing the broader `is_ncdm_family` / string-name-branch smell (orthogonal;
  Part A only renames the existing comparisons, it does not remove them).

## Decided questions

- **A vs B vs C for table location:** A — centralized in the input module.
- **Rename vs prefix-strip:** clean names via the strip rule (rename), not raw
  legacy keys.
- **UR included:** yes. UR's N-field is `ur.N` (strip rule), not `ur.N_ur`.
- **Fluid coverage:** full field set.
- **Bundling:** Parts A and B ship as one PR.
