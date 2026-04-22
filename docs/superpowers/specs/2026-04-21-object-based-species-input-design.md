# Object-based species input (issue #246) — Design

**Issue:** [AarhusCosmology/CLASSpp#246](https://github.com/AarhusCosmology/CLASSpp/issues/246)

## Goal

Introduce a dot-syntax input format for cosmological species
(`nu1.type = ncdm_standard; nu1.m = 0.06`) using the existing per-class
`CreateAll()` factories. NCDM is the first and only mover in this iteration:
dot-syntax is fully wired for the three existing NCDM variants, and legacy flat
syntax (`N_ncdm`, `m_ncdm`, `T_ncdm`, …) continues to work through the existing
legacy readers. Dot-syntax instances are discovered first and then synthesised
back into the flat legacy keys each factory already understands. Other existing
species (CDM, Lambda, UR, Fluid, IDM_DR_IDR, IDM_DRMD_IDR_DRMD, ScalarField,
DCDM_DR) are untouched — their construction path is unchanged and they are not
migrated to dot-syntax in this PR.

## Scope

**In scope.** `FileContent::instances_with()`, `SpeciesInput`, per-class
dot-syntax discovery inside the three NCDM `CreateAll()` factories
(`NCDMSpecies`, `DNCDMSpecies`/`DNCDM_DR_Species`, `NCDMInteractingSpecies`),
dot-syntax-to-legacy synthesis for the flat keys these classes already consume,
named `{key, species}` return values for insertion into `SpeciesCollection`,
tests and documentation for the new syntax.

**Out of scope.**
- Migrating non-NCDM species to dot-syntax. The framework is ready for them,
  but this PR does not touch their construction.
- Removing `pba->N_ncdm` or other legacy scalar background state. Where a
  counter or guard becomes trivially unused as a side effect of this refactor,
  delete it; otherwise leave it.
- Any input-format change besides adding the dot-syntax path.

## User-facing syntax

### New (dot-syntax)

```ini
nu1.type  = ncdm_standard
nu1.m     = 0.06
nu1.T     = 0.71611
nu1.deg   = 1.0

nu2.type  = ncdm_self_interacting
nu2.m     = 0.07
nu2.log10Geff = -5.5
```

Instance names (`nu1`, `nu2`, …) are user-chosen, must match
`[A-Za-z_][A-Za-z0-9_]*`, and become the `SpeciesCollection::Entry.key` for
that species. They must be unique; duplicates raise `std::invalid_argument` at
insert time.

### Legacy (still supported for standard NCDM)

```ini
N_ncdm  = 1
m_ncdm  = 0.06
T_ncdm  = 0.71611
```

Internally, the legacy path continues to feed the existing flat-key readers.
For standard NCDM the inserted species keys still default to `ncdm__1`,
`ncdm__2`, …, while dot-syntax keeps the user-provided instance names as the
`SpeciesCollection` keys.

### Mixing is rejected

The factories reject mixing dot-syntax `*.type = ncdm_*` discovery with the
legacy NCDM count keys (`N_ncdm`, `N_ncdm_standard`, `N_ncdm_decay_dr`,
`N_ncdm_interacting`). Users pick one entry path per NCDM family in a given
input file.

## Canonical `.type` names

| Class                       | `.type` string           |
|-----------------------------|--------------------------|
| `NCDMSpecies`               | `ncdm_standard`          |
| `DNCDM_DR_Species`          | `ncdm_decay_dr`          |
| `NCDMInteractingSpecies`    | `ncdm_self_interacting`  |

These strings also become the value of `BaseSpecies::name_` for instances of
these classes (reinterpretation of `name_` as a class-level identifier, per
the issue). Instance identity moves to `SpeciesCollection::Entry.key` —
already the case since PR #258.

## Architecture

```text
.ini file
  │
  ▼
FileContent
  │
  ├─ input_read_precisions()
  │    └─ optional early synthesis of scalar precision knobs that dot-syntax
  │       shares with legacy input (currently `ncdm_fluid_approximation`)
  │
  ├─ input_read_parameters()
  │    └─ temporary NCDM CreateAll() passes compute counts/Omega totals
  │
  └─ ConstructSpecies()
       ├─ NCDMSpecies::CreateAll()
       ├─ DNCDM_DR_Species::CreateAll()
       └─ NCDMInteractingSpecies::CreateAll()
             └─ each factory finds `<instance>.type`, synthesises the legacy
                flat keys it needs, and returns named `{key, species}` entries
```

## Components

### `species/species_input.{h,cpp}` (new)

A thin per-instance wrapper over `FileContent`, plus small helper routines used
by the NCDM factories when synthesising legacy CSV keys from dot-syntax input.

```cpp
class SpeciesInput {
 public:
  SpeciesInput(FileContent* pfc, std::string instance_name);

  const std::string& instance_name() const;

  // Typed reads. Each returns true if the field was present (and sets out);
  // false otherwise. Present keys are marked as read in the underlying
  // FileContent so unread_parameters() stays accurate.
  bool read_double(const std::string& field, double& out);
  bool read_int(const std::string& field, int& out);
  bool read_string(const std::string& field, std::string& out);
  bool read_list_of_doubles(const std::string& field, std::vector<double>& out);

  // Throws std::invalid_argument with a message naming the instance and
  // field if the key is missing.
  double required_double(const std::string& field);
  int    required_int(const std::string& field);
  std::string required_string(const std::string& field);

 private:
  FileContent* pfc_;
  std::string  instance_name_;
};
```

Implementation is straightforward: build `instance_name_ + "." + field` and
delegate to the existing `FileContent` methods.

### `FileContent::instances_with()` (new)

This helper scans the parsed parameter table for keys of the form
`<instance>.<field>` whose value matches a requested string and returns the
instance names in insertion order. The NCDM factories use it to find
`*.type = ncdm_standard`, `*.type = ncdm_decay_dr`, and
`*.type = ncdm_self_interacting` without consuming read state.

### Per-class `CreateAll()` synthesis

Each NCDM-family factory now follows the same pattern:

1. Find dot-syntax instances for its `.type`.
2. Reject incompatible legacy count keys when both worlds are present.
3. Read per-instance dot fields with `SpeciesInput`.
4. Synthesis the legacy flat keys that the existing constructor path already
   understands, including default-padding and compressed PSD filename lists
   where required by the legacy readers.
5. Construct species through the unchanged legacy constructor and return
   `{key, species}` pairs for insertion into `SpeciesCollection`.

`InputModule` uses the same factories twice: once temporarily during
`input_read_parameters()` to determine `pba->N_ncdm`, `pba->N_decay_dr`, and
`pba->Omega0_ncdm_tot`, and once in `ConstructSpecies()` for the real species
objects. The other branches (`has_cdm`, `has_lambda`, …) are untouched.

### `BaseSpecies::name_` reinterpretation

Today each NCDM instance passes a per-instance name (e.g. `"NCDM_0"`) into the
`BaseSpecies` base constructor. Under this design each class passes its
canonical registry type instead (`"ncdm_standard"`, etc.). Any downstream use
of `sp->name()` that wanted a *per-instance* identifier must switch to the
`SpeciesCollection::Entry.key` — the PR #258 audit already labelled these
uses, so the surface is small and known.

## Error handling

All framework errors are thrown as `std::invalid_argument` (consistent with
the existing `SpeciesCollection` error contract):

- Unknown `.type`:
  `"species 'nu1': unknown type 'foo'. Known types: ncdm_standard, ncdm_decay_dr, ncdm_self_interacting"`
- Legacy + dot NCDM mixed:
  `"cannot mix legacy N_ncdm/m_ncdm keys with dot-syntax nuN.type = ncdm_*; use one or the other"`
- Required field missing:
  `"species 'nu1' (ncdm_standard): missing required field 'm'"`
- Duplicate instance name: thrown by `SpeciesCollection::insert`.
- Typos on optional fields: surface via the existing `unread_parameters()`
  warning at the end of `InputModule` — no new mechanism needed, because
  `SpeciesInput` correctly marks *read* keys as read, leaving unread
  (misspelled) keys visible to the diagnostic.

## Testing

**Unit tests** (whichever harness the existing `test/` directory uses — plan
task will check and match):

- `SpeciesInput`: key-prefixing, read/missing behaviour, `required_*` throws
  with the right message.
- `FileContent::instances_with()`: discovers the expected instance names without
  marking `*.type` as read.
- Per-class `CreateAll()` synthesis: partial dot fields pick up the same legacy
  defaults, compressed PSD filename lists follow the `use_psd_file` flags, and
  mixed interacting `G_eff`/`log10G_eff` representations are rejected early.

**Integration tests** at `InputModule` level: feed a minimal .ini with
dot-syntax only and a matching legacy-only .ini; assert both produce the same
`SpeciesCollection` contents (instance count, types, per-instance parameter
values via species accessors).

**Byte-equal regression** against master for three inputs:
- `explanatory.ini` (legacy path, `N_ncdm=1`)
- `base_2015_plikHM_TT_lowTEB_lensing.ini` (the standard regression baseline)
- a new `test/dotsyntax_ncdm.ini` with `nu1.type = ncdm_standard` parameterised
  to reproduce `explanatory.ini` physics; asserted byte-equal against the
  legacy run.

## Build systems

Every new source file must be registered in all three build systems. In this PR
that applies to `species/species_input.{h,cpp}`:

- `Makefile` — add `.opp` entries to `SPECIES_OPP` / `SOURCE_OPP`.
- `setup.py` — add `.cpp` entries to the species / source source_files lists.
- `CLASS.xcodeproj/project.pbxproj` — add PBXBuildFile + PBXFileReference +
  group children + Sources build phase entries (the user validates Xcode
  manually; plan records the expected four-entry-per-file pattern established
  by PR #258).

## Migration notes for future species

A new species type `FooSpecies` opts into dot-syntax by:

1. Teaching its existing factory (or adding a `CreateAll()`-style factory) to
   discover `foo.type` instances with `FileContent::instances_with()`.
2. Reading per-instance fields through `SpeciesInput`.
3. Synthesising whatever legacy keys its constructor path already consumes.
4. Returning named `{key, species}` entries for insertion into
   `SpeciesCollection`.
5. Adding any new `.cpp`/`.h` files to all three build systems.

No other framework code needs editing. Users then write
`my_foo.type = foo; my_foo.param_a = ...` in the .ini.
