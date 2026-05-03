# NCDM-family key-name-driven input — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refactor `NCDMSpecies`, `DNCDMSpecies`, and `NCDMInteractingSpecies` so their constructors take an instance-name string instead of an integer index, reading parameters via `SpeciesInput`. Standard-NCDM legacy CSV is normalised into per-instance dot keys before any species is constructed; DNCDM and NCDMInteracting drop legacy CSV.

**Architecture:** Parallel-pathway migration. We add the new constructor (instance-name based) alongside the old (int-index based) for each class, switch the factory to the new path, validate with smoke-runs against `.ini` fixtures, then delete the old path. The build stays green throughout. After all three families are migrated, the dot→CSV synthesis helpers (`synthesise_*_flat_keys`, `CsvWithDefaults`, `CsvForPsdFilenames`, `SynthesiseIdenticalScalarField`) are deleted.

**Tech Stack:** C++17, `std::map`-backed `FileContent`, the existing `SpeciesInput` per-instance reader at `species/species_input.h`, the unified factory loop (`kAllSpeciesFactories`) at `species/all_species.h`. Tests are integration-only — `.ini` fixtures fed to the compiled `./class` binary, plus `python -m pytest` for the full suite.

**Spec:** `docs/superpowers/specs/2026-04-30-ncdm-key-name-input-design.md`

## File map

### Files modified

- `species/species_build_context.h` — add `int* ncdm_id_next` field
- `source/input_module.cpp` — initialise `ncdm_id_next = 0` before factory loop
- `species/ncdm_base_species.h` — add new constructor, `SetNcdmId` setter; in-class default for `ncdm_id_`
- `species/ncdm_base_species.cpp` — new constructor body using `SpeciesInput`; delete old
- `species/ncdm_species.h` — new constructor signature; delete old
- `species/ncdm_species.cpp` — `TransmuteLegacyStandardNcdmToDot`; new flat `CreateAll`; delete old `synthesise_standard_ncdm_flat_keys`
- `species/dncdm_species.h` — new constructor signature; delete old
- `species/dncdm_species.cpp` — new constructor reads decay_dr fields directly via `SpeciesInput`; `RejectLegacyDecayDrKeys`; new flat `CreateAll`; delete old `synthesise_decay_dr_ncdm_flat_keys` and bulk Omega/deg/Omega_ini loop
- `species/dncdm_dr_species.cpp` — adjust delegation to new DNCDMSpecies CreateAll
- `species/ncdm_interacting_species.h` — new constructor signature; delete old
- `species/ncdm_interacting_species.cpp` — new constructor reads interacting fields directly; `RejectLegacyInteractingKeys`; new flat `CreateAll`; delete old `synthesise_self_interacting_ncdm_flat_keys`
- `species/species_input.h` — delete unused helpers (`CsvWithDefaults`, `CsvForPsdFilenames`, `SynthesiseIdenticalScalarField`, `CollectInstanceFieldValues`, `AnyInstanceFieldValue`) after final cleanup
- `species/species_input.cpp` — delete corresponding implementations

### Files added

- `test/legacy_ncdm_two_species.ini` — pure-legacy multi-species fixture (regression coverage)
- `test/dotsyntax_dncdm.ini` — DNCDM dot-only fixture
- `test/dotsyntax_ncdm_interacting.ini` — NCDMInteracting dot-only fixture

## Conventions used in steps

- **Build command:** `make -j class` (from repo root). Expected: builds `class` binary at the repo root.
- **Smoke run:** `./class <path/to/fixture.ini>`. Expected: exits 0; produces `output/<root>...` files; no diagnostic line containing "ERROR" or "Error" or "Aborted".
- **Full test suite:** `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py`. Slow; only run at end of plan and at major checkpoints.
- **Commits:** small and frequent, message format `<area>: <imperative>` matching existing style (e.g. `species: introduce ncdm_id_next counter on SpeciesBuildContext`).

---

## Phase 1 — Infrastructure

### Task 1: Add `ncdm_id_next` counter to `SpeciesBuildContext`

**Files:**
- Modify: `species/species_build_context.h`
- Modify: `source/input_module.cpp` (the call site that constructs and passes `SpeciesBuildContext`)

- [ ] **Step 1: Locate the orchestrator call site**

Run: `grep -n "SpeciesBuildContext" source/input_module.cpp`
Expected: one or more lines that construct a `SpeciesBuildContext` literal and pass it to `kAllSpeciesFactories[i].create_all(ctx)`. Read the surrounding ~30 lines to confirm where `ctx` is built and where the loop is.

- [ ] **Step 2: Add the counter field**

In `species/species_build_context.h`, modify the `SpeciesBuildContext` struct:

```cpp
struct SpeciesBuildContext {
  FileContent* pfc;
  const background* pba;
  const precision* ppr;
  const NcdmSettings* ncdm_settings;  // non-null
  const BackgroundModule* bgm;        // nullptr at species-construction time
  int* ncdm_id_next;                  // mutable counter, advanced by NCDM-family CreateAll
};
```

- [ ] **Step 3: Initialise the counter in `input_module.cpp`**

At the call site identified in Step 1, declare a local `int ncdm_id_next = 0;` immediately *before* the `SpeciesBuildContext` is built, and add `&ncdm_id_next` to the `ctx` initialiser. For example:

```cpp
int ncdm_id_next = 0;
SpeciesBuildContext ctx{pfc, pba, ppr, &ncdm_settings, /*bgm=*/nullptr, &ncdm_id_next};
```

(Match the actual field order and any other initialisation patterns already present at the call site.)

- [ ] **Step 4: Build**

Run: `make -j class`
Expected: clean build. No new warnings about `ncdm_id_next` (it is unused at this point, but field-init is fine).

- [ ] **Step 5: Smoke-run a baseline NCDM fixture**

Run: `./class test/dotsyntax_ncdm.ini`
Expected: exits 0; no error output. (The counter exists but is not consumed yet.)

- [ ] **Step 6: Commit**

```bash
git add species/species_build_context.h source/input_module.cpp
git commit -m "species: introduce ncdm_id_next counter on SpeciesBuildContext"
```

### Task 2: Add `SetNcdmId` setter and in-class default for `ncdm_id_`

**Files:**
- Modify: `species/ncdm_species.h` (move `ncdm_id_` to in-class init `= -1`)
- Modify: `species/dncdm_species.h` (same)
- Modify: `species/ncdm_base_species.h` (declare `SetNcdmId` on the base)
- Modify: `species/ncdm_base_species.cpp` (or inline the setter in the header)

- [ ] **Step 1: Move `ncdm_id_` to in-class default in `NCDMSpecies`**

In `species/ncdm_species.h`, line ~97, change:
```cpp
int ncdm_id_;  // species index (0-based), used for pv->index_ncdm_ etc.
```
to:
```cpp
int ncdm_id_ = -1;  // perturbation-array slot index; assigned by CreateAll via SetNcdmId
```

- [ ] **Step 2: Move `ncdm_id_` to in-class default in `DNCDMSpecies`**

In `species/dncdm_species.h`, find the `ncdm_id_` member (`grep -n "int ncdm_id_" species/dncdm_species.h`) and change its declaration to:
```cpp
int ncdm_id_ = -1;
```

- [ ] **Step 3: Add `SetNcdmId` to `NCDMBaseSpecies`**

The base class doesn't currently own `ncdm_id_` (it lives on the subclasses). Add a virtual setter on the base in `species/ncdm_base_species.h`:

```cpp
// Public section, near the other public methods:
virtual void SetNcdmId(int id) = 0;
```

In each subclass header (`ncdm_species.h`, `dncdm_species.h`), override it inline:

```cpp
void SetNcdmId(int id) override { ncdm_id_ = id; }
```

(`NCDMInteractingSpecies` inherits from `NCDMSpecies`, so it inherits the override.)

- [ ] **Step 4: Update existing constructors to no longer set `ncdm_id_(species_index)` in initialiser list**

In `species/ncdm_species.cpp` line ~24:
```cpp
ncdm_id_(species_index), pba_(pba) {
```
becomes:
```cpp
pba_(pba) {
```
plus the body adds:
```cpp
ncdm_id_ = species_index;  // legacy path; will be replaced with SetNcdmId in Task 7
```

In `species/dncdm_species.cpp` line ~115, do the same: remove `ncdm_id_(ncdm_index)` from the initialiser list and add `ncdm_id_ = ncdm_index;` as the first body line.

- [ ] **Step 5: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 6: Smoke-run a baseline NCDM fixture**

Run: `./class test/dotsyntax_ncdm.ini`
Expected: exits 0; same output behaviour as before.

- [ ] **Step 7: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp species/dncdm_species.h species/dncdm_species.cpp species/ncdm_base_species.h
git commit -m "species: add SetNcdmId setter and in-class default for ncdm_id_"
```

---

## Phase 2 — Standard NCDM new pathway (alongside old)

### Task 3: Add a new `NCDMBaseSpecies` constructor that reads via `SpeciesInput`

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_base_species.cpp`

We add a *third* constructor overload that takes an instance name and reads via `SpeciesInput`. The two existing CSV-based overloads stay; we'll delete them in a later task.

- [ ] **Step 1: Declare the new constructor and a private helper**

In `species/ncdm_base_species.h`, in the `protected:` section near the existing constructors, add:

```cpp
// New input path: reads parameters per-instance via SpeciesInput.
// Used by NCDMSpecies / DNCDMSpecies / NCDMInteractingSpecies once their
// CreateAll factories are migrated to the dot-syntax-only internal form.
NCDMBaseSpecies(std::string name,
                EnergyType energy_type,
                FileContent* pfc,
                const std::string& instance_name,
                const NcdmSettings& settings);
```

In the `private:` section, add:
```cpp
void ReadParametersByInstance(FileContent* pfc,
                              const std::string& instance_name,
                              const NcdmSettings& settings);
```

- [ ] **Step 2: Implement the new constructor**

In `species/ncdm_base_species.cpp`, after the two existing constructor definitions, add:

```cpp
NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 const std::string& instance_name,
                                 const NcdmSettings& settings)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadParametersByInstance(pfc, instance_name, settings);
}
```

- [ ] **Step 3: Implement `ReadParametersByInstance`**

In `species/ncdm_base_species.cpp`, append:

```cpp
void NCDMBaseSpecies::ReadParametersByInstance(FileContent* pfc,
                                               const std::string& instance_name,
                                               const NcdmSettings& settings) {
  SpeciesInput input(pfc, instance_name);

  // Common scalar fields. read_double / read_int leave the destination unchanged
  // on miss, so each member keeps its in-class default. (Defaults that were
  // previously fixed-up here, like T_ncdm_default = 0.71611, must already be set
  // as in-class member initialisers — see Task 4.)
  input.read_double("m", m_in_eV_);
  input.read_double("T", T_);
  input.read_double("deg", deg_);
  input.read_double("Omega", Omega0_);
  input.read_double("omega", omega0_);
  input.read_double("ksi", ksi_);
  input.read_int("quadrature_strategy", quadrature_strategy_);
  input.read_int("momenta_bins", input_q_size_);
  input.read_double("max_q", qmax_);

  // psd_parameters: variable-length list (single instance, comma-separated values).
  std::vector<double> psd_params;
  if (input.read_list_of_doubles("psd_parameters", psd_params)) {
    psd_parameters_ = std::move(psd_params);
  }

  // PSD file: a single-instance flag and an optional filename.
  int use_psd_file = 0;
  input.read_int("use_psd_file", use_psd_file);
  got_file_ = (use_psd_file != 0);
  if (got_file_) {
    std::string filename;
    if (!input.read_string("psd_filename", filename)) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': use_psd_file=1 but psd_filename is missing");
    }
    psd_file_ = std::move(filename);
  }

  // Resolve omega/Omega conflict (matches existing semantics).
  if (omega0_ != 0.0) {
    if (Omega0_ != 0.0) {
      throw std::invalid_argument("species '" + instance_name +
                                  "': both Omega and omega specified — choose one");
    }
    Omega0_ = omega0_ / settings.h / settings.h;
  }
  else {
    omega0_ = Omega0_ * settings.h * settings.h;
  }

  // Ultra-relativistic default: if both Omega and m are zero, give a tiny mass.
  if ((Omega0_ == 0.0) && (m_in_eV_ == 0.0)) {
    m_in_eV_ = 1.e-5;
  }

  // The rest of the construction matches the legacy path.
  InitQuadrature(settings);

  if (m_in_eV_ != 0.0) {
    M_ = m_in_eV_ / _k_B_ * _eV_ / T_ / T_cmb_;
    double rho_ncdm;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
  }
  else {
    M_ = 0.;
  }
}
```

- [ ] **Step 4: Build**

Run: `make -j class`
Expected: clean build. The new constructor is unused, so no behaviour changes.

- [ ] **Step 5: Smoke-run baseline**

Run: `./class test/dotsyntax_ncdm.ini` and `./class explanatory.ini`
Expected: both exit 0.

- [ ] **Step 6: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp
git commit -m "species: add SpeciesInput-based NCDMBaseSpecies constructor (alongside legacy)"
```

### Task 4: Move per-member defaults to in-class initialisers

**Files:**
- Modify: `species/ncdm_base_species.h`

The new constructor in Task 3 only writes a member when the field is found. Members that today rely on `T_ncdm_default = 0.71611` etc. inside `ReadParameters` need their default expressed at the class level. (Several already do — confirm and fill gaps.)

- [ ] **Step 1: Audit current in-class defaults**

Run: `grep -n "= [0-9.]" species/ncdm_base_species.h | head -20`
Look at each protected member from `m_in_eV_` to `omega0_`. Confirm each has an in-class default that matches the legacy default-value used by `ReadParameters`.

Today:
| Member | Current in-class default | Required default |
|---|---|---|
| `m_in_eV_` | `0.` | `0.` ✓ |
| `M_` | `0.` | `0.` ✓ |
| `deg_` | `1.` | `1.` ✓ |
| `T_` | `0.71611` | `0.71611` ✓ |
| `ksi_` | `0.` | `0.` ✓ |
| `quadrature_strategy_` | `3` (private) | `0` ⚠️ — change to `0` to match legacy default |
| `input_q_size_` | `-1` (private) | `5` ⚠️ — change to `5` |
| `qmax_` | `15.` | `15.` ✓ |
| `Omega0_` | `0.` (private) | `0.` ✓ |
| `omega0_` | `0.` (private) | `0.` ✓ |
| `got_file_` | `false` | `_FALSE_`/false ✓ |
| `psd_file_` | `""` | `""` ✓ |
| `psd_parameters_` | `{}` | `{}` ✓ |

- [ ] **Step 2: Adjust quadrature defaults**

In `species/ncdm_base_species.h`, change:
```cpp
int quadrature_strategy_ = 3;
int input_q_size_        = -1;
```
to:
```cpp
int quadrature_strategy_ = 0;
int input_q_size_        = 5;
```

(These match the defaults the legacy `ReadParameters` was supplying via `read_list_of_ints("quadrature_strategy_ncdm", ..., 0)` and `read_list_of_ints("N_momentum_bins_ncdm", ..., 5)` at lines 247 and 251.)

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 4: Run baseline fixtures and full pytest run for safety**

```bash
./class test/dotsyntax_ncdm.ini
./class explanatory.ini
```
Expected: both exit 0. The legacy `ReadParameters` overwrites these members anyway (it always writes the field), so legacy fixtures behaviour is unchanged.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_base_species.h
git commit -m "species: align in-class quadrature defaults with legacy fallbacks"
```

### Task 5: Add a new `NCDMSpecies` constructor that takes an instance name

**Files:**
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp`

- [ ] **Step 1: Declare the new constructor**

In `species/ncdm_species.h`, in the `public:` section just below the existing constructor declaration, add:

```cpp
// New input path: parameters are read from PFC under <instance_name>.<field>.
NCDMSpecies(FileContent* pfc,
            const std::string& instance_name,
            const NcdmSettings& settings,
            const background* pba,
            const BackgroundModule* bgm);
```

- [ ] **Step 2: Implement the new constructor**

In `species/ncdm_species.cpp`, after the existing constructor (line ~26), add:

```cpp
NCDMSpecies::NCDMSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name,
                      EnergyType::Other,
                      pfc,
                      instance_name,
                      settings),
      pba_(pba) {
  bgm_ = bgm;

  // Standard-NCDM closure: NCDMBaseSpecies::ReadParametersByInstance only
  // computes M_ from m_in_eV_ (matching legacy ReadDecayDrParameters). For
  // standard NCDM we additionally:
  //   - if m != 0 and Omega0 == 0: derive Omega0_ from rho_ncdm
  //   - if m != 0 and Omega0 != 0: rescale factor_ and deg_ by fnu_factor
  //   - if m == 0:                 derive M_ via MFromOmega and back-compute m_in_eV_
  // (Today's ReadParameters does this in lines 446–466 of ncdm_base_species.cpp.)
  const double H0 = settings.h * 1.e5 / _c_;
  if (m_in_eV_ != 0.0) {
    double rho_ncdm = 0.;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
    if (Omega0_ == 0.0) {
      SetOmega0(rho_ncdm / H0 / H0, settings.h);
    }
    else {
      const double fnu_factor = H0 * H0 * Omega0_ / rho_ncdm;
      factor_ *= fnu_factor;
      deg_    *= fnu_factor;
    }
  }
  else {
    M_       = MFromOmega(H0, Omega0_, settings.tol_M_ncdm);
    m_in_eV_ = _k_B_ / _eV_ * T_ * M_ * T_cmb_;
  }

  // ncdm_id_ stays at its in-class default (-1) until SetNcdmId is called.
}
```

(The `name_` of the species is the instance name itself — log messages and column headers that today read `NCDM_<index>` will read e.g. `nu1` after migration. This is intended.)

`SetOmega0`, `MFromOmega`, `factor_`, `deg_`, `M_`, `Omega0_`, `m_in_eV_`, `T_`, `T_cmb_` are all already accessible to a subclass: `SetOmega0` and `MFromOmega` are `protected`; the data members are also `protected` on the base. Confirm via `grep -n "SetOmega0\|MFromOmega\|protected:" species/ncdm_base_species.h` if uncertain.

`_c_`, `_k_B_`, `_eV_` come from `common.h` via the include chain. `ComputeMomenta` is a public member of `NCDMBaseSpecies`.

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build. New constructor is unused.

- [ ] **Step 4: Smoke-run baseline**

Run: `./class test/dotsyntax_ncdm.ini`
Expected: exits 0. Behaviour unchanged.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp
git commit -m "species: add instance-name NCDMSpecies constructor (alongside legacy)"
```

### Task 6: Implement `TransmuteLegacyStandardNcdmToDot`

**Files:**
- Modify: `species/ncdm_species.cpp`

This helper rewrites legacy `N_ncdm` / `N_ncdm_standard` + per-field CSVs into per-instance dot keys (`ncdm__1.m`, `ncdm__1.T`, ..., `ncdm__N.<field>` plus `ncdm__i.type = "ncdm_standard"`), then marks the legacy keys as read.

- [ ] **Step 1: Add the helper inside the anonymous namespace**

In `species/ncdm_species.cpp`, inside the existing anonymous namespace (after `has_unconsumed_dot_type_keys`), add:

```cpp
struct LegacyFieldMap {
  const char* legacy;       // legacy CSV key, e.g. "m_ncdm"
  const char* dot;          // per-instance dot field, e.g. "m"
};

constexpr LegacyFieldMap kLegacyStandardFields[] = {
    {"m_ncdm",                            "m"},
    {"T_ncdm",                            "T"},
    {"deg_ncdm",                          "deg"},
    {"Omega_ncdm",                        "Omega"},
    {"omega_ncdm",                        "omega"},
    {"ksi_ncdm",                          "ksi"},
    {"ncdm_psd_parameters",               "psd_parameters"},
    {"use_ncdm_psd_files",                "use_psd_file"},
    {"ncdm_psd_filenames",                "psd_filename"},
    {"quadrature_strategy_ncdm_standard", "quadrature_strategy"},
    {"N_momentum_bins_ncdm_standard",     "momenta_bins"},
    {"maximum_q_ncdm_standard",           "max_q"},
};

// Returns true if any key of the form "<prefix>.*" exists in pfc.
bool any_dot_key_with_prefix(const FileContent& pfc, const std::string& prefix) {
  bool found = false;
  pfc.for_each([&](const std::string& name, const std::string&, bool) {
    if (name.rfind(prefix + ".", 0) == 0) {
      found = true;
    }
  });
  return found;
}

void TransmuteLegacyStandardNcdmToDot(FileContent& pfc) {
  // 1. Resolve N from N_ncdm / N_ncdm_standard.
  int N1 = 0, N2 = 0;
  bool has1 = pfc.read_int("N_ncdm_standard", N1);
  bool has2 = pfc.read_int("N_ncdm",          N2);

  if (has1 && has2) {
    throw std::invalid_argument(
        "In input file, you can only enter one of N_ncdm_standard and N_ncdm, choose one");
  }
  int N = has1 ? N1 : (has2 ? N2 : 0);
  if (N <= 0) {
    return;  // nothing to transmute
  }

  // 2. Synthesised prefix collision check.
  for (int i = 1; i <= N; ++i) {
    const std::string prefix = "ncdm__" + std::to_string(i);
    if (any_dot_key_with_prefix(pfc, prefix)) {
      throw std::invalid_argument(
          "legacy NCDM transmutation collides with existing dot-instance '" + prefix +
          "': rename your dot instance or remove the legacy N_ncdm keys");
    }
  }

  // 3. For each legacy CSV field, broadcast or distribute per instance.
  for (const auto& fm : kLegacyStandardFields) {
    std::string raw_value;
    if (!pfc.read_string(fm.legacy, raw_value)) {
      continue;
    }

    std::vector<std::string> values = FileContent::split_csv(raw_value);
    if (values.empty()) {
      continue;
    }

    if (values.size() == 1u) {
      // broadcast scalar to all N instances
      for (int i = 1; i <= N; ++i) {
        pfc.set("ncdm__" + std::to_string(i) + "." + fm.dot, values[0]);
      }
    }
    else if (static_cast<int>(values.size()) == N) {
      for (int i = 0; i < N; ++i) {
        pfc.set("ncdm__" + std::to_string(i + 1) + "." + fm.dot, values[i]);
      }
    }
    else {
      throw std::invalid_argument(
          std::string(fm.legacy) + " has " + std::to_string(values.size()) +
          " entries but N_ncdm_standard = " + std::to_string(N) +
          "; expected 1 (broadcast) or " + std::to_string(N));
    }
  }

  // 4. Stamp the type field on each synthesised instance.
  for (int i = 1; i <= N; ++i) {
    pfc.set("ncdm__" + std::to_string(i) + ".type", "ncdm_standard");
  }
}
```

(The `psd_filename` field will be set on every instance even when only some have `use_psd_file=1`. The legacy semantics mapped a filename per `_TRUE_` flag in the same order — the dot-syntax form gives every instance its own filename slot, so for legacy inputs we have to resequence: only instances with `use_psd_file=1` should receive the filenames in order. Handle this in the `psd_filename` loop iteration instead of the generic broadcast above.)

- [ ] **Step 2: Special-case `psd_filename` distribution to match legacy semantics**

Replace the single body of the `for (const auto& fm : kLegacyStandardFields)` loop with this branch for the `psd_filename` field. After the loop's per-field handling, in the `if (raw_value)` branch, add:

```cpp
if (std::string(fm.dot) == "psd_filename") {
  // Legacy: ncdm_psd_filenames lists one filename per instance whose
  // use_ncdm_psd_files flag is _TRUE_, in order. We've already written the
  // dot-syntax use_psd_file flags, so we now read them back.
  std::vector<int> use_flags(N, 0);
  for (int i = 1; i <= N; ++i) {
    int flag = 0;
    pfc.read_int("ncdm__" + std::to_string(i) + ".use_psd_file", flag);
    use_flags[i - 1] = flag;
  }
  size_t fn_idx = 0;
  for (int i = 0; i < N; ++i) {
    if (use_flags[i] != 0) {
      if (fn_idx >= values.size()) {
        throw std::invalid_argument(
            "ncdm_psd_filenames has fewer entries than enabled use_ncdm_psd_files flags");
      }
      pfc.set("ncdm__" + std::to_string(i + 1) + ".psd_filename", values[fn_idx]);
      ++fn_idx;
    }
  }
  continue;  // skip the generic broadcast/distribute path
}
```

Place this branch *before* the broadcast/distribute logic so the filename never goes through that path.

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build. The transmuter is unused — still no behaviour change.

- [ ] **Step 4: Commit**

```bash
git add species/ncdm_species.cpp
git commit -m "species: add TransmuteLegacyStandardNcdmToDot helper"
```

---

## Phase 3 — Switch standard-NCDM CreateAll over to the new pathway

### Task 7: Rewrite `NCDMSpecies::CreateAll` to use transmuter + new constructor

**Files:**
- Modify: `species/ncdm_species.cpp`

- [ ] **Step 1: Replace the `CreateAll` body**

In `species/ncdm_species.cpp`, replace the entire body of `NCDMSpecies::CreateAll(const SpeciesBuildContext& ctx)` (today lines ~99–143) with:

```cpp
std::vector<Named> NCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  TransmuteLegacyStandardNcdmToDot(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_standard");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    // Mark <instance>.type as consumed (instances_with does not mark anything as read)
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);

    auto sp = std::make_unique<NCDMSpecies>(
        ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
```

- [ ] **Step 2: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 3: Smoke-run dot-syntax fixture**

Run: `./class test/dotsyntax_ncdm.ini`
Expected: exits 0. Compare key output values to a baseline taken before this commit (`./class test/dotsyntax_ncdm.ini` against the same fixture run on the previous commit if available; otherwise just verify exit code and that `output/dotsyntax_ncdm_*` files are produced).

- [ ] **Step 4: Smoke-run a legacy fixture**

Run:
```bash
./class explanatory.ini
```
Expected: exits 0. (`explanatory.ini` typically contains `N_ncdm = ...`-style legacy inputs — confirm this is exercising the legacy path by checking it has either `N_ncdm` or `N_ncdm_standard` set: `grep -n "^N_ncdm" explanatory.ini`. If not, skip this step and use one of the `test_*.ini` fixtures with legacy NCDM, e.g. `grep -l "^N_ncdm" *.ini test/*.ini | head -3`.)

- [ ] **Step 5: Add a multi-species legacy fixture for regression coverage**

Create `test/legacy_ncdm_two_species.ini`. Copy `test/dotsyntax_ncdm.ini` as a starting point, then replace the `nu1.*` block with:

```
N_ncdm = 2
m_ncdm = 0.06, 0.10
T_ncdm = 0.71611, 0.71611
deg_ncdm = 1.0, 1.0
quadrature_strategy_ncdm_standard = 2, 2
N_momentum_bins_ncdm_standard = 7, 7
maximum_q_ncdm_standard = 18.0, 18.0
root = output/legacy_ncdm_two_species_
```

- [ ] **Step 6: Run the new fixture**

Run: `./class test/legacy_ncdm_two_species.ini`
Expected: exits 0. Verify two ncdm species were created by inspecting the output: `head -1 output/legacy_ncdm_two_species_background.dat` should contain `(.)number_ncdm[0]` and `(.)number_ncdm[1]`.

- [ ] **Step 7: Run the existing mixed-mode fixture**

The fixture `test/dotsyntax_ncdm_mixed.ini` mixes `N_ncdm_standard = 1` with `nu1.type = ncdm_standard`. Today's code throws "cannot mix..." on this; after this task it should succeed and produce *two* ncdm species (one synthesised as `ncdm__1` from the legacy block — but with no per-field CSV, all fields fall through to defaults — and one as `nu1`). 

Run: `./class test/dotsyntax_ncdm_mixed.ini`
Expected: exits 0; output contains both `(.)number_ncdm[0]` and `(.)number_ncdm[1]` in the background table.

- [ ] **Step 8: Commit**

```bash
git add species/ncdm_species.cpp test/legacy_ncdm_two_species.ini
git commit -m "species: rewrite NCDMSpecies::CreateAll to transmute legacy CSV into dot-syntax"
```

### Task 8: Delete `synthesise_standard_ncdm_flat_keys` and `kStandardFields`

**Files:**
- Modify: `species/ncdm_species.cpp`

The int-index `NCDMSpecies` constructor stays alive for now — `NCDMInteractingSpecies::NCDMInteractingSpecies` still chains to it (`NCDMSpecies(pfc, species_index, settings, pba, bgm, "_interacting")` at `ncdm_interacting_species.cpp:99`). It will be deleted in Task 16 once the interacting path is migrated. Only the dead synthesis helpers are removed in this task.

- [ ] **Step 1: Audit remaining callers of the int-index NCDMSpecies constructor**

Run: `grep -rn "make_unique<NCDMSpecies>\|new NCDMSpecies\|: NCDMSpecies(" species/ source/`
Expected: only `species/ncdm_interacting_species.cpp` (the chained constructor call) plus the new instance-name CreateAll. If a third caller exists, list it before continuing — it likely needs to be migrated to the instance-name form.

- [ ] **Step 2: Delete the dead helpers**

In `species/ncdm_species.cpp`, delete the following from the anonymous namespace:
- `struct FieldMap` and the `kStandardFields[]` constant array
- the `synthesise_standard_ncdm_flat_keys` function
- the `has_unconsumed_dot_type_keys` function (no longer referenced now that `NCDMSpecies::CreateAll` doesn't branch on dot-vs-legacy)

Confirm `has_unconsumed_dot_type_keys` is unused first: `grep -n "has_unconsumed_dot_type_keys" species/ncdm_species.cpp` — there should be no remaining references after Task 7.

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 4: Smoke-run NCDM fixtures**

```bash
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
```
Expected: all exit 0.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_species.cpp
git commit -m "species: drop synthesise_standard_ncdm_flat_keys and kStandardFields"
```

---

## Phase 4 — DNCDM (drop legacy entirely)

### Task 9: Add a new `DNCDMSpecies` constructor that takes an instance name

**Files:**
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`

This constructor reads everything DNCDM-specific (`Gamma`/`log10Gamma`/`lifetime`/`log10lifetime`, `Omega`/`omega`/`Omega_ini`/`omega_ini`/`Neff_ini`) directly via `SpeciesInput`. The Omega/deg closure logic that today lives in `DNCDMSpecies::CreateAll`'s bulk loop moves into the constructor (operating on the single instance).

- [ ] **Step 1: Declare the new constructor**

In `species/dncdm_species.h`, in the `public:` section, add:

```cpp
DNCDMSpecies(FileContent* pfc,
             const std::string& instance_name,
             const NcdmSettings& settings,
             const background* pba,
             const BackgroundModule* bgm);
```

- [ ] **Step 2: Implement the new constructor**

In `species/dncdm_species.cpp`, after the existing constructor (line ~115), add:

```cpp
DNCDMSpecies::DNCDMSpecies(FileContent* pfc,
                           const std::string& instance_name,
                           const NcdmSettings& settings,
                           const background* pba,
                           const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name,
                      EnergyType::Other,
                      pfc,
                      instance_name,
                      settings),
      pba_(pba) {
  bgm_ = bgm;

  SpeciesInput input(pfc, instance_name);

  // Read decay rate (exactly one of these must be specified)
  double Gamma_value      = 0.;
  double log10Gamma_value = 0.;
  double lifetime_value   = 0.;
  double log10lifetime_value = 0.;
  bool has_G    = input.read_double("Gamma",          Gamma_value);
  bool has_lG   = input.read_double("log10Gamma",     log10Gamma_value);
  bool has_lt   = input.read_double("lifetime",       lifetime_value);
  bool has_llt  = input.read_double("log10lifetime",  log10lifetime_value);

  int n_provided = (int)has_G + (int)has_lG + (int)has_lt + (int)has_llt;
  if (n_provided != 1) {
    throw std::invalid_argument(
        "species '" + instance_name +
        "': specify exactly one of Gamma, log10Gamma, lifetime, log10lifetime");
  }

  double Gamma_raw = 0.;
  if (has_G) {
    Gamma_raw = Gamma_value;
  }
  else if (has_lG) {
    Gamma_raw = std::pow(10., log10Gamma_value);
  }
  else if (has_lt) {
    Gamma_raw = 1. / lifetime_value / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  else {  // has_llt
    double lifetime = std::pow(10., log10lifetime_value);
    Gamma_raw       = 1. / lifetime / (365 * 24 * 60 * 60) * _Mpc_over_m_ * 1e-3;
  }
  Gamma_ = Gamma_raw * (1.e3 / _c_);

  // Read DNCDM-specific Omega / deg / Omega_ini values from the same instance.
  // These were previously bulk-loaded in CreateAll from CSV lists.
  // Semantics: at most one of {Omega, omega} and at most one of
  // {deg, Omega_ini, omega_ini, Neff_ini} (the latter group resolved later
  // by ConstructSpecies for shooting).
  double Omega_local = 0., omega_local = 0.;
  bool has_Omega  = input.read_double("Omega",  Omega_local);
  bool has_omega  = input.read_double("omega",  omega_local);
  if (has_Omega && has_omega) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of Omega, omega");
  }
  if (has_omega) {
    Omega_local = omega_local / settings.h / settings.h;
    has_Omega   = true;
  }
  if (has_Omega) {
    SetOmega0(Omega_local, settings.h);
  }

  double deg_local       = 0.;
  double Omega_ini_local = 0.;
  double omega_ini_local = 0.;
  double Neff_ini_local  = 0.;
  bool has_deg       = input.read_double("deg",       deg_local);
  bool has_Omega_ini = input.read_double("Omega_ini", Omega_ini_local);
  bool has_omega_ini = input.read_double("omega_ini", omega_ini_local);
  bool has_Neff_ini  = input.read_double("Neff_ini",  Neff_ini_local);

  int n_deg_options = (int)has_deg + (int)has_Omega_ini + (int)has_omega_ini + (int)has_Neff_ini;
  if (n_deg_options > 1) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify at most one of deg, Omega_ini, omega_ini, Neff_ini");
  }

  if (has_deg) {
    SetDegAndFactor(deg_local);
  }
  // The Omega_ini / omega_ini / Neff_ini closure runs in CreateAll after the
  // species exists (it needs pba->a_today and ppr->tol_ncdm). Stash the chosen
  // option here; CreateAll calls SetDeg_from_Omega_ini.
  if (has_Omega_ini) {
    Omega_ini_pending_ = Omega_ini_local;  // already in Omega units
  }
  else if (has_omega_ini) {
    Omega_ini_pending_ = omega_ini_local / settings.h / settings.h;
  }
  else if (has_Neff_ini) {
    Neff_ini_pending_  = Neff_ini_local;
  }

  // Compute dq_[i] = w_bg_[i] / f0(q_bg_[i]) (matches existing semantics)
  dq_ = ComputeDq();
}
```

- [ ] **Step 3: Add the pending-stash members and accessors**

In `species/dncdm_species.h`, add `<optional>` to the includes (or use sentinel-based flags if the codebase avoids `std::optional`).

In the `private:` section, add:
```cpp
std::optional<double> Omega_ini_pending_;  // for SetDeg_from_Omega_ini in CreateAll
std::optional<double> Neff_ini_pending_;
```

In the `public:` section, add accessors:
```cpp
const std::optional<double>& Omega_ini_pending() const { return Omega_ini_pending_; }
const std::optional<double>& Neff_ini_pending() const { return Neff_ini_pending_; }
```

Update Step 2's body to assign via `Omega_ini_pending_ = Omega_ini_local;` (the `optional<double>::operator=(double)` in-place sets the value). Same for `Neff_ini_pending_`.

- [ ] **Step 4: Build**

Run: `make -j class`
Expected: clean build. New constructor unused.

- [ ] **Step 5: Smoke-run baseline**

Run: `./class test/dotsyntax_ncdm.ini`
Expected: exits 0.

- [ ] **Step 6: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp
git commit -m "species: add instance-name DNCDMSpecies constructor"
```

### Task 10: Add a DNCDM dot-only fixture

**Files:**
- Create: `test/dotsyntax_dncdm.ini`

- [ ] **Step 1: Locate an existing DNCDM legacy fixture for reference**

Run: `grep -l "N_ncdm_decay_dr\|Gamma_ncdm_decay_dr" *.ini test/*.ini 2>/dev/null | head -5`
Pick one as a template. If none found, use `explanatory.ini` and search for the DNCDM section to read what valid legacy values look like.

- [ ] **Step 2: Write the dot-only fixture**

Create `test/dotsyntax_dncdm.ini` with cosmology copied from `test/dotsyntax_ncdm.ini` and a single decay_dr species block:

```
# (cosmology block copied from test/dotsyntax_ncdm.ini up to the species block)
h = 0.67556
T_cmb = 2.7255
omega_b = 0.022032
N_ur = 3.046
omega_cdm = 0.12038
Omega_dcdmdr = 0.0
Omega_k = 0.
Omega_fld = 0
Omega_scf = 0
Omega_idm_dr = 0.
stat_f_idr = 0.875
YHe = BBN
recombination = RECFAST
output = tCl,pCl,lCl
modes = s
lensing = yes
ic = ad
gauge = synchronous
P_k_ini type = analytic_Pk
k_pivot = 0.05
A_s = 2.215e-9
n_s = 0.9619
l_max_scalars = 2500
input_verbose = 1
background_verbose = 1

dncdm1.type = ncdm_decay_dr
dncdm1.m = 1.0
dncdm1.T = 0.71611
dncdm1.Gamma = 1e-3
dncdm1.Omega = 0.001

root = output/dotsyntax_dncdm_
```

(Fine-tune `Gamma` / `Omega` to match a known-stable existing legacy DNCDM test if available.)

- [ ] **Step 3: Commit**

```bash
git add test/dotsyntax_dncdm.ini
git commit -m "test: add dot-syntax DNCDM fixture"
```

### Task 11: Rewrite `DNCDMSpecies::CreateAll` and add `RejectLegacyDecayDrKeys`

**Files:**
- Modify: `species/dncdm_species.cpp`
- Modify: `species/dncdm_dr_species.cpp` (delegation)

- [ ] **Step 1: Add `RejectLegacyDecayDrKeys` to the anonymous namespace**

In `species/dncdm_species.cpp`, in the existing anonymous namespace at the top, add:

```cpp
constexpr const char* kLegacyDecayDrKeys[] = {
    "N_ncdm_decay_dr",
    "m_ncdm_decay_dr",
    "T_ncdm_decay_dr",
    "ksi_ncdm_decay_dr",
    "deg_ncdm_decay_dr",
    "quadrature_strategy_ncdm_decay_dr",
    "N_momentum_bins_ncdm_decay_dr",
    "maximum_q_ncdm_decay_dr",
    "Gamma_ncdm_decay_dr",
    "log10Gamma_ncdm_decay_dr",
    "lifetime_ncdm_decay_dr",
    "log10lifetime_ncdm_decay_dr",
    "Omega_dncdmdr",
    "omega_dncdmdr",
    "Omega_ini_dncdm",
    "omega_ini_dncdm",
    "Neff_ini_dncdm",
};

void RejectLegacyDecayDrKeys(FileContent& pfc) {
  for (const char* key : kLegacyDecayDrKeys) {
    std::string unused;
    if (pfc.read_string(key, unused)) {
      throw std::invalid_argument(
          std::string("'") + key + "' is no longer supported. Use dot-syntax: " +
          "'<instance>.<dot-name> = ...' with '<instance>.type = ncdm_decay_dr'.");
    }
  }
}
```

- [ ] **Step 2: Replace the body of `DNCDMSpecies::CreateAll`**

Replace the entire body of `DNCDMSpecies::CreateAll` (today lines ~179–337) with:

```cpp
std::vector<DNCDMSpecies::Named> DNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyDecayDrKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_decay_dr");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);  // mark consumed

    auto sp = std::make_unique<DNCDMSpecies>(
        ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    ApplyDncdmInitialClosure(*sp, ctx);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
```

The `ApplyDncdmInitialClosure` helper (added in Step 3) runs the `SetDeg_from_Omega_ini` closure for species whose `Omega_ini_pending_` or `Neff_ini_pending_` is set. This preserves the existing semantics (today at `dncdm_species.cpp:296–330`) but operates per species instead of in a bulk indexed loop.

- [ ] **Step 3: Implement `ApplyDncdmInitialClosure` in the anonymous namespace**

In the same `species/dncdm_species.cpp` anonymous namespace, add (after `RejectLegacyDecayDrKeys`):

```cpp
void ApplyDncdmInitialClosure(DNCDMSpecies& sp, const SpeciesBuildContext& ctx) {
  const auto& Omega_ini = sp.Omega_ini_pending();
  const auto& Neff_ini  = sp.Neff_ini_pending();
  if (!Omega_ini.has_value() && !Neff_ini.has_value()) {
    return;
  }

  const double a_ini_seed = ctx.pba ? ctx.pba->a_today * 1e-14 : 1e-14;
  const double a_today    = ctx.pba ? ctx.pba->a_today : 1.;
  const double a_ini      = sp.GetIni(a_ini_seed, a_today, ctx.ncdm_settings->tol_ncdm);
  const double z_ini      = 1.0 / a_ini - 1.0;
  const double H0         = ctx.pba ? ctx.pba->H0 : ctx.ncdm_settings->h * 1.e5 / _c_;

  if (Omega_ini.has_value()) {
    sp.SetDeg_from_Omega_ini(z_ini, H0, *Omega_ini);
  }
  else {
    // Neff_ini → Omega_ini using photon energy density at the same epoch.
    const double Omega0_g = ctx.pba ? ctx.pba->Omega0_g : 0.;
    const double Omega_ini_value =
        *Neff_ini * 7. / 8. * std::pow(4. / 11., 4. / 3.) * Omega0_g;
    sp.SetDeg_from_Omega_ini(z_ini, H0, Omega_ini_value);
  }
}
```

This is the per-species version of today's bulk loop at `dncdm_species.cpp:296–330`. The `SetDeg_from_Omega_ini` member already exists on `NCDMBaseSpecies`. Since `ApplyDncdmInitialClosure` calls public accessors and `SetDeg_from_Omega_ini` (which is `protected` on the base — verify with `grep -n "SetDeg_from_Omega_ini" species/ncdm_base_species.h`), if access is blocked, expose `SetDeg_from_Omega_ini` as `public` on `NCDMBaseSpecies` (move it from `protected:` to `public:` in `ncdm_base_species.h`) and add `species: expose SetDeg_from_Omega_ini for DNCDM closure` as a separate small commit before this one.

- [ ] **Step 4: Update `dncdm_dr_species.cpp` delegation**

In `species/dncdm_dr_species.cpp`, the `DNCDM_DR_Species::CreateAll` delegates to `DNCDMSpecies::CreateAll`. Confirm it still compiles. If `DNCDMSpecies::CreateAll`'s return type or wrapping has changed, adjust accordingly.

Run: `grep -n "DNCDMSpecies::CreateAll" species/dncdm_dr_species.cpp`
Read the surrounding code; adjust the call shape as needed.

- [ ] **Step 5: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 6: Smoke-run dot-only DNCDM fixture**

Run: `./class test/dotsyntax_dncdm.ini`
Expected: exits 0; produces `output/dotsyntax_dncdm_*` files. Inspect the background table header: `head -1 output/dotsyntax_dncdm_background.dat` should mention `dncdm1` (or similar).

- [ ] **Step 7: Smoke-run a legacy DNCDM fixture, expect rejection**

If you found a legacy DNCDM fixture in Task 10 Step 1, run it now:
Run: `./class <legacy_dncdm_fixture>.ini`
Expected: exits non-zero with an error message starting with `"'N_ncdm_decay_dr' is no longer supported. Use dot-syntax..."`.

- [ ] **Step 8: Migrate any tracked legacy DNCDM fixtures**

Run: `git ls-files | xargs grep -l "N_ncdm_decay_dr\|Gamma_ncdm_decay_dr\|Omega_dncdmdr" 2>/dev/null`
For each file in the tracked repo (not output/), rewrite the legacy DNCDM block into dot-syntax. Verify each runs.

- [ ] **Step 9: Smoke-run NCDM fixtures (regression)**

```bash
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
```
Expected: all exit 0.

- [ ] **Step 10: Commit**

```bash
git add species/dncdm_species.cpp species/dncdm_dr_species.cpp [any migrated fixtures]
git commit -m "species: rewrite DNCDMSpecies::CreateAll as dot-only, reject legacy keys"
```

### Task 12: Delete the legacy `DNCDMSpecies` constructor and `synthesise_decay_dr_ncdm_flat_keys`

**Files:**
- Modify: `species/dncdm_species.h`
- Modify: `species/dncdm_species.cpp`

- [ ] **Step 1: Audit remaining callers of the old DNCDM constructor**

Run: `grep -rn "make_unique<DNCDMSpecies>\|new DNCDMSpecies" species/ source/`
Expected: only `dncdm_dr_species.cpp` and `dncdm_species.cpp` itself. If something else creates DNCDMSpecies with the int-index form, address before deleting.

- [ ] **Step 2: Delete the old constructor**

In `species/dncdm_species.h`, remove the declaration with `int ncdm_index, int dncdm_index` parameters.

In `species/dncdm_species.cpp`, remove the implementation (today lines ~103–173).

- [ ] **Step 3: Delete `synthesise_decay_dr_ncdm_flat_keys`**

In `species/dncdm_species.cpp`, remove `synthesise_decay_dr_ncdm_flat_keys` (today line ~69) and its associated field-mapping table if any.

- [ ] **Step 4: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 5: Smoke-run all DNCDM and NCDM fixtures**

```bash
./class test/dotsyntax_dncdm.ini
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
```
Expected: all exit 0.

- [ ] **Step 6: Commit**

```bash
git add species/dncdm_species.h species/dncdm_species.cpp
git commit -m "species: drop legacy DNCDMSpecies constructor and synthesise_decay_dr_ncdm_flat_keys"
```

---

## Phase 5 — NCDMInteracting (drop legacy entirely)

### Task 13: Add a new `NCDMInteractingSpecies` constructor

**Files:**
- Modify: `species/ncdm_interacting_species.h`
- Modify: `species/ncdm_interacting_species.cpp`

- [ ] **Step 1: Declare the new constructor**

In `species/ncdm_interacting_species.h`, in the `public:` section, add:

```cpp
NCDMInteractingSpecies(FileContent* pfc,
                       const std::string& instance_name,
                       const NcdmSettings& settings,
                       const background* pba,
                       const BackgroundModule* bgm);
```

- [ ] **Step 2: Implement the new constructor**

In `species/ncdm_interacting_species.cpp`, add:

```cpp
NCDMInteractingSpecies::NCDMInteractingSpecies(FileContent* pfc,
                                               const std::string& instance_name,
                                               const NcdmSettings& settings,
                                               const background* pba,
                                               const BackgroundModule* bgm)
    : NCDMSpecies(pfc, instance_name, settings, pba, bgm) {
  SpeciesInput input(pfc, instance_name);

  double G_eff_value      = 0.;
  double log10G_eff_value = 0.;
  bool has_G  = input.read_double("G_eff",      G_eff_value);
  bool has_lG = input.read_double("log10G_eff", log10G_eff_value);

  if (has_G && has_lG) {
    throw std::invalid_argument("species '" + instance_name +
                                "': specify exactly one of G_eff or log10G_eff");
  }
  if (has_G) {
    G_eff_ = G_eff_value;
  }
  else if (has_lG) {
    G_eff_ = std::pow(10.0, log10G_eff_value);
  }
}
```

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 4: Smoke-run**

Run: `./class test/dotsyntax_ncdm.ini`
Expected: exits 0.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_interacting_species.h species/ncdm_interacting_species.cpp
git commit -m "species: add instance-name NCDMInteractingSpecies constructor"
```

### Task 14: Add an NCDMInteracting dot-only fixture

**Files:**
- Create: `test/dotsyntax_ncdm_interacting.ini`

- [ ] **Step 1: Find an existing legacy interacting fixture for reference**

Run: `grep -l "N_ncdm_interacting\|G_eff_ncdm_interacting" *.ini test/*.ini 2>/dev/null | head -3`

- [ ] **Step 2: Write the fixture**

Create `test/dotsyntax_ncdm_interacting.ini` based on `test/dotsyntax_ncdm.ini`, replacing the species block with:

```
nui1.type = ncdm_self_interacting
nui1.m = 0.06
nui1.T = 0.71611
nui1.deg = 1.0
nui1.G_eff = 1e-3

root = output/dotsyntax_ncdm_interacting_
```

- [ ] **Step 3: Commit**

```bash
git add test/dotsyntax_ncdm_interacting.ini
git commit -m "test: add dot-syntax NCDMInteracting fixture"
```

### Task 15: Rewrite `NCDMInteractingSpecies::CreateAll` and add `RejectLegacyInteractingKeys`

**Files:**
- Modify: `species/ncdm_interacting_species.cpp`

- [ ] **Step 1: Add `RejectLegacyInteractingKeys`**

In `species/ncdm_interacting_species.cpp`, in the anonymous namespace at the top, add:

```cpp
constexpr const char* kLegacyInteractingKeys[] = {
    "N_ncdm_interacting",
    "G_eff_ncdm_interacting",
    "log10G_eff_ncdm_interacting",
    "quadrature_strategy_ncdm_interacting",
    "N_momentum_bins_ncdm_interacting",
    "maximum_q_ncdm_interacting",
};

void RejectLegacyInteractingKeys(FileContent& pfc) {
  for (const char* key : kLegacyInteractingKeys) {
    std::string unused;
    if (pfc.read_string(key, unused)) {
      throw std::invalid_argument(
          std::string("'") + key + "' is no longer supported. Use dot-syntax: " +
          "'<instance>.<dot-name> = ...' with '<instance>.type = ncdm_self_interacting'.");
    }
  }
}
```

- [ ] **Step 2: Replace the `CreateAll` body**

Replace the body of `NCDMInteractingSpecies::CreateAll` (today lines ~126–163) with:

```cpp
std::vector<Named> NCDMInteractingSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  RejectLegacyInteractingKeys(*ctx.pfc);

  const auto instances = ctx.pfc->instances_with("type", "ncdm_self_interacting");

  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);  // mark consumed

    auto sp = std::make_unique<NCDMInteractingSpecies>(
        ctx.pfc, name, *ctx.ncdm_settings, ctx.pba, ctx.bgm);
    sp->SetNcdmId((*ctx.ncdm_id_next)++);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
```

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 4: Smoke-run**

```bash
./class test/dotsyntax_ncdm_interacting.ini
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
./class test/dotsyntax_dncdm.ini
```
Expected: all exit 0.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_interacting_species.cpp
git commit -m "species: rewrite NCDMInteractingSpecies::CreateAll as dot-only, reject legacy keys"
```

### Task 16: Delete the legacy `NCDMInteractingSpecies` and `NCDMSpecies` int-index constructors and `synthesise_self_interacting_ncdm_flat_keys`

**Files:**
- Modify: `species/ncdm_interacting_species.h`
- Modify: `species/ncdm_interacting_species.cpp`
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp`

- [ ] **Step 1: Audit remaining callers**

```bash
grep -rn "make_unique<NCDMInteractingSpecies>\|new NCDMInteractingSpecies" species/ source/
grep -rn "make_unique<NCDMSpecies>\|new NCDMSpecies" species/ source/
```
Expected: only the new instance-name forms remain. If any int-index call survives, address before deleting.

- [ ] **Step 2: Delete the old `NCDMInteractingSpecies` int-index constructor**

In `species/ncdm_interacting_species.h`, remove the int-index constructor declaration (today line ~95).

In `species/ncdm_interacting_species.cpp`, remove the implementation (today lines ~94–123).

- [ ] **Step 3: Delete `synthesise_self_interacting_ncdm_flat_keys`**

In `species/ncdm_interacting_species.cpp`, remove the function definition.

- [ ] **Step 4: Delete the old `NCDMSpecies` int-index constructor**

In `species/ncdm_species.h`, remove the int-index constructor declaration (today lines ~15–20).

In `species/ncdm_species.cpp`, remove its implementation (today lines ~12–26).

- [ ] **Step 5: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 6: Smoke-run all NCDM-family fixtures**

```bash
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
./class test/dotsyntax_dncdm.ini
./class test/dotsyntax_ncdm_interacting.ini
./class explanatory.ini
```
Expected: all exit 0.

- [ ] **Step 7: Commit**

```bash
git add species/ncdm_interacting_species.h species/ncdm_interacting_species.cpp species/ncdm_species.h species/ncdm_species.cpp
git commit -m "species: drop legacy int-index NCDM/NCDMInteracting constructors and synthesise helpers"
```

---

## Phase 6 — Base-class cleanup and final tests

### Task 17: Delete the legacy `NCDMBaseSpecies::ReadParameters` / `ReadDecayDrParameters`

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_base_species.cpp`

- [ ] **Step 1: Audit which constructors of `NCDMBaseSpecies` are still used**

Run: `grep -rn "NCDMBaseSpecies(" species/ | grep -v ".worktrees"`
Expected: only the new `(name, energy_type, pfc, instance_name, settings)` constructor is invoked. The two legacy overloads are dead.

- [ ] **Step 2: Delete the two legacy constructor declarations and definitions**

In `species/ncdm_base_species.h`, remove:
```cpp
NCDMBaseSpecies(std::string name, EnergyType energy_type, FileContent* pfc,
                int species_index, const NcdmSettings& settings,
                const std::string& suffix = "_standard");
NCDMBaseSpecies(std::string name, EnergyType energy_type, FileContent* pfc,
                int dncdm_index, const NcdmSettings& settings, bool is_decay_dr);
```

In `species/ncdm_base_species.cpp`, delete both implementations and the bodies of `ReadParameters` and `ReadDecayDrParameters`.

- [ ] **Step 3: Remove now-unused declarations and helpers**

In `species/ncdm_base_species.h`, remove:
```cpp
void ReadParameters(...);
void ReadDecayDrParameters(...);
```

In `species/ncdm_base_species.cpp`, remove the anonymous-namespace `readDoubleList`, `readIntegerList`, `readStringList` if no longer used.

Run: `grep -n "readDoubleList\|readIntegerList\|readStringList" species/ncdm_base_species.cpp` to confirm they're dead.

- [ ] **Step 4: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 5: Smoke-run all NCDM-family fixtures**

```bash
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
./class test/dotsyntax_dncdm.ini
./class test/dotsyntax_ncdm_interacting.ini
./class explanatory.ini
```
Expected: all exit 0.

- [ ] **Step 6: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp
git commit -m "species: drop legacy NCDMBaseSpecies CSV-based constructors and ReadParameters"
```

### Task 18: Delete unused `species_input` helpers

**Files:**
- Modify: `species/species_input.h`
- Modify: `species/species_input.cpp`

- [ ] **Step 1: Audit which helpers still have callers**

```bash
grep -rn "CsvWithDefaults\|CsvForPsdFilenames\|SynthesiseIdenticalScalarField\|CollectInstanceFieldValues\|AnyInstanceFieldValue" species/ source/ tools/ | grep -v ".worktrees"
```

For each helper, if there are zero non-self matches (i.e. only the definition itself), it's dead and removable.

- [ ] **Step 2: Delete the dead helpers**

For each dead helper found in Step 1, remove its declaration from `species/species_input.h` and its definition from `species/species_input.cpp`.

The `SpeciesInput` class itself stays — it is the new read API and is now the primary user-facing interface in this header.

- [ ] **Step 3: Build**

Run: `make -j class`
Expected: clean build.

- [ ] **Step 4: Smoke-run all NCDM fixtures**

```bash
./class test/dotsyntax_ncdm.ini
./class test/legacy_ncdm_two_species.ini
./class test/dotsyntax_ncdm_mixed.ini
./class test/dotsyntax_dncdm.ini
./class test/dotsyntax_ncdm_interacting.ini
```
Expected: all exit 0.

- [ ] **Step 5: Commit**

```bash
git add species/species_input.h species/species_input.cpp
git commit -m "species: drop unused dot→CSV synthesis helpers from species_input"
```

### Task 19: Run the full integration test suite

**Files:**
- (no changes; verification only)

- [ ] **Step 1: Run the pytest scenario suite**

```bash
cd python
TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py
cd ..
```
Expected: all tests pass. Investigate any failures — likely candidates are scenarios that use legacy DNCDM or legacy NCDMInteracting keys (now rejected). Migrate those scenarios to dot-syntax in `python/test_class.py`'s `CLASS_INPUT` dictionary if any exist.

- [ ] **Step 2: Run wrapper smoke**

```bash
pip install -e .
python -c "import classy; c = classy.Class(); c.set({'output': 'tCl', 'N_ncdm': 1, 'm_ncdm': 0.06}); c.compute(); print(c.lensed_cl(2)['ell'][:5])"
```
Expected: prints a small array; no exception.

- [ ] **Step 3: Commit any test-file migrations**

If `python/test_class.py` was modified to drop legacy DNCDM/NCDMInteracting input forms:
```bash
git add python/test_class.py
git commit -m "test: migrate legacy DNCDM/NCDMInteracting test scenarios to dot-syntax"
```

---

## Self-review checklist

After completing all tasks, before declaring done:

- [ ] No constructor in `NCDMSpecies` / `DNCDMSpecies` / `NCDMInteractingSpecies` / `NCDMBaseSpecies` takes an `int species_index` or `int ncdm_index` or `int dncdm_index`.
- [ ] `NCDMSpecies::CreateAll`, `DNCDMSpecies::CreateAll`, and `NCDMInteractingSpecies::CreateAll` are each ≤ ~15 lines (the discover-and-loop shape).
- [ ] `synthesise_standard_ncdm_flat_keys`, `synthesise_decay_dr_ncdm_flat_keys`, `synthesise_self_interacting_ncdm_flat_keys`, `kStandardFields`, `has_unconsumed_dot_type_keys` are gone.
- [ ] `test/legacy_ncdm_two_species.ini` (multi-species legacy CSV) passes.
- [ ] `test/dotsyntax_ncdm_mixed.ini` (legacy + dot mixed) passes.
- [ ] `test/dotsyntax_dncdm.ini` (DNCDM dot-only) passes.
- [ ] `test/dotsyntax_ncdm_interacting.ini` (NCDMInteracting dot-only) passes.
- [ ] Running an `.ini` file with `N_ncdm_decay_dr = ...` produces an error message starting with `"'N_ncdm_decay_dr' is no longer supported"`.
- [ ] Running an `.ini` file with `N_ncdm_interacting = ...` produces an analogous error.
- [ ] `python -m pytest -v -m test_scenario test_class.py` passes at TEST_LEVEL=1.
- [ ] `unread_parameters()` no longer reports false negatives for legacy NCDM keys after transmutation.

If any check fails, fix and re-run before declaring complete.
