# Species Index Encapsulation Design

## Goal

Remove all species-dependent `index_bg_rho_X_` mirror fields from `BackgroundModule`, eliminate `bg_rho_index()` / `bg_p_index()` from `BaseSpecies`, move background output to per-species virtuals, and replace manual `rho_r/rho_m/rho_nu` accumulation in perturbation initial conditions with dispatch. The result: every species owns its background index fully privately, and no external code can construct a pvecback offset by hand.

## Motivation

The PR that removed `pba->has_*` dispatch guards revealed a latent class of bug: `CDMSpecies`, `LambdaSpecies`, `DCDMSpecies`, and `FluidSpecies` each stored their background density index in a species-private field (`index_bg_rho_cdm_`, etc.) but never assigned `BaseSpecies::index_bg_rho_`. `bg_rho_index()` returned -1, so `BackgroundModule::index_bg_rho_cdm_` was always -1, and `pvecback[-1]` was accessed on every perturbation run — a heap buffer overflow caught only by ASan. The one-liner fix (`index_bg_rho_ = index_bg_rho_cdm_`) closes that specific gap, but the underlying system remains fragile: two parallel representations of the same integer can drift. This design eliminates the duplication entirely.

## Architecture

Two goals, executed sequentially in one branch:

1. **Background index encapsulation** — remove mirror fields from `BackgroundModule`, delete `bg_rho_index()` / `bg_p_index()` from `BaseSpecies`, migrate all call sites to `sp->Rho(pvecback)` / `sp->P(pvecback)`, move background output columns into per-species virtuals.
2. **IC accumulation via dispatch** — replace the manual `rho_r / rho_m / rho_nu` accumulation block in the perturbation IC setup with a species dispatch loop using `EnergyType` plus a free-streaming radiation contribution that aggregates correctly across composites.

**Migration strategy: Approach 2 (compiler-guided).** Remove the interface (fields + accessor methods) in a single commit; the compiler enumerates every call site. Fix module by module. No call site can be missed.

---

## Section 1: `BaseSpecies` interface changes

### Removed

- `int bg_rho_index() const` — returns `index_bg_rho_`
- `int bg_p_index() const` — returns `index_bg_p_`

The fields `index_bg_rho_` and `index_bg_p_` **stay `protected`** — subclasses still write to them directly in `RegisterBackgroundIndices` (e.g. `index_bg_rho_ = index_bg++`). What is removed is the *public getter*. External callers can no longer read the raw integer; `Rho(pvecback)` is the only external interface.

### Added

```cpp
/** True for free-streaming massless species (UR, DR, massless NCDM).
 *  Used by FreestreamingRho() for simple species. */
virtual bool IsFreestreaming() const { return false; }

/** Free-streaming radiation density contributed by this species at the
 *  current background state. Simple species derive this from
 *  IsFreestreaming(); CompositeSpecies sums over children. */
virtual double FreestreamingRho(const double* pvecback) const {
    return IsFreestreaming() ? Rho(pvecback) : 0.;
}

/** Write this species' background column titles into writer.
 *  Called by BackgroundModule during background_output_titles().
 *  Default: no-op. */
virtual void WriteBackgroundColumnTitles(BackgroundColumnWriter& writer) const {}

/** Write this species' background data values into writer.
 *  pvecback is the full background vector at the current time step.
 *  Called by BackgroundModule during background_output_data().
 *  Default: no-op. */
virtual void WriteBackgroundData(const double* pvecback,
                                 BackgroundColumnWriter& writer) const {}
```

`BackgroundColumnWriter` is a new thin struct in a new `source/background_column_writer.h` (keeping `background_module.h` from growing further), matching the pattern of `PerturbColumnWriter`:

```cpp
struct BackgroundColumnWriter {
  char*  titles;   // non-null in title mode
  double* data;    // non-null in data mode
  int&   index;
  bool   title_mode;
  void Add(const char* name, double value, bool condition = true);
};
```

### Unchanged

`Rho(pvecback)`, `P(pvecback)`, `DpDloga(pvecback)`, `EnergyType`, `RegisterBackgroundIndices`. These are the dispatch interface and remain as-is.

### `IsFreestreaming()` overrides

| Species | Returns |
|---------|---------|
| `UltraRelativisticSpecies` | `true` |
| `DarkRadiationSpecies` | `true` |
| `NCDMSpecies` | `true` (always free-streaming at IC time, deep in radiation domination) |
| `DNCDM_DecayRadiationSpecies` | `true` if it contributes free-streaming DR at IC time — implementer to verify |
| All others | `false` (default) |

---

## Section 2: `BackgroundModule` changes

### Fields removed from `background_module.h`

All species-dependent index mirrors:

```
index_bg_rho_g_          index_bg_rho_b_
index_bg_rho_cdm_        index_bg_rho_lambda_
index_bg_rho_fld_        index_bg_w_fld_       index_bg_dw_over_da_fld_
index_bg_rho_ur_
index_bg_rho_idm_dr_     index_bg_rho_idr_
index_bg_rho_dcdm_
index_bg_rho_idm_drmd_   index_bg_rho_idr_drmd_
index_bg_rho_dr_species_ index_bg_rho_dr_
index_bg_rho_scf_        index_bg_p_scf_        index_bg_p_prime_scf_
index_bg_rho_ncdm1_      index_bg_p_ncdm1_
```

**Not removed:** aggregated/physical module-owned fields — `index_bg_rho_tot_`, `index_bg_p_tot_`, `index_bg_p_tot_prime_`, `index_bg_rho_crit_`, `index_bg_H_`, `index_bg_a_`, `index_bg_Omega_r_`, `index_bg_conf_distance_`, etc. These belong to the module, not to any species.

### Index registration (`background_indices()`)

The manual per-species blocks:
```cpp
index_bg_rho_cdm_ = -1;
if (all_species_.count("CDM")) {
    all_species_.at("CDM")->RegisterBackgroundIndices(index_bg);
    index_bg_rho_cdm_ = all_species_.at("CDM")->bg_rho_index();  // deleted
}
```
become a uniform loop:
```cpp
for (auto& [name, sp] : all_species_)
    sp->RegisterBackgroundIndices(index_bg);
```

Composite sub-components register their own indices when the composite's `RegisterBackgroundIndices` is called (already the case via `CompositeSpecies::RegisterBackgroundIndices`).

### Background output

Species-specific `class_store_columntitle` / `class_store_double` calls in `background_output_titles()` and `background_output_data()` are removed. After the module writes its aggregate columns (a, H, rho_tot, etc.) it dispatches:

```cpp
for (auto& [name, sp] : all_species_)
    sp->WriteBackgroundColumnTitles(writer);  // title pass

for (auto& [name, sp] : all_species_)
    sp->WriteBackgroundData(pvecback, writer);  // data pass
```

Each species implements these to write its own density (and any additional quantities such as w_fld, phi_scf, etc.).

---

## Section 3: Call site migration

Removing the fields from `background_module.h` produces one compile error per call site. Fix in module order:

1. `source/background_module.cpp`
2. `source/perturbations_module.cpp` (18 sites)
3. `source/thermodynamics_module.cpp` (3 sites)
4. `source/input_module.cpp` (4 sites — background table lookups)

### Migration patterns

**Single species, current state:**
```cpp
// before
pvecback[background_module_->index_bg_rho_cdm_]
// after
all_species_.at("CDM")->Rho(pvecback)
```

**NCDM indexed by flavor:**
```cpp
// before
pvecback[background_module_->index_bg_rho_ncdm1_ + n]
// after
ncdm_species_sorted_[n]->Rho(pvecback)
```

**Composite sub-component:**
```cpp
// before
pvecback[background_module_->index_bg_rho_idm_dr_]
// after
static_cast<IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).idm_dr().Rho(pvecback)
```

**Background table lookup (`input_module.cpp`):**
```cpp
// before
background_table_[row * bg_size + background_module_->index_bg_rho_dcdm_]
// after
all_species_.at("DCDM_DR")->/* dcdm child */.Rho(background_table_ + row * bg_size)
```

**Absent-species guard (index `>= 0` sentinel disappears):**
```cpp
// before
if (background_module_->index_bg_rho_cdm_ >= 0)
    val += pvecback[background_module_->index_bg_rho_cdm_];
// after
if (all_species_.count("CDM"))
    val += all_species_.at("CDM")->Rho(pvecback);
```

---

## Section 4: IC accumulation via dispatch

Replaces the manual `rho_r / rho_m / rho_nu` accumulation block in `perturbations_module.cpp` (~lines 4700–4760) with a dispatch loop modelled on the `accumulate` lambda in `background_module.cpp:320`.

### New accumulation loop

```cpp
double rho_r  = all_species_.at("Photons")->Rho(pvecback);  // photons seed rho_r
double rho_m  = all_species_.at("Baryons")->Rho(pvecback);  // baryons seed rho_m
double rho_nu = 0.;

for (auto& [name, sp] : all_species_) {
    if (name == "Photons" || name == "Baryons") continue;  // already seeded
    const double rho = sp->Rho(pvecback);
    switch (sp->energy_type()) {
        case BaseSpecies::EnergyType::Matter:
            rho_m += rho; break;
        case BaseSpecies::EnergyType::Radiation:
            rho_r += rho; break;
        case BaseSpecies::EnergyType::Other:
            rho_r += 3. * sp->P(pvecback);
            rho_m += rho - 3. * sp->P(pvecback);
            break;
        default: break;  // DarkEnergy: no rho_r/rho_m contribution
    }
    rho_nu += sp->FreestreamingRho(pvecback);
}
```

`fracnu`, `fracg`, `fracb`, `rho_m_over_rho_r` are derived from these totals exactly as today. The `ic_ctx` structure and all per-species `ApplyInitialConditions` calls are untouched.

### Derived fractions

```cpp
const double rho_g   = all_species_.at("Photons")->Rho(pvecback);
const double rho_b   = all_species_.at("Baryons")->Rho(pvecback);
const double fracg   = rho_g / rho_r;
const double fracnu  = rho_nu / rho_r;
const double fracb   = rho_b / rho_m;

double fraccdm = 0.;
if (all_species_.count("CDM"))
    fraccdm = all_species_.at("CDM")->Rho(pvecback) / rho_m;

double fracidm_drmd = 0.;
if (all_species_.count("IDM_DRMD_IDR_DRMD")) {
    auto& drmd = static_cast<IDM_DRMD_IDR_DRMD_Species&>(*all_species_.at("IDM_DRMD_IDR_DRMD"));
    // idm_drmd() sub-component is only present when has_idm_drmd was true.
    // Requires a small IsPresent() (or equivalent) on the sub-component species —
    // this replaces the old index_bg_rho_idm_drmd_ >= 0 sentinel.
    if (drmd.idm_drmd().IsPresent())
        fracidm_drmd = drmd.idm_drmd().Rho(pvecback) / rho_m;
}
```

---

## Execution order

1. Add `BackgroundColumnWriter`, `IsFreestreaming()`, `WriteBackgroundColumnTitles()`, `WriteBackgroundData()` to `BaseSpecies` and implement in each species
2. Delete `bg_rho_index()`, `bg_p_index()`, and all mirror fields from `BackgroundModule` — let compiler enumerate call sites
3. Fix `background_module.cpp` (registration loop + output dispatch)
4. Fix `perturbations_module.cpp` (18 density access sites)
5. Fix `thermodynamics_module.cpp` (3 sites)
6. Fix `input_module.cpp` (4 background table sites)
7. Replace IC accumulation block with dispatch loop; add `IsFreestreaming()` overrides
8. Build, run 84 tests + ASan on `explanatory.ini`

## Testing

- `make class -j` must compile with zero warnings
- `cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -v -m test_scenario test_class.py` — all 84 must pass
- `./class explanatory.ini` under ASan — zero heap errors
- Background output columns for all species must appear correctly in output files
