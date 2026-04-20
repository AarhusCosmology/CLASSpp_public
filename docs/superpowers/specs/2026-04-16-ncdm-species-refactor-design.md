# NCDM Species Refactor Design

**Date:** 2026-04-16
**Context:** Dissolve the `NonColdDarkMatter` and `DarkRadiation` tool classes entirely into self-contained species classes, introduce `NCDMBaseSpecies` as a shared base for NCDM flavors, and give each species class a `CreateAll` factory so that `InputModule` becomes ignorant of per-species parameter details.

## 1. Motivation

The current design has two structural problems:

1. **`NonColdDarkMatter`** manages *all* NCDM species together in a single shared object. Each `NCDMSpecies` / `DNCDMSpecies` holds a `shared_ptr<NonColdDarkMatter>` and addresses its own data via an `n_ncdm` integer index. This means per-species data (mass, temperature, quadrature grid, distribution function) is scattered across vector-of-vectors indexed by species number rather than owned by each species.

2. **`input_module.cpp` lines 1022–1193** play orchestrator: after constructing `NonColdDarkMatter`, the input module reads DNCDM-specific parameters (`Omega_dncdmdr`, `deg_ncdm_decay_dr`, `Gamma_ncdm_decay_dr`, etc.) and calls `SetOmega0` / `SetDegAndFactor` / `SetDeg_from_Omega_ini` on the shared tool object. This logic belongs in the species themselves.

The goal is species that construct any number of themselves from `FileContent` via a `CreateAll` factory, own all their data, and leave `InputModule` with only a handful of calls.

A secondary benefit: when issue #246 (object-based input syntax) is implemented, only the `CreateAll` methods need to change — no ripple into `InputModule` or the rest of the species machinery.

## 2. Class Hierarchy

```
BaseSpecies
├── NCDMBaseSpecies               ← NEW: owns quadrature/distribution-function data
│   ├── NCDMSpecies               ← standard massive neutrinos / WDM
│   └── DNCDMSpecies              ← decaying NCDM (adds lnf, Gamma, decay props)
├── DNCDM_DecayRadiationSpecies   ← already self-contained, no changes needed
├── DarkRadiationSpecies          ← absorbs DarkRadiation tool (DCDM-only DR)
└── ... (CDM, baryons, photons, IDM, etc. — untouched)

CompositeSpecies
└── DNCDM_DR_Species              ← structurally unchanged; constructor updated
```

## 3. `NCDMBaseSpecies`

Abstract class. Does not implement `RegisterPerturbationIndices`, `PerturbDerivs`, `Delta`, `Theta`, `DeltaP`, `RhoPlusPShear` — those remain pure virtual and are implemented by concrete subclasses.

### 3.1 Data members

All single-species (no `n_ncdm` index):

```cpp
// Quadrature — perturbation sampling
std::vector<double> q_;
std::vector<double> w_;
std::vector<double> dlnf0_dlnq_;

// Quadrature — background sampling
std::vector<double> q_bg_;
std::vector<double> w_bg_;

double factor_;    // normalization factor for momentum integrals

// Species parameters
double m_in_eV_;
double M_;         // dimensionless mass = m / T_ncdm
double deg_;
double T_;         // T_ncdm / T_cmb
double ksi_;       // chemical potential / T_ncdm
```

### 3.2 Methods absorbed from `NonColdDarkMatter`

- `ComputeMomenta(double z, double* n, double* rho, double* p, double* drho_dM, double* pseudo_p)` — the old `background_ncdm_momenta(int n_ncdm, ...)` without the index parameter.
- `GetOmega0()`, `GetNeff(double z)`, `GetMassInElectronvolt()`
- `GetRescalingFactor(...)`, `GetRescaledParameters(...)`
- Private: `InitDistribution(FileContent*, int species_index)`, `InitQuadrature(...)`, and the static distribution function callbacks.

`background_parameters_for_distributions` becomes a private implementation detail inside `NCDMBaseSpecies`.

### 3.3 Construction

`NCDMBaseSpecies` has a `protected` constructor:

```cpp
NCDMBaseSpecies(FileContent* pfc, int species_index, const NcdmSettings&);
```

It reads the `species_index`-th entry from comma-separated input lists (`m_ncdm`, `T_ncdm`, `deg_ncdm`, `ksi_ncdm`, `ncdm_psd_parameters`, etc.), sets up quadrature, and initializes the distribution function. Subclasses call this via their own constructors.

`NcdmSettings` moves from `tools/non_cold_dark_matter.h` to `species/ncdm_base_species.h`.

## 4. `NCDMSpecies`

Thin subclass — inherits everything from `NCDMBaseSpecies`, adds only:

- Background/perturbation index members (`index_bg_number_`, `index_bg_pseudo_p_`, `index_pt_psi0_`).
- The full `BaseSpecies` virtual interface already implemented today (physics unchanged).
- Factory:

```cpp
static std::vector<std::unique_ptr<NCDMSpecies>>
    CreateAll(FileContent*, const NcdmSettings&);
```

Reads `N_ncdm` and the per-entry type list, skips `decay_dr` entries, constructs one instance per standard entry.

## 5. `DNCDMSpecies`

Inherits `NCDMBaseSpecies`, adds:

### 5.1 Decay-specific data (absorbed from `DecayDRProperties` + `NonColdDarkMatter`)

```cpp
double Gamma_;
std::vector<double> dq_;
int q_offset_;
int quadrature_strategy_decay_;
```

### 5.2 Background integration indices

Unchanged from today: `index_bi_lnf_decay_dr1_`, `index_bi_dlnfdlnq_separate_decay_`, `index_bg_lnf_decay_dr1_`, `index_bg_dlnfdlnq_decay_`, `index_bg_dlnfdlnq_sep_`.

### 5.3 Factory

```cpp
static std::vector<std::unique_ptr<DNCDMSpecies>>
    CreateAll(FileContent*, const NcdmSettings&);
```

Reads `N_ncdm`, skips standard entries, then reads the per-DNCDM parameters currently in `input_module.cpp` lines 1080–1190:
- `Omega_dncdmdr` / `omega_dncdmdr` → calls `SetOmega0` internally
- `deg_ncdm_decay_dr` / `Omega_ini_dncdm` / `omega_ini_dncdm` / `Neff_ini_dncdm` → calls deg-setting logic internally
- `Gamma_ncdm_decay_dr`

All validation (`class_test` checks) moves inside this factory.

## 6. `DNCDM_DR_Species`

Structurally unchanged. Its `CreateAll` factory:

```cpp
static std::vector<std::unique_ptr<DNCDM_DR_Species>>
    CreateAll(FileContent*, const NcdmSettings&);
```

Calls `DNCDMSpecies::CreateAll` internally and wraps each result in a composite with its `DNCDM_DecayRadiationSpecies`. `InputModule` calls only `DNCDM_DR_Species::CreateAll` — knowing nothing about the DNCDM/DR split.

## 7. `DarkRadiationSpecies` — Absorbing `DarkRadiation`

Since DNCDM-sourced DR is handled by `DNCDM_DecayRadiationSpecies`, `DarkRadiationSpecies` becomes the **DCDM-only DR species**. The multi-channel `N_species_` / `w_species_` / `cumulative_q_index_` design in `DarkRadiation` disappears.

### 7.1 Data members absorbed

```cpp
std::vector<double> q_;
std::vector<double> dq_;
std::vector<double> w_;
double deg_;
double factor_;
int N_q_;
DRType dr_type_;   // fermion / boson
```

### 7.2 Methods absorbed

- `IntegrateDistribution(double z, double* number, double* rho, double* p)` — single-channel, no `index_dr` parameter.
- `ComputeMeanMomentum()`
- Private: `InitialDistribution(...)`, quadrature setup.

### 7.3 Factory

```cpp
static std::unique_ptr<DarkRadiationSpecies>
    Create(FileContent*, const DCDMSpecies* dcdm);
```

Returns `nullptr` if `dcdm == nullptr` or no DR is requested (replaces `NoDrRequested` exception). The `DCDMSpecies` pointer replaces the current raw `dcdm_` pointer for coupling in background derivatives.

## 8. `InputModule` After Refactor

Lines 1022–1193 collapse to:

```cpp
// Standard NCDM
auto ncdm_species = NCDMSpecies::CreateAll(pfc, ncdm_settings);

// Decaying NCDM + decay radiation composites
auto dncdm_dr_species = DNCDM_DR_Species::CreateAll(pfc, ncdm_settings);

// DCDM dark radiation (nullptr if DCDM absent)
auto dr_species = DarkRadiationSpecies::Create(pfc, dcdm_species.get());
```

No more `SetOmega0`, `SetDegAndFactor`, `SetDeg_from_Omega_ini` calls. No more `dr_sources` / `dr_types` / `dr_deg` assembly. No more `ncdm_->N_ncdm_decay_dr_` checks.

`pba->Omega0_ncdm_tot` and `pba->N_ncdm` are derived from the species list in `BackgroundModule` (as they already are via the species loop), so those assignments in `InputModule` are removed.

## 9. Deletion Targets

| Item | Action |
|---|---|
| `tools/non_cold_dark_matter.h` / `.cpp` | Deleted |
| `tools/dark_radiation.h` / `.cpp` | Deleted |
| `NCDMType` enum | Deleted |
| `NcdmSettings` struct | Moves to `species/ncdm_base_species.h` |
| `DecayDRProperties` struct | Deleted (absorbed into `DNCDMSpecies`) |
| `NoNcdmRequested` exception | Deleted (`CreateAll` returns empty vector) |
| `NoDrRequested` exception | Deleted (`Create` returns `nullptr`) |
| `background_parameters_for_distributions` | Private detail in `NCDMBaseSpecies` |
| `decay_dr_map_` | Deleted (each `DNCDMSpecies` is its own instance) |
| `SourceType` enum | Deleted (only one DR source type per species class) |

## 10. Out of Scope

- Issue #246 (object-based input syntax): `CreateAll` is the exact seam where that change lands. No other code needs to change when #246 is implemented.
- Self-interacting NCDM branch: not merged to master, not considered here.
- IDM/DR species: untouched.
- Perturbation physics: unchanged throughout. This is a structural refactor only.
