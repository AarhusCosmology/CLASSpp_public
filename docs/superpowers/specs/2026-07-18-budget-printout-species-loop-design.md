# Budget printout: species-owned loop instead of per-type string branches

**Date:** 2026-07-18
**Status:** Approved
**Scope:** One PR. Resolves #377 (both defects).
**Issue:** #377 — "Budget printout: dcdm_wdm sector missing (TOTAL != 1), all
NCDM-family labeled 'Neutrino Species'".

## Problem

`BackgroundModule::background_output_budget()` (`source/background_module.cpp`,
guarded by `background_verbose > 1`) is the last holdout in the background module
still using the per-type "name-the-species" anti-pattern: a chain of ~10
`all_species_.count("CDM")` / `static_cast<const DCDM_DR_Species&>(...)` branches.
The rest of the module already migrated to the #308 convention (neutral accessors
looped over `all_species_`, e.g. `Omega0_dr_` is computed at
`background_module.cpp:693` by summing `DarkRadiationRhoToday()` over all species).

Two concrete defects fall out of the string-branch approach:

1. **The dcdm_wdm sector is absent, so TOTAL visibly misses ~1.** There is no
   `count("DCDM_WDM")` branch in the matter section, and `GetNcdmSpecies()`
   special-cases `DNCDM_DR_Species`' child but not `DCDM_WDM_Species`'. Running
   `test/dotsyntax_dcdm_wdm.ini` at `background_verbose 2` prints
   `TOTAL Omega = 0.974` — the missing ~0.026 is the dcdm_wdm sector.

2. **Every NCDM-family species prints the label "Neutrino Species"**
   (`NCDMBaseSpecies::PrintOmegaInfo`, `species/ncdm_base_species.cpp:637`), so
   grey-body dark radiation, thermal axions, and decay products are all called
   neutrinos.

A latent third symptom: Type3's (`CDM_SCF_Momentum`) `scf` child is omitted from
the budget entirely — only its `cdm` part is printed — so that sector does not
close today either.

## Key facts that shape the design

- **The budget already prints today densities.** `Omega0_dcdm_`/`Omega0_dr_` are
  integration-derived *today* values; `Omega0_g` and NCDM `Omega0_` are today
  values. Closure to 1 is defined by the today background-table row
  (`Omega0_de_ = 1 - Omega0_m_ - Omega0_r_ - Omega0_k`, `background_module.cpp:813`).
- **`Rho(pvecback)/pvecback[index_bg_rho_crit_]` *is* a species' Omega** (used at
  `background_module.cpp:794`). Reading the today row
  (`background_table_.data() + (bt_size_-1)*bg_size_`) and dividing each species'
  `Rho()` by `rho_crit` reproduces every current budget number exactly, with no
  integration-vector plumbing — and it is precisely the value that closes to 1.
  `background_output_budget()` runs after `background_solve_evolver()`
  (`background_module.cpp:478` then `:483`), so the today row is fully populated.
- **Every composite is `energy_type = Other`** (`dcdm_dr_species.cpp:14`, etc.),
  by design, so a composite cannot be bucketed as a single line — its children
  carry the real matter/radiation split. The composite must produce its own
  sector breakdown.
- **NCDM species are also `energy_type = Other`** (`ncdm_base_species.h:74`: they
  tally into both `rho_r` and `rho_m`), so `energy_type()` alone cannot separate
  NCDM from dark energy. `dynamic_cast<NCDMBaseSpecies*>` is the exact,
  already-blessed family test (the same one `GetNcdmSpecies` uses; the repo rule
  permits `HasNcdm`-style presence detection while forbidding concrete-type
  branches).
- **`energy_type` classification is complete:** Lambda/Fluid = `DarkEnergy`,
  ScalarField = `Other`, UR/IDR = `Radiation`, CDM/Baryons/DCDM = `Matter`,
  Photons = `Radiation`, NCDM-family/composites = `Other`.

## Design constraints (from review)

- **No new virtuals on `BaseSpecies`.** The vtable and species interface must not
  grow to serve a `verbose>1` printout.
- **No general `children()` accessor.** A composite is responsible for its whole
  sector; a raw child getter is a general-purpose hole that invites
  sector-bypassing use-cases. A *purpose-specific method on `CompositeSpecies`*
  is acceptable — children never leave the composite.

## Design

### 1. Shared POD types + free classifier (`species/base_species.h` / `.cpp`)

No vtable impact — plain types and one free function:

```cpp
enum class BudgetBucket { Radiation, NonRelativistic, Ncdm, Other };
struct BudgetLine { std::string label; double omega; BudgetBucket bucket; };

// Family-based classification; mirrors the GetNcdmSpecies idiom.
BudgetBucket BudgetBucketOf(const BaseSpecies& s);
```

`BudgetBucketOf` is defined in `base_species.cpp` (which includes
`ncdm_base_species.h` for the family cast):

```cpp
BudgetBucket BudgetBucketOf(const BaseSpecies& s) {
  if (dynamic_cast<const NCDMBaseSpecies*>(&s)) return BudgetBucket::Ncdm;
  switch (s.energy_type()) {
    case BaseSpecies::EnergyType::Radiation: return BudgetBucket::Radiation;
    case BaseSpecies::EnergyType::Matter:    return BudgetBucket::NonRelativistic;
    default:                                 return BudgetBucket::Other; // DarkEnergy/Other
  }
}
```

### 2. Sector-owned method on `CompositeSpecies` (`species/composite_species.h` / `.cpp`)

The composite produces its sector's lines; `children_` stays encapsulated:

```cpp
void AppendBudgetLines(const double* pvecback_today, double rho_crit,
                       std::vector<BudgetLine>& out) const {
  for (const auto& c : children_)
    out.push_back({c->name(), c->Rho(pvecback_today) / rho_crit, BudgetBucketOf(*c)});
}
```

### 3. Rewrite `background_output_budget()` (`source/background_module.cpp`)

Replace the ~10 `count(...)`/`static_cast` branches with a single loop plus a
local pretty-name map (defaulting to the terse `name()`), then print grouped by
bucket in a fixed section order with per-bucket subtotals and a grand total:

```cpp
const double* bg_today = background_table_.data() + (bt_size_ - 1) * bg_size_;
const double rho_crit  = bg_today[index_bg_rho_crit_];

// Display-only: terse internal names -> readable labels; unknown names pass through.
static const std::map<std::string, const char*> kPretty = {
    {"Baryons", "Baryons"}, {"Photons", "Photons"}, {"CDM", "Cold Dark Matter"},
    {"UR", "Ultra-relativistic relics"}, {"Lambda", "Cosmological Constant"},
    {"Fluid", "Dark Energy Fluid"}, {"ScalarField", "Scalar Field"},
    {"DCDM", "Decaying Cold Dark Matter"}, {"DR", "Dark Radiation"},
    {"IDR", "Interacting Dark Radiation"}, {"IDM_DR", "Interacting Dark Matter (DR)"}};
// (illustrative subset; the implementation seeds every known standard/sub-species
//  name. Any name absent from the map falls through to the terse name() verbatim.)

std::vector<BudgetLine> lines;
for (const auto& [key, sp] : all_species_) {
  if (auto* comp = dynamic_cast<const CompositeSpecies*>(sp.get()))
    comp->AppendBudgetLines(bg_today, rho_crit, lines);
  else
    lines.push_back({sp->name(), sp->Rho(bg_today) / rho_crit, BudgetBucketOf(*sp)});
}
if (pba->sgnK != 0)
  lines.push_back({"Spatial Curvature", pba->Omega0_k, BudgetBucket::Other});

// Prettify labels, then print sections [Radiation, NonRelativistic, Ncdm, Other]
// with subtotals, and a grand total (sum over all buckets).
```

The NCDM section header/total is renamed from "Massive Neutrino Species" /
"Neutrinos" to a neutral "Non-Cold Dark Matter (NCDM)" / "NCDM", since axions and
grey-body DR are not neutrinos.

### 4. Remove dead code (`species/ncdm_base_species.{h,cpp}`)

`PrintOmegaInfo`'s only caller is the budget loop (`background_module.cpp:1183`);
it is removed (declaration + definition). This deletes the `"Neutrino Species"`
hardcode that was defect 2.

`GetNcdmSpecies()` stays — it has other callers (`background_module.cpp:425, 430,
467, 694`). Concrete composite `#include`s in `background_module.cpp` that become
unused after the branch deletion are pruned; those still referenced elsewhere in
the file (e.g. `DCDM_DR_Species` at `:697`) are kept.

## How the design resolves each defect

- **Defect 1 (dcdm_wdm / TOTAL≠1):** the `DCDM_WDM` composite is now iterated via
  `AppendBudgetLines` → parent (`Matter` → Non-rel) + WDM daughter
  (`NCDMBaseSpecies` → NCDM) both appear from their today `Rho`; the sector
  density enters the grand total → **TOTAL → 1**.
- **Defect 2 (labels):** every line's label is the species' own `name()` (then
  optionally prettified), so the `"Neutrino Species"` hardcode is gone and each
  NCDM-family species prints its own identity.
- **Type3 latent gap:** iterating the sector now prints both `cdm` (Non-rel) and
  `scf` (Other) children, so that sector closes too.

## Numerical behaviour

`Rho(today)/rho_crit` equals every value the current printout uses
(`Omega0_dcdm_`, `Omega0_dr_`, `Omega0_g`, NCDM `Omega0_`, Lambda `GetOmega0`),
so existing single-sector cases print identical numbers. Cases previously broken
(dcdm_wdm, Type3) now close to 1. For `ScalarField`/`Fluid` the printout switches
from the "wished" `GetOmega0()` to the realized `Rho(today)/rho_crit`; these agree
to solver tolerance and the realized value is what actually closes the budget.

## Cost accounting

- 0 new virtuals.
- 1 purpose-specific method on `CompositeSpecies` (`AppendBudgetLines`).
- 1 free classifier (`BudgetBucketOf`) + 2 POD types.
- Deletes ~10 per-type string/`static_cast` branches and 1 dead method
  (`PrintOmegaInfo`).
- 2 `dynamic_cast`s, both to structural family bases (`CompositeSpecies`,
  `NCDMBaseSpecies`) and both already idiomatic in this file.

## Verification

Per the repo convention (no bit-identical requirement; budget is stdout-only, so
no output-golden impact):

1. Build (CMake + `make`), zero warnings on touched files.
2. `test/dotsyntax_dcdm_wdm.ini` at `background_verbose 2`: TOTAL ≈ 1.0, with
   dcdm_wdm parent + daughter lines present and correctly labeled.
3. A DCDM_DR case: "Decaying Cold Dark Matter" + "Dark Radiation" numbers match
   master to the printed precision (they read the same today densities).
4. A massive-ν / thermal-axion case: each NCDM-family species prints its own
   label (no "Neutrino Species"); NCDM subtotal unchanged.
5. A plain ΛCDM case: budget identical to master.
6. A Type3 (`CDM_SCF_Momentum`) case: both cdm and scf now appear; TOTAL ≈ 1.

## Non-goals

- Removing the remaining ~130 string-name branches elsewhere in the codebase
  (tracked separately; this PR only retires the budget-printout cluster).
- Changing any physics, output-file columns, or the Python wrapper.
- Per-channel refactor of `Omega0_dr_`/`Omega0_dcdm_` computation (already
  species-owned enough; untouched).
