# Remove species-specific physics params from pba — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move all 20 species-specific physics parameters off `pba` (the `background` struct) and onto the owning species. After this PR, `pba` retains only universal numerical state (H0, T_cmb, Omega0_g/b/k, h, K, sgnK, a_today, info flags, verbose). Closes issue #275.

**Architecture:** Each species gains private member fields for its physics params with per-field public const getters; its constructor takes one positional argument per field; `CreateAll(ctx)` parses every field from `ctx.pfc` and constructs the species. All cross-module readers migrate from `pba->X` to `bgm_->all_species_.find(...)` + accessor. Coupling-time reads (in `ReadCoupledOmegaBudget`) re-parse the keys inline — one duplicate `parser_read_double` per key, lookup cost is negligible.

**Tech Stack:** C++17 (species classes, BackgroundModule), Cython (classy.pyx — only `T_idr` needs a wrapper accessor), python pytest scenario harness for verification.

**Design spec:** `docs/superpowers/specs/2026-05-29-remove-species-physics-params-design.md` (head 02b49cce). Read before starting.

**Predecessor:** PR #276 (`275-remove-has-flags`, merged into master `b0494e56`) — established the `SpeciesCollection` / `SpeciesOmegaBudget` / `BackgroundModule::GetOmega0Species(key)` patterns this plan uses.

**Verification harness:** every commit is verified by running the 84-scenario test set; the final commit also runs the 1260-scenario set before pushing. See Task 0 for the exact commands. Diagnostic noise from clangd (e.g. `'common.h' file not found`, `Use of undeclared identifier '_TRUE_'`) is LSP false positives — ignore. Only trust `g++` errors from `make`.

**Working directory:** `/Users/au192734/Projects/class_claude/.claude/worktrees/remove-species-physics-params` (worktree on branch `worktree-remove-species-physics-params`, based on master HEAD `b0494e56`).

---

## Task 0: Verify baseline + survey

Confirm the worktree starts clean and the test harness works before any edits.

**Files:** none modified.

- [ ] **Step 1: Verify build clean**

```bash
make -j8 class 2>&1 | tail -3
```

Expected output ends with `g++ -O3 ... -o class` (final link line), no error markers.

- [ ] **Step 2: Install classy wrapper**

```bash
pip install . 2>&1 | tail -3
```

Expected: `Successfully installed classy-community-23.2.0`.

- [ ] **Step 3: Run baseline scenarios**

```bash
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -5
```

Expected: `84 passed, 98 deselected in ~40s`.

- [ ] **Step 4: No commit** — baseline check only.

---

## Task 1: Migrate Fluid (7 fields)

Move `fluid_equation_of_state`, `w0_fld`, `wa_fld`, `cs2_fld`, `Omega_EDE`, `use_ppf`, `c_gamma_over_c_fld` from `pba` onto `FluidSpecies`.

**Files:**
- Modify: `species/fluid.h` — add private fields + getters + expanded constructor
- Modify: `species/fluid.cpp` — implement expanded constructor, replace `pba_.X` with `X_` in self-reads, update `CreateAll` parsing
- Modify: `source/input_module.cpp` — delete fluid parser block lines 1055–1097
- Modify: `source/perturbations_module.cpp` — 5 cross-module reads → species lookup
- Modify: `source/background_module.cpp` — 4 cross-module reads → species lookup
- Modify: `source/nonlinear_module.cpp` — 1 halofit `wa_fld != 0.` check → species lookup
- Modify: `source/background.h` — delete the 7 fluid fields from struct

- [ ] **Step 1: Update `species/fluid.h`**

Add private members (with defaults matching today's pba), public getters, expanded constructor signature. Replace the existing single-Omega constructor declaration:

```cpp
// In class FluidSpecies, public section:
FluidSpecies(const background& pba,
             double omega0_fld,
             equation_of_state fluid_eos,
             double w0_fld,
             double wa_fld,
             double cs2_fld,
             double Omega_EDE,
             short  use_ppf,
             double c_gamma_over_c_fld);

double GetOmega0() const override { return Omega0_fld_; }
equation_of_state fluid_eos() const { return fluid_eos_; }
double w0_fld() const { return w0_fld_; }
double wa_fld() const { return wa_fld_; }
double cs2_fld() const { return cs2_fld_; }
double Omega_EDE() const { return Omega_EDE_; }
short  use_ppf() const { return use_ppf_; }
double c_gamma_over_c_fld() const { return c_gamma_over_c_fld_; }
```

In the private section (after existing `Omega0_fld_`):

```cpp
equation_of_state fluid_eos_ = CLP;
double w0_fld_              = -1.;
double wa_fld_              = 0.;
double cs2_fld_             = 1.;
double Omega_EDE_           = 0.;
short  use_ppf_             = _TRUE_;
double c_gamma_over_c_fld_  = 0.4;
```

- [ ] **Step 2: Update `species/fluid.cpp` constructor**

Replace the single-arg implementation:

```cpp
FluidSpecies::FluidSpecies(const background& pba,
                           double omega0_fld,
                           equation_of_state fluid_eos,
                           double w0_fld,
                           double wa_fld,
                           double cs2_fld,
                           double Omega_EDE,
                           short  use_ppf,
                           double c_gamma_over_c_fld)
    : BaseSpecies("Fluid", EnergyType::DarkEnergy),
      pba_(pba),
      Omega0_fld_(omega0_fld),
      fluid_eos_(fluid_eos),
      w0_fld_(w0_fld),
      wa_fld_(wa_fld),
      cs2_fld_(cs2_fld),
      Omega_EDE_(Omega_EDE),
      use_ppf_(use_ppf),
      c_gamma_over_c_fld_(c_gamma_over_c_fld) {}
```

- [ ] **Step 3: Replace all `pba_.X` reads in `species/fluid.cpp` with member-field reads**

In `species/fluid.cpp`, every site that today reads `pba_.fluid_equation_of_state`, `pba_.w0_fld`, `pba_.wa_fld`, `pba_.cs2_fld`, `pba_.Omega_EDE`, `pba_.use_ppf`, or `pba_.c_gamma_over_c_fld` becomes the member access (`fluid_eos_`, `w0_fld_`, etc.). Search:

```bash
grep -n "pba_\.\(fluid_equation_of_state\|w0_fld\|wa_fld\|cs2_fld\|Omega_EDE\|use_ppf\|c_gamma_over_c_fld\)" species/fluid.cpp
```

Replace each occurrence with the corresponding `<name>_` member. The EDE branch in `ComputeWFld` is the heaviest hitter.

- [ ] **Step 4: Update `species/fluid.cpp` `CreateAll`**

Current `CreateAll` parses `Omega_fld` only. Extend it to also parse the 7 physics fields, then construct with them. The parsing logic mirrors `input_module.cpp:1055–1097` (which we delete in Step 7). New body:

```cpp
std::vector<Named> FluidSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  // Existing Omega_fld parse + closure-override resolution — keep it.
  double omega0_fld = 0.;
  if (ctx.omega0_closure_override) {
    omega0_fld = *ctx.omega0_closure_override;
  } else {
    ctx.pfc->read_double("Omega_fld", omega0_fld);
  }
  if (omega0_fld == 0.)
    return {};

  // ── fluid_equation_of_state (string-keyed enum) ────────────────────────────
  equation_of_state fluid_eos = CLP;
  std::string eos_str;
  if (ctx.pfc->read_string("fluid_equation_of_state", eos_str)) {
    if (eos_str.find("CLP") != std::string::npos ||
        eos_str.find("clp") != std::string::npos) {
      fluid_eos = CLP;
    } else if (eos_str.find("EDE") != std::string::npos ||
               eos_str.find("ede") != std::string::npos) {
      fluid_eos = EDE;
    } else {
      throw std::invalid_argument(
          "incomprehensible input '" + eos_str +
          "' for the field 'fluid_equation_of_state'");
    }
  }

  // ── numeric params ────────────────────────────────────────────────────────
  double w0_fld    = -1.;
  double wa_fld    = 0.;
  double cs2_fld   = 1.;
  double Omega_EDE = 0.;
  if (fluid_eos == CLP) {
    ctx.pfc->read_double("w0_fld",  w0_fld);
    ctx.pfc->read_double("wa_fld",  wa_fld);
    ctx.pfc->read_double("cs2_fld", cs2_fld);
  } else {  // EDE
    ctx.pfc->read_double("w0_fld",    w0_fld);
    ctx.pfc->read_double("Omega_EDE", Omega_EDE);
    ctx.pfc->read_double("cs2_fld",   cs2_fld);
  }

  // ── PPF flag + sound-speed param ──────────────────────────────────────────
  short  use_ppf            = _TRUE_;
  double c_gamma_over_c_fld = 0.4;
  std::string ppf_str;
  if (ctx.pfc->read_string("use_ppf", ppf_str)) {
    use_ppf = (ppf_str.find("y") != std::string::npos ||
               ppf_str.find("Y") != std::string::npos) ? _TRUE_ : _FALSE_;
  }
  if (use_ppf == _TRUE_)
    ctx.pfc->read_double("c_gamma_over_c_fld", c_gamma_over_c_fld);

  std::vector<Named> result;
  result.push_back({"Fluid",
                    std::make_unique<FluidSpecies>(*ctx.pba, omega0_fld,
                                                   fluid_eos, w0_fld, wa_fld,
                                                   cs2_fld, Omega_EDE,
                                                   use_ppf, c_gamma_over_c_fld)});
  return result;
}
```

Add `#include <stdexcept>` and `#include <string>` at the top of `species/fluid.cpp` if not already present.

- [ ] **Step 5: Migrate cross-module readers**

5a. `source/perturbations_module.cpp` — 5 reads. Find them:

```bash
grep -n "pba->\(fluid_equation_of_state\|w0_fld\|wa_fld\|cs2_fld\|Omega_EDE\|use_ppf\|c_gamma_over_c_fld\)" source/perturbations_module.cpp
```

Replace each with the same pattern PR #276 used:

```cpp
const FluidSpecies* fluid = nullptr;
if (auto* p = all_species_.find("Fluid"))
  fluid = static_cast<const FluidSpecies*>(p->get());
// then: fluid ? fluid->w0_fld() : <default>
```

Define `fluid` once per function scope (locally) and use its accessor. If a function reads multiple fluid fields, hoist the lookup once at the top. `#include "../species/fluid.h"` if not already there.

5b. `source/background_module.cpp` — 4 reads. Same pattern. Both `background_functions` and the verbose print sites are touched.

5c. `source/nonlinear_module.cpp` — line ~3438 (`pba->wa_fld != 0.` Pk_equal guard). Same pattern: lookup Fluid via `all_species_.find("Fluid")` and read `wa_fld()`.

5d. `source/input_module.cpp` — line 2793 (`pba->wa_fld != 0.` halofit check, inside the same function where `fluid_present_pfc` is in scope). Replace with `parser_read_double("wa_fld", ...)` peek (since the species isn't built yet at this point in input_init), OR, simpler: move the halofit gate down past `ConstructSpecies` and do the species lookup there. Pick whichever is structurally cleaner; both preserve behavior.

- [ ] **Step 6: Delete fluid parser block in `source/input_module.cpp`**

Delete lines 1055–1097 (the block that parses `use_ppf`, `fluid_equation_of_state`, and the CLP/EDE branches into `pba->X`). Leave the `Omega_fld` parser_read_double (it's used by the closure decision and survives — see PR #276 design).

Actually those parses run in the same `if (fluid_present_pfc) { … }` block; delete the body, leave the `if` and replace with a comment:

```cpp
if (fluid_present_pfc) {
  // Fluid physics params (use_ppf, fluid_equation_of_state, w0_fld, wa_fld,
  // cs2_fld, Omega_EDE, c_gamma_over_c_fld) are parsed inside
  // FluidSpecies::CreateAll directly from pfc; no per-key writes to pba here.
}
```

Or just delete the whole `if (fluid_present_pfc) { … }` block if it's now empty after Step 6.

- [ ] **Step 7: Delete fluid fields from `source/background.h`**

Delete lines 76–94 (the `fluid_equation_of_state` enum + `w0_fld` + `wa_fld` + `Omega_EDE` + `cs2_fld` + `use_ppf` + `c_gamma_over_c_fld` declarations + their interleaved comments).

Leave the `enum equation_of_state { CLP, EDE };` declaration at line 21 in place (it's used by `species/fluid.h`'s new constructor signature and the Cython wrapper enums).

- [ ] **Step 8: Build**

```bash
make -j8 class 2>&1 | grep -E "error:" | head -20
```

Expected: no `error:` lines. If errors remain, they're consumer reads we missed in Step 5 — grep them out:

```bash
grep -rnE "pba(_)?(->|\.)w0_fld|pba(_)?(->|\.)wa_fld|pba(_)?(->|\.)cs2_fld|pba(_)?(->|\.)Omega_EDE|pba(_)?(->|\.)use_ppf|pba(_)?(->|\.)c_gamma_over_c_fld|pba(_)?(->|\.)fluid_equation_of_state" --include='*.cpp' --include='*.h' source species
```

Fix each and rebuild.

- [ ] **Step 9: Reinstall + run scenarios**

```bash
pip install . 2>&1 | tail -3
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -5
```

Expected: `84 passed`.

- [ ] **Step 10: Commit**

```bash
git add -- species/fluid.h species/fluid.cpp source/input_module.cpp \
            source/perturbations_module.cpp source/background_module.cpp \
            source/nonlinear_module.cpp source/background.h
git commit -m "$(cat <<'EOF'
fluid: own fluid_eos / w0 / wa / cs2 / Omega_EDE / use_ppf / c_gamma

FluidSpecies stores the 7 fluid physics params as private members with
per-field getters; CreateAll parses them from pfc and feeds the new
9-arg constructor.  Removes pba->fluid_equation_of_state /w0_fld /wa_fld
/cs2_fld /Omega_EDE /use_ppf /c_gamma_over_c_fld and the matching
parser block in input_module.cpp.  Cross-module readers (perturbations,
background, nonlinear, input halofit gate) migrate to bgm/input lookups
of the Fluid species.

84/84 scenarios pass at TEST_LEVEL=1.
EOF
)"
```

---

## Task 2: Migrate ScalarField (5 fields)

Move `attractor_ic_scf`, `phi_ini_scf`, `phi_prime_ini_scf`, `scf_parameters` (vector), `scf_tuning_index` from `pba` onto `ScalarFieldSpecies`.

**Files:**
- Modify: `species/scalar_field.h` — add private fields + getters + expanded constructor
- Modify: `species/scalar_field.cpp` — constructor body + replace `pba_.X` self-reads + extend `CreateAll`
- Modify: `source/input_module.cpp` — delete scf parser block lines 1100–1146
- Modify: `source/background_module.cpp` — 36 cross-module reads → species lookup (V_scf eval, IC seed, verbose print)
- Modify: `source/background.h` — delete the 5 scf fields

Same step-by-step shape as Task 1.

- [ ] **Step 1: Update `species/scalar_field.h`**

Public:

```cpp
ScalarFieldSpecies(const background& pba,
                   double omega0_scf,
                   std::vector<double> scf_parameters,
                   int    scf_tuning_index,
                   short  attractor_ic_scf,
                   double phi_ini_scf,
                   double phi_prime_ini_scf);

double GetOmega0() const override { return Omega0_scf_; }
const std::vector<double>& scf_parameters() const { return scf_parameters_; }
int    scf_tuning_index() const { return scf_tuning_index_; }
short  attractor_ic_scf() const { return attractor_ic_scf_; }
double phi_ini_scf() const { return phi_ini_scf_; }
double phi_prime_ini_scf() const { return phi_prime_ini_scf_; }
```

Private (after existing `Omega0_scf_`):

```cpp
std::vector<double> scf_parameters_;
int    scf_tuning_index_   = 0;
short  attractor_ic_scf_   = _TRUE_;
double phi_ini_scf_        = 1.;
double phi_prime_ini_scf_  = 1.;
```

Add `#include <vector>` if not already.

- [ ] **Step 2: Update `species/scalar_field.cpp` constructor**

```cpp
ScalarFieldSpecies::ScalarFieldSpecies(const background& pba,
                                       double omega0_scf,
                                       std::vector<double> scf_parameters,
                                       int    scf_tuning_index,
                                       short  attractor_ic_scf,
                                       double phi_ini_scf,
                                       double phi_prime_ini_scf)
    : BaseSpecies("ScalarField", EnergyType::Other),
      pba_(pba),
      Omega0_scf_(omega0_scf),
      scf_parameters_(std::move(scf_parameters)),
      scf_tuning_index_(scf_tuning_index),
      attractor_ic_scf_(attractor_ic_scf),
      phi_ini_scf_(phi_ini_scf),
      phi_prime_ini_scf_(phi_prime_ini_scf) {}
```

- [ ] **Step 3: Replace `pba_.scf_parameters` etc. in `species/scalar_field.cpp`**

```bash
grep -n "pba_\.\(scf_parameters\|scf_tuning_index\|attractor_ic_scf\|phi_ini_scf\|phi_prime_ini_scf\)" species/scalar_field.cpp
```

`pba_.scf_parameters` → `scf_parameters_`; `pba_.scf_parameters[0]` etc. → `scf_parameters_[0]`. Same for the rest.

- [ ] **Step 4: Update `species/scalar_field.cpp` `CreateAll`**

The current `CreateAll` parses `Omega_scf` + handles closure + has a shooter target path. Extend it to also parse the 5 physics fields. Code copied from `input_module.cpp:1100–1146`:

```cpp
// Inside CreateAll, after Omega_scf has been resolved (current code).
// If omega0_scf == 0., return empty (current code).

// ── scf_parameters (variable-length list) ────────────────────────────────
std::vector<double> scf_parameters;
{
  int flag = _FALSE_;
  // readDoubleList lives in input_module.h — replicate its semantics inline
  // using ctx.pfc directly:
  std::string raw;
  if (ctx.pfc->read_string("scf_parameters", raw)) {
    // parse the space-separated double list; readDoubleList already does
    // this — call it (it's static-free) or inline its logic. Easiest:
    // forward-declare extern: `int readDoubleList(FileContent*, const char*,
    // std::vector<double>&, int*, ErrorMsg);` and call it.
    char errmsg[2048] = {0};
    if (readDoubleList(ctx.pfc, "scf_parameters", scf_parameters, &flag, errmsg) != _SUCCESS_)
      throw std::invalid_argument(std::string("scf_parameters parse error: ") + errmsg);
  }
}

int scf_tuning_index = 0;
ctx.pfc->read_int("scf_tuning_index", scf_tuning_index);
if (scf_tuning_index >= static_cast<int>(scf_parameters.size()))
  throw std::invalid_argument(
      "Tuning index scf_tuning_index = " + std::to_string(scf_tuning_index) +
      " is larger than the number of entries " +
      std::to_string(scf_parameters.size()) + " in scf_parameters.");

// scf_shooting_parameter overwrites the slot at tuning_index
double scf_shooting_param = 0.;
if (ctx.pfc->read_double("scf_shooting_parameter", scf_shooting_param))
  scf_parameters[scf_tuning_index] = scf_shooting_param;

// ── attractor_ic_scf (y/n string) ────────────────────────────────────────
short  attractor_ic_scf   = _TRUE_;
double phi_ini_scf        = 1.;
double phi_prime_ini_scf  = 1.;
std::string attr_str;
if (ctx.pfc->read_string("attractor_ic_scf", attr_str)) {
  if (attr_str.find("y") != std::string::npos ||
      attr_str.find("Y") != std::string::npos) {
    attractor_ic_scf = _TRUE_;
  } else {
    attractor_ic_scf = _FALSE_;
    if (scf_parameters.size() < 2)
      throw std::invalid_argument(
          "Since you are not using attractor initial conditions, you must "
          "specify phi and its derivative phi' as the last two entries in "
          "scf_parameters. See explanatory.ini for more details.");
    phi_ini_scf       = scf_parameters[scf_parameters.size() - 2];
    phi_prime_ini_scf = scf_parameters[scf_parameters.size() - 1];
  }
}

// Construct with the parsed params.
result.push_back({"ScalarField",
                  std::make_unique<ScalarFieldSpecies>(*ctx.pba, omega0_scf,
                                                       std::move(scf_parameters),
                                                       scf_tuning_index,
                                                       attractor_ic_scf,
                                                       phi_ini_scf,
                                                       phi_prime_ini_scf)});
```

`readDoubleList` is declared in `source/input_module.h`. Add `#include "../source/input_module.h"` to `species/scalar_field.cpp` (it's already pulled in indirectly via other headers in the species file; verify with the build).

Also add `#include "parser.h"` (already there from PR #276's Phase 2B; verify).

- [ ] **Step 5: Migrate `source/background_module.cpp` SCF readers**

36 sites. Hoist a single lookup at function scope:

```cpp
const ScalarFieldSpecies* scf = nullptr;
if (auto* p = all_species_.find("ScalarField"))
  scf = static_cast<const ScalarFieldSpecies*>(p->get());
```

Then replace:
- `pba->scf_parameters` → `scf->scf_parameters()`
- `pba->scf_parameters[i]` → `scf->scf_parameters()[i]`
- `pba->attractor_ic_scf` → `scf->attractor_ic_scf()`
- `pba->phi_ini_scf` → `scf->phi_ini_scf()`
- `pba->phi_prime_ini_scf` → `scf->phi_prime_ini_scf()`

Search to confirm coverage:

```bash
grep -n "pba->\(scf_parameters\|scf_tuning_index\|attractor_ic_scf\|phi_ini_scf\|phi_prime_ini_scf\)" source/background_module.cpp
```

Add `#include "../species/scalar_field.h"` if not present.

- [ ] **Step 6: Delete the scf parser block in `source/input_module.cpp`**

Delete lines 1100–1146 (the `if (scf_present_pfc) { … }` body that does `readDoubleList`, sets `pba->scf_parameters`, `pba->scf_tuning_index`, the `scf_shooting_parameter` override, and the attractor IC branch).

Leave the `Omega_scf` parser_read_double at line 1003 (closure decision still needs it).

- [ ] **Step 7: Delete scf fields from `source/background.h`**

Delete lines 105–111 (`attractor_ic_scf`, `phi_ini_scf`, `phi_prime_ini_scf`, `scf_parameters`, `scf_tuning_index`).

- [ ] **Step 8: Build + test + commit**

Same shape as Task 1 Steps 8–10. Commit message:

```
scalar_field: own scf_parameters / phi_ini / phi_prime_ini / attractor_ic / scf_tuning_index

ScalarFieldSpecies stores the 5 SCF physics params as private members
with per-field getters; CreateAll parses them from pfc (including the
scf_parameters list, tuning index validation, scf_shooting_parameter
override, and attractor-IC y/n branch).  background_module's 36 reads
(V_scf eval, BackgroundDerivs, verbose print) migrate to species lookup.
Deletes the matching pba fields and the scf parser block in input_module.

84/84 scenarios pass at TEST_LEVEL=1.
```

---

## Task 3: Migrate DCDM (Gamma_dcdm, Omega_ini_dcdm)

Move `Gamma_dcdm` and `Omega_ini_dcdm` from `pba` onto `DCDMSpecies` (the DCDM sub-species inside `DCDM_DR_Species`). The composite forwards both to the sub-species at construction.

**Files:**
- Modify: `species/dcdm.h` — add private fields + getters + expanded constructor
- Modify: `species/dcdm.cpp` — constructor body + replace `pba_.Gamma_dcdm` / `pba_.Omega_ini_dcdm` self-reads
- Modify: `species/dcdm_dr_species.h` — composite constructor takes Gamma_dcdm + Omega_ini_dcdm, forwards
- Modify: `species/dcdm_dr_species.cpp` — composite constructor forwards; replace `pba_->Gamma_dcdm` reads with `dcdm_->Gamma_dcdm()`; update `CreateAll` to parse and pass through
- Modify: `source/input_module.cpp` — delete the Gamma_dcdm + Omega_ini_dcdm reads inside the DCDM parser block in `ReadCoupledOmegaBudget`; only the inline parses for budget-time use remain
- Modify: `source/background_module.cpp` — 2 reads → species lookup
- Modify: `source/background.h` — delete `Gamma_dcdm` and `Omega_ini_dcdm` fields

- [ ] **Step 1: Update `species/dcdm.h`**

```cpp
DCDMSpecies(const background& pba,
            double omega0_dcdmdr,
            double Gamma_dcdm,
            double Omega_ini_dcdm);

double Gamma_dcdm() const { return Gamma_dcdm_; }
double Omega_ini_dcdm() const { return Omega_ini_dcdm_; }
```

Private:

```cpp
double Gamma_dcdm_     = 0.;
double Omega_ini_dcdm_ = 0.;
```

- [ ] **Step 2: Update `species/dcdm.cpp` constructor body + self-reads**

```cpp
DCDMSpecies::DCDMSpecies(const background& pba,
                         double omega0_dcdmdr,
                         double Gamma_dcdm,
                         double Omega_ini_dcdm)
    : BaseSpecies("DCDM", EnergyType::Matter),
      pba_(pba),
      Omega0_dcdmdr_(omega0_dcdmdr),
      Gamma_dcdm_(Gamma_dcdm),
      Omega_ini_dcdm_(Omega_ini_dcdm) {}
```

Replace 5 self-reads in `species/dcdm.cpp` (`pba_.Omega_ini_dcdm` at line 21, `pba_.Gamma_dcdm` at lines 37, 54, 82, 208 in the current file — re-grep to confirm) with `Omega_ini_dcdm_` / `Gamma_dcdm_`.

- [ ] **Step 3: Update `species/dcdm_dr_species.h` constructor**

```cpp
DCDM_DR_Species(const background* pba,
                const BackgroundModule* bgm,
                double omega0_dcdmdr,
                double Gamma_dcdm,
                double Omega_ini_dcdm);
```

- [ ] **Step 4: Update `species/dcdm_dr_species.cpp` constructor body**

```cpp
DCDM_DR_Species::DCDM_DR_Species(const background* pba,
                                 const BackgroundModule* bgm,
                                 double omega0_dcdmdr,
                                 double Gamma_dcdm,
                                 double Omega_ini_dcdm)
    : CompositeSpecies("DCDM_DR", BaseSpecies::EnergyType::Other),
      pba_(pba),
      bgm_(bgm) {
  auto dcdm  = std::make_unique<DCDMSpecies>(*pba, omega0_dcdmdr,
                                              Gamma_dcdm, Omega_ini_dcdm);
  auto dr_sp = std::make_unique<DarkRadiationSpecies>("DR", pba, bgm);
  dcdm_      = dcdm.get();
  dr_sp_     = dr_sp.get();
  children_.push_back(std::move(dcdm));
  children_.push_back(std::move(dr_sp));
}
```

Replace 4 `pba_->Gamma_dcdm` reads in `species/dcdm_dr_species.cpp` (lines 46, 60, 102, 206, 284 — re-grep) with `dcdm_->Gamma_dcdm()`.

- [ ] **Step 5: Update `species/dcdm_dr_species.cpp` `CreateAll`**

Today's `CreateAll` reads `Omega_dcdmdr` from `ctx.omega_budget` and constructs the composite. Add parses for `Gamma_dcdm` and `Omega_ini_dcdm`:

```cpp
// Inside CreateAll, after existence determination:
double Gamma_dcdm = 0.;
if (ctx.pfc->read_double("Gamma_dcdm", Gamma_dcdm))
  Gamma_dcdm *= 1.e3 / _c_;  // convert km/s/Mpc → Mpc^-1 (same as input_module)

double Omega_ini_dcdm = 0.;
double omega_ini_dcdm = 0.;
if (ctx.pfc->read_double("Omega_ini_dcdm", Omega_ini_dcdm)) {
  /* use as-is */
} else if (ctx.pfc->read_double("omega_ini_dcdm", omega_ini_dcdm)) {
  Omega_ini_dcdm = omega_ini_dcdm / (ctx.pba->h * ctx.pba->h);
}

// pass to constructor
result.push_back({"DCDM_DR",
                  std::make_unique<DCDM_DR_Species>(ctx.pba, ctx.bgm,
                                                    omega0_dcdmdr,
                                                    Gamma_dcdm,
                                                    Omega_ini_dcdm)});
```

The shooter-discovery path (around line 436 in the current file, see PR #276) constructs a `tmp_composite` — update both call sites.

Add `#include "common.h"` if `_c_` isn't already visible.

- [ ] **Step 6: Delete DCDM physics writes in `source/input_module.cpp` ReadCoupledOmegaBudget**

In the DCDM block of `ReadCoupledOmegaBudget` (around lines 988–1010 in current `input_module.cpp`), delete the `pba->Omega_ini_dcdm = …` and `pba->Gamma_dcdm = …` writes (parsing for the budget computation stays as local reads — the budget only needs Omega_dcdmdr's value, not Gamma or ini, but verify by reading the block).

If the budget does need Gamma_dcdm for the shooter guess (it does — `ComputeShootingGuess` uses `ba.Gamma_dcdm`), parse it to a local variable and pass via context or accept the duplicate parse in the species's `CreateAll`. Cleanest: in `ReadCoupledOmegaBudget`, parse `Gamma_dcdm` into a local for the budget's `omega_budget_.dcdmdr` resolution, then DCDM_DR's `CreateAll` re-parses it for storage. Per the design, duplicate parse is fine.

- [ ] **Step 7: Migrate `source/background_module.cpp` DCDM readers**

```bash
grep -n "pba->\(Gamma_dcdm\|Omega_ini_dcdm\)" source/background_module.cpp
```

Hoist a DCDM_DR lookup once at function scope, read via `comp.dcdm().Gamma_dcdm()` etc.:

```cpp
const DCDM_DR_Species* dcdm_dr = nullptr;
if (auto* p = all_species_.find("DCDM_DR"))
  dcdm_dr = static_cast<const DCDM_DR_Species*>(p->get());
// then: dcdm_dr ? dcdm_dr->dcdm().Gamma_dcdm() : 0.
```

Also update DCDM_DR_Species's ComputeShootingGuess to read Gamma_dcdm from `dcdm_->Gamma_dcdm()` instead of `ba.Gamma_dcdm`.

- [ ] **Step 8: Delete `Gamma_dcdm` and `Omega_ini_dcdm` from `source/background.h`**

Delete lines 99 (Gamma_dcdm) and 103 (Omega_ini_dcdm).

- [ ] **Step 9: Build + test + commit**

Same shape. Commit message:

```
dcdm: DCDMSpecies owns Gamma_dcdm and Omega_ini_dcdm

The DCDM sub-species stores its two physics params; the DCDM_DR composite
constructor takes them and forwards.  DCDM_DR's ComputeShootingGuess and
SetBackgroundInitialConditions read via dcdm_->Gamma_dcdm() etc.
ReadCoupledOmegaBudget parses Gamma_dcdm inline for the shooter-guess
matter sum; DCDM_DR::CreateAll re-parses for storage (one duplicate
map lookup, negligible).

84/84 scenarios pass at TEST_LEVEL=1.
```

---

## Task 4: Migrate IDR (T_idr, l_max_idr)

Move `T_idr` and `l_max_idr` from `pba` onto `IDRSpecies` (sub-species of `IDM_DR_IDR_Species`). The composite forwards both to IDR's constructor.

**Files:**
- Modify: `species/idr.h` — add private fields + getters + expanded constructor
- Modify: `species/interacting_species.cpp` — replace `pba_.l_max_idr` with `l_max_idr_`
- Modify: `species/idm_dr_idr_species.h` — composite constructor takes T_idr + l_max_idr
- Modify: `species/idm_dr_idr_species.cpp` — forward to IDRSpecies, update `CreateAll`
- Modify: `source/input_module.cpp` — delete `pba->T_idr = …` writes in `ReadCoupledOmegaBudget`; delete `pba->l_max_idr = ppr->l_max_idr` write in `input_read_parameters`
- Modify: `source/thermodynamics_module.cpp` — 10 cross-module reads → composite/IDR accessor
- Modify: `source/background.h` — delete `T_idr` and `l_max_idr`

- [ ] **Step 1: Update `species/idr.h`**

Expand the existing constructor (PR #276 added `has_sibling_idm_dr`):

```cpp
IDRSpecies(const background& pba,
           double omega0_idr,
           bool has_sibling_idm_dr,
           double T_idr,
           int l_max_idr)
    : BaseSpecies("IDR", EnergyType::Radiation),
      pba_(pba),
      Omega0_idr_(omega0_idr),
      has_sibling_idm_dr_(has_sibling_idm_dr),
      T_idr_(T_idr),
      l_max_idr_(l_max_idr) {}

double T_idr() const { return T_idr_; }
int    l_max_idr() const { return l_max_idr_; }
```

Private:

```cpp
double T_idr_     = 0.;
int    l_max_idr_ = 0;
```

- [ ] **Step 2: Replace `pba_.l_max_idr` reads in `species/interacting_species.cpp`**

Three reads at lines ~140, 142, 144 inside `IDRSpecies::RegisterPerturbationIndices` — replace with `l_max_idr_`.

- [ ] **Step 3: Update `species/idm_dr_idr_species.h` + .cpp composite constructor**

```cpp
// h
IDM_DR_IDR_Species(const background& pba,
                   double omega0_idm_dr,
                   double omega0_idr,
                   double T_idr,
                   int    l_max_idr);

// cpp
IDM_DR_IDR_Species::IDM_DR_IDR_Species(const background& pba,
                                       double omega0_idm_dr,
                                       double omega0_idr,
                                       double T_idr,
                                       int    l_max_idr)
    : CompositeSpecies("IDM_DR_IDR", BaseSpecies::EnergyType::Other), pba_(pba) {
  has_idm_dr_ = (omega0_idm_dr != 0.);
  has_idr_    = (omega0_idr != 0.);
  auto idm = std::make_unique<IDM_DRSpecies>(pba, omega0_idm_dr);
  auto idr = std::make_unique<IDRSpecies>(pba, omega0_idr, has_idm_dr_,
                                          T_idr, l_max_idr);
  idm_dr_  = idm.get();
  idr_     = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}
```

- [ ] **Step 4: Update `species/idm_dr_idr_species.cpp` `CreateAll`**

Currently `CreateAll` reads Omega values from the budget. Add T_idr derivation + l_max_idr parsing. Note: the T_idr value comes from `N_idr`/`N_dg`/`xi_idr` parsing (same logic that lives in `ReadCoupledOmegaBudget`). Duplicate the parse here (one-off cost per call):

```cpp
// Inside CreateAll, after existence determined:
double T_idr_local = 0.;
double xi_idr;
double N_idr;
double N_dg;
double stat_f_idr = 7. / 8.;
ctx.pfc->read_double("stat_f_idr", stat_f_idr);

bool flag_N_idr = ctx.pfc->read_double("N_idr",  N_idr);
bool flag_N_dg  = ctx.pfc->read_double("N_dg",   N_dg);
bool flag_xi    = ctx.pfc->read_double("xi_idr", xi_idr);

if (flag_N_idr)
  T_idr_local = std::pow(N_idr / stat_f_idr * (7./8.) / std::pow(11./4., 4./3.),
                         1./4.) * ctx.pba->T_cmb;
else if (flag_N_dg)
  T_idr_local = std::pow(N_dg / stat_f_idr * (7./8.) / std::pow(11./4., 4./3.),
                         1./4.) * ctx.pba->T_cmb;
else if (flag_xi)
  T_idr_local = xi_idr * ctx.pba->T_cmb;

int l_max_idr_local = ctx.ppr->l_max_idr;

// Pass to constructor:
result.push_back({"IDM_DR_IDR",
                  std::make_unique<IDM_DR_IDR_Species>(*ctx.pba,
                                                       omega0_idm_dr,
                                                       omega0_idr,
                                                       T_idr_local,
                                                       l_max_idr_local)});
```

- [ ] **Step 5: Delete `pba->T_idr = …` and `pba->l_max_idr = …` writes**

In `source/input_module.cpp`:
- `ReadCoupledOmegaBudget`: replace `pba->T_idr = …` with a local `T_idr_local` used only to compute `omega_budget_.idr`. Delete the writes.
- `input_read_parameters`: delete `pba->l_max_idr = ppr->l_max_idr;`.

- [ ] **Step 6: Migrate thermodynamics_module.cpp**

10 `pba->T_idr` reads at lines 241, 244, 882, 884, 891, 924, 939, 950, 990, 1019 (re-grep). Hoist once per function scope:

```cpp
const IDM_DR_IDR_Species* comp = nullptr;
if (auto* p = all_species_.find("IDM_DR_IDR"))
  comp = static_cast<const IDM_DR_IDR_Species*>(p->get());
const double T_idr = comp ? comp->idr().T_idr() : 0.;
```

Then replace `pba->T_idr` with the local `T_idr`.

`#include "../species/idm_dr_idr_species.h"` if not already.

- [ ] **Step 7: Delete from `source/background.h`**

Delete lines 96–101 region (`T_idr` and `l_max_idr` declarations).

- [ ] **Step 8: Build + test + commit**

```
idr: IDRSpecies owns T_idr and l_max_idr

IDRSpecies stores both as private fields; the IDM_DR_IDR composite
constructor takes them and forwards.  IDM_DR_IDR::CreateAll re-derives
T_idr from N_idr/N_dg/xi_idr inputs (duplicating the budget step's parse,
negligible cost).  IDRSpecies::RegisterPerturbationIndices reads its own
l_max_idr_.  Cross-module readers (thermodynamics_module, 10 sites for
the IDM_DR coupling math) migrate to comp.idr().T_idr() via all_species_.

84/84 scenarios pass at TEST_LEVEL=1.
```

---

## Task 5: Migrate IDM_DRMD (f_idm_drmd, G_over_aH_drmd, delta_Neff_drmd, z_stop)

Move the 4 DRMD coupling fields from `pba` onto `IDM_DRMD_IDR_DRMD_Species` (the composite, since these are coupling params governing the IDM_DRMD ↔ IDR_DRMD relationship).

**Files:**
- Modify: `species/idm_drmd_idr_drmd_species.h` — composite constructor takes 4 new args + getters
- Modify: `species/idm_drmd_idr_drmd_species.cpp` — constructor body + update `CreateAll`
- Modify: `source/input_module.cpp` — delete the 4 `pba->X_drmd = …` writes in `ReadCoupledOmegaBudget`
- Modify: `source/background_module.cpp` — 6 cross-module reads → composite accessor
- Modify: `source/background.h` — delete the 4 DRMD fields

- [ ] **Step 1: Update `species/idm_drmd_idr_drmd_species.h`**

```cpp
IDM_DRMD_IDR_DRMD_Species(const background& pba,
                          double omega0_idm_drmd,
                          double omega0_idr_drmd,
                          double f_idm_drmd,
                          double G_over_aH_drmd,
                          double delta_Neff_drmd,
                          double z_stop);

double f_idm_drmd() const { return f_idm_drmd_; }
double G_over_aH_drmd() const { return G_over_aH_drmd_; }
double delta_Neff_drmd() const { return delta_Neff_drmd_; }
double z_stop() const { return z_stop_; }
```

Private:

```cpp
double f_idm_drmd_      = 0.;
double G_over_aH_drmd_  = 0.;
double delta_Neff_drmd_ = 0.;
double z_stop_          = 0.;
```

- [ ] **Step 2: Update `species/idm_drmd_idr_drmd_species.cpp` constructor**

```cpp
IDM_DRMD_IDR_DRMD_Species::IDM_DRMD_IDR_DRMD_Species(
    const background& pba,
    double omega0_idm_drmd,
    double omega0_idr_drmd,
    double f_idm_drmd,
    double G_over_aH_drmd,
    double delta_Neff_drmd,
    double z_stop)
    : CompositeSpecies(/* args */),
      pba_(pba),
      f_idm_drmd_(f_idm_drmd),
      G_over_aH_drmd_(G_over_aH_drmd),
      delta_Neff_drmd_(delta_Neff_drmd),
      z_stop_(z_stop) {
  has_idm_drmd_ = (omega0_idm_drmd != 0.);
  has_idr_drmd_ = (omega0_idr_drmd != 0.);
  auto idm = std::make_unique<IDM_DRMDSpecies>(pba, omega0_idm_drmd);
  auto idr = std::make_unique<IDR_DRMDSpecies>(pba, omega0_idr_drmd);
  idm_drmd_ = idm.get();
  idr_drmd_ = idr.get();
  children_.push_back(std::move(idm));
  children_.push_back(std::move(idr));
}
```

- [ ] **Step 3: Update `species/idm_drmd_idr_drmd_species.cpp` `CreateAll`**

Today `CreateAll` reads omegas from budget. Add parses for the 4 fields:

```cpp
// Inside CreateAll, after existence:
double f_idm_drmd = 0.;
double G_over_aH_drmd = 0.;
double delta_Neff_drmd = 0.;
double z_stop = 0.;
ctx.pfc->read_double("f_idm_drmd", f_idm_drmd);
ctx.pfc->read_double("G_over_aH_drmd_ini", G_over_aH_drmd);
ctx.pfc->read_double("delta_Neff_drmd", delta_Neff_drmd);
ctx.pfc->read_double("z_stop", z_stop);

result.push_back({"IDM_DRMD_IDR_DRMD",
                  std::make_unique<IDM_DRMD_IDR_DRMD_Species>(
                      *ctx.pba, omega0_idm_drmd, omega0_idr_drmd,
                      f_idm_drmd, G_over_aH_drmd, delta_Neff_drmd, z_stop)});
```

- [ ] **Step 4: Delete the 4 DRMD writes in `source/input_module.cpp` ReadCoupledOmegaBudget**

The block around the DRMD parsing (search for `pba->f_idm_drmd =`). The budget computation needs `f_idm_drmd` (for CDM subtraction) and `delta_Neff_drmd` (for Omega0_idr_drmd derivation); parse them to locals, use them for budget math, but don't write to pba.

- [ ] **Step 5: Migrate `source/background_module.cpp` readers**

6 sites (lines 509, 512, 1361, 1363, 1364, 1540 per the earlier survey). Hoist composite lookup:

```cpp
const IDM_DRMD_IDR_DRMD_Species* drmd = nullptr;
if (auto* p = all_species_.find("IDM_DRMD_IDR_DRMD"))
  drmd = static_cast<const IDM_DRMD_IDR_DRMD_Species*>(p->get());
```

Replace `pba->z_stop` with `drmd->z_stop()`, etc.

- [ ] **Step 6: Delete from `source/background.h`**

Delete lines 141–144 (`f_idm_drmd`, `G_over_aH_drmd`, `delta_Neff_drmd`, `z_stop`).

- [ ] **Step 7: Build + test + commit**

```
idm_drmd: composite owns f_idm_drmd / G_over_aH_drmd / delta_Neff_drmd / z_stop

IDM_DRMD_IDR_DRMD composite stores the 4 DRMD coupling fields as private
members with getters; CreateAll parses them from pfc.
ReadCoupledOmegaBudget keeps inline parses for f_idm_drmd (CDM subtraction)
and delta_Neff_drmd (Omega0_idr_drmd derivation) but no longer writes to pba.
background_module's 6 reads migrate to composite accessor.

84/84 scenarios pass at TEST_LEVEL=1.
```

---

## Task 6: Cython wrapper accessor for T_idr

`classy.pyx:1481` reads `self.ba.T_idr` to compute `xi_idr`. After Task 4, `pba->T_idr` is gone; add a `BackgroundModule::GetTIdr()` lookup and use it.

**Files:**
- Modify: `source/background_module.h` — add `double GetTIdr() const;` next to `GetOmega0Species`
- Modify: `source/background_module.cpp` — implement (delegate to `all_species_.find("IDM_DR_IDR")->idr().T_idr()`)
- Modify: `classy.pyx` — replace `self.ba.T_idr` with `deref(bam).GetTIdr()`
- (`cclassy.pxd` is auto-generated — do NOT edit; generate_wrapper.py picks up the new method)

- [ ] **Step 1: Add `GetTIdr` declaration to `source/background_module.h`**

```cpp
// Inside class BackgroundModule, in the public section near GetOmega0Species:
double GetTIdr() const;
```

- [ ] **Step 2: Implement in `source/background_module.cpp`**

```cpp
double BackgroundModule::GetTIdr() const {
  if (auto* p = all_species_.find("IDM_DR_IDR")) {
    const auto& comp = static_cast<const IDM_DR_IDR_Species&>(**p);
    return comp.idr().T_idr();
  }
  return 0.;
}
```

- [ ] **Step 3: Update `classy.pyx`**

Find line 1481:

```python
elif name == 'xi_idr':
    value = self.ba.T_idr/self.ba.T_cmb
```

Replace with:

```python
elif name == 'xi_idr':
    bam_tidr = deref(self._thisptr).GetBackgroundModule()
    value = deref(bam_tidr).GetTIdr()/self.ba.T_cmb
```

(Match the two-step pattern used elsewhere in the file for module accessors.)

- [ ] **Step 4: Rebuild Cython side + run scenarios**

```bash
pip install . 2>&1 | tail -3
TEST_LEVEL=1 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -5
```

Expected: `84 passed`. The `cclassy.pxd` auto-regenerates and includes `GetTIdr`.

- [ ] **Step 5: Commit**

```bash
git add -- source/background_module.h source/background_module.cpp classy.pyx
git commit -m "classy: route xi_idr through BackgroundModule::GetTIdr()

T_idr is no longer on pba; classy.pyx reads it via a new GetTIdr()
accessor on BackgroundModule that resolves through IDM_DR_IDR_Species
to its IDR sub-species.  cclassy.pxd auto-picked-up by generate_wrapper.py."
```

---

## Task 7: Sweep + level-2 verification + PR

Run one final sweep for any residual `pba->X` references to the deleted fields (we may have missed an obscure caller), then `TEST_LEVEL=2`.

- [ ] **Step 1: Residual reference sweep**

```bash
grep -rnE "pba(_)?(->|\.)(fluid_equation_of_state|w0_fld|wa_fld|cs2_fld|Omega_EDE|use_ppf|c_gamma_over_c_fld|attractor_ic_scf|phi_ini_scf|phi_prime_ini_scf|scf_parameters|scf_tuning_index|Gamma_dcdm|Omega_ini_dcdm|T_idr|l_max_idr|f_idm_drmd|G_over_aH_drmd|delta_Neff_drmd|z_stop)\b" --include='*.cpp' --include='*.h' --include='*.c' source species include tools main
```

Expected: zero matches (or only matches inside comments). If anything turns up, migrate it using the appropriate species lookup pattern, build, scenarios, commit.

- [ ] **Step 2: clang-format the diff**

```bash
git diff --name-only master...HEAD | grep -E '\.(cpp|c|h|hpp)$' | grep -v '^hyrec/' | xargs clang-format -i
git diff --stat
```

If clang-format changed anything, commit:

```bash
git add -u
git commit -m "chore: clang-format all modified files"
```

- [ ] **Step 3: TEST_LEVEL=2 scenarios**

```bash
TEST_LEVEL=2 python -m pytest -q -m test_scenario python/test_class.py 2>&1 | tail -10
```

Expected: `1260 passed, 1274 deselected in ~4.5 min`.

- [ ] **Step 4: Push branch**

```bash
git push origin worktree-remove-species-physics-params:275-remove-species-physics-params --set-upstream 2>&1 | tail -5
```

(Branch name follows the `<issue>-<topic>` convention established by PR #276 — same issue #275.)

- [ ] **Step 5: Open PR**

```bash
gh pr create --repo AarhusCosmology/CLASSpp --base master \
  --head 275-remove-species-physics-params \
  --title "Remove species-specific physics params from pba (#275, part 2)" \
  --body "$(cat <<'EOF'
## Summary

Closes #275.  Follow-up to PR #276 (which removed `pba->has_<species>` and
`pba->Omega0_<species>`).  Moves the remaining 20 species-specific physics
parameters off `pba` and onto the owning species:

- **Fluid (7):** `fluid_equation_of_state`, `w0_fld`, `wa_fld`, `cs2_fld`,
  `Omega_EDE`, `use_ppf`, `c_gamma_over_c_fld`.
- **ScalarField (5):** `attractor_ic_scf`, `phi_ini_scf`, `phi_prime_ini_scf`,
  `scf_parameters`, `scf_tuning_index`.
- **DCDM (2):** `Gamma_dcdm`, `Omega_ini_dcdm` — owned by the DCDM sub-species
  inside DCDM_DR composite.
- **IDR (2):** `T_idr`, `l_max_idr` — owned by the IDR sub-species inside
  IDM_DR_IDR composite.
- **IDM_DRMD (4):** `f_idm_drmd`, `G_over_aH_drmd`, `delta_Neff_drmd`,
  `z_stop` — owned by the IDM_DRMD_IDR_DRMD composite (coupling params).

Each species gains per-field private members with public const getters; the
constructor takes one positional arg per field; `CreateAll` parses every
field from `pfc` and feeds the constructor.  Cross-module readers
(`background_module`, `thermodynamics_module`, `perturbations_module`,
`nonlinear_module`) migrate to species lookups via `bgm_->all_species_.find(key)`.

`pba` retains only the universal state (H0, T_cmb, Omega0_g/b/k, h, K, sgnK,
a_today, info/inter flags, verbose, closure_species, background_method).

## Test Plan

- [x] `make -j8 class` clean
- [x] `pip install .` (Cython rebuild succeeds; `cclassy.pxd` auto-regenerated)
- [x] `TEST_LEVEL=1 pytest -m test_scenario` → 84/84 pass
- [x] `TEST_LEVEL=2 pytest -m test_scenario` → 1260/1260 pass
- [ ] reviewer: spot-check classy.pyx `xi_idr` derived param against a master build
- [ ] reviewer: confirm shooter convergence on a DCDM scenario with `Omega_ini_dcdm` target

Design spec: `docs/superpowers/specs/2026-05-29-remove-species-physics-params-design.md`
Implementation plan: `docs/superpowers/plans/2026-05-29-remove-species-physics-params.md`

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)" 2>&1 | tail -5
```

Report the PR URL.

---

## Notes for the executor

- **clangd diagnostics are LSP false positives.** You'll see `'common.h' file not found`, `Use of undeclared identifier '_TRUE_'`, etc. constantly. Ignore — only trust `g++` errors from `make`.
- **`cclassy.pxd` is auto-generated** by `generate_wrapper.py` during `pip install`. Do not edit it manually — generate from the C++ header instead.
- **Backstop pattern:** PR #276 deleted the temporary backstop. There is no equivalent here — once a `pba` field is gone, all writers AND readers must be migrated in the same commit.
- **Order matters within a commit:** add the species-side machinery (fields, constructor, accessor) and update the parser block FIRST, then migrate readers, then delete the pba field. The build will stay broken until all readers are migrated; that's expected.
- **Subagent budget:** Each Task in this plan (1–6) is roughly the size of a Phase from PR #276 (one subagent call, ~50-100 tool uses). Task 7 is one final sweep + push.
- **If a step's grep returns hits you didn't expect**, that's the build's way of telling you you missed a consumer. Migrate it and rebuild before committing.
