# Input Module Redesign Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure `source/input_module.cpp` into three explicit phases (read context → construct species → read derived), fixing the dump-order bug, single-sourcing the photon density formula, removing the coupled-cluster double-parses, and migrating the seven IDM_DR/IDR physics params off `thermo`/`perturbs` onto the species.

**Architecture:** The constructor runs `input_read_precisions → input_default_params → ReadContext → ConstructSpecies → ReadDerived → WriteParameterFiles`. Phase-i (`ReadContext`) reads only what building species needs; phase-iii (`ReadDerived`) reads everything else and asks the *built* species for densities/EoS. A new `CoupledClusterInputs` carrier holds raw intermediates parsed once and shared with the coupled factories. The IDM_DR/IDR interaction params move onto `IDM_DRSpecies`/`IDRSpecies`, read by `thermodynamics_module`/`perturbations_module` via typed accessors on the sanctioned `IDM_DR_IDR` downcast.

**Tech Stack:** C++17, CLASS Makefile build (`make class`), Python pytest harness (`python/test_class.py`), tolerance comparator (`test/scenarios/compare_tol.py`).

**Design spec:** `docs/superpowers/specs/2026-06-03-input-module-redesign-design.md`

**Verification model:** This is a behavior-preserving refactor. Verification is **characterization-based**: capture reference outputs before any change (Task 1), then after each task rebuild, re-run the fixtures, and compare with `compare_tol.py` at **~0.1% tolerance, zero-crossing-aware** (never bit-identical). Genuinely-new pure code (the photon formula helper) additionally gets a real unit test.

**Conventions:**
- Build the CLI binary: `make class` (from repo root). Build is incremental.
- Run a scenario: `./class <path/to.ini>` — writes `<root>*.dat` files (root set in the `.ini`).
- Commit after every green task. Branch is `input-module-redesign` (already created).
- Git commit trailer: `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`

---

## Phase 0 — Characterization baseline

### Task 1: Add fixtures and capture the golden baseline

**Files:**
- Create: `test/scenarios/redesign/lcdm.ini`, `ncdm.ini`, `idmdr_ethos.ini`, `idmdr_nadm.ini`, `idmdrmd.ini`, `dcdm.ini`, `scf.ini`, `fld_halofit.ini`, `s8.ini`, `thetas.ini`
- Create: `test/scenarios/redesign/run.sh`

- [ ] **Step 1: Create the fixture `.ini` files**

Most scenarios reuse known-good params from `test/scenarios/*.ini`. Each fixture sets a distinct `root` under `test/scenarios/redesign/out/`. Create each file exactly:

`test/scenarios/redesign/lcdm.ini`:
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
N_ur = 3.044
write parameters = yes
root = test/scenarios/redesign/out/lcdm_
```

`test/scenarios/redesign/ncdm.ini`:
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
N_ur = 2.0308
N_ncdm = 1
m_ncdm = 0.06
root = test/scenarios/redesign/out/ncdm_
```

`test/scenarios/redesign/idmdr_ethos.ini` (ETHOS path: `a_idm_dr`, free-streaming idr — exercises Blocker 5 params):
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
omega_b = 0.022032
omega_cdm = 0.10
N_ur = 2.0
Omega_idm_dr = 0.01
xi_idr = 0.5
stat_f_idr = 0.875
a_idm_dr = 1e3
idr_nature = free_streaming
input_verbose = 1
root = test/scenarios/redesign/out/idmdr_ethos_
```

`test/scenarios/redesign/idmdr_nadm.ini` (NADM path: `Gamma_0_nadm` — exercises the conversion + printf and the `nindex=0, idr_nature=fluid` defaulting):
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
omega_b = 0.022032
omega_cdm = 0.10
N_ur = 2.0
Omega_idm_dr = 0.01
N_idr = 0.5
Gamma_0_nadm = 1e-2
input_verbose = 2
root = test/scenarios/redesign/out/idmdr_nadm_
```

`test/scenarios/redesign/idmdrmd.ini` (DRMD double-parse — Blocker 4):
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
omega_b = 0.022032
omega_cdm = 0.10
N_ur = 2.0
z_stop = 1.0e4
G_over_aH_drmd_ini = 1.0
f_idm_drmd = 0.1
delta_Neff_drmd = 0.5
input_verbose = 1
root = test/scenarios/redesign/out/idmdrmd_
```

`test/scenarios/redesign/dcdm.ini`:
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_dcdmdr = 0.01
Gamma_dcdm = 100.0
root = test/scenarios/redesign/out/dcdm_
```

`test/scenarios/redesign/scf.ini` (scalar-field closure):
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
Omega_Lambda = 0
Omega_fld = 0
Omega_scf = -1
attractor_ic_scf = yes
scf_parameters = 10.0, 0.0, 0.0, 0.0, 0.0, 0.0
root = test/scenarios/redesign/out/scf_
```

`test/scenarios/redesign/fld_halofit.ini` (Halofit + CLP fluid with `wa_fld != 0` — exercises the `pk_eq` gate in Task 6):
```
output = mPk
P_k_max_1/Mpc = 1.0
non linear = halofit
Omega_Lambda = 0
w0_fld = -0.9
wa_fld = 0.1
pk_eq = yes
root = test/scenarios/redesign/out/fld_halofit_
```

`test/scenarios/redesign/s8.ini` (S8 → sigma8 conversion depends on CDM density — Task 5):
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
S8 = 0.82
root = test/scenarios/redesign/out/s8_
```

`test/scenarios/redesign/thetas.ini` (`100*theta_s` shooting):
```
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
100*theta_s = 1.041783
root = test/scenarios/redesign/out/thetas_
```

- [ ] **Step 2: Create the run helper `test/scenarios/redesign/run.sh`**

```bash
#!/usr/bin/env bash
# Run every redesign fixture, writing .dat outputs under out/. Usage: bash run.sh
set -e
cd "$(git rev-parse --show-toplevel)"
rm -rf test/scenarios/redesign/out
mkdir -p test/scenarios/redesign/out
for ini in test/scenarios/redesign/*.ini; do
  echo "=== $ini ==="
  ./class "$ini"
done
```

- [ ] **Step 3: Build the binary on the current (pre-change) tree**

Run: `make class`
Expected: compiles to `./class` with no errors.

- [ ] **Step 4: Capture the golden baseline**

Run:
```bash
bash test/scenarios/redesign/run.sh
cp -r test/scenarios/redesign/out /tmp/redesign_baseline
ls /tmp/redesign_baseline/*.dat | head
```
Expected: every fixture runs without error; `/tmp/redesign_baseline/` holds the `.dat` files. Confirm `idmdr_nadm` printed the ETHOS/NADM equivalence line and `lcdm` produced `test/scenarios/redesign/out/lcdm_parameters.ini`.

- [ ] **Step 5: Record the write-parameters bug baseline (Blocker 2 evidence)**

Run:
```bash
grep -c 'm_ncdm\|N_ur\|Omega_scf\|Omega_idm_dr' test/scenarios/redesign/out/*unused_parameters 2>/dev/null || echo "no unused file"
```
Then inspect: `cat test/scenarios/redesign/out/lcdm_unused_parameters`.
Expected (the bug): species params like `N_ur` appear in `*unused_parameters` even though they were used. Note this — Task 4 fixes it.

- [ ] **Step 6: Commit fixtures**

```bash
git add test/scenarios/redesign/*.ini test/scenarios/redesign/run.sh
git commit -m "Add characterization fixtures for input-module redesign

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```
(The `out/` dir and `/tmp/redesign_baseline` are not committed.)

---

## Phase 1 — Reorder (Blockers 1, 2, 3)

### Task 2: Add the photon `T_cmb ⇄ Omega0_g` formula to `PhotonsSpecies` (with unit test)

**Files:**
- Modify: `species/photons.h` (add two `static` declarations)
- Modify: `species/photons.cpp` (add definitions)
- Create: `species/photons_formula_test.cpp`
- Modify: `Makefile` (add `test-photons` target)

- [ ] **Step 1: Write the failing test `species/photons_formula_test.cpp`**

```cpp
#include "photons.h"

#include <cassert>
#include <cmath>
#include <cstdio>

int main() {
  const double h = 0.67556, T = 2.7255;

  // Round-trip identity in both directions.
  const double og = PhotonsSpecies::Omega0gFromTcmb(T, h);
  assert(std::fabs(PhotonsSpecies::TcmbFromOmega0g(og, h) - T) < 1e-10);
  const double og2 = 6.0e-5;
  assert(std::fabs(PhotonsSpecies::Omega0gFromTcmb(PhotonsSpecies::TcmbFromOmega0g(og2, h), h) -
                   og2) < 1e-15);

  // Sanity: photon Omega0 at default cosmology is ~5e-5.
  assert(og > 4.0e-5 && og < 7.0e-5);

  std::printf("photons formula tests passed (Omega0_g=%.6e)\n", og);
  return 0;
}
```

- [ ] **Step 2: Add the `test-photons` Makefile target**

Insert after the `test-bisection` target (`Makefile:115-116`):
```make
test-photons: photons.opp base_species.opp
	$(CXX) $(OPTFLAG) $(CXXFLAG) -Iinclude -Itools -Isource -Ispecies -I. species/photons_formula_test.cpp $(addprefix build/,$(notdir $^)) -o test-photons $(LIBRARIES)
```

- [ ] **Step 3: Run the test target to verify it fails (method not declared)**

Run: `make test-photons`
Expected: FAIL — `Omega0gFromTcmb`/`TcmbFromOmega0g` are not members of `PhotonsSpecies`.

- [ ] **Step 4: Declare the helpers in `species/photons.h`**

In the `public:` section of `PhotonsSpecies` (e.g. right after the `static std::vector<Named> CreateAll(...)` declaration):
```cpp
  /** Photon density parameter from CMB temperature.
   *  Omega0_g = (4 sigma_B / c) T^4 / (3 c^2 1e10 h^2 / Mpc_over_m^2 / 8 pi G),
   *  with sigma_B = 2 pi^5 k_B^4 / (15 h_P^3 c^2). Single source of truth for the
   *  T_cmb <-> Omega0_g conversion previously inlined in InputModule. */
  static double Omega0gFromTcmb(double T_cmb, double h);
  /** Inverse of Omega0gFromTcmb. */
  static double TcmbFromOmega0g(double Omega0_g, double h);
```

- [ ] **Step 5: Define the helpers in `species/photons.cpp`**

Add at the top of the file (after includes), or near `CreateAll`:
```cpp
double PhotonsSpecies::Omega0gFromTcmb(double T_cmb, double h) {
  const double sigma_B = 2. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 2);
  return (4. * sigma_B / _c_ * pow(T_cmb, 4.)) /
         (3. * _c_ * _c_ * 1.e10 * h * h / _Mpc_over_m_ / _Mpc_over_m_ / 8. / _PI_ / _G_);
}

double PhotonsSpecies::TcmbFromOmega0g(double Omega0_g, double h) {
  const double sigma_B = 2. * pow(_PI_, 5) * pow(_k_B_, 4) / 15. / pow(_h_P_, 3) / pow(_c_, 2);
  return pow(Omega0_g *
                 (3. * _c_ * _c_ * 1.e10 * h * h / _Mpc_over_m_ / _Mpc_over_m_ / 8. / _PI_ / _G_) /
                 (4. * sigma_B / _c_),
             0.25);
}
```
If `species/photons.cpp` does not already include `<cmath>`/the constants header, add `#include <cmath>` and ensure `background.h`/`common.h` (which define `_PI_`, `_c_`, `_G_`, `_k_B_`, `_h_P_`, `_Mpc_over_m_`) are included — they already are via `photons.h`.

- [ ] **Step 6: Run the test to verify it passes**

Run: `make test-photons && ./test-photons`
Expected: `photons formula tests passed (Omega0_g=5.xxxxxe-05)`

- [ ] **Step 7: Commit**

```bash
git add species/photons.h species/photons.cpp species/photons_formula_test.cpp Makefile
git commit -m "Add PhotonsSpecies T_cmb<->Omega0_g static formula + unit test

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 3: Use the helper in the InputModule photon block (behavior-identical)

**Files:**
- Modify: `source/input_module.cpp` (the `Omega_0_g` block, currently `:758-799`)

- [ ] **Step 1: Replace the inline `sigma_B` formula with the helper**

In `input_read_parameters()` at the `Omega_0_g` block, replace the three inline `Omega0_g`/`T_cmb` computations (`input_module.cpp:766-798`) so each branch calls the helper:
- "none of three" branch: `pba->Omega0_g = PhotonsSpecies::Omega0gFromTcmb(pba->T_cmb, pba->h);`
- `T_cmb` given (flag1): `pba->Omega0_g = PhotonsSpecies::Omega0gFromTcmb(param1, pba->h); pba->T_cmb = param1;`
- `Omega_g` given (flag2): `pba->Omega0_g = param2; pba->T_cmb = PhotonsSpecies::TcmbFromOmega0g(pba->Omega0_g, pba->h);`
- `omega_g` given (flag3): `pba->Omega0_g = param3 / pba->h / pba->h; pba->T_cmb = PhotonsSpecies::TcmbFromOmega0g(pba->Omega0_g, pba->h);`

The local `sigma_B` variable (`input_module.cpp:673,678`) is now unused here; leave it for now (Task 4 removes it with the relocation) or delete if no other use remains in the function. `#include "../species/photons.h"` is already pulled via `all_species.h` (`input_module.cpp:23`).

- [ ] **Step 2: Build**

Run: `make class`
Expected: compiles clean.

- [ ] **Step 3: Verify characterization unchanged**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: every line `OK`, exit 0. (Same formula, so identical to baseline.)

- [ ] **Step 4: Commit**

```bash
git add source/input_module.cpp
git commit -m "InputModule: use PhotonsSpecies formula helper for Omega0_g

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 4: Split into ReadContext / ConstructSpecies / ReadDerived; move the dump

**Files:**
- Modify: `source/input_module.h` (declare `ReadContext`, `ReadDerived`, `WriteParameterFiles`; remove `input_init`, `input_read_parameters`)
- Modify: `source/input_module.cpp` (restructure)

This task is a **behavior-preserving relocation**. The S8 and Halofit reads keep their current mechanism (`omega_budget_.cdm`, `pfc` peek) — Tasks 5 and 6 convert them afterward.

- [ ] **Step 1: Update declarations in `source/input_module.h`**

Replace the private method declarations `int input_init();` and `int input_read_parameters();` (`input_module.h:82-83`) with:
```cpp
  void ReadContext();           // phase i: inputs needed to build species
  void ReadDerived();           // phase iii: everything else + species-dependent reads
  void WriteParameterFiles();   // read/unread parameter dump (runs after ReadDerived)
```
Keep `input_read_precisions`, `input_default_params`, `input_default_precision`, `ConstructSpecies`, `ReadCoupledOmegaBudget`.

- [ ] **Step 2: Restructure the constructor (`input_module.cpp:170-179`)**

```cpp
InputModule::InputModule(FileContent& fc) : file_content_(fc) {
  file_content_.mark_all_unread();
  try {
    input_read_precisions();
    input_default_params();
    ReadContext();
    ConstructSpecies();
    ReadDerived();
    WriteParameterFiles();
  }
  catch (const std::runtime_error& e) {
    throw std::invalid_argument(e.what());
  }
}
```
Note `input_default_params()` was previously the first statement of `input_read_parameters()` (`input_module.cpp:682`); it is now called explicitly here. Remove that call from the relocated body.

- [ ] **Step 3: Define `ReadContext()` from the phase-i half of `input_read_parameters`**

Rename `int InputModule::input_read_parameters()` to `void InputModule::ReadContext()` and **keep only** the phase-i reads in it, in their current order:
- local-variable declarations needed by these reads (`errmsg`, `pfc`, `pba`, `ppt`, `pth`, flags/params, `input_verbose`)
- `input_verbose`, `threads` (`:688-708`)
- `gauge` (`:715-727`)
- `a_today` (`:732`)
- `h`/`H0`/`100*theta_s` (`:735-756`)
- `T_cmb`/`Omega_g`/`omega_g` → `Omega0_g` via helper (the Task 3 block, `:758-799`)
- `Omega_b`/`omega_b` (`:801-810`)
- `ceff2_ur`/`cvis2_ur`/`G_eff_ur` (`:815-832`)
- `ReadCoupledOmegaBudget()` call + the idm/idr thermo sub-block (`:839-962`) — *kept here verbatim for now*; it migrates in Phase 3
- `Omega_k` + curvature (`:971-978`)
- `Omega_Lambda`/`Omega_fld`/`Omega_scf` closure selection + the two `class_test`s + the `fluid_present_pfc` snapshot (`:986-1035`)

Replace `return _SUCCESS_;` style with plain returns (method is now `void`; the `class_call`/`class_test` macros still `throw` on failure — confirm they don't `return` an int; if a macro requires a return value, keep the method returning `int` and `return _SUCCESS_`, and have the constructor call it via `class_call`. **Check the macro:** if `class_call`/`class_test` expand to `return _FAILURE_`, keep these methods `int`-returning and call them with `class_call(ReadContext(), error_message_, error_message_)` inside the `try`. Adjust the `.h` signatures to `int` accordingly.)

- [ ] **Step 4: Define `ReadDerived()` from the phase-iii remainder**

Create `void InputModule::ReadDerived()` (or `int`, per the macro decision in Step 3) containing everything from the thermo block onward that was in `input_read_parameters` after closure selection:
- its own local-variable declarations (re-declare `errmsg`, `pfc`, the struct pointers, flags/params, `string1`, `sigma_B` if still referenced, `z_max`, etc. — copy the decl block from the old function head, dropping names not used here)
- `YHe`, recombination, reionization blocks, energy injection (`:1046-1253`)
- the `output`/perturbation/primordial/transfer/spectra/lensing blocks (`:1255` through the tensor/trigger/halofit/l_max blocks ending `:2761`)
- the S8 branch stays as-is for now (`:1686-1692`, still reads `omega_budget_.cdm`)
- the Halofit `pk_eq` peek stays as-is for now (`:2690-2726`)

- [ ] **Step 5: Define `WriteParameterFiles()` from the dump block**

Move the `write parameters` + `write warnings` block (`input_module.cpp:545-599`, currently inside `input_init`) into a new `void InputModule::WriteParameterFiles()`. Delete the now-empty `input_init()` (its `input_read_precisions` call and `input_verbose` "Reading input parameters" print move to the constructor / `ReadContext` respectively — put the "Reading input parameters" print at the top of `ReadContext`).

- [ ] **Step 6: Build**

Run: `make class`
Expected: compiles clean. Fix any local-variable omissions surfaced by the compiler (a phase-iii read referencing a local only declared in the phase-i head, or vice versa).

- [ ] **Step 7: Verify characterization unchanged**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`, exit 0.

- [ ] **Step 8: Verify the write-parameters bug is fixed (Blocker 2 proof)**

Run:
```bash
grep -c 'N_ur' test/scenarios/redesign/out/lcdm_unused_parameters || echo "N_ur not in unused (good)"
grep 'N_ur' test/scenarios/redesign/out/lcdm_parameters.ini
```
Expected: `N_ur` no longer appears in `*unused_parameters` (species reads now precede the dump) and **does** appear in `lcdm_parameters.ini`. Likewise spot-check `Omega_idm_dr` for `idmdr_ethos_` and `Omega_scf` for `scf_`.

- [ ] **Step 9: Commit**

```bash
git add source/input_module.h source/input_module.cpp
git commit -m "InputModule: split into ReadContext/ConstructSpecies/ReadDerived; fix dump order

Constructor now runs the three phases explicitly and emits the read/unread
parameter dump after construction, so species inputs are reported correctly
(fixes the write-parameters mis-report). Behavior otherwise identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 5: S8 → sigma8 asks the CDM species

**Files:**
- Modify: `source/input_module.cpp` (S8 branch in `ReadDerived`, was `:1686-1692`)

- [ ] **Step 1: Replace the budget read with a species query**

In the `S8` branch, replace:
```cpp
const double Omega0_cdm_for_S8 = omega_budget_.cdm.value_or(0.);
```
with a query of the built CDM species (now available — `ReadDerived` runs after `ConstructSpecies`):
```cpp
const double Omega0_cdm_for_S8 =
    all_species_.count("CDM") ? all_species_.at("CDM")->GetOmega0() : 0.0;
```
(Confirm `SpeciesCollection::at` returns a pointer-like usable with `->GetOmega0()`; from `cdm.cpp` the collection stores `std::unique_ptr<BaseSpecies>`, and `thermodynamics_module.cpp:204` dereferences `*all_species_.at(...)`. Use `all_species_.at("CDM")->GetOmega0()` if `at` returns the unique_ptr, else `(*all_species_.at("CDM")).GetOmega0()`.)

- [ ] **Step 2: Build**

Run: `make class`
Expected: clean.

- [ ] **Step 3: Verify**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK` — especially `s8_*.dat` (the resolved CDM Ω is identical to `omega_budget_.cdm`).

- [ ] **Step 4: Commit**

```bash
git add source/input_module.cpp
git commit -m "InputModule: S8->sigma8 reads CDM density from the built species

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 6: Halofit `pk_eq` gate asks the Fluid species; delete the peek

**Files:**
- Modify: `source/input_module.cpp` (Halofit block in `ReadDerived`, was `:2690-2726`; and the `fluid_present_pfc` snapshot from Task 4 / was `:1034`)
- Possibly modify: `species/fluid.h` (add accessors if absent)

- [ ] **Step 1: Check what the Fluid species exposes**

Run: `grep -n 'wa_fld\|w0_fld\|fluid_equation_of_state\|use_ppf\|EoS\|class FluidSpecies\|public:' species/fluid.h`
Determine whether `FluidSpecies` exposes `wa_fld` and an EoS-type query. If not, add minimal const accessors (e.g. `double wa_fld() const { return wa_fld_; }` and `bool eos_is_clp() const { return eos_ == fluid_eos_clp; }`) next to its other getters, returning the values it already parsed in `FluidSpecies::CreateAll`.

- [ ] **Step 2: Replace the pfc peek with a species query**

Replace the Halofit block (the `fluid_equation_of_state` probe + `wa_fld_peek` read + the gate, was `input_module.cpp:2700-2726`) with:
```cpp
if ((pnl->method == nl_halofit) && all_species_.count("Fluid")) {
  const auto& fluid = static_cast<const FluidSpecies&>(*all_species_.at("Fluid"));
  if (fluid.eos_is_clp() && fluid.wa_fld() != 0.) {
    class_call(parser_read_string(pfc, "pk_eq", &string1, &flag1, errmsg), errmsg, errmsg);
    if ((flag1 == _TRUE_) &&
        ((strstr(string1, "y") != nullptr) || (strstr(string1, "Y") != nullptr))) {
      pnl->has_pk_eq = _TRUE_;
    }
  }
}
```
`#include "../species/fluid.h"` if not already pulled via `all_species.h` (it is). Then **delete** the `fluid_present_pfc` snapshot line and its uses (it was only feeding this gate).

- [ ] **Step 3: Build**

Run: `make class`
Expected: clean.

- [ ] **Step 4: Verify**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK` — especially `fld_halofit_*.dat` (the `pk_eq` decision is unchanged).

- [ ] **Step 5: Commit**

```bash
git add source/input_module.cpp species/fluid.h
git commit -m "InputModule: Halofit pk_eq gate queries the built Fluid species; drop pfc peek

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Phase 2 — Blocker 4: `CoupledClusterInputs` carrier

### Task 7: Add the carrier struct, member, and context field; populate it

**Files:**
- Modify: `species/species_build_context.h` (add `CoupledClusterInputs`, add context field)
- Modify: `source/input_module.h` (add `coupled_inputs_` member)
- Modify: `source/input_module.cpp` (populate it in `ReadCoupledOmegaBudget`, rename to `ReadCoupledCluster`)

- [ ] **Step 1: Add `CoupledClusterInputs` to `species/species_build_context.h`**

After the `SpeciesOmegaBudget` struct, add:
```cpp
/**
 * Raw intermediates parsed once in phase i that the coupled factories need for
 * physics construction (not just budget math). Single source of truth for the
 * coupled cluster's reusable parsed values, so factories don't re-parse pfc
 * identically to the budget resolver. See spec Blocker 4.
 */
struct CoupledClusterInputs {
  std::optional<double> T_idr;       // resolved from N_idr | N_dg | xi_idr (+ stat_f_idr)
  double stat_f_idr = 7. / 8.;
  double f_idm_drmd = 0.;
  double G_over_aH_drmd_ini = 0.;
  double delta_Neff_drmd = 0.;
  double z_stop = 0.;
};
```
Add to `SpeciesBuildContext`:
```cpp
  const CoupledClusterInputs* coupled_inputs = nullptr;
```

- [ ] **Step 2: Add the member to `source/input_module.h`**

Next to `SpeciesOmegaBudget omega_budget_;` (`input_module.h:52`):
```cpp
  CoupledClusterInputs coupled_inputs_;
```

- [ ] **Step 3: Populate it in the budget resolver**

In `ReadCoupledOmegaBudget()` (rename the method and its `.h` declaration to `ReadCoupledCluster()` for accuracy), after `T_idr_local` is computed (`input_module.cpp:318-347`), store:
```cpp
  coupled_inputs_.stat_f_idr = stat_f_idr;
  if (flag1 == _TRUE_ || flag2 == _TRUE_ || flag3 == _TRUE_)
    coupled_inputs_.T_idr = T_idr_local;
```
And in the DRMD block, after the four `parser_read_double`s (`:444-453`), store the resolved locals:
```cpp
  coupled_inputs_.z_stop = z_stop_drmd;
  coupled_inputs_.f_idm_drmd = f_idm_drmd_local;
  coupled_inputs_.delta_Neff_drmd = delta_Neff_drmd_local;
  coupled_inputs_.G_over_aH_drmd_ini = (flag2 == _TRUE_) ? param2 : 0.;
```

- [ ] **Step 4: Thread `coupled_inputs_` into the build context**

In `ConstructSpecies()`, immediately after the `SpeciesBuildContext ctx{...}` aggregate is built (`input_module.cpp:195-203`), assign the new field:
```cpp
  ctx.coupled_inputs = &coupled_inputs_;
```
(`coupled_inputs` is a defaulted trailing member of the struct, so the existing aggregate initializer keeps compiling; set it by name afterward.) Do the same wherever `DoShooting` constructs a `SpeciesBuildContext` for its per-species guess loop — search `DoShooting` for `SpeciesBuildContext` and set `ctx.coupled_inputs = &<im>->coupled_inputs_;` pointing at the shooting `InputModule`'s member.

- [ ] **Step 5: Build (no behavior change yet — factories still self-parse)**

Run: `make class`
Expected: clean.

- [ ] **Step 6: Verify unchanged**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`.

- [ ] **Step 7: Commit**

```bash
git add species/species_build_context.h source/input_module.h source/input_module.cpp
git commit -m "Add CoupledClusterInputs carrier; populate from the budget resolver

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 8: `IDM_DR_IDR_Species::CreateAll` reads `T_idr` from the carrier

**Files:**
- Modify: `species/idm_dr_idr_species.cpp` (`CreateAll`, `:437-456`)

- [ ] **Step 1: Delete the duplicate `T_idr` parse, read from the carrier**

In `IDM_DR_IDR_Species::CreateAll`, replace the local `stat_f_idr`/`N_idr`/`N_dg`/`xi_idr` parse + `T_idr_local` derivation (`idm_dr_idr_species.cpp:437-456`) with:
```cpp
    const double T_idr_local = (ctx.coupled_inputs && ctx.coupled_inputs->T_idr)
                                   ? *ctx.coupled_inputs->T_idr
                                   : 0.;
```
Keep the rest (`l_max_idr_local`, the `make_unique`). If `ctx.coupled_inputs` is null (standalone/test callers), `T_idr_local` falls back to 0., matching the prior behavior when no IDR-temperature input was given.

- [ ] **Step 2: Build**

Run: `make class`
Expected: clean.

- [ ] **Step 3: Verify**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK` — especially `idmdr_ethos_*.dat` and `idmdr_nadm_*.dat` (same `T_idr`, one parse).

- [ ] **Step 4: Commit**

```bash
git add species/idm_dr_idr_species.cpp
git commit -m "IDM_DR_IDR: read T_idr from CoupledClusterInputs (kill double-parse)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 9: `IDM_DRMD_IDR_DRMD_Species::CreateAll` reads DRMD keys from the carrier

**Files:**
- Modify: `species/idm_drmd_idr_drmd_species.cpp` (`CreateAll`, `:392-395`)

- [ ] **Step 1: Delete the duplicate DRMD parse, read from the carrier**

Replace the four `ctx.pfc->read_double(...)` calls (`idm_drmd_idr_drmd_species.cpp:392-395`) and their local declarations with:
```cpp
    const auto* ci = ctx.coupled_inputs;
    const double f_idm_drmd      = ci ? ci->f_idm_drmd : 0.;
    const double G_over_aH_drmd  = ci ? ci->G_over_aH_drmd_ini : 0.;
    const double delta_Neff_drmd = ci ? ci->delta_Neff_drmd : 0.;
    const double z_stop          = ci ? ci->z_stop : 0.;
```
Leave the `make_unique<IDM_DRMD_IDR_DRMD_Species>(...)` call passing the same four values.

- [ ] **Step 2: Build**

Run: `make class`
Expected: clean.

- [ ] **Step 3: Verify**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK` — especially `idmdrmd_*.dat`.

- [ ] **Step 4: Commit**

```bash
git add species/idm_drmd_idr_drmd_species.cpp
git commit -m "IDM_DRMD_IDR_DRMD: read DRMD params from CoupledClusterInputs (kill double-parse)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Phase 3 — Blocker 5: IDM_DR/IDR param ownership migration

The seven params move onto the species incrementally. Tasks 10–11 add the storage and parse it on the species **while the monolith still parses to `pth`/`ppt`** (transient double-parse, behavior-identical). Tasks 12–14 repoint consumers one group at a time (thermo, then alpha/beta, then `idr_nature` alone). Task 15 deletes the monolith block and the global fields. Every task builds and is verified.

### Task 10: Add storage, ctor params, and accessors to `IDM_DRSpecies` / `IDRSpecies`

**Files:**
- Modify: `species/idm_dr.h` (members + accessors + ctor params on `IDM_DRSpecies`)
- Modify: `species/idr.h` (members + accessors + ctor params on `IDRSpecies`)

- [ ] **Step 1: Extend `IDM_DRSpecies` (`species/idm_dr.h`)**

Change the constructor and add members/accessors:
```cpp
  IDM_DRSpecies(const background& pba, double omega0_idm_dr,
                double a_idm_dr = 0., double nindex_idm_dr = 4., double m_idm = 1.e11)
      : BaseSpecies("IDM_DR", EnergyType::Matter), pba_(pba), Omega0_idm_dr_(omega0_idm_dr),
        a_idm_dr_(a_idm_dr), nindex_idm_dr_(nindex_idm_dr), m_idm_(m_idm) {}

  double a_idm_dr() const { return a_idm_dr_; }
  double nindex_idm_dr() const { return nindex_idm_dr_; }
  double m_idm() const { return m_idm_; }
```
Add to the `private:` section (defaults mirror the old `thermo` struct: `a_idm_dr=0`,
`nindex_idm_dr=4.`, `m_idm=1.e11`):
```cpp
  double a_idm_dr_ = 0.;
  double nindex_idm_dr_ = 4.;
  double m_idm_ = 1.e11;
```

- [ ] **Step 2: Extend `IDRSpecies` (`species/idr.h`)**

Change the constructor signature to also accept `b_idr`, `idr_nature`, and the `alpha`/`beta` vectors, and add accessors:
```cpp
  IDRSpecies(const background& pba,
             double omega0_idr,
             bool has_sibling_idm_dr = false,
             double T_idr            = 0.,
             int l_max_idr           = 0,
             double b_idr            = 0.,
             int idr_nature          = 0,
             std::vector<double> alpha_idm_dr = {},
             std::vector<double> beta_idr     = {})
      : BaseSpecies("IDR", EnergyType::Radiation), pba_(pba), Omega0_idr_(omega0_idr),
        has_sibling_idm_dr_(has_sibling_idm_dr), T_idr_(T_idr), l_max_idr_(l_max_idr),
        b_idr_(b_idr), idr_nature_(idr_nature),
        alpha_idm_dr_(std::move(alpha_idm_dr)), beta_idr_(std::move(beta_idr)) {}

  double b_idr() const { return b_idr_; }
  int idr_nature() const { return idr_nature_; }
  const std::vector<double>& alpha_idm_dr() const { return alpha_idm_dr_; }
  const std::vector<double>& beta_idr() const { return beta_idr_; }
```
Add to `private:`:
```cpp
  double b_idr_ = 0.;
  int idr_nature_ = 0;
  std::vector<double> alpha_idm_dr_;
  std::vector<double> beta_idr_;
```
Add `#include <vector>` to `idr.h` if not already present.

- [ ] **Step 3: Build (members unused — pure scaffolding)**

Run: `make class`
Expected: clean (default ctor args keep all existing call sites valid).

- [ ] **Step 4: Verify unchanged**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`.

- [ ] **Step 5: Commit**

```bash
git add species/idm_dr.h species/idr.h
git commit -m "IDM_DR/IDR species: add storage+accessors for migrating interaction params

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 11: Parse the interaction params in `IDM_DR_IDR_Species::CreateAll` and store on the children

**Files:**
- Modify: `species/idm_dr_idr_species.h` (composite ctor gains the new params)
- Modify: `species/idm_dr_idr_species.cpp` (`CreateAll` parses them; composite ctor forwards to children)
- Modify: `source/input_module.cpp` (remove the **printf** lines from the monolith block to avoid double-printing; keep the `pth`/`ppt` writes)

- [ ] **Step 1: Extend the composite constructor**

In `species/idm_dr_idr_species.h:31-33`, change to:
```cpp
  IDM_DR_IDR_Species(const background& pba, double omega0_idm_dr, double omega0_idr,
                     double T_idr, int l_max_idr,
                     double a_idm_dr, double nindex_idm_dr, double m_idm, double b_idr,
                     int idr_nature, std::vector<double> alpha_idm_dr,
                     std::vector<double> beta_idr);
```
In `species/idm_dr_idr_species.cpp:203-214`, forward the new args into the children:
```cpp
  auto idm = std::make_unique<IDM_DRSpecies>(pba, omega0_idm_dr, a_idm_dr, nindex_idm_dr, m_idm);
  auto idr = std::make_unique<IDRSpecies>(pba, omega0_idr, has_idm_dr_, T_idr, l_max_idr,
                                          b_idr, idr_nature, std::move(alpha_idm_dr),
                                          std::move(beta_idr));
```

- [ ] **Step 2: Parse the params in `CreateAll`**

In `IDM_DR_IDR_Species::CreateAll` (after `T_idr_local` from Task 8, before the `make_unique`), replicate the monolith parse logic (`input_module.cpp:848-961`) against `ctx.pfc`, producing locals: `a_idm_dr`, `nindex_idm_dr`, `m_idm`, `b_idr`, `idr_nature`, `alpha_idm_dr` (vector sized `ctx.ppr->l_max_idr - 1`), `beta_idr` (same). Use the resolved `omega0_idr` for the `Gamma_0_nadm → a_idm_dr` conversion and read `input_verbose` from `ctx.pfc` for the equivalence printfs. The exact branch structure (verbatim from the monolith, reading `ctx.pfc`):

```cpp
    int input_verbose = 0;
    ctx.pfc->read_int("input_verbose", input_verbose);
    const double h2 = ctx.pba->h * ctx.pba->h;

    double a_idm_dr = 0.; double nindex_idm_dr = 4.; double m_idm = 1.e11; double b_idr = 0.;
    int idr_nature = idr_free_streaming;  // == 0, the old perturbs struct default; overridden below
    double a_param = 0., a_dark = 0., gamma0 = 0.;
    bool f_a   = ctx.pfc->read_double("a_idm_dr", a_param);
    bool f_ad  = ctx.pfc->read_double("a_dark", a_dark);
    bool f_g   = ctx.pfc->read_double("Gamma_0_nadm", gamma0);
    // (replicate the at-least-two-of-three test as a throw)
    if (f_a)       a_idm_dr = a_param;
    else if (f_ad) a_idm_dr = a_dark;
    else if (f_g)  a_idm_dr = gamma0 * (3. / 4.) / (h2 * omega0_idr);
    if (f_g) {
      nindex_idm_dr = 0.; idr_nature = idr_fluid;
      if (input_verbose > 1) printf("NADM requested. Defaulting on nindex_idm_dr = 0 and idr_nature = fluid \n");
      if (input_verbose > 1) printf("You passed Gamma_0_nadm = %e, ... a_idm_dr = a_dark = %e ...\n", gamma0, a_idm_dr);
    } else {
      double nidx = 0.;
      if (ctx.pfc->read_double("nindex_dark", nidx) || ctx.pfc->read_double("nindex_idm_dr", nidx))
        nindex_idm_dr = nidx;  // double, default stays 4. (matches old thermo struct default)
      std::string nature;
      if (ctx.pfc->read_string("idr_nature", nature)) {
        if (nature.find("free_stream") != std::string::npos || nature.find("Free") != std::string::npos)
          idr_nature = idr_free_streaming;
        else if (nature.find("fluid") != std::string::npos || nature.find("Fluid") != std::string::npos)
          idr_nature = idr_fluid;
      }
    }
    { double v; if (ctx.pfc->read_double("m_idm", v) || ctx.pfc->read_double("m_dm", v)) m_idm = v; }
    { double v; if (ctx.pfc->read_double("b_dark", v) || ctx.pfc->read_double("b_idr", v)) b_idr = v; }

    std::vector<double> alpha_idm_dr, beta_idr;
    ctx.pfc->read_list_of_doubles("alpha_idm_dr", alpha_idm_dr) ||
        ctx.pfc->read_list_of_doubles("alpha_dark", alpha_idm_dr);
    if (alpha_idm_dr.empty()) alpha_idm_dr.assign(ctx.ppr->l_max_idr - 1, 1.5);
    else if ((int) alpha_idm_dr.size() != ctx.ppr->l_max_idr - 1)
      alpha_idm_dr.resize(ctx.ppr->l_max_idr - 1, alpha_idm_dr.back());
    ctx.pfc->read_list_of_doubles("beta_idr", beta_idr) ||
        ctx.pfc->read_list_of_doubles("beta_dark", beta_idr);
    if (beta_idr.empty()) beta_idr.assign(ctx.ppr->l_max_idr - 1, 1.5);
    else if ((int) beta_idr.size() != ctx.ppr->l_max_idr - 1)
      beta_idr.resize(ctx.ppr->l_max_idr - 1, beta_idr.back());
```
Pass these into the extended `make_unique<IDM_DR_IDR_Species>(...)`. **Cross-check** the branch logic line-by-line against `input_module.cpp:848-961` (alias precedence, the `nindex_idm_dr` cast, the `1.5` defaults, the resize-to-`l_max_idr-1` rule) so the species parse is byte-equivalent to the monolith's. Reuse the project's `readDoubleList`/`read_list_of_doubles` helper actually used at `:921` (match its name/signature). Include `<vector>`, `<string>`, and the enum header providing `idr_fluid`/`idr_free_streaming`.

- [ ] **Step 3: Remove the duplicate printfs from the monolith block**

In `input_module.cpp` (the idm/idr block still in `ReadContext`), delete only the `printf(...)` equivalence lines (`:858-862`, `:866-872`, `:876-881`, `:890-892`). Keep the `pth->...`/`ppt->...` assignments — consumers still read them until Task 15. This prevents double-printing during the transition.

- [ ] **Step 4: Build**

Run: `make class`
Expected: clean. The composite is now constructed with the new args; the species members hold the parsed values (still unused by consumers).

- [ ] **Step 5: Verify unchanged**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`, including identical printed equivalence lines for `idmdr_nadm` (now emitted from `CreateAll` instead of the monolith). Diff the captured stdout for `idmdr_nadm` against baseline if in doubt.

- [ ] **Step 6: Commit**

```bash
git add species/idm_dr_idr_species.h species/idm_dr_idr_species.cpp source/input_module.cpp
git commit -m "IDM_DR_IDR: parse interaction params in CreateAll onto the species

Transitional: monolith still writes pth/ppt for consumers; printfs moved to
CreateAll to avoid duplication. Behavior identical.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 12: Repoint the thermodynamics module to the species (`a_idm_dr`/`nindex`/`m_idm`/`b_idr`)

**Files:**
- Modify: `source/thermodynamics_module.cpp` (every `pth->a_idm_dr`/`nindex_idm_dr`/`b_idr`/`m_idm` read)

- [ ] **Step 1: Replace each global read with the typed accessor**

The module already binds `comp` (the `IDM_DR_IDR_Species&`) at `:203`, `:506`, `:830`. At each read site, source from the composite's children instead of `pth`:
- `pth->a_idm_dr` → `comp.idm_dr().a_idm_dr()`
- `pth->nindex_idm_dr` → `comp.idm_dr().nindex_idm_dr()`
- `pth->m_idm` → `comp.idm_dr().m_idm()`
- `pth->b_idr` → `comp.idr().b_idr()`

Sites (from the audit): `:209`, `:211-212`, `:218`, `:221`, `:242`, and the `Gamma_heat_idm_dr` block at `:878-879`, `:928-930`, `:954-956`, `:994-996`, `:1023-1025`, `:1041`, plus the `:1155` `nindex_idm_dr >= 2` test. Where a site is not already inside a scope that binds `comp`, bind it locally first:
```cpp
const auto& comp = static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR"));
```
(guarded by the existing `if (all_species_.count("IDM_DR_IDR") > 0)`).

- [ ] **Step 2: Build**

Run: `make class`
Expected: clean. (`pth->a_idm_dr` etc. still exist — just no longer read by the thermo module.)

- [ ] **Step 3: Verify**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK` — focus on `idmdr_ethos_*.dat` and `idmdr_nadm_*.dat`.

- [ ] **Step 4: Commit**

```bash
git add source/thermodynamics_module.cpp
git commit -m "thermo module: read a_idm_dr/nindex/m_idm/b_idr from the IDM_DR_IDR species

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 13: Repoint `alpha_idm_dr` / `beta_idr` to the species

**Files:**
- Modify: `species/idm_dr_idr_species.cpp` (`AddCouplingDerivs`, `:388-390`)
- Modify: any other `ppt->alpha_idm_dr`/`ppt->beta_idr` reader (audit found the composite's own `AddCouplingDerivs`)

- [ ] **Step 1: Read alpha/beta from the IDR child instead of `ppt`**

In `AddCouplingDerivs` (`idm_dr_idr_species.cpp:388-390`), replace `ppt->alpha_idm_dr[l - 2]` and `ppt->beta_idr[l - 2]` with the IDR child's vectors:
```cpp
        dy[idr_lay.idx_delta + l] -= (idr_->alpha_idm_dr()[l - 2] * dmu_idm_dr +
                                      idr_->beta_idr()[l - 2] * dmu_idr) *
                                     y[idr_lay.idx_delta + l];
```
Also update the TCA-shear use of `ppt->alpha_idm_dr[0]` (`idm_dr_idr_species.cpp:398`) to `idr_->alpha_idm_dr()[0]`. Confirm via `grep -n 'alpha_idm_dr\|beta_idr' source/perturbations_module.cpp species/*.cpp` that no other site reads these globals; repoint any that do.

- [ ] **Step 2: Build**

Run: `make class`
Expected: clean.

- [ ] **Step 3: Verify**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK` — `idmdr_ethos_*.dat` (free-streaming exercises the `l>=2` collision terms).

- [ ] **Step 4: Commit**

```bash
git add species/idm_dr_idr_species.cpp
git commit -m "IDM_DR_IDR: read alpha_idm_dr/beta_idr from the IDR child

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 14: Repoint `idr_nature` to the species (isolated, highest-risk)

**Files:**
- Modify: `source/perturbations_module.cpp` (`:2301`, `:3658`, `:4638`, `:5265`)
- Modify: `species/idm_dr_idr_species.cpp` (`:52`, `:166`, `:197`, `:382`, `:484`)

- [ ] **Step 1: Re-source the `scalar_ctx` feed and the module's other reads**

At `perturbations_module.cpp:5265`, change `ppw->scalar_ctx.idr_nature = ppt->idr_nature;` to read from the composite:
```cpp
ppw->scalar_ctx.idr_nature =
    static_cast<const IDM_DR_IDR_Species&>(*all_species_.at("IDM_DR_IDR")).idr().idr_nature();
```
(guarded by the existing `if (all_species_.count("IDM_DR_IDR"))` — if `:5265` is unconditional, wrap it so the feed only applies when the species exists, preserving the prior default for the no-IDR case). At `:2301`, `:3658`, `:4638`, replace `ppt->idr_nature` similarly, binding `comp` locally under the existing `IDM_DR_IDR` count guards (`:3659` already binds it). The deep consumers reading `ppw->scalar_ctx.idr_nature` / `ctx.idr_nature` (`interacting_species.cpp:135/178/264/326`, `idm_dr_idr_species.cpp:382`) need **no change** — they read the re-sourced copy.

- [ ] **Step 2: Repoint the composite's own direct `ppt->idr_nature` reads**

In `idm_dr_idr_species.cpp`, replace `ppt->idr_nature` / `mod.GetPerturbs()->idr_nature` / `ppt_->idr_nature` at `:52`, `:166`, `:197`, `:484` with `idr_->idr_nature()` (the composite owns the IDR child via `idr_`). For `:197` (inside `WriteOutputColumns`, which has `mod`), use `idr().idr_nature()`.

- [ ] **Step 3: Build**

Run: `make class`
Expected: clean.

- [ ] **Step 4: Verify (scrutinize this one)**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`. Pay special attention to `idmdr_ethos_*` (free-streaming) **and** `idmdr_nadm_*` (NADM forces `idr_nature=fluid`) — these exercise both branches of `idr_nature`. If any FAIL, the `scalar_ctx` feed guard or a missed read site is wrong; do not proceed.

- [ ] **Step 5: Commit**

```bash
git add source/perturbations_module.cpp species/idm_dr_idr_species.cpp
git commit -m "perturbations: source idr_nature from the IDR species (scalar_ctx feed + direct reads)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

### Task 15: Delete the monolith parse block and the global fields

**Files:**
- Modify: `source/input_module.cpp` (delete the idm/idr thermo parse block in `ReadContext`, was `:841-962` minus the printfs already removed)
- Modify: `source/thermodynamics.h` (delete `a_idm_dr`, `nindex_idm_dr`, `m_idm`, `b_idr`; `:169-174`)
- Modify: `source/perturbations.h` (delete `idr_nature`, `alpha_idm_dr`, `beta_idr`, `alpha_idm_dr_storage`, `beta_idr_storage`; `:200-205`)

- [ ] **Step 1: Delete the monolith idm/idr parse block**

Remove the entire `if (omega_budget_.idm_dr.value_or(0.) > 0.) { ... }` block that parsed `a_idm_dr`/`nindex`/`idr_nature`/`m_idm`/`b_idr`/`alpha`/`beta` into `pth`/`ppt` (`input_module.cpp:841-962`, printfs already gone). Keep the preceding `class_test` that requires non-zero IDR density when IDM_DR is requested (`:841-846`) **only if** it doesn't depend on the deleted writes — re-home that validation into `IDM_DR_IDR_Species::CreateAll` (throw if `omega0_idm_dr > 0 && omega0_idr == 0`) and delete it from the monolith.

- [ ] **Step 2: Delete the now-unread global fields**

In `source/thermodynamics.h` delete the declarations of `a_idm_dr` (`:169`), `b_idr` (`:171`), `nindex_idm_dr` (`:172`), `m_idm` (`:174`). In `source/perturbations.h` delete `alpha_idm_dr` (`:200`), `beta_idr` (`:201`), `alpha_idm_dr_storage` (`:202`), `beta_idr_storage` (`:203`), `idr_nature` (`:205`).

- [ ] **Step 3: Build — fix any stragglers the compiler finds**

Run: `make class`
Expected: it may fail if a forgotten reader of these fields remains. For each error, repoint that reader to the species accessor (it has, or can bind, the `IDM_DR_IDR` composite). Re-run until clean. Run `grep -rn 'a_idm_dr\|nindex_idm_dr\|->m_idm\|->b_idr\|idr_nature\|alpha_idm_dr\|beta_idr' source/ species/` and confirm the only remaining hits are the species members/accessors and `scalar_ctx.idr_nature`.

- [ ] **Step 4: Full characterization verification**

Run:
```bash
make class
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`, exit 0.

- [ ] **Step 5: Run the existing Python test suite (gauge + budget regression)**

Run:
```bash
make classy
cd python && python -m pytest test_class.py -q
```
Expected: passes (cross-gauge comparisons and scenario runs). Investigate any failure before committing — compare against a run of the suite on `master` if a pre-existing failure is suspected.

- [ ] **Step 6: Commit**

```bash
git add source/input_module.cpp source/thermodynamics.h source/perturbations.h
git commit -m "Delete IDM/IDR globals + monolith parse block; ownership now on the species

Completes Blocker 5: a_idm_dr/nindex/m_idm/b_idr off thermo and
idr_nature/alpha_idm_dr/beta_idr off perturbs. The IDM_DR_IDR species owns and
parses them; thermo/perturbation modules read via typed accessors.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Final verification

### Task 16: Whole-suite green + spec cross-check

- [ ] **Step 1: Clean rebuild**

Run: `make clean && make class && make classy`
Expected: clean build of both targets.

- [ ] **Step 2: Full characterization pass**

Run:
```bash
bash test/scenarios/redesign/run.sh
python test/scenarios/compare_tol.py /tmp/redesign_baseline test/scenarios/redesign/out
```
Expected: all `OK`.

- [ ] **Step 3: Existing scenario + gauge suite**

Run: `cd python && python -m pytest test_class.py -q`
Expected: green (or only pre-existing, documented failures — cf. the three historical cross-gauge items in project memory; none of which this change touches).

- [ ] **Step 4: Confirm Blocker fixes landed**
- Blocker 1: `grep -n 'omega_budget_.cdm\|wa_fld_peek\|fluid_present_pfc' source/input_module.cpp` → no S8/Halofit peeks remain.
- Blocker 2: `idmdr_ethos_parameters.ini` contains `Omega_idm_dr`/`a_idm_dr`; `*unused_parameters` does not.
- Blocker 3: the `sigma_B` inline formula is gone from `input_module.cpp` (only the `PhotonsSpecies` helper remains).
- Blocker 4: `grep -n 'read_double("N_idr"\|read_double("f_idm_drmd"' species/` → only the budget resolver parses these.
- Blocker 5: `grep -rn 'a_idm_dr_\|idr_nature_' species/idr.h species/idm_dr.h` → fields live on the species; `grep -n 'a_idm_dr' source/thermodynamics.h` → gone.

- [ ] **Step 5: Final commit (if any cleanup remained)**

```bash
git add -A && git commit -m "Input-module redesign: final cleanup + verification

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>" || echo "nothing to commit"
```

---

## Notes for the implementer

- **Macro return convention:** `class_call`/`class_test` historically `return _FAILURE_` on error. If `ReadContext`/`ReadDerived` use these macros, keep them `int`-returning and call them via `class_call(ReadContext(), error_message_, error_message_)` from the constructor. Decide this in Task 4 Step 3 by inspecting the macro in `include/`/`source/common.h`; keep the whole plan consistent with that choice.
- **`SpeciesCollection::at` return type:** check whether it returns `std::unique_ptr<BaseSpecies>&` or `BaseSpecies*`; the plan assumes `*all_species_.at(key)` yields a `BaseSpecies&` (as used in `thermodynamics_module.cpp:204`). Use the same dereference idiom already in the codebase.
- **`pfc` read helpers:** the species side uses `ctx.pfc->read_double/read_int/read_string/read_list_of_doubles` (see `ultra_relativistic.cpp`, `idm_dr_idr_species.cpp`). The monolith uses `parser_read_double`/`readDoubleList`. When moving parse logic into `CreateAll`, use the species-side helpers and confirm the exact list-reader name at `input_module.cpp:921`.
- **Never claim a task done on a partial build.** Each task ends green (`make class` clean + `compare_tol.py` exit 0) before its commit.
