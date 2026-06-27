# Type-3 momentum-transfer coupling (CDM ↔ scalar field) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the Pourtsidou-Tram Type-3 pure-momentum-transfer coupling between cold dark matter and a quintessence scalar field in the CLASSpp species system (synchronous gauge), as a stress-test of the composite-species coupling machinery.

**Architecture:** Three parametrized species pieces plus a composite, mirroring the `IDM_DR_IDR` precedent: (1) the scalar field's potential becomes injectable; (2) `CDMSpecies` gains an opt-in synchronous-gauge velocity; (3) `ScalarFieldSpecies` gains a coupling parameter β that modifies only its own background and free-streaming equations; (4) a dedicated `Type3Species : CompositeSpecies` owns the bidirectional perturbation coupling via `AddCouplingDerivs` and a `StressEnergy` cross-term. β=0 / non-coupled runs stay byte-identical.

**Tech Stack:** C++17, CMake (build dir `build/cmake`), CTest hand-rolled unit tests (assert-`main`), Cython `classy` wrapper, pytest harness `python/test_class.py` with the `classyref` reference workflow.

**Reference:** Design spec `docs/superpowers/specs/2026-06-26-type3-momentum-transfer-design.md`. Equation source of truth = the old-C `github.com/ThomasTram/class` branch `CoupledQuintessence` (β is `scf_veta`; `factorveta`≡1.0 vestigial).

## Global Constraints

- **Synchronous gauge only.** β≠0 in Newtonian gauge must be **rejected** with a clear error (not silently run uncoupled). The CDM velocity flag itself works in both gauges; only the coupling *source* is synchronous-only.
- **β=0 / non-coupled byte-identical.** Every refactor/parametrization defaults to the current behavior; existing runs and `classyref` comparisons must be unchanged.
- **No bit-identical requirement for coupled physics** — validate with ~0.1% tolerance, handle `Cl^TE` zero-crossings; don't gate on committed goldens that are stale under ffast-math (use HEAD-vs-master A/B with a scale-relative metric).
- **Never `git add -A`** in this repo (in-source CMake/Xcode artifacts get swept in). Stage explicit paths only.
- **Physicality:** `class_test` / reject β ≥ 1/2 (ghost/strong-coupling). Guard the denominator `D = 3·ρ_cdm·a² − 2β·Z̄²` against `D→0` (relevant only for the β>0 branch).
- **Input key:** the coupling parameter β is read from `.ini` key `scf_veta` (absent or 0 ⟹ no coupling).
- **C++ is C++-only** (no `extern "C"` / `#ifdef __cplusplus`); headers are plain C++.
- **cclassy.pxd is auto-generated** — never hand-edit it; `generate_wrapper.py` rebuilds it from the C++ headers at build.

**Build / test commands used throughout:**
- Configure + build a target: `cmake --build build/cmake --target <t> --parallel` (CMake auto-reconfigures when `CMakeLists.txt` changes; if `build/cmake` does not exist yet, run `cmake -S . -B build/cmake` once).
- Run a C++ unit test: `ctest --test-dir build/cmake -R <test-name> --output-on-failure`.
- Build the CLI: `make class` (binary lands at repo-root `./class`).
- Build the wrapper: `make classy-pip-dev` (reuses the scikit-build tree).
- Python suite (reference gate): `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest test_class.py -k <pattern> -v` (requires a `classyref` = fresh master build installed under module name `classyref`).

---

### Task 1: Injectable scalar-field potential

Extract `V/dV/ddV` into an injectable bundle so a composite (or any client) can supply its own potential without touching `ScalarFieldSpecies`. The bundle's callables take `(phi, params)`; the species keeps owning `scf_parameters_`/`scf_tuning_index_`, so shooting is unchanged. Default = today's built-in → byte-identical.

**Files:**
- Create: `species/scalar_field_potential.h`
- Create: `species/scalar_field_potential.cpp`
- Create: `species/scalar_field_potential_test.cpp`
- Modify: `species/scalar_field.h` (add member + ctor param; declare wrappers)
- Modify: `species/scalar_field.cpp:14-58` (ctor stores potential; `V_scf/dV_scf/ddV_scf` delegate)
- Modify: `CMakeLists.txt:74-92` (add `species/scalar_field_potential.cpp` to `classpp` sources) and `CMakeLists.txt:210-217` (add `test-scf-potential`)

**Interfaces:**
- Produces:
  - `struct ScalarFieldPotential { std::function<double(double, const std::vector<double>&)> V, dV, ddV; };`
  - `ScalarFieldPotential DefaultScalarFieldPotential();`
  - `ScalarFieldSpecies` ctor gains a trailing defaulted param `ScalarFieldPotential potential = DefaultScalarFieldPotential()`.

- [ ] **Step 1: Write the failing unit test**

Create `species/scalar_field_potential_test.cpp`:
```cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "scalar_field_potential.h"

int main() {
  const std::vector<double> p = {1.22, 2.0, 0.5, 0.1};  // lambda, alpha, A, B
  const double phi = 0.7;
  const ScalarFieldPotential pot = DefaultScalarFieldPotential();

  // Default bundle reproduces V = exp(-lambda*phi)*((phi-B)^alpha + A).
  const double lambda = p[0], alpha = p[1], A = p[2], B = p[3];
  const double Ve = std::exp(-lambda * phi);
  const double V_expected = Ve * (std::pow(phi - B, alpha) + A);
  const double dV_expected =
      -lambda * Ve * (std::pow(phi - B, alpha) + A) + Ve * alpha * std::pow(phi - B, alpha - 1);
  assert(std::fabs(pot.V(phi, p) - V_expected) < 1e-12);
  assert(std::fabs(pot.dV(phi, p) - dV_expected) < 1e-12);

  // An injected pure-1EXP bundle (params = [V0, lambda]) is honored.
  ScalarFieldPotential exp1{
      [](double f, const std::vector<double>& q) { return q[0] * std::exp(-q[1] * f); },
      [](double f, const std::vector<double>& q) { return -q[1] * q[0] * std::exp(-q[1] * f); },
      [](double f, const std::vector<double>& q) {
        return q[1] * q[1] * q[0] * std::exp(-q[1] * f);
      }};
  const std::vector<double> q = {3.0, 1.22};
  assert(std::fabs(exp1.V(phi, q) - 3.0 * std::exp(-1.22 * phi)) < 1e-12);

  std::printf("scalar_field_potential tests passed\n");
  return 0;
}
```

- [ ] **Step 2: Create the header**

Create `species/scalar_field_potential.h`:
```cpp
#pragma once
#include <functional>
#include <vector>

/**
 * A scalar-field potential bundled with its first two derivatives. Each callable
 * takes (phi, params); the owning ScalarFieldSpecies passes its scf_parameters_
 * as `params`, so the shooting machinery (which tunes params[tuning_index])
 * keeps working unchanged. The potential is evaluated only on the background ODE
 * (never in the per-k perturbation loop), so std::function indirection is free.
 */
struct ScalarFieldPotential {
  std::function<double(double phi, const std::vector<double>& params)> V;
  std::function<double(double phi, const std::vector<double>& params)> dV;
  std::function<double(double phi, const std::vector<double>& params)> ddV;
};

/** The historical built-in: V = exp(-lambda*phi) * ((phi-B)^alpha + A),
 *  params = [lambda, alpha, A, B]. */
ScalarFieldPotential DefaultScalarFieldPotential();
```

- [ ] **Step 3: Create the implementation**

Create `species/scalar_field_potential.cpp`:
```cpp
#include "scalar_field_potential.h"

#include <cmath>

ScalarFieldPotential DefaultScalarFieldPotential() {
  ScalarFieldPotential p;
  p.V = [](double phi, const std::vector<double>& params) {
    const double lambda = params[0], alpha = params[1], A = params[2], B = params[3];
    return std::exp(-lambda * phi) * (std::pow(phi - B, alpha) + A);
  };
  p.dV = [](double phi, const std::vector<double>& params) {
    const double lambda = params[0], alpha = params[1], A = params[2], B = params[3];
    const double Ve = std::exp(-lambda * phi);
    const double Vp = std::pow(phi - B, alpha) + A;
    return -lambda * Ve * Vp + Ve * alpha * std::pow(phi - B, alpha - 1);
  };
  p.ddV = [](double phi, const std::vector<double>& params) {
    const double lambda = params[0], alpha = params[1], A = params[2], B = params[3];
    const double Ve = std::exp(-lambda * phi);
    const double Vp = std::pow(phi - B, alpha) + A;
    const double dVe = -lambda * Ve;
    const double dVp = alpha * std::pow(phi - B, alpha - 1);
    const double ddVe = lambda * lambda * Ve;
    const double ddVp = alpha * (alpha - 1.) * std::pow(phi - B, alpha - 2);
    return ddVe * Vp + 2 * dVe * dVp + Ve * ddVp;
  };
  return p;
}
```

- [ ] **Step 4: Wire CMake (library source + test target)**

In `CMakeLists.txt`, add to the `classpp` source list (after `species/scalar_field.cpp`'s neighbors, keep alphabetical-ish):
```cmake
  species/scalar_field_potential.cpp
```
In the test block (`CMakeLists.txt:210-213`), add:
```cmake
  add_executable(test-scf-potential species/scalar_field_potential_test.cpp)
```
and add `test-scf-potential` to the `foreach(_t IN ITEMS ...)` list at `CMakeLists.txt:214`.

- [ ] **Step 5: Run the test to verify it fails to build/link**

Run: `cmake --build build/cmake --target test-scf-potential --parallel`
Expected: FAIL — `scalar_field_potential.h`/`DefaultScalarFieldPotential` not found until Steps 2-4 are saved. (After Steps 2-4 it builds.)

- [ ] **Step 6: Run the test to verify it passes**

Run: `ctest --test-dir build/cmake -R test-scf-potential --output-on-failure`
Expected: PASS — `scalar_field_potential tests passed`.

- [ ] **Step 7: Make `ScalarFieldSpecies` delegate to the bundle**

In `species/scalar_field.h`: add `#include "scalar_field_potential.h"`, add a trailing ctor parameter `ScalarFieldPotential potential = DefaultScalarFieldPotential()`, and add a private member `ScalarFieldPotential potential_;`. Keep the private `V_scf/dV_scf/ddV_scf` declarations.

In `species/scalar_field.cpp`, change the ctor (`:14-24`) to store it (add `, potential_(std::move(potential))` to the initializer list and the matching parameter), and replace the three method bodies (`:26-58`) with delegations:
```cpp
double ScalarFieldSpecies::V_scf(double phi) const { return potential_.V(phi, scf_parameters_); }
double ScalarFieldSpecies::dV_scf(double phi) const { return potential_.dV(phi, scf_parameters_); }
double ScalarFieldSpecies::ddV_scf(double phi) const { return potential_.ddV(phi, scf_parameters_); }
```

- [ ] **Step 8: Rebuild and run the byte-identical gates**

Run: `cmake --build build/cmake --target test-scf-potential test-species-types --parallel && ctest --test-dir build/cmake --output-on-failure`
Expected: PASS (all C++ unit tests).
Run: `make class` then `./class explanatory.ini` — Expected: completes without error.
Run (reference gate, if `classyref` installed): `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest test_class.py -k scalar_field -v`
Expected: PASS — the default potential reproduces master byte-for-byte.

- [ ] **Step 9: Commit**

```bash
git add species/scalar_field_potential.h species/scalar_field_potential.cpp \
        species/scalar_field_potential_test.cpp species/scalar_field.h species/scalar_field.cpp \
        CMakeLists.txt
git commit -m "scf: make the scalar-field potential injectable (default = built-in)"
```

---

### Task 2: Opt-in coupled CDM velocity (synchronous gauge)

Give `CDMSpecies` a construction flag that, when set, makes `θ_cdm` a dynamical variable in synchronous gauge (free-streaming only — the momentum-transfer source is added later by the composite). Default off ⟹ byte-identical. This is reusable infrastructure for any X-coupled-to-CDM model.

**Files:**
- Create: `species/cdm_coupled_test.cpp`
- Modify: `species/cdm.h:18-31,91-98` (ctor flag + member)
- Modify: `species/cdm.cpp:8-9,36-43,45-63,65-82,101-123,178-213` (registration, derivs, sources, IC, print)
- Modify: `CMakeLists.txt` test block (add `test-cdm-coupled`)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `CDMSpecies(const background& pba, double omega0_cdm, bool coupled = false)` and accessor `bool coupled() const;`. When `coupled` and gauge is synchronous, `PerturbLayout::idx_theta >= 0` and `θ_cdm` evolves as `−(a'/a)·θ_cdm`.

- [ ] **Step 1: Write the failing unit test**

Create `species/cdm_coupled_test.cpp`:
```cpp
#include <cassert>
#include <cstdio>
#include <memory>

#include "background.h"
#include "cdm.h"
#include "perturbations.h"  // possible_gauges

int main() {
  background pba{};
  pba.H0 = 1e-4;
  const int sync = static_cast<int>(possible_gauges::synchronous);

  // Uncoupled CDM in synchronous gauge: no theta variable.
  {
    CDMSpecies cdm(pba, 0.12, /*coupled=*/false);
    auto layout = cdm.CreatePerturbLayout();
    int index_pt = 0;
    cdm.RegisterPerturbationIndices(*layout, nullptr, nullptr, index_pt, nullptr, sync);
    const auto& l = static_cast<const CDMSpecies::PerturbLayout&>(*layout);
    assert(l.idx_delta == 0);
    assert(l.idx_theta == -1);  // theta_cdm = 0 by gauge choice
    assert(index_pt == 1);
  }

  // Coupled CDM in synchronous gauge: theta is registered.
  {
    CDMSpecies cdm(pba, 0.12, /*coupled=*/true);
    auto layout = cdm.CreatePerturbLayout();
    int index_pt = 0;
    cdm.RegisterPerturbationIndices(*layout, nullptr, nullptr, index_pt, nullptr, sync);
    const auto& l = static_cast<const CDMSpecies::PerturbLayout&>(*layout);
    assert(l.idx_delta == 0);
    assert(l.idx_theta == 1);
    assert(index_pt == 2);
  }

  std::printf("cdm coupled-velocity tests passed\n");
  return 0;
}
```

- [ ] **Step 2: Wire the test target and run to verify it fails**

In `CMakeLists.txt` test block add `add_executable(test-cdm-coupled species/cdm_coupled_test.cpp)` and add `test-cdm-coupled` to the `foreach` list.
Run: `cmake --build build/cmake --target test-cdm-coupled --parallel`
Expected: FAIL — `CDMSpecies` has no 3-arg constructor.

- [ ] **Step 3: Add the flag to the header**

In `species/cdm.h`, change the ctor to `explicit CDMSpecies(const background& pba, double omega0_cdm, bool coupled = false);`, add `bool coupled() const { return coupled_; }`, and add private member `bool coupled_ = false;`.

- [ ] **Step 4: Store the flag and flip the synchronous branches**

In `species/cdm.cpp`:

Ctor (`:8-9`):
```cpp
CDMSpecies::CDMSpecies(const background& pba, double omega0_cdm, bool coupled)
    : BaseSpecies("CDM", EnergyType::Matter), Omega0_cdm_(omega0_cdm), H0_(pba.H0),
      coupled_(coupled) {}
```

`RegisterPerturbationIndices` (`:45-63`) — allocate θ in synchronous too when coupled:
```cpp
  layout.idx_delta = index_pt;
  ++index_pt;

  layout.idx_theta = -1;
  const bool newtonian = (gauge == static_cast<int>(possible_gauges::newtonian));
  if (newtonian || coupled_) {
    layout.idx_theta = index_pt;
    ++index_pt;
  }
```

`PerturbDerivs` (`:65-82`) — coupled synchronous CDM carries θ (free-streaming only):
```cpp
  if (gauge == possible_gauges::newtonian) {
    dy[layout.idx_delta] = -(y[layout.idx_theta] + ctx.metric_continuity);
    dy[layout.idx_theta] = -ctx.a_prime_over_a * y[layout.idx_theta] + ctx.metric_euler;
  }
  else if (coupled_) { /* synchronous, coupled: theta sourced by composite */
    dy[layout.idx_delta] = -(y[layout.idx_theta] + ctx.metric_continuity);
    dy[layout.idx_theta] = -ctx.a_prime_over_a * y[layout.idx_theta] + ctx.metric_euler;
  }
  else { /* synchronous: theta_cdm = 0 by gauge choice */
    dy[layout.idx_delta] = -ctx.metric_continuity;
  }
```

`RegisterTransferSourceIndices` (`:36-43`) — register the θ source in synchronous too when coupled:
```cpp
  class_define_index(index_tp_delta_, ctx.wants_density, index_tp, 1);
  class_define_index(index_tp_theta_,
                     ctx.wants_velocity &&
                         (coupled_ ||
                          ctx.gauge != static_cast<int>(possible_gauges::synchronous)),
                     index_tp,
                     1);
```

`ApplyInitialConditions` (`:101-123`) — initialize θ_cdm to 0 (adiabatic) when present. After the existing `idx_delta` writes, add:
```cpp
  if (layout.idx_theta >= 0)
    y[layout.idx_theta] = 0.;
```

`PrintVariables` (`:194-199`) — output the real θ when coupled. Replace the forced-zero synchronous branch:
```cpp
    if (ppt->gauge == possible_gauges::synchronous && !coupled_) {
      theta_cdm = 0.;
    }
    else {
      theta_cdm = (layout.idx_theta >= 0) ? y[layout.idx_theta] : 0.;
    }
```

- [ ] **Step 5: Run the unit test to verify it passes**

Run: `cmake --build build/cmake --target test-cdm-coupled --parallel && ctest --test-dir build/cmake -R test-cdm-coupled --output-on-failure`
Expected: PASS — `cdm coupled-velocity tests passed`.

- [ ] **Step 6: Run the byte-identical gate (default off)**

Run: `cmake --build build/cmake --target class --parallel && ./class explanatory.ini`
Expected: completes; the `coupled=false` default leaves all paths unchanged.
Run (if `classyref` installed): `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest test_class.py -v`
Expected: PASS — no regression for any existing scenario.

- [ ] **Step 7: Commit**

```bash
git add species/cdm.h species/cdm.cpp species/cdm_coupled_test.cpp CMakeLists.txt
git commit -m "cdm: opt-in coupled velocity (theta_cdm in synchronous gauge, default off)"
```

---

### Task 3: Scalar-field coupling parameter β (own equations only)

Add β to `ScalarFieldSpecies` (default 0 → byte-identical). β modifies only the field's *own* equations: background `ρ_φ/P_φ/P'_φ/KG`, its free-streaming KG (`k² → k²/(1−2β)`), and the kinetic part of its `δρ_φ/δp_φ`. The momentum cross-coupling is **not** here (that is the composite, Task 5).

**Files:**
- Modify: `species/scalar_field.h` (ctor param + member + `beta()` accessor)
- Modify: `species/scalar_field.cpp:101-141` (background), `:206-210` (KG), `:405-440` (StressEnergy)
- Modify: `species/scalar_field_potential_test.cpp` (extend with a β background assertion) — or add `species/scalar_field_beta_test.cpp`

**Interfaces:**
- Consumes: `ScalarFieldPotential` ctor param from Task 1.
- Produces: `ScalarFieldSpecies` ctor gains a trailing `double beta = 0.` (after `potential`); accessor `double beta() const;`. The composite (Task 4/5) reads β via this accessor.

- [ ] **Step 1: Write the failing unit test**

Create `species/scalar_field_beta_test.cpp`:
```cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "background.h"
#include "scalar_field.h"

int main() {
  background pba{};
  const double beta = -2.0;
  const std::vector<double> params = {1.22, 0.0, 0.0, 0.0};  // pure exp: V = exp(-1.22 phi)
  ScalarFieldSpecies scf(pba, /*omega0=*/0.7, params, /*tuning=*/0, /*attractor=*/true,
                         /*phi_ini=*/1., /*phi_prime_ini=*/1.,
                         DefaultScalarFieldPotential(), beta);
  assert(std::fabs(scf.beta() - beta) < 1e-15);

  // Register the background/integration indices into fresh counters.
  int index_bg = 0, index_bi = 0;
  scf.RegisterBackgroundIndices(index_bg);
  scf.RegisterIntegrationIndices(index_bi);
  std::vector<double> pvecback(index_bg, 0.), pvecback_B(index_bi, 0.);

  const double a = 0.5, phi = 0.3, phi_prime = 0.4;
  pvecback_B[scf.bi_phi_index()] = phi;
  pvecback_B[scf.bi_phi_prime_index()] = phi_prime;
  scf.ComputeBackground(a, pvecback_B.data(), pvecback.data());

  const double V = std::exp(-1.22 * phi);
  const double rho_expected = ((1. - 2. * beta) * phi_prime * phi_prime / (2. * a * a) + V) / 3.;
  // rho is stored at the scalar field's background rho slot; Rho() reads it.
  assert(std::fabs(scf.Rho(pvecback.data()) - rho_expected) < 1e-12);

  std::printf("scalar_field beta background test passed\n");
  return 0;
}
```

- [ ] **Step 2: Wire the test target and run to verify it fails**

In `CMakeLists.txt` test block add `add_executable(test-scf-beta species/scalar_field_beta_test.cpp)` and add `test-scf-beta` to the `foreach` list.
Run: `cmake --build build/cmake --target test-scf-beta --parallel`
Expected: FAIL — ctor has no β parameter / no `beta()` accessor.

- [ ] **Step 3: Add β to the header**

In `species/scalar_field.h`: add a trailing ctor parameter `double beta = 0.` (after `ScalarFieldPotential potential = DefaultScalarFieldPotential()`), add `double beta() const { return beta_; }`, and add private member `double beta_ = 0.;`.

- [ ] **Step 4: Apply β to the background**

In `species/scalar_field.cpp`:

Ctor: store `, beta_(beta)` and add the parameter.

`ComputeBackground` (`:109-110`) — kinetic term gets `(1−2β)`:
```cpp
  pvecback[index_bg_rho_] = ((1. - 2. * beta_) * phi_prime * phi_prime / (2. * a * a) + V_scf(phi)) / 3.;
  pvecback[index_bg_p_]   = ((1. - 2. * beta_) * phi_prime * phi_prime / (2. * a * a) - V_scf(phi)) / 3.;
```

`BackgroundDerivs` (`:124`) — `1/(1−2β)` on `dV`:
```cpp
  dy[index_bi_phi_prime_scf_] = -a * (2. * H * phi_prime + a * dV_scf(phi) / (1. - 2. * beta_));
```

`PPrime` (`:133`) — `(1−2β)` on the kinetic term only:
```cpp
  return phi_prime * (-(1. - 2. * beta_) * phi_prime * H / a - 2. / 3. * dV_scf(phi));
```

- [ ] **Step 5: Apply β to the scalar field's own perturbations**

`PerturbDerivs` synchronous branch (`:208-210`) — `k² → k²/(1−2β)`:
```cpp
    dy[layout.idx_phi]       = y[layout.idx_phi_prime];
    dy[layout.idx_phi_prime] = -2. * a_prime_over_a * y[layout.idx_phi_prime] -
                               metric_continuity * phi_prime_bg -
                               (k2 / (1. - 2. * beta_) + a2 * ddV_bg) * y[layout.idx_phi];
```
(The Newtonian branch stays as-is; β≠0 in Newtonian gauge is rejected upstream by the composite, Task 4.)

`StressEnergy` (`:433`) — kinetic part of `δρ_φ` gets `(1−2β)`. Replace the `delta_rho` line:
```cpp
  double delta_rho = (1. / 3.) * ((1. - 2. * beta_) / a2 * phi_prime * delta_phi_prime + dV * y[layout.idx_phi]);
```
and the corresponding `delta_p` kinetic term (the `+ (1/3)(1/a2) phi_prime delta_phi_prime` part) likewise gets `(1−2β)`. (Leave `rho_plus_p_theta` at `:440` unchanged — the θ_cdm cross-term is the composite's job.)

- [ ] **Step 6: Run the unit test to verify it passes**

Run: `cmake --build build/cmake --target test-scf-beta --parallel && ctest --test-dir build/cmake -R test-scf-beta --output-on-failure`
Expected: PASS — `scalar_field beta background test passed`.

- [ ] **Step 7: Run the byte-identical gate (β=0 default)**

Run: `cmake --build build/cmake --target test-scf-beta test-scf-potential --parallel && ctest --test-dir build/cmake --output-on-failure`
Expected: PASS.
Run (if `classyref` installed): `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest test_class.py -k scalar_field -v`
Expected: PASS — β defaults to 0, so existing scalar-field runs are byte-identical.

- [ ] **Step 8: Commit**

```bash
git add species/scalar_field.h species/scalar_field.cpp species/scalar_field_beta_test.cpp CMakeLists.txt
git commit -m "scf: add coupling parameter beta to the field's own background + free-streaming equations"
```

---

### Task 4: `Type3Species` composite — scaffold, factory, input wiring, gauge guard

Create the composite that owns the two children (coupled CDM + β scalar field), register it in the factory, route input (`scf_veta`, `Omega_scf`, `Omega_cdm`), suppress the standalone CDM/scalar-field when the coupling is active, and reject Newtonian gauge. **No coupling terms yet** — `AddCouplingDerivs`/`StressEnergy` cross-term land in Task 5. After this task, a β≠0 run computes (modified-background scalar field + uncoupled CDM whose θ stays 0).

**Files:**
- Create: `species/type3_species.h`, `species/type3_species.cpp`
- Modify: `species/all_species.h:31-59` (include + factory row)
- Modify: `species/species_type_name_test.cpp:25-48` (pin the new `kTypeName`)
- Modify: `species/scalar_field.cpp` `CreateAll` (skip when `scf_veta` active)
- Modify: `species/cdm.cpp` `CreateAll` (skip when `scf_veta` active)
- Modify: `CMakeLists.txt:74-92` (add `species/type3_species.cpp`)
- Modify: `python/test_class.py` (add three integration tests)

**Interfaces:**
- Consumes: `CDMSpecies(pba, omega0, coupled=true)` (Task 2); `ScalarFieldSpecies(..., potential, beta)` (Tasks 1+3).
- Produces:
  - `class Type3Species : public CompositeSpecies` with `static constexpr std::string_view kTypeName = "cdm_scf_momentum";`
  - nested `struct PerturbLayout : BaseSpecies::PerturbLayout { CDMSpecies::PerturbLayout cdm; ScalarFieldSpecies::PerturbLayout scf; };`
  - `static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);`
  - accessors `CDMSpecies& cdm(); ScalarFieldSpecies& scf();` and `double beta() const { return scf_->beta(); }`
  - protected `void AddCouplingDerivs(...) const override;` (empty in this task)

- [ ] **Step 1: Write the header**

Create `species/type3_species.h`:
```cpp
#pragma once
#include <string_view>
#include <vector>

#include "background.h"
#include "cdm.h"
#include "composite_species.h"
#include "scalar_field.h"
#include "species_build_context.h"

/**
 * Type3Species: composite for the Pourtsidou-Tram Type-3 pure-momentum-transfer
 * coupling between cold dark matter and a quintessence scalar field
 * (arXiv:1604.04222). Children handle free-streaming; AddCouplingDerivs (Task 5)
 * adds the synchronous-gauge momentum exchange, and StressEnergy adds the
 * -(2 beta/3) Zbar^2 theta_cdm cross-term. Synchronous gauge only.
 */
class Type3Species : public CompositeSpecies {
 public:
  static constexpr std::string_view kTypeName = "cdm_scf_momentum";

  struct PerturbLayout : BaseSpecies::PerturbLayout {
    CDMSpecies::PerturbLayout cdm;
    ScalarFieldSpecies::PerturbLayout scf;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    return std::make_unique<PerturbLayout>();
  }

  Type3Species(const background& pba,
               double omega0_cdm,
               std::unique_ptr<ScalarFieldSpecies> scf);

  CDMSpecies& cdm() { return *cdm_; }
  ScalarFieldSpecies& scf() { return *scf_; }
  double beta() const { return scf_->beta(); }

  void RegisterPerturbationIndices(BaseSpecies::PerturbLayout& layout,
                                   perturb_vector* pv,
                                   const precision* ppr,
                                   int& index_pt,
                                   const perturb_workspace* ppw,
                                   int gauge) override;
  void PerturbDerivs(const BaseSpecies::PerturbLayout& layout,
                     double tau,
                     const double* y,
                     double* dy,
                     const perturb_parameters_and_workspace& ppaw) const override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  StressEnergyContribution StressEnergy(const BaseSpecies::PerturbLayout& layout,
                                        const perturb_vector* pv,
                                        const double* y,
                                        const double* pvecback,
                                        const perturb_workspace* ppw) const override;
  void FillSources(const BaseSpecies::PerturbLayout& layout,
                   const double* y,
                   const double* dy,
                   PerturbSourceContext& ctx) const override;
  void RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) override;
  void PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& layout,
                                     double* y,
                                     const PerturbIcContext& ctx) override;
  void CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_layout,
                                     const BaseSpecies::PerturbLayout& new_layout,
                                     const double* old_y,
                                     double* new_y,
                                     const PerturbSwitchContext& ctx) const override;
  void WriteOutputColumns(PerturbColumnWriter& writer,
                          const PerturbationsModule& mod,
                          file_format fmt,
                          TransferColumnSection section = TransferColumnSection::all) const override;
  void PrintVariables(PerturbColumnWriter& writer,
                      double tau,
                      const double* y,
                      const PerturbationsModule& mod,
                      const perturb_workspace* ppw) const override;

  std::vector<ShootingTarget> GetShootingTargets() const override { return scf_->GetShootingTargets(); }
  void ComputeShootingGuess(const SpeciesBuildContext& ctx,
                            std::vector<double>& guess,
                            std::vector<double>& dxdy) const override {
    scf_->ComputeShootingGuess(ctx, guess, dxdy);
  }
  double ComputeShootingResidual(const ShootingResidualContext& ctx,
                                 const ShootingTarget& target) const override {
    return scf_->ComputeShootingResidual(ctx, target);
  }

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

 protected:
  void AddCouplingDerivs(double tau,
                         const double* y,
                         double* dy,
                         const perturb_parameters_and_workspace& ppaw) const override;

 private:
  CDMSpecies* cdm_ = nullptr;
  ScalarFieldSpecies* scf_ = nullptr;
};
```

- [ ] **Step 2: Write the implementation (scaffold; AddCouplingDerivs empty)**

Create `species/type3_species.cpp`:
```cpp
#include "type3_species.h"

#include <stdexcept>

#include "perturbations.h"
#include "perturbations_module.h"

Type3Species::Type3Species(const background& pba,
                           double omega0_cdm,
                           std::unique_ptr<ScalarFieldSpecies> scf)
    : CompositeSpecies("CDM_SCF_Momentum", EnergyType::Other) {
  auto cdm = std::make_unique<CDMSpecies>(pba, omega0_cdm, /*coupled=*/true);
  cdm_ = cdm.get();
  scf_ = scf.get();
  children_.push_back(std::move(cdm));
  children_.push_back(std::move(scf));
}

void Type3Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                               perturb_vector* pv,
                                               const precision* ppr,
                                               int& index_pt,
                                               const perturb_workspace* ppw,
                                               int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  cdm_->RegisterPerturbationIndices(my.cdm, pv, ppr, index_pt, ppw, gauge);
  scf_->RegisterPerturbationIndices(my.scf, pv, ppr, index_pt, ppw, gauge);
}

void Type3Species::PerturbDerivs(const BaseSpecies::PerturbLayout& base,
                                 double tau,
                                 const double* y,
                                 double* dy,
                                 const perturb_parameters_and_workspace& ppaw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  cdm_->PerturbDerivs(my.cdm, tau, y, dy, ppaw);
  scf_->PerturbDerivs(my.scf, tau, y, dy, ppaw);
  AddCouplingDerivs(tau, y, dy, ppaw);
}

void Type3Species::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                          double* y,
                                          const PerturbIcContext& ctx) {
  auto& my = static_cast<PerturbLayout&>(base);
  cdm_->ApplyInitialConditions(my.cdm, y, ctx);
  scf_->ApplyInitialConditions(my.scf, y, ctx);
}

BaseSpecies::StressEnergyContribution Type3Species::StressEnergy(
    const BaseSpecies::PerturbLayout& base,
    const perturb_vector* pv,
    const double* y,
    const double* pvecback,
    const perturb_workspace* ppw) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  StressEnergyContribution se = cdm_->StressEnergy(my.cdm, pv, y, pvecback, ppw);
  se += scf_->StressEnergy(my.scf, pv, y, pvecback, ppw);
  return se;  // cross-term added in Task 5
}

void Type3Species::FillSources(const BaseSpecies::PerturbLayout& base,
                               const double* y,
                               const double* dy,
                               PerturbSourceContext& ctx) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  cdm_->FillSources(my.cdm, y, dy, ctx);
  scf_->FillSources(my.scf, y, dy, ctx);
}

void Type3Species::RegisterTransferSourceIndices(int& index_tp, const SourceRequestContext& ctx) {
  cdm_->RegisterTransferSourceIndices(index_tp, ctx);
  scf_->RegisterTransferSourceIndices(index_tp, ctx);
}

void Type3Species::PerturbSynchronousToNewtonian(const BaseSpecies::PerturbLayout& base,
                                                 double* y,
                                                 const PerturbIcContext& ctx) {
  auto& my = static_cast<PerturbLayout&>(base);
  cdm_->PerturbSynchronousToNewtonian(my.cdm, y, ctx);
  scf_->PerturbSynchronousToNewtonian(my.scf, y, ctx);
}

void Type3Species::CopyPerturbationsAcrossSwitch(const BaseSpecies::PerturbLayout& old_base,
                                                 const BaseSpecies::PerturbLayout& new_base,
                                                 const double* old_y,
                                                 double* new_y,
                                                 const PerturbSwitchContext& ctx) const {
  const auto& old_l = static_cast<const PerturbLayout&>(old_base);
  const auto& new_l = static_cast<const PerturbLayout&>(new_base);
  cdm_->CopyPerturbationsAcrossSwitch(old_l.cdm, new_l.cdm, old_y, new_y, ctx);
  scf_->CopyPerturbationsAcrossSwitch(old_l.scf, new_l.scf, old_y, new_y, ctx);
}

void Type3Species::WriteOutputColumns(PerturbColumnWriter& writer,
                                      const PerturbationsModule& mod,
                                      file_format fmt,
                                      TransferColumnSection section) const {
  cdm_->WriteOutputColumns(writer, mod, fmt, section);
  scf_->WriteOutputColumns(writer, mod, fmt, section);
}

void Type3Species::PrintVariables(PerturbColumnWriter& writer,
                                  double tau,
                                  const double* y,
                                  const PerturbationsModule& mod,
                                  const perturb_workspace* ppw) const {
  cdm_->PrintVariables(writer, tau, y, mod, ppw);
  scf_->PrintVariables(writer, tau, y, mod, ppw);
}

void Type3Species::AddCouplingDerivs(double /*tau*/,
                                     const double* /*y*/,
                                     double* /*dy*/,
                                     const perturb_parameters_and_workspace& /*ppaw*/) const {
  // Task 5 fills this in.
}

std::vector<Named> Type3Species::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // Coupling active iff scf_veta is present and nonzero.
  auto veta_opt = ctx.pfc->get<double>("scf_veta");
  if (!veta_opt.has_value() || *veta_opt == 0.)
    return result;
  const double beta = *veta_opt;

  if (beta >= 0.5)
    throw std::invalid_argument(
        "scf_veta (beta) must be < 1/2: beta >= 1/2 has a ghost / strong-coupling pathology.");

  // Synchronous gauge only.
  if (auto gauge = ctx.pfc->get<std::string>("gauge")) {
    if (gauge->find("newton") != std::string::npos || gauge->find("Newton") != std::string::npos)
      throw std::invalid_argument(
          "Type-3 (scf_veta) coupling is implemented in synchronous gauge only.");
  }

  // CDM Omega from the resolved coupled-species budget (same source CDMSpecies uses).
  const double omega0_cdm =
      (ctx.omega_budget && ctx.omega_budget->cdm.has_value()) ? *ctx.omega_budget->cdm : 0.;

  // Build the scalar-field child via the shared scalar-field input path, then
  // hand its single instance to the composite. ScalarFieldSpecies::CreateAll
  // returns at most one "ScalarField" entry; it must be present here.
  std::vector<Named> scf_built = ScalarFieldSpecies::CreateAllForComposite(ctx, beta);
  if (scf_built.empty())
    throw std::invalid_argument(
        "scf_veta is set but no scalar field was configured (set Omega_scf and scf_parameters).");
  auto scf_ptr = std::unique_ptr<ScalarFieldSpecies>(
      static_cast<ScalarFieldSpecies*>(scf_built.front().species.release()));

  result.push_back({"CDM_SCF_Momentum",
                    std::make_unique<Type3Species>(*ctx.pba, omega0_cdm, std::move(scf_ptr))});
  return result;
}
```

*Note:* `ScalarFieldSpecies::CreateAllForComposite(ctx, beta)` is a thin refactor of the existing `ScalarFieldSpecies::CreateAll` body that (a) takes β and threads it into the constructed `ScalarFieldSpecies`, and (b) returns the built species (with its shooting state set) for the composite to adopt. Extract the shared body in Step 4.

- [ ] **Step 3: Register the composite in the factory + pin its type name**

In `species/all_species.h`: add `#include "type3_species.h"` and add a factory row
`SpeciesFactoryEntry{Type3Species::kTypeName, &Type3Species::CreateAll},` to `kAllSpeciesFactories`.

In `species/species_type_name_test.cpp`: add `assert(Type3Species::kTypeName == "cdm_scf_momentum");` and add `"cdm_scf_momentum"` to the `expected` set.

In `CMakeLists.txt`: add `species/type3_species.cpp` to the `classpp` source list.

- [ ] **Step 4: Suppress standalone CDM / scalar field when coupling is active**

In `species/scalar_field.cpp` `CreateAll` (top of the function): after parsing, add
```cpp
  if (auto veta = ctx.pfc->get<double>("scf_veta"); veta.has_value() && *veta != 0.)
    return result;  // the Type3 composite owns the scalar field
```
Then extract the existing build body into a new
`std::vector<Named> ScalarFieldSpecies::CreateAllForComposite(const SpeciesBuildContext& ctx, double beta)` that constructs the `ScalarFieldSpecies` with the trailing `beta` argument (and the 1EXP/default potential), preserving the shooting setup (`needs_shooting_`, `shooting_target_`, `ComputeShootingGuess`). Declare it in `species/scalar_field.h`. `CreateAll` (standalone) calls it with `beta = 0.`.

In `species/cdm.cpp` `CreateAll`: after the budget read, add
```cpp
  if (auto veta = ctx.pfc->get<double>("scf_veta"); veta.has_value() && *veta != 0.)
    return result;  // the Type3 composite owns the CDM
```

- [ ] **Step 5: Build and run the type-name unit test**

Run: `cmake --build build/cmake --target test-species-types --parallel && ctest --test-dir build/cmake -R test-species-types --output-on-failure`
Expected: PASS — the registry now contains `cdm_scf_momentum`.

- [ ] **Step 6: Write the failing Python integration tests**

In `python/test_class.py`, add to the `TestClass` body:
```python
    def test_type3_synchronous_computes(self):
        scenario = {
            'output': 'tCl mPk',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '1.22, 0, 0, 0',
            'scf_veta': -0.5,
        }
        self._assert_compute_succeeds(scenario)

    def test_type3_newtonian_rejected(self):
        scenario = {
            'output': 'tCl',
            'gauge': 'newtonian',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '1.22, 0, 0, 0',
            'scf_veta': -0.5,
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        with self.assertRaises(Exception):
            self.cosmo.compute()

    def test_type3_beta_zero_unset_is_plain_scf(self):
        # scf_veta absent => no composite => plain CDM + scalar field.
        scenario = {
            'output': 'tCl',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '1.22, 0, 0, 0',
        }
        self._assert_compute_succeeds(scenario)
```

- [ ] **Step 7: Build the wrapper and run the integration tests**

Run: `make classy-pip-dev`
Run: `cd python && python -m pytest test_class.py -k type3 -v`
Expected: PASS — `test_type3_synchronous_computes` and `test_type3_beta_zero_unset_is_plain_scf` compute; `test_type3_newtonian_rejected` raises.

- [ ] **Step 8: Byte-identical gate for non-coupled runs**

Run: `cd python && COMPARE_OUTPUT_REF=1 TEST_LEVEL=2 python -m pytest test_class.py -v`
Expected: PASS — adding the (inactive) factory does not change any existing scenario.

- [ ] **Step 9: Commit**

```bash
git add species/type3_species.h species/type3_species.cpp species/all_species.h \
        species/species_type_name_test.cpp species/scalar_field.h species/scalar_field.cpp \
        species/cdm.cpp CMakeLists.txt python/test_class.py
git commit -m "type3: composite scaffold + factory/input wiring + synchronous-only guard"
```

---

> **Sequencing note (added 2026-06-27):** Task 4 surfaced finding **I-2** (a matter+DE
> composite pollutes the cold `delta_m` → P(k) ~10× low). Its fix is a separate plan —
> `docs/superpowers/plans/2026-06-27-composite-child-iteration-and-matter-tally.md` —
> whose tasks run as ledger **Tasks 5-7** (cached classification → generic composite
> child iteration → per-child `TallyStressEnergy` + dropping the `ρ−3P` proxy) **before**
> the coupling below. So this "Task 5" is **ledger Task 8** and "Task 6" (Fig-2
> validation) is **ledger Task 9**; the I-2 plan's verification sweep is **ledger Task 10
> (last)**. See `.superpowers/sdd/progress.md`. The `StressEnergy` cross-term below builds
> on the per-child structure from ledger Tasks 6-7 (Type3 re-overrides `StressEnergy` =
> generic child sum + the `−(2β/3)Z̄²θ_cdm` term).

### Task 5: The coupling — `AddCouplingDerivs` + `StressEnergy` cross-term

Fill in the synchronous-gauge momentum exchange: the φ-KG `θ_cdm` source, the CDM-Euler momentum-transfer source, and the `−(2β/3)Z̄²θ_cdm` cross-term in the composite's `StressEnergy`. Conventions: `δφ = y[scf.idx_phi]`, `δφ' = y[scf.idx_phi_prime]`, `θc = y[cdm.idx_theta]`, `Z̄ = −φ'_bg/a`, `D = 3·ρ_cdm·a² − 2β·Z̄²`.

**Files:**
- Modify: `species/type3_species.cpp` (`AddCouplingDerivs`, `StressEnergy`)
- Create: `species/type3_coupling_test.cpp` (β=0 ⟹ coupling vanishes)
- Modify: `CMakeLists.txt` test block (`test-type3-coupling`)
- Modify: `python/test_class.py` (β≠0 drives θ_cdm; suppression sanity)

**Interfaces:**
- Consumes: `Type3Species::PerturbLayout` (Task 4), `ScalarFieldSpecies` background indices (`index_bg_phi_prime_scf_`, `index_bg_dV_scf_`) via the scalar-field child, `CDMSpecies::index_bg_rho_cdm_` via `cdm_->Rho()`.
- Produces: a coupling-complete `Type3Species`.

- [ ] **Step 1: Write the failing unit test (coupling vanishes at β=0)**

Create `species/type3_coupling_test.cpp`. It builds a `Type3Species` with β=0, lays out a small `y/dy`, runs `PerturbDerivs`, and asserts the scalar-field-`phi_prime` and `theta_cdm` derivatives equal what the children alone produce (i.e. `AddCouplingDerivs` added nothing):
```cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "background.h"
#include "scalar_field.h"
#include "type3_species.h"
// NOTE: PerturbDerivs needs a populated perturb_parameters_and_workspace; this
// test exercises the pure-arithmetic guarantee that AddCouplingDerivs is a no-op
// at beta=0 by calling it directly through a minimal workspace built by the
// helper make_min_ppaw() (see species/test_support.h, added in this step).
```
Add a tiny shared test helper `species/test_support.h` exposing `make_min_ppaw(a, H, rho_cdm, phi_prime_bg, dV_bg, metric_continuity)` that fills only the `PerturbScalarContext` fields and `pvecback` slots the coupling reads. Assert: with β=0, `dy[scf.idx_phi_prime]` and `dy[cdm.idx_theta]` are unchanged by `AddCouplingDerivs` (compare a run with β=0 against the children-only derivatives); with β=−0.5, `dy[cdm.idx_theta]` changes by exactly the §5.3 expression evaluated by hand for the chosen inputs.

*(If wiring a full `perturb_parameters_and_workspace` in a unit test proves too heavy, downgrade this to: assert the closed-form coupling expression — factored into a free function `Type3CouplingDeltaThetaCdm(beta, k2, a_prime_over_a, rho_cdm, a2, Zbar, dV, phi, phi_prime, theta_cdm)` and `Type3CouplingDeltaPhiPrime(beta, phi_prime_bg, theta_cdm)` — returns 0 at β=0 and the hand-computed value at β=−0.5. Implement the coupling in those free functions and call them from `AddCouplingDerivs`; this keeps the physics unit-testable without a workspace.)*

- [ ] **Step 2: Wire the test target and run to verify it fails**

In `CMakeLists.txt` add `add_executable(test-type3-coupling species/type3_coupling_test.cpp)` and add it to the `foreach` list.
Run: `cmake --build build/cmake --target test-type3-coupling --parallel`
Expected: FAIL — the coupling free functions / behavior do not exist yet.

- [ ] **Step 3: Implement the coupling free functions**

Add to `species/type3_species.cpp` (and declare both in `species/type3_species.h`). `Zbar = -phi_prime_bg / a`; `D = 3*rho_cdm*a2 - 2*beta*Zbar^2`. Both return 0 at β=0 (every term carries an explicit `beta` factor):
```cpp
// phi-KG momentum source (added to dy[scf.idx_phi_prime]).
double Type3CouplingDeltaPhiPrime(double beta, double phi_prime_bg, double theta_cdm) {
  return 2. * (beta / (1. - 2. * beta)) * phi_prime_bg * theta_cdm;
}

// CDM-Euler momentum-transfer source (added to dy[cdm.idx_theta]); spec §5.3.
double Type3CouplingDeltaThetaCdm(double beta, double k2, double a_prime_over_a, double rho_cdm,
                                  double a2, double Zbar, double dV, double phi, double phi_prime,
                                  double phi_prime_bg, double theta_cdm) {
  const double one_m_2b = 1. - 2. * beta;
  const double D = 3. * rho_cdm * a2 - 2. * beta * Zbar * Zbar;
  const double term1 = -k2 * 2. * beta * (one_m_2b * Zbar * Zbar * phi_prime + dV * phi) /
                       (one_m_2b * D);
  const double term2 = 4. * beta * dV * phi_prime_bg *
                       (k2 * phi / phi_prime_bg - 2. * beta * theta_cdm) / (one_m_2b * D);
  const double term3 = -(6. * beta * a_prime_over_a * Zbar * Zbar + 4. * beta * dV * phi_prime_bg) *
                       theta_cdm / D;
  return term1 + term2 + term3;
}
```

- [ ] **Step 4: Call them from `AddCouplingDerivs`**

Replace the empty `AddCouplingDerivs` body with (mirrors `IDM_DR_IDR_Species::AddCouplingDerivs`'s workspace access):
```cpp
void Type3Species::AddCouplingDerivs(double /*tau*/,
                                     const double* y,
                                     double* dy,
                                     const perturb_parameters_and_workspace& ppaw) const {
  const double beta = scf_->beta();
  if (beta == 0.)
    return;

  const perturb_workspace* ppw = ppaw.ppw;
  const perturb_vector* pv = ppw->pv.get();
  const PerturbScalarContext& ctx = ppw->scalar_ctx;
  const auto& lay = static_cast<const PerturbLayout&>(*pv->species_layouts[collection_index_]);

  const double* pvecback = ppw->pvecback.data();
  const double a = ctx.a, a2 = ctx.a2, a_prime_over_a = ctx.a_prime_over_a, k2 = ctx.k2;
  const double phi_prime_bg = pvecback[scf_->index_bg_phi_prime_scf()];
  const double dV = pvecback[scf_->index_bg_dV_scf()];
  const double rho_cdm = cdm_->Rho(pvecback);
  const double Zbar = -phi_prime_bg / a;

  const double phi = y[lay.scf.idx_phi];
  const double phi_prime = y[lay.scf.idx_phi_prime];
  const double theta_cdm = y[lay.cdm.idx_theta];

  dy[lay.scf.idx_phi_prime] += Type3CouplingDeltaPhiPrime(beta, phi_prime_bg, theta_cdm);
  dy[lay.cdm.idx_theta] += Type3CouplingDeltaThetaCdm(beta, k2, a_prime_over_a, rho_cdm, a2, Zbar,
                                                      dV, phi, phi_prime, phi_prime_bg, theta_cdm);
}
```
Add public accessors `int index_bg_phi_prime_scf() const;` and `int index_bg_dV_scf() const;` to `ScalarFieldSpecies` (they currently exist only as private members) so the composite can read the background slots.

- [ ] **Step 5: Add the `StressEnergy` cross-term**

In `Type3Species::StressEnergy`, after summing the children, add:
```cpp
  const double beta = scf_->beta();
  if (beta != 0.) {
    const double a = ppw->scalar_ctx.a;
    const double phi_prime_bg = pvecback[scf_->index_bg_phi_prime_scf()];
    const double Zbar = -phi_prime_bg / a;
    const double theta_cdm = y[my.cdm.idx_theta];
    se.rho_plus_p_theta += -(2. * beta / 3.) * Zbar * Zbar * theta_cdm;
  }
  return se;
```

- [ ] **Step 6: Run the coupling unit test**

Run: `cmake --build build/cmake --target test-type3-coupling --parallel && ctest --test-dir build/cmake -R test-type3-coupling --output-on-failure`
Expected: PASS — coupling vanishes at β=0, matches the hand-computed value at β=−0.5.

- [ ] **Step 7: Write Python tests that the coupling actually bites**

In `python/test_class.py`, add:
```python
    def test_type3_drives_theta_cdm(self):
        scenario = {
            'output': 'mPk',
            'gauge': 'synchronous',
            'Omega_fld': 0,
            'Omega_scf': 0.7,
            'attractor_ic_scf': 'yes',
            'scf_parameters': '1.22, 0, 0, 0',
            'scf_veta': -0.5,
            'k_output_values': '0.1',
        }
        self.scenario = dict(scenario)
        self.cosmo.set(dict(self.verbose, **scenario))
        self.cosmo.compute()
        scalar = self.cosmo.get_perturbations()['scalar'][0]
        # theta_cdm is forced to 0 in uncoupled synchronous CDM; the coupling
        # must move it away from zero at late times.
        assert np.max(np.abs(scalar['theta_cdm'])) > 0.0

    def test_type3_suppresses_growth(self):
        base = {
            'output': 'mPk', 'P_k_max_1/Mpc': 1.0, 'z_pk': 0,
            'gauge': 'synchronous', 'Omega_fld': 0, 'Omega_scf': 0.7,
            'attractor_ic_scf': 'yes', 'scf_parameters': '1.22, 0, 0, 0',
        }
        import classy
        uncoupled = classy.Class(); uncoupled.set(dict(self.verbose, **base)); uncoupled.compute()
        coupled = classy.Class(); coupled.set(dict(self.verbose, scf_veta=-0.5, **base)); coupled.compute()
        k = 0.1
        try:
            ratio = coupled.pk(k, 0) / uncoupled.pk(k, 0)
            assert ratio < 1.0  # negative beta suppresses growth (paper Fig. 1/2)
        finally:
            for c in (uncoupled, coupled):
                c.struct_cleanup(); c.empty()
```

- [ ] **Step 8: Build the wrapper and run the coupling integration tests**

Run: `make classy-pip-dev`
Run: `cd python && python -m pytest test_class.py -k "type3" -v`
Expected: PASS — θ_cdm is driven nonzero and P(k) is suppressed for β=−0.5.

- [ ] **Step 9: Commit**

```bash
git add species/type3_species.h species/type3_species.cpp species/type3_coupling_test.cpp \
        species/scalar_field.h species/test_support.h CMakeLists.txt python/test_class.py
git commit -m "type3: implement the synchronous momentum-transfer coupling + stress-energy cross-term"
```

---

### Task 6: Figure-2 validation (large negative β)

Reproduce the qualitative features of the paper's Figure 2 (`cl_ratio_beta.pdf`, `pk_ratio_beta.pdf`) over a wide range of negative β, and lock a coarse regression check. This is a validation task, not new production code.

**Files:**
- Create: `notebooks/type3_figure2_validation.py` (a script producing the ratio plots)
- Modify: `python/test_class.py` (one coarse `T3_ph` regression assertion)

**Interfaces:**
- Consumes: the complete `Type3Species` (Tasks 1-5) via the `classy` wrapper.

- [ ] **Step 1: Write the ratio-plot script**

Create `notebooks/type3_figure2_validation.py` that, for `beta in [0, -1, -1e1, -1e2, -1e3, -1e4]`, runs `classy` with fixed `omega_b, omega_cdm, 100*theta_s, A_s, n_s, tau_reio` and `Omega_scf=0.7, scf_parameters='1.22,0,0,0', gauge='synchronous'`, and plots `Cl^TT(beta)/Cl^TT(0)` and `P(k,beta)/P(k,0)`. Save figures to `notebooks/type3_cl_ratio.pdf` and `notebooks/type3_pk_ratio.pdf`.
```python
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import classy

COMMON = dict(output='tCl mPk', gauge='synchronous', Omega_fld=0, Omega_scf=0.7,
              attractor_ic_scf='yes', scf_parameters='1.22, 0, 0, 0',
              omega_b=0.0223, omega_cdm=0.119, n_s=0.967, ln10^{10}A_s=3.06,
              tau_reio=0.07, l_max_scalars=2500, **{'P_k_max_1/Mpc': 1.0, 'z_pk': 0,
              '100*theta_s': 1.0420})

def run(beta):
    c = classy.Class()
    pars = dict(COMMON)
    if beta != 0:
        pars['scf_veta'] = beta
    c.set(pars); c.compute()
    return c

ref = run(0.0)
cl0 = ref.raw_cl(2500)['tt'][2:]
ks = np.logspace(-3, 0, 200)
pk0 = np.array([ref.pk(k, 0) for k in ks])

fig_cl, ax_cl = plt.subplots()
fig_pk, ax_pk = plt.subplots()
for beta in [-1, -1e1, -1e2, -1e3, -1e4]:
    c = run(beta)
    cl = c.raw_cl(2500)['tt'][2:]
    pk = np.array([c.pk(k, 0) for k in ks])
    ax_cl.semilogx(np.arange(2, 2501), cl / cl0, label=f'beta={beta:g}')
    ax_pk.semilogx(ks, pk / pk0, label=f'beta={beta:g}')
    c.struct_cleanup(); c.empty()
ax_cl.set_xlabel('l'); ax_cl.set_ylabel('Cl_TT(b)/Cl_TT(0)'); ax_cl.legend()
ax_pk.set_xlabel('k [1/Mpc]'); ax_pk.set_ylabel('P(k,b)/P(k,0)'); ax_pk.legend()
fig_cl.savefig('notebooks/type3_cl_ratio.pdf')
fig_pk.savefig('notebooks/type3_pk_ratio.pdf')
ref.struct_cleanup(); ref.empty()
print('saved type3_cl_ratio.pdf and type3_pk_ratio.pdf')
```
*(Fix the `ln10^{10}A_s` / `100*theta_s` keys to the exact `classy` spellings used elsewhere in the repo when implementing — grep `explanatory.ini`.)*

- [ ] **Step 2: Run the script and eyeball against the paper**

Run: `make classy-pip-dev && python notebooks/type3_figure2_validation.py`
Expected: produces the two PDFs. Compare visually with `arXiv-1604.04222v2/cl_ratio_beta.pdf` and `pk_ratio_beta.pdf`: `Cl^TT` ratio departs from 1 only at low `l` (ISW); `P(k)` ratio is suppressed (<1) over a range and turns up again for `beta < -1e2`.

- [ ] **Step 3: Add a coarse regression assertion**

In `python/test_class.py` add a test asserting the monotonic-then-turnup feature is present (a cheap proxy that survives across builds without committed goldens):
```python
    def test_type3_pk_ratio_turnup_feature(self):
        import classy
        base = {'output': 'mPk', 'P_k_max_1/Mpc': 1.0, 'z_pk': 0, 'gauge': 'synchronous',
                'Omega_fld': 0, 'Omega_scf': 0.7, 'attractor_ic_scf': 'yes',
                'scf_parameters': '1.22, 0, 0, 0'}
        ref = classy.Class(); ref.set(dict(self.verbose, **base)); ref.compute()
        k = 0.2
        pk0 = ref.pk(k, 0)
        ratios = {}
        for beta in (-0.5, -1e2, -1e4):
            c = classy.Class(); c.set(dict(self.verbose, scf_veta=beta, **base)); c.compute()
            ratios[beta] = c.pk(k, 0) / pk0
            c.struct_cleanup(); c.empty()
        ref.struct_cleanup(); ref.empty()
        # Suppression at moderate |beta|, recovering/enhancing at very large |beta|.
        assert ratios[-0.5] < 1.0
        assert ratios[-1e4] > ratios[-1e2]
```

- [ ] **Step 4: Run the regression assertion**

Run: `cd python && python -m pytest test_class.py -k type3_pk_ratio_turnup -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add notebooks/type3_figure2_validation.py python/test_class.py
git commit -m "type3: Figure-2 validation script + coarse P(k) ratio regression"
```

---

## Self-Review

**Spec coverage:**
- §2 ① synchronous CDM θ → Task 2. ② StressEnergy cross-term → Task 5 Step 5. ③ background-coupling hook → explicitly out of scope (spec §9), no task — correct.
- §4.1 injectable potential → Task 1. §4.2 CDM flag → Task 2. §4.3 scalar field β → Task 3. §4.4 composite → Tasks 4-5.
- §5.1 background → Task 3 Step 4. §5.2 child free-streaming → Task 2 (CDM) + Task 3 Step 5 (scf KG). §5.3 composite coupling → Task 5 Steps 3-4. §5.4 stress-energy → Task 3 Step 5 (scf δρ/δp) + Task 5 Step 5 (cross-term). §5.5 IC/guards → Task 2 Step 4 (θ=0 IC), Task 4 Step 2 (β<1/2, gauge), Task 5 (D in the denominator).
- §6 input/factory/budget/shooting/gauge guard → Task 4. §7 testing → woven into every task + Task 6. §8 build order → Tasks 1-5 in order. §9 follow-ups → not built (correct).

**Placeholder scan:** Task 6 Step 1 flags two `classy` key spellings (`ln10^{10}A_s`, `100*theta_s`) to confirm against `explanatory.ini` at implementation time — a genuine lookup, not a logic gap. Task 4 Step 2 references `ScalarFieldSpecies::CreateAllForComposite`, defined in Task 4 Step 4. Task 5 Step 1 offers a documented fallback (free-function coupling) if a full `perturb_parameters_and_workspace` is too heavy to fixture — the fallback is the path the implementation actually takes (Steps 3-4). No bare TODOs.

**Type consistency:** `Type3Species::PerturbLayout{cdm, scf}` used consistently (Tasks 4-5). `scf_->beta()` / `cdm_->Rho()` / `index_bg_phi_prime_scf()` / `index_bg_dV_scf()` accessors introduced in Task 5 Step 4 and used there. `CDMSpecies(pba, omega0, coupled)` signature consistent (Tasks 2, 4). `ScalarFieldSpecies(..., potential, beta)` trailing-arg order consistent (Tasks 1, 3, 4).
