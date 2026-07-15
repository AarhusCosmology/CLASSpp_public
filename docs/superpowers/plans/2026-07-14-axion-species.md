# Axion Species Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add axion physics per `docs/superpowers/specs/2026-07-14-axion-species-design.md`: (1) a (1−cos)ⁿ axion potential for `ScalarFieldSpecies`, (2) AxiCLASS's `pheno_axion` effective EDE fluid as a `FluidSpecies` subclass.

**Architecture:** Both components are species-plugin work only — no module changes. Component 1 injects a new `ScalarFieldPotential` bundle plus a factory branch in `ScalarFieldSpecies::CreateAll`. Component 2 adds a `PhenoAxion` equation-of-state enum value, virtualizes `FluidSpecies::ComputeWFld` and the sound speed (`Cs2(k2, a)` hook), and implements `AxionEDEFluid : FluidSpecies` built by `FluidSpecies::CreateAll`.

**Tech Stack:** C++17 (repo style), CMake, plain assert-based test executables, `./class <ini>` integration runs, python3+numpy for output verification.

## Global Constraints

- Work on branch `axion-species` (created in Task 1). **Subagents must never `git checkout`, `git reset`, or switch branches** — all work happens on the current checkout of this branch.
- **Never `git add -A` or `git add .`** — stage explicit paths only. Test .ini fixtures at repo root are gitignored: stage with `git add -f`.
- CLASSpp is C++-only: no `extern "C"`, no `#ifdef __cplusplus` in headers.
- Never hand-edit `cclassy.pxd` (auto-generated). No `struct background` members are added by this plan, so the wrapper is untouched.
- Numerical comparisons: never blind max-rel-diff; use scale-aware tolerances (~0.1% typical); avoid Cl^TE near zero-crossings (use TT/EE only).
- Build: `cmake --build build -j8` from repo root (if `build/` is stale/missing: `cmake -B build -DCMAKE_BUILD_TYPE=Release` first). Unit tests: `ctest --test-dir build -R <name> --output-on-failure`. Full suite: `ctest --test-dir build --output-on-failure`. Integration binary: `./class <file>.ini` from repo root.
- Reference formulas: AxiCLASS clone at
  `/private/tmp/claude-497981304/-Users-au192734-Projects-class-claude/235a779f-51ca-4b58-a9b6-fec0301d4959/scratchpad/AxiCLASS`
  (referenced below as `$AXICLASS`). Port physics, not bugs — the spec lists two known AxiCLASS bugs (comma-operator dw/da; n=1 ddV instability) that must NOT be ported.
- Commit after every task with the trailer:
  ```
  Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01EWPTcYe8AfdxcSqmx4b9ye
  ```

---

### Task 1: AxionScalarFieldPotential bundle

**Files:**
- Modify: `species/scalar_field_potential.h` (declare factory)
- Modify: `species/scalar_field_potential.cpp` (implement)
- Test: `species/scalar_field_potential_test.cpp` (extend existing main)

**Interfaces:**
- Produces: `ScalarFieldPotential AxionScalarFieldPotential();` with params layout `[m, f, n, Theta_ini]` (m in 1/Mpc, f and φ in reduced-Planck units, n ≥ 1, Θ_ini ∈ (0,π)). Task 2 consumes this factory and the params layout.

- [ ] **Step 1: Create the branch and verify baseline**

```bash
cd /Users/au192734/Projects/class_claude
git checkout -b axion-species
cmake --build build -j8 && ctest --test-dir build --output-on-failure
```
Expected: build succeeds, all existing tests pass. (This is the ONLY permitted checkout — it creates the feature branch at the start of the whole plan. If branch `axion-species` already exists, you are mid-plan: do NOT checkout, just verify `git branch --show-current` prints `axion-species`.)

- [ ] **Step 2: Write the failing test**

Append to the existing `main()` in `species/scalar_field_potential_test.cpp`, before the final `printf`/`return`:

```cpp
  // ── Axion bundle: V = m^2 f^2 (1 - cos(phi/f))^n, params = [m, f, n, Theta_ini]
  {
    const ScalarFieldPotential ax = AxionScalarFieldPotential();
    for (double n : {1.0, 2.0, 2.5, 3.0}) {
      const std::vector<double> ap = {2.0e-3, 0.5, n, 2.0};  // m, f, n, Theta
      const double m = ap[0], f = ap[1];
      // V matches the closed form away from the minimum.
      for (double ph : {0.05, 0.4, 1.0, 1.4}) {
        const double u = 1.0 - std::cos(ph / f);
        const double V_exp = m * m * f * f * std::pow(u, n);
        assert(std::fabs(ax.V(ph, ap) - V_exp) < 1e-12 * std::fabs(V_exp) + 1e-300);
        // dV, ddV against central finite differences.
        const double h = 1e-6;
        const double dV_fd  = (ax.V(ph + h, ap) - ax.V(ph - h, ap)) / (2 * h);
        const double ddV_fd = (ax.V(ph + h, ap) - 2 * ax.V(ph, ap) + ax.V(ph - h, ap)) / (h * h);
        assert(std::fabs(ax.dV(ph, ap) - dV_fd) < 1e-5 * std::fabs(dV_fd) + 1e-12);
        assert(std::fabs(ax.ddV(ph, ap) - ddV_fd) < 1e-3 * std::fabs(ddV_fd) + 1e-8);
      }
      // Regularity at the minimum: no NaN/inf, and for n=1 ddV(0) = m^2 exactly.
      assert(std::isfinite(ax.V(0.0, ap)));
      assert(std::isfinite(ax.dV(0.0, ap)));
      assert(std::isfinite(ax.ddV(0.0, ap)));
    }
    const std::vector<double> a1 = {2.0e-3, 0.5, 1.0, 2.0};
    assert(std::fabs(ax.ddV(0.0, a1) - a1[0] * a1[0]) < 1e-15);

    // Frozen-field shooting guess: V(Theta*f) = 3 H0^2 Omega at the guessed m.
    const double H0 = 2.2e-4, omega = 0.05;
    const auto [mg, dmdo] = ax.shooting_guess(omega, H0, a1, 0);
    std::vector<double> ag = a1;
    ag[0] = mg;
    assert(std::fabs(ax.V(a1[3] * a1[1], ag) - 3 * H0 * H0 * omega) <
           1e-10 * 3 * H0 * H0 * omega);
    assert(std::fabs(dmdo - mg / (2 * omega)) < 1e-12 * mg / omega);
  }
```

- [ ] **Step 3: Run test to verify it fails**

```bash
cmake --build build -j8 --target test-scf-potential 2>&1 | tail -5
```
Expected: compile FAILURE — `AxionScalarFieldPotential` not declared.

- [ ] **Step 4: Implement**

In `species/scalar_field_potential.h`, after `DefaultScalarFieldPotential()`:

```cpp
/** Axion bundle: V = m^2 f^2 (1 - cos(phi/f))^n, params = [m, f, n, Theta_ini]
 *  (m in 1/Mpc; f, phi in reduced-Planck units; n >= 1; Theta_ini in (0, pi)).
 *  ddV is written with sin^2 = u(2-u) so it stays regular at the minimum u -> 0.
 *  The shooting guess assumes a frozen field today: V(Theta_ini*f) = 3 H0^2 Omega. */
ScalarFieldPotential AxionScalarFieldPotential();
```

In `species/scalar_field_potential.cpp`, append:

```cpp
ScalarFieldPotential AxionScalarFieldPotential() {
  ScalarFieldPotential p;
  p.V = [](double phi, const std::vector<double>& q) {
    const double m = q[0], f = q[1], n = q[2];
    const double u = 1. - std::cos(phi / f);
    return m * m * f * f * std::pow(u, n);
  };
  p.dV = [](double phi, const std::vector<double>& q) {
    const double m = q[0], f = q[1], n = q[2];
    const double u = 1. - std::cos(phi / f);
    return m * m * f * n * std::pow(u, n - 1.) * std::sin(phi / f);
  };
  p.ddV = [](double phi, const std::vector<double>& q) {
    const double m = q[0], f = q[1], n = q[2];
    const double u = 1. - std::cos(phi / f);
    // ddV = m^2 n u^(n-1) [ (n-1)(2-u) + (1-u) ]: sin^2 = u(2-u) removes the
    // u^(n-2) singularity, so this is regular at u -> 0 for all n >= 1.
    return m * m * n * std::pow(u, n - 1.) * ((n - 1.) * (2. - u) + (1. - u));
  };
  p.shooting_guess = [](double omega,
                        double H0,
                        const std::vector<double>& q,
                        int /*tuning_index*/) -> std::pair<double, double> {
    const double f = q[1], n = q[2], theta = q[3];
    const double u0 = 1. - std::cos(theta);
    const double m  = std::sqrt(3. * H0 * H0 * omega / (f * f * std::pow(u0, n)));
    return {m, m / (2. * omega)};  // dm/dOmega for m = sqrt(c*Omega)
  };
  return p;
}
```

- [ ] **Step 5: Run test to verify it passes**

```bash
cmake --build build -j8 --target test-scf-potential && ctest --test-dir build -R test-scf-potential --output-on-failure
```
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add species/scalar_field_potential.h species/scalar_field_potential.cpp species/scalar_field_potential_test.cpp
git commit -m "Add axion (1-cos)^n scalar-field potential bundle"
```
(Include the Global Constraints trailer.)

---

### Task 2: Axion branch in ScalarFieldSpecies::CreateAll

**Files:**
- Modify: `species/scalar_field.h` (declare `CreateAxion`)
- Modify: `species/scalar_field.cpp` (branch in `CreateAll` ~line 518, new `CreateAxion`)
- Create: `species/axion_scf_factory_test.cpp`
- Modify: `CMakeLists.txt` (~line 215 test block: add `test-axion-scf-factory`)

**Interfaces:**
- Consumes: `AxionScalarFieldPotential()` from Task 1; existing `ScalarFieldSpecies` constructor and shooting members (`shooting_target_`, `needs_shooting_`, `scf_parameters_`, `ComputeShootingGuess`).
- Produces: input vocabulary `scf_potential = axion`, `m_axion`, `f_axion`, `n_axion`, `Theta_initial_scf` wired end-to-end (Task 3 runs it via .ini).

Semantics being implemented (from the spec, with one amendment): the axion scalar field always resolves m by **shooting** against Omega_scf (explicit `Omega_scf > 0` or budget-closure override); `m_axion` is an optional Newton seed that replaces the frozen-field guess. There is no "specify m, derive Omega" mode (the budget architecture has no place for a floating Omega — closure mode IS the floating-Omega mode, and it shoots).

- [ ] **Step 1: Write the failing test**

Create `species/axion_scf_factory_test.cpp`:

```cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

#include "background.h"
#include "parser.h"
#include "scalar_field.h"
#include "species_build_context.h"

static SpeciesBuildContext MakeCtx(FileContent& fc, const background& pba) {
  SpeciesBuildContext ctx{};
  ctx.pfc = &fc;
  ctx.pba = &pba;
  return ctx;
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  // ── Valid axion config: shooting on m, frozen ICs, params = [m, f, n, Theta].
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("n_axion", "2");
    fc.set("Theta_initial_scf", "2.0");
    auto ctx    = MakeCtx(fc, pba);
    auto result = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(!scf.attractor_ic_scf());
    assert(std::fabs(scf.phi_ini_scf() - 2.0 * 0.5) < 1e-15);
    assert(scf.phi_prime_ini_scf() == 0.);
    assert(scf.scf_tuning_index() == 0);
    assert(scf.scf_parameters().size() == 4);
    assert(scf.scf_parameters()[1] == 0.5);
    assert(scf.scf_parameters()[2] == 2.0);
    assert(scf.scf_parameters()[3] == 2.0);
    assert(scf.scf_parameters()[0] > 0.);  // frozen-field guess seeded
    assert(scf.GetShootingTargets().size() == 1);
    assert(scf.GetShootingTargets()[0].target_name == "Omega_scf");
    assert(scf.GetShootingTargets()[0].unknown_param == "scf_shooting_parameter");
  }

  // ── m_axion, when given, is used as the Newton seed.
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("m_axion", "3.3e-3");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    auto ctx    = MakeCtx(fc, pba);
    auto result = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(std::fabs(scf.scf_parameters()[0] - 3.3e-3) < 1e-18);
    assert(scf.scf_parameters()[2] == 1.0);  // n_axion defaults to 1
    assert(scf.GetShootingTargets().size() == 1);
  }

  // ── scf_shooting_parameter (resolved value from DoShooting) disables shooting.
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    fc.set("scf_shooting_parameter", "4.4e-3");
    auto ctx    = MakeCtx(fc, pba);
    auto result = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(std::fabs(scf.scf_parameters()[0] - 4.4e-3) < 1e-18);
    assert(scf.GetShootingTargets().empty());
  }

  // ── Closure mode: override supplies Omega, shooting still armed.
  {
    FileContent fc;
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    auto ctx                     = MakeCtx(fc, pba);
    ctx.omega0_closure_override  = 0.68;
    auto result                  = ScalarFieldSpecies::CreateAll(ctx);
    assert(result.size() == 1);
    auto& scf = static_cast<ScalarFieldSpecies&>(*result[0].species);
    assert(scf.GetShootingTargets().size() == 1);
    assert(std::fabs(scf.GetShootingTargets()[0].target_value - 0.68) < 1e-15);
  }

  // ── Rejections: attractor requested; scf_parameters given; missing keys; bad ranges.
  auto expect_throw = [&pba](void (*mut)(FileContent&)) {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    fc.set("f_axion", "0.5");
    fc.set("Theta_initial_scf", "2.0");
    mut(fc);
    auto ctx = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  };
  expect_throw([](FileContent& fc) { fc.set("attractor_ic_scf", "yes"); });
  expect_throw([](FileContent& fc) { fc.set("scf_parameters", "1.0, 2.0"); });
  expect_throw([](FileContent& fc) { fc.set("f_axion", "-0.5"); });
  expect_throw([](FileContent& fc) { fc.set("n_axion", "0.5"); });
  expect_throw([](FileContent& fc) { fc.set("Theta_initial_scf", "3.5"); });

  // Missing f_axion / Theta_initial_scf.
  for (const char* missing : {"f_axion", "Theta_initial_scf"}) {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "axion");
    if (std::string(missing) != "f_axion")
      fc.set("f_axion", "0.5");
    if (std::string(missing) != "Theta_initial_scf")
      fc.set("Theta_initial_scf", "2.0");
    auto ctx   = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  // ── Unknown scf_potential value still errors.
  {
    FileContent fc;
    fc.set("Omega_scf", "0.05");
    fc.set("scf_potential", "banana");
    auto ctx   = MakeCtx(fc, pba);
    bool threw = false;
    try {
      ScalarFieldSpecies::CreateAll(ctx);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  }

  std::printf("axion scf factory tests passed\n");
  return 0;
}
```

Register in `CMakeLists.txt`: add `add_executable(test-axion-scf-factory species/axion_scf_factory_test.cpp)` next to the other species tests (~line 217) and add `test-axion-scf-factory` to the `foreach(_t IN ITEMS ...)` list.

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build -j8 --target test-axion-scf-factory 2>&1 | tail -5
```
Expected: compile or assert FAILURE (the axion branch doesn't exist; `CreateAll` returns a default-potential species or throws nothing).

- [ ] **Step 3: Implement**

In `species/scalar_field.h`, add to the public section of `ScalarFieldSpecies` (near `CreateAll`):

```cpp
  // Axion factory branch: scf_potential = axion. params = [m, f, n, Theta_ini],
  // frozen ICs (phi = Theta*f, phi' = 0), always shoots m (tuning index 0)
  // against Omega_scf / the closure value; m_axion is an optional Newton seed.
  static std::vector<Named> CreateAxion(const SpeciesBuildContext& ctx, double omega0_scf);
```

In `species/scalar_field.cpp`, inside `CreateAll` right after the `if (omega0_scf == 0.) return result;` check (line ~518), insert:

```cpp
  // ── scf_potential dispatch ────────────────────────────────────────────────
  if (auto pot_opt = ctx.pfc->get<std::string>("scf_potential")) {
    if (*pot_opt == "axion")
      return CreateAxion(ctx, omega0_scf);
    throw std::invalid_argument("unknown scf_potential '" + *pot_opt +
                                "': supported values are 'axion' (absent = historical default)");
  }
```

Then add the new factory (before `CreateAllForComposite`):

```cpp
std::vector<Named> ScalarFieldSpecies::CreateAxion(const SpeciesBuildContext& ctx,
                                                   double omega0_scf) {
  std::vector<Named> result;

  const auto f_opt     = ctx.pfc->get<double>("f_axion");
  const auto theta_opt = ctx.pfc->get<double>("Theta_initial_scf");
  const double n       = ctx.pfc->get_or("n_axion", 1.);
  if (!f_opt || *f_opt <= 0.)
    throw std::invalid_argument(
        "scf_potential = axion requires f_axion > 0 (decay constant, reduced-Planck units)");
  if (!theta_opt || *theta_opt <= 0. || *theta_opt >= M_PI)
    throw std::invalid_argument(
        "scf_potential = axion requires Theta_initial_scf in (0, pi): phi_ini = Theta*f");
  if (n < 1.)
    throw std::invalid_argument("n_axion must be >= 1 (V = m^2 f^2 (1-cos(phi/f))^n)");
  const double f     = *f_opt;
  const double theta = *theta_opt;

  // A (1-cos)^n potential has no exponential attractor: frozen ICs are forced.
  if (auto attr = ctx.pfc->get<std::string>("attractor_ic_scf")) {
    if (attr->find("y") != std::string::npos || attr->find("Y") != std::string::npos)
      throw std::invalid_argument(
          "attractor_ic_scf = yes is incompatible with scf_potential = axion "
          "(frozen ICs phi = Theta_initial_scf * f_axion, phi' = 0 are set automatically)");
  }
  if (ctx.pfc->get<std::vector<double>>("scf_parameters").has_value())
    throw std::invalid_argument(
        "scf_parameters is not used with scf_potential = axion; "
        "give m_axion, f_axion, n_axion, Theta_initial_scf instead");
  ctx.pfc->get<int>("scf_tuning_index");  // consume; the axion branch forces 0

  std::vector<double> params = {ctx.pfc->get_or("m_axion", 0.), f, n, theta};

  // Resolved shooting value (written back by DoShooting) wins over any seed.
  bool shooting_param_present = false;
  if (auto sp = ctx.pfc->get<double>("scf_shooting_parameter")) {
    params[0]              = *sp;
    shooting_param_present = true;
  }

  auto species = std::make_unique<ScalarFieldSpecies>(*ctx.pba,
                                                      omega0_scf,
                                                      std::move(params),
                                                      /*scf_tuning_index=*/0,
                                                      /*attractor_ic_scf=*/false,
                                                      /*phi_ini_scf=*/theta * f,
                                                      /*phi_prime_ini_scf=*/0.,
                                                      AxionScalarFieldPotential());
  if (!shooting_param_present) {
    species->shooting_target_ = {"Omega_scf", "scf_shooting_parameter", omega0_scf};
    species->needs_shooting_  = true;
    if (species->scf_parameters_[0] <= 0.) {
      std::vector<double> g, d;
      species->ComputeShootingGuess(ctx, g, d);
      species->scf_parameters_[0] = g[0];
    }
  }
  result.push_back({"ScalarField", std::move(species)});
  return result;
}
```

Notes for the implementer:
- `#include "scalar_field_potential.h"` is already included via `scalar_field.h`. `M_PI` needs `<cmath>` (already included).
- The closure path: `CreateAll` computes `omega0_scf` from `ctx.omega0_closure_override` before the dispatch, so `CreateAxion` receives the closure value and arms shooting against it — this deliberately differs from the default potential's closure path (which self-normalizes via attractor ICs). `target_value = omega0_scf` covers both the explicit-Omega and closure cases.
- `FileContent::get<int>` is a supported instantiation (include/parser.h line ~58 lists int); the bare call consumes the key so the unread-parameter check stays quiet.

- [ ] **Step 4: Run tests**

```bash
cmake --build build -j8 && ctest --test-dir build -R "test-axion-scf-factory|test-scf" --output-on-failure
```
Expected: PASS (both new and existing scf tests).

- [ ] **Step 5: Full regression**

```bash
ctest --test-dir build --output-on-failure
```
Expected: all tests pass (no existing behavior touched — the new branch only triggers on `scf_potential`, a previously-absent key; run `./class explanatory.ini` as an extra smoke check if in doubt).

- [ ] **Step 6: Commit**

```bash
git add species/scalar_field.h species/scalar_field.cpp species/axion_scf_factory_test.cpp CMakeLists.txt
git commit -m "Wire scf_potential=axion into ScalarFieldSpecies::CreateAll"
```
(Include the Global Constraints trailer.)

---

### Task 3: scf-axion integration run (.ini fixture)

**Files:**
- Create: `test_axion_scf.ini` (repo root, gitignored → `git add -f`)

**Interfaces:**
- Consumes: Task 2's input vocabulary end-to-end (input → shooting → background → Cls).

- [ ] **Step 1: Write the fixture**

Create `test_axion_scf.ini`:

```ini
# Axion scalar field (exact Klein-Gordon), dark-energy regime: the field is
# frozen deep in the past and thaws around a ~ 0.5. Shooting tunes m_axion.
h = 0.67
omega_b = 0.02238
omega_cdm = 0.1201
output = tCl
Omega_scf = 0.05
scf_potential = axion
f_axion = 0.5
n_axion = 1
Theta_initial_scf = 2.0
write background = yes
root = output/test_axion_scf_
background_verbose = 1
input_verbose = 1
```

- [ ] **Step 2: Run it**

```bash
./class test_axion_scf.ini
```
Expected: exit 0; shooting converges (verbose output mentions the scf shooting); no NaN warnings.

- [ ] **Step 3: Verify Omega_scf was hit and the field history is sane**

First inspect the header (`head -8 output/test_axion_scf_background.dat`). CLASS background files carry a final comment line of `N:column_name` tokens. Then run:

```bash
python3 - <<'EOF'
import re
import numpy as np

path = 'output/test_axion_scf_background.dat'
names = None
with open(path) as fh:
    for line in fh:
        if line.startswith('#'):
            toks = re.findall(r'(\d+):(\S+)', line)
            if toks:
                names = [t[1] for t in toks]
        else:
            break
data = np.loadtxt(path)
col = {n: i for i, n in enumerate(names)}
# Column names to adjust after inspecting the header if they differ:
z, rho_scf, H = (data[:, col['z']], data[:, col['(.)rho_scf']], data[:, col['H']])

H0 = H[z.argmin()]                       # today = smallest z row
omega_scf = rho_scf[z.argmin()] / H0**2
assert abs(omega_scf - 0.05) < 1e-3, omega_scf

early = rho_scf[z > 10]
assert early.max() / early.min() - 1 < 0.01, "field not frozen for z > 10"
assert (rho_scf > 0).all()
print("scf-axion background OK: Omega_scf =", omega_scf)
EOF
```
Expected: prints `scf-axion background OK` with Omega_scf ≈ 0.05. Adjust the three column-name strings to the actual header tokens (e.g. `H` may print as `H[1/Mpc]`) — do not change the assertions. If the frozen-field check fails because the chosen (f, Θ, n) thaws earlier than expected, loosen only the redshift threshold (this is a smoke test, not physics validation) and note it in the commit message.

- [ ] **Step 4: Commit**

```bash
git add -f test_axion_scf.ini
git commit -m "Add scf-axion integration fixture (dark-energy regime smoke test)"
```
(Include the Global Constraints trailer.)

---

### Task 4: PhenoAxion enum + Cs2/ComputeWFld virtualization (mechanical, no behavior change)

**Files:**
- Modify: `source/background.h:20` (enum)
- Modify: `species/fluid.h` (virtualize `ComputeWFld`, add `Cs2` hook)
- Modify: `species/fluid.cpp:138,200-204,284-285` (route through `Cs2`)

**Interfaces:**
- Produces: `enum equation_of_state { CLP, EDE, PhenoAxion };`
  `virtual void ComputeWFld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) const;`
  `virtual double Cs2(double k2, double a) const;` (base returns `cs2_fld_`). Task 5 overrides both.

- [ ] **Step 1: Make the changes**

`source/background.h:20`:
```cpp
enum equation_of_state { CLP, EDE, PhenoAxion };
```

`species/fluid.h`: change the `ComputeWFld` declaration (line ~92) to `virtual`, and add below it:
```cpp
  virtual void ComputeWFld(double a, double* w_fld, double* dw_over_da_fld, double* integral_fld) const;

  /** Effective sound speed in the fluid rest frame, cs2(k^2, a). The base fluid
   *  uses the constant input cs2_fld; AxionEDEFluid overrides with the
   *  k- and a-dependent GDM formula. Takes k^2 to avoid a sqrt in PerturbDerivs. */
  virtual double Cs2(double /*k2*/, double /*a*/) const {
    return cs2_fld_;
  }
```

`species/fluid.cpp`, three replacements:
- Line 138: `const double cs2 = cs2_fld_;` → `const double cs2 = Cs2(k2, a);`
- Lines 200–204 in `ApplyInitialConditions`: insert `const double cs2 = Cs2(ctx.k * ctx.k, ctx.a);` after the `background_w_fld` call and replace all four `cs2_fld_` occurrences with `cs2`.
- Lines 284–285 in `StressEnergy`: insert `const double cs2 = Cs2(k2, a);` after `const double a_prime_over_a = ...` (both `k2` and `a` are already in scope) and replace both `cs2_fld_` occurrences with `cs2`.

Leave `ppf_fluid.cpp:114` alone (PPF's use of `cs2_fld_` is a transition-scale parameter, and PPF is never pheno_axion).

- [ ] **Step 2: Build and run the full suite (behavior-change check)**

```bash
cmake --build build -j8 && ctest --test-dir build --output-on-failure
```
Expected: everything passes — this is a pure mechanical refactor; the base `Cs2` returns exactly the old constant.

- [ ] **Step 3: Integration spot-check**

```bash
./class base_2018_plikHM_TTTEEE_lowl_lowE_lensing.ini 2>&1 | tail -3
```
Expected: runs to completion (this .ini exercises the CLP/PPF fluid path if configured; if it has no fluid, run `./class explanatory.ini` instead — the point is a smoke check that nothing crashed).

- [ ] **Step 4: Commit**

```bash
git add source/background.h species/fluid.h species/fluid.cpp
git commit -m "Add PhenoAxion eos value; virtualize fluid ComputeWFld and sound speed (no behavior change)"
```
(Include the Global Constraints trailer.)

---

### Task 5: AxionEDEFluid class + unit tests

**Files:**
- Create: `species/axion_ede_fluid.h`
- Create: `species/axion_ede_fluid.cpp`
- Create: `species/axion_ede_fluid_test.cpp`
- Modify: `CMakeLists.txt` (add `species/axion_ede_fluid.cpp` to the classpp source list ~line 83 next to `species/fluid.cpp`; add `test-axion-ede-fluid` executable + foreach entry ~line 215)

**Interfaces:**
- Consumes: Task 4's virtual hooks; `FluidSpecies` constructor (`fluid.h:35`); `BaseSpecies::energy_type()` / `GetOmega0()`; `PerturbIcContext` fields (`a`, `k`, `ktau_two`, `ktau_three`, `ppr->curvature_ini`, `s2_squared`, `index_ic`, `p_mod`).
- Produces (Task 6 consumes the constructor; Task 7 exercises the physics):
```cpp
AxionEDEFluid(const background& pba, double omega0_fld, double a_c, double n_axion,
              double nu, double w_i, double w_f, double theta_i);
static double WFinal(double n);
static double Integral3OnePlusWOverA(double a, double a_c, double nu, double w_i, double w_f);
static double OmegaZeroFromOmegaAc(double omega_ac, double a_c, double nu, double w_i, double w_f);
static double OmegaAcFromOmegaZero(double omega0, double a_c, double nu, double w_i, double w_f);
```

- [ ] **Step 1: Write the failing test**

Create `species/axion_ede_fluid_test.cpp`:

```cpp
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

#include "axion_ede_fluid.h"
#include "background.h"

// Log-spaced trapezoid quadrature of 3(1+w)/a on [a, 1] for the reference integral.
static double QuadIntegral(const AxionEDEFluid& fld, double a) {
  const int N  = 200000;  // log-spaced, plenty for 1e-14 -> 1
  const double la = std::log(a);
  const double h  = -la / N;
  auto f = [&](double lna) {
    double w, dw, integ;
    fld.ComputeWFld(std::exp(lna), &w, &dw, &integ);
    return 3. * (1. + w);
  };
  double s = 0.5 * (f(la) + f(0.));
  for (int i = 1; i < N; ++i)
    s += f(la + i * h);
  return s * h;
}

int main() {
  background pba{};
  pba.H0 = 2.2e-4;

  const double a_c = std::pow(10., -3.5), nu = 1., w_i = -1., n = 3.;
  const double w_f = AxionEDEFluid::WFinal(n);
  assert(std::fabs(w_f - 0.5) < 1e-15);

  AxionEDEFluid fld(pba, /*omega0=*/1e-6, a_c, n, nu, w_i, w_f, /*theta_i=*/2.8);

  // ── w(a): sigmoid limits and midpoint.
  double w, dw, integ;
  fld.ComputeWFld(1e-8, &w, &dw, &integ);
  assert(std::fabs(w - w_i) < 1e-10);
  fld.ComputeWFld(1.0, &w, &dw, &integ);
  assert(std::fabs(w - w_f) < 1e-3);  // a >> a_c
  fld.ComputeWFld(a_c, &w, &dw, &integ);
  assert(std::fabs(w - 0.5 * (w_i + w_f)) < 1e-12);  // x = 1 at a = a_c

  // ── a = 0 exactly must not NaN (BackgroundModule::background_checks calls it).
  fld.ComputeWFld(0., &w, &dw, &integ);
  assert(w == w_i);
  assert(dw == 0.);
  assert(std::isfinite(integ) || integ == 0.);

  // ── dw/da vs central finite difference over the transition.
  for (double a : {a_c / 10., a_c / 2., a_c, 2. * a_c, 10. * a_c}) {
    const double h = a * 1e-6;
    double wp, wm, d_;
    fld.ComputeWFld(a + h, &wp, &d_, &integ);
    fld.ComputeWFld(a - h, &wm, &d_, &integ);
    const double fd = (wp - wm) / (2. * h);
    fld.ComputeWFld(a, &w, &dw, &integ);
    assert(std::fabs(dw - fd) < 1e-4 * std::fabs(fd) + 1e-12);
  }

  // ── Closed-form integral vs quadrature at several epochs.
  for (double a : {1e-10, 1e-5, a_c, 1e-2, 0.5}) {
    fld.ComputeWFld(a, &w, &dw, &integ);
    const double q = QuadIntegral(fld, a);
    assert(std::fabs(integ - q) < 1e-4 * std::fabs(q) + 1e-8);
  }

  // ── Omega0 <-> Omega_ac round trip.
  const double om_ac = 0.05;
  const double om0   = AxionEDEFluid::OmegaZeroFromOmegaAc(om_ac, a_c, nu, w_i, w_f);
  assert(std::fabs(AxionEDEFluid::OmegaAcFromOmegaZero(om0, a_c, nu, w_i, w_f) - om_ac) <
         1e-12 * om_ac);
  // And it matches direct density evolution: rho(a_c)/rho(1) = exp(integral(a_c)).
  fld.ComputeWFld(a_c, &w, &dw, &integ);
  assert(std::fabs(om0 * std::exp(integ) - om_ac) < 1e-10 * om_ac);

  // ── cs2 limits: k -> infinity gives 1; k -> 0 gives w_f = (n-1)/(n+1);
  //    monotonic in k^2. (omega_axion_ must be set: use the test hook.)
  fld.SetOmegaAxionForTest(1e2 * pba.H0);
  const double a_test = 1e-2;
  assert(std::fabs(fld.Cs2(1e30, a_test) - 1.) < 1e-6);
  assert(std::fabs(fld.Cs2(1e-30, a_test) - (n - 1.) / (n + 1.)) < 1e-6);
  double prev = 0.;
  for (double k2 : {1e-8, 1e-4, 1., 1e4, 1e8}) {
    const double c = fld.Cs2(k2, a_test);
    assert(c >= prev - 1e-15 && c <= 1. + 1e-15);
    prev = c;
  }

  std::printf("axion EDE fluid tests passed\n");
  return 0;
}
```

Register in `CMakeLists.txt` (library source list + test executable + foreach), as in Task 2.

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build -j8 --target test-axion-ede-fluid 2>&1 | tail -5
```
Expected: compile FAILURE — `axion_ede_fluid.h` doesn't exist.

- [ ] **Step 3: Implement**

Create `species/axion_ede_fluid.h`:

```cpp
#pragma once

#include <string_view>

#include "fluid.h"

/**
 * Pheno-axion early-dark-energy fluid (AxiCLASS's ede_parametrization = pheno_axion;
 * Poulin et al. 1806.10608, 1811.04083, 1905.12618). Background: sigmoid equation of
 * state w(a) = w_i + (w_f - w_i)/(1 + (a_c/a)^r), r = 3(w_f - w_i)/nu, with a
 * closed-form density integral. Perturbations: true-fluid delta/theta with the
 * k- and a-dependent GDM effective sound speed
 *   cs2(k,a) = [2 a^2 (n-1) wbar^2 + k^2] / [2 a^2 (n+1) wbar^2 + k^2],
 *   wbar(a) = omega_axion * a^(-3(n-1)/(n+1)),
 * where omega_axion is derived from (a_c, Omega_fld_ac, Theta_i, n) at
 * SetBackgroundModule time (needs the full species budget for E(a_c)).
 * PPF does not apply (w > -1 for all a > 0, and the GDM cs2 must enter the true
 * fluid equations), so this class is always built as a non-PPF fluid.
 */
class AxionEDEFluid : public FluidSpecies {
 public:
  AxionEDEFluid(const background& pba,
                double omega0_fld,
                double a_c,
                double n_axion,
                double nu,
                double w_i,
                double w_f,
                double theta_i);

  // ── closed-form pieces (unit-testable without a background pipeline) ──────
  /** w_f = (n-1)/(n+1). */
  static double WFinal(double n);
  /** int_a^1 da' 3(1+w(a'))/a' for the sigmoid w(a). Overflow-safe (log1p form). */
  static double Integral3OnePlusWOverA(double a, double a_c, double nu, double w_i, double w_f);
  /** Omega0_fld from Omega_fld_ac = rho_fld(a_c)/rho_crit0 (and inverse). */
  static double OmegaZeroFromOmegaAc(double omega_ac, double a_c, double nu, double w_i, double w_f);
  static double OmegaAcFromOmegaZero(double omega0, double a_c, double nu, double w_i, double w_f);

  void ComputeWFld(double a,
                   double* w_fld,
                   double* dw_over_da_fld,
                   double* integral_fld) const override;
  double Cs2(double k2, double a) const override;
  void ApplyInitialConditions(const BaseSpecies::PerturbLayout& layout,
                              double* y,
                              const PerturbIcContext& ctx) override;
  void SetBackgroundModule(const BackgroundModule* bgm) override;

  double a_c() const {
    return a_c_;
  }
  double n_axion() const {
    return n_axion_;
  }
  double omega_axion() const {
    return omega_axion_;
  }
  /** Unit-test hook: Cs2 needs omega_axion_, normally derived in SetBackgroundModule. */
  void SetOmegaAxionForTest(double omega_axion) {
    omega_axion_ = omega_axion;
  }

 private:
  /** Port of AxiCLASS background.c 982-1023 (Eqs. 27/28/30 of 1806.10608):
   *  derive m_fld, alpha_fld, omega_axion from (a_c, Omega_fld_ac, Theta_i, n)
   *  and the species budget at a_c. Called from SetBackgroundModule. */
  void DeriveAxionScales();

  double a_c_         = 0.;
  double n_axion_     = 3.;
  double nu_          = 1.;
  double w_i_         = -1.;
  double w_f_         = 0.5;
  double theta_i_     = 0.;
  const background& pba_ref_;  // FluidSpecies::pba_ is private; keep our own ref for H0
  double omega_ac_    = 0.;  // rho_fld(a_c)/rho_crit0
  double m_fld_       = 0.;  // axion mass in units of H0 (derived, diagnostic)
  double alpha_fld_   = 0.;  // decay constant in reduced-Planck units (derived, diagnostic)
  double omega_axion_ = 0.;  // characteristic oscillation frequency today [1/Mpc]
};
```

Create `species/axion_ede_fluid.cpp`:

```cpp
#include "axion_ede_fluid.h"

#include <cmath>
#include <limits>

#include "background.h"
#include "background_module.h"
#include "perturbations.h"
#include "perturbations_module.h"

namespace {
// ln(1 + e^t), overflow-safe for large positive t.
double LnOnePlusExp(double t) {
  return (t > 0.) ? t + std::log1p(std::exp(-t)) : std::log1p(std::exp(t));
}
}  // namespace

AxionEDEFluid::AxionEDEFluid(const background& pba,
                             double omega0_fld,
                             double a_c,
                             double n_axion,
                             double nu,
                             double w_i,
                             double w_f,
                             double theta_i)
    : FluidSpecies(pba,
                   omega0_fld,
                   PhenoAxion,
                   /*w0_fld=*/w_i,
                   /*wa_fld=*/0.,
                   /*cs2_fld=*/1.,
                   /*Omega_EDE=*/0.),
      a_c_(a_c), n_axion_(n_axion), nu_(nu), w_i_(w_i), w_f_(w_f), theta_i_(theta_i),
      pba_ref_(pba), omega_ac_(OmegaAcFromOmegaZero(omega0_fld, a_c, nu, w_i, w_f)) {}

double AxionEDEFluid::WFinal(double n) {
  return (n - 1.) / (n + 1.);
}

double AxionEDEFluid::Integral3OnePlusWOverA(double a,
                                             double a_c,
                                             double nu,
                                             double w_i,
                                             double w_f) {
  // int_a^1 3(1+w)/a' da' = 3(1+w_i) ln(1/a) + nu [ ln(1+a_c^-r) - ln(1+(a/a_c)^r) ],
  // r = 3(w_f - w_i)/nu. (Derived fresh; agrees with AxiCLASS input.c:4154 for w_i = -1.)
  const double r = 3. * (w_f - w_i) / nu;
  return 3. * (1. + w_i) * std::log(1. / a) +
         nu * (LnOnePlusExp(-r * std::log(a_c)) - LnOnePlusExp(r * std::log(a / a_c)));
}

double AxionEDEFluid::OmegaZeroFromOmegaAc(double omega_ac,
                                           double a_c,
                                           double nu,
                                           double w_i,
                                           double w_f) {
  return omega_ac / std::exp(Integral3OnePlusWOverA(a_c, a_c, nu, w_i, w_f));
}

double AxionEDEFluid::OmegaAcFromOmegaZero(double omega0,
                                           double a_c,
                                           double nu,
                                           double w_i,
                                           double w_f) {
  return omega0 * std::exp(Integral3OnePlusWOverA(a_c, a_c, nu, w_i, w_f));
}

void AxionEDEFluid::ComputeWFld(double a,
                                double* w_fld,
                                double* dw_over_da_fld,
                                double* integral_fld) const {
  const double dw_tot = w_f_ - w_i_;
  const double r      = 3. * dw_tot / nu_;

  // a = 0 is a real call (background_checks probes w(a->0)); the sigmoid limit
  // is exact there: w = w_i, dw/da = 0. Large t covers underflow near a = 0 too.
  const double t = (a > 0.) ? r * std::log(a_c_ / a) : std::numeric_limits<double>::infinity();
  if (t > 700.) {
    *w_fld          = w_i_;
    *dw_over_da_fld = 0.;
    *integral_fld   = (a > 0.) ? Integral3OnePlusWOverA(a, a_c_, nu_, w_i_, w_f_) : 0.;
    return;
  }
  const double x = std::exp(t);  // (a_c/a)^r
  *w_fld          = w_i_ + dw_tot / (1. + x);
  *dw_over_da_fld = dw_tot * (r / a) * x / ((1. + x) * (1. + x));
  *integral_fld   = Integral3OnePlusWOverA(a, a_c_, nu_, w_i_, w_f_);
}

double AxionEDEFluid::Cs2(double k2, double a) const {
  // GDM effective sound speed, Poulin et al. 1806.10608 (AxiCLASS perturbations.c:7802).
  const double wbar = omega_axion_ * std::pow(a, -3. * (n_axion_ - 1.) / (n_axion_ + 1.));
  const double t    = 2. * a * a * wbar * wbar;
  return ((n_axion_ - 1.) * t + k2) / ((n_axion_ + 1.) * t + k2);
}

void AxionEDEFluid::ApplyInitialConditions(const BaseSpecies::PerturbLayout& base,
                                           double* y,
                                           const PerturbIcContext& ctx) {
  // EDE-adapted adiabatic ICs (AxiCLASS perturbations.c:5784-5795): regular as
  // w -> -1 because delta, theta carry explicit (1+w) suppression... delta does;
  // theta's (1+w) lives in rho_plus_p_theta downstream.
  const auto& layout = static_cast<const FluidSpecies::PerturbLayout&>(base);
  if (ctx.index_ic != ctx.p_mod->index_ic_ad_)
    return;
  if (layout.idx_delta < 0 || layout.idx_theta < 0)
    return;

  double w_fld, dw_over_da, integral;
  ComputeWFld(ctx.a, &w_fld, &dw_over_da, &integral);
  const double cs2   = Cs2(ctx.k * ctx.k, ctx.a);
  const double denom = 32. + 6. * cs2 + 12. * w_f_;

  y[layout.idx_delta] = 0.5 * ctx.ktau_two * (1. + w_fld) * (-4. + 3. * cs2) / denom *
                        ctx.ppr->curvature_ini * ctx.s2_squared;
  y[layout.idx_theta] = -0.5 * ctx.k * ctx.ktau_three * cs2 / denom *
                        ctx.ppr->curvature_ini * ctx.s2_squared;
}

void AxionEDEFluid::SetBackgroundModule(const BackgroundModule* bgm) {
  FluidSpecies::SetBackgroundModule(bgm);
  DeriveAxionScales();
}

void AxionEDEFluid::DeriveAxionScales() {
  // Port of AxiCLASS background.c 982-1023 (Eqs. 27, 28, 30 of 1806.10608).
  // E(a_c)^2 = sum over the other species scaled to a_c by EnergyType
  // (Radiation a^-4, Matter a^-3, DarkEnergy flat) + our own Omega_fld_ac.
  // Composites and EnergyType::Other are skipped, matching the approximate
  // nature of AxiCLASS's E(a_c) (this only calibrates omega_axion, not rho_fld).
  double omega_r = 0., omega_m = 0., omega_de = 0.;
  for (const auto& sp : bgm_->all_species_) {
    if (sp.get() == static_cast<const BaseSpecies*>(this))
      continue;
    switch (sp->energy_type()) {
      case EnergyType::Radiation:
        omega_r += sp->GetOmega0();
        break;
      case EnergyType::Matter:
        omega_m += sp->GetOmega0();
        break;
      case EnergyType::DarkEnergy:
        omega_de += sp->GetOmega0();
        break;
      case EnergyType::Other:
        break;
    }
  }

  const double n = n_axion_;
  if (n > 50.) {  // AxiCLASS guard: the Gamma-function factor degenerates
    m_fld_       = 0.;
    alpha_fld_   = 0.;
    omega_axion_ = 0.;
    return;
  }

  const double a_eq = omega_r / omega_m;
  const double p    = (a_c_ < a_eq) ? 0.5 : 2. / 3.;
  const double Eac  = std::sqrt(omega_r * std::pow(a_c_, -4.) + omega_m * std::pow(a_c_, -3.) +
                                omega_de + omega_ac_);
  const double xc          = p / Eac;
  const double f           = 7. / 8.;
  const double cos_initial = std::cos(theta_i_);
  const double sin_initial = std::sin(theta_i_);

  m_fld_ = std::pow(1. - cos_initial, (1. - n) / 2.) *
           std::sqrt((1. - f) * (6. * p + 2.) * theta_i_ / (n * sin_initial)) / xc;
  alpha_fld_ = std::sqrt(6. * omega_ac_) / m_fld_ / std::pow(1. - cos_initial, n / 2.);
  const double Gac = std::sqrt(M_PI) * std::tgamma((n + 1.) / (2. * n)) /
                     std::tgamma(1. + 1. / (2. * n)) *
                     std::pow(2., -(n * n + 1.) / (2. * n)) *
                     std::pow(3., 0.5 * (1. / n - 1.)) * std::pow(a_c_, 3. - 6. / (1. + n)) *
                     std::pow(std::pow(a_c_, 6. * n / (1. + n)) + 1., 0.5 * (1. / n - 1.));
  omega_axion_ = pba_ref_.H0 * m_fld_ * std::pow(1. - cos_initial, 0.5 * (n - 1.)) * Gac;
}
```

Implementation notes:
- `pba_ref_`: `FluidSpecies::pba_` is private, so `AxionEDEFluid` keeps its own reference. Add `const background& pba_ref_;` to the private members in `axion_ede_fluid.h` and initialize it with `pba_ref_(pba)` in the constructor's init list (after the base-class call).
- `theta_i_ = pi` would make `sin_initial = 0` → division by zero in `m_fld_`: Task 6 validates Θ ∈ (0, π) at input time.
- `bgm_->all_species_` iteration and `energy_type()` follow the exact pattern of `source/background_module.cpp:468` and `species/base_species.h:110`.

- [ ] **Step 4: Run tests**

```bash
cmake --build build -j8 && ctest --test-dir build -R test-axion-ede-fluid --output-on-failure
```
Expected: PASS.

- [ ] **Step 5: Full regression + commit**

```bash
ctest --test-dir build --output-on-failure
git add species/axion_ede_fluid.h species/axion_ede_fluid.cpp species/axion_ede_fluid_test.cpp CMakeLists.txt
git commit -m "Add AxionEDEFluid: pheno-axion EDE fluid with GDM sound speed"
```
(Include the Global Constraints trailer.)

---

### Task 6: pheno_axion wiring in FluidSpecies::CreateAll

**Files:**
- Modify: `species/fluid.h` (declare `CreatePhenoAxion`)
- Modify: `species/fluid.cpp:424-504` (`CreateAll` dispatch + new factory)
- Test: `species/axion_ede_fluid_test.cpp` (extend)

**Interfaces:**
- Consumes: `AxionEDEFluid` constructor + statics from Task 5.
- Produces: input vocabulary `fluid_equation_of_state = pheno_axion`, `n_pheno_axion` XOR `w_fld_f`, `a_c` XOR `log10_axion_ac`, exactly one of `Omega_fld`/`Omega_fld_ac`/`fraction_fld_ac`, `Theta_initial_fld` (required), `nu_fld` (default 1), `w_fld_i` (default −1).

- [ ] **Step 1: Write the failing tests**

Append to `main()` in `species/axion_ede_fluid_test.cpp` (add includes `#include "parser.h"`, `#include "species_build_context.h"`, `#include <stdexcept>`, `#include <string>`):

```cpp
  // ════ CreateAll wiring ════
  auto make_base_fc = [] {
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("n_pheno_axion", "3");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Theta_initial_fld", "2.8");
    return fc;
  };
  auto run_create = [&pba](FileContent& fc) {
    SpeciesBuildContext ctx{};
    ctx.pfc = &fc;
    ctx.pba = &pba;
    return FluidSpecies::CreateAll(ctx);
  };

  // Valid: Omega_fld_ac path.
  {
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    auto result = run_create(fc);
    assert(result.size() == 1);
    auto* ax = dynamic_cast<AxionEDEFluid*>(result[0].species.get());
    assert(ax != nullptr);
    const double ac  = std::pow(10., -3.5);
    const double om0 = AxionEDEFluid::OmegaZeroFromOmegaAc(0.05, ac, 1., -1., 0.5);
    assert(std::fabs(ax->GetOmega0() - om0) < 1e-12 * om0);
    assert(std::fabs(ax->a_c() - ac) < 1e-15);
    assert(std::fabs(ax->n_axion() - 3.) < 1e-15);
  }

  // Valid: direct Omega_fld path.
  {
    FileContent fc = make_base_fc();
    fc.set("Omega_fld", "1e-6");
    auto result = run_create(fc);
    assert(result.size() == 1);
    assert(std::fabs(result[0].species->GetOmega0() - 1e-6) < 1e-18);
  }

  // Valid: w_fld_f instead of n (n derived as (1+w_f)/(1-w_f)).
  {
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("w_fld_f", "0.5");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Theta_initial_fld", "2.8");
    fc.set("Omega_fld_ac", "0.05");
    auto result = run_create(fc);
    assert(result.size() == 1);
    auto* ax = dynamic_cast<AxionEDEFluid*>(result[0].species.get());
    assert(std::fabs(ax->n_axion() - 3.) < 1e-12);
  }

  // Rejections.
  auto expect_create_throw = [&](FileContent fc) {
    bool threw = false;
    try {
      run_create(fc);
    }
    catch (const std::invalid_argument&) {
      threw = true;
    }
    assert(threw);
  };
  {  // both n and w_fld_f
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("w_fld_f", "0.5");
    expect_create_throw(std::move(fc));
  }
  {  // missing Theta_initial_fld
    FileContent fc;
    fc.set("fluid_equation_of_state", "pheno_axion");
    fc.set("n_pheno_axion", "3");
    fc.set("log10_axion_ac", "-3.5");
    fc.set("Omega_fld_ac", "0.05");
    expect_create_throw(std::move(fc));
  }
  {  // both a_c and log10_axion_ac
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("a_c", "3e-4");
    expect_create_throw(std::move(fc));
  }
  {  // two density keys
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("Omega_fld", "1e-6");
    expect_create_throw(std::move(fc));
  }
  {  // no density key
    FileContent fc = make_base_fc();
    expect_create_throw(std::move(fc));
  }
  {  // explicit use_ppf = yes
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("use_ppf", "yes");
    expect_create_throw(std::move(fc));
  }
  {  // cs2_fld makes no sense here
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("cs2_fld", "1.0");
    expect_create_throw(std::move(fc));
  }
  {  // Theta out of range
    FileContent fc = make_base_fc();
    fc.set("Omega_fld_ac", "0.05");
    fc.set("Theta_initial_fld", "3.15");
    expect_create_throw(std::move(fc));
  }
```

Add `#include "axion_ede_fluid.h"` is already present. Note the `fraction_fld_ac` path is validated at .ini level in Task 7 (it needs a populated species collection).

- [ ] **Step 2: Run test to verify it fails**

```bash
cmake --build build -j8 --target test-axion-ede-fluid 2>&1 | tail -5
```
Expected: FAILURE — `CreateAll` throws `incomprehensible input 'pheno_axion'`.

- [ ] **Step 3: Implement**

In `species/fluid.h`, add to `FluidSpecies`'s public section (near `CreateAll`):

```cpp
  // Factory branch for fluid_equation_of_state = pheno_axion (AxiCLASS EDE fluid).
  static std::vector<Named> CreatePhenoAxion(const SpeciesBuildContext& ctx);
```

In `species/fluid.cpp` `CreateAll`: the eos string must be inspected BEFORE the Omega_fld/closure gate (pheno_axion runs may specify density via `Omega_fld_ac`/`fraction_fld_ac` with no `Omega_fld` key). Restructure the top of `CreateAll` (lines 424-452) to:

```cpp
std::vector<Named> FluidSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  // ── fluid_equation_of_state (string-keyed enum) ────────────────────────────
  equation_of_state fluid_eos = CLP;
  if (auto eos_opt = ctx.pfc->get<std::string>("fluid_equation_of_state")) {
    const std::string& eos_str = *eos_opt;
    if (eos_str.find("pheno_axion") != std::string::npos ||
        eos_str.find("PhenoAxion") != std::string::npos) {
      return CreatePhenoAxion(ctx);
    }
    else if (eos_str.find("CLP") != std::string::npos ||
             eos_str.find("clp") != std::string::npos) {
      fluid_eos = CLP;
    }
    else if (eos_str.find("EDE") != std::string::npos || eos_str.find("ede") != std::string::npos) {
      fluid_eos = EDE;
    }
    else {
      throw std::invalid_argument("incomprehensible input '" + eos_str +
                                  "' for the field 'fluid_equation_of_state'");
    }
  }

  // Resolve Omega0_fld: either the closure-budget override or pfc.
  double omega0_fld = 0.;
  if (ctx.omega0_closure_override.has_value()) {
    omega0_fld = *ctx.omega0_closure_override;
  }
  else {
    omega0_fld = ctx.pfc->get_or("Omega_fld", omega0_fld);
  }
  if (omega0_fld == 0.)
    return result;
  ...   // rest of the existing body unchanged from the "── numeric params" comment on
```

CAREFUL: the eos parse must not consume `fluid_equation_of_state` when the fluid ends up absent — but `FileContent::get` marks keys read regardless, and a user who sets `fluid_equation_of_state` without any density key gets a silent no-fluid... That is the existing behavior for `w0_fld` etc. (keys read, fluid absent) — preserve it, except pheno_axion which throws on missing density (explicitly better).

Then add the factory (after `CreateAll`):

```cpp
std::vector<Named> FluidSpecies::CreatePhenoAxion(const SpeciesBuildContext& ctx) {
  std::vector<Named> result;

  if (ctx.omega0_closure_override.has_value())
    throw std::invalid_argument(
        "the pheno_axion fluid cannot be the budget-closure species (its density is set "
        "by a_c and the EDE fraction); let Lambda close the budget");

  // ── n XOR w_fld_f ──────────────────────────────────────────────────────────
  const auto n_opt  = ctx.pfc->get<double>("n_pheno_axion");
  const auto wf_opt = ctx.pfc->get<double>("w_fld_f");
  if (n_opt && wf_opt)
    throw std::invalid_argument("give only one of 'n_pheno_axion' and 'w_fld_f'");
  if (!n_opt && !wf_opt)
    throw std::invalid_argument("pheno_axion needs 'n_pheno_axion' (or 'w_fld_f')");
  double n, w_f;
  if (n_opt) {
    n = *n_opt;
    if (n < 1.)
      throw std::invalid_argument("n_pheno_axion must be >= 1");
    w_f = AxionEDEFluid::WFinal(n);
  }
  else {
    w_f = *wf_opt;
    if (w_f <= -1. || w_f >= 1.)
      throw std::invalid_argument("w_fld_f must lie in (-1, 1)");
    n = (1. + w_f) / (1. - w_f);
  }
  const double w_i = ctx.pfc->get_or("w_fld_i", -1.);
  const double nu  = ctx.pfc->get_or("nu_fld", 1.);
  if (nu <= 0.)
    throw std::invalid_argument("nu_fld must be > 0");
  if (w_f <= w_i)
    throw std::invalid_argument("pheno_axion needs w_fld_f > w_fld_i");

  // ── a_c XOR log10_axion_ac ─────────────────────────────────────────────────
  const auto ac_opt    = ctx.pfc->get<double>("a_c");
  const auto log_ac_opt = ctx.pfc->get<double>("log10_axion_ac");
  if (ac_opt && log_ac_opt)
    throw std::invalid_argument("give only one of 'a_c' and 'log10_axion_ac'");
  if (!ac_opt && !log_ac_opt)
    throw std::invalid_argument("pheno_axion needs 'a_c' or 'log10_axion_ac'");
  const double a_c = ac_opt ? *ac_opt : std::pow(10., *log_ac_opt);
  if (a_c <= 0. || a_c >= 1.)
    throw std::invalid_argument("a_c must lie in (0, 1)");

  // ── Theta_initial_fld (required; enters only the omega_axion calibration) ──
  const auto theta_opt = ctx.pfc->get<double>("Theta_initial_fld");
  if (!theta_opt || *theta_opt <= 0. || *theta_opt >= M_PI)
    throw std::invalid_argument("pheno_axion requires Theta_initial_fld in (0, pi)");

  // ── density: exactly one of Omega_fld / Omega_fld_ac / fraction_fld_ac ────
  const auto om0_opt  = ctx.pfc->get<double>("Omega_fld");
  const auto omac_opt = ctx.pfc->get<double>("Omega_fld_ac");
  const auto frac_opt = ctx.pfc->get<double>("fraction_fld_ac");
  const int n_density = int(om0_opt.has_value()) + int(omac_opt.has_value()) +
                        int(frac_opt.has_value());
  if (n_density != 1)
    throw std::invalid_argument(
        "pheno_axion needs exactly one of 'Omega_fld', 'Omega_fld_ac', 'fraction_fld_ac'");

  double omega0_fld = 0.;
  if (om0_opt) {
    omega0_fld = *om0_opt;
  }
  else if (omac_opt) {
    omega0_fld = AxionEDEFluid::OmegaZeroFromOmegaAc(*omac_opt, a_c, nu, w_i, w_f);
  }
  else {
    // fraction f = rho_fld(a_c)/rho_tot(a_c): Omega_fld_ac = Omega_tot_others(a_c) f/(1-f),
    // with the other species scaled to a_c by EnergyType (AxiCLASS input.c:4157-4174;
    // Lambda is excluded there too — negligible at EDE-era a_c).
    const double f_ac = *frac_opt;
    if (f_ac <= 0. || f_ac >= 1.)
      throw std::invalid_argument("fraction_fld_ac must lie in (0, 1)");
    if (ctx.all_species == nullptr)
      throw std::invalid_argument(
          "fraction_fld_ac needs the species collection (internal: all_species missing)");
    double omega_tot_ac = 0.;
    for (const auto& sp : *ctx.all_species) {
      switch (sp->energy_type()) {
        case BaseSpecies::EnergyType::Radiation:
          omega_tot_ac += sp->GetOmega0() * std::pow(a_c, -4.);
          break;
        case BaseSpecies::EnergyType::Matter:
          omega_tot_ac += sp->GetOmega0() * std::pow(a_c, -3.);
          break;
        default:
          break;
      }
    }
    const double omega_ac = omega_tot_ac * f_ac / (1. - f_ac);
    omega0_fld            = AxionEDEFluid::OmegaZeroFromOmegaAc(omega_ac, a_c, nu, w_i, w_f);
  }
  if (omega0_fld <= 0.)
    return result;

  // ── PPF is meaningless here; cs2 is derived, not an input ──────────────────
  if (auto ppf_opt = ctx.pfc->get<std::string>("use_ppf")) {
    if (ppf_opt->find("y") != std::string::npos || ppf_opt->find("Y") != std::string::npos)
      throw std::invalid_argument(
          "use_ppf = yes is incompatible with pheno_axion (w > -1 always; the GDM "
          "sound speed requires the true fluid equations)");
  }
  if (ctx.pfc->get<double>("cs2_fld").has_value())
    throw std::invalid_argument(
        "cs2_fld is not an input for pheno_axion (the sound speed is the GDM formula "
        "derived from n, a_c and Theta_initial_fld)");

  result.push_back({"Fluid",
                    std::make_unique<AxionEDEFluid>(*ctx.pba, omega0_fld, a_c, n, nu, w_i, w_f,
                                                    *theta_opt)});
  return result;
}
```

Add `#include "axion_ede_fluid.h"` and `#include <cmath>` (already there) to `fluid.cpp`. Note the exact `BaseSpecies::EnergyType` qualification — check how `energy_type()`'s return type is spelled at `species/base_species.h:79,110` (nested enum → qualify as `BaseSpecies::EnergyType::Radiation`; from inside species code the shorter `EnergyType::` may resolve — match whatever `background_module.cpp` uses).
`SpeciesCollection` iteration: confirm range-for over `*ctx.all_species` yields `std::unique_ptr<BaseSpecies>` refs (same as `bgm_->all_species_` in `background_module.cpp:468`); adapt the loop variable if the collection yields a different element type.

- [ ] **Step 4: Run tests**

```bash
cmake --build build -j8 && ctest --test-dir build -R test-axion-ede-fluid --output-on-failure
```
Expected: PASS.

- [ ] **Step 5: Full regression + commit**

```bash
ctest --test-dir build --output-on-failure
git add species/fluid.h species/fluid.cpp species/axion_ede_fluid_test.cpp
git commit -m "Wire fluid_equation_of_state=pheno_axion into FluidSpecies::CreateAll"
```
(Include the Global Constraints trailer.)

---

### Task 7: EDE integration fixture + physics validation

**Files:**
- Create: `test_axion_ede.ini` (repo root, `git add -f`)
- Scratch scripts in the session scratchpad (not committed)

**Interfaces:**
- Consumes: everything from Tasks 4-6, end-to-end.

- [ ] **Step 1: Write the fixture**

Create `test_axion_ede.ini`:

```ini
# Pheno-axion early dark energy: f_EDE ~ 0.1 at z_c ~ 10^3.5 with n = 3 (w_f = 1/2).
# Canonical EDE configuration from Poulin et al. 1811.04083.
h = 0.67
omega_b = 0.02238
omega_cdm = 0.1201
output = tCl, mPk
fluid_equation_of_state = pheno_axion
n_pheno_axion = 3
log10_axion_ac = -3.5
fraction_fld_ac = 0.1
Theta_initial_fld = 2.8
write background = yes
root = output/test_axion_ede_
background_verbose = 1
input_verbose = 1
```

- [ ] **Step 2: Run + verify the EDE fraction**

```bash
./class test_axion_ede.ini
```
Expected: exit 0. Then verify (reuse the header-parsing prologue from Task 3 step 3, same `N:name` format, file `output/test_axion_ede_background.dat`):

```bash
python3 - <<'EOF'
import re
import numpy as np

path = 'output/test_axion_ede_background.dat'
names = None
with open(path) as fh:
    for line in fh:
        if line.startswith('#'):
            toks = re.findall(r'(\d+):(\S+)', line)
            if toks:
                names = [t[1] for t in toks]
        else:
            break
data = np.loadtxt(path)
col = {n: i for i, n in enumerate(names)}
z, rho_fld, w_fld, H = (data[:, col['z']], data[:, col['(.)rho_fld']],
                        data[:, col['(.)w_fld']], data[:, col['H']])
a = 1.0 / (1.0 + z)
a_c = 10**-3.5

# 1-2: EDE fraction at a_c. In CLASS units rho_crit(a) = H(a)^2.
f_ede = rho_fld / H**2
f_at_ac = np.interp(a_c, a[np.argsort(a)], f_ede[np.argsort(a)])
assert abs(f_at_ac - 0.1) < 0.015, f_at_ac

# 3: w crosses (w_i+w_f)/2 = -0.25 at a_c (within 5% in a).
order = np.argsort(a)
a_s, w_s = a[order], w_fld[order]
i_cross = np.argmax(w_s > -0.25)
assert abs(a_s[i_cross] / a_c - 1) < 0.05, a_s[i_cross]

# 4: w = -1 plateau below a_c; rho ~ a^-4.5 tail one decade above a_c.
plateau = rho_fld[a < a_c / 30]
assert plateau.max() / plateau.min() - 1 < 0.01
tail = (a_s > 10 * a_c) & (a_s < 30 * a_c)
slope = np.polyfit(np.log(a_s[tail]), np.log(np.interp(a_s[tail], a_s, rho_fld[order])), 1)[0]
assert abs(slope + 4.5) < 0.25, slope
print("EDE background OK: f_EDE(a_c) =", f_at_ac, "tail slope =", slope)
EOF
```
Expected: prints `EDE background OK` with f_EDE(a_c) ≈ 0.1 (the fraction→Omega conversion is exact for the fld itself; the residual comes from the neutrino-scaling approximation in Ω_tot(a_c)). Adjust column-name strings to the actual header; do not change the assertions.

- [ ] **Step 3: ΛCDM decoupling check**

Copy the fixture to the scratchpad with `fraction_fld_ac = 1e-8` and root `output/test_axion_tiny_`; also create a no-fluid version (delete all fluid keys) with root `output/test_axion_none_`. Run both. Compare Cl^TT (and P(k) if written): assert max over l≥2 of |Cl_tiny/Cl_none − 1| < 1e-4.
Expected: passes — the fluid decouples smoothly. (These two .ini files are scratch, not committed.)

- [ ] **Step 4: Recombination-code cross-check**

Run `test_axion_ede.ini` once with each recombination code (find the exact key and value spellings by searching `explanatory.ini` for `recombination`; set it explicitly in two scratch copies of the fixture). Compare Cl^TT: agreement < 0.5% for l ≥ 2. This guards against the HyRec w0-coupling worry (`thermodynamics_module.cpp:2901` reads `background_w_fld(1., ...)`): pheno_axion gives w(a=1) ≈ w_f ≈ 0.5, and if HyRec misused it the two codes would diverge visibly.
Expected: sub-0.5% agreement. If it fails, STOP and report — do not patch around it; the fix likely belongs in how thermodynamics reads w_fld and needs a human-visible decision. Document the outcome either way for the PR description.

- [ ] **Step 5: KG-vs-fluid cross-validation (advisory, goes in the PR text)**

In the scratchpad, configure a **late** transition where the exact KG run is affordable:
- scf run: `scf_potential = axion`, `n_axion = 2`, `f_axion = 0.1`, `Theta_initial_scf = 2.0`, `Omega_scf = 0.05`, `write background = yes`.
- From its background.dat, measure: a_half = scale factor where w_scf (= p_scf/rho_scf) first crosses (w_i+w_f)/2 sustainably, and the density fraction at that point.
- fluid run: pheno_axion with `a_c = a_half`, `fraction_fld_ac` = measured fraction, `n_pheno_axion = 2`, `Theta_initial_fld = 2.0`.
- Compare rho_fld(a) vs rho_scf(a) (cycle-averaged for the scf after a_c): report the max deviation over a ∈ [1e-3, 1]; ~10-20% is expected (the sigmoid is an approximation to the true thaw).
Record the numbers; no hard assert. If the scf run is too slow or oscillations dominate, reduce to a qualitative statement (frozen plateau matches; post-thaw slopes match) — this is an advisory check.

Optional, timeboxed to ~20 minutes: attempt `make -j8` inside `$AXICLASS` (old C CLASS 2.x — likely fails on modern clang; if so, skip immediately) and, if it builds, run one matched pheno-axion configuration and compare Cl^TT at ~1% against ours. Record the outcome (built or skipped) for the PR body.

- [ ] **Step 6: Commit**

```bash
git add -f test_axion_ede.ini
git commit -m "Add pheno-axion EDE integration fixture with f_EDE(a_c) validation"
```
(Include the Global Constraints trailer; summarize the step 3-5 results in the commit body.)

---

### Task 8: Documentation, follow-up issues, PR

**Files:**
- Modify: `explanatory.ini` (document both key sets, near the existing scf/fld blocks)
- No other source changes.

- [ ] **Step 1: Document input keys in explanatory.ini**

Find the scalar-field parameter block (search `scf_parameters`) and the fluid block (search `fluid_equation_of_state`) in `explanatory.ini`; append commented documentation for: `scf_potential = axion`, `m_axion`, `f_axion`, `n_axion`, `Theta_initial_scf`; and `fluid_equation_of_state = pheno_axion`, `n_pheno_axion`/`w_fld_f`, `a_c`/`log10_axion_ac`, `Omega_fld`/`Omega_fld_ac`/`fraction_fld_ac`, `Theta_initial_fld`, `nu_fld`, `w_fld_i`. State units (m in 1/Mpc, f in reduced-Planck units) and the shooting semantics (m tuned to hit Omega_scf). Follow the surrounding comment style exactly.

- [ ] **Step 2: Final full verification**

```bash
cmake --build build -j8 && ctest --test-dir build --output-on-failure
./class test_axion_scf.ini && ./class test_axion_ede.ini
./class explanatory.ini
```
Expected: all pass; explanatory.ini (no axion keys active) unchanged behavior.

- [ ] **Step 3: Commit + file follow-up issues**

```bash
git add explanatory.ini
git commit -m "Document axion input keys in explanatory.ini"
```
(Include the Global Constraints trailer.) Then file follow-up issues with `gh issue create` (one each):
1. "Axion scf: KG→fluid runtime switching (AxiCLASS scf_evolve_as_fluid parity)" — body: needs a background approximation-switch seam; perturbation side can reuse CopyPerturbationsAcrossSwitch; reference the spec's Approach B analysis.
2. "Axion scf inside the Type-3 composite" — CreateAllForComposite hardwires the 1EXP potential.
3. "pheno_axion: derived-parameter output (f_EDE peak, a_peak; m_fld/alpha_fld) via wrapper".
4. "Axion potentials phi_2n / axionquad via the ScalarFieldPotential seam".

- [ ] **Step 4: Push and open the PR**

```bash
git push -u origin axion-species
gh pr create --title "Axion support: (1-cos)^n scf potential + pheno-axion EDE fluid" --body "<summary per below>"
```
PR body: link the spec and plan docs; summarize both components; paste the Task 7 validation numbers (f_EDE(a_c) accuracy, decoupling check, recfast/hyrec agreement, KG-vs-fluid comparison); list the follow-up issues; note the two AxiCLASS bugs intentionally not ported. End with the standard PR trailer:

```
🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_01EWPTcYe8AfdxcSqmx4b9ye
```
