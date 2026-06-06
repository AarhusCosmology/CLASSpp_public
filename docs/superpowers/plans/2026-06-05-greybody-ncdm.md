# Grey-body NCDM Species Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a grey-body phase-space-distribution NCDM flavor (`type = ncdm_greybody`) that re-implements the physics of the unmergeable PR #94 against the post-refactor species system.

**Architecture:** A `GreyBodyNCDMSpecies : public NCDMSpecies` overrides only the analytic PSD and the quadrature sampling. The PSD is made overridable via a small virtual hook in `NCDMBaseSpecies`; construction is split into two phases (a deferred-init constructor) so the grey-body PSD drives its own sampling in a single clean pass. The moment-based inverse solve (with a vendored Riemann ζ used uniformly across compilers) lives in a dedicated, unit-testable `species/greybody_moments.{h,cpp}`. Per-species M₂/M₃/M₄ are exposed through the existing `GetSpeciesParam` path into `classy`.

**Tech Stack:** C++17 (Apple clang / GCC), CLASS++ species system, plain-`assert` C++ unit tests built via Makefile targets, Python `classy` integration tests.

**Reference spec:** `docs/superpowers/specs/2026-06-05-greybody-ncdm-design.md`

---

## File Structure

| File | Responsibility |
|------|----------------|
| `species/greybody_moments.h` / `.cpp` (create) | Vendored `greybody::riemann_zeta`; `greybody::GreyBodyParams` moment→(α,x,q0) inverse solve + closed-form moments. Namespace `greybody`. Unit-testable, no CLASS deps beyond `tools/bisection.h`. |
| `species/greybody_moments_test.cpp` (create) | Unit tests for ζ and the round-trip inverse solve. |
| `species/greybody_ncdm_species.h` / `.cpp` (create) | `GreyBodyNCDMSpecies`: input parsing (both modes), PSD override, quadrature-param fill, `GreyBodyMoments()`, `GetParam` override, `CreateAll`. |
| `species/ncdm_base_species.h` / `.cpp` (modify) | Extract `EvaluatePsdAnalytic` virtual; split construction into `ReadParametersByInstance` + `BuildQuadratureAndMass`; add `DeferInit` ctor; add `FillQuadratureParams` / `DefaultQuadratureStrategy` virtuals; make `InitQuadrature` protected. |
| `species/ncdm_species.h` / `.cpp` (modify) | Extract `ResolveMassOmegaClosure`; add `DeferInit` ctor. |
| `species/all_species.h` (modify) | Register `NCDMGreyBody` factory. |
| `include/quadrature.h` (modify) | Add `qm_GB_Laguerre`, `qm_trapz_log`; add `GBQuadParams` struct; thread it through `get_qsampling_manual`. |
| `tools/quadrature.cpp` (modify) | Implement the two grey-body quadrature cases. |
| `classy.pyx` (modify) | Parse `<instance>.M_2`/`.M_3`/`.M_4`/`.alpha`/`.q0`/`.x` derived params. |
| `Makefile` (modify) | `test-greybody-moments` target. |
| `test/scenarios/ncdm_greybody.ini`, `test/dotsyntax_greybody.ini` (create) | Scenario + dot-syntax examples. |
| `python/test_greybody.py` (create) | FD-limit physics test + moments round-trip via `classy`. |

---

## Phase 1 — Base refactor: overridable PSD + deferred-init construction

This phase changes only `ncdm_base_species.*` and `ncdm_species.*`, and must be **byte-stable** for all existing NCDM/DNCDM/interacting scenarios. Verify that before moving on.

### Task 1.1: Capture a regression baseline

**Files:** none (baseline artifacts only)

- [ ] **Step 1: Build current `class` and capture an NCDM baseline**

Run:
```bash
make class -j4
mkdir -p output /tmp/gb_baseline
./class test/dotsyntax_ncdm.ini
cp output/dotsyntax_ncdm_* /tmp/gb_baseline/
```
Expected: `class` builds and runs; baseline files in `/tmp/gb_baseline/`.

- [ ] **Step 2: Record a checksum**

Run: `md5 /tmp/gb_baseline/* 2>/dev/null || md5sum /tmp/gb_baseline/*`
Expected: a list of checksums — keep this output to compare after Phase 1.

### Task 1.2: Extract `EvaluatePsdAnalytic` virtual

**Files:**
- Modify: `species/ncdm_base_species.h` (add protected virtual)
- Modify: `species/ncdm_base_species.cpp:253-290` (`DistributionFunction`)

- [ ] **Step 1: Declare the virtual in the header**

In `species/ncdm_base_species.h`, in the `protected:` section (near `GetW0ForGwSource`), add:
```cpp
  /** Analytic PSD f0(q) for the non-file case. Base: Fermi-Dirac with chemical
   *  potential ksi_. GreyBodyNCDMSpecies overrides this. */
  virtual double EvaluatePsdAnalytic(double q) const;
```

- [ ] **Step 2: Implement it and route `DistributionFunction` through it**

In `species/ncdm_base_species.cpp`, replace the body of `DistributionFunction` (lines ~253-290) so the analytic branch dispatches to the virtual, and add the new method just above it:
```cpp
double NCDMBaseSpecies::EvaluatePsdAnalytic(double q) const {
  // Fermi-Dirac with chemical potential.
  return 1.0 / pow(2 * _PI_, 3) * (1. / (exp(q - ksi_) + 1.) + 1. / (exp(q + ksi_) + 1.));
}

int NCDMBaseSpecies::DistributionFunction(void* params, double q, double* f0) {
  auto* p                   = static_cast<DistributionParams*>(params);
  const NCDMBaseSpecies* sp = p->sp;

  if (p->tablesize > 0) {
    int lastidx = p->tablesize - 1;
    if (q < p->q[0]) {
      *f0 = p->f0[0];
    }
    else if (q > p->q[lastidx]) {
      double qlast   = p->q[lastidx];
      double f0last  = p->f0[lastidx];
      double dqlast  = qlast - p->q[lastidx - 1];
      double df0last = f0last - p->f0[lastidx - 1];
      *f0            = f0last * exp(-(qlast - q) * df0last / f0last / dqlast);
    }
    else {
      class_call(array_interpolate_spline(p->q.data(),
                                          p->tablesize,
                                          p->f0.data(),
                                          p->d2f0.data(),
                                          1,
                                          q,
                                          &p->last_index,
                                          f0,
                                          1,
                                          const_cast<char*>(sp->error_message_)),
                 const_cast<char*>(sp->error_message_),
                 const_cast<char*>(sp->error_message_));
    }
  }
  else {
    *f0 = sp->EvaluatePsdAnalytic(q);
  }
  return _SUCCESS_;
}
```

- [ ] **Step 3: Build**

Run: `make class -j4`
Expected: compiles with no errors.

- [ ] **Step 4: Verify byte-stable vs baseline**

Run:
```bash
./class test/dotsyntax_ncdm.ini
md5 output/dotsyntax_ncdm_* 2>/dev/null || md5sum output/dotsyntax_ncdm_*
```
Expected: checksums identical to Task 1.1 Step 2.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp
git commit -m "refactor(ncdm): extract EvaluatePsdAnalytic virtual for PSD override"
```

### Task 1.3: Split construction — `BuildQuadratureAndMass` + `DeferInit` ctor

**Files:**
- Modify: `species/ncdm_base_species.h`
- Modify: `species/ncdm_base_species.cpp:16-95,101` (ctor, `ReadParametersByInstance`, make `InitQuadrature` protected)

- [ ] **Step 1: Header — add tag, deferred ctor, protected build method; move `InitQuadrature`**

In `species/ncdm_base_species.h`:
- Add to the `protected:` section:
```cpp
  /** Tag type selecting the deferred-init constructor: reads parameters but
   *  does NOT build quadrature or mass, so a subclass can set up an overridden
   *  PSD first and then call BuildQuadratureAndMass itself. */
  struct DeferInit {};

  NCDMBaseSpecies(std::string name,
                  EnergyType energy_type,
                  FileContent* pfc,
                  const std::string& instance_name,
                  const NcdmSettings& settings,
                  DeferInit);

  /** Build perturbation/background quadrature and the dimensionless mass M_.
   *  Calls InitQuadrature (which dispatches the virtual PSD) then computes M_
   *  from m_in_eV_. Safe to call once, after the object is fully constructed. */
  void BuildQuadratureAndMass(const NcdmSettings& settings);
```
- Move the declaration `void InitQuadrature(const NcdmSettings& settings);` from the `private:` section to the `protected:` section.
- Add to `protected:` the two quadrature-customisation virtuals (defined inline):
```cpp
  /** Default quadrature strategy when the user does not override it.
   *  Base: qm_auto. GreyBodyNCDMSpecies returns a GB-specific method. */
  virtual int DefaultQuadratureStrategy() const { return qm_auto; }

  /** Populate grey-body quadrature parameters. Base: leaves p.active == false. */
  virtual void FillQuadratureParams(GBQuadParams& /*p*/) const {}
```
- Ensure `#include "quadrature.h"` is present in this header (it is, transitively via the existing include of `quadrature.h`; confirm `GBQuadParams` is visible — it will be after Phase 5 Task 5.1 adds it. **For this task, temporarily add a forward-compatible include**: `#include "quadrature.h"` already exists at the top, so just ensure Task 5.1 lands the struct. To keep Phase 1 self-compiling, add the struct stub now in `quadrature.h` — see Step 1b.)

- [ ] **Step 1b: Add the `GBQuadParams` struct + enum values now (used in Phase 5)**

In `include/quadrature.h`, change the enum line and add the struct directly below it:
```cpp
enum quadrature_method {
  qm_auto,
  qm_Laguerre,
  qm_trapz_indefinite,
  qm_trapz,
  qm_GB_Laguerre,
  qm_trapz_log
};

/** Optional grey-body parameters threaded into manual quadrature. When
 *  `active` is false the GB-specific methods are not used. Defined here (not in
 *  a species header) so quadrature.cpp stays species-agnostic. */
struct GBQuadParams {
  bool active        = false;
  double alpha       = 0.;
  double x_times_alpha = 0.;
  double alpham1_logq0 = 0.;
};
```

- [ ] **Step 2: cpp — split the constructor tail into `BuildQuadratureAndMass`**

In `species/ncdm_base_species.cpp`:
- Change `ReadParametersByInstance` so it ends at line ~82 (after the ultra-relativistic default) and **remove** the `InitQuadrature(settings);` call and the `if (m_in_eV_ != 0.0) { ... } else { M_ = 0.; }` block (lines ~84-94).
- Add the normal constructor body to call both, and add the deferred constructor:
```cpp
NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 const std::string& instance_name,
                                 const NcdmSettings& settings)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadParametersByInstance(pfc, instance_name, settings);
  BuildQuadratureAndMass(settings);
}

NCDMBaseSpecies::NCDMBaseSpecies(std::string name,
                                 EnergyType energy_type,
                                 FileContent* pfc,
                                 const std::string& instance_name,
                                 const NcdmSettings& settings,
                                 DeferInit)
    : BaseSpecies(std::move(name), energy_type), T_cmb_(settings.T_cmb), h_(settings.h) {
  ReadParametersByInstance(pfc, instance_name, settings);
  // Caller (subclass) is responsible for BuildQuadratureAndMass after setting up its PSD.
}

void NCDMBaseSpecies::BuildQuadratureAndMass(const NcdmSettings& settings) {
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

- [ ] **Step 3: Build**

Run: `make class -j4`
Expected: compiles. (`InitQuadrature` still uses `DistributionFunction`; `FillQuadratureParams` is a no-op here.)

- [ ] **Step 4: Verify byte-stable**

Run:
```bash
./class test/dotsyntax_ncdm.ini
md5 output/dotsyntax_ncdm_* 2>/dev/null || md5sum output/dotsyntax_ncdm_*
```
Expected: checksums identical to Task 1.1 Step 2 (construction order is preserved for the normal ctor).

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_base_species.h species/ncdm_base_species.cpp include/quadrature.h
git commit -m "refactor(ncdm): two-phase construction (DeferInit ctor + BuildQuadratureAndMass)"
```

### Task 1.4: Extract `ResolveMassOmegaClosure` + `DeferInit` ctor on `NCDMSpecies`

**Files:**
- Modify: `species/ncdm_species.h`
- Modify: `species/ncdm_species.cpp:12-48` (constructor)

- [ ] **Step 1: Header — declarations**

In `species/ncdm_species.h`, in the `protected:` section add:
```cpp
  // Deferred-init constructor for subclasses (e.g. GreyBodyNCDMSpecies) that
  // must configure an overridden PSD before quadrature is built.
  NCDMSpecies(FileContent* pfc,
              const std::string& instance_name,
              const NcdmSettings& settings,
              const background* pba,
              const BackgroundModule* bgm,
              NCDMBaseSpecies::DeferInit);

  // Standard-NCDM mass/Omega closure (factored out of the public constructor so
  // GreyBodyNCDMSpecies can re-run it after rebuilding quadrature).
  void ResolveMassOmegaClosure(const NcdmSettings& settings);
```

- [ ] **Step 2: cpp — refactor the constructor**

In `species/ncdm_species.cpp`, replace the constructor (lines ~12-48) with:
```cpp
NCDMSpecies::NCDMSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings), pba_(pba) {
  bgm_ = bgm;
  ResolveMassOmegaClosure(settings);
}

NCDMSpecies::NCDMSpecies(FileContent* pfc,
                         const std::string& instance_name,
                         const NcdmSettings& settings,
                         const background* pba,
                         const BackgroundModule* bgm,
                         NCDMBaseSpecies::DeferInit defer)
    : NCDMBaseSpecies(instance_name, EnergyType::Other, pfc, instance_name, settings, defer),
      pba_(pba) {
  bgm_ = bgm;
  // No closure: the subclass calls BuildQuadratureAndMass + ResolveMassOmegaClosure.
}

void NCDMSpecies::ResolveMassOmegaClosure(const NcdmSettings& settings) {
  const double H0 = settings.h * 1.e5 / _c_;
  if (m_in_eV_ != 0.0) {
    double rho_ncdm = 0.;
    ComputeMomenta(0., nullptr, &rho_ncdm, nullptr, nullptr, nullptr);
    if (GetOmega0() == 0.0) {
      SetOmega0(rho_ncdm / H0 / H0, settings.h);
    }
    else {
      const double fnu_factor  = H0 * H0 * GetOmega0() / rho_ncdm;
      factor_                 *= fnu_factor;
      deg_                    *= fnu_factor;
    }
  }
  else {
    M_       = MFromOmega(H0, GetOmega0(), settings.tol_M_ncdm);
    m_in_eV_ = _k_B_ / _eV_ * T_ * M_ * T_cmb_;
  }
}
```
Note: `m_in_eV_`, `factor_`, `deg_`, `M_` are protected members of `NCDMBaseSpecies` (already accessible). `MFromOmega`, `SetOmega0`, `ComputeMomenta`, `GetOmega0` are inherited.

- [ ] **Step 3: Build + verify byte-stable**

Run:
```bash
make class -j4
./class test/dotsyntax_ncdm.ini
md5 output/dotsyntax_ncdm_* 2>/dev/null || md5sum output/dotsyntax_ncdm_*
```
Expected: compiles; checksums identical to Task 1.1 Step 2.

- [ ] **Step 4: Run the broader scenario set to confirm DNCDM/interacting unaffected**

Run:
```bash
./class test/dotsyntax_dncdm.ini && ./class test/dotsyntax_ncdm_interacting.ini && ./class test/dotsyntax_ncdm_mixed.ini
```
Expected: all three run to completion with no error.

- [ ] **Step 5: Commit**

```bash
git add species/ncdm_species.h species/ncdm_species.cpp
git commit -m "refactor(ncdm): extract ResolveMassOmegaClosure + add DeferInit ctor"
```

---

## Phase 2 — Vendored Riemann ζ (unit-tested)

### Task 2.1: `greybody::riemann_zeta`

**Files:**
- Create: `species/greybody_moments.h`
- Create: `species/greybody_moments.cpp`
- Create: `species/greybody_moments_test.cpp`
- Modify: `Makefile` (add `test-greybody-moments`)

- [ ] **Step 1: Write the failing test**

Create `species/greybody_moments_test.cpp`:
```cpp
#include "greybody_moments.h"

#include <cassert>
#include <cmath>
#include <cstdio>

int main() {
  using greybody::riemann_zeta;
  const double PI = 3.14159265358979323846;
  // Known closed forms.
  assert(std::fabs(riemann_zeta(2.0) - PI * PI / 6.0) < 1e-9);
  assert(std::fabs(riemann_zeta(4.0) - std::pow(PI, 4) / 90.0) < 1e-9);
  // Apery's constant.
  assert(std::fabs(riemann_zeta(3.0) - 1.2020569031595942854) < 1e-9);
  // Non-integer argument vs reference value zeta(2.5) = 1.3414872573...
  assert(std::fabs(riemann_zeta(2.5) - 1.3414872573009741) < 1e-7);
  std::printf("greybody zeta tests passed\n");
  return 0;
}
```

- [ ] **Step 2: Create header + (empty-ish) implementation so it fails on assert, not link**

Create `species/greybody_moments.h`:
```cpp
#pragma once

namespace greybody {

/** Riemann zeta via the Borwein eta (alternating-series) algorithm. Vendored
 *  because libc++ lacks std::riemann_zeta; used uniformly on all compilers so
 *  results are toolchain-independent. Valid for s != 1. */
double riemann_zeta(double s);

}  // namespace greybody
```

Create `species/greybody_moments.cpp`:
```cpp
#include "greybody_moments.h"

#include <cmath>
#include <cstdlib>

namespace greybody {

double riemann_zeta(double s) {
  int n = (s > 0.) ? static_cast<int>(21 - s) : 21;
  if (n < 6) n = 6;

  // Borwein d_k coefficients via the e_j auxiliary array.
  double* e = static_cast<double*>(std::malloc(sizeof(double) * (n + 1)));
  double bnj = 1.0;
  e[n] = bnj;
  for (int j = n - 1; j >= 0; --j) {
    bnj *= (j + 1.) / (n - j);
    e[j] = e[j + 1] + bnj;
  }

  double S1 = 0., S2 = 0.;
  for (int k = 1; k <= n; ++k) {
    if ((k - 1) % 2 == 0) S1 += 1. / std::pow(k, s);
    else                  S1 -= 1. / std::pow(k, s);
  }
  for (int k = n + 1; k <= 2 * n; ++k) {
    if ((k - 1) % 2 == 0) S2 += e[k - n] / std::pow(k, s);
    else                  S2 -= e[k - n] / std::pow(k, s);
  }
  std::free(e);
  return (S1 + S2 / std::pow(2.0, n)) / (1. - std::pow(2.0, 1. - s));
}

}  // namespace greybody
```

- [ ] **Step 3: Add the Makefile target**

In `Makefile`, after the `test-bisection` target add:
```make
test-greybody-moments:
	$(CXX) $(OPTFLAG) $(CXXFLAG) -Iinclude -Itools -Isource -Ispecies -I. species/greybody_moments_test.cpp species/greybody_moments.cpp -o test-greybody-moments $(LIBRARIES)
```

- [ ] **Step 4: Build and run the test**

Run: `make test-greybody-moments && ./test-greybody-moments`
Expected: `greybody zeta tests passed`

- [ ] **Step 5: Commit**

```bash
git add species/greybody_moments.h species/greybody_moments.cpp species/greybody_moments_test.cpp Makefile
git commit -m "feat(greybody): vendored Riemann zeta (uniform across compilers) + test"
```

---

## Phase 3 — Moment inverse solve `GreyBodyParams` (unit-tested)

The grey-body PSD is `f0(q) = 2/(2π)³ (q/q0)^(α-1) / (exp(αxq)+1)`. Its dimensionless moments are
`M_n = q0^(1-α) (1 - 2^(1-(n+α))) Γ(n+α) ζ(n+α) / (αx)^(n+α)` (derived from the same building blocks as the solver; verified that `M2·M4/M3² = r_M(α)`). The solve inverts (r, M2, M3) → (α, x, q0).

### Task 3.1: `GreyBodyParams`

**Files:**
- Modify: `species/greybody_moments.h`
- Modify: `species/greybody_moments.cpp`
- Modify: `species/greybody_moments_test.cpp` (append round-trip test)

- [ ] **Step 1: Write the failing round-trip test**

Append to `main()` in `species/greybody_moments_test.cpp`, before the final `printf`/`return`:
```cpp
  // Round-trip: choose physical (alpha,x,q0), compute its moments, feed (r,M2,M3)
  // back into the solver, and confirm it recovers consistent moments.
  {
    using greybody::GreyBodyParams;
    // Forward: pick parameters and compute moments directly.
    GreyBodyParams forward = GreyBodyParams::FromDirect(2.5, 0.8, 1.3);  // alpha, x, q0
    double M2, M3, M4;
    forward.moments(M2, M3, M4);
    double r = M2 * M4 / (M3 * M3);

    // Inverse: solve from (r, M2, M3).
    GreyBodyParams solved = GreyBodyParams::FromMoments(r, M2, M3);
    double m2b, m3b, m4b;
    solved.moments(m2b, m3b, m4b);

    assert(std::fabs(solved.alpha() - 2.5) < 1e-6);
    assert(std::fabs(m2b - M2) / M2 < 1e-8);
    assert(std::fabs(m3b - M3) / M3 < 1e-8);
    assert(std::fabs(m4b - M4) / M4 < 1e-6);
  }
  std::printf("greybody moment round-trip passed\n");
```
(Place the existing `printf("greybody zeta tests passed\n");` before this block; keep one final `return 0;`.)

- [ ] **Step 2: Run to verify it fails to compile (no `GreyBodyParams` yet)**

Run: `make test-greybody-moments`
Expected: FAIL — `no member named 'GreyBodyParams' in namespace 'greybody'`.

- [ ] **Step 3: Declare `GreyBodyParams` in the header**

Append to `namespace greybody` in `species/greybody_moments.h`:
```cpp
/** Grey-body PSD parameters. Construct either directly from (alpha, x, q0) or
 *  by inverting the velocity moments (r, M2, M3). Stores the log-space
 *  combinations used by the PSD and quadrature. */
class GreyBodyParams {
 public:
  static GreyBodyParams FromDirect(double alpha, double x, double q0);
  static GreyBodyParams FromMoments(double r, double M2, double M3);

  double alpha() const { return alpha_; }
  double x() const { return x_; }
  double q0() const { return q0_; }
  double x_times_alpha() const { return x_times_alpha_; }
  double alpham1_logq0() const { return alpham1_logq0_; }

  /** Closed-form dimensionless moments of this distribution. */
  void moments(double& M2, double& M3, double& M4) const;

 private:
  GreyBodyParams() = default;
  double alpha_         = 1.;
  double x_             = 1.;
  double q0_            = 1.;
  double x_times_alpha_ = 1.;
  double alpham1_logq0_ = 0.;
};
```

- [ ] **Step 4: Implement in the cpp**

Append to `namespace greybody` in `species/greybody_moments.cpp` (add `#include "bisection.h"` and `#include <stdexcept>` at the top):
```cpp
namespace {

// J_alpha = (1 - 2^-(alpha+1)) zeta(alpha+2), with Taylor fallback near alpha=-1.
double J_alpha(double alpha) {
  if (std::fabs(alpha + 1.) < 1e-3) {
    double j = alpha + 1., j2 = j * j, j3 = j2 * j, j4 = j3 * j, j5 = j4 * j, j6 = j5 * j;
    return 0.693147180559945309 + 0.159868903742430972 * j - 0.0326862962794492996 * j2 +
           0.00156899170541551496 * j3 + 0.000749872421120475325 * j4 -
           0.000204290897136748191 * j5 + 0.0000231747166672581420 * j6;
  }
  double A = std::pow(2., alpha + 1.);
  return (1. - 1. / A) * riemann_zeta(alpha + 2.);
}

// H_alpha = (2*2^(alpha+1) - 1) zeta(alpha+3), with Taylor fallback near alpha=-2.
double H_alpha(double alpha) {
  if (std::fabs(alpha + 2.) < 1e-3) {
    double h = 2. + alpha, h2 = h * h, h3 = h2 * h, h4 = h3 * h, h5 = h4 * h, h6 = h5 * h;
    return 0.693147180559945309 + 0.640321917660632396 * h + 0.244638709603290757 * h2 +
           0.0557898423443070071 * h3 + 0.00952545989606883663 * h4 +
           0.00134002229372986064 * h5 + 0.000154346027153208059 * h6;
  }
  double A = std::pow(2., alpha + 1.);
  return (2. * A - 1.) * riemann_zeta(alpha + 3.);
}

// r_M(alpha) = M2*M4/M3^2 as a function of alpha only.
double r_M(double alpha) {
  if (alpha == -1.) {
    return 216. * std::log(2.) * riemann_zeta(3.) / std::pow(M_PI, 4);
  }
  double A = std::pow(2., alpha + 1.);
  return ((4. * A * A - 5. * A + 1.) * (alpha + 3.) * riemann_zeta(alpha + 2.) *
          riemann_zeta(alpha + 4.)) /
         (std::pow(2. * A - 1., 2) * (alpha + 2.) * std::pow(riemann_zeta(alpha + 3.), 2));
}

// Bracket the root of r_M(alpha) = r, then bisect with tools/bisection.h.
double solve_alpha(double r) {
  auto f = [&](double a) { return r_M(a) - r; };
  // Initial guess from the small-residual expansion (Q(0) constant from the PR).
  double a    = 2. * (1.43748132827993497652 - r) / (r - 1.);
  double left = a, right = a, lim = 1.;
  while (f(left) < 0.) {        // walk left toward alpha = -2
    lim *= 0.1;
    right = left;
    left  = -2. + lim;
  }
  while (f(right) > 0.) {       // walk right
    left  = right;
    right = std::fabs(right) * 10.;
  }
  // f(left) >= 0, f(right) <= 0; bisect_value expects pred(mid) true on the hi side.
  if (f(left) * f(right) > 0.) {
    throw std::runtime_error("grey-body: failed to bracket alpha for r=" + std::to_string(r));
  }
  return bisect_value(left, right, 1e-13, [&](double mid) { return f(mid) < 0.; });
}

}  // namespace

GreyBodyParams GreyBodyParams::FromDirect(double alpha, double x, double q0) {
  GreyBodyParams p;
  p.alpha_         = alpha;
  p.x_             = x;
  p.q0_            = q0;
  p.x_times_alpha_ = x * alpha;
  p.alpham1_logq0_ = (alpha - 1.) * std::log(q0);
  return p;
}

GreyBodyParams GreyBodyParams::FromMoments(double r, double M2, double M3) {
  double alpha = solve_alpha(r);
  double A     = std::pow(2., alpha + 1.);
  // x from M2/M3 closure (consistent with the moment definitions).
  double x = ((2. * A - 1.) * M2 * riemann_zeta(alpha + 3.) * std::tgamma(alpha + 3.)) /
             (2. * (A - 1.) * alpha * M3 * riemann_zeta(alpha + 2.) * std::tgamma(alpha + 2.));
  double x_times_alpha = x * alpha;
  // alpham1_logq0 from the M2 normalization (log form to avoid overflow).
  double alpham1_logq0 =
      std::log(J_alpha(alpha)) + std::lgamma(alpha + 2.) -
      (alpha + 2.) * std::log(x_times_alpha) - std::log(M2);
  double q0 = std::min(1e100, std::exp(alpham1_logq0 / (alpha - 1.)));

  GreyBodyParams p;
  p.alpha_         = alpha;
  p.x_             = x;
  p.q0_            = q0;
  p.x_times_alpha_ = x_times_alpha;
  p.alpham1_logq0_ = alpham1_logq0;
  return p;
}

void GreyBodyParams::moments(double& M2, double& M3, double& M4) const {
  auto Mn = [&](int n) {
    double s = n + alpha_;
    return std::pow(q0_, 1. - alpha_) * (1. - std::pow(2., 1. - s)) * std::tgamma(s) *
           riemann_zeta(s) / std::pow(x_times_alpha_, s);
  };
  M2 = Mn(2);
  M3 = Mn(3);
  M4 = Mn(4);
}
```
(Suppress the unused `H_alpha`/`x_` warnings only if the compiler complains; `H_alpha` is retained for the quadrature `qm_trapz_log` path in Phase 5 — if a `-Werror=unused-function` build flags it before then, mark it `[[maybe_unused]]`.)

- [ ] **Step 5: Build and run**

Run: `make test-greybody-moments && ./test-greybody-moments`
Expected:
```
greybody zeta tests passed
greybody moment round-trip passed
```

- [ ] **Step 6: Commit**

```bash
git add species/greybody_moments.h species/greybody_moments.cpp species/greybody_moments_test.cpp
git commit -m "feat(greybody): moment -> (alpha,x,q0) inverse solve + round-trip test"
```

---

## Phase 4 — `GreyBodyNCDMSpecies` skeleton + PSD (FD-limit physics test)

### Task 4.1: Create the species and register it

**Files:**
- Create: `species/greybody_ncdm_species.h`
- Create: `species/greybody_ncdm_species.cpp`
- Modify: `species/all_species.h`

- [ ] **Step 1: Header**

Create `species/greybody_ncdm_species.h`:
```cpp
#pragma once
#include <optional>
#include <string>
#include <vector>

#include "../species/greybody_moments.h"
#include "../species/ncdm_species.h"
#include "../species/species_build_context.h"

/** NCDM with a grey-body phase-space distribution:
 *  f0(q) = 2/(2π)³ (q/q0)^(α-1) / (exp(α x q) + 1).
 *  Parameterized directly by (alpha, q0, x) or by inverting velocity moments
 *  (r, M2, M3). All perturbation/background behavior is inherited from NCDMSpecies. */
class GreyBodyNCDMSpecies : public NCDMSpecies {
 public:
  GreyBodyNCDMSpecies(FileContent* pfc,
                      const std::string& instance_name,
                      const NcdmSettings& settings,
                      const background* pba,
                      const BackgroundModule* bgm);

  static std::vector<Named> CreateAll(const SpeciesBuildContext& ctx);

  // Per-species derived parameters surfaced via BackgroundModule::GetSpeciesParam.
  std::optional<double> GetParam(const std::string& name) const override;

 protected:
  double EvaluatePsdAnalytic(double q) const override;
  quadrature_method DefaultQuadratureStrategy() const override;
  void FillQuadratureParams(GBQuadParams& p) const override;

 private:
  // Computes the species' dimensionless moments from its own quadrature.
  void GreyBodyMoments(double& M2, double& M3, double& M4) const;

  greybody::GreyBodyParams gb_ = greybody::GreyBodyParams::FromDirect(1., 1., 1.);
  bool gb_ready_ = false;
};
```

- [ ] **Step 2: Implementation**

Create `species/greybody_ncdm_species.cpp`:
```cpp
#include "greybody_ncdm_species.h"

#include <cmath>
#include <stdexcept>

#include "background_module.h"
#include "species/species_input.h"

namespace {
constexpr double kGreyBodyLeadingFactor = 2.0;  // the leading 2 in f0
}

GreyBodyNCDMSpecies::GreyBodyNCDMSpecies(FileContent* pfc,
                                         const std::string& instance_name,
                                         const NcdmSettings& settings,
                                         const background* pba,
                                         const BackgroundModule* bgm)
    : NCDMSpecies(pfc, instance_name, settings, pba, bgm, NCDMBaseSpecies::DeferInit{}) {
  SpeciesInput input(pfc, instance_name);

  std::string mode;
  input.read_string("parameterization", mode);
  if (mode != "direct" && mode != "moments") {
    throw std::invalid_argument("grey-body species '" + instance_name +
                                "': parameterization must be 'direct' or 'moments'");
  }

  if (mode == "direct") {
    double alpha = 0., q0 = 0., x = 0.;
    bool ok = input.read_double("alpha", alpha);
    ok &= input.read_double("q0", q0);
    ok &= input.read_double("x", x);
    if (!ok) {
      throw std::invalid_argument("grey-body species '" + instance_name +
                                  "': direct mode requires alpha, q0, x");
    }
    gb_ = greybody::GreyBodyParams::FromDirect(alpha, x, q0);
  }
  else {
    double r = 0., M2 = 0., M3 = 0.;
    bool ok = input.read_double("r", r);
    ok &= input.read_double("M2", M2);
    ok &= input.read_double("M3", M3);
    if (!ok) {
      throw std::invalid_argument("grey-body species '" + instance_name +
                                  "': moments mode requires r, M2, M3");
    }
    gb_ = greybody::GreyBodyParams::FromMoments(r, M2, M3);

    // Optional mass-via-M2 convenience: m = m_M2 / M2.
    double m_M2 = 0.;
    if (input.read_double("m_M2", m_M2)) {
      m_in_eV_ = m_M2 / M2;
    }
    if (m_in_eV_ != 0.0 && m_in_eV_ < 0.01) {
      throw std::invalid_argument(
          "grey-body species '" + instance_name +
          "': inferred mass < 0.01 eV; moment inversion is unreliable in the low-mass limit");
    }
  }
  gb_ready_ = true;

  // If the user did not pin a quadrature strategy, use the grey-body default.
  // (quadrature_strategy_ defaults to 0 == qm_auto when unspecified.)
  // BuildQuadratureAndMass reads quadrature_strategy_; set it here if still auto.
  // NOTE: quadrature_strategy_ is private to NCDMBaseSpecies; it is consulted via
  // DefaultQuadratureStrategy() inside InitQuadrature (see Phase 5 Task 5.2).

  BuildQuadratureAndMass(settings);
  ResolveMassOmegaClosure(settings);
}

double GreyBodyNCDMSpecies::EvaluatePsdAnalytic(double q) const {
  if (!gb_ready_) {
    // During base construction (before gb_ is configured) we should never reach
    // here, because BuildQuadratureAndMass is deferred. Guard defensively.
    return NCDMSpecies::EvaluatePsdAnalytic(q);
  }
  // Log-space evaluation for numerical stability at large parameters.
  double logf = (gb_.alpha() - 1.) * std::log(q) - gb_.alpham1_logq0() -
                std::log(std::exp(gb_.x_times_alpha() * q) + 1.);
  return kGreyBodyLeadingFactor / std::pow(2. * _PI_, 3) * std::exp(logf);
}

quadrature_method GreyBodyNCDMSpecies::DefaultQuadratureStrategy() const {
  return qm_GB_Laguerre;
}

void GreyBodyNCDMSpecies::FillQuadratureParams(GBQuadParams& p) const {
  p.active        = true;
  p.alpha         = gb_.alpha();
  p.x_times_alpha = gb_.x_times_alpha();
  p.alpham1_logq0 = gb_.alpham1_logq0();
}

void GreyBodyNCDMSpecies::GreyBodyMoments(double& M2, double& M3, double& M4) const {
  const std::vector<double>& qv = q_bg_;
  const std::vector<double>& wv = w_bg_;
  M2 = M3 = M4 = 0.;
  for (size_t i = 0; i < qv.size(); ++i) {
    double q2 = qv[i] * qv[i];
    M2 += q2 * wv[i];
    M3 += q2 * qv[i] * wv[i];
    M4 += q2 * q2 * wv[i];
  }
  double pref = 0.5 * std::pow(2. * M_PI, 3);
  M2 *= pref;
  M3 *= pref;
  M4 *= pref;
}

std::optional<double> GreyBodyNCDMSpecies::GetParam(const std::string& name) const {
  if (name == "alpha") return gb_.alpha();
  if (name == "q0") return gb_.q0();
  if (name == "x") return gb_.x();
  if (name == "M_2" || name == "M_3" || name == "M_4") {
    double M2, M3, M4;
    GreyBodyMoments(M2, M3, M4);
    if (name == "M_2") return M2;
    if (name == "M_3") return M3;
    return M4;
  }
  return NCDMSpecies::GetParam(name);
}

std::vector<Named> GreyBodyNCDMSpecies::CreateAll(const SpeciesBuildContext& ctx) {
  const auto instances = ctx.pfc->instances_with("type", "ncdm_greybody");
  std::vector<Named> result;
  result.reserve(instances.size());
  for (const auto& name : instances) {
    std::string unused_type;
    ctx.pfc->read_string(name + ".type", unused_type);
    auto sp = std::make_unique<GreyBodyNCDMSpecies>(ctx.pfc, name, *ctx.ncdm_settings,
                                                    ctx.pba, ctx.bgm);
    result.push_back({name, std::move(sp)});
  }
  return result;
}
```
Note: `q_bg_`, `w_bg_`, `m_in_eV_` are protected `NCDMBaseSpecies` members. `_PI_` is defined in `background.h` (pulled via the headers). If `q0`/`x`/`m_M2` etc. need different `SpeciesInput` accessors, mirror the calls in `NCDMBaseSpecies::ReadParametersByInstance` (`input.read_double(...)`).

- [ ] **Step 3: Register the factory**

In `species/all_species.h`:
- Add near the other includes: `#include "greybody_ncdm_species.h"`
- Add to the factory table after the `NCDM` entry:
```cpp
    SpeciesFactoryEntry{"NCDMGreyBody", &GreyBodyNCDMSpecies::CreateAll},
```

- [ ] **Step 3b: Add the new objects to the build**

In `Makefile`, the `SPECIES_OPP` variable (line ~82) is an **explicit** list. Append both new objects so they link into `class`:
```make
SPECIES_OPP = ... perturb_column_writer.opp greybody_moments.opp greybody_ncdm_species.opp
```
(The generic `species/%.cpp -> build/%.opp` rule handles compilation; only the list needs the two names.)

- [ ] **Step 4: Make `GetParam` overridable confirm**

Confirm `species/base_species.h:137` declares `virtual std::optional<double> GetParam(const std::string&) const`. If `NCDMSpecies` does not already declare an override, calling `NCDMSpecies::GetParam` from the grey-body override resolves to `BaseSpecies::GetParam` (returns `std::nullopt`) — which is correct.

- [ ] **Step 5: Build**

Run: `make class -j4`
Expected: compiles and links. (Quadrature still falls through to the manual path; `qm_GB_Laguerre` is added in Phase 5 — until then, set the scenario to a working strategy, see Task 4.2. For now `DefaultQuadratureStrategy` returning `qm_GB_Laguerre` is consulted only once Phase 5 wires it; in Phase 4 `BuildQuadratureAndMass` still uses the user/auto strategy.)

- [ ] **Step 6: Commit**

```bash
git add species/greybody_ncdm_species.h species/greybody_ncdm_species.cpp species/all_species.h
git commit -m "feat(greybody): GreyBodyNCDMSpecies (PSD override, both input modes, GetParam)"
```

### Task 4.2: FD-limit physics test

With α=1, x=1, q0=1 the grey-body PSD reduces to `2/(2π)³ · 1/(e^q+1)` — i.e. Fermi-Dirac with ksi=0 and deg accounting for the factor of 2. The test asserts Ω_ncdm and ρ_ncdm match an equivalent `ncdm_standard` species.

**Files:**
- Create: `test/scenarios/ncdm_greybody.ini`
- Create: `python/test_greybody.py`

- [ ] **Step 1: Write the failing physics test**

Create `python/test_greybody.py`:
```python
"""Grey-body NCDM: FD limit and moment round-trip via classy."""
import numpy as np
from classy import Class

COMMON = {
    'output': 'mPk',
    'P_k_max_1/Mpc': 1.0,
    'N_ur': 2.0308,
    'h': 0.67556,
    'omega_b': 0.022032,
    'omega_cdm': 0.12038,
}


def _run(extra):
    cosmo = Class()
    cosmo.set({**COMMON, **extra})
    cosmo.compute()
    return cosmo


def test_greybody_fd_limit():
    # Standard NCDM neutrino.
    std = _run({
        'nu.type': 'ncdm_standard', 'nu.m': 0.06, 'nu.deg': 1.0,
        'nu.momenta_bins': 7, 'nu.quadrature_strategy': 2, 'nu.max_q': 18.0,
    })
    # Grey-body in the FD limit (alpha=1, x=1, q0=1). With ksi=0 the standard FD
    # formula already equals 2/(2pi)^3 / (e^q + 1) (particle + antiparticle), which
    # is exactly the grey-body alpha=1,x=1,q0=1 PSD -> SAME deg (1.0).
    gb = _run({
        'gb.type': 'ncdm_greybody', 'gb.parameterization': 'direct',
        'gb.alpha': 1.0, 'gb.x': 1.0, 'gb.q0': 1.0,
        'gb.m': 0.06, 'gb.deg': 1.0,
        'gb.momenta_bins': 7, 'gb.quadrature_strategy': 2, 'gb.max_q': 18.0,
    })
    om_std = std.get_current_derived_parameters(['Omega_m'])['Omega_m']
    om_gb = gb.get_current_derived_parameters(['Omega_m'])['Omega_m']
    assert abs(om_gb - om_std) / om_std < 1e-3, (om_gb, om_std)
    std.struct_cleanup(); gb.struct_cleanup()


def test_greybody_moment_roundtrip():
    # Direct params; read back the moments derived parameters and confirm they are
    # finite and positive, and that M2*M4/M3^2 is a sensible ratio (> 1).
    gb = _run({
        'gb.type': 'ncdm_greybody', 'gb.parameterization': 'direct',
        'gb.alpha': 2.5, 'gb.x': 0.8, 'gb.q0': 1.3, 'gb.m': 0.1,
        'gb.momenta_bins': 15, 'gb.quadrature_strategy': 2, 'gb.max_q': 30.0,
    })
    d = gb.get_current_derived_parameters(['gb.M_2', 'gb.M_3', 'gb.M_4'])
    M2, M3, M4 = d['gb.M_2'], d['gb.M_3'], d['gb.M_4']
    assert M2 > 0 and M3 > 0 and M4 > 0, d
    assert M2 * M4 / (M3 * M3) > 1.0, d
    gb.struct_cleanup()


if __name__ == '__main__':
    test_greybody_fd_limit()
    test_greybody_moment_roundtrip()
    print('grey-body classy tests passed')
```

- [ ] **Step 2: Wire the moments derived params into classy**

In `classy.pyx`, in the derived-parameter `get_current_derived_parameters` dispatch (the `elif name == ...` chain near line 1690), add a generic branch BEFORE the final `else: raise CosmoSevereError`:
```cython
            elif name.endswith('.M_2') or name.endswith('.M_3') or name.endswith('.M_4') \
                 or name.endswith('.alpha') or name.endswith('.q0') or name.endswith('.x'):
                instance, _, field = name.rpartition('.')
                background_module = deref(self._thisptr).GetBackgroundModule()
                value = deref(background_module).GetSpeciesParam(
                    instance.encode('utf-8'), field.encode('utf-8'))
```
(`rpartition` splits on the last dot so multi-dot instance names still work.)

- [ ] **Step 3: Build classy**

Run: `make -j4` (builds `class` and the python extension; or `pip install -e .` per project convention)
Expected: builds the `classy` extension with no errors.

- [ ] **Step 4: Run the physics test**

Run: `cd python && python test_greybody.py`
Expected: `grey-body classy tests passed`

- [ ] **Step 5: Commit**

```bash
git add classy.pyx python/test_greybody.py test/scenarios/ncdm_greybody.ini
git commit -m "feat(greybody): classy per-species moments + FD-limit physics test"
```

---

## Phase 5 — Grey-body quadrature

### Task 5.1: Thread `GBQuadParams` through `get_qsampling_manual`

**Files:**
- Modify: `include/quadrature.h` (signature)
- Modify: `tools/quadrature.cpp` (signature + new cases)
- Modify: `species/ncdm_base_species.cpp` (`InitQuadrature` call site)

- [ ] **Step 1: Extend the signature**

In `include/quadrature.h`, change `get_qsampling_manual` to take the GB params (default-constructed = inactive):
```cpp
int get_qsampling_manual(double* x,
                         double* w,
                         double* dq,
                         int N,
                         double qmax,
                         enum quadrature_method method,
                         double* qvec,
                         int qsiz,
                         int (*function)(void* params_for_function, double q, double* f0),
                         void* params_for_function,
                         const GBQuadParams& gb,
                         ErrorMsg errmsg);
```

- [ ] **Step 2: Implement the two cases in `tools/quadrature.cpp`**

**Port the exact bodies from the `grey-body` branch** (`git show origin/grey-body:tools/quadrature.cpp`, the `qm_GB_Laguerre` case ~lines 45-62 and `qm_trapz_log` case ~lines 64-125), changing only the parameter source: where the branch reads `pbadist->alpha_` / `pbadist->x_times_alpha_` / `pbadist->alpham1_logq0_`, read `gb.alpha` / `gb.x_times_alpha` / `gb.alpham1_logq0`. The sketch below shows the expected structure — verify each line against the branch rather than trusting the sketch verbatim:
```cpp
  case (qm_GB_Laguerre): {
    // Generalized Gauss-Laguerre tuned to q^(alpha-1) exp(-alpha x q).
    double* b = (double*)malloc(N * sizeof(double));
    double* c = (double*)malloc(N * sizeof(double));
    compute_Laguerre(x, w, N, gb.alpha + 1., b, c, _TRUE_);
    double scaling = 1. / gb.x_times_alpha;
    for (int i = 0; i < N; i++) {
      w[i] *= pow(x[i], -(gb.alpha + 1.));   // strip the built-in weight power
      x[i] *= scaling;                        // rescale node to alpha*x grid
      w[i] *= scaling;                        // and the measure
      double y;
      (*function)(params_for_function, x[i], &y);
      w[i] *= y;
    }
    free(b);
    free(c);
    return _SUCCESS_;
  }
  case (qm_trapz_log): {
    // Log-spaced trapezoid with analytic qmin/qmax for the non-thermal tail.
    double GB_alpha = gb.alpha;
    double GB_x     = gb.x_times_alpha / GB_alpha;
    double GB_q0    = exp(gb.alpham1_logq0 / (GB_alpha - 1.));
    double beta     = GB_alpha + 2.;
    double param4   = 1e-12;  // target tail fraction
    double qmin = pow(param4 * 2 * beta * (1 - pow(2., 1. - beta)) *
                          pow(GB_x * GB_alpha, -beta) * tgamma(beta) *
                          greybody::riemann_zeta(beta),
                      1. / beta);
    qmin = MIN(qmin, 0.1 / (GB_x * GB_alpha));
    // Solve (q/15)^4 (q/q0)^(alpha-1) exp(-alpha x q) = exp(-qmax) for the upper limit.
    double xjp1 = qmax / (GB_x * GB_alpha), xj;
    do {
      xj   = xjp1;
      xjp1 = (qmax + 4 * log(xj / qmax) + (GB_alpha - 1.) * log(xj / GB_q0)) /
             (GB_x * GB_alpha);
    } while (fabs(xjp1 - xj) > 1e-10);
    double qmax_log = xjp1;
    double lqmin = log(qmin), lqmax = log(qmax_log);
    for (int i = 0; i < N; i++) {
      double h = (lqmax - lqmin) / (N - 1);
      x[i]     = exp(lqmin + i * h);
      double y;
      (*function)(params_for_function, x[i], &y);
      w[i] = y * x[i] * h;             // dq = q d(ln q)
      if (i == 0 || i == N - 1) w[i] *= 0.5;
    }
    return _SUCCESS_;
  }
```
**ζ dependency:** `qm_trapz_log` needs `riemann_zeta` for its analytic qmin. Add `#include "../species/greybody_moments.h"` to `tools/quadrature.cpp` and link `species/greybody_moments.cpp` into the `class` build (it already compiles as part of `$(SPECIES_OPP)` once the file exists — confirm the Makefile's species object glob picks it up; if `SPECIES_OPP` is an explicit list, add `greybody_moments.opp` and `greybody_ncdm_species.opp` to it). This is the one cross-dependency from `tools/` into the grey-body unit; acceptable and isolated to this single call.

- [ ] **Step 3: Update the `InitQuadrature` call site**

In `species/ncdm_base_species.cpp` `InitQuadrature`, before the manual branch build the GB params and select the strategy via the virtual; pass `gb` into `get_qsampling_manual`:
```cpp
  GBQuadParams gb;
  FillQuadratureParams(gb);
  int strategy = quadrature_strategy_;
  if (strategy == qm_auto && gb.active) {
    strategy = DefaultQuadratureStrategy();
  }
```
Then change the `if (quadrature_strategy_ == qm_auto)` test to `if (strategy == qm_auto)`, and in the manual branch pass `(enum quadrature_method) strategy` and add `gb,` as the new argument before `error_message_`.

- [ ] **Step 4: Build**

Run: `make class -j4`
Expected: compiles and links (other `get_qsampling_manual` callers, if any, pass a default `GBQuadParams{}` — grep `get_qsampling_manual(` and update any call sites to add the `GBQuadParams{}` argument).

- [ ] **Step 5: Run the grey-body scenario with the GB strategy**

Update `test/scenarios/ncdm_greybody.ini` to NOT set `gb.quadrature_strategy` (so it uses the GB default), then:
Run: `./class test/scenarios/ncdm_greybody.ini`
Expected: runs to completion; background output written.

- [ ] **Step 6: Re-run the physics test with GB quadrature**

Edit `python/test_greybody.py` to remove the `gb.quadrature_strategy`/`gb.momenta_bins`/`gb.max_q` overrides from `test_greybody_moment_roundtrip` (let it use the GB default), then:
Run: `cd python && python test_greybody.py`
Expected: `grey-body classy tests passed`

- [ ] **Step 7: Commit**

```bash
git add include/quadrature.h tools/quadrature.cpp species/ncdm_base_species.cpp test/scenarios/ncdm_greybody.ini python/test_greybody.py
git commit -m "feat(greybody): GB-tuned quadrature (qm_GB_Laguerre, qm_trapz_log)"
```

---

## Phase 6 — Scenario + dot-syntax examples, docs

### Task 6.1: Example `.ini` files

**Files:**
- Create: `test/dotsyntax_greybody.ini`
- Ensure: `test/scenarios/ncdm_greybody.ini` (created in Phase 4/5)

- [ ] **Step 1: Direct-mode dot-syntax example**

Create `test/dotsyntax_greybody.ini` (model the header on `test/dotsyntax_ncdm.ini`, replacing the `nu1.*` block with):
```ini
output = mPk
P_k_max_1/Mpc = 1.0
h = 0.67556
omega_b = 0.022032
omega_cdm = 0.12038
N_ur = 2.0308
gb.type = ncdm_greybody
gb.parameterization = direct
gb.alpha = 2.5
gb.q0 = 1.3
gb.x = 0.8
gb.m = 0.1
root = output/dotsyntax_greybody_
```

- [ ] **Step 2: Moments-mode scenario**

Ensure `test/scenarios/ncdm_greybody.ini` includes a moments-mode block:
```ini
output = tCl,pCl,mPk
l_max_scalars = 1000
P_k_max_1/Mpc = 1.0
N_ur = 2.0308
gb.type = ncdm_greybody
gb.parameterization = moments
gb.r = 1.5
gb.M2 = 1.0
gb.M3 = 1.0
gb.m_M2 = 0.06
root = test/scenarios/out/ncdm_greybody_
```

- [ ] **Step 3: Run both**

Run:
```bash
mkdir -p output test/scenarios/out
./class test/dotsyntax_greybody.ini && ./class test/scenarios/ncdm_greybody.ini
```
Expected: both run to completion.

- [ ] **Step 4: Commit**

```bash
git add test/dotsyntax_greybody.ini test/scenarios/ncdm_greybody.ini
git commit -m "test(greybody): direct + moments example inis"
```

### Task 6.2: Documentation

**Files:**
- Modify: a user-facing parameter doc if one exists (e.g. `explanatory.ini` analogue or README section). If none, add a comment block to `test/dotsyntax_greybody.ini` documenting every field.

- [ ] **Step 1: Document the fields**

Add a header comment to `test/dotsyntax_greybody.ini` listing each `gb.*` field, both modes, and the `m_M2` convenience with the low-mass guard.

- [ ] **Step 2: Commit**

```bash
git add test/dotsyntax_greybody.ini
git commit -m "docs(greybody): document grey-body input fields"
```

---

## Phase 7 — CI wiring + final verification

### Task 7.1: Add the unit test to CI

**Files:**
- Modify: `.github/workflows/test_on_pull_request.yml`

- [ ] **Step 1: Inspect how existing C++ unit tests run in CI**

Run: `grep -n "test-bisection\|test-parser\|make test" .github/workflows/test_on_pull_request.yml`
Expected: find where `test-bisection`/`test-parser` are built+run.

- [ ] **Step 2: Add `test-greybody-moments`**

Mirror the existing pattern, adding:
```yaml
      - run: make test-greybody-moments && ./test-greybody-moments
```
in the same job/step group as the other unit tests.

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/test_on_pull_request.yml
git commit -m "ci(greybody): run greybody moments unit test"
```

### Task 7.2: Full verification sweep

- [ ] **Step 1: Unit tests**

Run: `make test-greybody-moments && ./test-greybody-moments`
Expected: both "passed" lines.

- [ ] **Step 2: Existing NCDM byte-stability (final guard)**

Run:
```bash
make class -j4 && ./class test/dotsyntax_ncdm.ini
md5 output/dotsyntax_ncdm_* 2>/dev/null || md5sum output/dotsyntax_ncdm_*
```
Expected: checksums identical to Task 1.1 Step 2.

- [ ] **Step 3: Grey-body physics**

Run: `cd python && python test_greybody.py`
Expected: `grey-body classy tests passed`

- [ ] **Step 4: Existing python test suite still green**

Run: `cd python && python -m pytest test_class.py -k ncdm -q` (or the project's standard invocation)
Expected: ncdm-related tests pass.

- [ ] **Step 5: Final commit / branch ready for PR**

```bash
git status
```
Expected: clean tree; branch `greybody-ncdm` ready to open a PR that references PR #94 as the superseded original.

---

## Notes for the implementer

- **`SpeciesInput` accessor names:** Phase 4 assumes `input.read_double("field", dst)` returns `bool` (true on hit). Confirm against `species/species_input.h` and `NCDMBaseSpecies::ReadParametersByInstance`; if the return type differs, adapt the `ok &=` pattern accordingly.
- **`_PI_` vs `M_PI`:** the codebase uses `_PI_` (from `background.h`) in PSD/normalization and `M_PI` (from `<cmath>`) in `tools/`. Match the surrounding file.
- **`riemann_zeta` cross-include:** Phase 5 introduces the only `tools/ -> species/greybody_moments.h` include. If the team prefers `tools/` stay independent of `species/`, the fallback is to promote `greybody_moments.{h,cpp}` to `tools/` at that point (the spec permits promotion once a second consumer appears — `quadrature.cpp` is that second consumer). Flag this in the PR for reviewer preference.
- **Unused `H_alpha`:** retained for parity with the PR's `x_times_alpha` derivation path; if a stricter build flags it, mark `[[maybe_unused]]` rather than deleting (it documents the closed form).
```
