# Dormand-Prince (rkdp45) Evolver Port — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port the already-implemented Dormand-Prince 4(5) explicit adaptive ODE solver from branch `PhD2024-EBH` to current master as an opt-in evolver (`evolver = 2`), and verify it reproduces the `ndf15` reference observables within the established `RTOL=1e-3` bar.

**Architecture:** Standalone solver file `tools/evolver_rkdp45.cpp` + header, conforming exactly to master's canonical evolver signature so it slots into the existing `generic_evolver` function-pointer dispatch. Selection via a new `rkdp45` value appended to `enum evolver_type` (read as integer `2`). No generic-RK-core refactor, no controller abstraction — that is the follow-on programme.

**Tech Stack:** C++17, master's `class_call`/`class_test` error macros, `MIN`/`MAX` from `common.h`, the existing `class` binary + `test/scenarios/compare_tol.py` (RTOL=1e-3, zero-crossing-aware) for verification, and `class_profiled` (PR #299, now on master) for benchmarking.

**Branch:** `feat/evolver-rkdp45` (already created, rebased on master @ `a53b251a`, spec committed).

---

## Reference: canonical evolver signature (must match exactly)

From `include/evolver_rkck.h` (`evolver_rk`) and `evolver_ndf15` — the dispatch uses
`auto generic_evolver = &evolver_ndf15;` then assigns `&evolver_rk`/`&evolver_rkdp45`,
so the ported function pointer type must be **identical**:

```cpp
int evolver_xxx(
    int (*derivs)(double x, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double minimum_variation,
    int (*evaluate_timescale)(double x, void* parameters_and_workspace, double* timescale, ErrorMsg error_message),
    double timestep_over_timescale,
    double* x_sampling,
    int x_size,
    int (*output)(double x, double y[], double dy[], int index_x, void* parameters_and_workspace, ErrorMsg error_message),
    int (*print_variables)(double x, double y[], double dy[], void* parameters_and_workspace, ErrorMsg error_message),
    ErrorMsg error_message);
```

The branch's `evolver_rkdp45` had one extra param (`std::vector<double>& t_vec_evolver`) inserted
before `error_message` — **dropped** in this port (it was only a diagnostic `push_back(t)`).

DP45 ignores `minimum_variation`, `evaluate_timescale`, `timestep_over_timescale` (it does its own
adaptive error control), exactly as a true adaptive solver should. They remain in the signature for
ABI compatibility with the dispatch.

## Step-controller bug fix (intentional deviation from the branch)

The branch's step-rejection block sets `nofailed = _FALSE_` only inside the `else` branch, which is
only entered when `nofailed` is already false — so it is **never** set false, the "halve the step on a
repeated failure" path is dead code, and every rejection uses the error-proportional formula. This is a
transcription bug versus standard `ode45`. **We fix it in this port** (per user decision): set
`nofailed = _FALSE_` on the first failure so consecutive failures fall back to halving. Task 3's code
already contains the corrected logic. Note the deviation in the PR body (Task 8).

## File structure

- **Create** `include/evolver_rkdp45.h` — declaration (mirrors `evolver_rkck.h`, `extern "C"` wrapper).
- **Create** `tools/evolver_rkdp45.cpp` — the ported solver.
- **Modify** `include/common.h` — append `rkdp45` to `enum evolver_type`.
- **Modify** `source/background_module.cpp` — add dispatch branch + include.
- **Modify** `source/perturbations_module.cpp` — add dispatch branch.
- **Modify** `source/perturbations.h` and/or `source/background.h` — add `#include "evolver_rkdp45.h"` beside existing evolver includes.
- **Modify** `Makefile` — add `evolver_rkdp45.opp` to `TOOLS`.
- **Create** `test/scenarios/compare_evolver.sh` — verification harness (ndf15 vs rkdp45 across scenarios).
- **Create** `docs/superpowers/specs/2026-06-07-evolver-rkdp45-results.md` — recorded correctness + benchmark results (written in the final task).

---

### Task 1: Add the `rkdp45` enum value

**Files:**
- Modify: `include/common.h` (the `enum evolver_type` near line 428)

- [ ] **Step 1: Read the enum**

Run: `grep -n -A4 "enum evolver_type" include/common.h`
Expected: `enum evolver_type { rk, ndf15 };` (with comments).

- [ ] **Step 2: Append `rkdp45`**

Change:
```c
enum evolver_type {
  rk,   /* Runge-Kutta integrator */
  ndf15 /* stiff integrator */
};
```
to:
```c
enum evolver_type {
  rk,    /* Runge-Kutta integrator */
  ndf15, /* stiff integrator */
  rkdp45 /* Dormand-Prince 4(5) explicit adaptive integrator */
};
```
(Integer mapping: `evolver = 0` → rk, `1` → ndf15, `2` → rkdp45. `read_enum` casts the integer directly; no string map exists.)

- [ ] **Step 3: Commit**

```bash
git add include/common.h
git commit -m "feat(evolver): add rkdp45 to evolver_type enum"
```

---

### Task 2: Create the header

**Files:**
- Create: `include/evolver_rkdp45.h`

- [ ] **Step 1: Write the header**

```cpp
#ifndef __EVORKDP45__
#define __EVORKDP45__

#include "common.h"

/**************************************************************/

/**
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

int evolver_rkdp45(
    int (*derivs)(
        double x, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double minimum_variation,
    int (*evaluate_timescale)(
        double x, void* parameters_and_workspace, double* timescale, ErrorMsg error_message),
    double timestep_over_timescale,
    double* x_sampling,
    int x_size,
    int (*output)(double x,
                  double y[],
                  double dy[],
                  int index_x,
                  void* parameters_and_workspace,
                  ErrorMsg error_message),
    int (*print_variables)(
        double x, double y[], double dy[], void* parameters_and_workspace, ErrorMsg error_message),
    ErrorMsg error_message);

#ifdef __cplusplus
}
#endif

/**************************************************************/

#endif
```

- [ ] **Step 2: Commit**

```bash
git add include/evolver_rkdp45.h
git commit -m "feat(evolver): add evolver_rkdp45 header"
```

---

### Task 3: Port the solver implementation

**Files:**
- Create: `tools/evolver_rkdp45.cpp`

Port notes applied: `malloc`→`std::vector`; VLAs (`double ci[s]` with runtime `s`)→fixed arrays;
`class_evolver_output` (branch-only macro, absent on master)→`class_call`; failure `printf`+`return`
→`class_test`; dropped unused `dy`, `output_return`, dead `done` branch; dropped `t_vec_evolver`.
Step-rejection logic corrected to standard ode45 (see "Step-controller bug fix" callout).

- [ ] **Step 1: Write the implementation**

```cpp
#include "evolver_rkdp45.h"

#include <cfloat>
#include <cmath>
#include <vector>

/**
 * Dormand-Prince 4(5) explicit adaptive Runge-Kutta integrator with 4th-order
 * dense output. Ported from branch PhD2024-EBH. Conforms to the canonical CLASS
 * evolver signature (shared with evolver_ndf15 / evolver_rk).
 *
 * minimum_variation, evaluate_timescale and timestep_over_timescale are part of
 * the shared signature but unused here: this solver performs its own embedded
 * error control.
 */
int evolver_rkdp45(
    int (*derivs)(
        double x, double* y, double* dy, void* parameters_and_workspace, ErrorMsg error_message),
    double x_ini,
    double x_end,
    double* y,
    int* used_in_output,
    int y_size,
    void* parameters_and_workspace_for_derivs,
    double tolerance,
    double minimum_variation,
    int (*evaluate_timescale)(
        double x, void* parameters_and_workspace, double* timescale, ErrorMsg error_message),
    double timestep_over_timescale,
    double* x_sampling,
    int x_size,
    int (*output)(double x,
                  double y[],
                  double dy[],
                  int index_x,
                  void* parameters_and_workspace,
                  ErrorMsg error_message),
    int (*print_variables)(
        double x, double y[], double dy[], void* parameters_and_workspace, ErrorMsg error_message),
    ErrorMsg error_message) {

  (void) minimum_variation;
  (void) evaluate_timescale;
  (void) timestep_over_timescale;

  const int neq          = y_size;
  const double rtol      = tolerance;
  const double abstol    = 1e-15; /* matches ndf15 */
  const double threshold = abstol / rtol;
  const double pow_grow  = 0.2;
  const int s            = 7; /* DP45 stages */

  /* Butcher tableau (Dormand-Prince 4(5)) */
  const double ci[s] = {0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0};
  const double bi[s] =
      {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0};
  const double bi_diff[s] = {71.0 / 57600.0,
                             0.0,
                             -71.0 / 16695.0,
                             71.0 / 1920.0,
                             -17253.0 / 339200.0,
                             22.0 / 525.0,
                             -1.0 / 40.0};
  double ai[s][s] = {{0.0}};
  ai[1][0] = 0.2;
  ai[2][0] = 3.0 / 40.0;        ai[2][1] = 9.0 / 40.0;
  ai[3][0] = 44.0 / 45.0;       ai[3][1] = -56.0 / 15.0;     ai[3][2] = 32.0 / 9.0;
  ai[4][0] = 19372.0 / 6561.0;  ai[4][1] = -25360.0 / 2187.0;
  ai[4][2] = 64448.0 / 6561.0;  ai[4][3] = -212.0 / 729.0;
  ai[5][0] = 9017.0 / 3168.0;   ai[5][1] = -355.0 / 33.0;    ai[5][2] = 46732.0 / 5247.0;
  ai[5][3] = 49.0 / 176.0;      ai[5][4] = -5103.0 / 18656.0;
  ai[6][0] = 35.0 / 384.0;      ai[6][1] = 0.0;              ai[6][2] = 500.0 / 1113.0;
  ai[6][3] = 125.0 / 192.0;     ai[6][4] = -2187.0 / 6784.0; ai[6][5] = 11.0 / 84.0;

  /* 4th-order dense-output interpolation coefficients */
  const double ixx[5][3] = {{1500.0 / 371.0, -1000.0 / 159.0, 1000.0 / 371.0},
                            {-125.0 / 32.0, 125.0 / 12.0, -375.0 / 64.0},
                            {9477.0 / 3392.0, -729.0 / 106.0, 25515.0 / 6784.0},
                            {-11.0 / 7.0, 11.0 / 3.0, -55.0 / 28.0},
                            {1.5, -4.0, 2.5}};
  const double i01 = -183.0 / 64.0, i02 = 37.0 / 12.0, i03 = -145.0 / 128.0;

  std::vector<double> ynew(neq), ytemp(neq), err(neq), yinterp(neq), dyinterp(neq);
  std::vector<double> ki(s * neq);
  double bi_vec_y[s], bi_vec_dy[s];

  double t = x_ini;

  /* initialise k0 = f(t, y) */
  class_call((*derivs)(t, y, ki.data(), parameters_and_workspace_for_derivs, error_message),
             error_message,
             error_message);

  const double hmax = fabs(x_end - x_ini) / 10.0;
  double absh;
  if (x_sampling != nullptr && x_size > 1)
    absh = MIN(hmax, fabs(x_sampling[1] - x_sampling[0]));
  else
    absh = hmax;
  if (absh == 0.0)
    absh = hmax;

  /* initial step from derivative scale */
  double rh = 0.0;
  for (int k = 0; k < neq; k++) {
    double maxtmp = MAX(fabs(y[k]), threshold);
    rh            = MAX(rh, fabs(ki[k]) / maxtmp);
  }
  rh /= 0.8 * pow(rtol, pow_grow);
  if (absh * rh > 1.0)
    absh = 1.0 / rh;

  const int tdir = (x_end > x_ini) ? 1 : -1;
  double hnew    = absh * tdir;
  int nofailed   = _TRUE_;

  int idx = 0;
  if ((t - x_end) * tdir < 0.0 && x_sampling != nullptr)
    for (idx = 0; (idx < x_size) && ((x_sampling[idx] - t) * tdir < 0.0); idx++)
      ;

  while ((t - x_end) * tdir < 0.0) {
    double h          = hnew;
    const double hmin = 100.0 * DBL_MIN * fabs(t);
    class_test(fabs(h) < hmin,
               error_message,
               "rkdp45: step size %e fell below minimum %e at x=%e",
               fabs(h),
               hmin,
               t);
    if (fabs(h) > 0.9 * fabs(x_end - t))
      h = x_end - t;

    /* stage 0 (FSAL: ki[0..neq) already holds k0) */
    for (int k = 0; k < neq; k++) {
      ynew[k] = y[k] + h * bi[0] * ki[k];
      err[k]  = h * bi_diff[0] * ki[k];
    }
    for (int i = 1; i < s; i++) {
      for (int k = 0; k < neq; k++)
        ytemp[k] = y[k];
      for (int j = 0; j < i; j++)
        for (int k = 0; k < neq; k++)
          ytemp[k] += h * ai[i][j] * ki[j * neq + k];
      class_call((*derivs)(t + ci[i] * h,
                           ytemp.data(),
                           ki.data() + i * neq,
                           parameters_and_workspace_for_derivs,
                           error_message),
                 error_message,
                 error_message);
      for (int k = 0; k < neq; k++) {
        ynew[k] += h * bi[i] * ki[i * neq + k];
        err[k] += h * bi_diff[i] * ki[i * neq + k];
      }
    }

    double errmax = 0.0;
    for (int k = 0; k < neq; k++) {
      double errtemp = fabs(err[k] / MAX(threshold, fabs(ynew[k])));
      if (errtemp > errmax)
        errmax = errtemp;
    }

    if (errmax > rtol) {
      /* step rejected. Standard ode45 control: first failure shrinks by the
       * error-proportional factor and marks nofailed; consecutive failures
       * halve. (This corrects a transcription bug in the PhD2024-EBH branch
       * where nofailed was never set, making the halving path dead code.) */
      if (nofailed == _TRUE_) {
        nofailed = _FALSE_;
        hnew     = tdir * MAX(hmin, fabs(h) * MAX(0.1, 0.8 * pow(rtol / errmax, pow_grow)));
      }
      else {
        hnew = tdir * MAX(hmin, 0.5 * fabs(h));
      }
      continue;
    }

    /* step accepted */
    if (print_variables != nullptr)
      class_call((*print_variables)(t + h,
                                    ynew.data(),
                                    ki.data() + 6 * neq,
                                    parameters_and_workspace_for_derivs,
                                    error_message),
                 error_message,
                 error_message);

    nofailed          = _TRUE_;
    hnew              = tdir * MAX(hmin, fabs(h) * MAX(0.1, 0.8 * pow(rtol / errmax, pow_grow)));
    const double tnew = t + h;

    /* emit output at all sampling points within (t, tnew] */
    for (; (idx < x_size) && ((tnew - x_sampling[idx]) * tdir >= 0.0); idx++) {
      if (tnew == x_sampling[idx]) {
        class_call((*output)(tnew,
                             ynew.data(),
                             ki.data() + 6 * neq,
                             idx,
                             parameters_and_workspace_for_derivs,
                             error_message),
                   error_message,
                   error_message);
      }
      else {
        const double ti  = x_sampling[idx];
        const double ss1 = (ti - t) / h, ss2 = ss1 * ss1, ss3 = ss2 * ss1, ss4 = ss2 * ss2;
        bi_vec_y[0]  = ss1 + i01 * ss2 + i02 * ss3 + i03 * ss4;
        bi_vec_dy[0] = 1.0 + i01 * 2.0 * ss1 + i02 * 3.0 * ss2 + i03 * 4.0 * ss3;
        for (int i = 2; i < 7; i++) {
          bi_vec_y[i]  = ixx[i - 2][0] * ss2 + ixx[i - 2][1] * ss3 + ixx[i - 2][2] * ss4;
          bi_vec_dy[i] = ixx[i - 2][0] * 2 * ss1 + ixx[i - 2][1] * 3 * ss2 + ixx[i - 2][2] * 4 * ss3;
        }
        for (int k = 0; k < neq; k++) {
          if (used_in_output[k] == _TRUE_) {
            yinterp[k]  = y[k];
            dyinterp[k] = 0.0;
            for (int i = 0; i < 7; i++) {
              if (i != 1) {
                yinterp[k] += h * bi_vec_y[i] * ki[i * neq + k];
                dyinterp[k] += bi_vec_dy[i] * ki[i * neq + k];
              }
            }
          }
        }
        class_call((*output)(ti,
                             yinterp.data(),
                             dyinterp.data(),
                             idx,
                             parameters_and_workspace_for_derivs,
                             error_message),
                   error_message,
                   error_message);
      }
    }

    for (int k = 0; k < neq; k++) {
      y[k]  = ynew[k];
      ki[k] = ki[6 * neq + k]; /* FSAL: last stage becomes next k0 */
    }
    t = tnew;
  }

  return _SUCCESS_;
}
```

- [ ] **Step 2: Commit**

```bash
git add tools/evolver_rkdp45.cpp
git commit -m "feat(evolver): port Dormand-Prince 4(5) solver (rkdp45)"
```

---

### Task 4: Wire dispatch + build

**Files:**
- Modify: `Makefile` (the `TOOLS = ...` line, ~line 80)
- Modify: `source/perturbations.h` (evolver includes, ~lines 10-11)
- Modify: `source/background_module.cpp` (~lines 892-895) and its include block
- Modify: `source/perturbations_module.cpp` (~lines 2794-2797)

- [ ] **Step 1: Add to Makefile TOOLS**

In the `TOOLS = ...` line, append `evolver_rkdp45.opp`:
```make
TOOLS = growTable.opp dei_rkck.opp sparse.opp evolver_rkck.opp arrays.opp parser.opp quadrature.opp hyperspherical.opp common.opp trigonometric_integrals.opp exceptions.opp evolver_ndf15.opp evolver_rkdp45.opp
```

- [ ] **Step 2: Add the include next to the other evolver includes**

Run: `grep -n "evolver_ndf15.h\|evolver_rkck.h" source/perturbations.h source/background.h`
In `source/perturbations.h`, after the existing two evolver includes, add:
```cpp
#include "evolver_rkdp45.h"
```
Then check where `background_module.cpp` gets `evolver_ndf15`:
Run: `grep -n "evolver_ndf15\|#include" source/background_module.cpp | grep -i evolver`
If `background_module.cpp` does not transitively see the header, add `#include "evolver_rkdp45.h"` to its include block (next to where `evolver_ndf15.h` is reachable, e.g. via `background.h`). Add to `background.h` beside `dei_rkck.h` if that's where evolver decls live.

- [ ] **Step 3: Add dispatch branch in background_module.cpp**

Change (~line 892):
```cpp
  auto generic_evolver = &evolver_ndf15;
  if (ppr->evolver == rk) {
    generic_evolver = &evolver_rk;
  }
```
to:
```cpp
  auto generic_evolver = &evolver_ndf15;
  if (ppr->evolver == rk) {
    generic_evolver = &evolver_rk;
  }
  else if (ppr->evolver == rkdp45) {
    generic_evolver = &evolver_rkdp45;
  }
```

- [ ] **Step 4: Add dispatch branch in perturbations_module.cpp**

Change (~line 2794):
```cpp
    auto generic_evolver = &evolver_ndf15;
    if (ppr->evolver == rk) {
      generic_evolver = &evolver_rk;
    }
```
to:
```cpp
    auto generic_evolver = &evolver_ndf15;
    if (ppr->evolver == rk) {
      generic_evolver = &evolver_rk;
    }
    else if (ppr->evolver == rkdp45) {
      generic_evolver = &evolver_rkdp45;
    }
```

- [ ] **Step 5: Build**

Run: `make class -j4`
Expected: clean link, no errors. (If `generic_evolver = &evolver_rkdp45` fails to compile, the signature does not match `evolver_ndf15` exactly — fix the signature in Task 2/3, do not cast.)

- [ ] **Step 6: Commit**

```bash
git add Makefile source/perturbations.h source/background.h source/background_module.cpp source/perturbations_module.cpp
git commit -m "feat(evolver): wire rkdp45 into background/perturbations dispatch + build"
```

---

### Task 5: Smoke test — prove rkdp45 actually runs and is correct on one case

The two-sided check: with `evolver = 2` the raw output must **differ byte-wise** from the ndf15 run
(proving rkdp45 was actually dispatched, not silently falling back to ndf15) yet **pass compare_tol**
(proving correctness).

**Files:**
- Uses: `explanatory.ini`, `test/scenarios/compare_tol.py`, the `class` binary.

- [ ] **Step 1: Run ndf15 reference for explanatory**

```bash
mkdir -p /tmp/ev/ref /tmp/ev/new
printf "root = /tmp/ev/ref/out_\n" | cat explanatory.ini - > /tmp/ev/ref.ini
./class /tmp/ev/ref.ini
ls /tmp/ev/ref/
```
Expected: `out_*.dat` files written.

- [ ] **Step 2: Run rkdp45 (evolver = 2)**

```bash
{ cat explanatory.ini; printf "\nroot = /tmp/ev/new/out_\nevolver = 2\n"; } > /tmp/ev/new.ini
./class /tmp/ev/new.ini
ls /tmp/ev/new/
```
Expected: completes without error; `out_*.dat` files written.

- [ ] **Step 3: Confirm rkdp45 was actually used (outputs not byte-identical)**

```bash
diff -q /tmp/ev/ref/out_cl.dat /tmp/ev/new/out_cl.dat && echo "IDENTICAL (rkdp45 NOT used!)" || echo "DIFFERENT (rkdp45 ran)"
```
Expected: `DIFFERENT (rkdp45 ran)`. If IDENTICAL, dispatch wiring is wrong (Task 4) — fix before continuing.

- [ ] **Step 4: Confirm correctness within tolerance**

```bash
python3 test/scenarios/compare_tol.py /tmp/ev/ref /tmp/ev/new 'out_*.dat'
```
Expected: every file `OK`, exit 0. If FAIL, the port has a numerical bug — debug (use superpowers:systematic-debugging) before continuing.

- [ ] **Step 5: Commit (nothing to commit if all in /tmp; this is a verification gate)**

No commit. Record the smoke-test outcome for the results doc (Task 8).

---

### Task 6: Verification harness across all scenarios

**Files:**
- Create: `test/scenarios/compare_evolver.sh`

- [ ] **Step 1: Write the harness**

```bash
#!/usr/bin/env bash
# Compare rkdp45 (evolver=2) against the ndf15 default across all test scenarios.
# Pass = every scenario's outputs agree within test/scenarios/compare_tol.py (RTOL=1e-3).
# Usage: test/scenarios/compare_evolver.sh [path-to-class-binary]
set -u
CLASS="${1:-./class}"
HERE="$(cd "$(dirname "$0")" && pwd)"
WORK="$(mktemp -d)"
overall=0
for ini in "$HERE"/*.ini; do
  name="$(basename "$ini" .ini)"
  refd="$WORK/$name/ref"; newd="$WORK/$name/new"
  mkdir -p "$refd" "$newd"
  # strip any existing root=, then append our own + evolver setting
  grep -vE '^\s*root\s*=' "$ini" > "$WORK/$name/base.ini"
  { cat "$WORK/$name/base.ini"; printf "\nroot = %s/out_\n" "$refd"; }                       > "$WORK/$name/ref.ini"
  { cat "$WORK/$name/base.ini"; printf "\nroot = %s/out_\nevolver = 2\n" "$newd"; }          > "$WORK/$name/new.ini"
  if ! "$CLASS" "$WORK/$name/ref.ini" > "$WORK/$name/ref.log" 2>&1; then
    echo "RUNFAIL(ndf15) $name"; overall=1; continue
  fi
  if ! "$CLASS" "$WORK/$name/new.ini" > "$WORK/$name/new.log" 2>&1; then
    echo "RUNFAIL(rkdp45) $name (see $WORK/$name/new.log)"; overall=1; continue
  fi
  if ! ls "$refd"/out_*.dat >/dev/null 2>&1; then
    echo "NOOUT $name (scenario produces no .dat)"; continue
  fi
  echo "== $name =="
  if ! python3 "$HERE/compare_tol.py" "$refd" "$newd" 'out_*.dat'; then
    overall=1
  fi
done
echo "WORKDIR=$WORK"
[ "$overall" -eq 0 ] && echo "ALL SCENARIOS PASS" || echo "SOME SCENARIOS FAILED"
exit $overall
```

- [ ] **Step 2: Make executable & run**

```bash
chmod +x test/scenarios/compare_evolver.sh
make class -j4
test/scenarios/compare_evolver.sh ./class | tee /tmp/ev/scenarios.txt
```
Expected: a per-scenario block of `OK` lines. Note any `FAIL` or `RUNFAIL(rkdp45)` — these are the legitimate research findings (e.g. a stiff scenario the explicit solver can't meet at RTOL=1e-3). Do **not** mask them.

- [ ] **Step 3: Commit the harness**

```bash
git add test/scenarios/compare_evolver.sh
git commit -m "test(evolver): scenario comparison harness for rkdp45 vs ndf15"
```

---

### Task 7: Benchmark ndf15 vs rkdp45

**Files:**
- Uses: `class_profiled` (built from `make class_profiled`), `explanatory.ini`, one representative scenario.

- [ ] **Step 1: Build the profiler**

```bash
make class_profiled -j4
```
Expected: builds `class_profiled`.

- [ ] **Step 2: Profile ndf15 (default) on explanatory**

```bash
./class_profiled explanatory.ini 7 2>/dev/null | tail -15
```
Record the Perturbations + Background + TOTAL medians.

- [ ] **Step 3: Profile rkdp45 on explanatory**

```bash
{ cat explanatory.ini; printf "\nevolver = 2\n"; } > /tmp/ev/expl_rkdp45.ini
./class_profiled /tmp/ev/expl_rkdp45.ini 7 2>/dev/null | tail -15
```
Record the same medians.

- [ ] **Step 4: Repeat for one representative scenario (e.g. `test/scenarios/gauge_lcdm.ini` and `ncdm_single.ini`)**

```bash
for s in gauge_lcdm ncdm_single; do
  echo "=== $s : ndf15 ==="; ./class_profiled test/scenarios/$s.ini 5 2>/dev/null | tail -13
  { cat test/scenarios/$s.ini; printf "\nevolver = 2\n"; } > /tmp/ev/$s.ini
  echo "=== $s : rkdp45 ==="; ./class_profiled /tmp/ev/$s.ini 5 2>/dev/null | tail -13
done
```
Record results.

- [ ] **Step 5: No commit (data captured for Task 8).**

---

### Task 8: Record results + open PR

**Files:**
- Create: `docs/superpowers/specs/2026-06-07-evolver-rkdp45-results.md`

- [ ] **Step 1: Write the results doc**

Include: smoke-test outcome (Task 5); the full scenario pass/fail table from `/tmp/ev/scenarios.txt` (Task 6); the benchmark table ndf15 vs rkdp45 (Task 7) with per-module medians and the wall-time delta on the perturbations path; and a conclusion answering "does rkdp45 reproduce observables within RTOL=1e-3, and is it competitive on speed?" Note the step-controller bug fix vs the PhD2024-EBH branch (plan callout) as an intentional deviation.

- [ ] **Step 2: Commit**

```bash
git add docs/superpowers/specs/2026-06-07-evolver-rkdp45-results.md
git commit -m "docs(evolver): rkdp45 correctness + benchmark results"
```

- [ ] **Step 3: Run the existing classy test suite as a regression gate (default evolver unchanged)**

```bash
make classy
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python -m pytest -q -m test_scenario test_class.py; cd ..
```
Expected: pass (we did not change the default evolver; this confirms no collateral breakage).

- [ ] **Step 4: Push and open the PR**

```bash
git push -u origin feat/evolver-rkdp45
gh pr create --title "feat(evolver): port Dormand-Prince 4(5) solver (rkdp45), opt-in & verified" --body "<summary: what, opt-in via evolver=2, scenario results table, benchmark table, intentional step-controller bug fix vs PhD2024-EBH branch, scope boundaries (no Tsit5/DOP853/PI-PID yet)>"
```

---

## Self-review notes (addressed)

- **Spec coverage:** enum+wiring (T1,T2,T4) ✓; faithful standalone port (T3) ✓; RTOL=1e-3 scenario verification (T5,T6) ✓; class_profiled benchmark (T7) ✓; scope boundaries respected (no Tsit5/DOP853/controllers); default unchanged + regression gate (T8 S3) ✓.
- **Placeholders:** none — full source and exact commands inline. (Only the PR body text and results doc prose are authored at execution time from captured data, which is correct.)
- **Type/signature consistency:** ported `evolver_rkdp45` signature is byte-for-byte the canonical evolver signature; verified by the compile gate in T4 S5 (the `auto generic_evolver` assignment will not compile otherwise).
- **Intentional deviation from branch:** the dead step-rejection branch is corrected to standard ode45 (T3 code + callout), flagged in the PR (T8).
