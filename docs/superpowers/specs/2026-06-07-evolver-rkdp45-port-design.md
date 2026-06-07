# Design: Port the Dormand-Prince 4(5) evolver (`rkdp45`) to master

**Date:** 2026-06-07
**Status:** Approved (design phase)
**Scope:** Port the already-implemented DP45 ODE solver from branch `PhD2024-EBH`
to current master as an opt-in evolver, and rigorously verify that observables
agree with the `ndf15` reference within the established tolerance bar.

## Background & motivation

CLASS++ integrates the background and perturbation ODE systems through a
`generic_evolver` function pointer selected by `enum evolver_type` (in
`include/common.h`), currently `{ rk, ndf15 }` with `ndf15` (stiff, implicit)
the default. `rk` (`evolver_rkck.cpp`) is a thin wrapper over the legacy
Cash-Karp `generic_integrator`.

A Dormand-Prince 4(5) explicit adaptive solver was implemented on the research
branch `PhD2024-EBH` as `tools/evolver_rkdp45.cpp`, but never merged because it
was never rigorously verified to reproduce observables within tolerance. That
branch has diverged heavily from master (different `.c`/`.cpp` file layout,
predates several refactors), so a merge is not viable — the solver must be
**ported** to current master.

This is the first increment of a larger programme (later: Tsit5, DOP853,
PI/PID step controllers). This increment deliberately lands **one** verified
solver and the evidence needed to decide whether the family is worth building
on.

## Decisions (locked during brainstorming)

1. **Structure:** Standalone file, minimal. Port DP45 as
   `tools/evolver_rkdp45.cpp` mirroring `evolver_rkck.cpp`, with a hardcoded
   Butcher tableau and a fixed I-controller. No generic-RK-core refactor in
   this increment.
2. **Acceptance bar:** The established `test/scenarios/compare_tol.py` bar
   (zero-crossing-aware, `RTOL=1e-3`) that master's CI already enforces.
   `rkdp45` output must pass it against the `ndf15` reference across the
   scenario set.

## Architecture & components

### New solver
- `tools/evolver_rkdp45.cpp` + `include/evolver_rkdp45.h`.
- `int evolver_rkdp45(...)` conforms **exactly** to master's canonical evolver
  signature — the one shared by `evolver_ndf15` and `evolver_rk` — so it slots
  into the existing `generic_evolver` function-pointer dispatch with no
  signature changes.
- The branch's extra diagnostic parameter `std::vector<double>& t_vec_evolver`
  (only used to `push_back(t)` recording step locations) is **dropped**.

### Selection / wiring
- Extend the enum: `enum evolver_type { rk, ndf15, rkdp45 }` in
  `include/common.h`. Default stays `ndf15`.
- Map the string `"rkdp45"` in the `read_enum` path used by
  `input_module.cpp` (`read_enum(fc, "evolver", evolver)`), so
  `evolver = rkdp45` selects it from an `.ini`.
- Add a third dispatch branch in **both** `source/background_module.cpp` and
  `source/perturbations_module.cpp`:
  `else if (ppr->evolver == rkdp45) generic_evolver = &evolver_rkdp45;`

### Build
- Add `evolver_rkdp45.opp` to `TOOLS` in the `Makefile`.
- Add the source to any other build manifests that enumerate tool sources
  (CMake / `setup.py` / Xcode project) **only if** they list evolvers
  explicitly. Verify during implementation; do not assume.

## The port itself

Faithful 1:1 port of the branch algorithm — **no algorithm changes**:
- Same DP45 Butcher tableau (`ci`, `bi`, `ai`).
- Same embedded 4th/5th-order error estimate (`bi_diff`).
- Same 4th-order dense-output interpolant (`ixx`, `i01..i03`) used to emit
  values at the requested `x_sampling` output points.
- Same fixed I-controller step logic: `pow_grow = 0.2`, 0.8 safety factor,
  `abstol = 1e-15`, `threshold = abstol/rtol`, step grow/shrink and
  reject-and-retry on error.

Modernized lightly to match master conventions:
- `malloc`/`free` → `std::vector<double>` (RAII, as in `evolver_rkck.cpp`).
- `_TRUE_`/`_FALSE_`, `class_call`/`class_test` error handling.

## Verification (the reason it never merged)

1. Build master `class` (unchanged) and generate **ndf15 reference** outputs by
   running the `test/scenarios/*.ini` set.
2. Re-run the same scenarios with `evolver = rkdp45` and compare via
   `test/scenarios/compare_tol.py` at `RTOL=1e-3` (zero-crossing-aware) — the
   exact bar CI enforces. **Pass = every scenario within tolerance.**
3. The scenario set spans diverse physics already present in
   `test/scenarios/` (lcdm, ncdm, idmdr, dcdm, scalar field, curved, …),
   including stiff regimes — the real test for a non-stiff explicit solver.

## Benchmark

- Compare **ndf15 vs rkdp45** per-module wall time using the `class_profiled`
  harness (Perturbations and Background are the affected stages), on
  `explanatory.ini` plus 1–2 representative scenarios, median of N loops.
- **Dependency (resolved):** `class_profiled` (PR #299) is merged into master
  and this branch is rebased onto it, so the harness is available here.
- Report the wall-time delta. This tells us whether DP45 is competitive on the
  ~57%-of-runtime perturbations path or is purely a correctness exercise.

## Expected-outcome note

DP45 is a **non-stiff** explicit method, while the perturbation system is stiff
during tight coupling. It is plausible that DP45 either (a) requires many tiny
steps (slow) or (b) fails some stiff scenarios at `RTOL=1e-3`. **Either outcome
is a legitimate, useful result** — part of what this increment measures is
whether the explicit-RK family is viable at all on these systems. The design
does not assume success.

## Scope boundaries (explicitly NOT in this increment)

- No Tsit5 / DOP853 solvers.
- No PI/PID step controllers.
- No generic-RK-core refactor (tableau-parameterized stepper).
- No change to the default evolver — `rkdp45` is opt-in.

These constitute the follow-on programme, each its own spec → plan → build
cycle, informed by the correctness/benchmark evidence from this increment.
