# Perturbations RHS: collapse the four stress-energy virtuals into a single `StressEnergy` struct method

**Date:** 2026-06-21
**Status:** Design approved, ready for implementation plan
**Area:** `species/` (BaseSpecies + all concrete species) and `source/perturbations_module.cpp`

## Motivation

The perturbations RHS (`perturb_derivs` → `perturb_einstein` → `perturb_total_stress_energy`)
is the only reducible part of the C++/C performance gap. After PR #328, perturbations is
~1.09× fair-flags C (`-O3`, same Apple-clang, no `-ffast-math`); the evolver (~76% of the
stage) is shared with C and equal, so the entire reducible gap lives in the **~24% RHS** —
non-inlinable virtual dispatch per species per ODE step.

A session of paired/rotated benchmarking (2026-06-21, recorded in the
`project-perturbations-regression-vs-c` memory) established:

- The hottest perturbations-*module* leaf after the solver was the naive base `StressEnergy`
  wrapper (533 `sample` hits) — because a base wrapper that delegates to N inner virtuals
  makes cheap species pay *more* dispatches, not fewer.
- A **direct per-species `StressEnergy` returning a small struct** is the keeper: bit-identical,
  and it moved the perturbations stage from ~6% above master to **~3.4% above fair-flags C**
  (master → delegate ≈ −6%, paired mean/se 9.18, 18/18 rounds).
- A **by-reference accumulator** variant (each species adds only its nonzero terms into a
  passed-in accumulator) was measured to give **no gain** (paired delegate→byref mean/se −1.07,
  statistically zero) and it **broke bit-identity** (ULP drift). At `-O3` the 6-double struct
  return is already RVO'd into registers, so there is no round-trip to eliminate. **Rejected.**

This spec captures the keeper design, refined so that the four existing per-species stress-energy
virtuals are **removed entirely** and replaced by the one struct-returning method — a net
reduction in virtual surface area with no duplicated code.

## Goals

- Replace the four pure virtuals `DeltaRho`, `RhoPlusPTheta`, `DeltaP`, `RhoPlusPShear` on
  `BaseSpecies` with a single pure virtual `StressEnergy(...) -> StressEnergyContribution`.
- Each species computes all its scalar stress-energy perturbations in **one** dispatch.
- NCDM computes them in **one fused pass over the q-grid** (`epsilon = √(q² + (M·a)²)` once),
  reused by both the hot loop and NCDM's own output paths.
- Fold the anisotropic-stress (`rho_plus_p_shear`) pre-pass into the single main loop.
- Stay within **~0.1% physics agreement** (`Cl^TT`, `P(k)`) of current master — NOT bit-identical.
  This is a perf refactor: prefer clean, fast code over byte-for-byte reproduction. Most of the
  change happens to come out bit-identical, but NCDM's single fused q-pass deliberately reorders
  the per-bin reductions and drifts at the ULP level (~1e-5 on `Cl^TT`, well under the bar) — that
  is accepted, not worked around. Baselines/`classyref` are regenerated at the end.

## Non-goals (explicitly out of scope)

- The `all_species_.count("IDM_DR_IDR")` per-step string-scan caching in `perturb_einstein`
  (the "einstein-bool" + count-harvest, ~3%, separately bit-identical). **Separate follow-up PR.**
- Any change to the ODE evolver, background, or the vector/tensor stress-energy blocks.
- Devirtualizing `Rho`/`P` or the per-species `PerturbDerivs`.
- The by-reference accumulator (measured to give nothing — see Motivation).
- A real physics fix for scalar-field-in-Newtonian-gauge (issue #285); this design only
  preserves the existing hard `class_test` error for that path.

## Design

### The struct and the single virtual

In `species/base_species.h`, **remove** the four pure virtuals `DeltaRho`, `RhoPlusPTheta`,
`DeltaP`, `RhoPlusPShear`. Add:

```cpp
/** The scalar stress-energy contribution of one species at one ODE step
 *  (Ma & Bertschinger). Returned by value: six doubles, RVO'd into registers
 *  at -O3 (the by-reference accumulator variant was measured to give no gain). */
struct StressEnergyContribution {
  double rho              = 0.;  // background ρ
  double p                = 0.;  // background P
  double delta_rho        = 0.;  // δρ
  double rho_plus_p_theta = 0.;  // (ρ+P)θ
  double delta_p          = 0.;  // δp
  double rho_plus_p_shear = 0.;  // (ρ+P)σ

  // Element-wise accumulation, used by composite species to sum children.
  StressEnergyContribution& operator+=(const StressEnergyContribution& o) {
    rho += o.rho; p += o.p; delta_rho += o.delta_rho;
    rho_plus_p_theta += o.rho_plus_p_theta; delta_p += o.delta_p;
    rho_plus_p_shear += o.rho_plus_p_shear;
    return *this;
  }
};

/** All scalar stress-energy perturbations in one call. Pure virtual: every
 *  concrete species computes its quantities directly (and, for NCDM, in a
 *  single fused pass over the q grid). */
virtual StressEnergyContribution StressEnergy(const PerturbLayout& layout,
                                              const perturb_vector* pv,
                                              const double* y,
                                              const double* pvecback,
                                              const perturb_workspace* ppw) const = 0;
```

`StressEnergy` is **pure virtual — no base default** (user-confirmed). Every concrete species
already implements all four removed virtuals, so the conversion is a uniform 4-bodies → 1-body
mechanical change with no species left behind. A base default that delegated to inner virtuals
is exactly the wrapper that the profile flagged as the 2nd-hottest leaf, so we deliberately
avoid it.

`Rho` and `P` remain their own virtuals — they are used throughout the code for background
quantities (vector/tensor sources, `perturb_einstein`'s idr terms, etc.). The struct merely
carries copies so the hot loop need not separately dispatch them.

### Per-species implementations

| Species | Implementation |
|---|---|
| `PhotonsSpecies`, `BaryonsSpecies`, `CDMSpecies`, `UltraRelativisticSpecies` | Fold the four existing method bodies into one `StressEnergy`; compute `Rho`/`P` inline once. Cold species naturally return zeros for the fields they don't source. |
| NCDM (placed at the **`NCDMBaseSpecies`** level so all ncdm variants share it) | **Single fused pass over the q grid**: compute `epsilon = √(q² + (M·a)²)` and the `(a₀/a)⁴` power once; accumulate `delta_rho`, `rho_plus_p_theta`, `delta_p`, `rho_plus_p_shear` in one loop. Keep each per-term expression and its left-associative grouping identical to the four current methods, and handle the `ncdmfa_on` branch as today. This is the **sole** q-pass implementation; it also feeds NCDM's own output paths. |
| Composites (`DCDM_DR`, `DNCDM_DR`, `IDM_DR_IDR`, `IDM_DRMD`) | `return child_a.StressEnergy(...) += child_b.StressEnergy(...)` — i.e. start from the first child's struct (`child_a`, matching today's left operand) and `+=` the second child's. Replaces four delegating methods with one. |
| `FluidSpecies` | One `StressEnergy`; `delta_rho`/`rho_plus_p_theta`/`delta_p` as in the current methods, `rho_plus_p_shear = 0`. PPF self-zeroing still falls out of the existing layout-index guards (`idx_delta >= 0 ? … : 0.`). |
| Any other concrete species with the four methods (Λ, scalar field, etc.) | Same mechanical 4 → 1 conversion. |

### Loop: `perturb_total_stress_energy` (scalar block, ~lines 4811-4885)

- **Delete the shear pre-pass** (the `for (... ) ppw->rho_plus_p_shear += ...RhoPlusPShear(...)`
  loop, current lines ~4838-4839).
- In the single species loop, call `e.species->StressEnergy(...)` **once** and read the six
  fields. Accumulate `delta_rho`, `rho_plus_p_theta`, `delta_p`, `rho_plus_p_tot (= rho + p)`,
  **and `rho_plus_p_shear`** (folded in here now). Accumulation order over `active_species` is
  unchanged, so the totals are bit-identical.
- The `(ρ−3P)` cold/warm matter tally stays **loop-local**, reading `rho`/`p`/`delta_rho`/
  `delta_p`/`rho_plus_p_theta` straight off the returned struct (unchanged arithmetic). **No
  accumulator object** — that was the rejected by-reference design.
- Add a code comment at the deleted pre-pass site documenting that a future issue-#285
  Newtonian-scalar-field fix must restore "complete shear before the scalar field's
  `DeltaRho`/`DeltaP`" — either by reinstating a pre-pass or by ordering the scalar field last.
- Everything after the loop is **unchanged**: the `delta_cb`/`theta_cb`/`delta_m`/`theta_m`
  derivations and the post-loop fluid `ComputePpf` block (which adds `delta_rho_fld` etc. and
  `f->Rho() + f->P()`) stay exactly as today.

### Other callers (all cold paths)

- **IC path** `source/perturbations_module.cpp:3984-3985` (downcast-free IC analogue, runs once
  per k): replace the separate `DeltaRho` + `RhoPlusPTheta` calls with one `StressEnergy(...)`,
  reading `.delta_rho` and `.rho_plus_p_theta`.
- **NCDM output** `species/ncdm_species.cpp` — `FillSources` (~285/297, needs δρ and (ρ+P)θ) and
  `PrintVariables` (~479-482, needs **all four**): one `StressEnergy(...)` call each, read the
  needed fields. Collapses 2-4 separate q-passes into one in each output function.
- **Fluid output** `species/fluid.cpp:256-258` (`PrintVariables`, non-PPF branch): one
  `StressEnergy(...)` call, read `.delta_rho`/`.rho_plus_p_theta`/`.delta_p`.

Documentation comments that *mention* the old method names (`base_species.h:70`,
`composite_species.h:21`, `idr.h:122`, `ultra_relativistic.cpp:73`, `photons.cpp:164`) are
updated to refer to `StressEnergy`. These are comments, not call sites.

## Why this is bit-identical

- **Hot loop:** the struct delegation was measured byte-identical (`cl_lensed.dat` + `pk.dat`)
  on the 1-massive-ν benchmark; field accumulation order over `active_species` is preserved.
- **Shear fold:** the only mid-loop reader of the running `rho_plus_p_shear` is scalar-field-in-
  Newtonian-gauge, which is a hard `class_test` error at perturbation init (the #285 path), so
  no runnable config observes a partial sum. Post-loop readers (`perturb_einstein` psi/α′,
  fluid `ComputePpf`) see the same complete sum. Previously verified byte-identical.
- **NCDM fused pass:** per-term expressions and grouping are preserved, so each accumulator is
  summed in the same order; `epsilon` and the `(a₀/a)⁴` power are deterministic. The prototype
  was verified byte-identical to the four separate methods.
- **Composites via `+=`:** for fields where both children currently contribute, `a.field +
  b.field` matches today's `child_a->Method + child_b->Method` order. For the fields the current
  code selectively omits (a cold child's `delta_p`/`rho_plus_p_shear`), that child's struct field
  is exactly `0.0`, and `0.0 + x == x`. Confirm during implementation that each cold child's
  `StressEnergy` returns literal `0.0` for the omitted fields (the existing code comments —
  "IDM_DR has DeltaP == 0", "both have RhoPlusPShear == 0" — assert this by physics).

## Testing

1. **Physics agreement (primary bar, ~0.1%):** run against the clean-master reference and require
   `Cl^TT` and `P(k)` within `1e-3` relative (`.superpowers/sdd/tol_check.py`) on:
   - the 1-massive-ν benchmark (`base_2018_plikHM_TTTEEE_lowl_lowE_lensing`-style), and
   - at least one **composite** scenario (DNCDM_DR or IDM_DR) to exercise the struct-`+=` path,
   - optionally a multi-ncdm run to exercise the shared `NCDMBaseSpecies` pass.
   Most paths land bit-identical; NCDM drifts at ~1e-5 (ULP). Do not chase byte-identity.
2. Full `TEST_LEVEL=2` suite green (with `COMPARE_OUTPUT_REF=1` against `classyref`). Since NCDM's
   fused pass intentionally shifts output at the ULP level, **regenerate `classyref` / golden refs
   at the end** as the deliberate baseline update for this change.
3. **Performance:** re-run the paired/rotated single-thread Perturbations benchmark
   (`OMP_NUM_THREADS=1`, `-O3` no-ffast-math) vs clean master; expect ≈ −6% on the stage,
   landing ~3.4% above fair-flags C-O3. Confirm no regression vs the prototype numbers.

## Risks / notes

- **#285 coupling:** folding the shear pre-pass is safe only while scalar-field-Newtonian stays
  a hard error. The code comment at the deleted pre-pass site records the requirement for any
  future fix. (See the `scf-newtonian-missing-physics` memory / issue #285.)
- **Diff size:** larger and more invasive than an additive prototype (touches every species and
  the composites), but it *removes* net code and surface area. The change is mechanical and each
  species is independently verifiable.
- **`Rho`/`P` retained:** intentionally not collapsed into the struct as the only source, because
  they have many background callers outside stress-energy.
- **Bench tooling:** `build_prof/` CMake build of `./class_profiled`; fair-flags C harness in the
  `809e706` worktree; paired-rotated harness from the prior session (bash-3.2-safe).
