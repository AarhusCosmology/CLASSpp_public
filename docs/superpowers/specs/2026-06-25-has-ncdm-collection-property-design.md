# Design: `has_ncdm` as a cached property of `SpeciesCollection`

Date: 2026-06-25

## Goal

Replace the scattered, O(n) "does the model contain any NCDM-family species?"
checks with a single boolean cached once on `SpeciesCollection` and exposed as
`all_species_.has_ncdm()`.

## Motivation

Today the NCDM-presence predicate is recomputed all over `perturbations_module`:

- `source/perturbations_module.cpp:53` — `HasNcdm(all_species_)` free function
  (anonymous namespace) that loops over `all_species_` doing
  `dynamic_cast<NCDMBaseSpecies*> || dynamic_cast<DNCDM_DR_Species*>` per species.
  Called from 19 sites (per-k, per-interval, and setup paths), all boolean
  presence (`> 0`) checks.

`source/output_module.cpp:668` uses the **same per-species `dynamic_cast` pair**,
but as a `count` feeding `class_test(n_ncdm_family > 1, ...)` (CAMB format permits
at most one ncdm species). That is a `> 1` check, not presence, so a boolean
`has_ncdm()` does not apply there — it is out of scope (see below).

`SpeciesCollection` is frozen after construction (`freeze()`), and already caches
derived hot-path state there (`photons_`/`baryons_` pointers and indices). The
NCDM-presence flag is the same kind of frozen-once derived state and belongs in
the same place. Since the collection is immutable for the whole run, the cached
bool is always valid post-`freeze()`.

A previous performance pass added `ResolvedSpecies::has_ncdm` (a per-module copy
of the same bool) so the per-ODE-step hot path could avoid the scan. Once the
collection owns the flag, that module-level copy is redundant.

## Design

### 1. `SpeciesCollection` owns the flag (`species/species_collection.{h,cpp}`)

- Add a private member `bool has_ncdm_ = false;`.
- Add a hot-path accessor alongside `photons()`/`baryons()`:
  ```cpp
  bool has_ncdm() const { assert(frozen_); return has_ncdm_; }
  ```
- In `freeze()`, compute it inside the existing index-resolving loop (no second
  pass):
  ```cpp
  if (dynamic_cast<NCDMBaseSpecies*>(species_[i].get()) ||
      dynamic_cast<DNCDM_DR_Species*>(species_[i].get()))
    has_ncdm_ = true;
  ```
- `species_collection.cpp` gains `#include "ncdm_base_species.h"` and
  `#include "dncdm_dr_species.h"`.

This couples the otherwise-generic container to two concrete species types. That
is the accepted tradeoff of keeping the diff minimal (vs. a virtual predicate on
`BaseSpecies`); it is consistent with existing repo practice that sanctions
`dynamic_cast` for NCDM presence-detection.

### 2. Caller cleanup

- `source/perturbations_module.cpp`:
  - Delete the `HasNcdm` free function and its enclosing anonymous namespace
    (it is the only member).
  - Replace every `HasNcdm(all_species_)` call (19 sites) with
    `all_species_.has_ncdm()`.
  - Leave the existing ncdm `#include`s — `NCDMBaseSpecies`/`DNCDM_DR_Species`
    are still used elsewhere in the file (e.g. lines 575, 5006, 5627).

### 3. Remove the now-redundant `ResolvedSpecies::has_ncdm`

- `source/perturbations_module.h`: drop `bool has_ncdm` from the
  `ResolvedSpecies` struct.
- `source/perturbations_module.cpp`: drop the assignment in `ResolveSpecies()`
  and repoint the three per-step readers (`perturb_approximations`,
  `perturb_timescale_member`, `perturb_total_stress_energy`) at
  `all_species_.has_ncdm()`. The accessor returns a cached bool and is
  inlinable, so the hot path keeps O(1) reads.

## Out of scope

- `source/output_module.cpp:668` counts ncdm-family species and `class_test`s on
  `n_ncdm_family > 1` (CAMB-format single-ncdm constraint). That is a `> 1`
  check, not presence; a boolean `has_ncdm()` cannot express it, so it is left
  unchanged. (A future `ncdm_count()` could fold it in, but that is not in scope.)
- `source/background_module.cpp` `GetNcdmSpecies()` returns the actual
  `NCDMBaseSpecies*` list (only `NCDMBaseSpecies`, **not** the `DNCDM_DR_Species`
  composite). It answers a different question and is left unchanged.
- No `ncdm_count()` / generalization to other species — `has_ncdm()` (the `> 0`
  check) is all that's requested (YAGNI).

## Verification

The returned boolean is identical to what `HasNcdm(all_species_)` produced today,
so this is a pure dedup with no intended change in output.

- Build cleanly (C++) into `build_prof/` (`./class`).
- Capture a master baseline before editing, then A/B with
  `test/scenarios/compare_tol.py` (RTOL=1e-3, zero-crossing-aware): every
  output file must report `OK`. Scenarios: `base_2018...` (`N_ncdm=1`, the
  `true`/`NCDMBaseSpecies` path) and `gauge_lcdm` (no ncdm, the `false` path).
  Byte-identity is **not** required — under ffast-math (#338) ULP-level drift is
  expected and acceptable; only a meaningfully-off column fails
  (`feedback_no_bit_identical_requirement`).
- Final gate: `TEST_LEVEL=2 COMPARE_OUTPUT_REF=1` full suite (covers the
  DNCDM_DR-composite perturbation path, both gauges). No baseline regeneration
  expected.

## Risks

- Low. The only behavioral risk is an incorrect `freeze()` predicate; it mirrors
  the existing one exactly. Touching the `ResolvedSpecies` struct is mechanical.
