# Issue #325 — Perturbations hot-path dispatch cache

Date: 2026-06-19
Issue: [#325](https://github.com/AarhusCosmology/CLASSpp/issues/325)

## Problem

After the species-layout refactor, the perturbations module regressed ~1.7 s
(`3966 ms` → `5633 ms` perturbations median) versus the old fork on the
Planck-2018 single-NCDM benchmark. Issue #325 isolates the cause to the way the
RHS hot path accesses per-species perturbation layouts. The dominant costs, all
paid **per RHS evaluation** (millions of calls):

1. **String lookups in the hot path.** Photon/baryon layouts are reached via
   `all_species_.index_of("Photons")` / `index_of("Baryons")` — a binary search
   over keyed `std::string` — repeated at ~15 hot sites
   (`perturb_einstein`, `perturb_total_stress_energy`, `perturb_derivs`,
   `perturb_sources`, `perturb_tca_slip_and_shear`, `perturb_rsa`).
2. **Two full passes over `all_species_` in the scalar RHS.** Pass 1 skips
   photons/baryons by `sp->name()` **string compare** and skips "deferred"
   species; pass 2 runs only the deferred species. Both manually carry a
   parallel index into `pv->species_layouts`.
3. **No-op calls.** Lambda (and, under RSA, photons) are visited every pass even
   though their `PerturbDerivs` does nothing.

A probe that cached photon/baryon perturbation indices recovered ~0.63 s of the
~1.7 s and produced **bitwise-identical** `cl`/`pk`/nonlinear output, confirming
this is pure dispatch overhead, not a physics change.

### Out of scope

The issue's later comments also document that most of the *public* `class_public`
speed advantage is a **workload/precision** difference (halofit `k_max`
handling; the New-Limber `k_scalar_max_tau0_over_l_max` 2.4→1.8 reduction), not
RHS implementation. Those are not addressed here. `HasNcdm` cleanup is also out
of scope — #325 measured it as not a source of this regression.

## Key findings (drive the design)

- **The "deferred second loop" is spurious.** The only deferred species is the
  PPF fluid, and its deferred `PerturbDerivs` does exactly one thing:
  `dy[idx_Gamma] = ppw->Gamma_prime_fld`. `Gamma_prime_fld` is published by
  `FluidSpecies::ComputePpf`, which already runs **earlier in the same RHS**
  (`perturb_derivs` → `perturb_einstein` (5850) → `perturb_total_stress_energy`
  → `ComputePpf` (5002)), well before the species loop (6061). So the fluid
  derivs can run in-loop at any position. The deferral mechanism can be deleted.
- **PPF is already module-owned.** The module caches `ppf_fluid_` and calls
  `ComputePpf` explicitly (PR #322). Nothing more is needed module-side.
- **`RequiresDeferredPerturbDerivs`** is overridden only by `FluidSpecies`
  (unconditionally `true`, even non-PPF fluid, whose derivs were always
  self-contained). It is deleted entirely.
- **One real ordering constraint.** The perturbed-recombination block
  (`perturbations_module.cpp:6042–6056`) reads `dy[baryon.idx_delta]`, so it must
  run *after* baryon derivs. Nothing reads its output, so it is safe to move it
  below the dispatch loop.
- **Lambda** is the canonical no-op: `lambda.h:44` is an inline empty
  `PerturbDerivs`. It registers zero perturbation indices in every mode.
- **A `perturb_vector` is single-mode.** `perturb_solve` builds a `pv` per
  `index_md`; `perturb_vector_init` populates `species_layouts` through exactly
  one of the `_scalars_` / `_vectors_` / `_tensors_` branches. So a `pv` never
  needs more than one mode's dispatch list. Vectors and tensors are **not** cold
  — tensors become as hot as scalars when B-mode output is requested — so all
  three modes receive identical treatment.

## Goals

- Remove per-RHS string lookups, name comparisons, the second pass, and no-op
  derivs calls from all three (scalar/vector/tensor) RHS dispatch loops.
- Keep `all_species_` as the canonical container everywhere else; the dispatch
  cache is a hot-path-only derived view rebuilt with each `pv`.
- Byte-identical output expected (same calls, same order, no reordered
  reductions); verified, with a ~0.1 % fallback bar.

## Design

### 1. `SpeciesCollection`: cache photon/baryon indices (issue suggestion #1)

`freeze()` already resolves the cached `photons_`/`baryons_` pointers. Add the
parallel indices and accessors:

```cpp
std::size_t photons_index_ = 0, baryons_index_ = 0;   // set in freeze()
std::size_t photons_index() const { assert(frozen_); return photons_index_; }
std::size_t baryons_index() const { assert(frozen_); return baryons_index_; }
```

These are used only to resolve the typed `pv` pointers once per `pv`-build (not
per RHS).

### 2. `perturb_vector`: per-mode dispatch list + typed photon/baryon pointers

```cpp
struct ActiveSpecies {
  BaseSpecies*                species;
  BaseSpecies::PerturbLayout* layout;
};

// Species with ≥1 registered perturbation variable in THIS pv's mode, in
// lex-key order, no-ops (e.g. Lambda; photons under RSA) excluded.
std::vector<ActiveSpecies> active_species;

// Always-present species; typed, resolved once per pv. Used by the ~15 hot
// sites that read photon/baryon layout fields directly.
PhotonsSpecies::PerturbLayout* photon_layout = nullptr;
BaryonsSpecies::PerturbLayout* baryon_layout = nullptr;
```

Resolved immediately after `species_layouts` is populated
(`perturbations_module.cpp:2939`), guarded by `assert` that both are non-null:

```cpp
ppv->photon_layout = static_cast<PhotonsSpecies::PerturbLayout*>(
    ppv->species_layouts[all_species_.photons_index()].get());
ppv->baryon_layout = static_cast<BaryonsSpecies::PerturbLayout*>(
    ppv->species_layouts[all_species_.baryons_index()].get());
```

**Ownership — non-owning views.** `species_layouts` remains the *sole owner*
(`vector<unique_ptr<PerturbLayout>>`, full, parallel to `all_species_` including
Lambda's empty slot). `active_species`, `photon_layout`, and `baryon_layout` are
**non-owning** views into it. Ownership must stay collection-parallel because
cold paths index `species_layouts` by collection position — the registration
loop, `Mark*UsedInSources`, `CopyPerturbationsAcrossSwitch` (all iterate every
`i`, incl. Lambda), and per-species self-reads via
`species_layouts[collection_index_]` (e.g. `Delta`/`Theta`/`FillSources`/
`ApplyInitialConditions`, `fluid.cpp:255`). Lifetime-safe: `species_layouts` is
`reserve`+filled once (2936–2939) *before* any `.get()` is taken, never resized
or reset afterward, and the heap `PerturbLayout` objects are address-stable;
the `pv` is rebuilt wholesale per approximation switch, so owner and views share
one lifetime.

**Header coupling — pointers are base-typed.** `photons.h` and `baryons.h` both
`#include "perturbations.h"` (line 4), so `perturbations.h` *cannot* include them
back (cycle), and a nested `PhotonsSpecies::PerturbLayout` cannot be
forward-declared. Therefore `photon_layout` / `baryon_layout` are stored as base
`BaseSpecies::PerturbLayout*` and `static_cast` to the concrete type at the read
sites (in `perturbations_module.cpp`, which already includes the concrete
headers). This is minimal churn: the read sites already write
`static_cast<PhotonsSpecies::PerturbLayout&>(*pv->species_layouts[g_i])`; only the
operand changes to `*pv->photon_layout`. The expensive part — the binary-search
string lookup — is what's removed.

### 3. Building `active_species` — data-driven, during registration

Each mode branch already threads one `index_pt` accumulator through
`entry->Register*PerturbationIndices(layout[i], …, index_pt, …)`. A species
advances `index_pt` by exactly the number of variables it registers, and its
`Perturb*Derivs` can only write `dy[layout.idx_*]` for slots it registered.
Therefore *index_pt advanced ⟺ the derivs is non-trivial in this mode* — an
exact test, not a heuristic, needing no per-species virtual flag:

```cpp
size_t i = 0;
for (auto& entry : all_species_) {
  const int before = index_pt;
  entry->RegisterPerturbationIndices(*ppv->species_layouts[i], …, index_pt, …);
  if (index_pt > before)
    ppv->active_species.push_back({entry.get(), ppv->species_layouts[i].get()});
  ++i;
}
```

Applied identically in the `_scalars_`, `_vectors_`, and `_tensors_` branches
(with `RegisterVector*` / `RegisterTensor*`). Properties this buys for free:
Lambda is dropped in every mode; a species that perturbs scalars but not tensors
is in the scalar list and absent from the tensor list; and because the list is
rebuilt with each `pv` (each approximation switch), RSA-disabled photons drop
out of the scalar list automatically — eliminating the issue's "empty calls".

### 4. Unified scalar RHS loop

Replace the explicit baryon/photon calls **and** both passes
(`perturbations_module.cpp:6022–6075`) with a single loop, and move the
perturbed-recombination block to *after* it (it reads `dy[baryon.idx_delta]`):

```cpp
// scalar_ctx populated; perturb_tca_slip_and_shear already ran above.
for (const auto& [species, layout] : pv->active_species)
  species->PerturbDerivs(*layout, tau, y, dy, *pppaw);

// perturbed recombination (moved down: reads dy[baryon.idx_delta])
if (ppt->has_perturbed_recombination && tca_off) { … }
```

Loop order is free: every species reads `scalar_ctx`/`y` and writes only its own
`dy`. Lex order keeps baryons before photons, matching today.

### 5. Vector and tensor RHS loops

Replace the all-species loops at `6119` and `6138` with the same form over
`pv->active_species` (calling `PerturbVectorDerivs` / `PerturbTensorDerivs`).
The typed `photon_layout` / `baryon_layout` pointers serve the vector/tensor hot
sites (`g_vec`, `g_tens`, `b_vec`) too.

### 6. Scattered hot-site replacement

Convert the ~15 **hot** `index_of("Photons")` / `index_of("Baryons")` read sites
to `pv->photon_layout` / `pv->baryon_layout` (in `perturb_einstein`,
`perturb_total_stress_energy`, `perturb_derivs` setup, `perturb_sources`,
`perturb_tca_slip_and_shear`, `perturb_rsa`, and the vector RHS). Leave the
**cold** per-approximation-switch sites in the `perturb_vector` construction
region (`~3209–3527`, `CopyPerturbationsAcrossSwitch`) on `index_of` — they
juggle old and new `pv`s and are not per-RHS.

### 7. Deletions

- `BaseSpecies::RequiresDeferredPerturbDerivs` (`base_species.h:260`).
- `FluidSpecies::RequiresDeferredPerturbDerivs` override (`fluid.h:83`).
- The deferred second pass and the `sp->name()` string filtering in the scalar
  loop.

## Correctness argument (why byte-identical)

- Same set of `PerturbDerivs` calls, in the same lex order (baryons before
  photons), with the same arguments. The explicit photon/baryon calls become
  ordinary list entries at the same relative position.
- Fluid moves from the second pass to its lex position in the single loop, but
  its `dy[idx_Gamma] = Gamma_prime_fld` depends only on `Gamma_prime_fld`
  (already published) — value-independent of loop position.
- The recombination block moves below the loop but reads the same
  `dy[baryon.idx_delta]` (written by baryon, now in the loop) — identical value;
  its outputs feed nothing in the RHS.
- No reductions are reordered (the dispatch loop performs independent `dy`
  writes, no accumulation), so the NCDM ULP-drift risk from
  [[feedback_vectorization_reduction_drift]] does not apply.

## Verification

1. Build clean (CMake) with asserts (Debug) and Release.
2. A few bit-exact scenarios (incl. one with PPF fluid, one with massive NCDM,
   one with tensors/B-modes) — expect byte-identical `cl`/lensed `cl`/`pk`.
3. Full `TEST_LEVEL=2` against regenerated `classyref`.
4. Re-profile with `class_profiled` on the issue's Planck-2018 single-NCDM
   benchmark (single-threaded) and report the perturbations median vs. the
   `5633 ms` baseline and the `5075 ms` cached-index experiment.

## Risks

- **Include cycle** from typed pointers in `perturbations.h` — mitigated by the
  base-pointer + cast fallback (§2).
- **A species whose derivs is non-trivial but registers no indices** would be
  wrongly filtered — but this cannot happen: derivs writes only registered
  slots. Asserted indirectly by byte-identical output.
- **PPF/fluid edge cases** — covered by a PPF-fluid bit-exact scenario.
