# Clean up `perturb_vector_init` — design

**Date:** 2026-05-30
**Predecessors:**
- `2026-04-13-perturbations-species-dispatch-design.md` (started the layout migration)
- `2026-05-07-encapsulate-species-layouts-design.md` (introduced `ppv->species_layouts`)

## Goal

Eliminate three accidental complexities in `PerturbationsModule::perturb_vector_init`
(`source/perturbations_module.cpp:3269-3747`):

1. Precision-consistency tests that live on the hot path and that are gated on
   species presence with hand-rolled scans (including a `dynamic_cast` to detect
   any DR-emitting species).
2. A duplicated registration call (`sp->RegisterPerturbationIndices(layout, …)`
   followed by `sp->RegisterPerturbationIndices(ppv, …)`) where the second call
   is now a stub on every species and exists only as a leftover from the layout
   migration.
3. A long sequence of `if (all_species_.count("X"))` blocks — one per species —
   that each grab `index_of("X")` and call exactly the same method. Plus a
   second cluster of the same shape inside the "used in sources" block.

After this PR, the registration of scalar / vector / tensor perturbation slots
is a single dispatch loop over `all_species_`; the "used in sources" masking is
a single dispatch loop; and all precision-consistency tests live once in the
input module.

## Scope

Single PR. Files touched:

| Area | Files |
|---|---|
| Precision tests upstream | `source/input_module.cpp` |
| Dead overload removed | `species/base_species.h`, `species/composite_species.{h,cpp}`, every species `.h` that declares `RegisterPerturbationIndices(perturb_vector*, …)` (≈20 files) |
| DNCDM_DR composite override | `species/dncdm_dr_species.{h,cpp}` |
| Species register methods read their own ppr | `species/ultra_relativistic.{h,cpp}`, `species/ncdm_species.{h,cpp}` |
| Baryons vector mode | `species/baryons.{h,cpp}` |
| `MarkUsedInSources` signature | `species/base_species.h`, `species/photons.{h,cpp}`, `species/ultra_relativistic.{h,cpp}`, `species/idm_dr_idr_species.{h,cpp}`, `species/ncdm_species.{h,cpp}` |
| Call sites | `source/perturbations_module.cpp` (`perturb_vector_init`, lines 3269-3747 region) |

Out of scope:

- `perturbed_recombination`, `tensor_ur_layout`, and the metric slots
  (`eta`, `phi`, `hv_prime`, `V`, `gw`, `gwdot`) stay as module-side blocks.
  They are not species; promoting them to a species API just to be uniform is
  over-abstraction.
- The rest of `perturbations_module.cpp` (initial conditions, sources,
  approximation transitions, …) is not touched.
- The legacy `pv->index_pt_<species>_*` member fields stay on
  `perturb_vector`. Removing them is a separate cleanup — many derivs/source
  call sites still read them through `pv->`.

## A. Precision tests → `input_read_parameters` end

The block of `class_test(ppr->l_max_X < 4, …)` calls scattered through
`perturb_vector_init` (lines 3294-3336, 3492-3500, 3533-3541, 3571-3576) moves
to a single block at the end of `InputModule::input_read_parameters`. Each
test runs **unconditionally** — no `all_species_.count(...)` gating, no
`dynamic_cast<DNCDM_DR_Species*>` scan to detect DR-emitting species. A
violated minimum is a configuration error even when the species would not
have been instantiated; flagging it always is simpler and more user-friendly
than silently letting the bad value sit until somebody adds the species.

Tests moved:

```cpp
class_test(ppr->l_max_g       < 4, …);  // photon temperature
class_test(ppr->l_max_pol_g   < 4, …);  // photon polarisation
class_test(ppr->l_max_ur      < 4, …);  // ultra-relativistic
class_test(ppr->l_max_dr      < 4, …);  // decay radiation (DCDM_DR + DNCDM_DR share)
class_test((ppr->l_max_idr < 4) && (ppt->idr_nature == idr_free_streaming), …);
class_test(ppr->l_max_g_ten     < 4, …);
class_test(ppr->l_max_pol_g_ten < 4, …);
class_test(ppr->l_max_ncdm      < 4, …);
```

Location: end of `input_read_parameters`, because the IDR test reads
`ppt->idr_nature` which is parsed there (not in `input_read_precisions`). All
other tests could live earlier in `input_read_precisions`, but splitting them
into two locations buys nothing.

## B. Single dispatch loop for registration

### B.1 Delete the legacy `RegisterPerturbationIndices(perturb_vector*, …)` overload

This overload (`species/base_species.h:227-231`) is now an empty `{}` body on
every species (e.g. `species/photons.h:69-73`, `species/cdm.h:49-53`). The
only non-stub implementation is `CompositeSpecies` (`species/composite_species.cpp:13-20`),
which just iterates its children and calls each child's empty stub. The
overload is pure transitional fossil from the layout migration.

Deletions:

- `BaseSpecies::RegisterPerturbationIndices(perturb_vector*, …)` — drop the
  pure-virtual declaration. Same for `RegisterVectorPerturbationIndices(perturb_vector*, …)`
  and `RegisterTensorPerturbationIndices(perturb_vector*, …)`.
- All ≈20 species-level stub overrides of these three methods.
- `CompositeSpecies::RegisterPerturbationIndices` and its vector/tensor twins.
- All `sp->RegisterPerturbationIndices(ppv, ppr, index_pt, ppw, ppt->gauge);`
  call sites in `perturb_vector_init` (one per species block, plus the inline
  NCDM dispatch loop at line 3465).

### B.2 `DNCDM_DR_Species` overrides `RegisterPerturbationIndices`

Today `DNCDM_DR_Species` deliberately does **not** override
`RegisterPerturbationIndices` (`species/dncdm_dr_species.h:72` comment). The
module instead pulls the composite apart and registers DR children with the
photon-side block (`perturbations_module.cpp:3396-3411`) and DNCDM children
with the NCDM-side block (`perturbations_module.cpp:3467-3473`). This is
what blocks a single dispatch loop.

After this PR, `DNCDM_DR_Species` overrides `RegisterPerturbationIndices` the
same way `DCDM_DR_Species` already does
(`species/dcdm_dr_species.cpp:68-77`):

```cpp
void DNCDM_DR_Species::RegisterPerturbationIndices(BaseSpecies::PerturbLayout& base,
                                                   perturb_vector* pv,
                                                   const precision* ppr,
                                                   int& index_pt,
                                                   const perturb_workspace* ppw,
                                                   int gauge) {
  auto& my = static_cast<PerturbLayout&>(base);
  dr_->RegisterPerturbationIndices(my.dr, pv, ppr, index_pt, ppw, gauge);     // DR child first
  dncdm_->RegisterPerturbationIndices(my.dncdm, pv, ppr, index_pt, ppw, gauge); // DNCDM child second
}
```

**Order chosen:** DR child first, DNCDM child second — matches the legacy
sequence (DR registered with photon-side, DNCDM with NCDM-side).

### B.3 Y-vector layout shift

The single-loop dispatch reorders the integration vector for setups with one
or more DNCDM_DR composites. Today:

```
… photons, baryons, …, [DR of all DNCDM_DR composites], …, [DNCDM of all DNCDM_DR composites], …
```

After:

```
… photons, baryons, …, [DR of composite 0, DNCDM of composite 0], [DR of composite 1, DNCDM of composite 1], …
```

Standalone NCDM species (no DR child) and all other species keep their
relative order. The `perturbed_recombination` insertion point also shifts
from "between ScalarField and UR" to "after the species loop, before metric
slots".

This is **not bit-identical** with master. Verification target: ~0.1% on Cl
spectra across a representative `.ini` suite that exercises DNCDM_DR (per
established project tolerance; see `feedback_no_bit_identical_requirement`).

### B.4 Per-species layout self-setup

Two pre-loop assignments in `perturb_vector_init` move inside the species
register methods:

- `ur_lay.l_max = ppr->l_max_ur;` (line 3444) moves into
  `UltraRelativisticSpecies::RegisterPerturbationIndices`. Photons already
  read `ppr->l_max_g` themselves (`species/photons.cpp:25`), so this just
  evens the pattern.
- The NCDM pre-loop block setting `ncdm_lay.l_max = ppr->l_max_ncdm;
  ncdm_lay.q_size = nsp->q_size();` (lines 3577-3584) moves into
  `NCDMSpecies::RegisterTensorPerturbationIndices` (the species knows its own
  `q_size`). The accompanying `dynamic_cast<const NCDMSpecies*>` scan goes
  away.

### B.5 Baryons vector mode

The inline `b_vec_lay.idx_theta = index_pt; ++index_pt;` block at lines
3502-3509 moves into `BaryonsSpecies::RegisterVectorPerturbationIndices`
(currently absent on Baryons; it gets a non-stub body).

### B.6 The collapsed call sites

Scalar mode (replacing lines 3338-3476):

```cpp
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->RegisterPerturbationIndices(
      *ppv->species_layouts[i], ppv, ppr, index_pt, ppw, ppt->gauge);
}

// perturbed_recombination — module-owned, not a species
if (ppt->has_perturbed_recombination == _TRUE_ &&
    ppw->approx[ppw->index_ap_tca] == (int) tca_off) {
  class_define_index(ppv->index_pt_perturbed_recombination_delta_temp, _TRUE_, index_pt, 1);
  class_define_index(ppv->index_pt_perturbed_recombination_delta_chi,  _TRUE_, index_pt, 1);
}
```

Vector mode (replacing lines 3502-3520):

```cpp
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->RegisterVectorPerturbationIndices(
      *ppv->species_layouts[i], ppv, index_pt, ppw, ppt->gauge);
}
```

Tensor mode (replacing lines 3543-3594, minus `tensor_ur_layout` which stays
module-side):

```cpp
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->RegisterTensorPerturbationIndices(
      *ppv->species_layouts[i], ppv, index_pt, ppw, ppt->gauge);
}
```

## C. "Used in sources" → species own their masks

Today `MarkUsedInSources` on `BaseSpecies` is a no-op default with the
signature `(const PerturbLayout&, int* used_in_sources)`, called once in a
dispatch loop (`perturbations_module.cpp:3628-3630`). After this loop, the
module runs four bespoke approximation-conditional masking blocks: photons
(scalar, lines 3637-3655), UR (3658-3670), IDM_DR_IDR (3672-3692), NCDM
(3694-3713), and photons (tensor, 3722-3741).

The single method splits into three, matching the
`Register{,Vector,Tensor}PerturbationIndices` pattern:

```cpp
virtual void MarkUsedInSources(const PerturbLayout& /*layout*/,
                               const perturb_workspace* /*ppw*/,
                               int* /*used_in_sources*/) const {}

virtual void MarkVectorUsedInSources(const PerturbLayout& /*layout*/,
                                     const perturb_workspace* /*ppw*/,
                                     int* /*used_in_sources*/) const {}

virtual void MarkTensorUsedInSources(const PerturbLayout& /*layout*/,
                                     const perturb_workspace* /*ppw*/,
                                     int* /*used_in_sources*/) const {}
```

(`MarkUsedInSources` keeps the existing scalar name but gains the `ppw`
argument it needs for approximation-conditional masking.)

Each affected species overrides the appropriate method:

- `PhotonsSpecies::MarkUsedInSources` (scalar): mask l3..lmax and pol1,
  pol3..lmax_pol when `rsa_off && tca_off`.
- `PhotonsSpecies::MarkTensorUsedInSources`: mask theta, l3, l5..lmax,
  pol1, pol3, pol5..lmax_pol when `rsa_off && tca_off`.
- `UltraRelativisticSpecies::MarkUsedInSources` (scalar): mask l3..lmax when
  `rsa_off && ufa_off`.
- `IDM_DR_IDR_Species::MarkUsedInSources` (scalar): mask idr l3..lmax when
  `rsa_idr_off && tca_idm_dr_off && idr_nature == free_streaming`.
- `NCDMSpecies::MarkUsedInSources` (scalar): mask per-q l>2 when
  `ncdmfa_off`.

The existing default-`_TRUE_`-then-mask split (set everything true, then
species mark slots false) is preserved.

Module side, each mode-block in `perturb_vector_init` calls the matching
loop:

```cpp
// scalar block
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->MarkUsedInSources(
      *ppv->species_layouts[i], ppw, ppv->used_in_sources);
}

// tensor block
for (size_t i = 0; i < all_species_.size(); ++i) {
  all_species_[i]->MarkTensorUsedInSources(
      *ppv->species_layouts[i], ppw, ppv->used_in_sources);
}
```

(Vector mode has no current approximation-conditional masking in
`perturb_vector_init`; `MarkVectorUsedInSources` defaults to no-op for
every species. The method exists to keep the trio symmetric with
`Register{,Vector,Tensor}PerturbationIndices`.)

The standalone `ppv->used_in_sources[ppv->index_pt_gw] = _FALSE_;`
(line 3746) stays in the module (metric slot, not a species).

## Validation

1. `make` must succeed.
2. Run the `.ini` files in `test/` and `notebooks/` that exercise:
   - Plain LCDM (regression on most-common path).
   - UR + standard NCDM (regression on photon/UR/NCDM masks).
   - IDR + IDM_DR_IDR (regression on IDR mask + IDR composite).
   - DCDM_DR (regression on single-composite DR).
   - **DNCDM_DR ≥ 1** (this is where the y-vector reorder bites).
3. Compare Cl spectra with the master baseline at ~0.1% tolerance, with
   Cl^TE zero-crossings handled per `feedback_no_bit_identical_requirement`.

## Out-of-scope follow-ups (don't do in this PR)

- Remove the legacy `pv->index_pt_<species>_*` members on `perturb_vector`
  in favour of layout-only access. Many derivs / sources / initial-condition
  call sites still read these. Separate, larger cleanup.
- Promote `perturbed_recombination` slots onto `BaryonsSpecies`. Plausible
  long-term but unrelated to the three issues this PR addresses.
- The `tensor_ur_layout` block (`perturbations_module.cpp:3554-3567`) stays
  module-owned; promoting it to a species requires resolving whether the
  tensor UR hierarchy belongs to `UltraRelativisticSpecies` or is a
  mode-only construct shared by massless NCDM.

## Open questions

None at design time. Verification at ~0.1% is the only thing left to confirm
empirically.
