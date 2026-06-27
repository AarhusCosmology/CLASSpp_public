# Composite generic child iteration + per-child matter tally (`TallyStressEnergy`)

Date: 2026-06-27
Branch: `type3-momentum-transfer`
Follows: `2026-06-20-perturb-stress-energy-deltarho-matter-tally-design.md` (the
tally being amended), `2026-06-21-perturb-stressenergy-delegation-design.md`
(the single-call `StressEnergy()` struct), `2026-06-26-type3-momentum-transfer-design.md`.

## Problem

The matter tally introduced in the 2026-06-20 design classifies **each
`active_species` entry's whole `StressEnergyContribution`** using two booleans
cached once per `pv` on the entry (`perturbations.h:263-264`):

```cpp
struct ActiveSpecies {
  const BaseSpecies* species;
  const BaseSpecies::PerturbLayout* layout;
  bool clusters_as_matter;  // ClustersAsMatter()
  bool is_cold;             // IsColdMatterSpecies()
};
```

The hot loop (`perturbations_module.cpp:4853-4894`) sums each entry's six
stress-energy quantities into `ppw`, then — if the entry is tagged — folds its
`ρ−3P` / `δρ−3δP` matter contribution into a `cold` or `warm` bucket.

For a **composite** these two booleans are OR-reductions over children
(`composite_species.cpp:96-108`): `ClustersAsMatter()` is true iff *any* child
clusters, `IsColdMatterSpecies()` true iff *any* child is cold. That is
meaningless for a sector that **straddles the matter / non-matter partition**.

**Type3** (`cdm_scf_momentum`) = coupled CDM (`EnergyType::Matter`, cold,
clusters) + quintessence scalar field (`EnergyType::Other`, not cold, does not
cluster). The composite-level booleans collapse to `clusters_as_matter = true`,
`is_cold = true`, so the **entire lumped** contribution (CDM **+** scalar field)
is dumped into the `cold` bucket.

Because the tally uses the background-consistent `ρ−3P` / `δρ−3δP` identity, the
pollution is exactly the scalar field's `ρ−3P` and `δρ−3δP`. Those are **zero
only for radiation** (`w = 1/3`): the scalar field has `w ≠ 1/3`, so its
contribution is nonzero and wrong-signed. The result is corrupted `delta_cb`,
`delta_m`, `theta_cb`, `theta_m`, and therefore corrupted `P_m(k)` / `P_cb(k)`
sources, **whenever the Type3 coupling is on and transfers/P(k) are requested.**

### Root cause

The tally granularity is the **top-level `active_species` entry**, but the
matter / cold partition is inherently **per child**. A composite's lumped
`StressEnergyContribution` **cannot be re-split after summation** — the per-child
information is gone the moment `Type3Species::StressEnergy` returns
`cdm + scf` (`type3_species.cpp:68-69`).

### Why the existing composites are unaffected

DCDM_DR (cold dcdm + dark radiation) and DNCDM_DR (warm dncdm + dark radiation)
*also* straddle, but their non-matter child is **radiation**, whose `ρ−3P` and
`δρ−3δP` are identically zero. So the whole-composite tally already equals the
matter-child-only contribution. Per-child tally reproduces their results
**exactly** (verified component-by-component in "Correctness" below). Only Type3,
whose second child is not radiation, is broken — and only it changes.

## Constraints (from the discussion)

1. **Hot path.** `perturb_total_stress_energy` runs inside every RHS evaluation.
   No added virtual dispatch per *plain-species* entry.
2. **`StressEnergy()` stays a pure function** — it must not mutate `ppw`. It is
   reused un-tallied at the IC `delta_tot` loop (`perturbations_module.cpp:4016`).
   Mutating `ppw` inside `StressEnergy` would push the tally burden onto every
   species author and destroy the interface.
3. **No double-evaluation** of `StressEnergy`.
4. **The composite stays the entity that represents the sector to the module**
   (no flattening of children into module-visible `active_species` entries).

## Design

Three parts. Part A is reusable infrastructure the whole composite pipeline
wants; Part B is the tally fix built on it; Part C retires the `ρ−3P` / `δρ−3δP`
matter proxy, which the per-child gate of Part B makes unnecessary.

### Part A — generic composite child + layout iteration

Today every `Type3Species` forwarding method is identical boilerplate:
`static_cast` the layout, then call the same method on `cdm_`/`scf_` with
`my.cdm`/`my.scf` (`type3_species.cpp:31-80`). The piece that blocks a generic
loop is that the child **layouts** are concrete-named value members
(`Type3Species::PerturbLayout { cdm; scf; }`, `type3_species.h:22-25`).

Introduce a layout-aligned handle on the children:

```cpp
// composite_species.h
class CompositeSpecies : public BaseSpecies {
 public:
  struct PerturbLayout : BaseSpecies::PerturbLayout {
    // Aligned 1:1 with children_ (owning). Built by CreatePerturbLayout below.
    std::vector<std::unique_ptr<BaseSpecies::PerturbLayout>> child_layouts;
  };

  std::unique_ptr<BaseSpecies::PerturbLayout> CreatePerturbLayout() const override {
    auto l = std::make_unique<PerturbLayout>();
    l->child_layouts.reserve(children_.size());
    for (const auto& c : children_)
      l->child_layouts.push_back(c->CreatePerturbLayout());
    return l;
  }
  // ...
};
```

**Storage decision (settled):** owning `vector<unique_ptr<PerturbLayout>>`, not
inline members + non-owning views. The per-child pointer indirection in the hot
loop is dwarfed by the per-child virtual `StressEnergy` call already there; the
heap allocations happen once at `pv` construction (cold path), never per step.
This keeps `CreatePerturbLayout` and the forwarders **fully generic** — zero
per-composite boilerplate.

**Alignment invariant:** `children_` order == `child_layouts` order ==
construction order. Children are pushed in the concrete constructor (Type3: cdm
then scf, `type3_species.cpp:14-18`); `CreatePerturbLayout` runs afterward, so
`children_` is populated.

**Generic forwarding defaults** on `CompositeSpecies`, each a loop over
`children_[i]` with `*child_layouts[i]`:

- `RegisterPerturbationIndices`
- `ApplyInitialConditions`
- `StressEnergy` (sum of children's contributions)
- `PerturbDerivs` (free-streaming loop, **then** the existing `AddCouplingDerivs`
  hook — preserves the two-phase contract in `composite_species.h:22-24`)
- `FillSources`, `CopyPerturbationsAcrossSwitch`,
  `PerturbSynchronousToNewtonian`, `MarkUsedInSources`
- `DelegateTally` (Part B)

A concrete composite **only** overrides what is genuinely coupling- or
type-specific. For Type3 that is: `AddCouplingDerivs`; the `StressEnergy`
cross-term once Task 5 lands (override = generic child sum **+** the
`-(2β/3) Z̄² θ_cdm` term, so it reuses the base sum and adds one line); the
shooting forwarders to `scf_`; `kTypeName`; `CreateAll`; the background column
writers. Type3's `PerturbLayout` loses its `cdm`/`scf` members; coupling code
that needs typed fields reaches them by index with a cast helper:

```cpp
enum ChildIndex { kCdm = 0, kScf = 1 };
static CDMSpecies::PerturbLayout& cdm_layout(const PerturbLayout& my) {
  return static_cast<CDMSpecies::PerturbLayout&>(*my.child_layouts[kCdm]);
}
```

Casting `*child_layouts[k]` to the concrete child layout is the composite
operating on **its own** children — not module code downcasting a species — so it
respects [[feedback_no_species_picking_in_modules]].

### Part B — `TallyStressEnergy` (the tally fix)

A **non-virtual** wrapper on `BaseSpecies` fuses the totals accumulation and the
matter tally into one call that invokes the virtual `StressEnergy()` **exactly
once**. Composites delegate the *whole* call to each child, so each child is
evaluated once (no double-eval) and tallied with its **own** booleans (no
lump-splitting).

```cpp
// base_species.h — declared AND defined inline here. The body only PASSES ppw
// through to StressEnergy (never dereferences it) and accumulates into caller-
// owned StressEnergyContribution structs, so it compiles with perturb_workspace
// merely forward-declared — no include of perturbations.h, no layering inversion.
class BaseSpecies {
 protected:
  bool clusters_as_matter_ = false;  // own classification, cached once (see "Cached own-booleans")
  bool is_cold_            = false;
  bool delegates_tally_    = false;  // CompositeSpecies sets true in its ctor

 public:
  // Non-virtual: a direct, inlinable call at the loop site (e.species is BaseSpecies*).
  // Only the inner StressEnergy() is virtual -> one indirect call for plain species,
  // identical to today. Composites take the delegates_tally_ branch instead.
  // The three accumulators are owned by the module loop; we reuse
  // StressEnergyContribution (and its operator+=, base_species.h:340-348) as the
  // bucket type, so there is no bespoke MatterTally struct and no ppw bucket state.
  void TallyStressEnergy(const PerturbLayout& layout, const perturb_vector* pv,
                         const double* y, const double* pvecback, const perturb_workspace* ppw,
                         StressEnergyContribution& total,
                         StressEnergyContribution& total_cold,
                         StressEnergyContribution& total_warm) const {
    if (delegates_tally_) {
      DelegateTally(layout, pv, y, pvecback, ppw, total, total_cold, total_warm);
      return;
    }
    const StressEnergyContribution se = StressEnergy(layout, pv, y, pvecback, ppw);
    total += se;
    // Actual ρ/δρ/(ρ+P)θ — no ρ−3P proxy (Part C). Radiation never reaches here
    // (it does not cluster), so there is nothing to "zero out".
    if (clusters_as_matter_)
      (is_cold_ ? total_cold : total_warm) += se;
  }

  // Seam: only composites (delegates_tally_ == true) override; base is unreachable.
  virtual void DelegateTally(const PerturbLayout& /*layout*/, const perturb_vector* /*pv*/,
                             const double* /*y*/, const double* /*pvecback*/,
                             const perturb_workspace* /*ppw*/,
                             StressEnergyContribution& /*total*/,
                             StressEnergyContribution& /*total_cold*/,
                             StressEnergyContribution& /*total_warm*/) const {}
};
```

```cpp
// composite_species.cpp — the one generic delegation, reused by every composite.
void CompositeSpecies::DelegateTally(const PerturbLayout& base, const perturb_vector* pv,
                                     const double* y, const double* pvecback,
                                     const perturb_workspace* ppw,
                                     StressEnergyContribution& total,
                                     StressEnergyContribution& total_cold,
                                     StressEnergyContribution& total_warm) const {
  const auto& my = static_cast<const PerturbLayout&>(base);
  for (size_t i = 0; i < children_.size(); ++i)
    children_[i]->TallyStressEnergy(*my.child_layouts[i], pv, y, pvecback, ppw,
                                    total, total_cold, total_warm);
}
```

**Virtual-call accounting (the whole point):**

| entry kind | indirect calls / entry | vs today |
|------------|------------------------|----------|
| plain species | `TallyStressEnergy` is **direct/inlined** → 1 (`StressEnergy`) | unchanged |
| composite (Type3) | 1 (`DelegateTally`) + N children × 1 (`StressEnergy`) = 3 | same total work as today's lumped path; no double-eval |

The non-virtual wrapper is what keeps plain species at one indirect call. A
*virtual* `TallyStressEnergy` would cost two (outer dispatch + inner
`StressEnergy`), because a base-class method is compiled once for all derived
types and cannot devirtualize `this->StressEnergy()`.

**Cached own-booleans (`clusters_as_matter_`, `is_cold_`).** Set once,
post-construction (NOT in the constructor — `ClustersAsMatter()` /
`IsColdMatterSpecies()` are virtual and composites scan children, so they must
run after the object and its children are fully built). A small finalize pass
over `all_species_` (recursing into composite `children_`) calls the virtuals
once and stores the result. The species objects are shared read-only across
threads thereafter, so this is thread-safe.

**Module hot loop** (`perturb_total_stress_energy`, scalars block) becomes:

```cpp
StressEnergyContribution total, total_cold, total_warm;   // module-owned accumulators

// A canonical scalar field's Newtonian StressEnergy reconstructs psi from the
// running shear in ppw->rho_plus_p_shear (active_species orders scalar fields
// last). Expose that running shear from `total` after each species.
ppw->rho_plus_p_shear = 0.;
for (const auto& e : pv->active_species) {
  e.species->TallyStressEnergy(*e.layout, pv, y, pb, ppw, total, total_cold, total_warm);
  ppw->rho_plus_p_shear = total.rho_plus_p_shear;
}

ppw->delta_rho        = total.delta_rho;
ppw->rho_plus_p_theta = total.rho_plus_p_theta;
ppw->delta_p          = total.delta_p;
ppw->rho_plus_p_tot   = total.rho + total.p;
// ppw->rho_plus_p_shear already holds total.rho_plus_p_shear (last sync above).

// matter outputs (only consumed when transfers/P(k) requested):
//   delta_cb = total_cold.delta_rho / total_cold.rho
//   theta_cb = total_cold.rho_plus_p_theta / (total_cold.rho + total_cold.p)
//   delta_m  = (total_cold.delta_rho + total_warm.delta_rho)
//                  / (total_cold.rho + total_warm.rho)
//   theta_m  = (total_cold.rho_plus_p_theta + total_warm.rho_plus_p_theta)
//                  / ((total_cold.rho + total_cold.p) + (total_warm.rho + total_warm.p))
// + the PPF fluid closure (unchanged).
```

**No workspace additions, no bespoke bucket type.** The three accumulators are
**module-loop locals**, and we reuse `StressEnergyContribution` (with its existing
`operator+=`, base_species.h:340-348) as the bucket type — so `total_cold` /
`total_warm` carry ρ, P, δρ, (ρ+P)θ in named fields and we read `delta_rho`,
`rho`, `rho_plus_p_theta`, and `rho+p` straight off them. The two buckets that are
function-locals today (`4838-4840`) stay locals; they just change type from the
ad-hoc `struct Tally` to `StressEnergyContribution`. There is **no** `MatterTally`
struct and **no** `do_matter_tally` flag: the buckets are always accumulated (a
handful of adds per clustering species), and the post-loop code ignores them
unless `delta_m`/`theta_m` are requested.

**Why this matters for layering.** Because `TallyStressEnergy` only *passes* `ppw`
through to `StressEnergy` and never dereferences it, the body needs only
`StressEnergyContribution` (defined in `base_species.h`). So it is a normal inline
member of `base_species.h`, with `perturb_workspace` left forward-declared. An
earlier sketch had the body write `ppw->delta_rho += …` directly, which needs the
*complete* `perturb_workspace` and would have forced the definition into
`perturbations.h` — an unacceptable layering inversion. Passing the accumulators by
reference (the three-struct design) avoids that and is also less state to manage.

**Scalar-field shear read-back.** The one quantity a species reads back from `ppw`
mid-loop is `rho_plus_p_shear` (a canonical scalar field's Newtonian `StressEnergy`
reconstructs ψ from it; that is why `active_species` orders scalar fields last,
`4843-4852`). The module keeps `ppw->rho_plus_p_shear` synced from the running
`total` after each species (see the loop above), so a last-ordered scalar field
still sees the full accumulated shear — behaviour identical to today.

**`ActiveSpecies` shrinks.** The `clusters_as_matter` / `is_cold` fields move
onto the species (cached) and are removed from the entry struct. They are
currently only consumed by this tally; the derivs loops destructure but do not
use them (`perturbations_module.cpp:5871, 5931, 5947`) — verify, then drop.

### Part C — remove the `ρ−3P` / `δρ−3δP` matter proxy

The 2026-06-20 design defined each species' matter contribution as the
background-consistent `ρ−3P` (density) and `δρ−3δP` (perturbation), so that
radiation contributes 0 (`ρ−3·ρ/3`) without an `IsMatterSpecies` hook. With the
per-child gate (Part B), **radiation is already excluded** — it does not cluster
— so the proxy's only remaining effect is to **pressure-subtract warm species**,
which is wrong: structure formation is sourced by the full `δρ`, and `δP` for a
warm species is an independent variable, not noise to be cancelled.

So the matter tally uses the **actual** stress-energy of each clustering species:
`ρ`, `δρ`, `(ρ+P)θ`, `(ρ+P)` — all already computed for the Einstein totals, so
this is *fewer* operations than the proxy (no subtraction, no `δP` term in the
tally). `se.delta_p` is still summed into `ppw->delta_p` for the Einstein totals;
it simply no longer enters the matter tally.

This is purely the **perturbation** matter tally (`delta_m/cb`, `theta_m/cb`).
The **background** `rho_m` budget split (which adds `ρ−3P` for `EnergyType::Other`
at `background_module.cpp:344`) is a separate axis-1 concern and is **left
untouched** here.

### Ordering note (scalar fields last)

`active_species` keeps `ScalarFieldSpecies` entries last so their Newtonian
`StressEnergy` sees the fully accumulated shear (`perturbations_module.cpp:4847`,
ordering at `2997-3006`). Type3 is a `CompositeSpecies`, not a
`ScalarFieldSpecies`, so it is **not** sorted last — but Type3 is
**synchronous-gauge only** (it throws in Newtonian, `type3_species.cpp:213-217`),
and the shear-reconstruction path is Newtonian-only. So no ordering change is
needed. A future Newtonian-capable composite containing a scalar field would need
the ordering predicate generalized to "contains a scalar field" — **out of
scope** here.

## Rejected alternatives (recorded)

- **Flatten composite children into `active_species`.** Destroys the composite
  as the sector's representative, and `active_species` is shared with the
  derivs/coupling dispatch (`PerturbDerivs` runs children **then**
  `AddCouplingDerivs` as a unit, `5871`); flattening would drop the coupling.
- **Tally off the composite lump via a precomputed-`se` helper.** Impossible —
  the per-child split is gone after summation.
- **Re-run children's `StressEnergy` in a separate tally hook.** Double-evaluates
  the most expensive call on the hot path.
- **All-virtual `TallyStressEnergy`.** Correct, but +1 indirect call for every
  plain-species entry; the non-virtual wrapper + `delegates_tally_` removes it.
- **Inline-member layouts + non-owning view vector.** Faster hot path, but
  per-composite boilerplate and no generic `CreatePerturbLayout`; rejected in
  favor of uniform genericity (the indirection is immaterial here).

## Correctness (invariants — must NOT change beyond ULP)

- **Einstein totals** (`delta_rho / rho_plus_p_theta / rho_plus_p_shear /
  delta_p / rho_plus_p_tot`) are the same sums. A composite that previously added
  `(cdm + scf)` in one `+=` now adds `cdm` then `scf` in two; addition is
  associative, so this differs only at the ULP level (amplified by the ODE),
  accepted under the ≤0.1 % bar ([[feedback_no_bit_identical_requirement]],
  [[feedback_vectorization_reduction_drift]]). Part C does not touch the totals.
- **`delta_cb` / `theta_cb`** (cold bucket) — cold species are pressureless
  (`P=0`), so `ρ−3P=ρ` and `δρ−3δP=δρ` already: dropping the proxy is a **no-op**
  for the cold bucket. Unchanged.
- **Purely-cold cosmologies** (ΛCDM, ΛCDM+baryons, DCDM_DR — dcdm is cold, DR is
  excluded by the cluster gate): `delta_m`/`theta_m` unchanged. The loop body for
  plain cold species is logically identical (same single `StressEnergy` call,
  same accumulation order), so these runs are byte-stable modulo ULP reordering.
- **Type3 fix:** the scalar-field child no longer enters any matter bucket (it
  neither clusters nor is cold), so the cold tally is CDM-only — and the CDM
  child being cold, its contribution is the full `ρ_cdm`/`δρ_cdm`, as required.

## Behavior changes (intentional — regenerate reference, document)

Removing the `ρ−3P` / `δρ−3δP` proxy (Part C) changes the matter tally **only for
warm species**, where it previously pressure-subtracted:

- **Standalone NCDM / DNCDM, and the DNCDM_DR composite:** `delta_m` and
  `theta_m` now use the species' **full** `δρ` / `(ρ+P)θ` (its true clustering
  contribution) instead of the `ρ−3P`-weighted proxy. Magnitude ~`3w·f_warm`
  (sub-percent at z=0 for realistic neutrino masses, larger at high z and in the
  high-`k` matter power). This is the "ncdm slightly, to the better" change: the
  matter power spectrum now sees the full neutrino over/under-density rather than
  a pressure-cancelled fraction.
- **`delta_cb`/`theta_cb` (cold) and all non-`mPk`/non-transfer runs: unchanged**
  (the tally is gated by `has_source_delta_m_ || has_source_theta_m_`).
- The **IDM_DR P(k) fix** from the 2026-06-20 design is preserved (IDM_DR is
  `EnergyType::Matter`, so it clusters and enters the tally with its full
  density).

## Verification

1. Build clean (CMake) Debug (asserts) + Release.
2. **Type3 fix:** a `scf_veta`-on run with `mPk`/transfer output. Confirm the
   cold/matter tally excludes the scalar field — e.g. cross-check `delta_cb`
   against the CDM child's `δρ_cdm/ρ_cdm` directly, and confirm `P_m(k)` no
   longer carries the scalar field's `ρ−3P`. Compare coupling-off (`scf_veta=0`)
   to the master CDM+quintessence baseline (must match within ULP).
3. **Cold cosmologies (no regression):** ΛCDM, ΛCDM+baryons, and DCDM_DR (cold
   dcdm + excluded DR) — TT **and** mP(k) byte-stable modulo ULP vs master
   (TE zero-crossing aware, never blind max-rel-diff). Cold buckets are
   pressureless, so Part C is a no-op for them.
4. **Warm species (intentional change, Part C):** standalone massive-NCDM and
   DNCDM_DR — **TT ≤0.1 %** (the matter tally feeds only transfer/P(k) sources,
   not the ODE, so the CMB is unaffected beyond ULP), while **mP(k) / `delta_m`
   / `theta_m` change**; confirm the shift is the expected direction and
   magnitude (full `δρ` > pressure-subtracted proxy, ~`3w·f_warm`), then
   regenerate the reference.
5. Full `TEST_LEVEL=2` against regenerated `classyref`; regenerate goldens for
   any intentional Type3-config changes ([[reference_classyref_testing]]).
6. Re-profile `class_profiled` (Planck-2018 single-NCDM) to confirm the
   plain-species hot path is unchanged (one indirect call per entry).

## Risks

- **`delegates_tally_` set but `DelegateTally` not overridden** → base
  `DelegateTally` asserts/throws; a generic `CompositeSpecies::DelegateTally`
  default covers every composite, so only a non-composite mis-setting the flag
  could trip it.
- **`child_layouts` / `children_` misalignment** → single construction-order
  invariant; assert `child_layouts.size() == children_.size()` in
  `CreatePerturbLayout`.
- **Finalize pass misses a species** (stale `clusters_as_matter_`/`is_cold_`) →
  drive it from the same `all_species_` walk that builds `active_species`, and
  recurse composites; assert it ran before the first RHS.
- **`ActiveSpecies` bool still read somewhere** → grep before removing the
  fields (the derivs structured bindings name but don't use them).
- **Type3 typed-layout access after dropping `cdm`/`scf` members** → the
  `cdm_layout()`/`scf_layout()` cast helpers are the single chokepoint; the
  coupling code (Task 5) uses only those.
