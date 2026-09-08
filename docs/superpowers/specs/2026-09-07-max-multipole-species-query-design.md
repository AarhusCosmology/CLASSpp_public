# s_l is sized by asking the species, not by naming them

**Status (2026-09-07).** Shipped in full: `BaseSpecies::MaxMultipole(const
precision*, bool tensors)`, overrides on the six species families that own
hierarchies, a `CompositeSpecies` override that scans children, and
`perturb_workspace_init` reduced to a loop. Fixes #421. Output-neutral: verified
bit-identical on `dcdm_dr` at `l_max_dr` = 17 and 200 **[M]**.

## The bug

`ppw->max_l_max` sized `ppw->s_l` from a hard-coded list of species **keys**:

```cpp
if (all_species_.count("DCDM_DR"))
  ppw->max_l_max = std::max(ppw->max_l_max, ppr->l_max_dr);
```

A DNCDM composite also owns a `DarkRadiationSpecies` daughter running to
`ppr->l_max_dr` (`species/dncdm_dr_species.cpp:81`), but it is keyed by its
instance name (`dncdm1`). The test missed it, `s_l` was allocated 18 entries, and
the daughter's free-streaming loop reads `s_l[l+1]` for `l < l_max_dr`, so it indexed
up to `s_l[l_max_dr]`.

Any dncdm run with `l_max_dr > 17` read heap garbage: silently wrong spectra, and
nondeterministic `sp_ludcmp` / "step size too small" aborts at roughly 50% on
identical input which **survived `OMP_NUM_THREADS=1`** — so not a race. ASan named
it at `dark_radiation_species.cpp:107` against the allocation at
`perturbations_module.cpp:1883` — a 144-byte region, i.e. 18 doubles. The report
reads "0 bytes after" that region because the *first* out-of-bounds access is
`s_l[18]`, at `l = 17`, not the eventual `s_l[l_max_dr]`. **[M]**

The nondeterminism is the expensive part of this bug. It presents as integrator
fragility, which is a plausible-looking dead end for a stiff decaying species, and
it silently poisoned a truncation-convergence campaign before ASan was reached
for. **Run-to-run variation on identical input is a memory bug until proven
otherwise** — one thread, then ASan, before any numerical hypothesis.

## The fix

```cpp
virtual int MaxMultipole(const precision* ppr, bool tensors) const { return 0; }
```

Overridden by `PhotonsSpecies` (the only mode-dependent one),
`UltraRelativisticSpecies`, `NCDMBaseSpecies` (covering ncdm, dncdm, `dr_psd` and
the wdm daughters in one place), `DarkRadiationSpecies`, `IDRSpecies` (non-zero
only when free-streaming), and `CompositeSpecies`, which maxes over children —
the same shape as `SupportsExplicitPerturbationEvolver`. The module loops.

Over-reporting is safe: `s_l` is merely longer and the extra entries are never
read. Under-reporting is a heap overflow. The contract says so, so a species
author facing doubt has a correct default.

Two behaviour differences from the old key list, both widenings and both
deliberate: a standalone free-streaming `IDR` now contributes `l_max_idr` (before,
only one nested in an `IDM_DR_IDR` composite did), and any `NCDMBaseSpecies`
contributes `l_max_ncdm` rather than only those the collection's `has_ncdm()`
gate admits. **[C]**

## Why not the one-line patch

Adding `l_max_dr` unconditionally fixes this instance and leaves the mechanism.
The tensor branch immediately below already carries a comment about a *previous*
instance of the same sizing bug (`l_max_ncdm` vs `l_max_ur` when `N_ur = 0` with a
massive neutrino), which is the second time this exact class of defect has been
paid for. A module cannot know what storage a species needs; asking removes the
question.

## What is not covered

The tensor relativistic-neutrino hierarchy is owned by `perturb_vector`, not by a
species, and still contributes through the `evolve_tensor_ur_` flag in the module.
Moving it would mean giving that hierarchy an owner, which is a larger change than
this one. **[C]**
