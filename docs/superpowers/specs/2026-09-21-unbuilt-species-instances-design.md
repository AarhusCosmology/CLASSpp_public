# A dot-syntax species instance nobody builds is an input error

**Status (2026-09-21).** Shipped: `RejectUnbuiltSpeciesInstances`
(`species/species_input.h`), called once at the end of
`InputModule::ConstructSpecies`. It rejects an unrecognised `<instance>.type` (#430
fix (2)), and also the case where an instance has fields but no `.type` line at all.
**Not shipped:** fix (1) of #430, making the CLI reject *every* unread key the way
classy does. That stays a separate decision; the measurement below is there for it.
Unread fields of an instance that *was* built (`nuM.mass` for `nuM.m`) are also not
checked here, because they belong to (1).

## The bug

```ini
nuM.type = ncdm          # the type is `ncdm_standard`
nuM.m    = 0.06
```

`./class` ran to completion, wrote output, and printed nothing. The background table
had no `rho_nuM` column: the massive neutrino was simply absent. **[M]** (reproduced
on `7cc41a51`: exit 0, no message). Each factory in `kAllSpeciesFactories` asks
`instances_with("type", kTypeName)` for its own type. A value that matches none of
them matches nothing, and every `nuM.*` key is left unread. The CLI reports unread
keys only under `write warnings = yes`, which is off by default. classy raised
`Class did not read input parameter(s)`, so the two front ends disagreed.

A wrong cosmology reported as success is the worst failure mode this code has. This
one cost about an hour in the H0 Olympics work: two species that should have been
identical differed by 22% in lensed TT, and the cause was that one of them had never
been built.

## The check

After every factory has run, look at every dot-syntax key:

| input | verdict |
|---|---|
| `N.type` read | a factory built `N`: fine |
| `N.type = V` unread, `V` not in `kAllSpeciesFactories` | not a species type; lists the species types |
| `N.type = V` unread, `V` a species type | `V` cannot be declared with `.type` yet; use its legacy keys |
| `N.type` with `N` not `[A-Za-z_][A-Za-z0-9_]*` | illegal instance name (rejected even if read) |
| `N.<field>` and no `N.type` anywhere | no species is built from these keys |

All offenders go into one `class_stop_severe`. The check is structural: it is decided
by the text of the input, never by a value a sampler varies. So it goes in the severe
channel (`include/errors.h`), which becomes `CosmoSevereError` in classy.

**The contract it relies on** is written at `kAllSpeciesFactories`: a factory reads
`N.type` for exactly the instances it builds. Every dot-syntax factory already did.
`CreateAllNcdmInstances`, `DNCDMSpecies`, `DCDM_WDM_Species` and `DrPsdSpecies` read
it inside their per-instance loops, and `TranslateSingleInstanceDotSyntax` reads it
for photons, baryons, cdm, ur, lambda and fluid. **[C]** The legacy `N_ncdm` path
stamps `ncdm__N.type` on every instance it synthesises, so it passes. **[C]**, and
**[M]** in the test.

## Why read status, not a name check

Checking `V` against the registry is the obvious one-liner, and it is not enough.
Five registered types never read `N.type`: `dcdm_dr`, `idm_dr_idr`,
`idm_drmd_idr_drmd`, `scalar_field` and `cdm_scf_momentum` are built from legacy keys
only. **[C]** `x.type = scalar_field` would pass a name check and still build
nothing. That is the same silent absence, with a correctly spelled name. Whether the
key was read is the direct statement of "a species was built from this instance".
It needs no list that could drift from the factories, and a future species is covered
without touching the check.

The registry is still used for the message, to tell "not a type" apart from "a type
that is not declarable". The unknown-type message lists every registered type,
including the five legacy-only ones. Listing only the declarable ones would need a
second list kept in step with the factories by hand. Choosing a legacy-only type
from the list gets its own message, so the listing is never a dead end.

**Rejected: a "did you mean".** Edit distance ranks `cdm` (distance 1) above
`ncdm_standard` (distance 9) for the very input that motivated this. The full list
is short and does not mislead.

**Rejected: checking only `N.type`.** `nuM.typ = ncdm_standard`, or a forgotten type
line, fails in exactly the same way. Nothing reads a dotted key whose instance has no
`.type` line: there are no non-species dotted keys, and the only key-pattern scan
(`.fluid_approximation`) leaves such keys unread. **[C]** So rejecting them changes
nothing for any input that currently works.

## Front ends

* **CLI:** such an input now exits 1 with the message above, where before it
  "succeeded".
* **classy:** the same inputs were already refused by the unread-key check in
  `reset()`. The exception now comes from `Cosmology` construction, before that
  check runs, so the error names the cause instead of listing unread keys. No input
  changes from accepted to refused.

## For fix (1): what "reject every unread key" would cost the CLI

Building the `InputModule` alone is where classy's unread check looks. Doing that for
every tracked `.ini` **[M]** (2026-09-21):

* 45 of the 47 files that build read every key. (`python/pytest.ini` is not a
  CLASS input; `test/scenarios/dncdm_dr_bare_omega.ini` is a negative scenario that
  stops on its own guard.)
* `test/scenarios/type3_scf_veta.ini` leaves `l_max_scalars` unread.
* `explanatory.ini` leaves 74 keys unread. Most of them are features that are not
  enabled there (`w0_fld`, `scf_parameters`, `binned_reio_z`, ...). The file
  documents keys; it is not a run.

So (1) would need a way for `explanatory.ini` to stay runnable (an opt-out, or
commenting out the inactive blocks), and one scenario fixed. That is small, but it is
a change of behaviour for every user `.ini`, which is why it is not bundled here.

## Verification

* `species/species_instance_check_test.cpp` (`test-species-instances`). Unit cases
  on a bare `FileContent` cover all five rows of the table, grouping, key order and
  one message for several offenders. `InputModule` cases run the real factories: the
  #430 input verbatim now throws; with `ncdm_standard` it builds `nuM`; single-instance
  `photons`, the legacy `N_ncdm` path, `scalar_field` and `nuM.typ` all behave as the
  table says. With the call in `ConstructSpecies` commented out, the #430 case fails.
  **[M]**
* All 49 tracked `.ini` files and 18 untracked local ones were run through
  `InputModule` construction with the check live. None tripped it. 47 tracked files
  built; the other two fail as they did before, for the reasons given above. **[M]**
* `make test` passed 44/44. The CI quick scenario run (`TEST_LEVEL=1`,
  `-m test_scenario`) passed 84/84. The non-scenario `test_class.py` run gave 96 passed and
  5 failed. The five are reference comparisons against the local `classyref`, and
  master `7cc41a51` fails the same five with identical numbers (rs_drag
  147.38426900 vs 147.38405632), so the local `classyref` is stale. **[M]**
