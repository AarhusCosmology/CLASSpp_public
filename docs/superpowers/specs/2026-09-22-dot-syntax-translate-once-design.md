# Single-instance dot syntax is translated once per input

**Status (2026-09-22).** Shipped: `TranslateSingleInstanceDotSyntax` copies dot
fields into their legacy keys only on the first pass over a `FileContent`, marked
by `FileContent::single_instance_translated()`. Later passes still consume the dot
keys. This fixes `pk_eq = yes` with a dot-syntax fluid, found while fixing #430 and
noted in PR #437. It changes output for nothing else **[C]**, below.

## The bug

```ini
pk_eq    = yes
non_linear = halofit
Omega_Lambda = 0
f.type = fluid
f.w0   = -0.9
f.wa   = 0.1
```

This aborted with `input sets both legacy key 'w0_fld' and dot-syntax 'f.w0' with
different values`, while the same run spelled `w0_fld = -0.9, wa_fld = 0.1`
completed. **[M]** (CLI on master `1a8c5775`, classy on `7cc41a51`)

pk_eq finds, for each z, a constant-w model with the same distance to recombination.
It copies the input, overrides `w0_fld`/`wa_fld` with the trial effective values,
and builds a fresh `InputModule` from the copy. The copy already carries the first
translation's output, `w0_fld = -0.9` taken from `f.w0`. The new `InputModule` runs
the translation again. `f.w0` still says -0.9, `w0_fld` now holds the override, and
the check meant for a user who gave both spellings fires on the code's own rewrite.

## The fix

The translation rewrites the user's input in place, so after the first pass the
legacy key is the parameter and the dot key is the record of how it was spelled. A
reparse must not re-copy the dot value, because that would silently undo the
override. It must not re-check either, because the pair can only disagree through a
deliberate override. So a reparse only consumes `N.type` and the dot fields. They
still have to count as read, both for `RejectUnbuiltSpeciesInstances` (#430) and for
the unread-key warnings.

The once-only state is a flag on the `FileContent`, beside
`legacy_ncdm_transmuted()`. The NCDM transmutation hit the same reparse problem with
shooting (which made `100*theta_s` unusable with massive neutrinos) and solved it the
same way. The flag is copied with the `FileContent`, so pk_eq's copies and shooting's
copies both see it. It is not a parameter, so no input file can set it.

**Output-neutral elsewhere [C].** A reparse whose legacy keys were not overridden
used to re-copy the same values: a no-op. The only behaviour change is on a reparse
with an override. Today that is pk_eq alone: shooting's unknowns (`h`,
`Omega_ini_dcdm`, `scf_shooting_parameter`, `<instance>.deg`, ...) include no key
that has a single-instance alias.

**The trade-off.** After the first pass, a programmatic `set()` of a *dot* key on
that `FileContent` is ignored; to change a translated parameter, set the legacy key,
as pk_eq does. Nothing sets these dot keys programmatically **[C]**. classy builds a
fresh `FileContent` on every `reset()` **[C]**, so a user changing `f.w0` between
computes is translated afresh.

## Rejected

* **An alias-aware setter that pk_eq calls instead of `set()`**, updating the
  legacy key and any dot spelling of it. It works, but it leaves the trap for the
  next override site, and plain `set()` stays the obvious call.
* **Erasing the dot keys after translation.** This needs a `FileContent::erase`
  that nothing else wants, and `write parameters` would stop recording what the user
  wrote. Erasing `N.type` as well is worse: the reparse would see the instance's
  other fields with no `N.type` line, and the #430 check would reject as orphans keys
  that the first pass accepted.
* **Teaching pk_eq the fluid's instance name.** That is a module picking out a
  species by type, which modules here do not do.

**Not changed:** `SynthesiseIdenticalScalarField` (the family-wide
`ncdm_fluid_approximation` from `N.fluid_approximation`) has the same copy-then-check
shape, but nothing overrides `ncdm_fluid_approximation` programmatically **[C]**,
so nothing triggers it. Give it the same treatment if something ever does.

## Verification

* `tools/parser_test.cpp`: translate, copy, override `w0_fld`/`wa_fld`, then reparse
  the copy. The override survives and `f.type`, `f.w0`, `f.wa` are read. A fresh
  input giving both spellings still throws. Without the fix, the reparse throws the
  reported error. **[M]**
* CLI on the input above: exits 0. `pk.dat` and `pk_nl.dat` are byte-identical to
  the legacy-spelled run. pk_eq is live in that comparison: `pk_eq = no` moves
  P_nl(z=0) by up to 0.70%. **[M]**
* `python/test_class.py::test_dot_syntax_fluid_with_pk_eq_matches_legacy` passes on
  this branch and fails on a `7cc41a51` build with the reported error. **[M]**
* `make test` 44/44. Quick scenarios (`TEST_LEVEL=1`) 84/84. Non-scenario
  `test_class.py`: 97 passed and the same 5 stale-local-`classyref` failures as
  master, with identical numbers. **[M]**
