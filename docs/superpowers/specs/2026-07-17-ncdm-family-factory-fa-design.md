# NCDM-family shared factory loop + family-wide fluid_approximation gathering

**Issues:** #376 (per-instance `fluid_approximation` dot-key only gathered from
`ncdm_standard` instances) and #378 (factor the duplicated per-instance
`CreateAll` loop shared by NCDM-family factories). One PR; the two share the
"NCDM family enumeration" concern.

## Problem

1. `InputModule::input_read_precisions()` synthesizes the family-wide
   `ncdm_fluid_approximation` precision parameter from `<name>.fluid_approximation`
   dot keys, but only collects `ncdm_standard` instances. The parameter gates
   `index_ap_ncdmfa` for the *whole* family, and the flag is consumed by
   `NCDMSpecies` (hence also `ncdm_axion`, `ncdm_greybody` which inherit its
   perturbation code), `NCDMInteractingSpecies`, and the `ncdm_decay_dr`
   composite (`dncdm_dr_species.cpp:175`). So the dot key is silently ignored
   on four of the five types that honour it, and the "identical for all"
   guard misses them.
2. Four `CreateAll` factories repeat the same per-instance loop byte-for-byte
   (`ncdm_species.cpp`, `axion_ncdm_species.cpp`, `greybody_ncdm_species.cpp`,
   `ncdm_interacting_species.cpp`).

## Design decisions (issue #376 asked for these explicitly)

- **Contributing types** = exactly the types whose perturbation code consults
  `index_ap_ncdmfa`: `ncdm_standard`, `ncdm_axion`, `ncdm_greybody`,
  `ncdm_self_interacting`, `ncdm_decay_dr`.
- **`dcdm_wdm` (and any other species type) with a `fluid_approximation` dot
  key → hard error**, not silent acceptance. The daughter has no fluid
  approximation; silently accepting the key is the same bug class as #376.
  Rule: any `<name>.fluid_approximation` key whose instance has a `<name>.type`
  outside the consumer list throws `std::invalid_argument`. Keys without a
  matching `.type` stay in the generic unused-parameter warning path.
- The identity guard now spans the family: if any contributing instance sets
  the key, all must, and all values must agree (existing
  `SynthesiseIdenticalScalarField` semantics, unchanged).

## Approach

New species-owned pair `species/ncdm_family.{h,cpp}` — the single place that
knows the NCDM family roster:

1. `template <class T> std::vector<Named> CreateAllNcdmInstances(const SpeciesBuildContext&)`
   (header) — the shared loop: `instances_with("type", T::kTypeName)`, mark
   `.type` consumed, construct `T(pfc, name, settings, pba, bgm)`. The four
   factories become one-line bodies (NCDM keeps its legacy transmute pre-step,
   Interacting keeps `RejectLegacyInteractingKeys`). `DNCDMSpecies::CreateAll`
   is left alone: different `Named` type (concrete `unique_ptr<DNCDMSpecies>`)
   plus a per-instance closure post-step.
2. `void SynthesiseNcdmFluidApproximation(FileContent*)` (cpp) — gathers
   instances of the five consumer types, runs the rejection scan (via
   `for_each`), then delegates to `SynthesiseIdenticalScalarField` with
   description "dot-syntax NCDM-family species".
   `input_module.cpp:587-593` collapses to one call — the module no longer
   names a species type (per the no-species-picking rule).

Alternatives considered: (a) inline type-name list in `input_module.cpp` —
rejected, adds another string-name branch in module code; (b) a
`uses_ncdm_fluid_approximation` flag on `SpeciesFactoryEntry` — rejected,
couples the generic factory table to one precision parameter.

## Tests

New `species/ncdm_family_test.cpp` (executable `test-ncdm-family`, registered
in **both** CMakeLists.txt and Makefile `TEST_TARGETS`):

- #376 repro: axion-only instance with `fluid_approximation` → key synthesized.
- Family-wide identity guard: standard+axion disagreeing → throw;
  one-sets-one-doesn't → throw; agreeing → synthesized once.
- Rejection: `dcdm_wdm` instance with the dot key → throw.
- Conflict with an explicit global `ncdm_fluid_approximation` → throw.
- No dot keys anywhere → no-op.
- `CreateAllNcdmInstances<NCDMSpecies>` / `<AxionNCDMSpecies>`: instance count,
  keys, `.type` consumption (via `unread_parameters`).

## Verification

- `make test` (all existing + new).
- #378 output-neutrality: A/B `./class test/dotsyntax_ncdm_mixed.ini` and
  `test/dotsyntax_dncdm.ini` HEAD vs master — same build flags, expect
  identical output files (pure factory refactor, no numeric change).
- #376 behavior demo: `myaxion.fluid_approximation = 3` (ncdmfa_none) now
  changes Cl output vs the default CLASS method; on master it is a silent
  no-op.
