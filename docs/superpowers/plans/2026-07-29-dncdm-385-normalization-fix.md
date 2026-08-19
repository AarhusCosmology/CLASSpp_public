# Plan: #385 normalization reconciliation (DNCDM inverse decays)

Fixes issue #385 on branch `dncdm-inverse-qs` (PR #387). Design decision made by
Thomas Tram 2026-07-29: **occupation-suppression convention** (option 1 below).

## Context (read once)

`DecayTransitionKernel` (tools/decay_transition_kernel.{h,cpp}) physically requires
**bare per-dof occupation numbers** for all three species: its band factor
`Λ = f_l f_φ − f_H + f_H(f_l − f_φ)` comes from the QS collision algebra
(−f_H(1−f_l)(1+f_φ) + f_l f_φ(1−f_H), cubic cancels), which only holds in that
convention, and its exact discrete conservation identities are g-weighted with the
model's dof structure (g_H, g_l, g_φ) = (2, 2, 1):
`N_H + N_l = 0`, `2N_H + N_φ = 0`, `2E_H + 2E_l + E_φ = 0`
(the boson ×2 is folded into df_φ at deposit; the kernel is CORRECT and is
**frozen — do not modify the two kernel files**).

Two confirmed convention bugs (issue #385):

- **Bug A**: the parent (`DNCDMSpecies` on `NCDMBaseSpecies`) stores
  `f_stored = 1/(2π)³ · [FD(q−ξ) + FD(q+ξ)]` (ncdm_base_species.cpp:271-274), but
  `DNCDMInvSpecies` passes it to the kernel raw (dncdm_inv_species.cpp:190-198
  background; :252-254 perturbations exponentiate `ln f̄_H` with no conversion).
  Daughters (`DrPsdSpecies::ThermalF0`) store bare occupation. Ratio ≈ 1/124 at ξ=0.
- **Bug B**: the daughters' `factor_` is built from `SetDegAndFactor(2 or 1)` with
  `T = dr_T` (dr_psd_species.cpp:37), i.e. deg·T⁴ semantics **without** the 1/(2π)³
  that pairs with bare-f storage, and with T untied to the parent's. Physical
  (factor_-weighted) conservation then fails whenever shooting fits `deg_H`.

**Chosen fix (occupation suppression):** the shooting-fit `deg_H` means occupation
suppression. Conversion at the kernel boundary:

    κ = deg_H · (2π)³ / 2          (g_H = 2)
    f_bare = κ · f_stored          (in);   df_stored = df_bare / κ   (out)

Sanity anchors: deg_H = 1 ⇒ f_bare = FD exactly (a standard-abundance neutrino is
fully occupied); deg_H = 0.015 (Neff_ini fit) ⇒ f_bare = 0.015·FD, an under-occupied
FD that inverse decays drive toward the full-FD detailed-balance endpoint.

Daughters: factor semantics become `g/(2π)³` with **T pinned to the parent's** —
`SetDegAndFactor(2/(2π)³)` fermion, `SetDegAndFactor(1/(2π)³)` boson. Then each
species' `factor_ ∝ deg·T⁴` weights make physical conservation follow from the
kernel's exact (2,2,1)-weighted bare identities:
parent weight `deg_H·T_H⁴·(1/κ) = 2·T_H⁴/(2π)³ ∝ g_H`, fermion `2T_H⁴/(2π)³ ∝ g_l`,
boson `T_H⁴/(2π)³ ∝ g_φ`. Daughter factors are constants (no shooting coupling).

**Why decay-only regression cannot regress:** `DNCDMInvSpecies` is only constructed
when `inverse_decays=yes` (dncdm_dr_species.cpp:21-22). The `--decay-only` /
`--pert-regression` A/B paths never build the composite. Additionally, in a
decay-only composite (test-built) Λ = −f_H is linear, so κ cancels structurally in
the parent's Ψ-space cancellation.

## Global Constraints

- **Never modify** `tools/decay_transition_kernel.h` / `.cpp`.
- Never run `git checkout`, `git reset`, `git stash`, `git clean` — shared working tree.
- Never `git add -A` / `git add .`; stage explicit paths only.
- Build: `cmake --build build/cmake --target <targets> --parallel 8`; run tests via
  `ctest --test-dir build/cmake -R <regex> --output-on-failure`.
- C++ only, no `extern "C"` / `__cplusplus` guards.
- Error severity (#382): `class_test_severe` ONLY for structural checks decidable
  from which keys are set / flag combinations; any check on a numeric value a
  sampler could vary uses `class_test`.
- `OMP_NUM_THREADS=1` for any `class` executable run with `inverse_decays=yes`
  (kernel perturbation operator is not thread-safe, #386).
- Comment style: state constraints the code can't show; no change-narration.

## Task 1 — TDD: physical-conservation tests (red), then the fix (green), then adapt existing tests

**Step 1 (red).** Two new tests; run both BEFORE the fix and record the failing
magnitudes in the report (they must fail for the *predicted* reason):

1a. In `species/dncdm_inv_test.cpp`, new section "physical conservation (#385)":
  - `BuildComposite` via the factory with `dncdm1.inverse_decays=yes`,
    `dncdm1.quantum_statistics=yes`, plus `fc.set("dncdm1.deg", "0.015")`
    (a shooting-realistic amplitude; `deg` is a real input key,
    ncdm_base_species.cpp:66).
  - Seed random positive f's in the y-slots exactly like the existing
    "Conservation identities" section; `a = 1e-3`;
    `comp->ApplyKernelBackgroundDerivs(a, y, dy)`.
  - Compute per-species raw grid moments of dy (parent: bg grid `GetQ()/dq()`,
    `ε = sqrt(q² + a²·GetMass()²)`; daughters: `q_bg()/dq_bg()`, ε = q):
    `dN_i = Σ dq q² dy_i`, `dE_i = Σ dq q² ε dy_i`.
  - Assert the PHYSICAL identities using only public per-species accessors
    (`factor()`, `GetT()`), no κ and no (2π)³ anywhere in the test:
    - energy: `factor_H·dE_H + factor_l·dE_l + factor_φ·dE_φ = 0` (rel 1e-10,
      scale-relative like the existing `Scale3` pattern)
    - number vs fermion: `(factor_H/T_H)·dN_H + (factor_l/T_l)·dN_l = 0`
    - number vs boson: `(factor_H/T_H)·dN_H + (factor_φ/T_φ)·dN_φ = 0`
    - non-triviality: `|dN_H| > 0`.
  - `GetT()` does not exist yet: add `double GetT() const { return T_; }` next to
    `GetDeg()` in `species/ncdm_base_species.h`.

1b. In `species/dr_psd_test.cpp`, new section "absolute normalization":
  - Standalone fermion `DrPsdSpecies` with `dr_T = pow(4./11., 1./3.)`, thermal
    seed (set `dr_T` so `has_initial_abundance_` triggers), fine grid
    (e.g. `dr_q_min=1e-3, dr_q_max=50, dr_N_q=400`), `a = 1`:
    `Rho(bg)` must equal `rho_nu_rel()` — one standard massless neutrino —
    to quadrature tolerance (use rtol 1e-3; tighten if observed convergence
    allows). Add accessor `double rho_nu_rel() const { return rho_nu_rel_; }` in
    `ncdm_base_species.h`. Pre-fix this fails by ≈ (2π)³/2 ≈ 124×.

**Step 2 (green).** The fix:

2a. `species/dr_psd_species.cpp:37`:
    `SetDegAndFactor((stat_ == Statistics::Fermion ? 2.0 : 1.0) / pow(2 * _PI_, 3));`
    with a comment: bare-occupation storage ⇒ factor carries g/(2π)³, so deg means
    physical dof (g=2 fermion = particle+antiparticle, matching a deg=1 standard
    NCDM neutrino; g=1 boson).

2b. `species/dr_psd_species.{h,cpp}`: add `void PinTemperature(double T)` that sets
    `T_ = T` and re-runs `SetDegAndFactor(GetDeg())`. Call it from the
    **`DNCDMInvSpecies` constructor** (not Create — direct-built test composites
    must get it too) for both daughters with `parent_->GetT()`, before anything
    reads their `factor_`.

2c. `species/dncdm_inv_species.{h,cpp}` κ boundary. Single source of truth:
    `double KappaStoredToBare() const { return parent_->GetDeg() * pow(2 * _PI_, 3) / 2.; }`
    (public — the tests use it). Read `GetDeg()` live at every call (shooting
    updates deg between iterations).
  - `ApplyKernelBackgroundDerivs`: new ctor-sized scratch `fH_bare_` (parent
    q_size); fill `fH_bare_[i] = κ·y[parent slot i]`; pass to the kernel instead
    of the raw slot pointer; write back `dy[parent slot i] = df_H_[i] / κ`;
    daughter legs unchanged.
  - `ConservationMoments`: scale the parent's dy leg by κ (scratch copy) before
    handing to `kernel_->ComputeMoments`, so the reported moments are all-bare and
    the existing bare-identity test stays valid as written. Document in a comment.
  - `ApplyKernelPerturbDerivs` line ~254: `fH_gather_[i] = κ · std::exp(pvecback[pf + i]);`
    — everything downstream (F_H = f̄·Ψ, the /f̄ write-back, the floor guard) is
    already expressed in `fH_gather_` and needs no further change.

2d. Guards in `DNCDMInvSpecies::Create` (species/dncdm_inv_species.cpp):
  - parent uses a PSD file → `class_test_severe` (structural: key presence).
    Detect via the parent (add a minimal `bool UsesPsdFile() const` accessor to
    `ncdm_base_species.h` if none exists — check for one first).
  - parent `ksi ≠ 0` → `class_test` (numeric, sampler-varyable): the stored
    FD(q−ξ)+FD(q+ξ) sum has no single bare-occupation reading. Needs a
    `GetKsi()` accessor (add like `GetT()`).
  - `dncdm1.dr_T` explicitly set and ≠ parent `GetT()` (beyond 1e-12 rel) →
    `class_test` explaining T is pinned to the parent's for kernel-unit
    consistency.

**Step 3 (adapt + green).** Existing-test updates, then the full suite:

3a. `species/dncdm_inv_test.cpp` perturbation-conservation section (~line 537):
    the reconstruction `fi = std::exp(bg[...])` must become
    `fi = comp->KappaStoredToBare() * std::exp(bg[...])` (the kernel's cached
    `df_bg_H` is now bare-space). The decay-only cancellation section and the
    background bare-identity section need **no** change (κ cancels structurally /
    ConservationMoments converts) — if they fail, that's a bug in your fix, not
    the test.
3b. `species/dr_psd_test.cpp`: rescale any absolute-density expectations by
    1/(2π)³ (grep for factor/Rho assertions; re-derive rather than fudge).
3c. Build all five targets + run the FULL ctest suite (not just the five):
    `test-decay-kernel test-dr-psd test-dncdm-switch-copy test-dncdm-inv
    test-species-types` must pass; report any other target that breaks with output.
    (`test-decay-kernel` must pass UNTOUCHED — it pins the frozen kernel.)

**Step 4.** Run both new tests post-fix (green), confirm the 1a/1b failures are
gone with the same seeds. Commit everything as ONE commit with explicit paths
(source + tests + this plan file), message style:
`dncdm_inv: reconcile PSD conventions at the kernel boundary (#385)` with a body
explaining the occupation-suppression convention, ending with the repo's standard
Co-Authored-By/Claude-Session trailer block (ask the controller if you don't have
the session URL).

## Task 2 — Physics validation reruns (after Task 1 merges to the branch)

All via `python/tests/dncdm_inv_ab.py` (read its argparse/docstring first),
`OMP_NUM_THREADS=1`, scratch under the session scratchpad:

1. `--decay-only` and `--pert-regression`: must remain byte-identical to
   origin/master (composite never built on that path — any diff means Task 1
   leaked outside the composite).
2. `--tiny-gamma`: sane Γ→0 limit.
3. `--lasing`: the q_min sweep. Expect: sector energy conservation now holds in
   PHYSICAL (factor-weighted) units, not just bare moments; report drift numbers.
   If CPU budget allows, add the fixed-log-step control `dr_q_min=1e-5,
   dr_N_q=166` (the handoff's unfinished confound-disentangling run, >30 min).
4. `--anisotropic-stress`: report the envelope ratios; note the parent Ψ_H
   blow-up at a≈3.2e-5 (#386) is NOT fixed by this task — record where it occurs
   now, don't chase it.
5. Report all numbers verbatim; flag anything that contradicts the expectations
   above rather than rationalizing it.
