# Decay channel multiplicity: scenario B for the decaying-neutrino sector

**Status — 2026-09-11.** Design, written to scope the code change that lets the
`ncdm_decay_dr` sector represent more than one parent or more than one daughter mass
eigenstate, so that Chen, Oldengott, Pierobon & Wong's "scenario B" (arXiv:2203.09075,
Table 1) can be run *and validated against the exact q-resolved scheme*, rather than
asserted from their analytic factor of two.

**Implementation status — 2026-09-11, same day.** §4 SHIPPED in full, for BOTH
representations, and §5 step 1 passes:

* `DecayTransitionKernel::Config::n_parent` / `n_daughter`, with the three leg factors
  applied in the background sweep, both perturbation operators, `CollisionDiagonal` and
  `PositivityDecomposition`, plus the split/clamp energy bookkeeping.
* `dncdm<i>.dr_n_parent` / `dr_n_daughter` on **both** composites. **`dr_n_parent > 1`
  is REJECTED** — see the correction below; only B1 is runnable.
* The daughter's momentum-integration weight, by a different mechanism in each
  representation and idiomatic to each: the proxy's `DarkRadiationSpecies` has no `deg`,
  so the weight is applied in `DaughterRho` and the ρ source; `DrPsdSpecies` already has
  one, so the exact scheme uses `SetSpeciesMultiplicity` and its ρ, being computed from
  the PSD through `factor()`, needs no source scaling.
* The `nu_H` / `nu_l` / `nu_phi` exchange rates turned out to carry the same three leg
  factors, which is why `gross` is now split into a per-channel quantity. `eps_ne` is
  invariant under multiplicity, as a dimensionless departure from balance must be.
* **`dr_rate_cap` resolved**, and it was NOT a free choice. The cap documents itself as
  bounding the parent's per-bin rate in units of the expansion rate; capping `K` alone
  bounds the rate *per channel*, so `dr_rate_cap = 1e3` would permit a 2x stiffer system
  in scenario B than in scenario A and an A/B comparison would inherit the difference.
  The cap now divides by the parent's leg factor. Still one global scale on `K`, so the
  conservation identities are untouched.

Tests: `decay_kernel_test` gains `test_channel_multiplicity` at (1,1), (1,2), (2,1),
(2,3) and (1,2)+balanced_gather, plus `test_rate_cap_bounds_stiffness` (with the
disarmed-cap control, without which the claim would also pass if multiplicity did
nothing); `test_perturbation_jacobian` and `test_collision_diagonal` now also run at
non-unit multiplicity, so the operator is held to being the Jacobian of the background
RHS at every multiplicity. `dncdm_proxy_test` and `dncdm_inv_test` gain the weight and
guard checks. Full suite 43/43.

Measured end-to-end, not only in unit tests:

| check | proxy | exact |
|---|---|---|
| sector `a⁴ρ` drift across the decay epoch, (1,2) | 0.003% | 0.0028% |
| same at (1,1) — i.e. multiplicity adds no leak | 0.003% | 0.0028% |
| background at (1,1) vs the pre-change binary | **bit-identical** | — |
| lensed `C_l` at (1,1) vs pre-change | 1.7e-4% rms, 4.8e-4% max | — |

The `C_l` difference is reassociation in the perturbation operator amplified by the ODE —
the drift `decay_transition_kernel.cpp`'s own deposit-restructuring comment documents,
~600x inside the repo's 0.1% verification tolerance and ~1000x below the m = 0.06 eV
signal. Not chased to bit-equality, per that comment's advice.

### Correction, 2026-09-11 (PR #431 review)

**`deg` cannot carry a parent species count, and §4 item 3 was wrong to say it could.**
Caught in review; the error was mine, and the measurement that settles it is:

| config | sector `a⁴ρ` drift across the decay epoch | ρ_φ at a = 1e-4 |
|---|---|---|
| deg 1, n_p 1 | 0.0043% | 1.12 |
| deg 2, n_p 1 | 0.0207% | 8.34 (7.4×) |
| deg 2, n_p 2 | **3.64%** | 25.06 (22×) |

`deg` is **overloaded**. In the `Omega_dncdmdr` shooting path it is a PSD *amplitude* —
a genuinely diluted population, for which `f_bare = deg·f0` IS the physical occupation
and `KappaStoredToBare()` is right to carry `GetDeg()`. A species *count* is the
opposite: n identical species share ONE per-dof occupation. Setting `deg = n` to weight
the parent therefore also scales the occupation the kernel sees, and the collision is
nonlinear in it (Λ carries `f_l f_φ` and `f_H(f_l − f_φ)`), so the rate comes out far
more than n times too large. ρ_φ starts empty and is pure decay product: 7.4× at deg = 2
where a pure weight would give 2×.

The earlier observation that "deg = 2 doubles `rho_dncdm1` exactly" is true and was
misread — it says nothing about what the kernel boundary receives.

Supporting B2 properly means giving the parent a species count **separate from `deg`**,
multiplying `factor()` (hence ρ/n/Π) without touching the kernel boundary. That is a
change to `NCDMBaseSpecies`, deliberately not bundled here.

**B1 is unaffected and is supported**: `dr_n_daughter` is applied to ρ/n outside the
kernel and never reaches the boundary (measured drift 0.003%, i.e. the baseline). The
falsification test of §5 step 3 is a B1 test, so the campaign is not blocked.

The kernel's `n_parent` leg factor is kept, and tested: it is correct where it lives,
and it is what a proper parent species count will drive.

Still open, and it is a measurement rather than code: **§5 steps 2 and 3 have not been
run.** Nothing here says the A-calibrated `Gamma_T` formula times two reproduces the
measured scenario-B rate; that is what step 3 is for, and §6 says what is expected to
survive it and what is not.

Evidence tags below: **[M]** measured in this session, **[C]** stated by the code,
**[P]** earlier campaign, not re-measured, **[?]** no evidence — do not act on it
without checking.

---

## 1. Why

The MCMC campaign for the RTA representation reproduces Chen et al.'s **scenario A**
exactly as the code stands: one free-streaming spectator, one parent, one daughter,
`N_ur ≈ 1.0176` [C, the campaign inis]. Their **scenario B** — no spectator, all three
mass eigenstates in the interacting sector — is where their tightest lifetime bounds
live, and it is the more natural model in that nothing singles out one eigenstate for
exclusion from the coupling.

The cheap way to get scenario B is a multiplier on the RTA transport rate, matching
their `C = 2`. **That is rejected** (§6): it would leave scenario B with no exact
counterpart anywhere in this code, so scenario A would be anchored in the q-resolved
solve while scenario B was this campaign's fitted amplitudes carrying someone else's
analytic factor. The distinction is the whole value of the exercise.

Because the exact scheme (`DNCDMInvSpecies`) and the proxy (`DNCDMProxySpecies`) share
one `DecayTransitionKernel` [C, `dncdm_inv_species.cpp:382` and the proxy's ctor], a
change made in the kernel buys both an exact scenario-B reference and an RTA
scenario B from one implementation.

## 2. The counting

Let `n_H` be the number of parent mass eigenstates and `n_l` the number of daughter
mass eigenstates, all identical within each group. The physical cases:

| case | ordering | decay | `n_H` | `n_l` | spectator |
|---|---|---|---|---|---|
| A  | either | one channel | 1 | 1 | 1 (`N_ur`) |
| B1 | NO | ν₃ → ν₁, ν₂ | 1 | 2 | none |
| B2 | IO | ν₁, ν₂ → ν₃ | 2 | 1 | none |

Treating the members of each group as identical is consistent with the oscillation
data: |Δm²₃₁| and |Δm²₃₂| differ by Δm²₂₁/|Δm²₃₁| ≈ 3%, and Chen et al. explicitly do
not distinguish them (their Table 1 caption, "less than 5%").

The channel count is `n_H · n_l`. Every leg of the collision network then scales as:

* **parent legs** (forward loss *and* inverse-decay gain) × `n_l` — each parent dof has
  `n_l` channels open to it;
* **daughter legs** × `n_H` — each daughter dof is fed by, and consumed by, `n_H` parent
  species;
* **boson legs** × `n_H · n_l` — φ is one species and counts every channel.

Detailed balance survives this automatically, and that is the point of writing it per
leg rather than per species: forward and inverse on the same leg carry the same factor,
so the kernel's equilibrium condition is untouched at any `(n_H, n_l)`.

**Consistency check, derived independently.** Chen et al. §4.3 give the recipe as: for
NO, multiply by two the ν_H and φ collision integrals and all momentum-integrated ν_l
quantities; for IO, multiply by two the ν_l and φ collision integrals and all
momentum-integrated ν_H quantities. Both reduce to the table above — the collision-integral
factors are the per-leg rate scalings, the "momentum-integrated quantities" are the
integration weights of §4. The two derivations agree exactly.

## 3. What the code already gets right

The concern that motivated this investigation — degeneracy leaking into the `(1±f)`
Pauli/Bose blocking terms, where it would be wrong in a way indistinguishable from
physics — **is already prevented architecturally** and needs no new guard.

The kernel works entirely in **bare per-dof occupation** for all three species
[C, `decay_transition_kernel.h:422-425`]. The parent's stored PSD is deg-suppressed, and
`KappaStoredToBare()` converts it at the kernel boundary; the daughters already store
bare occupation, so their legs are untouched [C, `dncdm_inv_species.cpp:429-439`].
Degeneracy is therefore already a normalisation applied *outside* the kernel, which is
exactly where a channel multiplicity must not live and a momentum weight must.

This forces the design to keep **two concepts separate**, and they are separate today
only because both happen to equal one:

* **Integration weight** — `deg` → `factor()` [C, `ncdm_base_species.cpp:586`], and for
  the daughters `BareFactor() = parent_->factor() / KappaStoredToBare()`
  [C, `dncdm_proxy_species.h:310`]. The fermion daughter rides on it at an implicit
  weight of 1 and the boson at 0.5 — the latter being `g_φ/g_l` [C, `:305-309`].
* **Channel multiplicity** — currently hard-coded. The `-2.0` at
  `decay_transition_kernel.cpp:659` is `g_H/g_φ = 2`; the fermion leg carries no factor
  because `g_H/g_l = 1`.

Conflating them is the failure mode to design against: an ini that sets the parent's
`deg = 2` without the matching multiplicity would run, produce no error, and violate
energy conservation in the sector. That is the same class of silent trap as a
misspelled species type.

## 4. Change list

1. **`DecayTransitionKernel::Config`** gains `n_parent` and `n_daughter`, both defaulting
   to 1. The hard-coded `-2.0` becomes `-2.0 · n_parent · n_daughter`; parent legs scale
   by `n_daughter`; daughter legs by `n_parent`. One global scale per leg, never per bin
   — the same constraint the existing `max_rate` cap documents
   [C, `decay_transition_kernel.h:157-160`].
2. **A momentum-integration weight on the daughters**, which does not exist today:
   `DarkRadiationSpecies` and `DrPsdSpecies` carry no degeneracy, and
   `DaughterRho`/`DaughterNumber` hard-code the 1 and 0.5
   [C, `dncdm_proxy_species.cpp:720-738`].
3. **One input key pair on the composite** that sets the parent's `deg` and the daughter
   weight *together* with the kernel multiplicities, so an inconsistent sector cannot be
   spelled. An unread or half-set key here is silent, which is why this is one key pair
   and not four independent knobs.
4. **Kernel conservation tests** at `(1,2)` and `(2,1)`. Number and energy conservation
   on the grid must hold exactly as they do at `(1,1)`; those identities are the kernel's
   and are already pinned in `decay_kernel_test.cpp` [C].

`N_ur` is the ini's business, not the code's: scenario B sets it to zero.

## 5. Validation ladder

**Step 1 — local, minutes.** `(n_H, n_l) = (1,1)` through the new code path reproduces
current scenario-A results bit-for-bit; kernel conservation tests pass at `(1,2)` and
`(2,1)`. The second half is a real physics check, not plumbing, and it costs seconds.

**Step 2 — local, ~1 h.** Background-only cells across the whole Γ ladder for B1 and B2:
number and energy conservation, ρ_H/ρ_sec landing near 1/3 (B1) and 2/3 (B2), and the
ε_ne curve. ε_ne is the input the C₃ term reads and is the quantity expected to move,
since two daughters shift the detailed-balance point. Background-only cells are cheap
[M: 0.7 s at Γ=1e5 to 29 s at Γ=1e11, m=0.06, proxy grid defaults].

**Step 3 — cluster, ≈ 80–100 core-h.** The falsification test. `iso` cells at
m = 0.3 eV, Γ = 1e8 … 1e9.5, `x` plus refine twins, for B1; plus two or three scenario-A
cells on the new binary as a lineage control. Measure Γ_T as `proxynote/tools/isodata.py`
does and ask whether **2× the A-calibrated formula reproduces the measured B rate**.

B2 is deliberately not measured at this step. It shares the code path and the counting,
and its background is checked at step 2; if B1's ×2 holds, B2 is taken on the strength of
it. If B1's ×2 fails, B2 must be measured too, because the refit would not transfer
between two cases with different ρ_H/ρ_sec.

m = 0.3 eV and Γ ≥ 1e8 are not free choices: below that rung the ladder has no power to
distinguish closures at all, and at m = 0.06 eV nothing about the rate is measurable
[P, `rtanote/note.tex` §2].

Measured cost of the exact `iso` cells, from `hpc_prod2/logs/*.out` `ELAPSED_S`
[M], 8 cores, `n_k = 8`, `momenta_bins = 32`, m = 0.3 eV:

| Γ | wall | core-h | refine twin |
|---|---|---|---|
| 1e8 | 10 min | 1.3 | 0.6 |
| 1e8.5 | 28 min | 3.7 | 1.7 |
| 1e9 | 1.5 h | 12.0 | 5.1 |
| 1e9.5 | 4.7 h | 37.3 | 15.2 |
| 1e10 | 13.5 h | 108 | 44 |
| 1e11 | 166 h | 1329 | 465 |

Cost rises ×3 per half decade above 1e8. **`hpc_prod2/manifest.csv`'s own
`est_work_core_h` is ~2× low at the top** (637 estimated against 1329 measured at 1e11):
the cost model still has no background term, so do not size this round from it.

## 6. What is expected to survive, and what is not

**α_ℓ^eff should survive untouched, by argument rather than by hope.** If the daughters
are identical, their hierarchies are identical at every ℓ and every q for all time, so
the sector's Π is a reweighted sum over the same three fields and the collision
operator's ℓ-structure is unchanged. Only the weights and the overall amplitude move,
and both are already explicit in the rate formula.

**ρ_H/ρ_sec re-derives itself**, since it is read off this code's own background rather
than assumed (Chen et al. take it as 1/3 throughout).

**C₃, C₅ and n₃ are what the test is for.** They are fitted, and what they absorb is
precisely the gap between the analytic structure above and the exact kernel. Two
specific reasons the ×2 could fail: ε_ne traces a different curve when the
detailed-balance point moves, and n₃ is already the softest constant in the calibration
(estimators spread it over 0.28–0.51 [P, `proxynote`]). C₅ = 0.53 against COPW's 1 is
itself the standing evidence that "the derivation says so" has been 1.9× wrong before
in exactly this sector.

If the ×2 fails, the outcome is a refit of C₃/C₅ on scenario B, stated as such. That is
a result, and finding it from ten cells is the reason step 3 precedes any chain.

## 7. Rejected alternative

**`dr_rta_C = 2`, a multiplier on the RTA transport rate.** Costs nothing, reproduces
Chen et al.'s scenario B on their own terms, and cannot be validated: it changes only
the closure, so the exact scheme has nothing to say about it and there is no reference
to measure the closure against. It would put scenario A and scenario B on different
evidentiary footings within one paper, which is the specific failure this design exists
to avoid.
