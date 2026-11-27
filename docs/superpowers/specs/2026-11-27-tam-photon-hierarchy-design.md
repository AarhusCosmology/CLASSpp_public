# Total angular momentum (TAM) hierarchy for photons

**Status (2026-11-27):** built as designed on branch `tam-photon-hierarchy` (on master after
#502), in four commits:
1. the tensor gravitational-wave fix;
2. the abstract base plus the optimal leaf (byte-identical);
3. a fix found on the way: the scalar tight-coupling hand-over lacked s₂;
4. `PhotonsTamSpecies`.

Where the build departs from the text below:
- **Dot spelling.** It is `<name>.hierarchy` with `<name>.type = photons`, which is how
  single-instance dot syntax works. "`photons.hierarchy`" below is wrong.
- **No separate input check for TAM with vectors.** `perturb_init` stops on vector modes for
  every basis. The TAM vector hooks throw, in case that changes.
- **E₃ (my call).** The TAM scalar hand-over sets E₃ at second order, which class_public does
  not. The optimal basis hands over the corresponding G₁ and G₃, so both bases now hand over
  the same state. `test-scalar-tca` checks both against their own equations.
- **The +2 offset is kept, as in class_public.** The reason is that Θ₂ … Θ_{L+2} has as many
  multipoles as F₀ … F_L, so each hierarchy is as long as the optimal one and B is the extra
  one. It is not precision parity: that would be +1 (see Measured). "Switches basis without
  retuning" under Multipole counts overstates it.
- **CI regressions** are C++ tests (`test-tam-hierarchy`, `test-photons-tam`,
  `test-scalar-tca`, `test-tensor-tight-coupling`), not `TestReviewRegressions`.
- **s₂ fix (commit 3).** The optimal scalar hand-over (and its output columns) lacked the s₂
  factors. The hierarchy's P then jumped by 12 % at the switch for K/k² = 0.1. Spectra move by
  ≤ 0.007 % even at |Ω_k| = 0.1 **[M]**.

Measured, against the verification list below **[M]**:
- **Optimal path.** The refactor and TAM leave it byte-identical (FAST_MATH=OFF, one build
  directory, 28 output files).
- **Flat space.** With long hierarchies TAM equals optimal to ≤ 7e-7.
  - Scalars at defaults: ≤ 0.013 %.
  - Tensors, truncation alone (tight trigger, `benchmarks/tam/lmax_scan.py`): TAM at internal
    l_max L + 1 matches optimal at L, so the +2 offset puts TAM between optimal at L + 1 and
    L + 2. Tensor EE at l ≤ 30 with input 5: TAM 0.6–0.8 %, optimal 2.3–2.9 %. Variables:
    3(L − 1) against 2(L + 1).
  - At the default trigger TAM looks closer still to converged (0.16–0.41 % against
    1.9–2.4 %), but only because its truncation error partly cancels the trigger's.
- **Curved space** (converged runs):
  - Scalar EE at l = 2–10: 0.63 |Ω_k| %.
  - Tensor EE and BB: about 50 |Ω_k| % (4.94 % and 5.12 % at Ω_k = 0.1).
  - TT: ≤ 2e-4.
- **class_public v3.4.3**, its own TAM/optimal ratio, same settings: agrees to ≤ 0.009
  percentage points.
  - class_public v3.4.1–v3.4.3 fail on every tensor run in a closed universe
    (`index_start_spline outside of range`, harmonic module), with either basis. That failure
    is theirs, not ours.
- **Cost.** One thread, AC power: 0.498 s against 0.495 s (lensed scalars), 0.696 s against
  0.694 s (with tensors).
- Harnesses: `benchmarks/tam/` (`verify.py`, `cross.py`, `timing.py`, `tca_trigger.py`,
  `ab/`).

**Credit.** The TAM hierarchy in CLASS is Cyril Pitrou's work: he wrote it for class_public
around 2019–2020, and it was released in class_public v3.4.0 (2026). The equations, the coefficients,
the closure and the tight-coupling relations ported here are his, as is the photon
tight-coupling term in the tensor gravitational-wave source. CLASS++ contributes the
species structure that holds them. The code, the commits and the PR credit him the same
way.

Evidence tags: **[M]** measured here, **[C]** the code says so (CLASS++ or class_public
v3.4.x), **[P]** stated in a paper or in class_public's documentation, not re-measured,
**[?]** no evidence yet.

## Problem

class_public v3.4 added the TAM hierarchy for photons and ultra-relativistic species
(class_public v3.4.0; equations from Hu & White
astro-ph/9702170 for flat space, Hu, Seljak, White & Zaldarriaga astro-ph/9709066 for
curved space, and the closure of Pitrou, Pereira & Lesgourgues 2005.12119). One input key,
`hierarchy = optimal | tam`, chooses the basis.

The optimal hierarchy (1305.3261) is exact in flat space, and for temperature and for
unpolarised radiation also in curved space. For polarisation in curved space it is an
approximation, with errors at large angular scales of order |Ω_k| per cent for scalars and
50 |Ω_k| per cent for tensors **[P: class_public explanatory.ini, 2005.12119]**. TAM is exact
there. CLASS++ should offer it under the same key, so that curved-space polarisation can be
computed exactly and the two codes can be compared directly.

## Scope

**In:** photons, scalar and tensor modes. Also one tensor fix from the same release that is independent
of the basis (see "Tensor tight coupling feeds the gravitational-wave source").

**Out:**
- **Vector modes.** class_public v3.4.0 also reworks them in both bases (signs, vector Bessel
  functions, the `icv_iso`/`icv_oct` initial conditions, the primordial vector spectrum,
  the `wTk` output). That is a separate PR. Until then, `hierarchy = tam` with vector modes
  requested is an input error.
- **ur, dr, idr, ncdm, and the module-owned tensor ur hierarchy.** They are unpolarised, so
  the optimal hierarchy is already exact for them in curved space. TAM would be a change of
  basis with nothing gained and, for tensors, two more multipoles. This departs from
  class_public, where `hierarchy = tam` also switches ur. The comparison with class_public
  below has to allow for that.

## Approach: two sibling photon classes under an abstract base

`PhotonsSpecies` becomes an abstract base with two concrete leaves:
`PhotonsOptimalSpecies` (today's photon code, moved) and `PhotonsTamSpecies` (new file).
`PhotonsSpecies::CreateAll` reads `hierarchy` and builds exactly one of them. The two
therefore can never both be present, with no check needed.

**What the base keeps.** class_public's TAM keeps δ_γ, θ_γ and σ_γ in the same
normalisation as the optimal basis **[C]**. Everything that reads only those is therefore
shared and stays in the base:
- the background;
- initial conditions (photons start tightly coupled and set only δ_γ, θ_γ **[C]**) and the
  synchronous-to-Newtonian shift;
- scalar tight coupling and radiation streaming, the baryon equation, perturbed
  recombination, the Einstein equations;
- the polarisation source P while tight coupling is on (5 s₂ σ_γ^tca / 8; P is a physical
  quantity and does not depend on the basis);
- `StressEnergy`, which reads only δ_γ, θ_γ, σ_γ **[C]**;
- the transfer outputs of δ_γ and θ_γ.

The common part of the layout (`idx_delta`, `idx_theta`, `idx_shear`) lives in a base
layout struct. Each leaf derives its own layout from it.

**What each leaf must provide.** The base re-declares every basis-dependent `BaseSpecies`
hook as pure (`override = 0`) and adds the photon-only hooks the module needs. If a leaf
forgets one, the code does not compile. With "TAM inherits from optimal", the same mistake
would silently run optimal code on a TAM layout. The basis-dependent methods are:
- index registration and derivatives, scalar and tensor;
- `MaxMultipole`: TAM tensors need l_max + 2 (below);
- `CopyPerturbationsAcrossSwitch`, `MarkUsedInSources`, `MarkTensorUsedInSources`,
  `AddPerturbationOutput`;
- the vector hooks. `PhotonsOptimalSpecies` keeps today's vector code unchanged.
  `PhotonsTamSpecies` throws; the input check makes that unreachable.
- the new photon-only hooks, which replace the module sites listed next.

**Module sites that move into the photons.** In about nine places, the perturbations module
reads photon multipoles beyond δ, θ, σ. These become calls to photon methods:
- `perturb_vector_init`, when tight coupling ends:
  - scalars: fill l = 3 and the polarisation multipoles;
  - tensors: the hierarchy from `TensorTightCoupling`, including the drift term.
- `perturb_sources`: the polarisation source P^(0) with tight coupling off, and the tensor
  P^(2).
- The tensor gravitational-wave source: the photon anisotropic stress.
- The tensor rows of `perturb_print_variables`.

Photons already have a sanctioned exception to "no species picking in modules". This change
shrinks that exception: afterwards the module no longer knows which basis the photons use.

The vector-mode sites (vector sources, vector baryon coupling) still read the optimal
layout. They stay in place until the vector PR, and are reachable only with
`PhotonsOptimalSpecies`.

**TAM coefficients.** TAM streams with the k-dependent coefficients ₛκ_l^m
(class_public's `twokappam`, `zerokappam`, `q_m`). They depend on k, K and m.
- For m = 0, temperature needs no table: q ₀κ_l^0 = k l s_l, so it uses the shared `s_l`
  **[C]**.
- The other coefficients live in `PhotonsTamSpecies`'s `PerturbScratch`. Scratch is
  private to one k-mode integration and survives approximation switches. None of this goes
  into the shared workspace.

**Input.**
- `hierarchy = optimal | tam`, default `optimal`. class_public's spellings (`tam`, `TAM`)
  are accepted; any other value is an error.
- `photons.hierarchy` comes for free from the single-instance dot-notation table.
- `hierarchy = tam` together with vector modes is an error that names the key.
- `explanatory.ini` gets class_public's entry, with one added line: in CLASS++ only the
  photons switch basis, and why.

**Multipole counts.** `l_max_g_ten` and `l_max_pol_g_ten` keep meaning "the same expansion
order" in both bases: TAM tensors add 2 internally, as class_public does **[C]**. Scalars are
unchanged. A user then switches basis without retuning, and runs compare like for like with
class_public.

At defaults the tensor photon hierarchy grows from 12 variables (F₀…F₅, G₀…G₅) to 18
(Θ₂…Θ₇, E₂…E₇, B₂…B₇). The scalar one shrinks by 2: E₂…E₁₀ instead of G₀…G₁₀ **[C]**.

## Tensor tight coupling feeds the gravitational-wave source

During tensor tight coupling, the photon anisotropic stress currently contributes nothing to
the gravitational-wave source (`perturb_einstein`, tensor branch) **[C]**. class_public adds
8/5 a² ρ_γ Θ₂ with Θ₂ = −4/3 h′/κ′. In the optimal basis, the same stress is
−4√6 a² ρ_γ F₀/15 with F₀ = √6 · 4/3 h′/κ′. Both give −32/15 a² ρ_γ h′/κ′ **[C, by hand]**.

The fix takes the photon state from `TensorTightCoupling` with drift = 0. The drift needs
h″, and h″ is what the gravitational-wave source is computing. It goes in its own commit,
before the refactor, so its effect is measured alone. Like #502, it changes every tensor
lvl2 scenario by design.

**Measured 2026-11-27 [M]** (defaults against the tight-coupling trigger 0.002, and against
#502's tight reference):
- With the term, h no longer depends on when tight coupling ends. TT at every l, and EE and
  BB above l = 100, now agree with the tighter trigger to 4e-5; without it they differ by
  2e-3.
- The default spectra move by up to 0.2 % (TT and EE at l = 101–300).
- The term does **not** remove the −0.3 % in EE and BB at l = 30–40. EE and BB below
  l = 100 still depend on the trigger, by 3e-3 to 5e-3, through the polarisation source of
  modes that are still tightly coupled during recombination. That is a separate question,
  outside this spec.

## Approaches rejected

- **A `hierarchy` flag inside one Photons class** (class_public's way: 69
  `switch (ppt->hierarchy)` sites **[C]**). Every photon method would carry two interleaved
  versions. The default path would get harder to read for a feature that matters only when
  Ω_k ≠ 0.
- **`PhotonsTamSpecies` derived from today's `PhotonsSpecies`.** Any method TAM does not
  override would run the optimal code on a TAM layout. That code reads `idx_pol0` and
  `idx_pol1`, which TAM leaves at −1, so it would read `y[-1]` with no error. The abstract
  base turns "remember to override" into a compile error.
- **TAM only, dropping optimal.** In flat space it gains nothing. The tensor hierarchy is
  1.5 times larger, and every lvl2 reference would move.
- **TAM for ur and the other massless species.** Nothing gained (see Scope).

## Verification

1. **Refactor alone** (abstract base plus `PhotonsOptimalSpecies`, module sites moved): pure
   code motion.
   - Run a HEAD-vs-master A/B in one build directory with FAST_MATH=OFF. Expect identical
     output, and investigate any relative difference above 1e-10.
   - Committed goldens are stale under `-ffast-math`, so they are not the test.
2. **Flat space, scalars.**
   - Under Θ_l = (2l+1)/4 F_l, the TAM temperature equations and closure *are* the optimal
     ones **[C, by hand]**. Only polarisation differs: other variables, and a closure
     coefficient of (l+3) instead of (l+1).
   - TAM-vs-optimal C_l differences at defaults must be no larger than the optimal
     polarisation truncation error. Measure that error by raising `l_max_pol_g`.
   - Both bases must converge to the same answer as `l_max_pol_g` grows.
3. **Flat space, tensors.** With long hierarchies (`l_max_g_ten` = `l_max_pol_g_ten` = 20),
   TAM and optimal agree to 1e-4. At defaults, report the difference of each from that
   converged answer: TT, EE and BB, in l bands.
4. **Curved space.**
   - At Ω_k = ±0.01 and ±0.1, report the TAM-vs-optimal difference in TE, EE and BB at low
     l, and compare with the claimed |Ω_k| % and 50 |Ω_k| % **[P]**.
   - Cross-code check against class_public's TAM: same Ω_k, matched precision (see the
     class_public build recipe), and long ur hierarchies in both. The long ur hierarchies
     remove class_public's ur-in-TAM truncation difference.
   - Neither code is ground truth: a converged run (long hierarchies) settles any
     disagreement.
5. **Unit tests.**
   - The coefficients: the m = 0 identity q ₀κ_l^0 = k l s_l, the closed forms for
     ₂κ_l^m and ₀κ_l^2, and the K = 0 limit.
   - The TAM tensor tight-coupling state gives the same P^(2) as the optimal one, to
     rounding, drift term included.
   - The input errors: TAM with vectors, an unknown value, and the `photons.hierarchy`
     spelling.
6. **Cost.** Time TAM against optimal at defaults with scalars, and with scalars plus
   tensors, recording the power source.
7. **CI.**
   - PR CI compares every `CLASS_INPUT` scenario with a classyref built from master, and
     master rejects the unknown key `hierarchy`. So no `CLASS_INPUT` entry in this PR.
   - TAM runs (flat and curved, scalars and tensors) go into `TestReviewRegressions` instead,
     with numbers from items 2–4.
   - The existing lvl2 scenarios test the refactor and the tensor fix: scalar scenarios
     unchanged, tensor scenarios different by design.

## Constraints

- **Tensor accuracy baseline.** #502 (tensor initial conditions to (kτ)², Pitrou's
  correction of class_public v3.3.0) is on master, so tensor C_l no longer depend on the
  starting time. The remaining tensor errors at defaults **[P: #502]**:
  - EE and BB at l = 30–40: −0.3 %, from the tensor tight-coupling trigger;
  - EE at l = 10–20: +1.9 %, from `l_max_g_ten` and `l_max_pol_g_ten`;
  - TT at l = 2: −0.75 %, from the time and reionisation sampling.

  The tight-coupling fix may move the first; TAM's longer tensor hierarchy may move the
  second. Report both against #502's tight reference.

## Left to the builder

- **Second-order E₃ at the scalar tight-coupling exit.** Whether to fill it. Optimal fills
  G₁ and G₃ at second order; class_public's TAM fills only E₂ **[C]**. Add it if a scan of
  the tight-coupling trigger shows TAM more sensitive to the switch time than optimal.
- **The TAM tensor tight-coupling state.** Map `TensorTightCoupling`'s (F₀, G₀) onto
  (Θ₂, E₂) through the flat-space basis relation: at leading order Θ₂ = −F₀/√6 and
  E₂ = −G₀, which reproduces the optimal P^(2) exactly **[C, by hand]**. B₂ is zero at this
  order in both bases. Verification item 5 checks the mapping.
- **When the TAM coefficients are filled.** For example on first use for a given k, cached
  in scratch.
- **How the module reaches the photon hooks.** `photons()` returning `PhotonsSpecies&`, or
  one cast.
- **Output column names** for E_l and B_l in the perturbation output.
- Cite this spec from the declaration of `PhotonsSpecies`.
