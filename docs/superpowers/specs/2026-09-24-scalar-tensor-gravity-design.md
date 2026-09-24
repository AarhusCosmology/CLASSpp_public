# Modified gravity: a Jordan-frame scalar-tensor species behind a gravity seam

**Status (2026-09-24).** Shipped on branch `scalar-tensor-gravity`: the seam
(`species/modified_gravity.h`), `ScalarTensorSpecies`, and the source-sampling hook of
section 4. Section 4 was found during verification and is not part of the original
design. The closure unknown shipped as `Omega_Lambda_st`, in units of 3H0², not as
`V_Lambda_st` in 1/Mpc²; section 2 says why. Measured results are in *Verification*.
Everything above it is the design as reasoned, tagged where it was checked.

## Problem

The H0 Olympics (arXiv:2107.10291) ranks eighteen models. After #440-#445, CLASS++
runs seventeen of them. The one left is Early Modified Gravity (EMG, sec. 2.4.6;
Braglia, Ballardini, Finelli & Koyama, arXiv:2011.12934). EMG is a scalar σ with a
non-minimal coupling to the Ricci scalar,

    S = ∫ d^4x √-g [ (M_pl² + ξσ²) R/2 − (∂σ)²/2 − Λ − λσ⁴/4 ] + S_m ,

which the paper ran in a modified hi_class. Modified gravity does not fit the plugin
architecture. A species is a sector that couples to the rest only through gravity. EMG
changes gravity itself: the Friedmann equation and the Einstein equations.

The capability should be useful beyond EMG, leave a small footprint in the modules, and
build nothing that is not needed yet.

## What people test, and what this covers

Modified-gravity work with Boltzmann codes falls into four families:

1. **Non-minimally coupled scalars, F(φ)R** (CLASSig, hi_class). These are the early-MG
   H0 models: EMG; the massless non-minimally coupled scalar (Braglia et al. 2004.11161);
   induced gravity, which is Brans–Dicke with ω = 1/(4ξ) (Umiltà, Ballardini, Rossi et
   al.). Late-time ones are extended quintessence (Planck 2015 XIV) and non-minimally
   coupled thawing quintessence (the recent DESI phantom-crossing papers).
2. **The μ–Σ parametrisation** (MGCAMB). This is the quasi-static late-time
   phenomenology of DES, KiDS and the DESI full-shape papers.
3. **Horndeski α-basis** (hi_class, EFTCAMB): α_M, α_B, α_K, α_T.
4. **f(R) and nDGP**, mostly constrained on non-linear scales.

**This work implements family 1 only.** It is covariant, a single code path covers early
and late models, and it contains EMG. It is Horndeski's G4(φ) subclass: α_T = 0 and
α_B = −α_M. Families 2 and 3 are deliberately not built (see *Not built*). They would
attach to the same seam without touching the modules again.

## Approach

### 1. The seam: `ModifiedGravity`

CLASS assumes GR in exactly two places **[C]**:

- the Friedmann equations in `BackgroundModule::background_functions`, which says so in
  a comment ("this is the only place where the Friedmann equation is assumed");
- the scalar Einstein equations in `PerturbationsModule::perturb_einstein`.

A new interface, `species/modified_gravity.h`, has one method for each place where GR is
assumed:

    CloseFriedmann(a, K, pvecback_B, pvecback, rho_tot, p_tot, rho_r, rho_m)
    HPrime, EtaPrime, HPrimePrime, AlphaPrime   // synchronous-gauge scalar equations

A species that modifies gravity implements it next to `BaseSpecies`. Each module
resolves the role once with a `dynamic_cast` to the interface, as it already does for
the PPF fluid, the other module-owned closure of the Einstein equations. A second
modified-gravity species is an error. Modules ask for the role, never the type.

The footprint in the modules:

- **background:** one call after the species loop, before the unchanged
  `H = sqrt(rho_tot − K/a²)`;
- **perturbations:** four ternaries in the synchronous block of `perturb_einstein`, each
  choosing between the hook and GR's line, which stays as it is. The tca and rsa
  corrections interleaved with those lines are untouched. Each hook reads the module's
  running totals at the point where the GR line would.

With no modified-gravity species present, both modules take the GR branch. ΛCDM is
unchanged up to code-generation noise: under -ffast-math, a branch around GR's lines
changes how they are compiled. That is accepted rather than engineered away; the code
stays readable.

### 2. The background: an effective fluid that makes GR's lines exact

Use units with M_pl = 1 and CLASS densities ρ ≡ 8πGρ_phys/3, and write f(φ) ≡ F/M_pl².
The Jordan-frame equations are:

    f (H² + K/a²) + ḟ H = ρ_bare                          (00)
    −f (2Ḣ + 3H² + K/a²) = 3 p_bare + f̈ + 2Hḟ            (ij)

Here ρ_bare and p_bare are every other species plus the scalar's canonical
ρ_φ = (φ'²/2a² + V)/3. `CloseFriedmann` solves (00) for H with a cancellation-free form
of the quadratic root. It gets φ̈ from the Klein–Gordon equation in trace form (below),
which contains no Ḣ, then Ḣ from (ij). It then writes an **effective** density and
pressure into its own slots:

    ρ_eff = H² + K/a² − ρ_rest ,   p_eff = −(2Ḣ + 3H² + K/a²)/3 − p_rest ,

and adds them to the module's totals. With those totals, GR's two unchanged lines
reproduce H and H' = aḢ exactly (up to rounding). So `rho_crit`, `Omega_r`, the
thermodynamics and everything else downstream see a consistent GR-plus-dark-component
background. hi_class uses the same convention for `rho_smg` and `p_smg`.

**Energy buckets.** The species is `EnergyType::DarkEnergy`, so the generic bucket split
ignores it. `CloseFriedmann` instead rescales the radiation and matter buckets by
1/f − 1: each is a fraction of the Friedmann equation, and G_eff = G/f. Three things
follow:

- early times still pass the Ω_r ≈ 1 check;
- `thermodynamics_delta_neff_at_bbn`, which reads Ω_r · ρ_crit, picks up the BBN
  expansion-rate change as an equivalent ΔN_eff with no BBN code touched;
- z_eq stays G-independent.

For EMG today f = 1 to O(10⁻⁷), so nothing reported today changes.

**Initial conditions.** The perturbation IC code counts only Radiation, Matter and Other
species. It therefore skips this species, and its fractions (f_ν, f_b, …) are true ratios
of physical densities, which is what the G-independent leading-order adiabatic mode
needs. The next-to-leading `om` term assumes GR's a(τ). For f_ini ≠ 1 its coefficient is
off by 1/√f, on a term that is already O(ωτ) at τ_ini **[C]**.

**Closure.** The Lambda species keeps closing the GR budget. The scalar's potential
carries a constant V_Λ = 3H0² Ω_Λ,st, which the existing shooting machinery tunes so that
ρ_eff = 0 today, i.e. H(z=0) = H0. The species reports `GetOmega0() = 0`, and the residual
responds to Ω_Λ,st as 1/f₀. For EMG the untuned residual is 5e-6 **[M]**, inside the
shooter's tolerance, so shooting costs one extra background solve.

The unknown must be O(1). `fzero_Newton` stops on an absolute step of 1e-3. With the
unknown in 1/Mpc² (~1e-7), it declared convergence after a single step. That was right
for EMG but left H(0) 34% off for the massless coupling at ξ = 0.3, with exit status 0
**[M]**.

### 3. The perturbations (synchronous gauge)

Define δf = f_φ δφ, δf' = f_φ δφ' + f_φφ φ' δφ, and δf'' likewise. Start from
G_μν = [T_μν + ∇_μ∇_ν f − g_μν □f]/f, perturbed in Ma & Bertschinger's synchronous gauge
(derivation in the implementation comments). The CLASS-unit totals δρ, (ρ+p)θ, δp and
(ρ+p)σ include the scalar's canonical stress-energy. With ℋ = aH:

    h'  = [k²η + (1.5a²δρ − ½(3ℋδf' + k²δf + 3ℋ²δf))/f] / (½ℋ + ¼f'/f)
    η'  = [1.5a²(ρ+p)θ + ½k²(δf' − ℋδf)] / (f k²)
    h'' = −2ℋh' + 2k²η − [9a²δp + 3δf'' + 3ℋδf' + 2k²δf + f'h' + 3(2ℋ' + ℋ²)δf]/f
    α'  = −2ℋα + η − [4.5(a²/k²)(ρ+p)σ + δf + f'α]/f

Each reduces to CLASS's GR line at f = 1, δf = 0. The last one gives the known
gravitational slip of scalar-tensor theories, ψ − φ = −δf/f in Newtonian gauge (using
δf_N = δf_S + αf').

**Klein–Gordon in trace form.** The field equation □φ = V_φ − f_φR/2 contains R, and
through R it contains h''. The trace of the field equations, R = (3□f − T)/f, removes it:

    C □φ = S ,   C = 1 + 3f_φ²/(2f) ,   S = V_φ − (f_φ/2f)(3 f_φφ (∂φ)² − T)

T is the total canonical trace, −δρ + 3δp for the perturbations. Perturbing gives

    δφ'' = −2ℋδφ' − k²δφ − ½h'φ' − (a²/C)(δS − C_φ δφ S/C)

It needs only h' and the module's totals, which are ready before species derivatives
run **[C]**. The background uses the same form, which is why it needs no Ḣ.

### 4. Source sampling for an oscillating field (found in verification)

EMG's field keeps oscillating at late times, 20 to 70 half-periods per unit redshift
below z = 10, with frequency ω ≈ 50–100 aH. That makes the metric, and so the late-ISW
source η' − ℋ'α − ℋα', ring at ω. `perturb_timesampling_for_sources` samples the
sources at the visibility and late-ISW rates, about 0.1/aH. It therefore catches the
ringing at random phases, and the transfer integral turns it into noise.

The first EMG run showed TT at ℓ = 30 at 2.9× ΛCDM, erratic over 6 ≤ ℓ ≤ 70, with EE
clean **[M]**. Only the late-ISW contribution was off, by up to 133×. Refining the source
step converged it: at stepsize 0.01 and 0.003 the same model gives 1.00–1.05× ΛCDM at low
ℓ **[M]**. The minimally coupled field (ξ = 0) aliases too, less.

Two hypotheses were ruled out, both **[M]**:
- δf terms: switching them off below z = 300 still left 37×.
- Parametric resonance of λφ⁴: δφ decays relative to δ_cdm at every k.

The fix is a generic hook, `BaseSpecies::SourceSamplingRate(pvecback)` (default 0). The
sampling takes the fastest of the default rates and the species' rate. The scalar-tensor
species reports ω/π, from the quartic oscillator's period at the amplitude its energy
implies, which gives twenty samples per period at the default step.

The sampling loop now reads the background in normal rather than short format, so that
species slots exist. With that, EMG at default precision matches the converged spectrum to
3.5e-3 in TT and 1e-4 in EE for ℓ ≤ 200 **[M]**. Ten samples per period left 1.2%. The
run takes 0.82 s against ΛCDM's 0.33 s, most of it the ndf15 fallback.

The existing axion scalar field would benefit from the same hook but does not use it yet.

### 5. The species: `ScalarTensorSpecies` (`scalar_tensor`)

The model is kept to what EMG and the external check (below) need:

    f(φ) = N_st + ξ_st φ² ,   V(φ) = 3H0² Ω_Λ,st + ¼ λ_st M_pl² φ⁴ ,   canonical kinetic term

φ is in reduced-Planck units. λ_st is the dimensionless physical coupling, and the
species converts it with M_pl in 1/Mpc taken from CLASS's constants. This covers:

- EMG: N = 1;
- the massless non-minimally coupled scalar: λ = 0;
- induced gravity / Brans–Dicke: N = 0, ω_BD = 1/(4ξ);
- Rock 'n' Roll early dark energy: ξ = 0.

The conformal coupling is ξ = −1/6 in this sign convention.

| key | meaning | default |
|---|---|---|
| `xi_st` | ξ; its presence builds the species | — |
| `phi_ini_st` | φ at a_ini (M_pl units) | required |
| `N_st` | bare coefficient of R | 1 |
| `phi_prime_ini_st` | φ' at a_ini | 0 (frozen) |
| `lambda_st` | quartic coupling λ; EMG's papers quote V₀ with λ = 10^(2V₀)/3.516e109 (see *Verification*) | 0 |
| `Omega_Lambda_st` | the shooting unknown for the closure, in units of 3H0²; not for users | — |

Background outputs: `phi_st`, `phi'_st`, `f_st`, `(.)rho_st` and `(.)p_st`, where the
last two are effective. Derived parameters, today, for a massless field:

- `ScalarTensor.G_eff_today`: G_eff/G = (2f + 4f_φ²)/(f(2f + 3f_φ²));
- `ScalarTensor.gamma_PPN_minus_1` = −f_φ²/(f + 2f_φ²).

These are the Solar System quantities the induced-gravity literature constrains. The
generic `ScalarTensor.f_peak` is EMG's peak energy injection, as Braglia et al. define it.

**Guards, in the factory, where the input keys are.** These fail rather than skip:

- Newtonian gauge;
- tensor or vector modes;
- non-adiabatic initial conditions;
- Ω_k ≠ 0;
- `BackgroundDensityOverH0Sq`, which throws: NEDE's pre-solve query assumes GR.

The perturbations module refuses two more things alongside any modified-gravity species:
- a PPF fluid, in `ResolveSpecies`: PPF closes the Einstein equations too and reads
  `p_tot'`;
- the N-body-gauge sources, where they are switched on: they take δρ_com from GR's
  Poisson equation and c_a² from `p_tot'`. Two inputs reach them, `Nbody gauge transfer
  functions` and `extra metric transfer functions`; the second route was found in review.

Both modules resolve the species through one helper, `FindModifiedGravity`, which also
refuses a second modified-gravity species as a severe error.

`PPrime` returns the canonical scalar's p' only. The effective p' would need Ḧ, and its
only consumers are the two refused above.

## Alternatives rejected

- **A hi_class-style α-basis engine** (G2–G5 or α_i functions, a V_X field, stability
  checks). It covers everything, but the footprint is large, EMG would still need a
  covariant background solver, and hi_class's initial-condition machinery is a known pain
  point. Most of that generality is not needed.
- **Effective-fluid perturbations**, i.e. the scalar returning the δρ_eff that makes GR's
  equations exact. δρ_eff depends on h', the quantity being solved for, so it cannot be a
  stress-energy contribution without a module-side linear solve. Rewriting the totals in
  place would also feed dressed totals to the trace-form Klein–Gordon equation, which
  needs canonical ones.
- **Replacing the whole synchronous block** under a gravity branch. The tca and rsa
  corrections are interleaved with GR's four lines, so they would be duplicated in the
  species. Per-equation hooks leave them alone.
- **Virtuals on `BaseSpecies`.** Five methods that one species implements would pollute
  every species. A separate interface is found once.
- **The scalar as the closure species** (Ω_Λ = 0, the scalar fills the budget). That needs
  a new `ClosureSpecies` value and input-module changes, and makes shooting start far from
  the root. Tuning a small correction beside Lambda needs neither.
- **An `EnergyType::Other` bucket split** for the effective fluid. The IC code would then
  count the effective radiation in f_ν's denominator and set a wrong leading-order
  adiabatic mode (f_ν too large by a factor f).

## Not built (YAGNI), and what each would cost

- **Normalising G today** (G_eff,Cavendish = G_N), which induced-gravity and
  massless-coupling papers impose by shooting φ_i or N. Cost: one more `ShootingTarget`,
  and a rule for which unknown it tunes. EMG needs none: σ → 0.
- **μ–Σ.** A species with no variables whose four hooks are MGCAMB's synchronous
  equations and whose `CloseFriedmann` leaves GR alone. Needs its own external reference.
- **Newtonian gauge, tensors (α_M friction: h'' + (2ℋ + f'/f)h' + k²h = source/f),
  curvature, a general Z(φ), other potentials.** Each is local to the species, except
  tensors, which would need a fifth hook.
- **MG-aware growth factor D, halofit/HMcode.** These remain GR-calibrated. So is
  `delta_tot`, which sums canonical perturbations.

## Verification plan

1. **ΛCDM unchanged** from master to the numerical floor (no modified-gravity species, so
   the GR branches): rtol 1e-4 on TT, EE and P(k), and TE measured against √(TT·EE).
2. **GR limit:** ξ = 0, λ = 0, φ constant. The species is inert and the spectra equal
   ΛCDM to the shooting tolerance. This exercises every seam with trivial physics.
3. **Bianchi identities.** The four scalar equations are redundant through ∇_μG^μν = 0.
   h' from the (00) constraint must have the derivative h'' from (ij), and α from the
   constraints must have the derivative α' from the traceless equation, along an actual
   EMG solution. This tests the derivation and the trace-form Klein–Gordon equation
   together, with no external code.
4. **External: hi_class `brans_dicke`** (public hi_class, built locally). It solves
   G4 = φ_BD/2, G2 = −V + ωX/φ_BD, which is this species at N = 0, ξ = 1/(4ω),
   φ_BD = ξφ², V = const. It is a different formulation: α-functions and a V_X
   variable, written by other people. Compare H(z), φ(z) and the lensed C_ℓ and P(k).
5. **EMG itself:** H(z=0) = H0 after shooting, ΔH0 at fixed 100θ_s for the paper's
   posterior region (ξ ≲ 0.3, σ_i ≈ 0.5, V₀ ≈ 2.5), and the guards firing.

## Verification

All **[M]**, measured on 2026-09-24 on this branch.

1. **ΛCDM against master:** TT 4e-6, EE 3e-6, P(k) 3e-5. This is code-generation drift
   from the ternaries, not a physics change.
2. **GR limit** (ξ = 0, λ = 0, φ_ini = 0.3), both codes on ndf15 at
   `tol_perturb_integration` 1e-8: ≤ 7e-7 in C_ℓ (TT, EE, TE, φφ), 6e-8 in P(k). At the
   default tolerance the two extra ODE variables alone move the spectra by 1e-4. Nor is
   the default evolver the same: a run with the species present falls back from rkdp45 to
   ndf15, because the species does not declare itself explicit-safe.
3. **Bianchi identities.** EMG at ξ = 0.3, σ_i = 0.49, λ = 3.42e-108, with k = 0.01,
   0.1, 0.5 and tolerance 1e-10. The quantity is |Δh' − ∫h''dτ| (and the same for α)
   over the integrator's own steps, relative to max|X|, in τ bins from 1 to 300 Mpc.
   - EMG: 1e-8 to 1e-5.
   - GR control with the same machinery: the same values, bin by bin. Its largest,
     5e-5 near τ = 200–250, comes from tight coupling.
   - Deliberately wrong equations: 3ℋ²δf → 2ℋ²δf in h' gives up to 2e-3; flipping the
     sign of C_φ in δφ'' gives 1.4e-3.
4. **hi_class `brans_dicke`** (commit ece81a57), ω_BD = 50, φ_BD,ini = 1, both codes
   at tolerance 1e-6 and source step 0.02:

   | quantity | agreement |
   |---|---|
   | f against M*² | ≤ 2e-7 at all z |
   | H(z) | ≤ 3.5e-5 |
   | P(k) at z = 0 | ~1e-4 (single interpolation glitches ≤ 2e-3 at k = 0.1 h/Mpc and at the edge) |
   | lensed TT, BD/ΛCDM ratio | 4e-4 rms, max 1.8e-3 |
   | EE | 1.1e-3 rms, max 6.4e-3 in the EE trough near ℓ = 500 |
   | φφ | 7e-4 |

   The BD effects themselves are 20–110%. At CLASS++'s default precision, TT at
   ℓ ≤ 100 is off by up to 7e-3. This is the default source step: 0.02 removes it, while
   a tighter tolerance does not. ΛCDM loses 1e-3 at ℓ ≤ 30 to the same step. The pins in
   `test_scalar_tensor_matches_hi_class_brans_dicke` avoid ℓ < 100 in TT.
5. **Braglia et al. 2011.12934**, the EMG paper, whose code is not public.
   - **Fig. 1** (σ_i = 0.54, V₀ = 2, ξ = 0.1/0.3/0.5, their θ_s and densities):

     | | theirs | here |
     |---|---|---|
     | peak Ω_σ | 0.117 / ≈0.24 / ≈0.375 | 0.117 / 0.238 / 0.376 |
     | early plateau | −0.031 / −0.117 / −0.23 | −0.031 / −0.117 / −0.228 |
     | σ_max | 0.56 / 0.61 / 0.65 | 0.559 / 0.610 / 0.649 |
     | peak redshifts | ≈2700–3100 | 2729–3183 |

     These match **only** with λ 10³ smaller than their printed λ = 10^(2V₀)/3.516e109.
     With λ as printed, the field rolls at z ≈ 2×10⁴. The H0 Olympics prints 3.156e109, a
     transposition of Braglia's 3.516e109. Whether the H0 Olympics chains used the figure's
     normalisation is **[?]**; `explanatory.ini` states both.
   - **Fig. 3** (TT against ΛCDM). Their ξ-dependence is reproduced at fixed h = 0.7069,
     the h of the ξ = 0 point, not at fixed θ_s: the figure varies ξ at fixed H0. For
     ξ = 0.1 the ratio matches theirs to ~0.003 at every ℓ read (10–2500). For ξ = 0.3 the
     peaks match within ~15%, and the troughs are shallower here. Troughs are the most
     sensitive to their unreported h and to the λ normalisation.
6. **EMG at fixed 100θ_s = 1.04092**, Planck 2018 densities, one massive neutrino,
   ξ = 0.3, σ_i = 0.49, V₀ = 2.54 (the paper's posterior mean):
   - ΛCDM: H0 = 67.04;
   - λ as printed: H0 = 71.84 (f_peak = 0.148 at z = 3.4×10⁴);
   - the figure's normalisation: H0 = 73.68 (f_peak = 0.173 at z = 5.7×10³).

   G_eff_today = 1 to 1e-9 in both.
7. **Aliasing** (section 4). EMG at default precision matches its converged spectrum
   (source step 0.003) to 3.5e-3 in TT and 1e-4 in EE for ℓ ≤ 200, down from 1.76. The
   converged answer is itself stable to three digits between steps 0.01 and 0.003, and
   between tolerances 1e-6 and 1e-8.
8. **Guards.** Newtonian gauge, tensors, curvature, PPF, NEDE and both N-body-gauge
   routes each raise their own severe error. Closure holds to 1e-3 for EMG and for the massless ξ = 0.3 case, also
   with 100θ_s shooting h.
