# ETDRK4 on the perturbations: extending the diagonal-exponential evolver past the background

Date: 2026-08-06
Branch: `dncdm-inverse-qs`
Follows: `04110024` (ETDRK4, background only)

## Problem

`04110024` put a fourth-order exponential integrator (Cox–Matthews ETDRK4, step
doubling) on the background and removed the collision stiffness there: 8.6x at
Gamma=1e8, with the gain *growing* in Gamma, which is the signature of removing a
stiffness proportional to Gamma rather than winning a constant factor.

The background is 0.4% of runtime. The perturbations are the other 99.6%, and they
carry the same disease:

| Gamma | 1e5 | 3e5 | 1e6 | 3e6 | 1e7 |
|---|---|---|---|---|---|
| perturbations [s] | 9.7 | 13.6 | 24.7 | 51.4 | 155 |

Local log-slope 1.0 — cost is a straight line in Gamma. The step-limiting variable,
decoded with `CLASS_PERT_INDEX_MAP=1`, is overwhelmingly the *parent's own hierarchy*
(`parent bin0 l15`, `bin15 l15`, ...), i.e. the same relaxation rate `lambda = K/eps1
-> a*Gamma` that capped the background.

`perturbations_module.cpp:2296` currently refuses `evolver=etd` and falls back to
rkdp45 with a warning, because no species reports a Jacobian diagonal on the
perturbation path. That hook is the missing piece.

## Why the background result should transfer

Two structural facts, both already in the tree:

1. **`ApplyPerturbationOperator` is the exact Jacobian of the background collision
   network.** That was the point of the `440ec131` lumping. Its diagonal is therefore
   the same object `CollisionDiagonal` already computes, per (species, bin).

2. **Since `#386`, all three legs carry unnormalized `F = f_bar * Psi = delta_f`.**
   The parent no longer evolves the normalized `Psi_H`, so there is no
   `-(f_bar_dot/f_bar) Psi` dilution term to cancel against the collision diagonal —
   it cancels structurally. The parent's perturbation damping is the physical
   `-K/eps <= a*Gamma` directly, identical to the background `diag_H`.

The diagonal is expected to be l-independent on every leg: the parent's loss is a
same-bin local term with no angular factor, and the daughters' loss was lumped so it
couples a bin only to itself. **This is a claim to be measured, not assumed** — see
Verification (a).

## Non-goals

- Photon Thomson damping is also a stiff diagonal. Tight-coupling already handles it;
  it is not this session's problem.
- No re-derivation of the evolver's order. `04110024` verified order 4 by local error
  at fixed step (31.8x per halving, against 32x for h^5). Not re-litigated.
- No new embedded pair. Cox–Matthews has none; `04110024` chose step doubling
  deliberately, because inventing one risks controlling the wrong order silently.

## Design

### 1. Put the diagonal callback on the shared evolver signature

All four evolvers gain a trailing parameter:

```cpp
void (*derivs_diagonal)(double x, double* y, double* diag, void* parameters_and_workspace)
```

`evolver_etd` already has it. `evolver_rk`, `evolver_rkdp45` and `evolver_ndf15` accept
and ignore it — which is idiomatic in this codebase, not a wart: `evolver_rkdp45`
already documents `minimum_variation`, `evaluate_timescale` and
`timestep_over_timescale` as parts of the shared CLASS evolver signature that it does
not use.

Consequences:

- `background_module.cpp:702`'s direct `evolver_etd(...)` call collapses back into the
  `generic_evolver` pattern. The special case introduced to keep `04110024`'s blast
  radius at one site disappears rather than multiplying.
- The `evolver_type::etd` warn-and-fall-back block at `perturbations_module.cpp:2296`
  is deleted.
- Three modules (background, thermodynamics, perturbations), one dispatch pattern,
  zero special cases.
- It leaves the door open for `ndf15` to use an analytic diagonal as a Jacobian seed
  later, without another signature change.

Drive-by corrections: the `evolver_type` enum comment (`precision.h:22`) and the
`evolver_etd.h` header block both still describe the abandoned order-2 ETD2RK scheme.

### 2. Species hook

```cpp
/** Contribute this species' Jacobian DIAGONAL for the scalar perturbation ODE.
 *  Default: no-op (zero diagonal), at which the phi functions reduce ETDRK4 to its
 *  explicit counterpart. */
virtual void PerturbDerivsDiagonal(const PerturbLayout& layout,
                                   double tau,
                                   const double* y,
                                   double* diag,
                                   const perturb_parameters_and_workspace& ppaw) const {}
```

Mirrors `PerturbDerivs` (`base_species.h:349`) in shape and constness, and
`BackgroundDerivsDiagonal` (`base_species.h:223`) in contract. `CompositeSpecies`
loops its children, as it already does for the background diagonal
(`composite_species.cpp:54`).

Only `DNCDMInvSpecies` overrides it. Every other species reports zero, so the standard
sector is integrated exactly as before.

### 3. The `DNCDMInvSpecies` override

One kernel pass per step:

1. Rebuild `fH_gather` from `pvecback` (the coarse daughter-PSD gather).
2. `PrepareTransitions(a, m, Gamma, fH_gather, f_l, f_phi, a_prime_over_a)`.
3. `CollisionDiagonal(diag_H, diag_l, diag_phi)`.
4. Write each bin's value into **every l slot** of that bin:
   `diag[layout.index_per_q[i] + l] += diag_*[i]` for `l = 0..L`.

The parent needs no kappa conversion: `d(dy_i)/d(y_i)` is invariant under the
stored/bare rescaling, since the write-back divides the rate by the same kappa the
read multiplies the occupation by. This is the same argument
`BackgroundDerivsDiagonal` makes at `dncdm_inv_species.cpp:462`.

Two rules that are load-bearing:

- **Buffers live in the per-workspace `Scratch`** (`dncdm_inv_species.h:156`), never in
  members. The background override uses member `diag_H_` / `df_H_` and is safe because
  the background is single-threaded; the perturbation path is not. `CreatePerturbScratch`
  exists for exactly this class of bug.
- **Re-prepare; do not read the kernel cache.** The evolver calls derivs before the
  diagonal at the same state, so the cache would usually be right. `04110024` rejected
  that shortcut on the background for a reason worth repeating: a diagonal silently
  taken from a *different* state than the RHS it is paired with is the failure mode
  that stays stable while converging to the wrong attractor.

### 4. Module-level callback

The RHS uses a two-part pattern: a static thunk `perturb_derivs`
(`perturbations_module.cpp:94`) that casts `parameters_and_workspace` and forwards to
the member `perturb_derivs_member` (`:5525`). The diagonal mirrors it exactly —
`perturb_derivs_diagonal` static thunk beside `:94`, forwarding to
`perturb_derivs_diagonal_member`, which:

1. Zeroes `diag` over `pv->pt_size`.
2. Interpolates the background into `ppw->pvecback` at `tau`, with the *same* flags the
   RHS uses at `:5554` — `pba->normal_info`, `pba->inter_closeby`,
   `&ppw->last_index_back`. Sharing `last_index_back` with the RHS is intended: the
   two calls are at the same `tau`, so the closeby hint is warm.
3. Loops `species->PerturbDerivsDiagonal(*layout, tau, y, diag, *pppaw)` over the same
   species/layout list the RHS uses at `:5761`.

## Verification

Ordered. Each gates the next.

**(a) The diagonal is what we say it is.** Unit test: extract the true diagonal of
`ApplyPerturbationOperator` per l by unit vectors, and pin it against
`CollisionDiagonal`. `DumpCollisionSpectrum` (`dncdm_inv_species.cpp:663`) already
performs exactly this column extraction, so the machinery exists and only needs
lifting into a test. This converts the l-independence claim from an assumption into a
measurement, which is the specific caution carried over from the previous session.

Follows the precedent of `decay_kernel_test.cpp`, which pins `CollisionDiagonal`
against a finite-difference Jacobian at 3e-9.

**(b) The physics is unchanged.** P(k) and C_l, ETD vs rkdp45, matched settings,
~0.1% tolerance with C_l^TE zero-crossings handled — never a blind max-rel-diff.

**Magnitude before agreement.** Check `P[0] ~ 51.8` first, always. The recorded method
lesson is that two routes agreeing is worthless if both carry the same spurious growing
mode.

**Grid stays converged**: `dr_N_q >= 100` at Gamma=1e6, `>= 200` at Gamma=1e7. An
entire round of perturbation timings was voided once by measuring on `dr_N_q = 20`,
which is below the stability threshold.

**(c) The headline.** Perturbation cost vs Gamma, both evolvers, converged grid. The
falsifiable claim is that **the ~Gamma slope breaks** — not any single speed-up number.
A constant-factor win at fixed Gamma would not be this scheme's signature and would
suggest the diagonal is not what is limiting the step.

## Risks, carried explicitly

**Dense output may eat the win.** ETD's dense output is cubic Hermite, one order below
the step, so a step much larger than the sampling spacing makes the *table* less
accurate than the trajectory. On the background this was minor (capping at 1x moved
`rho_dr_phi` 9.7e-4 -> 6.1e-5 for 2.4x the cost, so it ships off by default). In the
perturbations the sampled sources *are* the deliverable — every C_l and P(k) comes
through them — so the same effect is far more consequential. This may force
`CLASS_ETD_HMAX_OVER_DX` on by default for the perturbation path. Measured, not
assumed.

**Per-step cost is ~1.8x DP45's.** Step doubling costs ~11 RHS per attempted step
against DP45's 6, so the step must grow by more than 1.8x merely to break even. The
background cleared this comfortably at high Gamma; at low Gamma the crossover may sit
above the Gamma range of interest, and ETD may simply lose there. That is an acceptable
outcome to measure, not a failure — `evolver` is a per-module setting.

**Diagonal dominance is a property of the trajectory, not the operator.** At a=0.3 the
off-diagonal sits 5-9 orders below the diagonal only because the emission band has
swept to q ~ a*m/2 ~ 267 while the daughters still occupy q ~ 1-6. With uniformly
populated daughters the ratio returns to ~1 at every a*m from 1 to 1e4. Do not assert
dominance in any unit test. (The scheme does not require it — an exponential method
integrates the diagonal exactly and treats the rest explicitly regardless — but any
claim about *why* it wins does.)

**Order discipline.** The ETD2RK detour cost a full build-and-measure cycle to learn
that an order-2 exponential method is worse than explicit here (10-25x slower, with
`errmax` pinned at the controller's own fixed point, i.e. accuracy-limited everywhere).
Nothing below order 4 goes in, and the constant-N exactness check cannot verify
coefficients on its own — it only pins the sum `c1 + 2*c2 + c4 = phi1`.

## Success criteria

1. The per-l diagonal test passes at finite-difference tolerance.
2. `evolver_perturbations=etd` produces P(k) and C_l agreeing with rkdp45 to ~0.1% on a
   converged grid, at sane magnitude.
3. A measured cost-vs-Gamma curve for the perturbations under both evolvers, with the
   slope stated. Either the ~Gamma law breaks, or we have a measured reason why it does
   not — both are publishable outcomes for the branch.
4. No change to any non-dncdm result: the standard sector under `evolver=etd` sees a
   zero diagonal by construction.
