# DNCDM: the minimal knob set for the three run types

**For the cleanup/simplification branch.** Written 2026-08-17, at the end of the
`hpc_cmb` / `hpc_defl` / `hpc_reduced` / `hpc_nodefl` / `hpc_lna` / `hpc_freeze`
campaign. Its job is to say which parameters have to survive a cleanup and which
were introduced to vary something that has since been fixed.

**Status (added 2026-08-22): executed.** The cleanup commit `9f369204` ("dncdm:
remove the diagnostics and collapse the knobs to the shipped set") applied §2's
dispositions. **Read §2 as a record of what was done, not as pending work** — it is
written in the imperative and no longer means it. Of the 31 knobs §2 names, 15 are
gone as input keys and 16 remain, and the deletions follow the recommendation almost
row for row. §4 was carried out in full and then some: `deflate_collision` is gone
*and so is the machinery* — zero occurrences of "deflat" anywhere under `species/`,
`include/`, `source/`, `tools/`.

Five departures, each checked against the code on 2026-08-22:

| knob | §2 said | what happened |
|---|---|---|
| `dr_rate_cap` | remove outright | **kept.** Tagged **[?]** ("I never touched it") and it is load-bearing at `species/dncdm_proxy_species.cpp:413`. The tag did its job |
| `fluid_approximation` | "the guard can go with the key" | **kept** — it is a shared CLASS input key, not DNCDM's to delete |
| `momenta_bins` | fix at 16, delete | **kept**, consistent with this note's own ⚠ in §5 (not adequate for a4Pi) |
| `quadrature_strategy` | fix at 3, delete | **kept**; reason not recorded |
| `dr_rta_exchange` | *"do not delete on my authority"* | **deleted**, in the same commit that spared `dr_rate_cap` — differentiated treatment, so it reads as a considered call rather than a sweep, but it was not this note's to authorise |

Everything in §2's "Keep as a knob" table survived unchanged.

---

**Evidence tags.** Every disposition below carries one:

| tag | meaning |
|---|---|
| **[M]** | measured on the cluster in this campaign, number quoted |
| **[C]** | stated by a code comment, a test assertion, or a shipped default |
| **[P]** | from an earlier campaign, recorded but NOT re-measured here |
| **[?]** | I have no evidence either way — **do not delete on my say-so** |

Anything marked **[?]** should be resolved by whoever does the cleanup, not
assumed. I was not working on the RTA/proxy path, so most of it is [P] or [?].

---

## 1. The three run types

`dr_representation` selects them (`species/dncdm_dr_species.cpp`). It defaults to
`psd` when `inverse_decays = yes` and `integrated` otherwise, so type 1 and 2 can
omit it. **Do not remove `dr_representation`** — it is the dispatch. [C]

### Type 1 — the full method (`psd`, resolved daughter PSD, no reduction)

```ini
dncdm1.type              = ncdm_decay_dr
dncdm1.m                 = 0.3            # eV
dncdm1.T                 = 0.71611
dncdm1.deg               = 1
dncdm1.Gamma             = 1e8            # 1/Mpc, code units
dncdm1.inverse_decays    = yes
dncdm1.quantum_statistics = yes
dncdm1.dr_f_ini_l        = 1
dncdm1.dr_f_ini_phi      = 0
dncdm1.momenta_bins      = 16
dncdm1.quadrature_strategy = 3
dncdm1.lump_limiter      = yes
dncdm1.dr_q_min          = 1e-2
dncdm1.dr_q_max          = 50
dncdm1.dr_N_q            = 149            # per the grid rule, see §3
dncdm1.dr_bg_refine      = 1
evolver_background       = 3
evolver_thermodynamics   = 1
evolver_perturbations    = 3
tol_perturb_integration  = 3e-8
back_integration_stepsize = 7e-3          # per Gamma, see §3
l_max_ncdm               = 17
```

Verified against `hpc_nodefl/out/G1e8_m0.3_m6t_nd_used.ini` — every key above is one
the campaign actually set, except `back_integration_stepsize`, which is new here and
appears in `hpc_lna`'s inis. `N_ur = 1.0176` is also required: CLASS does not subtract
a species when an NCDM is added, so leaving the default gives N_eff ~ 5.1.

### Type 2 — the moment method (`psd` + Galerkin reduction)

Type 1 **plus**:

```ini
dncdm1.dr_reduced_moments      = 6
dncdm1.dr_reduced_weight       = entropy
dncdm1.dr_reduced_table        = yes
dncdm1.dr_reduced_table_stride = 1
```

and `deflate_collision` **absent** (its default is already off). See §4.

### Type 3 — the cheap RTA method (`proxy`)

```ini
dncdm1.type              = ncdm_decay_dr
dncdm1.dr_representation = proxy
dncdm1.m, .T, .deg, .Gamma  as above
dncdm1.dr_f_ini_l        = 1
dncdm1.dr_f_ini_phi      = 0
dncdm1.dr_N_q            = <coarse; the background cost lever>   [P]
dncdm1.dr_q_min / dr_q_max                                        [?]
```

`inverse_decays`, `quantum_statistics` and `lump_limiter` all default to **`true`**
in `DNCDMProxySpecies::Create` (lines 321–323), where the `psd` path defaults
`inverse_decays` and `quantum_statistics` to **`false`**. [C] **That asymmetry is a
trap** — a cleanup that unifies the defaults will silently change one of the two
paths. Decide deliberately and write the choice down.

---

## 2. Disposition table

### Keep as a knob — these select the method or the convergence rung

| knob | value | why it stays |
|---|---|---|
| `dr_representation` | `psd` / `proxy` | the dispatch between the three types [C] |
| `dr_reduced_moments` | `0` (off) or `6` | selects type 1 vs type 2, **and** it is the reduction's own convergence axis. 4/6/8 give genuinely different answers (1.3–1.7e-4 in C_ℓ^TT at Γ=10⁸) so the rung is real, not decorative [M] |
| `dr_N_q` (+ `dr_q_min`, `dr_q_max`) | see §3 | the quadrature/grid convergence axis. ⚠ **it means different things per type** — see §5 |
| `dr_bg_refine` | `1` | see §5; keep only if removing it costs a lot of code |
| `Gamma`, `m`, `T`, `deg` | — | physics |

### Fix at a value and delete the knob

| knob | fix at | evidence |
|---|---|---|
| `dr_reduced_weight` | `entropy` | the only weight that is a congruence — it is what makes `P` and `R` adjoint, hence `max Re λ(M̃) ≤ max Re λ(M)`. `flat` diverges outright in a real run (ΔP/P = 5.9e4) and `occupation` is 3–4 orders worse on the dissipativity test. It is an input only so that prediction can be re-checked [C][P] |
| `dr_reduced_table` | `yes` | 11–26× faster than matrix-free and the accuracy difference is entirely the ln-a interpolation, which §3 fixes. Matrix-free is now a *diagnostic*, not an option [M] |
| `dr_reduced_table_stride` | `1` | stride 4 costs 4.1e-3 → 5.4e-1 in C_ℓ^TT across Γ=10⁸→10¹⁰ and **saves no wall time at all** (13.16 vs 13.67 core-h). It is a pure memory knob with a large accuracy price [M] |
| `lump_limiter` | `yes` | ⚠ its shipped default is **`false`** (`decay_transition_kernel.h:392`) while every production run in this campaign sets `yes`. Fixing it means changing the default, not just deleting the key [C][M] |
| `quadrature_strategy` | `3` (trapezoid) | Gauss-Laguerre measured 2× cost for nothing [P] |
| `momenta_bins` | `16` | ψ 4.8e-4, δ_cdm 1.0e-5 against a finer reference; 1.85× cheaper than 32. ⚠ **but see §5** — it is *not* adequate for a4Pi [P] |
| `stratified_quadrature` | `yes` (already the default) | cures a ~100%-per-bin sampled-emission error (wiggle 1.315 → 0.073 at 1.4× cost) [P] |
| `n_gauss`, `n_sub` | delete | **inert** under `stratified_quadrature`: the emission nodes become the daughter cell edges, so `dr_N_q` sets the quadrature. The code says so at the parse site [C][P] |
| `deflate_collision` | `no` | §4 |

### Remove outright — diagnostics that must never ship

| knob | why |
|---|---|
| `lumped_loss` | its own comment says "**Do not ship a run that sets it**". `no` restores the energy-exact two-bin split, which is *not* positivity-preserving. It exists only to size what lumping costs [C] |
| `balanced_gather`, `balanced_gather_jacobian` | off by default; the jacobian variant is "correct but numerically fragile" [C], and `balanced_gather` was recorded as breaking `CollisionDiagonal` [P] |
| `deposit_order` | fix at `2`. The order-4 B-spline was tried against the q-oscillation problem and does not work [P] |
| `dr_write_psd` | debug output [C] |
| `dr_bg_threads` | measured **not** to pay — the background is serial and accuracy-limited, not throughput-limited [P] |
| `dr_rate_cap` | [?] — I never touched it |
| `fluid_approximation` | already a hard error in both `psd` and `proxy` (`class_test_severe` on it being *set at all*). The guard can go with the key [C] |

### RTA / proxy knobs — **[P] and [?], do not delete on my authority**

`dr_rta_C3` (0.2802), `dr_rta_C5` (2.5228), `dr_rta_n3` (0.8714), `dr_rta_vshut`
(1.7995), `dr_rta_exchange` (1.0), `dr_rta_form` (`powers`|`paper`),
`dr_q_sampling` (`log`).

These are a **fitted calibration, not a derivation** — the transport rate
Γ_T = (ρ_H/ρ_sec)(1−1/γ²)^q[C₃ε^n₃γ⁻³ + C₅γ⁻⁵] was fitted to ±40% [P]. A later
result found the exponent `p` runs 2.98→5.03 tracking ε_H rather than Γ, and that
coarse grids fake a γ⁻⁵ tail [P]. **So the constants may not be final**, and
hard-coding them could freeze in a calibration that is still moving. Ask whoever
owns the proxy work before collapsing them.

---

## 3. Parameter values by Γ — the campaign's output

Recommended settings for types 1 and 2. `back_integration_stepsize` is the one
this campaign changed most and it is **not currently in any harness config**.

| Γ | m = 0.3 | | m = 0.06 | |
|---|---|---|---|---|
| | `back_integration_stepsize` | `dr_N_q` | `back_integration_stepsize` | `dr_N_q` |
| 10⁸ | 7e-3 (default) | 149 | 3.5e-3 | 297 |
| 10⁹ | 1.75e-3 | 441 | 1.75e-3 | 437 |
| 10¹⁰ | 8.75e-4 | 657 | 4.4e-4 | 657 |
| 10¹¹ | 2.2e-4 | 1169 | — | — |

Two laws behind it, both measured [M]:

* **ln-a interpolation error ∝ Γ·h²**, so `back_integration_stepsize ∝ Γ^-0.5`.
  `ReducedOperatorTable::Apply` blends two rows *linearly*; the ladder measured
  order 1.98–2.08 over four cells. Refining is **free** — b8/b1 = 0.97–1.08× across
  seven cells, because `Apply` touches two rows whatever the resolution.
* **`dr_N_q`** follows the existing rule
  (`max(28 bins/dec at Γ=10⁸ × Γ^0.25, 27 bins/dec at Γ=10⁷ × Γ^0.17) × mass
  factor`, over 3.699 decades), doubled for m=0.3 above Γ=10⁸.

⚠ **Every Γ=10¹¹ result produced before 2026-08-17 carries a ~44% low-ℓ error in
C_ℓ^TT** from the default ln-a grid [M]. Any golden file or reference spectrum at
that Γ is void and must be regenerated at `back_integration_stepsize = 2.2e-4`.

There is a **~5e-4 C_ℓ^TT floor at m=0.3** once the stepsize is ≤ 8.75e-4,
Γ-independent, absent at m=0.06 [M]. Ruled out: the table, the deflation, the
reduced basis (4/6/8 moments give it to 0.13%), the background sector density, and
generic CLASS. Not chased further because it sits ~15× below cosmic variance and at
the level of the reduction's own closure systematic. **It is not a reason to keep
any knob.**

---

## 4. The deflation — recommend removing it entirely

`deflate_collision` is a shift-invert inverse iteration on a descending σ ladder
with a positivity acceptance test, a two-sided rank-one update, and a
congruence-transformed version of that update for the reduced block. It is a large
code footprint.

What switching it **off** costs, against the reduction's own ~5e-4 closure
systematic [M]:

| | Γ=10⁸ | 10⁹ | 10¹⁰ | 10¹¹ |
|---|---|---|---|---|
| m=0.3 | 5.6e-4 | 1.0e-3 | 2.1e-3 | 4.3e-3 |
| m=0.06 | 5.7e-5 | 1.5e-4 | — | — |

It saves **42–65% of every reduced cell** (it is O(N³) in the *unreduced* state
size and the reduction does not touch it). The effect is confined to a few low ℓ —
medians are 1e-8 to 1e-5 — and even 4.3e-3 is ~100× below cosmic variance there.

Two further measurements that bear on removing the code rather than the knob [M]:

* it is **orthogonal to the table**: the tabulation drift is 2.60e-3 deflated and
  2.61e-3 undeflated; the deflation effect is 9.98e-4 tabulated and 9.88e-4
  matrix-free. No cross-term, so nothing else depends on it;
* it barely touches the operator the evolver actually integrates. On real
  background rows it cuts the reduced block's worst spectral abscissa by 1.48× at
  Γ=10⁸ and only **1.045×** at 10⁹ — its grip *weakens* as Γ rises, because most of
  what is positive in M̃ is the physical non-equilibrium positivity a linearised
  collision operator has away from equilibrium, not the spurious mode.

**Recommendation: delete the deflation and quote the residual as a known
systematic.** If a more conservative line is wanted, keep it behind the flag for
Γ ≥ 10¹⁰ only — but that keeps all the code for one corner of parameter space.

**Keep `CLASS_DNCDM_REDUCED_ABSCISSA`** (the env-var diagnostic in
`BuildReducedTable`) even if the deflation goes. It is ~60 lines, costs one dense
28×28 eigensolve per tabulated row, and it is the only instrument that reports the
abscissa of the block the evolver integrates rather than of the full-space operator.
It is what showed the deflation was correcting the wrong thing [M].

---

## 5. Traps

**`dr_bg_refine` means different things per run type.** In the exact scheme it
subsamples the *perturbation state* under a fine background — a state-size change.
In the reduced tabulated scheme the state is moments, so it moves the *quadrature
grid only* and is nearly free. A cleanup that documents it once, for one type, will
be wrong for the other. [M]

**`momenta_bins = 16` is right for C_ℓ and P(k) and wrong for a4Pi.** Every cell in
this campaign ran 16, and a4Pi differences between schemes were 10–100× larger than
C_ℓ differences. If a4Pi is ever a deliverable it needs `momenta_bins = 32` and
`l_max_ncdm = 25`, and the a4Pi numbers in `docs/reports/` bound *scheme agreement*,
not accuracy. Do not fix `momenta_bins` in a way that makes 32 unreachable. [P][M]

**`tol_perturb_integration = 3e-8`, not the 3e-6 of the older matrices.** At 3e-6 a
28-variable and a 418-variable system take different step sequences at the same
*relative* tolerance and that artefact swamps everything being compared. It is a
requirement of *comparing* type 1 against type 2; a production type-2 run alone may
be affordable looser, which was never measured. [M]

**The proxy path's defaults differ from the psd path's** for `inverse_decays`,
`quantum_statistics` and `lump_limiter` (see §1). [C]

**`l_max_ncdm = 17` bounds where a4Pi is physical** at kτ ≲ 30; beyond that it is
hierarchy-truncation reflection. This does not affect C_ℓ or P(k). [P]

---

## 6. What this note does not cover

* **The RTA/proxy calibration constants** — [P]/[?] throughout, see §2.
* **Anything about `ncdm_decay_dr`'s `integrated` representation**, which none of
  these three run types uses but which `dr_representation` still dispatches to.
* **Whether the grid rule itself should change.** The campaign found it too shallow
  (it pays Γ^0.167 below 10¹⁰ crossing to Γ^0.25, while holding the spurious
  eigenvalue fixed needs ≥Γ^0.27), and §3 patches it with per-Γ overrides rather
  than re-deriving it. Fixing the rule properly is a separate decision. [M]
* **Background-module knobs other than `back_integration_stepsize`.**

## Provenance

`~/dncdm-harness/xreport/report.pdf` is the cross-matrix analysis through
2026-08-17; the matrices are `hpc_cmb`, `hpc_defl`, `hpc_reduced` (the three it
covers) and `hpc_nodefl`, `hpc_lna`, `hpc_freeze` (the follow-ups, which produced
§3 and §4). Each folder's `README.md` states what it measured and how to read it;
`collect.py --status/--cl/--ladder/--abscissa/--verdict` regenerate the numbers.
