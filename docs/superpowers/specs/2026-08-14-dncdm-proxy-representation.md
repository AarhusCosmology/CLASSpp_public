# `dr_representation = proxy` — a core-minute stand-in for the DNCDM inverse-decay sector

2026-08-14. `species/dncdm_proxy_species.{h,cpp}`, branch `dncdm-inverse-qs`.

## What it replaces, and what it keeps

The exact scheme (`dr_representation = psd`, `DNCDMInvSpecies`) carries a resolved
PSD for both daughters in the perturbations — 7524 of a run's 7812 variables — and
applies the linearised collision operator to every one of them at every *k*. At high
Γ with a light parent that is thousands of core-hours per cosmology.

The proxy keeps the background essentially exact and replaces only the
perturbations:

| | exact (`psd`) | proxy |
|---|---|---|
| parent background | own q-grid | **same** |
| daughter background | PSD, `dr_N_q` ~ 600 (grows as Γ^0.25) | PSD, `dr_N_q` ~ 96–224 |
| collision (background) | `DecayTransitionKernel` | **same kernel, unchanged** |
| daughter perturbations | one hierarchy per momentum bin | one integrated hierarchy each (`DarkRadiationSpecies`) |
| collision (perturbations) | linearised operator, per (k, q, ℓ) | relaxation-time closure |

Because the background runs the shipped kernel, number conservation, energy
conservation and detailed balance are *inherited*, not re-derived. The daughter
grid is the only background approximation, and it is a convergence knob:

| `dr_N_q` | 32 | 48 | 96 | 192 | 384 | exact (101) |
|---|---|---|---|---|---|---|
| sector E_comov | 1.3738 | 1.2800 | 1.2540 | 1.2513 | 1.2509 | **1.2536** |
| background wall | 1.6 s | 2.1 s | 3.8 s | 8.1 s | 18.6 s | — |

(Γ = 10⁷, m = 0.3 eV.) 96 reproduces the reference to 3×10⁻⁴. Default is 96.

**But `dr_N_q` is not the main grid lever — `dr_q_min` is.** Because the proxy needs
only the daughters' *moments*, q²f̄ is negligible over the bottom decade of a log
grid while the emission band still sweeps every one of those fine bins per node.
Raising `dr_q_min` is cheaper *and* more accurate in both directions tested
(Γ = 10⁹ background, error in the sector's comoving energy):

| `dr_q_min` | 1e-3 | 1e-2 | 5e-2 | 1e-1 |
|---|---|---|---|---|
| core-s | 62.2 | 28.9 | 13.5 | 8.2 |
| error | +0.41% | +0.17% | +0.09% | +0.04% |

Two parent-side knobs measured on the same ladder: `momenta_bins` = 12 beats 16 on
both cost and accuracy (21.4 s/+0.13% vs 28.5 s/+0.17%), and `quadrature_strategy`
= 3 is the most accurate of the four available. `dr_q_sampling = linear` is ~10×
cheaper per bin than log but converges more slowly — a near-tie at matched accuracy.
`dr_q_max` is weak and inverted (larger is cheaper *and* less accurate).

Together, at Γ = 10¹¹ (background, m = 0.3): `dr_N_q` 224 / `dr_q_min` 1e-1 /
`momenta_bins` 12 costs **28.5 core-s** against **431 core-s** for 448 / 1e-2 / 16,
at the same E_comov to 0.05% — and is self-converged to +0.06% under 2× daughter
refinement and −0.05% under 2× parent refinement.

## The perturbation closure

Two relaxation terms, applied to the parent's per-q hierarchy and to both
daughters' integrated hierarchies:

**ℓ = 0, 1 — exchange.** Each species' contrast relaxes toward the sector's common
contrast at ν_i = (gross transition rate)/(comoving number of species i). This is
the boosted decay rate for the parent and goes to zero for the daughters once the
parent is extinct. The target is the ν·B-weighted mean, which makes the exchange
**exactly conservative in δρ and (ρ+P)θ** — not optional, since the Einstein
equations see the sector's total.

**ℓ ≥ 2 — damping** at the transport rate

    Γ_T^phys/Γ = (ρ_H/ρ_sec) (1 − 1/γ²)^q [ β_ℓ C₃ ε^n₃ γ⁻³ + α_ℓ C₅ γ⁻⁵ ]

with α_ℓ = (3ℓ⁴+2ℓ³−11ℓ²+6ℓ)/32 (arXiv:2203.09075 eq. 16) and β_ℓ = ℓ(ℓ+1)/6, both
normalised to 1 at ℓ = 2. `ε` is the background's own departure from detailed
balance (net/gross parent loss), published as a background column; it falls from
0.72 to 10⁻⁴ across Γ = 10⁵ → 10¹⁰, which is what turns the effective exponent from
3 into 5.

**The coefficients are fitted, not derived**: C₃ = 0.2802, n₃ = 0.8714,
C₅ = 2.5228, q = 1.7995, from 72 points over 6 cells (Γ = 10⁵…10¹⁰, m = 0.3),
restricted to the campaign's own convergence gate γ ≤ 20. They reproduce the
measured Γ_T with a 16/84% spread of 0.74–1.46. Both correction factors earn their
place: dropping the velocity shut-off costs ~2× at recombination (which sits at
γ ≈ 2.1 for m = 0.3 — inside the turnover) and ~17× by a = 4×10⁻³ at Γ = 10¹⁰;
dropping the ε dependence widens the spread to 0.48–1.72.

arXiv:2203.09075's own (1/12)X⁵𝓕(X) is available as `dr_rta_form = copw` but is
not the default: over γ = 3 → 6 it falls 10.0× where the runs fall 16.9×, because
X ~ 1 there and 𝓕's own decline eats two of the five powers.

## Validation

Full lensed C_ℓ + P(k), m = 0.3 eV, rkdp45 at `tol_perturb_integration = 3e-8`,
against the exact scheme's **converged** rung (r1):

| ℓ band | Γ=10⁷ | Γ=10⁸ | Γ=10⁹ |
|---|---|---|---|
| 2–10 | −0.00% | +0.01% | +0.07% |
| 30–100 | +0.15% | +0.16% | +0.17% |
| 200–500 | −0.04% | +0.05% | +0.08% |
| 1200–2000 | −0.03% | −0.00% | +0.03% |

and the Γ = 10¹⁰/10⁹ acoustic signal at ℓ = 200–800 comes out **1.0074** against the
exact **1.0065**.

**Cost, after the background work.** Γ = 10⁹ full lensed C_ℓ + P(k): 88 core-s /
17 s wall on the tuned grid (was 151 / 70), same C_ℓ rms to the reference's own
error bar. The expensive corner, Γ = 10¹¹ with m = 0.06 eV: **151 core-s / 48 s
wall**, agreeing to ≤0.08% in every band with the 2099 core-s run of the same cell.

⚠ **Compare against r1, not r2.** The campaign's own convergence pair disagrees at
low ℓ, and the gap grows with Γ: r2/r1 at ℓ = 2–10 is −0.91% / −1.03% / −2.58% for
Γ = 10⁷/10⁸/10⁹. At Γ = 10¹⁰ only r2 exists, so the proxy's apparent +5.4% there is
consistent with the reference's own convergence error and is **not established
either way**.

## Evolvers

Measured as delivered accuracy vs cost (equal `rtol` is not equal accuracy),
truth = ndf15 @ 3×10⁻¹⁰:

- **Perturbations: rkdp45 at `tol_perturb_integration = 3e-8`.** 2×10⁻⁴ rms for
  72 core-s at Γ = 10⁷. etd is dominated everywhere (~2× the cost at equal
  accuracy). ndf15 is 3–7× dearer but far more accurate, and its cost is nearly
  flat in Γ where rkdp45's grows ~1.5–2.4×/decade — so **expect ndf15 to overtake
  at high enough Γ; check rather than assume.**
- **Background: etd or rkdp45. ndf15 fails** (step collapse at a_ini).
  Implementing `BackgroundDerivsDiagonal` cut etd from 67 296 to 4 900 steps.

## Input keys

`dr_representation = proxy` plus: `dr_N_q` (96), `dr_q_min`, `dr_q_max`,
`dr_f_ini_l` (1), `dr_f_ini_phi` (0), `n_gauss`, `dr_rate_cap` (1e3),
`dr_rta_C3`, `dr_rta_C5`, `dr_rta_n3`, `dr_rta_vshut`, `dr_rta_exchange` (1),
`dr_rta_form` (`powers`), `dr_write_psd` (no).

## Known limits

- The rate coefficients are a calibration to one campaign, over a γ window 0.7
  decades wide, which that campaign itself says a single power fits as well.
- Γ ≥ 10¹⁰ has no converged exact reference to test against.
- Synchronous gauge only, scalars only, no fluid approximation — same guards as
  the exact composite.
