# Precision files

CLASS++ ships one set of compiled defaults (`include/precision.h`) and a small
number of `.pre` files that override them. Pass one as the second argument:

```
./class my_model.ini pre/precision_tuned.pre
```

Every number quoted in these files was measured, not estimated. The metric is

> **S** = the standard deviation of Δχ² against a converged reference, over 120
> draws from the Planck 2018 `base_mnu` posterior, through the Planck 2018
> `plik_lite` TTTEEE likelihood.

A constant offset in Δχ² cancels in the posterior and is not penalised; what
biases parameters is how Δχ² *varies* across the parameter space, and S measures
exactly that. S ≈ 2 × the induced parameter bias in σ, so S = 0.5 is roughly a
0.1σ shift. The reference is itself converged to S = 0.008.

| file | S | speed vs defaults | use it for |
|---|---|---|---|
| *(compiled defaults)* | 0.387 | 1.00× | everything |
| `precision_fast.pre` | 0.554 | 2.4× faster | ΛCDM only — see the warning in the file |
| `precision_tuned.pre` | 0.301 | 1.7× slower | when accuracy matters more than throughput |
| `precision_precise.pre` | 0.009 | ~60× slower | reference-grade, for convergence tests |
| `precision_pre_2026_08_21.pre` | 1.725 | 2.2× slower | reproducing results from before 2026-08-21 |

## Two things these files do not cover

**Low-ℓ EE.** `plik_lite` starts at ℓ = 30, so nothing below that entered the
objective these settings were tuned on. `fast`, `tuned` and `precise` all leave
`k_step_super` at a cheap rung, which puts the ℓ = 2–29 EE residual at 1.1–2.7×
the full-sky cosmic variance of a multipole. τ is measured there. Set
`k_step_super = 0.002` before pairing any of them with a low-ℓ EE likelihood
(Planck lowE, SRoll2); it costs about 10%. **The compiled defaults do not have
this problem** — they carry `k_step_super = 0.002` for exactly this reason.

**Nonlinear.** `non linear = none` throughout. The halofit and HMcode knobs were
never exercised and are at their own defaults. Linear P(k) *was* measured: the
defaults are accurate to ~0.2% rms over 10⁻³ < k < 10 h/Mpc. The *nonlinear*
spectrum is a different matter — `nonlinear_min_k_max = 5` truncates the integrals
halofit needs, costing 5–9% at k ≳ 0.5 /Mpc on defaults old and new alike. That is
issue #399 and is not something these files address.

**What the presets still carry.** The compiled defaults hold four knobs back from
the campaign's configuration (table above). The presets do not: `fast` and `tuned`
are the campaign's own settings, so they reintroduce the low-ℓ EE, P(k)-truncation,
lensing-headroom and tensor-BB issues. Each file's header names the affected
sector, the measured size, and the one line that fixes it. `precise` is the
exception — it no longer overrides `k_max_tau0_over_l_max`, because a
reference-grade preset that is 12% wrong on tensor BB is not a reference; that cost
S 0.0083 → 0.0087.

**Tensors.** The campaign scored scalars only, at ℓ ≥ 30, with no `mPk`. **Four**
knobs are therefore held at their pre-2026-08-21 values in the compiled defaults,
because that objective could not constrain them:

| knob | what it would have broken |
|---|---|
| `k_step_super` | low-ℓ EE, hence τ — the residual reached 1.6× cosmic variance |
| `k_min_tau0` | P(k) below k ≈ 3×10⁻⁴ h/Mpc, to 2.5% |
| `delta_l_max` | lensing headroom — the lensed spectrum degrades to 7% at the top of the requested range |
| `k_max_tau0_over_l_max` | tensors need a wider k range per ℓ than scalars; a tensor scenario went to 2.5% |

Each was caught by running the reference-comparison workflow and asking, for
every failing scenario, whether the new defaults are *closer to a converged
reference* than the old ones — not merely different from them.

A fifth, `tight_coupling_trigger_tau_c_over_tau_k`, was **split by mode** instead
of held back. It cost 14% on tensor BB over ℓ = 100–300 at the campaign's value,
because while the tight-coupling approximation is on the tensor equations set the
photon polarisation quadrupole — the entire B-mode source — identically to zero,
and unlike the scalar case there is no higher-order closure to fall back on
(every `tight_coupling_approximation` scheme gives bit-identical tensor output;
the selector is only ever consulted from the scalar branch). Scalars keep the
loose, measured-good 0.04; tensors get their own
`tight_coupling_trigger_tau_c_over_tau_k_ten`, still 0.01. Splitting rather than
reverting keeps 26% of the runtime, and follows the `_ten` convention CLASS
already uses for `l_max_g_ten` and `l_max_pol_g_ten`.

## Reproducing these numbers

The harness lives on the **`benchmarks/precision-campaign`** branch, under
`benchmarks/precision/` — it is deliberately kept off `master`, like
`benchmarks/class-public-comparison`. `s23_write_pre.py` regenerates this
directory; `paths.py` takes its locations from the environment. The campaign
report explains the method and the measurement behind every claim above:

<https://claude.ai/code/artifact/bc8ddc19-fd2d-4de7-a05c-58c974af45b9>
