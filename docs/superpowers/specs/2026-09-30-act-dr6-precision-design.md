# ACT DR6 precision: what the data need, and what CLASS++ was missing

**Status (2026-09-30).** Shipped on branch `act-dr6-precision`: the streamed
accurate-lensing transform (§4.1), `recfast_calibration` (§4.2) and
`pre/act_dr6.pre` (§4.3). **Not done:** HMcode 2020 (§5.1), a posterior-averaged
S for ACT (§5.2), the multi-frequency likelihood (§5.3), any change to the
compiled defaults (§5.4). The harness is on the not-for-merge branch
`benchmarks/act-dr6`.

Evidence tags: **[M]** measured here, **[C]** the code says so, **[L]** published,
quoted with its source, **[?]** no evidence, do not act on it.

## 1. The question

Someone heard that ACT DR6 needs (a) a more precise recombination history,
"perhaps forcing HyRec", and (b) an accurate lensing potential to ℓ ≈ 6000. How
accurate, exactly? How much of that is still true, and how much comes from old
class_public settings? CLASS++ has no users analysing ACT yet, so nobody has
complained. That is not the same as it working.

## 2. Sources

- **C25**: Calabrese et al., ACT DR6 extended models, arXiv:2503.14454, app. A
  "Theory specifications". It holds ACT's published CAMB and CLASS settings and
  the Δχ² behind each choice.
- **L25**: Louis et al., ACT DR6 power spectra and ΛCDM, arXiv:2503.14452. The
  multi-frequency spectra reach ℓ = 8500.
- **CAMB 2**: Lewis, arXiv:2607.14854 (2026-07-16), plus CAMB's `CHANGELOG.md`
  and `fortran/recfast.f90` at 56a95f7.
- The likelihood repositories `ACTCollaboration/DR6-ACT-lite`
  (`yamls/theories/camb.yaml`) and `ACT-DR6-parameters`.

## 3. What the literature says

**Multipole range [L].** L25 uses lensed TT/TE/EE to ℓ = 8500. ACT-lite, the
foreground-marginalised likelihood, has bandpowers over 600 ≤ ℓ ≤ 6500 whose
windows reach ℓ = 8502. ACT's CLASS recipe accordingly sets
`l_max_scalars = 9500` and `delta_l_max = 1800`, i.e. unlensed spectra and
C_L^φφ to ≈ 11300.

**Recombination [L].** C25: at the L25 best fit, Δχ²(RecFast − HyRec) = 0.37,
Δχ²(RecFast − CosmoRec) = 0.49 and Δχ²(HyRec − CosmoRec) = 0.12. These are
differences in χ² *to the data*. Away from the best fit they reach O(100),
mostly from TT at ℓ ~ 2000, yet "no significant difference in parameter
posteriors is found". ACT chose CosmoRec in CAMB, and HyRec in CLASS only
"since CosmoRec is not available with class".

CAMB 2 then refits RECFAST to a high-accuracy CosmoRec history, which is now
CAMB's default. The maximum x_e residual to CosmoRec falls:

| | Planck-era RECFAST | refit | HyRec-2 |
|---|---|---|---|
| hydrogen window | 9.8e-4 | 1.1e-4 | 3.5e-4 |
| helium window | 3.3e-3 | 7.8e-5 | 4.5e-3 |

On CAMB 2's SO-like Gaussian distance over 2 ≤ ℓ ≤ 4000, CosmoRec is 0.56 from
Planck-era RECFAST, 0.03 from the refit and 3.8 from HyRec-2. The HyRec-2 gap
is its simpler helium radiative transfer.

**So "forced to use HyRec" is dated.** HyRec was the CLASS stand-in for
CosmoRec. The RECFAST refit tracks CosmoRec better than HyRec-2 does.

**Lensing [L].** CAMB 2: with default lensing settings, v1.6.6 is 1.5% wrong in
TT by ℓ = 3000 (SO-like Δχ² ≈ 6000). With matched high-k lensing support it is
still 2.6% wrong in TT at ℓ = 4000. The need for accurate lensing at high ℓ is
real and not a class_public artefact.

**Nonlinear [L].** C25: Δχ²(Halofit − HMcode 2020) = 0.09. Baryonic feedback
adds 0.20, and 0.39 with stronger AGN feedback. ACT's CLASS recipe asks for
`hmcode_version: '2020'`.

## 4. What was measured in CLASS++, and what shipped

**Harness [M].** A standalone transcription of ACT-lite, with no Cobaya
dependency, reproduces the collaboration's own test value −395.48 with CAMB
1.6.6 (−395.482). Every spectrum difference is scored as the distance
dist = (t₁ − t₂)ᵀ C⁻¹ (t₁ − t₂) under ACT-lite's covariance. This is
data-independent. √dist bounds the parameter shift, in σ, along the worst
direction. It is not C25's χ²-to-data difference, which carries a cross term
2 rᵀ C⁻¹ δt with the data residual and so grows away from the best fit. Where
C25's quantity is quoted, it is labelled "ddata".

**Reference [M].** The reference has tight CLASS++ unlensed spectra and φφ,
then full-range Gauss–Legendre lensing by `camb.correlations` (independent code).

- REF_A vs REF_B (one rung tighter on every knob that moved high-ℓ spectra):
  0.0002, bands ≤ 3e-4 through ℓ = 8500.
- REF_C (REF_A with ncdm tolerances tightened) vs REF_A: 0.0012.

REF_C is the reference below, with a floor of ~0.001.

**Fiducial.** One ACT-like ΛCDM point: ω_b 0.0225, ω_cdm 0.1193, H0 67.6,
n_s 0.974, τ 0.056, one 0.06 eV neutrino.

### 4.1 accurate_lensing was unusable at ACT's ℓ

With ACT's published CLASS settings, the lensing module took 190 of a 198 s run
[M]. It tabulated every Wigner d^l_mm′(μ) at every quadrature node: 12 GB at
num_mu ≈ l_max ≈ 11300, transformed back with a row stride, on one thread while
ten idled [C].

The fast (default) mode is not a substitute. Against the full-sky reference it
is 2–5% wrong in lensed TT at 3500 ≤ ℓ < 6500 and ~10% above [M].

**Shipped.** The nodes are streamed in chunks of 32. Each task tabulates its
chunk, forms the correlation functions and adds its share of the quadrature
into per-chunk partial sums, which are added in chunk order.

- Output agrees with the previous code to ≤ 2e-10 and is bit-identical at 1, 7
  and 11 threads [M].
- ACT's settings: 198 s → 8.3 s, peak RSS 2.0 GB [M].

### 4.2 recfast_calibration

**Shipped.** `recfast_calibration` selects `planck` (the default) or `cosmorec`
(CAMB 2's seven constants plus the He I rate factor F_He, gated exactly as in
CAMB) [C]. Measured against CAMB 2.0 [M]:

- x_e agrees with CAMB 2's `recfast_cosmorec` to ~1e-5 over 600 < z < 1400.
- Residuals of up to 1.3e-3 near z ≈ 1500, 1900 and 2700 appear identically for
  planck, cosmorec and HyRec. They come from CLASS++'s phase-switching driver,
  not from the fit.
- The planck → cosmorec shift of the lensed spectra is 0.0275 on ACT-lite in
  CLASS++ and 0.0274 in CAMB 2. The distance between the two shifts is < 1e-5.

On ACT-lite at the fiducial, CLASS++ gives [M]:

| pair | dist | ddata |
|---|---|---|
| RECFAST planck vs cosmorec | 0.027 | −0.45 (C25, at its best fit: 0.49 in magnitude) |
| HyRec-2 vs RECFAST cosmorec | 0.017 | +0.80 |
| RECFAST planck vs HyRec-2 | 0.072 | −1.25 |

**For ACT DR6 alone the recombination choice is worth ≤ 0.17σ in the worst
direction.** That is consistent with C25 finding no posterior shift. It grows
with SO [L].

### 4.3 The numerical settings

All scored against REF_C on ACT-lite (dist; wall time on an M3 Pro, 11 threads,
minimum of 3 runs) [M]:

| configuration | dist | wall |
|---|---|---|
| compiled defaults, l_max_scalars 9500 | 0.394 | 1.2 s |
| ACT's published CLASS settings, translated | 0.059 | 8.3 s |
| ACT's settings + class_public's `neglect_CMB_sources_below_visibility` and `q_logstep_limber` | 0.003 | 8.7 s |
| **`pre/act_dr6.pre`** | **0.009** | **4.9 s** |

What the ladder and the knock-outs show [M]:

- **The defaults' 0.39 is not lensing.** It is EE/TE at 600 ≤ ℓ < 3500.
  Raising `radiation_streaming_trigger_tau_over_tau_k` from the campaign's 26 to
  240 removes 0.35 of it. Fixing lensing alone (accurate mode, headroom, high-k
  source) cuts high-ℓ TT errors from 8.5% to 0.46% but leaves the distance at
  0.40. ACT-lite puts little weight above ℓ ≈ 3500, where foreground
  marginalisation dominates.
- **ACT's recipe is sound, but not in CLASS++ as-is.** It never sets
  `neglect_CMB_sources_below_visibility` or `q_logstep_limber`, whose CLASS++
  defaults (0.008 and 1.1, chosen for Planck) are looser than class_public's
  (1e-3 and 1.025) [C]. With class_public's two values it reaches 0.003.
  `transfer_neglect_late_source` is redundant once the visibility cut is tight
  (0.0155 either way).
- **The preset** keeps ACT's damping-tail knobs that measurably matter and the
  visibility cut. It moves the high-k lensing source to `nonlinear_min_k_max`
  (1/Mpc), so it works without `mPk` in `output`. It adds `q_logstep_limber`
  from the reference and the CosmoRec calibration as physics. Lensed TT is within 1e-3 of the reference to
  ℓ = 6500 and within 5e-3 to ℓ = 8500. It was measured from the `.pre` file
  itself, as `./class my.ini pre/act_dr6.pre`, with RECFAST and HyRec, and runs
  with HMcode (7.9 s; not scored, since the reference uses halofit). Every key is
  read.

**Rejected.** Shipping ACT's own settings (0.059, 8.3 s): 7× worse and 1.7×
slower than the preset. Shipping ACT's settings plus the two class_public
values (0.003, 8.7 s): better, but 1.8× dearer for a gain of √0.009 − √0.003 ≈
0.04σ in the worst direction, below what one fiducial point can support
(§5.2).

## 5. Open

1. **HMcode 2020 [C].** CLASS++ implements HMcode 2016 only, and
   `hmcode_version` is not a key. Through `classy`, ACT's yaml fails loudly with
   "Class did not read input parameter(s)" rather than silently falling back
   [C]. Physics gap: 0.09 in C25's ddata, plus baryonic feedback [L].
2. **One point, not a posterior [?].** No ACT posterior draws were available,
   so there is no S like the Planck presets have. Extended models are
   unmeasured.
3. **The multi-frequency likelihood (MFLike, to ℓ = 8500) [?].** Not run.
   High-ℓ lensed CMB sits under the foreground power there, but that is not
   measured here.
4. **Compiled defaults [?].** `radiation_streaming_trigger_tau_over_tau_k = 26`
   and `neglect_CMB_sources_below_visibility = 0.008` were tuned for Planck
   `plik_lite` by the 2026-08-21 campaign and cost ACT 0.35 and ~0.05.
   Changing them moves the Planck S and cost, which were not re-measured, so
   they are untouched here.
5. **Silent mismatches that remain [C].** ACT's class_public key
   `perturbations_sampling_stepsize` is `perturb_sampling_stepsize` in CLASS++.
   `classy` refuses the unknown name, but the `class` binary only warns under
   `write warnings = yes`.
   BBN: CLASS++'s `YHe = BBN` reads `bbn/sBBN_2017.dat` (PArthENoPE 2017). CAMB
   2 defaults to PRIMAT 2024; C25 used PRIMAT 2021 [C]. The two PRIMAT tables
   differ by 2.0e-4 in Y_He [L]. The PArthENoPE–PRIMAT difference and its
   ACT-lite effect are not measured here [?].
6. **The RECFAST driver's He-window residual [M].** Up to 1.3e-3 in x_e near
   z ≈ 1500 against CAMB, shared by every model. The cross-code distance at
   matched calibration is 0.016, which bounds it together with every other
   code difference.

## 6. How to reproduce

The harness lives on branch `benchmarks/act-dr6`, directory `benchmarks/act_dr6/`
(scripts only; `fetch.sh` downloads the ACT-lite data, both CAMB versions and
the papers). `run_refs2.py` builds the references (REF_B: 18 min), `presets.py`
and `final_preset.py` the tables above, `recomb_score.py` and `xe_compare2.py`
§4.2, and `ab_lensing.py` the lensing A/B.
