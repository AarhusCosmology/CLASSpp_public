# act_dr6.pre for ACT DR6, SPT-3G D1 and their lensing, and for self-interacting neutrinos

**Status (2026-10-29).** Shipped: a new `pre/act_dr6.pre` (§5), measured on four likelihoods
at two points. **Not done:**
ACT's multi-frequency likelihood MFLike, a posterior-averaged S, other models (§6). The harness
is local, in `benchmarks/lmax_timing/` (§7).

Evidence tags: **[M]** measured here, **[C]** the code says so, **[L]** published,
**[?]** no evidence, do not act on it.

## 1. The question

The SIν + DRMD runs on the H0 World Cup data (ACT-lite, SPT-3G D1 lite, ACT DR6 + Planck PR4
lensing, SPT-3G MUSE lensing, with Planck at low ℓ, DESI and Pantheon+) used `act_dr6.pre`
without its `nonlinear_min_k_max = 70`, `nonlinear_min_k_max = 10` instead, and on top of it
`evolver_perturbations = 1` and `tol_perturb_integration = 1e-7` from earlier runs, at
`l_max_scalars = 9001`. One evaluation took 38 s on one thread (70 core seconds on the
cluster). Are these settings needed for those likelihoods; is CLASS++ still that slow after
the lensing and transfer work of 2026-10-26/27; and does the preset need changing?

Here SIν is three degenerate `ncdm_self_interacting` neutrinos (N_ur 0.00441, so that
N_eff = 3.044 without the dark radiation), and DRMD the `IDM_DRMD_IDR_DRMD` composite.

## 2. What the likelihoods use

| likelihood | asks CLASS for | uses |
|---|---|---|
| ACT-lite (`act_dr6_cmbonly` 1.0) | lensed TT/TE/EE to 9001 | ℓ ≤ 6325 **[M]** |
| SPT-3G D1 lite (candl 2.2.0, `spt_candl_data` 3.0.1) | lensed TT/TE/EE to 4095 | 350 ≤ ℓ ≤ 4095 **[M]** |
| ACT DR6 + PR4 lensing (`act_dr6_lenslike` 1.2.1, `actplanck_baseline`) | lensed TT/TE/EE/BB and φφ to 4000 | L ≤ 2998 **[C]** |
| SPT-3G MUSE, φφ (`muse3glike` 1.0.0) | φφ to 3000, unlensed EE to 5000 | 20 ≤ L ≤ 3073 **[M]** |

- **ACT-lite** keeps the bandpowers at 600 ≤ ℓ ≤ 6500 by default, and every window of those
  bandpowers is exactly zero above ℓ = 6325 **[M]**. Its code asks for
  max(`lmax_theory`, top of the window arrays + 1) + 1, and the window arrays run to 8501, so it
  never asks for less than 8503, whatever `lmax_theory` says **[C]**. Cobaya's classy theory sets
  `l_max_scalars` to the largest request of any likelihood (`set_cl_reqs`) **[C]**, so the
  request alone forces ℓ = 9001.
- The subclass below cuts the windows to their non-zero rows. Against the unmodified code, with
  the cobaya base class stubbed out, the χ² is identical to the last bit: on the reference
  spectra, with every C_ℓ above 6326 set to 1e30, and with spectra that stop at 6326 **[M]**.
  The request becomes 6327.
- MUSE asks for unlensed EE to 5000 also when only φφ is used. Its wrapper says this works
  around a Cobaya bug **[C]**.

```python
from types import SimpleNamespace
import numpy as np
from act_dr6_cmbonly import ACTDR6CMBonly

class ACTDR6CMBonlyTrimmed(ACTDR6CMBonly):
    lmax_theory: int = 0    # 0: the floor comes from the windows

    def initialize(self):
        super().initialize()
        top = 0
        for m in self.spec_meta:
            w = m["window"]
            last = int(np.flatnonzero(np.any(w.weight != 0, axis=1)).max())
            m["window"] = SimpleNamespace(values=w.values[:last + 1], weight=w.weight[:last + 1])
            top = max(top, int(w.values[last]))
        self.lmax_theory = top + 1
```

## 3. Harness

- **Metric [M].** For each likelihood, dist = (t − t_ref)ᵀ C⁻¹ (t − t_ref) in its own bandpower
  space, √dist being the bound on the parameter shift in σ along the worst direction (as in the
  2026-09-30 note). ACT-lite and SPT through their windows (SPT with Tcal = Ecal = 1); ACT lensing
  through its binning and its CMB-dependent norm and N1 corrections; MUSE in its transformed
  space, where it is Gaussian. The likelihoods are independent, so "all" is the sum.
- **Likelihood checks [M].** SPT-3G D1 lite reproduces its test value (relative difference 0.0),
  and act_dr6_lenslike passes its own 12 tests. MUSE has no Python test value; with the CLASS++
  ΛCDM φφ its χ² is 21.3 for 16 bandpowers. It and spt_candl_data were installed without their
  pins (numpy < 2; jax, gpjax), with numpy 2.5.
- **References [M].** REF_C's recipe of the 2026-09-30 note, run with master 7536377a, lensed by
  `camb.correlations` (full-range Gauss–Legendre). The ΛCDM reference is 0.05σ (√dist on
  ACT-lite) from the September one with the same recombination, about that note's floor. For
  SIν + DRMD, REF_C's neutrino hierarchy (`l_max_ncdm` 50, fluid trigger 200,
  `tol_ncdm_synchronous` 1e-6) under ndf15 did not finish in 30 min; the reference uses
  (35, 120, 3e-4). In the ladder of §4.2 that rung is 0.04σ from (25, 60, 1e-3) and 0.08σ from
  (50, 200, 1e-4), so the SIν numbers are not resolved below ~0.05–0.08σ on ACT-lite.
- **Points.** ΛCDM: the ACT-like point of 2026-09-30. SIν + DRMD: the best fit of the World Cup
  runs without SH0ES (H0 69.48, log10 G_eff −2.859, ΔN_eff 0.223, f_idm 0.0116,
  log10 z_stop 4.50, log10 (G/H)_ini 7.26). Halofit, RECFAST, `l_max_scalars` 6327 or 9001.
- **Timings [M].** Apple M3 Pro, mains power, one thread, classy: modules timed one by one
  through `compute(level=...)`, the faster of two warm runs (the Bessel table kept, as in Cobaya).

## 4. Results

### 4.1 The settings, on the four likelihoods (√dist)

| point | settings | ACT-lite | SPT | ACT lens | MUSE | all | 1 thread, ℓ 6327 |
|---|---|---|---|---|---|---|---|
| ΛCDM | new file | 0.088 | 0.083 | 0.022 | 0.012 | 0.123 | 4.8 s |
| | defaults | 0.156 | 0.109 | 0.020 | 0.013 | 0.192 | 1.5 s |
| | file until 10-29 | 0.103 | 0.081 | 0.042 | 0.040 | 0.143 | 7.6 s |
| SIν + DRMD | new file | 0.105 | 0.091 | 0.042 | 0.058 | 0.157 | 8.9 s |
| | defaults | 0.181 | 0.117 | 0.023 | 0.017 | 0.217 | 4.7 s |
| | file until 10-29 | 0.434 | 0.252 | 0.203 | 0.053 | 0.544 | 13.8 s |
| | the runs' settings, ℓ 9001 | 0.413 | 0.244 | 0.197 | 0.053 | 0.521 | 33.3 s |

At ℓ 9001 the new file gives 0.119 and 0.155 in 7.4 and 13.3 s **[M]**.

### 4.2 Why the old file fails for SIν: its neutrino lines [M]

The old file kept the neutrino defaults of 2026-08-21 (`l_max_ncdm` 10, fluid trigger 15,
`tol_ncdm_synchronous` 0.03). The latter also sets the momentum quadrature, which is adaptive.
In ΛCDM they act on one 0.06 eV state next to 2.03 massless neutrinos. In SIν every neutrino is
in that hierarchy. ACT-lite, every other line of the old file kept, against the top rung
(50, 200, 1e-4), at ℓ 6327:

| neutrino lines | √dist |
|---|---|
| the old file's (10, 15, 0.03) | 0.62 |
| defaults, but `l_max_ncdm` 10 | 0.73 |
| defaults, but trigger 15 | 0.44 |
| defaults, but tol 0.03 | 0.28 |
| defaults (17, 31, 0.0034) | 0.08 |
| (25, 60, 1e-3) | 0.08 |
| (35, 120, 3e-4) | 0.08 |

The three lines are also what puts SIν 0.20σ off on ACT lensing: without them, but with the
file's other pins, it is 0.020. The other pins (constant transfer cuts, Bessel sampling 8, the
k range) take SIν on ACT-lite from 0.116 to 0.201 against the reference (default neutrino lines
in both); for ΛCDM all the pins together were 0.103 against 0.094. Dropping them is faster.

### 4.3 Which lines matter [M]

Each line removed from the new file plus the old ur and visibility lines (√dist, ΛCDM / SIν):

| without | ACT-lite | SPT |
|---|---|---|
| (nothing) | 0.094 / 0.116 | 0.085 / 0.093 |
| `radiation_streaming_trigger_tau_over_tau_k` 240 | 0.120 / 0.149 | 0.088 / 0.097 |
| `perturb_sampling_stepsize` 0.05 | 0.132 / 0.149 | 0.106 / 0.112 |
| `l_linstep` 20 | 0.108 / 0.132 | 0.089 / 0.099 |
| `delta_l_max` 1800 | 0.118 / 0.142 | 0.084 / 0.093 |
| `accurate_lensing` | 0.094 / 0.117 | 0.084 / 0.093 |
| `neglect_CMB_sources_below_visibility` 1e-4 | 0.094 / 0.115 | 0.088 / 0.094 |
| `l_max_ur` 35, `ur_fluid_trigger_tau_over_tau_k` 130 | 0.087 / 0.107 | 0.080 / 0.090 |

- Neither lensing likelihood moves by more than 0.02 for any of these.
- `q_logstep_limber` 1.025 changes neither and costs nothing measurable (4.78 s either way).
- `accurate_lensing` changes none of the four by more than 0.001, and costs 1.8 s (ΛCDM) and
  1.9 s (SIν). The fast mode is 1.9% wrong in lensed TT at 3500 < ℓ < 5000, which the four do
  not weight.
- `nonlinear_min_k_max`: 2 (the default), 10 and 70 agree on the four within 0.04 for SIν and
  0.002 for ΛCDM. 70 costs 51% (ΛCDM) and 57% (SIν) more than 10. 10 stays, for the reason in the
  runs' draft: the continuation beyond the last mode has the high-k form of cold dark matter,
  and at points with dark acoustic oscillations the default moved C_L^φφ by up to 9% at L ≤ 400
  (the runs' authors' measurement, not repeated here).
- `tol_perturb_integration` 1e-7 against the default 1e-5 (SIν): 0.018–0.026 on ACT-lite, 0.014
  on SPT, ≤ 0.002 on the lensing likelihoods. None is closer to the reference (ACT-lite 0.182
  against 0.181), and it doubles the run time (33.3 against 17.8 s at ℓ 9001).
- `evolver_perturbations = 1` changes nothing. ndf15 is what master picks for SIν + DRMD
  (§6.3) **[C]**.
- At the best fit itself, the χ² to the ACT-lite data differs by 0.05 between the runs'
  settings and the new file. Both are 0.40–0.45 below the reference. A distance bounds the
  parameter shift; it does not have to show up as χ² at one point.

### 4.4 Timings, bfe6729e (the runs' build) against master 7536377a [M]

The runs' settings, one thread. bfe6729e reproduces the runs' own report to 1% on this machine.

| ℓ_max | ΛCDM bfe6729e | ΛCDM master | SIν bfe6729e | SIν master |
|---|---|---|---|---|
| 2500 | 7.1 | 9.0 | 12.4 | 18.3 |
| 4000 | 9.7 | 10.6 | 17.2 | 20.9 |
| 6500 | 15.2 | 14.1 | 25.9 | 27.1 |
| 9001 | 21.9 | 17.9 | 38.0 | 33.3 |

- At ℓ 9001, transfer went from 6.9 to 4.5 s and lensing from 5.9 to 3.3 s (SIν).
- The perturbations stay at 25 s, almost all of it the tolerance.
- Below ℓ 9500 the old file is slower on master than before, because its
  `k_max_tau0_above_l_max = 10840` (2026-10-27) only reproduces its old k range at 9500.
  At 2500 it solves to k τ0 = 13340, where it used to solve to 7740.

## 5. What shipped

`pre/act_dr6.pre`:
- **Kept:** `radiation_streaming_trigger_tau_over_tau_k` 240, `perturb_sampling_stepsize` 0.05,
  `l_linstep` 20, `delta_l_max` 1800, `accurate_lensing` 1, `q_logstep_limber` 1.025, and
  `nonlinear_min_k_max` 10 instead of 70.
- **Gone:** the three neutrino lines; the k range; the eight transfer cuts; `hyper_sampling_flat`;
  the ur lines; the visibility cut.
- **Header:** the table of §4.1, the ACT-lite cut, and the advice on the evolver and the
  tolerance.

The old file is not kept as a dated copy; results made with it are reproduced from the history
(2026-09-30 to 2026-10-29).

**Rejected.**
- *One file per likelihood:* the lensing likelihoods are satisfied by the defaults, SPT gains
  only from the sampling step, and the four are always fitted together through one theory
  call.
- *Dropping `accurate_lensing`:* 22–37% faster and the same on the four, but 2% wrong in lensed
  TT above ℓ ≈ 3500. The header says when to drop it.

## 6. Open

1. **MFLike [?].** ACT's multi-frequency likelihood uses ℓ to 8500 without marginalising the
   foregrounds into the bandpowers. It was not run.
2. **Two points, not a posterior [?].**
3. **Self-interacting neutrinos and etd [M][C].**
   - SIν alone gets etd. With DRMD the content gets ndf15, because DRMD's children report
     neither explicit safety nor a diagonal. The selector applies the rule as intended
     (`BaseSpecies::SupportsExponentialPerturbationEvolver`), and DRMD is right: its drag starts
     at G/aH ≈ 1.8e7 at this best fit.
   - A diagonal would not cure it. The drag has the Thomson structure
     G·[[−1, 1], [R, −R]], for which no diagonal splitting gets below √R·G (memory note
     2026-08-19, photons and baryons).
   - Forced etd is 3.8× slower than ndf15 here. With DRMD's dark tight coupling switched on
     (trigger 100), etd is 4–10% faster than ndf15.
   - At this coupling (log10 G_eff −2.86), etd and ndf15 cost the same on SIν alone (0.95
     against 0.93 s at ℓ 2500, at default tolerances).
4. **The runs [?].** The SIν + DRMD runs used settings that are 0.52σ from converged on these
   four likelihoods at the best fit, as a bound. Whether to redo them is for their authors.

## 7. How to reproduce

Local harness `benchmarks/lmax_timing/`, scripts only:
- `lmax_timing.py`: module timings; `run_timings.sh` for the A/B.
- `actlite_master.py ref` and `POINT=SIDR ... ref`: the references; its `score` builds the
  variants and scores ACT-lite.
- `score_all.py` and `final_table.py`: the four likelihoods.
- `ncdm_ladder.py`: §4.2.
- `act_dr6_cmbonly_trimmed.py` and `test_trimmed.py`: §2.
- `plots.py`: figures.

The likelihoods are installed into the `benchmarks/act_dr6/venv` of the 2026-09-30 note.
