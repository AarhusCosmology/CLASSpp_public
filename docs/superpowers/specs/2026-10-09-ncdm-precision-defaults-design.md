# Neutrino precision defaults: back to before #397

**Status (2026-10-09).** Implemented on branch `ncdm-precision-defaults`. Three
non-cold-relic knobs go back to their values from before 2026-08-21:
`ncdm_fluid_trigger_tau_over_tau_k` 15 → 31, `l_max_ncdm` 10 → 17 and
`tol_ncdm_synchronous` 0.03 → 0.0034. With massive neutrinos, the linear P(k)
error against a converged run drops from 3.8e-3 to 4.3e-4 for k ≤ 5/Mpc
(m_ν = 0.06 eV), for 2% of runtime on one thread. The presets keep the values
they were measured with. The Planck S of the defaults has not been re-measured;
a fuller comparison run is planned (the user's call).

Evidence tags: **[M]** measured here, **[C]** the code says so, **[?]** no
evidence.

## Problem

The precision campaign of #397 tuned the defaults on one objective: Planck
`plik_lite` TTTEEE, with `non linear = none` and no `mPk`. It loosened the three
knobs above, which that objective could not see at high k. The campaign's P(k)
check reported an rms over 10⁻³–10 h/Mpc, which hides the high-k end.

Found while validating HMcode 2020 against CAMB (PR #452). With m_ν = 0.06 eV
the HMcode boost was 3.3e-3 off CAMB at k ≤ 5/Mpc, against 1.1e-3 in ΛCDM, in
both HMcode versions. class devel's own HMcode gets 1.1e-3 there [M]. The halo
model parameters of the two CLASS codes agree to 1e-3; the linear spectra do
not. Relative to CAMB, normalised at k = 0.01/Mpc, CLASS++'s P_m at z = 0 falls
to −3.8e-3 at k = 5/Mpc, class devel's to −1.2e-3 [M].

## Cause [M]

Against a converged run (fluid trigger 200, `l_max_ncdm` 120,
`tol_ncdm_synchronous` 1e-5, `tol_perturb_integration` 1e-7; trigger 50 with
l_max 30 agrees with it to 6e-5 for k ≥ 0.1/Mpc), linear P_m(k) at z = 0,
m_ν = 0.06 eV:

| settings (trigger, l_max_ncdm, tol) | k = 0.5 | k = 1 | k = 2 | k = 5/Mpc |
|---|---|---|---|---|
| 15, 10, 0.03 (2026-08-21 defaults) | −3.7e-4 | −1.3e-3 | −2.4e-3 | −3.8e-3 |
| 15, 10, 1e-5 | −2.0e-3 | −3.1e-3 | −4.3e-3 | −5.9e-3 |
| 31, 17, 0.03 | +2.0e-3 | +2.2e-3 | +2.3e-3 | +2.5e-3 |
| 31, 17, 1e-3 | −7.3e-4 | −8.3e-4 | −9.5e-4 | −1.1e-3 |
| **31, 17, 0.0034** (before 2026-08-21, and now) | +1.6e-4 | +1.6e-4 | +1.3e-4 | +7.7e-5 |
| 25, 17, 0.03 | +1.5e-3 | +1.3e-3 | +1.1e-3 | +6.5e-4 |

* The error does not depend on the neutrino mass: m_ν = 0.06 and 0.3 eV give
  the same numbers at k ≥ 0.5. At kτ = 15 the neutrinos are still relativistic
  at those k, so the fluid approximation starts too early whatever their mass.
* The knobs interact. The loose q tolerance partly cancels the early switch
  (tightening it alone makes the defaults worse), and a later switch needs the
  longer hierarchy (trigger 31 with l_max 10 is 1.3e-2 off).
* The error does not fall smoothly with the tolerance: 1e-3 is further from
  converged than 0.0034. The momentum grid changes in steps [?]. 0.0034 is the
  pre-#397 value, so it is kept rather than re-tuned on one cosmology.
* A run without the fluid approximation is not a reference: truncated at
  `l_max_ncdm = 50`, it is 1e-2 above CAMB at k = 5.

## Change

* `include/precision.h`: the three defaults, with a comment citing this spec.
* `pre/precision_pre_2026_10_09.pre`: the three values of 2026-08-21, to
  reproduce results from 2026-08-21 to 2026-10-09.
* `pre/precision_fast.pre` and `pre/precision_tuned.pre` relied on the old
  defaults for some of these knobs. They now set them explicitly, so they stay
  the configurations whose S was measured; `fast` gets a header note on P(k).
  `precise` already set all three.
* `pre/README.md`: seven knobs held back from the campaign instead of four;
  the S of the defaults is marked as measured before this change.
* `python/test_class.py`: the l_max_dr_col test sets `l_max_ncdm = 10`, so its
  minimum rule still binds (the derived default would otherwise be 17).

`l_max_ncdm` also sizes the hierarchies of the other species of the NCDM family
(decaying NCDM, dcdm_wdm daughters, DR phase-space species) when they do not set
their own [C].

## Results [M]

Measured on the branch build; "before" passes the three 2026-08-21 values to
the same build.

**Linear P(k), z = 0**, maximum |P/P_converged − 1| for k ≤ 5/Mpc:

| | m_ν = 0.06 eV | m_ν = 0.3 eV |
|---|---|---|
| before | 3.8e-3 | 3.8e-3 |
| after | 4.3e-4 | 1.8e-3 |

* For m_ν = 0.3 eV what remains sits at k ≈ 0.017/Mpc, near the free-streaming
  scale. It is the fluid approximation itself: trigger 50 with l_max 30 still
  leaves 1.6e-3, and the tolerance and `tol_perturb_integration` do not move it.
* The converged run agrees with CAMB 2.0.4 to 2.5e-4 for m_ν = 0.06 eV,
  normalised at k = 0.01/Mpc. HMcode 2016 and 2020 against CAMB, which were
  3.3–4.4e-3 off with m_ν = 0.06 eV, come to 1.1e-3 and 7.4e-4, as in ΛCDM
  (PR #452 spec, §Massive neutrinos).

**Cost**, median wall time, before and after interleaved on one build:

| case | 1 thread | 11 threads |
|---|---|---|
| lensed CMB + P(k), m_ν = 0.06 eV | +2.3% | −0.1% |
| P(k) to k = 20/Mpc and z = 2, HMcode | +1.5% | |
| lensed CMB + P(k), three 0.02 eV neutrinos | +1.5% | |
| decaying NCDM (`ncdm_decay_dr`, Γ = 1e3, TT to ℓ = 200) | +29% | |

Decaying NCDM pays for the longer hierarchy in its kernel. Its recommended knob
set already used `l_max_ncdm = 17`
(`docs/superpowers/specs/2026-08-17-dncdm-minimal-knob-set.md`).

**CMB at the ACT-like point of #448** (ACT-lite distance to `refA_rc_hf`, with
`l_max_scalars = 9500`):

| | defaults | `act_dr6.pre` |
|---|---|---|
| master | 0.364 | 0.0096 |
| branch, before pinning the preset | 0.356 | 0.044 |
| branch, preset pinned (shipped) | 0.356 | 0.0096 |

With otherwise default settings, the neutrino sector alone is 0.051 at the old
values and 0.002 at the restored ones; TT and EE move by 8e-4 below ℓ = 2000.

**The ACT preset relied on the old values.** Its 0.0096 includes a cancellation:
the old neutrino settings shift TT at 600 < ℓ < 2000 by about the same amount
as a residual of the transfer integrals, in the opposite direction. With the
reference's own neutrino settings the preset measures 0.057. Adding one reference
knob at a time to the preset (with the restored defaults) finds the residual:
`transfer_neglect_delta_k_S_t0 = 1` takes it to 0.012, adding
`transfer_neglect_delta_k_S_t2 = 1` to 0.0046 (7.3 s against 4.6 s). The old
values with `t0 = 1` give 0.015. The preset keeps the old values for now, so
that its quoted number stays true; re-tuning it belongs to the planned
comparison run.

## Verification [M]

* After the rebase on master (HMcode 2020, #452): `ctest` 49/49,
  `python/test_class.py -m "not test_scenario"` 132 passed.
* Level 2 against a fresh master `classyref`: 9 of 1512 scenarios fail, all on
  P(k) in Newtonian gauge, all the same cosmology (three degenerate 0.06 eV
  neutrinos with `N_ur = 0`, scenarios 0288, 0396, 0402–0405, 0411–0413). There
  the restored defaults are the more accurate side. Against a converged run,
  P(k) up to k = 1.9/Mpc is off by:

  | | synchronous | Newtonian |
  |---|---|---|
  | master | 7.0e-3 | 1.4e-2 |
  | restored | 1.3e-3 | 3.4e-3 |

  With all neutrino radiation in the massive species, the early fluid switch
  costs more than with one massive neutrino, and Newtonian gauge more than
  synchronous.
