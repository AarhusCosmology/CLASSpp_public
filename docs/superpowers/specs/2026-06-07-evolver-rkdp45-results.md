# rkdp45 evolver — correctness & benchmark results

**Date:** 2026-06-07
**Branch:** `feat/evolver-rkdp45`
**Hardware:** Apple M3 Pro, 11 threads, `-O3` (no `-march=native`), `evolver = 2` selects rkdp45 (default remains `ndf15`).

Companion to the spec (`2026-06-07-evolver-rkdp45-port-design.md`) and plan
(`../plans/2026-06-07-evolver-rkdp45-port.md`).

## Summary

The Dormand-Prince 4(5) solver was ported from `PhD2024-EBH` and is **correct** —
it reproduces the `ndf15` reference observables within the established
`RTOL=1e-3` bar across the entire scenario suite — but it is **substantially
slower** than the stiff `ndf15` default on the perturbation integration, which
dominates runtime. Recommendation: merge as an opt-in solver (default unchanged)
and treat it as the validated baseline for the follow-on programme (Tsit5 /
DOP853 / PI-PID controllers); it is not a candidate to replace `ndf15` as the
default.

## Correctness (RTOL = 1e-3, `test/scenarios/compare_tol.py`)

### Smoke test (`explanatory.ini`)
rkdp45 output is **not** byte-identical to ndf15 (confirms it actually runs, not a
silent fallback) yet passes compare_tol:
- `out_cl.dat`: worst |Δ|/colpeak = 1.77e-05
- `out_cl_lensed.dat`: worst = 1.67e-05

### Full scenario suite (rkdp45 vs ndf15 reference)
Every scenario that produces spectra passes within RTOL=1e-3; worst deviation
across all of them is ~5.6e-05 (P(k), gauge_lcdm) — comfortably under the bar.

| Scenario | Result | worst |Δ|/colpeak |
|---|---|---|
| gauge_lcdm | PASS | 5.64e-05 (pk) |
| gauge_dcdm | PASS | 3.68e-05 (pk) |
| gauge_fluid | PASS | 3.01e-05 (pk) |
| gauge_idmdr | PASS | 8.37e-06 |
| gauge_ncdm | PASS | 1.13e-05 |
| gauge_scf | PASS | 9.61e-06 |
| idm_dr_full | PASS | 1.52e-05 |
| idm_drmd_full | PASS | 1.44e-05 |
| ncdm_single | PASS | 2.41e-05 |
| ncdm_self_interacting | PASS | 2.41e-05 |
| ncdm_multi_unsorted | PASS | 2.97e-05 |
| ncdm_greybody | PASS | 1.04e-05 |
| dncdm_dr | PASS | 1.28e-06 (background) |
| ncdm_dncdm_idmdr_combined | PASS | 2.19e-05 (background) |

Not solver-related (excluded): `dncdm_dr_bare_omega` is an intentionally-invalid
scenario that errors under **ndf15 too** (decaying species rejected by `Omega`
normalization); `dncdm_dr_combined` and `dncdm_dr_initial` are background-only
with no spectra output.

No stiffness-induced failures were observed — the explicit solver completed every
physics scenario, including stiff tight-coupling regimes.

## Benchmark (median of 5 loops, `class_profiled`, explicit-root inis)

The selected evolver drives **both** the Background and Perturbations
integrations (both go through the same `generic_evolver` dispatch). In practice
the measured Background cost is unchanged between evolvers (~0.7 ms either way)
because that system is small and smooth; the effect is concentrated in
**Perturbations**, which dominates runtime. Times in ms.

| Case | stage | ndf15 | rkdp45 | rkdp45/ndf15 |
|---|---|---|---|---|
| explanatory | Perturbations | 245.7 | 963.5 | **3.9×** |
| explanatory | TOTAL | 432.6 | 1141.5 | 2.6× |
| gauge_lcdm | Perturbations | 170.3 | 441.1 | **2.6×** |
| gauge_lcdm | TOTAL | 260.6 | 529.9 | 2.0× |
| ncdm_single | Perturbations | 923.5 | 1436.8 | **1.56×** |
| ncdm_single | TOTAL | 1146.4 | 1658.9 | 1.45× |

The slowdown is concentrated entirely in Perturbations and is scenario-dependent:
worst for the light/non-stiff LCDM cases (~3.9×), smallest for the NCDM case
(~1.56×) where ndf15's stiff variable-order BDF advantage narrows. This is the
expected behavior of a fixed-order **explicit** RK method on stiff perturbation
ODEs: it is forced to take many small steps through tight coupling.

## Intentional deviation from the source branch

The branch's step-rejection controller set `nofailed = _FALSE_` only in an
unreachable `else`, so the "halve step on consecutive failures" fallback was dead
code. The port corrects this to standard ode45 control (mark `nofailed` on the
first failure). This is a behavior change versus `PhD2024-EBH`, flagged here and
in the PR. The port also added a missing `idx < x_size` bound on the initial
output-index search (the original could read past the array) and replaced raw
`malloc`/`printf` with `std::vector` and `class_test`.

## Reproduction note (footgun)

`class_profiled`/`class` derive the default output `root` from the ini path; an
ini placed in `/tmp` with no `root=` yields an invalid path
(`output//tmp/...parameters.ini`) and fails in `WriteParameterFiles` during input
— independent of the evolver. Always set an explicit `root =` when benchmarking
from a scratch ini. (Minor usability gotcha in the input module, not specific to
this work.)

## Conclusion

rkdp45 is a correct, validated, opt-in addition. It is not faster than `ndf15`
and should not become the default. Its value is as a verified explicit-RK
baseline for the next increment (Tsit5, DOP853, and PI/PID step controllers),
where better tableaus/controllers may close — though for genuinely stiff regimes
are unlikely to fully overcome — the gap against the stiff integrator.
