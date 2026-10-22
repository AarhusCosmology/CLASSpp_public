# Follow-ups to the lensing-source continuation: threads, warm dark matter, HyRec

**Status (2026-10-22).** Implemented on `perf/nonlinear-followups`, on master after #470. Three
items that note left open or found:

1. the non-linear stage runs in parallel over time (24 -> 4.3-4.6 ms at 11 threads);
2. a species whose free-streaming cutoff lies in the continued range has its modes solved
   through it (thermal warm dark matter of 0.1 to 7 keV);
3. HyRec reports a missing or unreadable table under a long path as an error instead of
   dying with SIGTRAP.

On the way, ThreadSanitizer found one data race in a full run, a dead write in the NCDM
species, now deleted. Follows `2026-10-21-lensing-source-extrapolation-design.md`.

Evidence tags: **[M]** measured for this note, **[C]** read from the code, **[P]** an
earlier note, not re-measured, **[?]** no evidence.

## 1. The non-linear stage over threads

**Problem.** halofit and HMcode computed the correction R_NL one time after another. With the
continuation, halofit runs to z = 7.4 instead of 3.3, and the stage took 24 ms of 255 at 11
threads, unchanged from one thread **[M]**. 90% of it is halofit's bisection for R_NL: a spline
and an integral over the extended k grid at each step **[M]** (sampled profile).

**Design.** The times are independent **[C]**: each needs P_L(k) at that time, which it
computes itself, and writes only its own slice of `nl_corr_density_`, `nl_corr_density_extra_`
and `k_nl_`. One worker per thread takes the latest time not yet taken from an atomic counter.

- *The stop.* The serial loop ran from today backwards and stopped at the latest time at which
  the correction is not computable (k_NL beyond the trusted k_max), leaving R_NL = 1 and
  k_NL = 0 before it. A worker now stops at the first time before the latest such time found
  so far. Times taken before that was known are reset afterwards, so the result is the serial
  one. Every time after the latest stop was computed, since a time is only skipped below a
  stop **[C]**.
- *Errors.* An exception at a time the serial loop never reached is discarded for the same
  reason. HMcode 2020 can fail to converge for r_nl below its stop: a thread reached such a
  time in 2 of 12 runs of `test-nonlinear-threads` (z_max_pk = 10), which aborted the run
  before the rule **[M]**. Each time keeps its own exception (one `std::exception_ptr` per
  time, written only by the thread computing it, so no lock), and after the threads join, the
  first one from today down to the stop is thrown: the error the serial loop would have met.
  An error forced at every time (halofit_min_k_nonlinear = 30) gives master's message at 1 and
  8 threads **[M]**.
- *HMcode's tables.* Its sigma(R) table and wiggle spectrum change with time, so they move out
  of the shared `nonlinear_workspace` into `nonlinear_time_workspace`, one per worker. The
  shared part (growth tables, Eisenstein-Hu spectrum, the 512 x 512 smoothing weights) is
  filled before the loop and only read in it. Its `sigma_8`, `sigma_disp`, `sigma_disp_100`
  and `sigma_prime` tables were written at every time but read only by four interpolators that
  nothing calls; both are deleted **[C]**.

**Verification.**
- Byte-identical output to the serial loop, fast-math off, at 1, 4 and 11 threads: halofit with
  0.06 eV, 3 x 0.2 eV and w0wa; HMcode 2016, 2020, 2020 with feedback; halofit without the
  continuation, which stops at z = 2.2 and prints the same warning **[M]**
  (`benchmarks/nl_followups/ab.py`).
- ThreadSanitizer: nothing in the module, for halofit, HMcode 2016 and 2020 with feedback, and
  the stop at z = 2.2 **[M]**.
- `test-nonlinear-threads` pins 1 against 8 threads, bit for bit, for four cases that reach the
  stop, and checks that the stop happens; it fails if the reset of earlier times is removed
  **[M]**.

**Timing** (`class_profiled`, fiducial lensing run with halofit, l_max 3500):

| | 1 thread | 11 threads |
|---|---|---|
| non-linear stage | 23.8 -> 26.1 ms (code layout, see below) | 24.0 -> 4.3-4.6 ms |
| fiducial total (l_max 3500, halofit) | 1144-1202 -> 1150-1183 ms (noise) | 251-256 -> 235-240 ms (-6%) |
| CAMB-matched total (l_max 2500, halofit, mPk) | 857-899 -> 866-901 ms (noise) | 196-204 -> 180-185 ms (-7%) |

Ranges over three alternating runs, each a median of 9 loops (3 at one thread), #470 against
this branch **[M]** (`benchmarks/nl_followups/timing_final.log`). On one thread the loop does
the same work, 64 times x 2 spectra in the fiducial run, and the extra 2 ms is code layout:
with timers around the calls in both builds, the branch's halofit takes 21.0 ms against
master's 22.9 **[M]**. The source of halofit is unchanged.

## 2. Warm dark matter: the floor follows the free-streaming cutoff

**Problem.** Since #470 the modes are solved to `nonlinear_min_k_max` = 2/Mpc and the linear
P(k) that halofit and HMcode integrate is continued from there to `k_max_extrapolated` = 70/Mpc,
in the high-k form of cold matter. A free-streaming cutoff is not of that form. With 1 keV
thermal WDM as a quarter of the dark matter, P_NL at z = 3 was 5.3e-2 off and C_L^phiphi
4.4e-2, against modes solved to 70/Mpc **[M]**. The #470 note left it as a manual setting.

**Where it matters: only through the non-linear method.** Without halofit, the continued
linear lensing source (from the 1.5/Mpc of a linear run) is within 1.8e-4 of solved modes for
every model below, WDM included, and within 8.4e-5 for LCDM **[M]**
(`benchmarks/nl_followups/lin_scan.py`). The cutoff matters where halofit integrates the
continued P_L(k) for sigma(R), k_NL and the 1-halo term.

**The scale.** Each species now reports its free-streaming wavenumber
`BaseSpecies::FreeStreamingWavenumbers()`: k_fs = 1/lambda_fs, lambda_fs the distance travelled
at the rms velocity sqrt(3 P/rho) since the big bang. NCDM integrates it over the background
table; composites gather their children's; other species have none. It is dominated by the
time the species becomes non-relativistic, and it follows the WDM cutoff: m^1.2 at fixed
density, like the half-mode scale (m^1.11 in Viel et al. 2005), where the Jeans wavenumber at
the non-relativistic transition goes as m^0.66 **[M]**.

| model | k_fs (1/Mpc) |
|---|---|
| 0.06 eV, 3 x 0.2 eV, 1 eV, 3 eV neutrinos | 4.1e-4, 7.8e-4, 2.0e-3, 4.1e-3 |
| thermal WDM, a quarter of the DM: 0.1, 0.3, 1, 3, 10, 30 keV | 0.17, 0.60, 2.5, 9.5, 41, 161 |
| thermal WDM, all of the DM: 3, 10 keV | 6.2, 27 |

**The rule** (`perturb_get_k_list`, with a non-linear method only). The cutoff spans k_fs to
10 k_fs: at 0.3 keV, fitting above 11 k_fs (last mode 20/Mpc) gives 1.7e-4 in C_L^phiphi and
fitting above 5.5 k_fs (10/Mpc) 9.3e-3 **[M]**. When that span reaches the fitted modes (above
k_max / 3) and the cutoff starts before k_max_extrapolated / 2, the modes are solved up to
k_max_extrapolated. Raising the last mode only to 30 k_fs was not enough for the lightest
masses (0.1 keV: 1.1e-2 in C_L^phiphi at 5/Mpc), so it goes all the way **[M]**.

**Calibration** (halofit, defaults against modes solved to 70/Mpc without continuation,
q_logstep_limber = 1.02 in every run; `benchmarks/nl_followups/rule_scan.py`,
`verify_rule.py`):

| WDM, a quarter of the DM | C_L^phiphi before -> now | P_NL (z <= 3) before -> now | time |
|---|---|---|---|
| 0.1 keV | 1.3e-1 -> **0** | 5.7e-1 -> **0** | 0.17 -> 0.43 s |
| 0.3 keV | 4.8e-2 -> **0** | 2.8e-1 -> **0** | 0.17 -> 0.42 s |
| 1 keV | 4.4e-2 -> **0** | 5.3e-2 -> **0** | 0.17 -> 0.42 s |
| 3 keV | 8.1e-3 -> **0** | 2.2e-3 -> **0** | 0.18 -> 0.43 s |
| 7 keV | 1.1e-3 -> **0** | 6.2e-4 -> **0** | 0.17 -> 0.43 s |
| 10 keV (not triggered) | 3.9e-4 | 8.2e-4 | 0.17 s |
| 15, 20, 30 keV (not triggered) | 9.0e-5, 5.6e-5, 8.0e-5 | 9.0e-4 .. 9.4e-4 | 0.17 s |
| all of the DM, 3 and 10 keV | 1.1e-2, 6.6e-4 -> **0** | 5.3e-2, 4.9e-4 -> **0** | 0.17 -> 0.43 s |
| LCDM, 0.06 eV (not triggered) | 7.6e-5 | 5.2e-4 | 0.19 s |

"0" is identical: the run solves the same modes as the reference. HMcode behaves as halofit
(1 keV: 8.0e-3 and 7.3e-2, now 0) **[M]**. Figure:
`benchmarks/nl_followups/plots/free_streaming_cutoff.png`.

**Heavy neutrinos are not the continuation.** At 1 eV C_L^phiphi is 1.4e-3 to 1.6e-3 off
whatever the floor (2 to 60/Mpc), and P_NL at z = 3 1.5e-3 at the default, 6e-5 from 5/Mpc. At 3 eV it is
6e-3 to 1.7e-2 and not monotonic in the floor, and P_NL at z = 1 still 2.5e-3 at 60/Mpc
**[M]**: halofit's neutrino terms (Bird et al.) are calibrated for small f_nu, and at f_nu ~
0.2 its result moves with the k grid. Both are far better than before #470 (C_L^phiphi 6.6e-2
and 2.4e-1) **[M]**. The rule leaves them alone (k_fs < 0.005/Mpc).

**Comparing C_L^phiphi between runs with different last modes.** With q_logstep_limber = 1.1,
two runs whose Limber q grids differ scatter by about 1e-3 at L > 1400, alternating in sign,
and the spline between the sampled multipoles (every 40) magnifies it to 6e-3 in bands; LCDM
shows it between runs solved to 70 and to 120/Mpc **[M]**. The calibration above therefore used
1.02 everywhere. Separately, the solved reference itself moves by 1.3e-3 (LCDM) to 1.7e-2 (0.3
keV) at L ~ 3000 between 70 and 120/Mpc with halofit **[M]**: the source beyond 70/Mpc is not
negligible there, which is the choice of k_max_extrapolated = 70, not the continuation.

**Test.** `test-free-streaming-cutoff`: 1 keV is solved to 70/Mpc and its P_NL(z = 3) equals
the solved run; 10 keV, 0.06 eV neutrinos and a linear run keep their last mode near 2/Mpc. It
fails without the rule **[M]**.

## 3. HyRec under a long path

HYREC-2 reports a table it cannot open or read by printing the whole path into a 128-character
message (`hydrogen.c`) **[C]**. With any path longer than about 60 characters the fortified
`sprintf` aborted the run with SIGTRAP and no message: a missing `hyrec/` directory under a
deep site-packages, for one (exit 133 in a hand-made classy build directory, 2026-10-21) **[M]**.
Its file names are 512-character buffers, which a longer path overflows as well **[C]**.

- `HyrecModel` checks that the four tables open, with paths that fit, before HYREC-2 allocates
  anything. The error names the table and says what `hyrec_path` must be.
- `hydrogen.c`: the two messages that print paths are `SIZE_ErrorM` (2048) characters. That
  covers a table that opens but ends early, which no check outside HYREC-2 can see. It is the
  second local patch in `hyrec/PROVENANCE.md`; neither is CLASS-specific.
- `test-hyrec`: a missing directory at 200 and 600 characters and truncated tables at 200 each
  give an exception naming the path. Before, the first ended with exit 133, and so did the
  truncated tables with the check alone **[M]**.

## 4. A data race in the perturbations

`NCDMSpecies` and `DNCDMSpecies::RegisterPerturbationIndices` stored the first index of the
hierarchy in the species, which every perturbation thread shares, so all of them wrote it
concurrently **[M]** (ThreadSanitizer, the one race in a full run). Nothing read it; each thread's
`PerturbLayout` holds the indices **[C]**. Deleted.

## Checks on the whole branch

All **[M]**, 2026-10-22, on master 59064528 plus these commits:
- C++ unit tests 57 of 57 (two new, `test-hyrec` extended).
- classy against a classyref built from master (fast-math off): lvl2 1512 of 1512; the
  non-scenario tests of `test_class.py`, `test_greybody.py`, `test_hyperspherical.py` and
  `test_transfer_totals.py` pass; Cobaya 4 of 4.
- `./class` output byte-identical to master (fast-math off) for the seven configurations of
  section 1 at 1 and 11 threads; the background and transfer tables of the golden tests' cases
  bit-identical to master's. (The golden tests themselves fail on master too: their goldens
  predate later changes, and their `dncdm_dr` case is no longer valid input.)
- ThreadSanitizer, full runs with halofit, HMcode 2020 with feedback, and halofit with 1 keV WDM
  at 8 threads: no race.
- clang-format, clang-tidy (CI's checks) on the changed translation units, and doxygen
  (warnings as errors) are clean.

## Open

- The rule's two numbers (the cutoff ends at 10 k_fs; it matters below k_max_extrapolated / 2)
  are calibrated on thermal relics. A non-thermal distribution with a long high-momentum tail
  (a decay product) has the same lambda_fs integral but a different cutoff shape **[?]**.
- Halofit with heavy neutrinos (f_nu ~ 0.2) depends on the k grid at the 1e-3 to 1e-2 level
  **[M]**; not addressed.
- Dark acoustic oscillations (IDM-DRMD) also break the high-k form (#470 note); interacting
  species report no scale **[C]**.
