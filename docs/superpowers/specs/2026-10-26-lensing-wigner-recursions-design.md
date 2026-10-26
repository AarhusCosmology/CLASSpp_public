# Lensing: the Wigner recurrences run across the correlation-function nodes (2026-10-26)

**Status 2026-10-26:** implemented as described below (branch `lensing/wigner-recursions-vectorized`).
Evidence tags: [M] measured here, [C] stated by the code.

## Problem

CLASS++'s lensing module took 61 ms on one thread for LCDM lensed to l = 3000; CAMB's took 13.9 ms [M]
(`benchmarks/camb_anatomy/FINDINGS.md`). The module evaluates the lensed correlation functions at 437
nodes in theta (θ ≤ π/16 by default, `accurate_lensing = 0`). At every node it needs up to twelve Wigner
d^l_{mm'} for every l up to l_unlensed_max = 3500, by recurrence in l. An instrumented split [M]:

| piece | ms |
|---|---|
| the 12 recurrences (18.4M values) | 38.8 |
| ...of which the coefficient tables, recomputed per chunk | 1.1 |
| correlation sums (exp and the lensing terms at every l) | 16.8 |
| transform back to C_l (sampled l only) | 0.3 |

The recurrences ran one node at a time with l innermost. Each step waits for the previous one, about
7–8 cycles per value, and every value was written through a row pointer, `d[node][l]`. The 32 nodes of a
chunk are independent chains.

## Approach

* The d-functions move out of `LensingModule` into `include/wigner_d.h` / `tools/wigner_d.cpp`
  (namespace `wigner`). They are pure functions, and the classy wrapper generator turns every public
  member of a module into a Cython declaration, which cannot spell `std::array<double, kChunk>`.
* One kernel, `WignerRecurrence`, runs every recurrence for `kChunk` = 32 nodes at once, l outside and
  nodes inside, with the node state in `std::array`s. Its result is stored as `d[l * kChunk + node]`.
  Each d-function keeps its own coefficients and starting rows. A template parameter (`WignerShift`)
  spells mu, mu + b_l or mu - b_l exactly as before, so each recurrence does the arithmetic it did.
* The last chunk is padded with copies of its last node, whose results are not used.
* In `lensing_init` the sums over l (Cgl, Cgl2, the four correlation functions) and the transform back
  run over the nodes innermost too. Each node's sums still add their terms in increasing l, and each
  partial C_l adds the nodes in increasing order, as before. The expressions per element are unchanged.

## Alternatives not taken

* **CAMB's scheme** [C] (`lensing.f90`, `CorrFuncFullSkyApodized`): one Legendre recurrence per node, the
  other d-functions from it algebraically, and the lensed integrand summed at every 10th l only
  (`interp_fac`), weighted by 10. The algebraic forms need a small-angle Taylor patch for stability in CAMB.
  The subsampled sum is an approximation. With this change the 12 stable recurrences cost a few ms, so
  neither is needed.
* **Hoisting the coefficient tables** out of the chunk loop: 1.1 ms, not worth the extra state.

## Verification

* `unittests/tools/wigner_d_test.cpp`: all twelve d^l_{mm'} against d(0) = δ_{mm'} and
  d(π) = (-1)^(l-m) δ_{m,-m'}, and d^l_00 against P_l at 32 nodes, to l = 3000, with the special nodes
  at different positions of the chunk. With the sign of one recurrence's shift flipped, it fails (error
  1e6) [M].
* Against master bfe6729e [M], max |dC_l/C_l|: TT 4.5e-12, EE 8e-13, TE 8e-13, BB 7e-11 (LCDM to 3000,
  two draws); below 1.5e-13 (BB 3e-11) for Ω_k = ±0.01 and tensors; with `accurate_lensing = 1` below
  6e-9. These are summation-order effects.
* Timing [M], 1 thread, LCDM lensed to 3000 (`codecmp/stages.py`, median of 4 draws, best of 3):
  lensing 55–58 → 23 ms, whole run 0.911 → 0.878 s. One classy object over 16 draws, minimum of 3
  rounds: 1 thread −3.2% (−29.5 ms), 11 threads −2.2% (−4.4 ms; the chunks were already threaded).
  `benchmarks/camb_anatomy/dfunc_bench.cpp` puts the recurrences at 4.3 ms in this layout,
  against 39.1 ms in the old one.
* ctest 59/59, TEST_LEVEL=1 84/84, TEST_LEVEL=2 1512/1512 (as on master). The python non-scenario tests fail the same
  18 golden comparisons on master and on this branch, on this Mac.

What remains in the 23 ms is mostly the correlation sums: one `exp` per node and multipole, and the full
sum over l.
