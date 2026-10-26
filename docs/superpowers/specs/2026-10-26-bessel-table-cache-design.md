# Keeping the transfer module's Bessel table between cosmologies (2026-10-26)

**Status 2026-10-26:** implemented as described below (branch `transfer/bessel-table-cache`). Evidence tags:
[M] measured here, [C] stated by the code.

## Problem

The transfer module interpolates spherical Bessel functions j_l(x) from a table (`HyperInterpStruct`, the flat
"BIS") built at the start of every `TransferModule` construction: Phi and dPhi at every multipole of the l list
on a uniform x grid up to x_max = q_max tau0. For LCDM with lensed C_l to l = 3000 that is 155 multipoles times
~8000 nodes, built in 24 ms [M] (one thread, 0.85 s run). The build is serial and runs before the threaded
loop over wavenumbers, so at 11 threads it is a larger share: 23 ms of 196 ms [M].

A sampler calls `compute()` thousands of times on one classy object, and nothing in the table depends on the
cosmology except how far in x it must reach. CAMB keeps its Bessel tables in module variables for this reason
and rebuilds only when `Max_l` or `Max_eta_k` grow [C] (`cmbmain.f90`, `InitSpherBessels`). Found in the CAMB
comparison, `benchmarks/camb_anatomy/FINDINGS.md`.

## Approach

* **`BesselTableCache`** (`include/bessel_table_cache.h`) holds one `shared_ptr<const HyperInterpStruct>` and
  the parameters it was built for (l list, x_min, sampling, l_WKB, phi_min). `Get()` returns it when the
  parameters match and it reaches the requested x_max. Otherwise it builds a new one with 5% headroom in x and
  keeps that instead. A mutex guards it.
* **Ownership is explicit** (STYLE.md §1, §10: no hidden mutable state). `Cosmology` holds a
  `shared_ptr<BesselTableCache>` and hands it to `TransferModule`. classy's `PyCosmology` creates one cache and
  passes it to every `Cosmology` it builds; `reset()` replaces the `Cosmology` on every parameter change, the
  cache outlives it. A `Cosmology` constructed without one (the executable, the shooting trials) gets its own.
  When `GetInputModule()` adopts a converged shooting build, it keeps its own cache rather than the trial's.
* **A result must not depend on what was cached.** Two pieces make reuse exact:
  - *Fixed spacing.* The existing constructor fits the spacing to x_max (`(xmax - xmin)/(nx - 1)`), so every
    cosmology gets different nodes. `HyperInterpStruct::OnFixedSpacing` uses 2 pi / (beta sampling) from x_min,
    so the nodes of a table built for a larger x_max coincide with those of a smaller one. It also pads the
    node count to whole recurrence chunks (`_HYPER_CHUNK_`). Without the padding a node that falls in the
    last, partial chunk of one table and in a full chunk of another gets a dPhi that differs in the last bits.
    The unit test caught this.
  - *Decisions use the needed x_max.* The Limber switch at high q compared q with the table's last node
    (`q_max_bessel`). It now uses the x_max this cosmology needs (`bessel_x_max_`), which equals the old
    last node, so no decision changes and none depends on the table's size.
* The existing constructor is unchanged for its other callers (the per-q curved-space tables, the Python
  interface, `test_hyperspherical.py`).

## Alternatives rejected

* **Shared memory between processes (POSIX shm and the like).** Not portable to Windows. It would only save
  memory: each MPI process or walker builds its table once, in 24 ms.
* **A process-wide static cache.** Simpler, and it would cover every wrapper. But it is hidden mutable
  state, and two classy objects in different threads would contend for it.
* **Keeping the fitted spacing and rounding x_max up into buckets.** Results would then jump slightly between
  cosmologies that fall into different buckets, and would depend on which table happened to be cached.

## Verification

* `unittests/source/bessel_table_cache_test.cpp`: tables built for five larger x_max hold the reference
  table's nodes with bit-identical Phi and dPhi at every l. The cache serves covered requests and rebuilds for
  a wider x range, another l list or another sampling. With the chunk padding removed, the test fails [M].
* `test_kept_bessel_table_does_not_change_results` (python): lensed TT/EE/TE/BB/pp of a cosmology computed
  after another on the same object equal those computed fresh, exactly. The three predecessors exercise a
  reused larger table, a rebuild, and reuse through the headroom (confirmed with a temporary trace) [M].
* One-time change from the fixed spacing, vs master bfe6729e [M] (max |dC_l/C_l| / CV distance to l_max):

  | model | TT | EE | TE | BB | pp |
  |---|---|---|---|---|---|
  | LCDM + 0.06 eV, lensed to 3000 (2 draws) | 2.2e-5 / 0.001 | 2.4e-6 | 1.1e-5 | 7.6e-8 | 6.8e-7 |
  | Omega_k = +0.01 / -0.01 | 1.6e-7 | 5.9e-8 | 2.8e-8 | 2.6e-8 | 9.4e-9 |
  | tensors, r = 0.1 | 2.5e-5 / 0.001 | 3.0e-6 | 7.5e-6 | 1.5e-5 | 1.7e-6 |
  | number counts, 2 bins (density, lensing) | 4.1e-7 (all four) | | | | |

* Timing [M], one classy object over 16 LCDM draws (lensed to 3000, as cobaya runs it), master vs branch
  interleaved, minimum of 3 rounds, median of draws 1-15: 1 thread 878.5 -> 867.4 ms (-1.6%); 11 threads
  196.0 -> 171.9 ms (-12%). The first draw costs the same as before. Harness: `benchmarks/camb_anatomy/ab_cache.py`.
* ctest 59/59 and TEST_LEVEL=1 scenarios 84/84 [M]. The python non-scenario tests fail the same 18
  background/transfer golden comparisons on master and on this branch, on this Mac [M].
