# N-body gauge output: computing `H_T_Nb''` for the `k2gamma_Nb` transfer

**Date:** 2026-07-06
**Status:** Approved design, ready for implementation plan
**Author:** Thomas Tram (with Claude)

## Problem

N-body gauge output is currently blocked. Requesting either

- `extra metric transfer functions: yes` (sets `has_metricpotential_transfers`), or
- `Nbody gauge transfer functions: yes` (sets `has_Nbody_gauge_transfers`)

sets `has_source_k2gamma_Nb_ = true` (`perturbations_module.cpp:900-901`, `:909-911`), which trips a
`class_stop` in `perturb_sources` (`perturbations_module.cpp:4954-4959`):

> "We need to compute the derivative of `H_T_Nb_prime` numerically. Written by T. Tram but not yet
> propagated here."

So today *any* N-body-gauge output aborts `compute()` — including the density/velocity gauge
conversion (`delta_i`, `theta_i` in N-body gauge), which physically only needs `H_T_Nb'` and was
merely collateral damage from the stop.

The scalar N-body-gauge potential `γ` (emitted as the `k2gamma_Nb` transfer column) is defined
(paper arXiv:1505.04756v2, Eq. "gamma") as

```
γ ≡ Ḧ_T + (ȧ/a) Ḣ_T − 8πG a² p Π
```

In CLASS the stored `k2gamma_Nb` currently holds only the **partial** value
(`perturbations_module.cpp:4948-4951`):

```
k2gamma_Nb (partial) = −ℋ · H_T_Nb'  +  (9/2) a² (ρ+p) σ
```

The missing term is exactly `−H_T_Nb''` — the conformal-time second derivative of the N-body-gauge
`H_T`, i.e. the derivative of the already-computed `H_T_Nb'` source. That derivative is very hard to
form analytically inside the RHS (it would require `σ'`, the derivative of the total shear, which is
not available — the known "shear' blocker"). This spec computes it numerically instead.

The reference prototype is `ComputeGamma.ipynb`, which reproduces `arXiv-1505.04756v2/gammaPhiRatio.pdf`
by post-processing transfer functions: `k2gamma = get_transfer(z)['k2gamma_Nb'] − H_T_Nb''(z)`, where
`H_T_Nb''` is a numerical `τ`-derivative of `H_T_Nb'`.

## Goal

Compute `H_T_Nb''` inside CLASS so that the stored `k2gamma_Nb` transfer column is the full,
correct `k²γ`, and remove the `class_stop`. This unblocks all N-body-gauge output
(`get_transfer`, `*_tk.dat`) for both input flags above.

## Non-goals

- **No `perturbations_k*.dat` change.** `k2gamma_Nb` appears only in the mTk transfer output
  (`perturbations_module.cpp:269`, `:328`); the `k_output_values` time-series has no `γ` column.
  Correcting the `sources_` array is therefore sufficient and complete.
- **No committed golden/regression test.** The large-`k` values of `γ` are intrinsically noisy and
  grid-sensitive (see Validation), so a golden column would be brittle. Validation is the interactive
  notebook check only.
- **No new physics or analytic `H_T_Nb''`.** We take the numerical derivative of the existing
  `H_T_Nb'` source; the RHS and the partial `k2gamma_Nb` computation are unchanged.
- **No performance concern.** Output-only, guarded; runs only when N-body-gauge transfers are
  requested.

## Approach (chosen)

**Post-process per `k`: spline the stored `H_T_Nb'(τ)` time series and read `H_T_Nb''` off the
spline.** The derivative is the analytic derivative of the cubic that interpolates the data on the
native `tau_sampling_` grid — there is **no finite-difference step `dτ` to choose**, which dissolves
the timescale-of-variation problem (the grid is already tuned to resolve the sources, so it is
denser exactly where the physics varies fast). This was chosen over two alternatives:

- **On-demand finite difference at output time** (a direct C++ port of the notebook): rejected — it
  reintroduces the `dτ` knob and, since we already have a spline of `H_T_Nb'`, differencing it is
  strictly worse than reading the derivative off it; also requires intercepting the column at every
  output site.
- **In-RHS predictor** (`y+dτ·dy`, recompute `H_T_Nb'`, difference per sample): rejected — worst
  `dτ` behaviour (finite-difference truncation *plus* Euler-linearization error), and most invasive
  (must re-run the whole total-stress-energy assembly over all species at the shifted state).

The spline helpers already exist and are used identically elsewhere (e.g. `dκ/dτ` in
`thermodynamics_module.cpp:642-651`):

- `array_spline_table_line_to_line(x, n_lines, array, n_columns, index_y, index_ddy, spline_mode)`
- `array_derive_spline_table_line_to_line(x, n_lines, array, n_columns, index_y, index_ddy, index_dy)`

Both operate on strided tables, matching the `sources_` layout.

## Design

### Where

A new block in `perturb_solve` (the per-`k` routine), placed immediately after the evolver returns
(`perturbations_module.cpp:2282-2296`) and before the zero-fill tail loop (`:2307-2312`).
`tau_actual_size` and `index_k` are both in scope there. The block is guarded by
`has_source_k2gamma_Nb_` (only ever true for scalars; `has_source_H_T_Nb_prime_` is guaranteed true
alongside it, so the `H_T_Nb'` source column is always present).

### Algorithm (per `k`)

Let `T = tau_actual_size`, `S = k_size_[index_md]`, `base_HTp` / `base_k2g` the `sources_` vectors for
`index_tp_H_T_Nb_prime_` / `index_tp_k2gamma_Nb_` at this `(index_md, index_ic)`.

1. Allocate a scratch table `work` of shape `[T][3]`, columns `(y=0, ddy=1, dy=2)`.
2. Copy the `H_T_Nb'` time series for this `k`:
   `work[i*3 + 0] = base_HTp[i*S + index_k]` for `i ∈ [0, T)`.
3. `array_spline_table_line_to_line(tau_sampling_.data(), T, work.data(), 3, 0, 1, _SPLINE_EST_DERIV_)`
   → fills column 1 (second derivatives). `_SPLINE_EST_DERIV_` estimates the edge derivatives, which
   is the right endpoint behaviour for a subsequent derivative.
4. `array_derive_spline_table_line_to_line(tau_sampling_.data(), T, work.data(), 3, 0, 1, 2)`
   → fills column 2 with `H_T_Nb''(τ_i)`, the analytic derivative of the interpolating cubic.
5. Subtract into the stored partial value:
   `base_k2g[i*S + index_k] -= work[i*3 + 2]` for `i ∈ [0, T)`.

After this, `sources_[…][index_tp_k2gamma_Nb_]` for this `k` holds the full
`k²γ = −ℋ·H_T_Nb' + (9/2)a²(ρ+p)σ − H_T_Nb''`.

The scratch table is allocated per `k` call (cheap; performance is irrelevant here). It may be
hoisted to a single allocation outside the `k` loop if convenient, but this is not required.

### Remove the block

Delete the `class_stop` at `perturbations_module.cpp:4954-4959`. Keep the partial `k2gamma_Nb`
computation at `:4948-4951` unchanged. Update the now-stale comment near `:4947` to state that the
`−H_T_Nb''` term is applied as a post-processing spline derivative in `perturb_solve`, so future
readers know where the full value is assembled.

### Data flow summary

```
RHS / perturb_sources (per τ sample, per k):
    store  H_T_Nb'(τ_i, k)            → sources_[H_T_Nb_prime]
    store  partial k2gamma_Nb(τ_i, k) → sources_[k2gamma_Nb]   (= −ℋ H_T' + shear)

perturb_solve (per k, after evolver):        [NEW]
    spline H_T_Nb'(·, k) over tau_sampling_[0..T)
    H_T_Nb''(τ_i, k) = d/dτ of that spline
    sources_[k2gamma_Nb][τ_i, k] −= H_T_Nb''(τ_i, k)

output (get_transfer / *_tk.dat):
    interpolate sources_[k2gamma_Nb] at requested z → full k²γ
```

## Validation

`ComputeGamma.ipynb` is the validation harness. After the change:

1. Run the notebook (Newtonian gauge, `extra metric transfer functions: yes`) and confirm
   `|γ|/Φ` reproduces `arXiv-1505.04756v2/gammaPhiRatio.pdf`:
   - low-`k` plateaus per redshift (e.g. `z=1000 ≈ 0.1–0.2`; `z=49` crossing the 1% line near
     `k ≈ 5×10⁻³ h/Mpc`; lower-`z` curves progressively lower),
   - the 1% reference line,
   - the falloff toward large `k`.
   The notebook currently computes `H_T_Nb''` externally (via `get_transfer` at `z±0.1`); after the
   change, `get_transfer(z)['k2gamma_Nb']` alone should already be the full `γ`, so the notebook's
   external correction term should become redundant (a good consistency check: the in-code value and
   the notebook's externally-corrected value agree to plot tolerance).
2. Large-`k` oscillatory noise is expected and acceptable — it is the intrinsic error in `H_T_Nb''`.
   The spline-over-source-grid derivative is smoother than the original integrator-timestep finite
   difference used to make `gammaPhiRatio.pdf`, so we expect equal-or-less noise. Acceptance:
   **no worse than `gammaPhiRatio.pdf`**; low-`k` plateaus and the 1% crossings match.
3. Confirm `Nbody gauge transfer functions: yes` no longer aborts and produces sensible
   N-body-gauge density/velocity transfers.

No committed golden test (per decision (b)).

## Risks & mitigations

- **Source-grid density.** `H_T_Nb''` accuracy is capped by `tau_sampling_`. If that grid
  under-resolves the large-`k` oscillations relative to the original integrator-timestep method, the
  large-`k` tail could *differ* (aliasing) rather than merely look noisier. Expectation: the grid is
  tuned to resolve CMB sources at the same acoustic frequencies, so it is adequate; the notebook/plot
  comparison is what confirms it. Fallback if inadequate: increase source sampling for this output,
  or fall back to an integrator-timestep-based derivative.
- **Endpoints.** `array_derive_spline_table_line_to_line` uses one-sided spline formulas at `τ_ini`
  and `τ_0`; both ends are benign for `H_T_Nb'` (smooth, small deep in radiation era; `z=0` at the
  other end). `_SPLINE_EST_DERIV_` improves edge behaviour.
- **Integrated vs zero-filled region.** The spline must run over `[0, tau_actual_size)` only; the
  tail `[tau_actual_size, tau_size_)` is zero-filled afterwards and must not enter the derivative.

## Files touched

- `source/perturbations_module.cpp`
  - `perturb_solve`: add the per-`k` spline post-processing block (after the evolver, before the
    zero-fill tail).
  - `perturb_sources`: remove the `class_stop` (`:4954-4959`); update the comment near `:4947`.

## References

- `arXiv-1505.04756v2/finalpaper.tex` — Eq. "gamma" (line 170-173), γ definition and the large-scale
  cancellation note (line 326).
- `arXiv-1505.04756v2/gammaPhiRatio.pdf` — validation target.
- `ComputeGamma.ipynb` — prototype / validation harness.
- `tools/arrays.cpp:49` (`array_derive_spline_table_line_to_line`), `:223`
  (`array_spline_table_line_to_line`); `include/arrays.h:16-17` (`_SPLINE_NATURAL_`,
  `_SPLINE_EST_DERIV_`).
- `source/thermodynamics_module.cpp:642-651` — existing `dκ/dτ` spline-derivative usage to mirror.
```
