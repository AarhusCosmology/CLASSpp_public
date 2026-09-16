# The LSS methods Cobaya asks classy for

**Status (2026-09-16).** Shipped, and this note records where the work diverged
from the design it started as. Shipped: the six wrapper methods below,
`NonlinearModule::GetPkGrid`, the `d_m` transfer column, the
`perturb_output_data` loop hoist, and the helper de-duplication in `classy.pyx`.
**Deliberately not shipped: `get_sources`** (section *What is not here*).

Two things the design did not foresee, both written up in place below:

* `PerturbationsModule::perturb_output_data_at_index_tau` had to be added, with
  `perturb_store_columns` factored out of the two entry points. Converting the
  topmost redshift of the grid back to a conformal time lands marginally outside
  the tabulation range and trips the guard there, so the transfer grid has to be
  read at its nodes, not at redshifts. The loop hoist stayed — it is a real fix
  — but it is no longer what makes `get_transfer_and_k_and_z` viable.
* The C++ guards in `GetPkGrid` reached Python as bare `ValueError` /
  `RuntimeError`, because the module methods harvested into `cclassy.pxd` carry
  a plain `except +`. Cobaya keys on the wrapper's own exception types, so the
  boundary restates which one it was. See *The exception types stop at the
  boundary*.

## The problem

PR #434 made Cobaya's built-in `classy` theory drive CLASS++ for CMB,
BAO-distance and SN likelihoods. Everything touching the matter power spectrum,
growth, or CLASS internals still failed on the first evaluation with an
`AttributeError` or a `TypeError` from the wrapper, because the methods Cobaya
calls were either missing or had a different signature. Issue #435 lists them.

What Cobaya actually calls, read off `cobaya/theories/classy/classy.py` (3.6.2)
rather than inferred from the requirement names **[M]**:

| Cobaya requirement | Call |
|---|---|
| `Pk_grid`, `Pk_interpolator` | `get_pk_and_k_and_z(nonlinear=, only_clustering_species=, h_units=)` |
| `Pk_grid` for `("Weyl", "Weyl")` | `get_Weyl_pk_and_k_and_z(nonlinear=, h_units=)` |
| `sigma8_z` | `sigma(8, z, h_units=True)`, once per z |
| `sigma_R` | `sigma(R, z, h_units=False)` / `sigma_cb(...)` |
| `fsigma8` | `effective_f_sigma8(z, z_step=0.1)` |
| `angular_diameter_distance_2` | `angular_distance_from_to(z1, z2)` |
| `CLASS_sources` | `get_sources()` |

The `Pk_grid` collector post-processes with `lambda P, kk, z: (kk, z, P.T)`, so
the returned grid is indexed `pk[index_k, index_z]`. That fixes the shape
contract; it is not a free choice.

## The principle this follows

A wrapper method should be a rename and a reshape. Where the Python side would
otherwise have to know how a C++ module lays out its arrays, the work belongs in
C++ instead — the module owns the layout, and it also owns the guards that say
when the layout is meaningful.

That is why the C++ diff here is four small changes rather than a list of
newly-public members, and why every new method in `classy.pyx` is a plain `def`
doing numpy over primitives that were already typed. The only Cython left in them
is what genuinely crosses the boundary: the two `vector[double]` buffers
`GetPkGrid` fills, and the reads of already-public members. Move the boundary
crossings into one private helper each and what remains is ordinary Python — which
is the point, because that is what would move into a Python module if `classy` is
ever split into a package over a thin extension. That split is not done here; see
*Rejected alternatives*.

## The four C++ changes

### 1. `NonlinearModule::GetPkGrid`

```cpp
void GetPkGrid(enum pk_outputs pk_output,
               int index_pk,
               std::vector<double>& k,
               std::vector<double>& pk) const;
```

Fills `k[index_k]` in 1/Mpc and `pk[index_tau * k.size() + index_k]` in Mpc^3
from the stored table — `ln_pk_l_` or `ln_pk_nl_`, exponentiated, on CLASS's own
k sampling. Not an interpolation: these are the nodes.

It carries three guards, so no caller can reach the arrays in a state where they
mean nothing:

* `ppt->has_pk_matter` false — no P(k) was computed at all.
* `pk_output == pk_nonlinear` while `pnl->method == nl_none` — `ln_pk_nl_` is
  empty, and reading it would be a segfault rather than an error.
* `pk_output == pk_nonlinear` while the non-linear corrections do not reach the
  top of the stored z range — halofit/HMcode ran out of k_NL, so the early end
  of the table is linear. class_public's wrapper raises here and so do we;
  silently handing back a linear spectrum labelled non-linear is the worst of
  the three outcomes.

The third one is an index comparison, not a redshift comparison:
`index_tau_min_nl_` indexes the full `tau_` grid while the stored table covers
its last `ln_tau_size_` entries, so the corrections reach the top exactly when
`index_tau_min_nl_ <= tau_size_ - ln_tau_size_`. The redshifts are recovered
from the background only to write the error message.

The first two are `class_test_severe` (`std::invalid_argument`): they depend on
which keys are in `output` and `non_linear`, never on a sampled value, and they
are API-argument validation in the sense of STYLE §8. The third is `class_stop`
(`std::runtime_error`), because it depends on `nonlinear_min_k_max` and on the
computed background. That is a deliberate deviation from class_public, which
raises `CosmoSevereError` for all three; the repo's severity rule is that a check
depending on a parsed numeric value must reject the point rather than abort the
chain.

#### The exception types stop at the boundary

Those C++ types do not arrive in Python as the wrapper's own. Only the
`Cosmology` getters, hand-declared in `classy.pyx`, carry
`except +raise_my_py_error`; every method `generate_wrapper.py` harvests into
`cclassy.pxd` carries a plain `except +`, which Cython maps to `ValueError` and
`RuntimeError`. So a module method has never been able to tell a caller to abort
rather than reject — and Cobaya, which is the caller this work exists for, keys
on exactly that distinction.

`get_pk_and_k_and_z` therefore translates at the call site, through
`reraise_as_cosmo_error`. That is a patch on one boundary, not the fix: teaching
the generator to emit `except +raise_my_py_error` for every harvested method
would fix it everywhere. That was left out of this PR deliberately — the handler
would have to become visible to the generated `.pxd`, which means either a new
shipped `.pxd` carrying an inline copy or moving the exception classes into a
Python module, i.e. the packaging change rejected below. Worth doing; worth
doing on its own.

`k` is filled alongside `pk` rather than left to the caller, although `ln_k_` and
`k_size_` are already public. The pairing of the two axes is the invariant worth
protecting, and it costs one loop.

### 2. A `d_m` column in the density transfer output

CLASS++ emits `d_tot` in the `class_format` density section but not `d_m`;
class_public emits both (`perturbations.c`, `class_store_columntitle(titles,
"d_m", ppt->has_source_delta_m)`). `d_m` is the total *matter* density contrast —
the one P(k) is built from — and the Weyl spectrum is the matter spectrum
rescaled by `((phi+psi)/2/d_m)^2`. Without the column there is no way to form
that ratio from the wrapper.

This is a parity gap, not a Weyl-specific need, so it is fixed as a parity gap:
one title and one `class_store_double`, positioned before `d_tot` to match
class_public's column order.

Blast radius, checked rather than assumed **[M]**: `has_source_delta_m_` is set
only by `has_pk_matter`, `has_nl_corrections_based_on_delta_m` or number-count
density. `python/gen_transfer_golden.py` runs `output = mTk,vTk`, so every golden
in `python/transfer_golden/` is byte-unchanged. `test_class.py`'s
`compare_output` compares Cl and P(k) only, so `COMPARE_OUTPUT_REF` against a
master `classyref` is unaffected. A run asking for both `mPk` and `mTk` does gain
a column in its `.dat` output — intended, and the class_public-matching
behaviour.

### 3. Reading the transfer functions at a node, not at a redshift

`get_transfer_and_k_and_z` wants the transfer functions on the same z grid as the
P(k) grid, and that grid *is* the `ln_tau` sampling. Asking
`perturb_output_data(z)` for them converts each z back to a conformal time, and
for the topmost node that round trip lands just below `ln_tau_[0]` and trips

```cpp
class_test(log(tau) < ln_tau_[0],
           "Asking sources at a z bigger than z_max_pk, something probably went wrong\n");
```

So the grid is read at its nodes instead, through a second entry point:

```cpp
void perturb_output_data_at_index_tau(file_format output_format,
                                      int index_tau,
                                      int number_of_titles,
                                      double* data) const;
```

which indexes `late_sources_` directly. class_public has the same pair of
functions for the same reason. Besides being the only way to reach the top node,
it is exact — no interpolation at all — and cheaper.

The two entry points differ only in how they fill the table of `T_i(k)`, so the
forty-line block that lays that table out into output columns is now
`perturb_store_columns`, called by both, rather than a second copy.

### 4. Hoisting `perturb_sources_at_tau` out of the k loop

`perturb_output_data` had

```cpp
for (index_k ...) for (index_tp ...) for (index_ic ...) {
  perturb_sources_at_tau(index_md, index_ic, index_tp, tau, pvecsources.data());
  tkfull[...] = pvecsources[index_k];
}
```

`perturb_sources_at_tau` interpolates the source for **all** k in one call. The
loop therefore did `k_size * tp_size * ic_size` full-array interpolations and
threw away all but one value each time, where `tp_size * ic_size` calls suffice —
a factor of `k_size`, around 600 in a default run.

Nobody noticed because `get_transfer(z)` is called once per run. It is still
called that way after this work — the grid goes through the node reader above —
so this is an independent fix rather than a prerequisite, and it is kept on its
own merits: it is five lines, it costs nothing, and it makes repeated
`get_transfer(z)` calls from a script affordable. Inherited from upstream, where
the loop has the same shape.

The result is arithmetically identical — the same interpolant, read at the same
points. Verified as such: every column of the `mTk,vTk` output at z=1.0 and
z=0.5 is bit-identical before and after **[M]**.

## The six wrapper methods

All in `classy.pyx`, all plain `def`.

* **`sigma(R, z, h_units=False)`** and **`sigma_cb(R, z, h_units=False)`** —
  `R` in Mpc/h when set, converted with `ba().h` before the existing
  `nonlinear_sigmas_at_z` call. `sigma(8, z, h_units=True)` is then σ8(z). The
  only change to an existing signature; it is keyword-defaulted, so every
  existing caller is unaffected.
* **`effective_f_sigma8(z, z_step=0.1)`** — the two-sided finite difference
  `d sigma8 / d ln a = -(d sigma8 / dz)(1+z)`, with class_public's step handling:
  the step shrinks to `z` for `z_step/10 < z < z_step` and the derivative goes
  one-sided below that. fσ8 from a definition, not from a second growth-rate
  code path.
* **`angular_distance_from_to(z1, z2)`** — `sin_K(chi2-chi1)/(1+z2)` with the
  three curvature cases off `ba().K`; zero for `z1 >= z2`, as Cobaya expects for
  the unordered half of a z-pair grid.
* **`get_pk_and_k_and_z(nonlinear=True, only_clustering_species=False,
  h_units=False)`** — `GetPkGrid` plus a reshape to `pk[index_k, index_z]` and
  the k rescaling. `only_clustering_species` selects `index_pk_cb_` when
  `has_pk_cb_` and `index_pk_m_` otherwise, which is the definition of
  class_public's `index_pk_cluster`.
* **`get_transfer_and_k_and_z(output_format='class', h_units=False)`** — the
  existing `get_transfer(z)` evaluated over the z grid. Raises for more than one
  initial condition, as class_public does: the return shape has no ic axis.
* **`get_Weyl_pk_and_k_and_z(nonlinear=False, h_units=False)`** —
  `pk * ((phi+psi)/2/d_m)^2 * k^4`. The `k^4` is a convention that gives the
  result a matter-spectrum-like shape, and it uses the same `k` the caller
  asked for, so the h-units convention propagates consistently.

### The z grid has one owner

Both grids are the `ln_tau` table, converted with the wrapper's existing
`z_of_tau` and with the last entry pinned to exactly 0. It is built once, in
Python, from `PerturbationsModule::ln_tau_` — already public and already in
`cclassy.pxd`.

This is load-bearing, not tidiness. `get_Weyl_pk_and_k_and_z` multiplies a P(k)
grid by a transfer grid element-wise; if the two z vectors came from separate
implementations they could disagree at the edges and the product would be
quietly wrong. `NonlinearModule::ln_tau_` is a copy of the perturbations one
(`nonlinear_get_tau_list`), and `NonlinearModule::k_[0:k_size_]` is a copy of
`PerturbationsModule::k_[scalars][0:k_size_]` (`nonlinear_get_k_list`), so one
source for each axis is also the true one. **[C]**

### Two deliberate deviations from class_public

* **The full computed k range, not `k_size_pk`.** class_public truncates its
  returned k vector at `k_max_for_pk`; CLASS++ has no `k_size_pk` member, and the
  stored table is `ln_tau_size_ * k_size_` — the whole perturbation k list. We
  return all of it. The extra high-k values are genuinely computed, not
  extrapolated (`k_size_extra_`, which *is* extrapolated, is excluded), and
  handing Cobaya's interpolator a wider grid can only reduce extrapolation. Not
  adding a member to mirror a truncation nobody needs.
* **`class_stop` rather than `class_test_severe`** for the non-linear z-range
  guard — `std::runtime_error`, so the sampler rejects the point rather than
  aborting, per STYLE §8; see above.

## What is not here

**`get_sources()` is deliberately deferred.** class_public's version builds its
output from a hardcoded list of about 25 `has_source_*` / `index_tp_*` pairs. In
CLASS++ the source index set is assembled at runtime by the species themselves,
through `BaseSpecies::RegisterTransferSourceIndices`; a hardcoded list in the
wrapper would silently omit every species-registered source, and picking species
by name in a module is the thing this codebase most consistently refuses to do.

A faithful `get_sources` therefore needs a module-side registry mapping source
index to name, of the kind `perturb_output_titles` already has for the transfer
columns but which does not exist for sources. That is its own change to the
species source interface, not wrapper glue. No built-in Cobaya likelihood
requests `CLASS_sources`; `Pk_grid`, `sigma8_z`, `fsigma8` and
`angular_diameter_distance_2` all do.

## Rejected alternatives

* **Make the Fourier arrays public and do the indexing in Cython** — the shape
  issue #435 sketches. It needs `ln_pk_l_`, `ln_pk_nl_`, `ln_tau_`, `tau_` and
  `index_tau_min_nl_` public, puts `[index_tau * k_size_ + index_k]` in the
  wrapper, and leaves the three guards with no home. One const method is a
  smaller diff *and* a smaller public surface.
* **Rebuild the z grid in Python from `z_max_pk`** — also mentioned in #435. It
  avoids exposing anything, but it reimplements `nonlinear_get_tau_list`'s choice
  of nodes in a second language, and the sampling would no longer be the one the
  table was computed on.
* **Evaluate `nonlinear_pks_at_kvec_and_zvec` at the stored nodes** — needs no
  new array access and would return the node values to spline round-trip
  accuracy. Rejected because it still needs `ln_tau_` and a non-linear-range
  accessor exposed, so it exposes *more* than `GetPkGrid` while also paying for
  interpolation it does not need.
* **Deriving `d_m` in the wrapper** instead of adding the column — it is a sum
  over the matter species weighted by their densities, so the wrapper would have
  to decide which species are matter. Same objection as `get_sources`.
* **Splitting `classy` into a Python package over a thin extension now.** The
  right long-term shape, and the reason the new methods are written as plain
  `def`. It is a packaging change, not a wrapper change: `CMakeLists.txt` install
  rules, the configurable module name that the `classyref` reference build
  depends on, `wheel.packages`, and the failure mode where a stale
  `site-packages/classy/__init__.py` shadows a freshly built `.so`. Doing it
  inside a functional PR would make both halves unreviewable.

## Verification

* `python/test_cobaya.py` — the probe likelihood grows `Pk_interpolator`
  (`delta_tot`, linear and non-linear), `sigma8_z`, `fsigma8` and
  `angular_diameter_distance_2` requirements, and each is checked against the
  wrapper's own scalar calls on the same instance: `pk()` / `pk_lin()` at sample
  (k, z), the finite-difference definition of fσ8, and
  D_A(z1,z2) = (chi2-chi1)/(1+z2) in a flat cosmology. The test pins the
  plumbing; it is not a comparison against another code, and it needs no
  likelihood data.
* `python/test_class.py`, `TestLargeScaleStructureProducts` — seventeen wrapper
  tests that do not need Cobaya. Beyond the round-trips against the scalar calls,
  three of them are the ones that would catch a plausible wrong answer rather
  than a crash:
  * the non-linear grid must exceed the linear one by more than 1.5 at high k,
    so a table that is silently linear fails rather than passing every
    shape-and-axis check;
  * `P_lin(k,z) / d_m(k,z)^2` must not depend on z, since it is the primordial
    spectrum. That pins `d_m` to the contrast P(k) is actually built from,
    without committing to a normalisation convention;
  * with a 0.3 eV neutrino, `only_clustering_species` must reproduce
    `pk_cb_lin()` and exceed the total spectrum at high k — free-streaming
    suppression, in the right direction.

**Measured** **[M]**, on this branch:

* 17/17 new wrapper tests, 4/4 `test_cobaya.py` (including the new LSS probe),
  94 passed / 6 skipped for the rest of `test_class.py`, 84 `test_scenario` at
  `TEST_LEVEL=1`, 2 `test_greybody.py`, 43/43 ctest.
* The loop hoist: every column of the `mTk,vTk` output at z=1.0 and z=0.5 is
  bit-identical before and after.
* `test_transfer_columns.py` and `test_background_columns.py` fail 18/19 on this
  machine — **and fail identically on the unmodified tree**, checked by
  rebuilding from a stash. Values drift in the eighth digit, including the
  `k (h/Mpc)` column, which nothing here touches: the known staleness of
  committed goldens across builds under `-ffast-math`. What matters for this
  work is that the *column sets and order* assertions pass, which is the direct
  evidence that `d_m` stays out of an `mTk,vTk` run.
