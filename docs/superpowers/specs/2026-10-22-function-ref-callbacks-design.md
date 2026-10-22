# Callbacks carry their own context: FunctionRef in place of (function pointer, void*)

**Status — 2026-10-22. SHIPPED as described here.** `include/function_ref.h` is the
type. Every callback interface that took a `void*` context now takes a `FunctionRef`:
the five evolvers with `EvolverOptions`, `numjac`, `fzero_Newton`, `generic_integrator`,
the quadrature sampler and `fzero_ridder`. All their callers are converted. Evidence
tags: **[M]** measured here, **[C]** stated by the code, **[?]** no evidence.

## Problem

The C interfaces passed a callback as a function pointer plus a `void*` to its
context. A module method cannot be passed as a plain function pointer, so each module
kept a static wrapper per callback. The wrapper cast the `void*` to the module's
workspace struct, read a back-pointer to the module out of it, and called a
`*_member` method. The method then cast the same `void*` again to read the
workspace **[C]**:

    void PerturbationsModule::perturb_derivs(double tau, double* y, double* dy, void* p) {
      auto* pppaw = static_cast<perturb_parameters_and_workspace*>(p);
      pppaw->perturbations_module->perturb_derivs_member(tau, y, dy, p);
    }

There were twelve such wrappers in perturbations, background, thermodynamics and
primordial. In all, the code had 24 `static_cast`s from `void*`. Three workspace
structs held a module back-pointer that only the wrappers used. Nothing checked
types: passing the wrong workspace compiled and crashed at run time.

## Approach

A `FunctionRef<R(Args...)>` is that same (pointer, context) pair, with the compiler
writing the forwarding function. A call site names a lambda that captures the module
and its workspace, and passes the lambda:

    auto derivs = [this, &ppaw](double x, double* y, double* dy) { perturb_derivs(x, y, dy, &ppaw); };
    generic_evolver(derivs, tau_ini, tau_end, y, y_size, options);

The interfaces changed as follows:

* **Evolvers.** The `void* parameters_and_workspace` argument is gone from the evolvers
  and from every callback signature.
* **Module methods.** The former `*_member` methods keep their own names and take a typed
  workspace pointer. The wrappers, their casts and the back-pointers are deleted.
* **Optional callbacks.** A `FunctionRef` cannot be null, as `std::function_ref` cannot.
  An optional callback (`print_variables`, `evaluate_timescale`, `derivs_diagonal`) is a
  `std::optional<FunctionRef<...>>`, which is what it would be in C++26.
* **`output`.** The callback is required, but it is also an optional, checked by
  `EvolverOptionsCheck`. The struct is filled in after construction, and a `FunctionRef`
  member would make it impossible to default-construct.
* **Context-free callbacks.** Where a callback needs no context, a plain function still
  works: `NCDMBaseSpecies::TestFunction`, and the unit tests' free functions.

Three pieces of dead code went with it:
* `thermodynamics_recfast_timescale`: nothing referenced it.
* The no-op `thermodynamics_recfast_output_none`: a lambda replaces it.
* `HMcodeGrowthWorkspace`: the two HMcode lambdas capture what it carried.

## The contract with C++26

The interface is the subset of `std::function_ref` (P0792R14) that this code uses,
spelled the same way:
* a constructor from a function pointer, and one from a callable;
* defaulted copy;
* deleted assignment from a callable;
* `R operator()(Args...) const`;
* no default constructor and no `operator bool`.

Not provided, because nothing here needs them: the `nontype_t` constructors, and the
`const`/`noexcept` signature specialisations.

**One deliberate difference.** `FunctionRef` binds only lvalues. `std::function_ref`
also binds temporaries. That is safe for an argument, but a stored reference to a
temporary dangles, and `EvolverOptions` stores callbacks:

    options.output = [&](double x, double* y, double* dy, int i) { ... };  // dangles in C++26

That line does not compile here. Call sites name their lambdas first. The guard
goes with the swap.

**Verified against the real thing [M].** `unittests/tools/function_ref_test.cpp` builds
against `std::function_ref` when `CLASS_TEST_STD_FUNCTION_REF` is defined.
* Under GCC 16.2 (`gcc:16` image), `-std=c++26`, every check common to both passes.
* The four checks marked FunctionRef-only (temporaries rejected, two-word size) fail,
  as they should.
* With `include/function_ref.h` replaced by
  `template <class S> using FunctionRef = std::function_ref<S>;`, all 136 production
  and test translation units compile (`-fsyntax-only`).
* Only the syntax was checked: nothing was linked or run against `std::function_ref`
  beyond that unit test.

**The swap.** Replace the body of `function_ref.h` with that alias, or rename
`FunctionRef` to `std::function_ref` throughout.

## Alternatives rejected

* **`std::function`.** It owns a copy of the callable, may allocate when the capture is
  larger than its small buffer, and checks for null on every call. Ownership is the
  wrong meaning: an evolver calls its callbacks before it returns and keeps none of them.
  Measured cost: 1.20 ns per call against 1.09 ns **[M]**. That is harmless, but it buys
  nothing.
* **Templating the evolvers on the callable.** That would let the compiler inline the
  derivatives into the step loop. But the derivatives are far too large for inlining to
  pay, ndf15 (about 1500 lines) would move into a header and be compiled once per caller,
  and `auto generic_evolver = &evolver_ndf15;` (choosing the evolver at run time) would
  stop working.
* **Keeping the wrappers.** They are the cost this removes.

## Performance [M]

Microbenchmark of pure call overhead, with the evolver loop in its own translation unit
as in `tools/`. Apple clang 21, M3 Pro:

| | ns per call |
|---|---|
| function pointer + static wrapper (before) | 1.09 |
| `FunctionRef` + lambda | 1.09 |
| `std::function` | 1.20 |

The forwarder holds the callable's body inlined, and the evolver makes one indirect
call, exactly as before.

Whole runs, master against this branch, Release with fast-math. The input is LCDM with
one 0.06 eV neutrino, halofit, lensed C_l to l = 2500 and P(k) to 1 h/Mpc. Runs were
interleaved, 12 per binary:

| run | threads | master min / median (s) | branch min / median (s) |
|---|---|---|---|
| default evolvers | 1 | 0.821 / 0.838 | 0.811 / 0.844 |
| default evolvers | 11 | 0.178 / 0.191 (40 runs) | 0.176 / 0.187 (40 runs) |
| ndf15 throughout | 1 | 1.294 / 1.314 | 1.290 / 1.318 |
| ndf15 throughout | 11 | 0.243 / 0.255 | 0.243 / 0.251 |

The two binaries are indistinguishable.

## Verification [M]

* **Fast-math off: byte-identical.** Master and branch were built in one build
  directory, then run on 27 inputs. The inputs exercise:
  * every evolver on the perturbations, and etd on the background;
  * RECFAST with three thermodynamics evolvers;
  * all three inflation types (which go through `generic_integrator`);
  * HMcode 2020, shooting on `100*theta_s`, and open and closed curvature (`fzero_ridder`);
  * `k_output_values` (the print_variables path);
  * the solved BBN network;
  * automatic and manual ncdm quadrature;
  * ten `test/scenarios` inputs (greybody, self-interacting, decaying ncdm, interacting
    dark matter, scalar fields, fluid).

  All 209 output files are byte-identical.
* **Fast-math on, the default build.** 26 of the 27 inputs are byte-identical. The
  exception solves the BBN network. In master, `BbnDerivs` is a standalone function; in
  the branch, the compiler inlines it into the `FunctionRef` forwarder, and under
  `-ffast-math` the inlined copy rounds differently. The final abundances move by about
  3e-9, and the downstream adaptive solvers then take slightly different steps:
  * TT moves by 4.5e-6, EE by 3.4e-6, P(k) by 3.3e-5;
  * larger relative numbers appear only at zero crossings;
  * the same input with the tabulated BBN is byte-identical.

  This is the ULP drift that inlining under fast-math is known to produce.
* **Static checks.** Unit tests: 58 pass, one of them new. clang-tidy: no findings on
  the changed files, with the CI check set, on macOS and in the Linux CI image.
  `generate_wrapper.py` output: identical. classy: builds and runs.
* **Lifetimes: AddressSanitizer + UBSan.** A non-owning reference fails by outliving
  its callable, so the branch was also built with both sanitizers. All 27 inputs above
  run clean.
  * 57 of the 58 unit tests are clean.
  * `test-nonlinear-threads` stops on a UBSan report in `nonlinear_hmcode_compute`:
    `r_real[index_nl - 1]` with `index_nl = 0`. The same sanitized test on master
    stops at the same place, so it predates this change and is not part of it.
