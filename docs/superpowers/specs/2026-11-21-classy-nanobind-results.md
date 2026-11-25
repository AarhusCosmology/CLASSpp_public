# Results: classy's compiled core in nanobind, against Cython

**Status (2026-11-21): decided by Thomas the same day: switch to nanobind, as recommended
below. Done in PR #506 (branch `classy-nanobind`, on master after #505).** What shipped differs
from the prototype measured here in these points; everything below this header describes the
prototype:
- **No Cython left.** The `CLASS_PYTHON_BINDING` option is gone. `src/classy/_core_nb.cpp` is
  now `src/classy/_core.cpp`. `_core.pyx`, `generate_wrapper.py` and `tools/exceptions.{h,cpp}`
  are deleted, and the build requires `scikit-build-core` and `nanobind` only.
- **Docstrings:** the Cython core's docstrings were copied in full.
- **Constructors:** the default constructors that `OutputTable` and `OutputColumn` had for
  Cython are gone.
- **clang-tidy:** CI's checks pass on `_core.cpp` when run locally. That took one const fix:
  `NonlinearModule::nonlinear_pks_at_kvec_and_zvec` now takes its k and z as `const double*`.
  CI's clang-format job checks the file; its clang-tidy job configures without Python and does
  not see it.
- **From an independent review of the PR [M], all reproduced before fixing:**
  - **`NB_DOMAIN`:** each package (`classy`, `classyref`, ...) has its own nanobind internals.
    Without it, two nanobind builds in one process shared a type registry, and the second
    could not register its `PyCosmology`. CI's level-2 job would have failed as soon as master
    was nanobind. classy's exception translator also reached every other nanobind module in the
    process. Both now work: a nanobind `classy` and `classyref` side by side, and a foreign
    nanobind module that keeps its RuntimeError after `import classy`.
  - **Inputs as before:** None is accepted where Cython accepted it (`output_format`,
    `nonlinear`, `derivatives`, ...). The int arguments of the public methods take floats,
    truncated through `__int__` (`raw_cl(2.5e3)`). A non-str parameter key is an
    AttributeError, a negative `lmax` in `raw_cl_no_copy` a ValueError, and empty Bessel inputs
    an IndexError, not a segfault. `compute()` with an unknown level builds nothing.
  - **No special case for `nb::cast_error`,** which is `std::bad_cast`. It had also caught
    CLASS's own failed `dynamic_cast`s, which are now CosmoSevereError, like any other error.
  - **`get_pk_array` checks `k_size` and `z_size` against its arrays,** which Cython read past.
  - **`fc_` is destroyed after the Cosmology** that keeps a reference to it.
- **Known differences that stay:**
  - returned arrays are views (`owndata` False);
  - read-only, non-contiguous and other-dtype input arrays are accepted (converted where
    needed);
  - calling `__init__` a second time on a live object is refused;
  - `inspect.signature` does not work on the core's methods;
  - the messages of TypeErrors for wrong argument types differ.
- **2026-11-25, PR #509:** the parameter dict is in `_class.py`, and the core takes `(name, value)`
  strings. A non-str parameter name is now a TypeError from `set()`, before anything changes,
  not an AttributeError at the next build. `raw_cl_no_copy` and `cl_output_no_copy` are gone.
  The arguments are typed, which supersedes "Inputs as before" below. `lmax` is an int and the
  sizes are `size_t`, so a float or a negative size is a TypeError (`raw_cl(2.5e3)` included;
  MontePython never passes one), and NumPy integers are taken. A size beyond its array is a
  ValueError that says so, and `output_format` is 'class' or 'camb'. The public flags
  (`get_pk_array`'s `nonlinear`, ...) are truth values read in Python; the core's are `bool`.
- **Checks on the PR [M]:**
  - 1668 calls are identical to master's Cython core;
  - ctest, the level-1 scenarios, the other Python suites (the goldens included) and
    `test_cobaya.py` give the same results as on master.

Phase 2 of `2026-11-20-classy-package-design.md`. Every number below was measured on the
prototype [M], on an M3 Pro on AC power, unless tagged otherwise.

## What was built

- `src/classy/_core_nb.cpp`: the `_core` interface of `_core.pyx` in nanobind 3.1.0, with the
  same class name, method names, argument names and defaults, return types and exceptions.
  The plain-Python package (`_class.py`, `_bessel.py`, `_errors.py`) runs unchanged on it.
- A CMake option `CLASS_PYTHON_BINDING=cython|nanobind`, default `cython`, that builds one or
  the other as `<pkg>/_core`. `pyproject.toml` now also requires `nanobind>=3.1` at build
  time; Cython builds do not use it.
- **Exceptions** are mapped as `_core.pyx` maps them: by the name of the dynamic type,
  `std::invalid_argument` is CosmoSevereError, `std::runtime_error` is CosmoComputationError,
  and anything else is CosmoSevereError. The message is passed as bytes, as Cython passes a
  `std::string`. The Bessel functions keep Cython's plain `except +` mapping (ValueError,
  RuntimeError, ...), as `_core.pyx` declares them.
- **Arrays without a copy.** A table or a vector returned by the C++ code moves into a capsule
  that owns it, and the numpy arrays view it. Cython copied every column once.

Where the nanobind core behaves differently, none of it caught by the tests:
- The arrays it returns are views: `owndata` is False, `resize()` fails, and a table's memory
  is freed when the last of its columns is.
- It accepts read-only arrays, which Cython's `double[::1]` refused.
- Some docstrings are shorter than Cython's. They would have to be copied in full before a
  merge.

## Same output

- **A/B against the Cython core built from the same tree:** the 151 calls of
  `benchmarks/classy_glue/compare_tables.py` in 11 cosmologies, and 7 calls of the Bessel
  functions, error cases included. All 1668 are identical: values, key order and exception
  classes.
- **Tests, the same on both cores:**
  - ctest 73/73 (C++ only);
  - TEST_LEVEL=1 scenarios: 84 passed;
  - `test_class.py` without the scenarios, `test_package.py`, `test_greybody.py`,
    `test_hyperspherical.py`, `test_transfer_totals.py`: 240 passed and the same 3 failed on
    both. The 3 compare against a stale local `classyref`.
  - `test_cobaya.py` with `CLASSPP_REQUIRE_COBAYA=1`: 4 passed.
- **Wheels, side by side:**
  - a default Cython `classy` wheel and a nanobind wheel renamed as CI renames its reference
    (`classyref`, distribution `classy-community-ref`) were installed into one fresh venv;
  - in one process both import, run HyRec from their own data, and give identical `z_rec`
    and C_l;
  - each raises its own exception classes;
  - uninstalling `classyref` leaves `classy` working.
- **GCC:** GCC 13.3 on Ubuntu 24.04 (arm64, Docker) builds the core with no warning from
  `_core_nb.cpp`, and a run with ad + cdi, tables, C_l and an error behaves as on macOS.

## Measured

| | Cython | nanobind |
|---|---|---|
| lines of binding code | 1401: `_core.pyx` 948 + `generate_wrapper.py` 453 | 1066: one C++ file |
| the same, without blank and comment lines | 1075: 718 + 357 | 965 |
| generated on every build | `cclassy.pxi` (793 lines), `_core.cpp` (51,000 lines) | nothing |
| clean build of the binding, 1 job | 5.2 s | 10.8 s |
| clean build of the binding, 8 jobs | 5.2 s | 3.5 s |
| rebuild after editing the core | 4.8 s | 1.9 s |
| extension size, both stripped with `strip -x` | 1.92 MB | 1.66 MB |
| `import classy` | 75 ms | 78 ms |

- **Build times:** Cython's build is a chain (generate the declarations, cythonize, compile one
  file), so more jobs do not help it. nanobind compiles its own runtime, about ten files, in
  parallel with the core.
- **Line counts:** both count docstrings, and nanobind's are shorter. clang-format puts one
  argument per line when a call wraps, which lengthens the C++ file.
- **Import time:** most of it is numpy, in both.

**Timing**, best of 4 × 15, each build in its own process, alternating
(`benchmarks/classy_glue/timing_core.py`, a local harness outside the repository), ms:

| call | Cython | nanobind | ratio |
|---|---|---|---|
| `Hubble(z)`, 1700 redshifts | 0.440 | 0.443 | 1.01 |
| `angular_distance(z)`, 1700 redshifts | 0.439 | 0.443 | 1.01 |
| `luminosity_distance(z)`, 1700 redshifts | 0.438 | 0.443 | 1.01 |
| `z_of_r`, 1700 redshifts | 0.350 | 0.351 | 1.00 |
| `pk(0.1, z)`, 1700 redshifts | 4.22 | 4.22 | 1.00 |
| `get_pk`, 50 × 20 × 1 | 2.46 | 2.45 | 1.00 |
| `h()`, 1000 calls | 0.128 | 0.093 | 0.72 |
| `sigma8()`, 1000 calls | 0.190 | 0.142 | 0.75 |
| `get_current_derived_parameters`, 15 Cobaya names | 0.003 | 0.002 | 0.71 |
| `raw_cl()` | 0.444 | 0.436 | 0.98 |
| `get_background()` | 0.395 | 0.373 | 0.94 |
| `get_thermodynamics()` | 1.445 | 1.355 | 0.94 |
| `get_perturbations()` | 0.029 | 0.019 | 0.64 |
| `get_transfer(0)` | 0.040 | 0.037 | 0.92 |

- **The per-z calls are the same:** within 1%, about 2 ns in a call of 260 ns. That is the
  anti-goal of the spec, met.
- **A call through the core costs about 35 ns less**, which is the gain on `h()` and the
  derived parameters.
- **The tables gain the copy** that nanobind does not make.

## What changes besides the numbers

- **Nothing is generated from the headers.** `generate_wrapper.py` turns the C++ declarations
  into Cython ones, and needs hand-written declarations for what it cannot parse (the output
  table, the exception handler). With nanobind, a renamed or changed C++ method is a compile
  error in `_core_nb.cpp`, from the compiler that builds the rest, and the file can be checked
  by clang-format and clang-tidy like the other C++. Removing Cython would also remove
  `generate_wrapper.py`, `tools/exceptions.{h,cpp}` (used only by `_core.pyx`) and the struct
  harvesting.
- **A new build dependency:** nanobind, from pip; its runtime is compiled into the module.
- **The core is C++** where it was Cython. Students edit the plain-Python layer, which is
  unchanged; the core is one method per C++ call either way.
- **Not measured, and out of scope:**
  - nanobind can build against the stable ABI for Python 3.12 and later, one wheel per
    platform instead of one per Python version;
  - it supports free-threaded Python;
  - a method can release the GIL (`nb::call_guard<nb::gil_scoped_release>`), so that
    `compute()` could run beside other Python threads.

## Recommendation

**Switch to nanobind, in one PR that removes the Cython core.** Two cores would be the double
maintenance the package spec set out to avoid. The output is identical, and the per-z calls
are as fast. Small calls and tables are faster, the extension is smaller, and rebuilding the
core while working on it takes 1.9 s instead of 4.8 s. The case rests mostly on the toolchain:
no generated declarations, one language, and C++ checking of the binding. The costs are a pip
build dependency and a slower serial clean build. Before the switch:
- docstrings copied in full;
- Cython, `generate_wrapper.py` and `tools/exceptions.*` removed;
- CI's level-2 job keeps working while master is still Cython: its `classyref` is built from
  master, and a Cython and a nanobind package live side by side, as checked above.

This is Thomas's decision.
