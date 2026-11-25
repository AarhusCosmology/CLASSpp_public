# Design: classy as a package, with a thin compiled core and a plain-Python `Class`; then a nanobind prototype of the core

**Status (2026-11-20): approved by Thomas; phase 1 implemented on branch `classy-package`
(PR #504, after #503, which merged as 67e51a6b). Phase 2: nanobind replaced Cython, decided by
Thomas on 2026-11-21 (PR #506), with the numbers in `2026-11-21-classy-nanobind-results.md`.
Everything below about Cython, `cclassy.pxi` and `generate_wrapper.py` is superseded by it.**
**2026-11-25 (PR #509): the parameter dict moved from the core to `_class.py`.** `set`, `empty`, `pars`, `state`,
`struct_cleanup` and `compute` are Python; the core takes the parameters as `(name, value)` strings
(`_set_input`) and builds a level by name (`_build`). `get_pk_array` and `get_pk_cb_array` are Python
over the core's `_pk_arrays`. `raw_cl_no_copy` is gone. Section 2's table is superseded on these rows.
Where the
implementation differs from the text below:
- A module's derived-parameter name keeps the module's meaning when a species answers it too
  (section 3, found while building: NEDE's `z_star`).
- `get_perturbations` and `transfer_general` still branch on the input struct in the core.
  The output-table PR agreed on 2026-11-20 (an insertion-ordered table replacing the title
  strings and NamedColumns) reworks both, so they were left for it.
- `_z_grid_for_output` stays in the core: it is a single C++ call.
- `generate_wrapper.py` keeps harvesting the structs: the core still reads `perturbs` and
  `nonlinear` (section 4).
- Two crashes that master has too were fixed along the way: `nonlinear_scale` without a
  nonlinear method, and `raw_cl_no_copy` with several initial conditions. After the review,
  also `GetInputModule` (the shooting) throwing through Cython, which aborted the interpreter.
- The names read from the input (`h`, `H0`, `T_cmb`, `Omega_b`, `omega_b`, `Omega_g`,
  `Omega0_k`) are a table in `Cosmology` that builds no module, as classy read them from the
  input struct before. Found in review: the first version built the background for `h()`.
- The test that a species name needs no other code uses the pheno-axion fluid's `a_c`, which
  no table lists, rather than a test-only species.
- `Omega0_cdm`, a getter before, is answered by CDM and has a placeholder (0), as the getter had.
- The core still makes three small decisions: the default P(k) output (`_pk_output_default`),
  aliasing P_cb to P_m when there is no separate cdm+baryon spectrum, and `lmax = -1` meaning
  the largest l. `PyCosmology` stays importable from `classy`, as it was.

Builds on PR #503 (branch
`classy-glue-to-cpp`), which moved the table and value-at-z glue out of classy.pyx into the
modules. Phase 1 below is step 2 of that work; phase 2 is a prototype whose results go in a
companion `-results.md`.

## Why

The aim is the thinnest possible compiled layer between the C++ modules and Python, for two
reasons.

1. **Where code belongs.** Physics, checks and rules belong in C++, where they are compiled,
   linted (clang-tidy, clang-format), unit-tested and shared with the `class` binary.
   Composition, legacy shims and input normalisation belong in plain Python, which students
   read, which needs no rebuild, and which gives readable tracebacks. classy.pyx sits between
   the two and holds both kinds.
2. **To decide nanobind against Cython with numbers.** The comparison only means something
   once the compiled layer does nothing but conversion; otherwise it measures how much logic
   each tool happens to carry.

After #503, classy.pyx has 1,787 lines. About 30% of its code is plain Python that only lives
in a .pyx file:
- the derived-parameter table (`get_current_derived_parameters`, about 200 lines);
- the reshaping of the P(k,z) and T(k,z) grids, `effective_f_sigma8`, `angular_distance_from_to`;
- the CosmoHammer `__call__`, the slow-roll formulas and the Bessel input checks.

lesgourg/class#274 (reviewed November 2026) generates a Cython layer that exposes every struct
field and function, with a plain-Python classy on top. Its plain-Python half is taken over here.
The automatic exposure is not (see Anti-goals).

## Goals

- classy installs as a package: `classy/__init__.py`, the compiled `classy/_core`, and the
  data in `classy/bbn/` and `classy/hyrec/`.
- `_core` contains only:
  - (a) the `Cosmology` lifetime and the parameter gate (`_cosmo()`, `set`, `empty`);
  - (b) one method per C++ call, converting types (float, numpy array, dict);
  - (c) exception translation;
  - (d) loops that evaluate one C++ function over an input array.

  It holds no physics, no unit conversions, no checks, no name tables.
- Physics and checks move to C++, composition to plain Python.
- The public Python API is unchanged, except:
  - **`get_input_*` (10 methods) is removed.** Nobody uses it (Thomas, 2026-11-20). It is also
    the only thing that needs Cython's automatic struct-to-dict conversion.
  - **`omega_m` is fixed.** It returns Omega_m/h² instead of Omega_m·h² on master (0.695
    instead of 0.140 for Omega_m = 0.312, h = 0.67).
- Everything else gives the same output as master, bit for bit where the code path is unchanged.

## Non-goals

- No change to the numbers CLASS computes.
- No automatic exposure of C++ members.
- No change to the method names, signatures or exception classes that Cobaya, MontePython,
  CosmoHammer and the tests use.
- Phase 1 does not decide on nanobind. Phase 2 measures; Thomas decides.

## Anti-goals: failures even if everything else works

- A Python layer that reads C++ internals by name at run time, for example a dict of fields.
  A renamed C++ member would then break at run time instead of at build time.
- Per-z methods (`Hubble`, `angular_distance`, `luminosity_distance`, `z_of_r`, `pk`) slower
  than on master by more than the noise. Cobaya calls them about a thousand times per point.
- `classy` and `classyref` unable to live in one environment, or uninstalling one breaking the
  other.
- A source tree in which `import classy` picks up an uncompiled source package. This would
  happen to anyone who runs Python at the repository root, and to pytest in `python/`.

## Constraints

- **Build:** scikit-build-core drives CMake. The CI reference build renames the installed
  package with
  `--config-settings=cmake.define.CLASS_PYTHON_MODULE_NAME=classyref` (and renames the
  distribution in pyproject.toml with sed). That command line must keep working, because CI
  runs it on master.
- **Generated declarations:** `cclassy.pxi` stays generated by `generate_wrapper.py` for the
  Cython core.
- **MontePython ≥ 3.6** expects `python/build/lib.<version>.<platform>-<tag>/`
  (`scripts/montepython_layout.py`).
- **What callers import:**
  - the tests: `from classy import Class, CosmoSevereError, CosmoComputationError`, and the
    hyperspherical Bessel functions from `classy`;
  - CosmoHammer: `Class.__call__`;
  - Cobaya: `classy.__version__`.
- **Data paths:** the HyRec data must be found by `classy` and by `classyref`. A missing
  `hyrec/` directory ends in a SIGTRAP inside HyRec, not in an error.
- **STYLE.md error conventions:** `class_test` (std::runtime_error) becomes
  CosmoComputationError; `class_test_severe` (std::invalid_argument) becomes CosmoSevereError.
  Modules do not pick species by type name. No POSIX math constants.

## Approaches considered

**1. How the Python layer reads C++ numbers** (`z_rec`, `rs_d`, `n_s`, ...).

| | Approach | Trade-off |
|---|---|---|
| **A** | Each module answers named derived parameters; one forwarder | Names live in C++; the nanobind port stays small. **Chosen by Thomas, 2026-11-20.** |
| B | One hand-written Cython getter per number (about 45) | No C++ change, but it is getter boilerplate, and nanobind would have to repeat it. |
| C | The generator exposes every public member to Python | Least code, but #274's approach: a rename breaks at run time. |

*My call within A:* the lookup is per name, not one whole map per module. This has two
consequences:
- The module that refuses a name gives its own reason: `100*theta_d` without
  `compute damping scale`, BBN abundances without `YHe = network`, `sigma8` without a matter
  spectrum.
- Only the module that owns a name is computed. Asking for `h` must not compute the
  thermodynamics, as today.

**2. How the Python `Class` relates to the core.** *My call.*
- **(i) Chosen: `Class` subclasses the compiled core.** Forwarders that need no logic are
  inherited, so a per-z call stays a single call.
- **(ii) Rejected: `Class` wraps a core object.** Every call would go through a Python method,
  about 70 ns more each, which is +25% on `Hubble(z)` (now about 265 ns per call).

**3. Where the package source lives.** *My call.*
- **Rejected: `classy/` at the repository root.** It shadows the installed package for anyone
  running Python at the root; `scripts/montepython_layout.py` already strips the root from
  `sys.path` against a similar hack.
- **Rejected: `python/classy/`.** pytest prepends `python/` to `sys.path`, so the tests would
  import the uncompiled source.
- **Chosen: `src/classy/`.** It is never on `sys.path`.

## Design, phase 1

### 1. Layout and build

- **Sources:**
  - `src/classy/__init__.py` and any further plain-Python modules;
  - `src/classy/_core.pyx`, which replaces the root `classy.pyx`;
  - `generate_wrapper.py` stays at the root and still writes `cclassy.pxi` into the build
    directory.
- **Installed layout, with `<pkg>` = `CLASS_PYTHON_MODULE_NAME` (default `classy`):**
  - `<pkg>/__init__.py` (and any other .py files);
  - `<pkg>/_core.<soabi>.so`;
  - `<pkg>/data/bbn/*.dat` and `<pkg>/data/hyrec/*.dat`.
- **Module naming:** the Cython module is `_core` (`PyInit__core`) whatever the package is
  called, so the rename only changes the install destination. The CMake target needs a name
  that does not clash (for example `classy_core`, with OUTPUT_NAME `_core`).
- **Data:** each package carries its own data, so `classy` and `classyref` no longer share
  `site-packages/classy/` as they do today. The default `class_dir` becomes `<pkg>/data`.
  *My call:* the `data/` level keeps the new files off the paths that every old-layout build
  writes, whatever its import name (`site-packages/classy/{bbn,hyrec}`). A new `classy` and an
  old `classyref` or `classyold` in one environment then never own the same file, so
  uninstalling one cannot delete the other's data. That covers CI's level-2 job, which installs
  master's `classyref` next to the branch.
- **The default `class_dir` applies whenever the parameters lack one.** *My call*, a bundled
  fix. On master, `__init__` puts it into `_pars`, but `empty()` and the CosmoHammer `__call__`
  replace `_pars` wholesale. The next run then falls back to `__CLASSDIR__`, the source tree the
  build was compiled from. A wheel built elsewhere misses `hyrec/` there and traps. Cobaya's
  `close()` calls `empty()`, and the ClientCluster kit works around it by keeping the build tree
  in place. The gate adds the default when it builds the `FileContent`; `_pars` and `pars`
  behave as on master.
- **Exceptions** (`CosmoError`, `CosmoSevereError`, `CosmoComputationError`) are defined once
  and re-exported by `classy/__init__.py`. Where to define them is the builder's call, under
  one constraint: `except classy.CosmoSevereError` must catch what the core raises, also in
  phase 2's nanobind core. A plain-Python `_errors.py` imported by the core is the natural
  place, if the core's relative import works under the `classyref` rename. Otherwise the core
  defines them and the package re-exports.
- **`__version__`** stays importable as `classy.__version__`.
- **Upgrading an existing install:**
  - pip removes the old top-level `classy.<soabi>.so`, which is listed in the old RECORD.
  - A stale one left behind loses to the new package, because a directory with `__init__.py`
    is found before a module of the same name.
  - Check this once by hand; see Verification.
- **MontePython:** `scripts/montepython_layout.py` copies the package directory into
  `lib.*/classy` (today it copies the module file and the data directory).
- **ClientCluster kit.** *My call.* `grendel/install/09_classpp_guard.sh` and
  `10_classpp_build.sh` stop with FATAL unless `lib/classy*.so` exists, and a build with this
  layout installs `lib/classy/_core*.so` instead. A small ClientCluster PR adapts them when
  phase 1 merges; Thomas runs the sync on Grendel. The CLASS++ PR says so.
- **Dev builds:** a plain CMake build (`cmake --build build/dev --target ...`) should leave a
  runnable package in the build directory: the `.py` files and data copied or linked next to
  `_core`, so that `PYTHONPATH=build/dev` works. This replaces the manual symlinks of `hyrec/`
  and `bbn/`. The Makefile comments follow.
- **sdist:** `src/` is included by default. Check that the root `classy.pyx` is gone and
  nothing excludes `src/`.

### 2. What goes where

The rule: a method stays in `_core` if it is a single C++ call plus conversion, or a loop of
one C++ function over an input array. Logic moves to C++ if it is physics or a check, and to
Python if it is composition.

| Today in classy.pyx | Goes to | Note |
|---|---|---|
| `__init__`, `reset`, `_cosmo`, `set`, `empty`, `pars`, `state`, `struct_cleanup`, `compute` | core | the gate |
| C_l methods (`raw_cl`, `lensed_cl`, `raw_cl_no_copy`, `lensed_cl_computed`) | core | conversion |
| per-z methods: `Hubble`, `angular_distance`, `luminosity_distance`, growth factors, `Om_m`, `z_of_tau`, `ionization_fraction`, `baryon_temperature`, `z_of_r` | core | hot; inherited by `Class` |
| `pk`, `pk_cb`, `pk_lin`, `pk_cb_lin`, `get_pk*` loops, `get_pk_array*`, `pk_tilt`, `nonlinear_scale*` | core | single calls or array loops |
| tables: `get_background`, `get_thermodynamics`, `get_primordial`, `get_perturbations`, transfer tables | core | conversion |
| `sigma`, `sigma_cb`: the two checks (`has_pk_matter`; `k_max_for_pk` ≥ h) | C++ | a public NonlinearModule method that refuses; the internal `nonlinear_sigmas_at_z` stays unchecked, because the module calls it for sigma8 without mPk |
| `sigma`, `sigma_cb`: `h_units` | Python | unit convention |
| `get_current_derived_parameters` table, `get_slowroll_parameters`, `r_0002`, the `_h`-scaled distances, `exp_m_2_tau_As`, `omega_m`, the species dot syntax | C++ | see section 3 |
| `get_current_derived_parameters` itself | Python | type check, then one lookup per name |
| scalar getters: `age`, `h`, `n_s`, `A_s`, `tau_reio`, `Omega_m`, `Omega_r`, `theta_*_100`, `Omega_Lambda`, `Omega_g`, `Omega_b`, `omega_b`, `Neff`, `k_eq`, `sigma8*`, `rs_drag`, `z_reio`, `T_cmb`, `Omega0_k`, `Omega0_cdm`, `f_idr_drmd`, `z_dec_drmd`, `Omega_nu` | Python, via the derived-parameter lookup | one per point, not hot; each becomes a derived-parameter name if it is not one already |
| `angular_distance_from_to` | C++ (BackgroundModule) | physics; `np.sin`/`np.sinh` become `sin`/`sinh`, so ULP-level differences are allowed and listed in the PR |
| `get_pk_and_k_and_z`, `get_transfer_and_k_and_z`, `get_Weyl_pk_and_k_and_z`, `effective_f_sigma8`, `_z_grid_for_output` | Python | composition over core forwarders |
| `__call__` (CosmoHammer) | Python | through `empty()` and `set()`; same semantics |
| hyperspherical Bessel functions | Python wrappers, with a core forwarder per C++ entry point | `_normalize_bessel_inputs` and the shape handling are Python |
| `get_input_*` | removed | see Goals |

### 3. Derived parameters in C++

- **Entry point:** one call, `double Cosmology::GetDerivedParameter(const std::string& name)`.
  It is non-const, because it builds modules on demand.
- **Routing:**
  - BackgroundModule, ThermodynamicsModule, PrimordialModule and NonlinearModule each answer
    their own names, and can say statically which names they own, so that Cosmology builds
    only the owner.
  - Names that combine modules (`exp_m_2_tau_As`) belong to Cosmology.
  - So does the `_h` suffix on the six distances (`rs_rec_h`, ...), which multiplies by h.
- **Species answer their own names** (Thomas, 2026-11-20). This uses the virtual that already
  exists, `std::optional<double> BaseSpecies::GetParam(const std::string&)`. Today only the dot
  syntax calls it. AxionEDEFluid, IDM_DR_IDR, NEDE, scalar-tensor and greybody NCDM already
  override it. No new virtual is needed.
  - **Plain names** (`N_dg`): if a module owns the name, the module answers and the species
    are not asked. Otherwise BackgroundModule asks every species. A module's name keeps its
    meaning whatever the species say: NEDE answers `z_star` (its transition) as well as the
    thermodynamics (recombination), and a plain `z_star` stays the recombination value, as on
    master. NEDE's is `<instance>.z_star`. *Found while building, 2026-11-20; this replaces
    the earlier rule that a species reusing a module's name is an error, which would have
    broken `z_star` in every NEDE run.*
  - **Composites** answer their own coupling names and pass the rest to their children. The
    default `CompositeSpecies::GetParam` loops over `children_`, as the composite already does
    for Omega0. An IDR inside IDM_DR_IDR therefore answers `N_dg` exactly as a standalone IDR
    does.
  - **Exactly one answer** gives the value.
  - **Two or more species answering** is a `class_test_severe` that names them and suggests
    `instance.name`. This covers two instances of one type and two species types that claim
    the same name.
  - **No answer** gives the placeholder below if the name is a legacy one. Otherwise it is the
    "not recognized" severe error. A new name asked for in a run without its species is
    therefore an error, not 0.
  - **Dot syntax** (`axion1.f_peak`): BackgroundModule asks that instance only, first through
    `GetParam`, then through the generic `f_peak`/`a_peak`/`z_peak` from the table. The
    `_SPECIES_DERIVED_FIELDS` whitelist goes: any name the instance answers works.
  - **Adding a derived parameter to a new species** means overriding `GetParam` in that
    species. No edit in BackgroundModule, Cosmology, `_core` or the Python layer.
- **Legacy names that move to species,** with the placeholder that master returns when the
  species is absent. The placeholders are one table in Cosmology, with a comment saying why
  they exist: legacy callers ask for these names in runs without the species.

  | Names | Answered by | Placeholder |
  |---|---|---|
  | `Omega_Lambda`, `Omega0_lambda` | Lambda | 0 |
  | `Omega0_fld` | Fluid | 0 |
  | `N_dg`, `xi_idr` | IDR | 0 |
  | `a_dark`, `Gamma_0_nadm` | IDM_DR_IDR (coupling) | 0 |
  | `f_idr_drmd` | IDM_DRMD_IDR_DRMD | 0 |
  | `z_dec_drmd` | IDM_DRMD_IDR_DRMD | -1 |

  The `static_cast`s in `BackgroundModule::f_idr_drmd()` and `z_dec_drmd()` go. They are the
  species picking that STYLE.md forbids.

  *Checked while building:* AxionEDEFluid inherits `FluidSpecies::GetParam`, but the
  pheno-axion fluid is keyed `Fluid` like the plain fluid, so a run has one fluid instance and
  `Omega0_fld` one answer, as on master.
- **Sums over species stay in the modules** (`m_ncdm_tot`, `Neff`, `Omega_m`, `Omega_g`). A
  loop that adds up the species is allowed; a species cannot answer for the others.
- **Physics moved into C++:**
  - the slow-roll parameters (`epsilon_V`, `eta_V`, `ksi_V^2`) and `r_0002` go to PrimordialModule;
  - `omega_m` is computed as Omega_m·h².
- **Errors keep their classes:**
  - an unknown name is `class_test_severe` ("... was not recognized as a derived parameter");
  - a name this run cannot answer is severe, with its reason.

  Messages gain the usual C++ prefix. Tests that match on text keep the phrase they match.
- **The Python side** is the type check (`TypeError` for a non-list) and
  `{name: self._derived_parameter(name)}`.
- **Values** are what master gives, apart from `omega_m`. Where master combined two numbers in
  Python (theta_s, ln10^{10}A_s, ...), the same arithmetic in C++ may differ by an ULP under
  -ffast-math. Allowed, and listed if it happens.

### 4. Removed

- `get_input_*`.
- In `generate_wrapper.py`, the struct harvesting may become unnecessary once the core reads
  no input struct. Remove it only if nothing in the core still needs it. The builder checks.

### 5. Verification

- **A/B against master:** extend `benchmarks/classy_glue/compare.py` from the methods #503
  touched to every public method of `Class`, and compare with a master build (`classyold`).
  Everything must be identical, except:
  - `omega_m` (fixed);
  - `get_input_*` (removed);
  - ULP-level differences where arithmetic moved from Python to C++, each named.
- **Test suites:**
  - ctest;
  - `test_class.py` without the scenarios, and TEST_LEVEL=1 scenarios;
  - `test_greybody.py`, `test_hyperspherical.py`, `test_transfer_totals.py`, `test_cobaya.py`;
  - CI's level-2 comparison.
- **New tests:**
  - `omega_m` = Omega_m·h²;
  - a species name answered by two instances is severe, and `instance.name` resolves it;
  - a legacy name in a run without its species gives the placeholder;
  - a ctest with a test species that overrides `GetParam`: its name is reachable through
    `Cosmology::GetDerivedParameter` with no other code change;
  - every scalar getter equals the derived parameter of the same name;
  - severe errors for an unknown derived name and for names the run cannot answer;
  - `Class` is an instance of the core class (inheritance).
- **Packaging, checked by hand in a scratch venv and reported in the PR:**
  - `classy` and `classyref` installed side by side: both import, both run HyRec, and
    uninstalling one leaves the other working;
  - `make classy` gives a `lib.*` directory from which `import classy` works with that
    directory alone on `sys.path`;
  - a stale top-level `classy.<soabi>.so` next to the new package: the package wins.
- **Timing:** `benchmarks/classy_glue/timing.py` before and after. The per-z calls must stay
  within noise; also report `get_current_derived_parameters` for a typical Cobaya set of
  about 15 names.
- **Size:** report line counts of `_core.pyx` and of the Python modules.

## Design, phase 2: nanobind prototype

- **Branch:** `classy-nanobind`, on top of phase 1.
- **What gets built:** `src/classy/_core_nb.cpp` implements the same `_core` interface (class
  name, method names, return types, exceptions) with nanobind. A CMake option
  `CLASS_PYTHON_BINDING=cython|nanobind`, default `cython`, builds one or the other as
  `<pkg>/_core`. The plain-Python layer and every test run unchanged on either.
- **Exceptions:** a nanobind translator maps std::invalid_argument to CosmoSevereError and
  std::runtime_error to CosmoComputationError, the same classes the Python package exports.
- **Arrays:** returned as `nb::ndarray` owning the moved std::vector, where that keeps the
  output identical; this removes the extra copy #503 added to `get_background()`.
- **Measure, both bindings, same machine:**
  - lines of binding code, counting `generate_wrapper.py` for Cython, since nanobind needs none;
  - clean build time of the binding;
  - size of the extension;
  - `timing.py`;
  - A/B of the nanobind core against the Cython core (must be identical);
  - all test suites.
- **Deliverable:** `docs/superpowers/specs/YYYY-MM-DD-classy-nanobind-results.md` with the
  numbers and a recommendation, plus the branch pushed (a draft PR at most). Nothing merges
  before Thomas decides.
- **Out of scope:** stable-ABI wheels, free-threaded Python and releasing the GIL. The results
  doc mentions them as consequences, unmeasured.

## Left to the builder

- How the Python package is split into modules.
- The names of the new C++ methods.
- Forwarders not named in the table.
- The details of the dev-build layout.
- Whether the generator's struct harvesting can go.
- The order of commits.

Phase 1 is one PR stacked on #503 (base `classy-glue-to-cpp` until #503 merges).
