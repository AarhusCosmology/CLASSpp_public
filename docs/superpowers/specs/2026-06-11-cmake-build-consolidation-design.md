# CMake + scikit-build-core build consolidation — Design

**Date:** 2026-06-11
**Status:** Approved by Thomas (design dialogue, this session)

## Goal

One build description for everything. Today the source list is maintained in
three places — `Makefile`, `setup.py`, and `CLASS.xcodeproj` — and they drift.
After this change, `CMakeLists.txt` is the single source of truth; the CLI
executables, the C++ unit tests, the Xcode project, VS Code builds, and the
`classy` Python wheel (via scikit-build-core / cibuildwheel) all derive from it.
Must work on macOS, Linux, and Windows, and in the existing
`build_wheels.yml` cibuildwheel workflow.

## Decisions (made with Thomas)

1. **`cclassy.pxd` is generated at build time into the build directory**, not
   committed. `generate_wrapper.py` gains an output-directory argument with
   proper dependency tracking.
2. **Old build files:** `setup.py`, `MANIFEST.in`, and `CLASS.xcodeproj` are
   deleted. The `Makefile` is replaced by a **thin shim** (`class`, `classy`,
   `clean` targets that call cmake/pip) to preserve muscle memory and CI.
3. **Wheel layout stays byte-for-byte identical**: top-level single-file
   `classy` extension module + `classy/bbn/*.dat` + `classy/hyrec/*.dat`
   data dirs. No risk to MontePython or other consumers; modernizing the
   layout is a possible later PR.
4. **Architecture: one static library compiled once.** All shared sources
   build into static lib `classpp`; every executable and the Python extension
   link it.
5. **Explicit source lists** in the top-level `CMakeLists.txt` (no globbing) —
   stray scratch files in `source/` must not silently join the build, and the
   unit-test mains live inside the library directories.

## Design

### 1. Top-level `CMakeLists.txt`

- `cmake_minimum_required(VERSION 3.21)`, `project(CLASSpp LANGUAGES C CXX)`.
- Static library **`classpp`**: all `source/*.cpp`, `species/*.cpp`,
  `tools/*.cpp` (excluding the three test mains `tools/parser_test.cpp`,
  `tools/bisection_test.cpp`, `species/photons_formula_test.cpp`) plus the
  hyrec C sources (`helium.c`, `history.c`, `hydrogen.c`, `hyrectools.c` —
  matching the Makefile; `hyrec/hyrec.c` is the standalone HyRec driver with
  its own `main()` and stays out).
  Properties:
  - C++17 (`CXX_STANDARD 17`, `CXX_STANDARD_REQUIRED ON`)
  - `POSITION_INDEPENDENT_CODE ON`
  - `CXX_VISIBILITY_PRESET hidden`, `VISIBILITY_INLINES_HIDDEN ON`
    (preserves the classy/classyref symbol-isolation fix documented in
    setup.py; harmless for the executables)
  - compile definitions `HYREC` and `__CLASSDIR__="${PROJECT_SOURCE_DIR}"`
  - include dirs: `include`, `tools`, `source`, `species`, `hyrec`, `.`
  - link `m` and Threads on non-Windows
- Default `CMAKE_BUILD_TYPE Release` when unset (gives `-O3` / `/O2`,
  matching today's flags).
- Executables **`class`** (`main/class.cpp`) and **`class_profiled`**
  (`main/class_profiled.cpp`), linking `classpp`. For non-scikit-build
  configures, `RUNTIME_OUTPUT_DIRECTORY = ${PROJECT_SOURCE_DIR}` so
  `./class explanatory.ini` from the repo root keeps working exactly as now.
- **Unit tests** `test-parser`, `test-bisection`, `test-photons` as
  executables + `add_test` (CTest), built from the existing test mains with
  the same link lines as the current Makefile targets.

### 2. Python extension (gated on scikit-build / option)

Built when `SKBUILD` is set (pip build) or `-DCLASS_BUILD_PYTHON=ON`.

- `find_package(Python COMPONENTS Interpreter Development.Module NumPy)`.
- Custom command 1: run `generate_wrapper.py --output <builddir>` producing
  `cclassy.pxd` in the build dir. `DEPENDS` on the script and the C++ headers
  it parses, so it re-runs exactly when inputs change. The script is modified
  to (a) accept the output directory, (b) restrict its header walk to the five
  real header directories (deterministic, immune to in-tree build dirs), and
  (c) drop the `gcc -E` / `include/tmp` special case — investigation showed it
  is dead code (it compares a full path against the bare string `'common.h'`,
  which never matches), so the script has always parsed the raw header.
- Custom command 2: Cython transpile `classy.pyx` → `classy.cpp` in the build
  dir (`-3 --cplus`, include path = build dir for `cclassy.pxd`).
  `DEPENDS classy.pyx cclassy.pxd`.
- `python_add_library(classy MODULE ... WITH_SOABI)`, links `classpp` and
  NumPy headers, defines `NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION`.
- Installs reproduce today's wheel layout:
  - `classy.<soabi>.so|pyd` at wheel root
  - `bbn/*.dat` → `classy/bbn/`, `hyrec/*.dat` → `classy/hyrec/`
  (`classy.pyx` resolves `class_dir = dirname(__file__)/classy` at runtime —
  unchanged.)
- The stale generated files in the source tree (`cclassy.pxd`, `classy.cpp`,
  `include/tmp`) are deleted and gitignored.

### 3. `pyproject.toml`

- `[build-system] requires = ["scikit-build-core", "cython", "numpy"]`,
  `build-backend = "scikit_build_core.build"`.
- `[tool.scikit-build]`: minimum cmake/ninja versions; sdist configured so all
  build inputs (sources, headers, `generate_wrapper.py`, `classy.pyx`, data
  files) are included.
- Project name `classy-community` and version `23.2.0` unchanged.
- `requires-python` bumped to `>=3.10` (matches the wheel matrix
  cp310–cp314; scikit-build-core does not support 3.7 anyway).

### 4. Makefile shim

Replaces the current Makefile. Targets:

- `class` (default-ish, plus `all`): configure + build via CMake
  (`cmake -S . -B build/cmake && cmake --build build/cmake --target class -j`).
- `class_profiled`, `test-*`: same pattern.
- `classy`: `pip install .`, then recreate the MontePython compatibility
  layout `python/build/lib.<sys.version>.../classy*.so`. Since
  scikit-build-core has no `build/lib.*` directory, the shim copies the
  installed module out of site-packages (`python -c "import classy; ..."`)
  into the same destination the old hack produced. Observable result for
  MontePython is unchanged.
- `clean`: remove the cmake build dir.

### 5. Xcode, VS Code, CI

- **Xcode:** `CLASS.xcodeproj` deleted. Generate locally:
  `cmake -S . -B build-xcode -G Xcode`. With
  `CMAKE_XCODE_GENERATE_SCHEME ON` and, on the `class` target,
  `XCODE_SCHEME_WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}` and
  `XCODE_SCHEME_ARGUMENTS explanatory.ini`, the generated scheme runs CLASS
  from the project dir with `explanatory.ini` — the one Xcode setting worth
  preserving. `build-xcode/` gitignored. Brief README note.
- **VS Code:** committed `.vscode/settings.json` with `cmake.debugConfig`
  (`args: ["explanatory.ini"]`, `cwd: ${workspaceFolder}`) so CMake Tools
  launch matches.
- **CI:** `build.yml` keeps `make class` (now exercising CMake via the shim)
  and `./class explanatory.ini`; add a `ctest` step for the three unit tests.
  `build_wheels.yml` unchanged — cibuildwheel picks up scikit-build-core
  automatically on all three OSes. `test_nightly.yml` uses `pip install .`,
  which keeps working. `test_on_pull_request.yml` builds the reference
  wrapper by renaming the extension in master's `setup.py`; that file
  disappears once this merges, so the module name became a CMake cache
  variable (`CLASS_PYTHON_MODULE_NAME`) and the workflow branches on
  `[ -f setup.py ]` to support both pre- and post-CMake master.

## Error handling

- Missing Cython/numpy in a non-isolated build → scikit-build-core /
  find_package fail with clear messages at configure time.
- `generate_wrapper.py` failure aborts the build (custom command exit code).
- Configuring with `CLASS_BUILD_PYTHON=ON` without a usable Python dev
  environment fails at configure, not link.

## Verification

1. `make class && ./class explanatory.ini` runs to completion.
2. `ctest` green (parser, bisection, photons).
3. `pip install .` succeeds; `import classy`; run a subset of
   `python/test_class.py` scenarios.
4. **Wheel-layout diff:** build a wheel before and after; `unzip -l` listings
   must match (same entries, modulo timestamps).
5. Spot-check CLI spectra output old-vs-new at ~0.1% tolerance (per project
   convention; identical flags should make it effectively identical).
6. `cmake -G Xcode` generates and builds; scheme runs with explanatory.ini.
7. CI: all existing workflows green, including cibuildwheel on the three OS
   runners.

## Out of scope

- Modernizing the wheel/package layout (top-level module → proper package).
- Changing compiler flag policy (e.g. `-march=native`, SIMD options).
- The `CPU`/`CPU.py` helpers and doc build.
