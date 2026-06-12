# CMake + scikit-build-core Build Consolidation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the triple-maintained Makefile / setup.py / Xcode project with a single `CMakeLists.txt` consumed by CLI builds, IDEs, and `pip install` (scikit-build-core), on macOS/Linux/Windows + cibuildwheel.

**Architecture:** One static library `classpp` (all shared C++/C sources, compiled once with PIC + hidden visibility) linked by the `class`/`class_profiled` executables, three CTest unit tests, and the Cython `classy` extension. `cclassy.pxd` and `classy.cpp` are generated into the build dir by custom commands. A thin Makefile shim preserves `make class` / `make classy` muscle memory and the MontePython layout hack.

**Tech Stack:** CMake ≥3.21, scikit-build-core, Cython 3, `FindPython` (`Development.Module` + `NumPy`), cibuildwheel.

**Spec:** `docs/superpowers/specs/2026-06-11-cmake-build-consolidation-design.md`

**Branch:** `cmake-build-system` (already created off master).

**Working assumptions verified during design:**
- `hyrec/hyrec.c` has its own `main()` → excluded from the library (Makefile already excludes it; setup.py's inclusion was inert).
- `source/output_module.cpp` goes into `classpp` even though setup.py omitted it (static lib: unreferenced objects aren't pulled into the extension).
- The `gcc -E` block in `generate_wrapper.py` (lines ~190-193) is dead code: it compares a full path against the literal `'common.h'`, which is never true. Delete it, don't port it.
- Current wheel layout to replicate: top-level `classy.<soabi-tag>.so/.pyd` + `classy/bbn/*.dat` + `classy/hyrec/*.dat` (no `__init__.py` files anywhere).

---

### Task 1: Capture baseline artifacts (old build system)

Everything later is verified against these. Do this BEFORE touching any file.

**Files:** none modified.

- [ ] **Step 1.1: Baseline wheel listing (old setup.py)**

```bash
cd /Users/au192734/Projects/class_claude
python3 -m pip wheel . --no-deps -w /tmp/classy_baseline_wheel
unzip -l /tmp/classy_baseline_wheel/*.whl | sort > /tmp/wheel_baseline.txt
cat /tmp/wheel_baseline.txt
```

Expected: a wheel containing `classy.cpython-*-darwin.so` at root, `classy/bbn/sBBN*.dat` (3 files), `classy/hyrec/*.dat` (4 files), and a `*.dist-info/` dir. Save the listing — it is the layout contract.

- [ ] **Step 1.2: Baseline CLI spectra (old Makefile)**

```bash
make -j class
mkdir -p output /tmp/classy_baseline_cli
./class explanatory.ini
cp output/explanatory*_cl.dat /tmp/classy_baseline_cli/
ls /tmp/classy_baseline_cli
```

Expected: `class` builds and runs to completion; at least `explanatory00_cl.dat` copied. (If `output/explanatory*` already exists, `class` auto-increments the index — copy the newest.)

- [ ] **Step 1.3: Refresh the in-tree generated pxd as reference**

```bash
python3 generate_wrapper.py
cp cclassy.pxd /tmp/cclassy_baseline.pxd
```

Expected: exits 0, `cclassy.pxd` rewritten in the source tree (it is gitignored).

No commit (nothing changed).

---

### Task 2: Make `generate_wrapper.py` build-dir aware

**Files:**
- Modify: `generate_wrapper.py` (lines 1–11, 65–73, 189–193, 408)

- [ ] **Step 2.1: Add an `--output-dir` argument and restrict the header walk**

Replace lines 1–11 (imports + rootdir) with:

```python
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import os
import pathlib

rootdir = pathlib.Path(__file__).parent

_argparser = argparse.ArgumentParser(description='Generate cclassy.pxd from the C++ headers.')
_argparser.add_argument('--output-dir', type=pathlib.Path, default=rootdir,
                        help='Directory to write cclassy.pxd into (default: repo root).')
# parse_known_args: tolerate being run via runpy with foreign argv.
_args, _ = _argparser.parse_known_args()
outdir = _args.output_dir
```

(`import subprocess` is dropped — it was only used by the dead `gcc -E` block.)

Replace the repo-wide walk (lines 65–73, the `h_files`/`ignored_dirs`/`os.walk` block) with a deterministic walk of the real header directories, so a CMake build dir inside the source tree can never contaminate the generated file:

```python
h_files = []
header_dirs = ['include', 'main', 'source', 'species', 'tools']
for d in header_dirs:
    h_files.extend(sorted((rootdir / d).glob('*.h')))
```

- [ ] **Step 2.2: Delete the dead `gcc -E` special case**

In the structs loop (around line 189), delete these four lines entirely:

```python
    # Special treatment of common.h
    if file == 'common.h':
        subprocess.run(['gcc', '-E', rootdir+'/include/common.h','-o', rootdir+'/include/tmp'])
        file = rootdir+'/include/tmp'
```

(`file` holds a full path, so the comparison never matched; the script has always parsed the raw `common.h`.)

- [ ] **Step 2.3: Write to the output dir**

Change the final `open` (line ~408) from:

```python
with open(rootdir / 'cclassy.pxd', 'w', encoding='utf-8') as fid:
```

to:

```python
with open(outdir / 'cclassy.pxd', 'w', encoding='utf-8') as fid:
```

- [ ] **Step 2.4: Verify output is equivalent to baseline**

```bash
mkdir -p /tmp/pxd_new
python3 generate_wrapper.py --output-dir /tmp/pxd_new
diff /tmp/cclassy_baseline.pxd /tmp/pxd_new/cclassy.pxd && echo IDENTICAL
```

Expected: `IDENTICAL`, or differences that are pure re-ordering of `cdef extern` blocks (the walk is now sorted). If blocks are reordered, confirm equivalence with:

```bash
diff <(sort /tmp/cclassy_baseline.pxd) <(sort /tmp/pxd_new/cclassy.pxd) && echo SAME-CONTENT
```

Any *content* difference (missing struct members, missing classes) is a bug in the dir restriction — check which directory the missing header lives in and add it to `header_dirs`.

- [ ] **Step 2.5: Verify default behaviour unchanged (setup.py still works at this commit)**

```bash
python3 generate_wrapper.py && ls -la cclassy.pxd
```

Expected: writes to repo root as before.

- [ ] **Step 2.6: Commit**

```bash
git add generate_wrapper.py
git commit -m "generate_wrapper.py: --output-dir, deterministic header walk, drop dead gcc -E branch"
```

---

### Task 3: Core CMakeLists.txt (library, executables, tests, Xcode scheme)

**Files:**
- Create: `CMakeLists.txt`
- Modify: `.gitignore`

- [ ] **Step 3.1: Write `CMakeLists.txt`**

Complete file (the Python section comes in Task 4 — include the `if(CLASS_BUILD_PYTHON)` stub now):

```cmake
cmake_minimum_required(VERSION 3.21)
project(CLASSpp LANGUAGES C CXX)

if(SKBUILD)
  set(_class_python_default ON)
else()
  set(_class_python_default OFF)
endif()
option(CLASS_BUILD_PYTHON "Build the classy Python extension" ${_class_python_default})

# Single-config generators default to Release (-O3 / /O2, matching the old builds).
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _is_multi_config AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_XCODE_GENERATE_SCHEME ON)

# ------------------------------------------------------------------ sources --
set(CLASS_SOURCE_FILES
  source/background_module.cpp
  source/cosmology.cpp
  source/input_module.cpp
  source/lensing_module.cpp
  source/nonlinear_module.cpp
  source/output_module.cpp
  source/perturbations_module.cpp
  source/primordial_module.cpp
  source/spectra_module.cpp
  source/thermodynamics_module.cpp
  source/transfer_module.cpp
)

set(CLASS_SPECIES_FILES
  species/baryons.cpp
  species/base_species.cpp
  species/cdm.cpp
  species/composite_species.cpp
  species/dark_radiation_species.cpp
  species/dcdm.cpp
  species/dcdm_dr_species.cpp
  species/dncdm_dr_species.cpp
  species/dncdm_species.cpp
  species/fluid.cpp
  species/greybody_ncdm_species.cpp
  species/idm_dr_idr_species.cpp
  species/idm_drmd_idr_drmd_species.cpp
  species/interacting_species.cpp
  species/lambda.cpp
  species/ncdm_base_species.cpp
  species/ncdm_interacting_species.cpp
  species/ncdm_species.cpp
  species/perturb_column_writer.cpp
  species/photons.cpp
  species/scalar_field.cpp
  species/species_collection.cpp
  species/species_input.cpp
  species/ultra_relativistic.cpp
)

set(CLASS_TOOLS_FILES
  tools/arrays.cpp
  tools/common.cpp
  tools/dei_rkck.cpp
  tools/evolver_ndf15.cpp
  tools/evolver_rkck.cpp
  tools/evolver_rkdp45.cpp
  tools/exceptions.cpp
  tools/hyperspherical.cpp
  tools/parser.cpp
  tools/quadrature.cpp
  tools/sparse.cpp
  tools/trigonometric_integrals.cpp
)

# hyrec/hyrec.c is the standalone HyRec driver (has its own main) and stays out.
set(CLASS_HYREC_FILES
  hyrec/helium.c
  hyrec/history.c
  hyrec/hydrogen.c
  hyrec/hyrectools.c
)

# ------------------------------------------------------------------ library --
add_library(classpp STATIC
  ${CLASS_SOURCE_FILES}
  ${CLASS_SPECIES_FILES}
  ${CLASS_TOOLS_FILES}
  ${CLASS_HYREC_FILES}
)
target_compile_features(classpp PUBLIC cxx_std_17)
# Hidden visibility keeps a parallel build (e.g. classyref in the lvl2 nose
# tests) loadable in the same Python process without symbol collisions.
set_target_properties(classpp PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)
target_include_directories(classpp PUBLIC
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/tools
  ${PROJECT_SOURCE_DIR}/source
  ${PROJECT_SOURCE_DIR}/species
  ${PROJECT_SOURCE_DIR}/main
  ${PROJECT_SOURCE_DIR}/hyrec
  ${PROJECT_SOURCE_DIR}
)
target_compile_definitions(classpp PUBLIC
  HYREC
  "__CLASSDIR__=\"${PROJECT_SOURCE_DIR}\""
)
if(NOT WIN32)
  find_package(Threads REQUIRED)
  target_link_libraries(classpp PUBLIC Threads::Threads m)
endif()

# ----------------------------------------------------- executables + tests --
if(NOT SKBUILD)
  add_executable(class main/class.cpp)
  add_executable(class_profiled main/class_profiled.cpp)
  target_link_libraries(class PRIVATE classpp)
  target_link_libraries(class_profiled PRIVATE classpp)
  # Binaries land in the repo root so `./class explanatory.ini` keeps working.
  # ($<1:...> stops multi-config generators appending a per-config subdir.)
  set_target_properties(class class_profiled PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "$<1:${PROJECT_SOURCE_DIR}>"
  )
  set_target_properties(class PROPERTIES
    XCODE_SCHEME_ARGUMENTS "explanatory.ini"
    XCODE_SCHEME_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )

  enable_testing()
  add_executable(test-parser tools/parser_test.cpp)
  add_executable(test-bisection tools/bisection_test.cpp)
  add_executable(test-photons species/photons_formula_test.cpp)
  foreach(_t IN ITEMS test-parser test-bisection test-photons)
    target_link_libraries(${_t} PRIVATE classpp)
    add_test(NAME ${_t} COMMAND ${_t} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  endforeach()
endif()

# ------------------------------------------------------------------- python --
if(CLASS_BUILD_PYTHON)
  # Filled in by Task 4.
endif()
```

- [ ] **Step 3.2: Configure and build**

```bash
cmake -S . -B build/cmake
cmake --build build/cmake --parallel
```

Expected: configures with `CMAKE_BUILD_TYPE=Release`, all targets compile, `./class` and `./class_profiled` appear in the repo root.

- [ ] **Step 3.3: Run the CLI and CTest**

```bash
./class explanatory.ini
ctest --test-dir build/cmake --output-on-failure
```

Expected: CLASS runs to completion; 3/3 tests pass.

- [ ] **Step 3.4: Compare spectra against baseline**

```bash
ls -t output/explanatory*_cl.dat | head -3
# diff the newest 00_cl.dat against the Task-1 baseline:
python3 - <<'EOF'
import glob, os
import numpy as np
new = sorted(glob.glob('output/explanatory*00_cl.dat'), key=os.path.getmtime)[-1]
old = '/tmp/classy_baseline_cli/' + sorted(os.listdir('/tmp/classy_baseline_cli'))[0]
a, b = np.loadtxt(old), np.loadtxt(new)
rel = np.max(np.abs(a[:, 1] - b[:, 1]) / np.abs(a[:, 1]))  # TT column
print('max rel diff TT:', rel)
assert rel < 1e-3, 'spectra differ beyond 0.1% tolerance'
print('OK')
EOF
```

Expected: `OK` (same compiler + `-O3` should make it effectively identical; the project bar is 0.1%, never bit-identical).

- [ ] **Step 3.5: Update `.gitignore`**

Add these lines (the existing `build/` entry already covers `build/cmake`):

```
build-xcode/
```

- [ ] **Step 3.6: Commit**

```bash
git add CMakeLists.txt .gitignore
git commit -m "Add CMakeLists.txt: classpp static lib, CLI executables, CTest unit tests"
```

---

### Task 4: Python extension in CMake

**Files:**
- Modify: `CMakeLists.txt` (the `if(CLASS_BUILD_PYTHON)` block)

- [ ] **Step 4.1: Fill in the Python block**

Replace the `if(CLASS_BUILD_PYTHON)` stub with:

```cmake
if(CLASS_BUILD_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

  # 1) cclassy.pxd, generated from the C++ headers into the build dir.
  file(GLOB _class_wrapper_headers CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/include/*.h
    ${PROJECT_SOURCE_DIR}/main/*.h
    ${PROJECT_SOURCE_DIR}/source/*.h
    ${PROJECT_SOURCE_DIR}/species/*.h
    ${PROJECT_SOURCE_DIR}/tools/*.h
  )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cclassy.pxd
    COMMAND Python::Interpreter ${PROJECT_SOURCE_DIR}/generate_wrapper.py
            --output-dir ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${PROJECT_SOURCE_DIR}/generate_wrapper.py ${_class_wrapper_headers}
    COMMENT "Generating cclassy.pxd"
    VERBATIM
  )

  # 2) Cython transpile (cclassy.pxd is found via -I <build dir>).
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/classy.cpp
    COMMAND Python::Interpreter -m cython -3 --cplus
            -I ${CMAKE_CURRENT_BINARY_DIR}
            -o ${CMAKE_CURRENT_BINARY_DIR}/classy.cpp
            ${PROJECT_SOURCE_DIR}/classy.pyx
    DEPENDS ${PROJECT_SOURCE_DIR}/classy.pyx ${CMAKE_CURRENT_BINARY_DIR}/cclassy.pxd
    COMMENT "Cythonizing classy.pyx"
    VERBATIM
  )

  # 3) The extension module.
  python_add_library(classy MODULE ${CMAKE_CURRENT_BINARY_DIR}/classy.cpp WITH_SOABI)
  target_link_libraries(classy PRIVATE classpp Python::NumPy)
  target_compile_definitions(classy PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
  set_target_properties(classy PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
  )

  # 4) Wheel layout — identical to the old setup.py wheel:
  #    classy.<soabi>.so at the root, data under classy/{bbn,hyrec}.
  if(SKBUILD)
    install(TARGETS classy LIBRARY DESTINATION .)
    file(GLOB _class_bbn_data ${PROJECT_SOURCE_DIR}/bbn/*.dat)
    file(GLOB _class_hyrec_data ${PROJECT_SOURCE_DIR}/hyrec/*.dat)
    install(FILES ${_class_bbn_data} DESTINATION classy/bbn)
    install(FILES ${_class_hyrec_data} DESTINATION classy/hyrec)
  endif()
endif()
```

- [ ] **Step 4.2: Build the extension standalone (no pip yet)**

```bash
cmake -S . -B build/cmake -DCLASS_BUILD_PYTHON=ON
cmake --build build/cmake --parallel --target classy
ls build/cmake/classy.*
```

Expected: `generate_wrapper.py` and Cython run as custom commands; `build/cmake/classy.cpython-3XX-darwin.so` exists. (Cython and numpy must be importable by `python3` — `pip install cython numpy` if missing.)

- [ ] **Step 4.3: Smoke-test import and dependency tracking**

```bash
PYTHONPATH=build/cmake python3 -c "import classy; print(classy.__file__)"
touch include/common.h
cmake --build build/cmake --target classy 2>&1 | grep -E "Generating|Cythonizing"
```

Expected: import succeeds and prints the build-dir path; the touch retriggers BOTH "Generating cclassy.pxd" and "Cythonizing classy.pyx".

- [ ] **Step 4.4: Commit**

```bash
git add CMakeLists.txt
git commit -m "CMake: classy Cython extension with build-dir pxd generation"
```

---

### Task 5: pyproject.toml → scikit-build-core; delete setup.py/MANIFEST.in

**Files:**
- Modify: `pyproject.toml`
- Delete: `setup.py`, `MANIFEST.in`

- [ ] **Step 5.1: Rewrite `pyproject.toml`**

Complete new content:

```toml
[build-system]
requires = ["scikit-build-core>=0.11", "cython>=3.0", "numpy"]
build-backend = "scikit_build_core.build"

[project]
name = "classy-community"
version = "23.2.0"
authors = [
  { name="Thomas Tram", email="thomas.tram@phys.au.dk"},
]
description = "Python interface to the C++ community-edition of the cosmological Boltzmann code CLASS"
readme = "README.md"
requires-python = ">=3.10"
classifiers = [
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    "Operating System :: OS Independent",
]

[project.urls]
"Homepage" = "https://github.com/AarhusCosmology/CLASSpp_public"
"Bug Tracker" = "https://github.com/AarhusCosmology/CLASSpp_public/issues"

[tool.scikit-build]
cmake.version = ">=3.21"
wheel.packages = []
# sdist: scikit-build-core includes the source tree minus .gitignore matches.
# Exclude local clutter that is not gitignored; everything the build needs
# (sources, headers, generate_wrapper.py, classy.pyx, bbn/hyrec data,
# CMakeLists.txt) is included by default.
sdist.exclude = [
  "/doc",
  "/docs",
  "/notebooks",
  "/test",
  "/scripts",
  "/python",
  "*.dSYM",
  "*.ipynb",
]
```

- [ ] **Step 5.2: Delete the old files**

```bash
git rm setup.py MANIFEST.in
```

- [ ] **Step 5.3: pip install + import + greybody test**

```bash
pip install . -v 2>&1 | tail -20
cd python && python3 -m pytest -q test_greybody.py && cd ..
python3 -c "import classy; c = classy.Class(); print('classy OK')"
```

Expected: scikit-build-core drives CMake; install succeeds; greybody tests pass (they compute real cosmology, exercising the data-file lookup `class_dir = <module>/classy`).

- [ ] **Step 5.4: Wheel layout diff against baseline**

```bash
python3 -m pip wheel . --no-deps -w /tmp/classy_new_wheel
unzip -l /tmp/classy_new_wheel/*.whl | sort > /tmp/wheel_new.txt
diff <(grep -oE '[^ ]+$' /tmp/wheel_baseline.txt | grep -v dist-info | grep -v '^$' | sort) \
     <(grep -oE '[^ ]+$' /tmp/wheel_new.txt      | grep -v dist-info | grep -v '^$' | sort)
```

Expected: empty diff of non-dist-info entries (same `classy.*.so` name, same 7 `.dat` paths). If scikit-build-core adds extra entries (e.g. license files), confirm they're additive-only and harmless; any MISSING entry must be fixed in the `install()` rules.

- [ ] **Step 5.5: sdist sanity**

```bash
pipx run build --sdist
tar -tzf dist/classy_community-23.2.0.tar.gz | grep -cE "source/|species/|tools/|hyrec/.*\.(c|h|dat)|bbn/.*\.dat|classy\.pyx|generate_wrapper\.py|CMakeLists.txt"
pip install dist/classy_community-23.2.0.tar.gz --force-reinstall -q && python3 -c "import classy; print('sdist install OK')"
```

Expected: count > 100 (all build inputs present) and the install-from-sdist succeeds.

- [ ] **Step 5.6: Commit**

```bash
git add pyproject.toml
git commit -m "Switch Python build to scikit-build-core; drop setup.py and MANIFEST.in"
```

---

### Task 6: Makefile shim + MontePython layout script

**Files:**
- Rewrite: `Makefile`
- Create: `scripts/montepython_layout.py`

- [ ] **Step 6.1: Write `scripts/montepython_layout.py`**

```python
"""Recreate the python/build/lib.* layout that MontePython >= 3.6 expects.

setuptools used to leave the built extension in build/lib.<platform>-<tag>/,
which `make classy` renamed to lib.<sys.version>.<platform>-<tag>/ (see
https://github.com/brinckmann/montepython_public/issues/371). scikit-build-core
has no such directory, so this script reproduces the same end result from the
module installed in site-packages.
"""
import os
import shutil
import sys
import sysconfig

repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Never import a stray in-repo classy/ directory (dev symlink hack).
sys.path = [p for p in sys.path
            if os.path.abspath(p if p else os.getcwd()) != repo_root]
import classy  # noqa: E402

libdir = os.path.join(
    repo_root, 'python', 'build',
    'lib.{}.{}-{}'.format(sys.version, sysconfig.get_platform(),
                          sys.implementation.cache_tag))
os.makedirs(libdir, exist_ok=True)
shutil.copy2(classy.__file__, libdir)
data_dir = os.path.join(os.path.dirname(classy.__file__), 'classy')
shutil.copytree(data_dir, os.path.join(libdir, 'classy'), dirs_exist_ok=True)
print('MontePython layout created at', libdir)
```

- [ ] **Step 6.2: Replace the `Makefile` with the shim**

Complete new content:

```make
# Thin convenience shim. The build system is CMakeLists.txt; see README.md.
BUILD_DIR ?= build/cmake

.PHONY: all class class_profiled classy test test-parser test-bisection test-photons clean

all: class

$(BUILD_DIR)/CMakeCache.txt:
	cmake -S . -B $(BUILD_DIR)

class class_profiled test-parser test-bisection test-photons: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $@ --parallel

test: test-parser test-bisection test-photons
	ctest --test-dir $(BUILD_DIR) --output-on-failure

# Builds + installs classy via pip (scikit-build-core -> CMake), then recreates
# the python/build/lib.* layout MontePython expects.
classy:
	rm -rf python/build && mkdir -p python/build
	pip install .
	python scripts/montepython_layout.py

clean:
	rm -rf $(BUILD_DIR)
```

(Note: recipe lines must be TAB-indented.)

- [ ] **Step 6.3: Verify the shim**

```bash
make clean && make class && ./class explanatory.ini && make test
make classy && ls python/build/
```

Expected: `class` builds via CMake and runs; ctest 3/3; `python/build/lib.3.14.…` dir exists containing `classy.*.so` and `classy/{bbn,hyrec}/*.dat`.

- [ ] **Step 6.4: Commit**

```bash
git add Makefile scripts/montepython_layout.py
git commit -m "Replace Makefile with thin CMake shim; MontePython layout via script"
```

---

### Task 7: Delete Xcode project; add VS Code config; README

**Files:**
- Delete: `CLASS.xcodeproj/`
- Create: `.vscode/settings.json`
- Modify: `README.md` (build instructions section)

- [ ] **Step 7.1: Delete the Xcode project**

```bash
git rm -r CLASS.xcodeproj
rm -rf CLASS.xcodeproj   # clears leftover gitignored user files
```

- [ ] **Step 7.2: Verify the CMake-generated replacement (macOS only)**

```bash
cmake -S . -B build-xcode -G Xcode
xcodebuild -list -project build-xcode/CLASSpp.xcodeproj
grep -l "explanatory.ini" build-xcode/CLASSpp.xcodeproj/xcshareddata/xcschemes/class.xcscheme
grep -l "customWorkingDirectory" build-xcode/CLASSpp.xcodeproj/xcshareddata/xcschemes/class.xcscheme
```

Expected: schemes listed include `class`; the `class` scheme file contains the `explanatory.ini` argument and a custom working directory set to the repo root.

- [ ] **Step 7.3: Add `.vscode/settings.json`**

```json
{
  "cmake.buildDirectory": "${workspaceFolder}/build/cmake",
  "cmake.debugConfig": {
    "args": ["explanatory.ini"],
    "cwd": "${workspaceFolder}"
  }
}
```

- [ ] **Step 7.4: Update README build section**

Find the compilation/installation section in `README.md` and replace make/setup.py instructions with (adapt surrounding prose as needed):

```markdown
## Building

CLASS uses CMake. Quick start:

    make class          # CLI binary ./class (thin shim around CMake)
    make classy         # Python wrapper: pip install . + MontePython layout
    make test           # C++ unit tests via CTest

or directly:

    cmake -S . -B build/cmake && cmake --build build/cmake --parallel
    pip install .       # Python wrapper (scikit-build-core)

For faster repeated Python installs with build dependencies preinstalled
(cython, numpy, scikit-build-core): `pip install --no-build-isolation .`

IDE projects are generated, not committed:

    cmake -S . -B build-xcode -G Xcode    # Xcode (run schemes are preconfigured
                                          # to run from the repo root with
                                          # explanatory.ini)

VS Code: install the CMake Tools extension; .vscode/settings.json already
points it at build/cmake and sets the launch args/cwd.
```

- [ ] **Step 7.5: Commit**

```bash
git add -A .vscode README.md
git commit -m "Drop committed Xcode project (generate via cmake -G Xcode); VS Code + README build docs"
```

---

### Task 8: CI workflow updates

**Files:**
- Modify: `.github/workflows/build.yml`

(`build_wheels.yml` and `test_nightly.yml` need no changes: cibuildwheel and `pip install .` pick up scikit-build-core automatically. `test_on_pull_request.yml` DOES need one — caught in final review: its `make reference` step seds `Extension('classy'` → `Extension('classyref'` in the master checkout's `setup.py`, which disappears once this merges. Fixed via a `CLASS_PYTHON_MODULE_NAME` CMake cache variable plus an `if [ -f setup.py ]` branch in the workflow, so the step works against both pre- and post-CMake master.)

- [ ] **Step 8.1: Update `build.yml` build job**

In the `build` job, replace the `make`/`run class` steps with:

```yaml
    - name: cmake version
      run: cmake --version
    - name: make
      run: cd main_class && make class
    - name: run class
      run: cd main_class && ./class explanatory.ini
    - name: unit tests
      run: cd main_class && make test
```

(`make -j class` → `make class`: the shim always builds with `--parallel`. The `cmake --version` step gives an actionable failure if the self-hosted runner lacks CMake ≥3.21 — install it on the runner if this fails.)

- [ ] **Step 8.2: Commit**

```bash
git add .github/workflows/build.yml
git commit -m "CI: build job uses CMake shim + ctest"
```

---

### Task 9: Final verification + PR

- [ ] **Step 9.1: Clean-slate end-to-end**

```bash
make clean && rm -rf build-xcode
git status --short            # confirm no unexpected tracked changes
make class && ./class explanatory.ini
make test
pip install . --force-reinstall
cd python && python3 -m pytest -q test_greybody.py && cd ..
```

Expected: all green from scratch.

- [ ] **Step 9.2: Scenario smoke subset**

```bash
cd python && TEST_LEVEL=1 COMPARE_OUTPUT_REF=0 python3 -m pytest -v -m test_scenario test_class.py -k "0001 or 0002 or 0003" && cd ..
```

Expected: selected scenarios pass (full level-1 matrix runs in CI).

- [ ] **Step 9.3: Push and open PR**

```bash
git push -u origin cmake-build-system
gh pr create --title "Consolidate build system: CMake + scikit-build-core" --body "..."
```

PR body should call out: single source of truth in CMakeLists.txt; wheel layout unchanged (verified by listing diff); MontePython layout preserved via scripts/montepython_layout.py; Xcode project now generated (`cmake -G Xcode`) with run scheme = repo root + explanatory.ini; `requires-python` bumped to >=3.10; cibuildwheel workflow unchanged. Mention that **cibuildwheel on all three OSes is the one thing only CI can prove** — watch `build_wheels.yml` on the PR (or trigger it manually if it only runs on master pushes).

- [ ] **Step 9.4: Watch CI**

All workflows green, especially the Windows wheel (MSVC + generate_wrapper + Cython path).

---

## Self-review notes

- Spec coverage: spec §1→Task 3, §2→Tasks 2+4, §3→Task 5, §4→Task 6, §5→Tasks 7+8, verification→Tasks 1/3.4/5.4/5.5/9. `requires-python` bump and stale-generated-file cleanup are in Task 5 / gitignore already covers `cclassy.pxd`, `classy.cpp`, `include/tmp`.
- `build_wheels.yml` only triggers on pushes to master — Step 9.3 explicitly flags manual verification need.
- Naming consistent: `classpp`, `classy`, `CLASS_BUILD_PYTHON`, `build/cmake` used throughout.
