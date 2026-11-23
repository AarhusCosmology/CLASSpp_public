# Thin convenience shim. The build system is CMakeLists.txt; see README.md.
BUILD_DIR ?= build/cmake
XCODE_DIR ?= build/xcode
# The Python the classy targets build and install for, e.g. `PYTHON=python3.12 make xcode`.
PYTHON ?= python

.PHONY: all class class_profiled classy classy-pip-dev xcode test clean

# The real parallelism happens inside `cmake --build --parallel`; serialize
# the shim targets so `make -j class class_profiled` cannot run two CMake
# builds in the same build directory at once.
.NOTPARALLEL:

all: class

$(BUILD_DIR)/CMakeCache.txt:
	cmake -S . -B $(BUILD_DIR)

class class_profiled: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $@ --parallel

# One unit test, e.g. `make test-parser`: any target unittests/CMakeLists.txt declares.
test-%: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $@ --parallel

# Every unit test unittests/CMakeLists.txt declares (its unit-tests target), in a
# single cmake invocation so `make -j test` cannot race parallel builds of the
# same build directory.
test: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target unit-tests --parallel
	ctest --test-dir $(BUILD_DIR) --output-on-failure

# Builds + installs classy via pip (scikit-build-core -> CMake), then recreates
# the python/build/lib.* layout MontePython expects.
classy:
	rm -rf python/build && mkdir -p python/build
	$(PYTHON) -m pip install .
	$(PYTHON) scripts/montepython_layout.py

# Faster local wrapper rebuild: reuse scikit-build-core's CMake tree and build
# with 12 parallel jobs. This assumes build dependencies are already installed.
# Direct CMake can be faster still:
#   cmake -S . -B build/dev -DCLASS_BUILD_PYTHON=ON
#   cmake --build build/dev --target classy_core --parallel
# which leaves a runnable package in build/dev/classy/ (PYTHONPATH=build/dev),
# but neither installs it nor recreates the python/build/lib.* layout.
classy-pip-dev:
	rm -rf python/build && mkdir -p python/build
	$(PYTHON) -m pip install . --no-build-isolation -Cbuild-dir=build/skbuild -Cbuild.tool-args=-j12
	$(PYTHON) scripts/montepython_layout.py

# Xcode project with the classy core (classy_core: src/classy/_core.cpp), opened in Xcode.
# PYTHON gets nanobind (the bound in pyproject.toml) if it lacks it, and the core is built
# against it. Building classy_core in Xcode leaves an importable package in $(XCODE_DIR)/classy/
# (PYTHONPATH=$(XCODE_DIR)). A build tree of its own: generated in the repo root, that package
# would shadow the installed classy for every Python started there.
xcode:
	$(PYTHON) -m pip install -q "nanobind>=3.1"
	cmake -G Xcode -S . -B $(XCODE_DIR) -DCLASS_BUILD_PYTHON=ON \
	  -DPython_EXECUTABLE="$$($(PYTHON) -c 'import sys; print(sys.executable)')"
	open $(XCODE_DIR)/CLASSpp.xcodeproj

clean:
	rm -rf $(BUILD_DIR)
