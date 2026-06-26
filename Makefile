# Thin convenience shim. The build system is CMakeLists.txt; see README.md.
BUILD_DIR ?= build/cmake

.PHONY: all class class_profiled classy classy-pip-dev test test-parser test-bisection test-photons test-species-types clean

# The real parallelism happens inside `cmake --build --parallel`; serialize
# the shim targets so `make -j class class_profiled` cannot run two CMake
# builds in the same build directory at once.
.NOTPARALLEL:

all: class

$(BUILD_DIR)/CMakeCache.txt:
	cmake -S . -B $(BUILD_DIR)

class class_profiled test-parser test-bisection test-photons test-species-types: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $@ --parallel

# Single cmake invocation so `make -j test` cannot race parallel builds
# of the same build directory.
test: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target test-parser test-bisection test-photons test-species-types --parallel
	ctest --test-dir $(BUILD_DIR) --output-on-failure

# Builds + installs classy via pip (scikit-build-core -> CMake), then recreates
# the python/build/lib.* layout MontePython expects.
classy:
	rm -rf python/build && mkdir -p python/build
	pip install .
	python scripts/montepython_layout.py

# Faster local wrapper rebuild: reuse scikit-build-core's CMake tree and build
# with 12 parallel jobs. This assumes build dependencies are already installed.
# Direct CMake can be faster still:
#   cmake -S . -B build/dev -DCLASS_BUILD_PYTHON=ON
#   cmake --build build/dev --target classy --parallel
# but that only builds the extension in build/dev; it does not install the
# package or recreate the python/build/lib.* layout.
classy-pip-dev:
	rm -rf python/build && mkdir -p python/build
	pip install . --no-build-isolation -Cbuild-dir=build/skbuild -Cbuild.tool-args=-j12
	python scripts/montepython_layout.py

clean:
	rm -rf $(BUILD_DIR)
