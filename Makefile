# Thin convenience shim. The build system is CMakeLists.txt; see README.md.
BUILD_DIR ?= build/cmake

.PHONY: all class class_profiled classy test test-parser test-bisection test-photons clean

# The real parallelism happens inside `cmake --build --parallel`; serialize
# the shim targets so `make -j class class_profiled` cannot run two CMake
# builds in the same build directory at once.
.NOTPARALLEL:

all: class

$(BUILD_DIR)/CMakeCache.txt:
	cmake -S . -B $(BUILD_DIR)

class class_profiled test-parser test-bisection test-photons: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $@ --parallel

# Single cmake invocation so `make -j test` cannot race parallel builds
# of the same build directory.
test: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target test-parser --target test-bisection --target test-photons --parallel
	ctest --test-dir $(BUILD_DIR) --output-on-failure

# Builds + installs classy via pip (scikit-build-core -> CMake), then recreates
# the python/build/lib.* layout MontePython expects.
classy:
	rm -rf python/build && mkdir -p python/build
	pip install .
	python scripts/montepython_layout.py

clean:
	rm -rf $(BUILD_DIR)
