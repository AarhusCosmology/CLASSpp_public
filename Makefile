# Thin convenience shim. The build system is CMakeLists.txt; see README.md.
BUILD_DIR ?= build/cmake

TEST_TARGETS := \
	test-parser \
	test-bisection \
	test-quadrature \
	test-photons \
	test-species-types \
	test-scf-potential \
	test-cdm-coupled \
	test-scf-beta \
	test-composite-classification \
	test-type3-gauge-guard \
	test-composite-layout \
	test-type3-coupling \
	test-axion-scf-factory \
	test-axion-ede-fluid \
	test-dcdm-wdm \
	test-axion-ncdm \
	test-dncdm-switch-copy

.PHONY: all class class_profiled classy classy-pip-dev test clean $(TEST_TARGETS)

# The real parallelism happens inside `cmake --build --parallel`; serialize
# the shim targets so `make -j class class_profiled` cannot run two CMake
# builds in the same build directory at once.
.NOTPARALLEL:

all: class

$(BUILD_DIR)/CMakeCache.txt:
	cmake -S . -B $(BUILD_DIR)

class class_profiled $(TEST_TARGETS): $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $@ --parallel

# Single cmake invocation so `make -j test` cannot race parallel builds
# of the same build directory.
test: $(BUILD_DIR)/CMakeCache.txt
	cmake --build $(BUILD_DIR) --target $(TEST_TARGETS) --parallel
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
