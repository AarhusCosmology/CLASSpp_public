cmake_minimum_required(VERSION 3.21)
project(CLASSpp LANGUAGES C CXX)

# Emit build/<dir>/compile_commands.json so clangd (and any other libclang
# tooling) sees the real include paths and flags instead of guessing. .clangd
# points at BUILD_DIR from the Makefile. Honoured by the Makefile and Ninja
# generators only -- Xcode/Visual Studio ignore it, so those users need a
# throwaway Makefile-generator build dir to get a database.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Version: pyproject.toml's [project] version is the single source of truth.
# scikit-build-core passes it as SKBUILD_PROJECT_VERSION for wheel builds; for a
# standalone CMake build we read the `version = "x.y.z"` line out of pyproject.
# Either way we generate version.h so the C++ _VERSION_ macro is derived from
# that one value instead of being hand-edited.
if(DEFINED SKBUILD_PROJECT_VERSION)
  set(CLASS_VERSION "${SKBUILD_PROJECT_VERSION}")
else()
  file(STRINGS "${PROJECT_SOURCE_DIR}/pyproject.toml" _class_version_line
       REGEX "^version *=")
  string(REGEX REPLACE "^version *= *\"([^\"]+)\".*" "\\1"
         CLASS_VERSION "${_class_version_line}")
endif()
configure_file(
  "${PROJECT_SOURCE_DIR}/include/version.h.in"
  "${PROJECT_BINARY_DIR}/generated/version.h"
  @ONLY
)

if(SKBUILD)
  set(_class_python_default ON)
else()
  set(_class_python_default OFF)
endif()
option(CLASS_BUILD_PYTHON "Build the classy Python extension" ${_class_python_default})
option(CLASS_FAST_MATH "Enable relaxed floating-point optimizations" ON)
set(CLASS_ARCH_FLAGS "" CACHE STRING
  "Additional architecture flags, e.g. -march=x86-64-v3 or /arch:AVX2")

# The lvl2 nose tests install a second wrapper built from master under the
# name `classyref`; the package name must therefore be configurable.
set(CLASS_PYTHON_MODULE_NAME classy CACHE STRING
  "Import name of the Python package")

# Single-config generators default to Release (-O3 / /O2, matching the old builds).
get_property(_is_multi_config GLOBAL PROPERTY GENERATOR_IS_MULTI_CONFIG)
if(NOT _is_multi_config AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Match the historical Makefile flags (-O3 -g, no NDEBUG): assert() stays
# active — the unit tests check via assert, and the valgrind nightlies want
# line info in the Release binary.
foreach(_lang C CXX)
  if(MSVC)
    string(REPLACE "/DNDEBUG" "" CMAKE_${_lang}_FLAGS_RELEASE "${CMAKE_${_lang}_FLAGS_RELEASE}")
  else()
    string(REPLACE "-DNDEBUG" "" CMAKE_${_lang}_FLAGS_RELEASE "${CMAKE_${_lang}_FLAGS_RELEASE}")
    string(APPEND CMAKE_${_lang}_FLAGS_RELEASE " -g")
  endif()
endforeach()

set(CMAKE_XCODE_GENERATE_SCHEME ON)

# ------------------------------------------------------------------ sources --
set(CLASS_SOURCE_FILES
  source/background_module.cpp
  source/cosmology.cpp
  source/energy_deposition.cpp
  source/input_module.cpp
  source/lensing_module.cpp
  source/nonlinear_module.cpp
  source/output_module.cpp
  source/perturbations_module.cpp
  source/hyrec_model.cpp
  source/recfast_model.cpp
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
  species/dcdm_wdm_species.cpp
  species/dncdm_dr_species.cpp
  species/dncdm_sector.cpp
  species/dncdm_inv_species.cpp
  species/dncdm_proxy_species.cpp
  species/dncdm_species.cpp
  species/dr_psd_species.cpp
  species/fluid.cpp
  species/axion_ede_fluid.cpp
  species/ppf_fluid.cpp
  species/greybody_ncdm_species.cpp
  species/axion_ncdm_species.cpp
  species/idm_dr_idr_species.cpp
  species/idm_drmd_idr_drmd_species.cpp
  species/interacting_species.cpp
  species/lambda.cpp
  species/majoron_species.cpp
  species/ncdm_base_species.cpp
  species/ncdm_family.cpp
  species/ncdm_interacting_species.cpp
  species/ncdm_species.cpp
  species/nede_species.cpp
  species/photons.cpp
  species/scalar_field.cpp
  species/scalar_tensor_species.cpp
  species/scalar_field_potential.cpp
  species/species_collection.cpp
  species/species_input.cpp
  species/type3_species.cpp
  species/ultra_relativistic.cpp
  species/wdm_decay_product.cpp
)

set(CLASS_TOOLS_FILES
  tools/arrays.cpp
  tools/asymptotic_extrapolation.cpp
  tools/bbn_plasma.cpp
  tools/bbn_rates.cpp
  tools/bbn_solver.cpp
  tools/bbn_weak.cpp
  tools/common.cpp
  tools/output_table.cpp
  tools/decay_harness.cpp
  tools/decay_transition_kernel.cpp
  tools/dei_rkck.cpp
  tools/acoustic_spline.cpp
  tools/evolver_erk.cpp
  tools/evolver_etd.cpp
  tools/evolver_ndf15.cpp
  tools/evolver_options.cpp
  tools/nonfinite.cpp
  tools/hyperspherical.cpp
  tools/parser.cpp
  tools/quadrature.cpp
  tools/reduced_collision_operator.cpp
  tools/sparse.cpp
  tools/trigonometric_integrals.cpp
  tools/wigner_d.cpp
)

# hyrec/hyrec.c is the standalone HyRec driver (has its own main) and stays out.
set(CLASS_HYREC_FILES
  hyrec/energy_injection.c
  hyrec/helium.c
  hyrec/history.c
  hyrec/hydrogen.c
  hyrec/hyrectools.c
)

# Headers are listed only so IDE projects (Xcode/VS) display them; CMake marks
# .h files HEADER_FILE_ONLY and never compiles them. CONFIGURE_DEPENDS re-globs
# on build, so adding/removing a header refreshes the project automatically.
file(GLOB CLASS_HEADER_FILES CONFIGURE_DEPENDS
  ${PROJECT_SOURCE_DIR}/include/*.h
  ${PROJECT_SOURCE_DIR}/main/*.h
  ${PROJECT_SOURCE_DIR}/source/*.h
  ${PROJECT_SOURCE_DIR}/species/*.h
  ${PROJECT_SOURCE_DIR}/tools/*.h
  ${PROJECT_SOURCE_DIR}/hyrec/*.h
)

# ------------------------------------------------------------------ library --
add_library(classpp STATIC
  ${CLASS_SOURCE_FILES}
  ${CLASS_SPECIES_FILES}
  ${CLASS_TOOLS_FILES}
  ${CLASS_HYREC_FILES}
  ${CLASS_HEADER_FILES}
)

# Mirror the on-disk folder layout into IDE groups instead of one flat list.
source_group(TREE ${PROJECT_SOURCE_DIR} FILES
  ${CLASS_SOURCE_FILES}
  ${CLASS_SPECIES_FILES}
  ${CLASS_TOOLS_FILES}
  ${CLASS_HYREC_FILES}
  ${CLASS_HEADER_FILES}
)
target_compile_features(classpp PUBLIC cxx_std_17)
if(CLASS_FAST_MATH)
  if(MSVC)
    target_compile_options(classpp PRIVATE /fp:fast)
  else()
    target_compile_options(classpp PRIVATE -ffast-math)
  endif()
endif()
# tools/nonfinite.cpp holds the ONLY non-finite test in the build. Under
# -ffast-math (-ffinite-math-only) every such test folds to false -- measured, and
# true even for the memcpy bit-pattern form -- so this one file is compiled without
# it. Removing this line silently disables NaN/Inf detection everywhere.
if(CLASS_FAST_MATH)
  if(MSVC)
    set_source_files_properties(tools/nonfinite.cpp PROPERTIES COMPILE_OPTIONS "/fp:precise")
  else()
    set_source_files_properties(tools/nonfinite.cpp PROPERTIES COMPILE_OPTIONS "-fno-fast-math")
  endif()
endif()
if(CLASS_ARCH_FLAGS)
  separate_arguments(_class_arch_flags NATIVE_COMMAND "${CLASS_ARCH_FLAGS}")
  target_compile_options(classpp PRIVATE ${_class_arch_flags})
endif()
# Hidden visibility keeps a parallel build (e.g. classyref in the lvl2 nose
# tests) loadable in the same Python process without symbol collisions.
set_target_properties(classpp PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON
)
target_include_directories(classpp
  PUBLIC
    ${PROJECT_BINARY_DIR}/generated
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/tools
    ${PROJECT_SOURCE_DIR}/source
    ${PROJECT_SOURCE_DIR}/species
    ${PROJECT_SOURCE_DIR}/main
    ${PROJECT_SOURCE_DIR}
  PRIVATE
    ${PROJECT_SOURCE_DIR}/hyrec
)
target_compile_definitions(classpp PUBLIC
  HYREC
  "__CLASSDIR__=\"${PROJECT_SOURCE_DIR}\""
)
if(NOT WIN32)
  find_package(Threads REQUIRED)
  target_link_libraries(classpp PRIVATE Threads::Threads m)
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
  set_target_properties(class class_profiled PROPERTIES
    XCODE_SCHEME_ARGUMENTS "explanatory.ini"
    XCODE_SCHEME_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
  )
  # Standalone decay-kernel driver (tools/decay_harness.h, issue #389). Stays in the
  # build tree: it is a validation tool, not part of the `./class` workflow.
  add_executable(decay_harness main/decay_harness.cpp)
  target_link_libraries(decay_harness PRIVATE classpp)
  # Hannestad & Madsen's three-level cascade on the same engine (issue #389, step 3).
  add_executable(decay_cascade main/decay_cascade.cpp)
  target_link_libraries(decay_cascade PRIVATE classpp)

  enable_testing()
  add_subdirectory(unittests)
endif()

# ------------------------------------------------------------------- python --
if(CLASS_BUILD_PYTHON)
  find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)

  # 1) The compiled core, src/classy/_core.cpp, in nanobind
  #    (docs/superpowers/specs/2026-11-21-classy-nanobind-results.md). The module is always
  #    _core (PyInit__core); a renamed build (CLASS_PYTHON_MODULE_NAME) only installs the
  #    package under another name, and the module takes its full name from the import.
  #    nanobind's CMake package lives inside its Python package.
  set(_class_pkg_src ${PROJECT_SOURCE_DIR}/src/classy)
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT
    RESULT_VARIABLE _class_nanobind_found)
  if(NOT _class_nanobind_found EQUAL 0)
    message(FATAL_ERROR "classy needs nanobind in ${Python_EXECUTABLE} (pip install nanobind)")
  endif()
  find_package(nanobind CONFIG REQUIRED)

  # 2) The extension module, built straight into a package directory in the build tree, so
  #    that PYTHONPATH=<build dir> imports a runnable package (see 3). NB_STATIC: the nanobind
  #    runtime goes into the module. NB_DOMAIN: each package (classy, classyref, ...) has its own
  #    nanobind internals. Without it, two of them in one process share one type registry, so the
  #    second cannot register its PyCosmology, and classy's exception translator would apply to
  #    every other nanobind module in the process.
  set(_class_pkg_dir ${CMAKE_CURRENT_BINARY_DIR}/${CLASS_PYTHON_MODULE_NAME})
  nanobind_add_module(classy_core NB_STATIC NB_DOMAIN ${CLASS_PYTHON_MODULE_NAME}
                      ${_class_pkg_src}/_core.cpp)
  target_link_libraries(classy_core PRIVATE classpp)
  set_target_properties(classy_core PROPERTIES
    OUTPUT_NAME _core
    # $<1:...>: multi-config generators (Xcode, Visual Studio) would add a per-config
    # subdirectory to a plain path, and the package would not be importable from the build tree.
    LIBRARY_OUTPUT_DIRECTORY $<1:${_class_pkg_dir}>
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
  )

  # 3) The rest of the package next to it in the build tree: the plain-Python modules, copied
  #    on every build so that edits are picked up, and the data, linked (copied where links
  #    are not possible).
  file(GLOB _class_pkg_py CONFIGURE_DEPENDS ${_class_pkg_src}/*.py)
  set(_class_pkg_py_out)
  foreach(_py ${_class_pkg_py})
    get_filename_component(_py_name ${_py} NAME)
    add_custom_command(
      OUTPUT ${_class_pkg_dir}/${_py_name}
      COMMAND ${CMAKE_COMMAND} -E copy ${_py} ${_class_pkg_dir}/${_py_name}
      DEPENDS ${_py}
      VERBATIM
    )
    list(APPEND _class_pkg_py_out ${_class_pkg_dir}/${_py_name})
  endforeach()
  add_custom_target(classy_python_files ALL DEPENDS ${_class_pkg_py_out})
  add_dependencies(classy_core classy_python_files)
  file(GLOB _class_bbn_data CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/bbn/*.dat)
  file(GLOB _class_hyrec_data CONFIGURE_DEPENDS ${PROJECT_SOURCE_DIR}/hyrec/*.dat)
  foreach(_subdir bbn hyrec)
    file(MAKE_DIRECTORY ${_class_pkg_dir}/data/${_subdir})
  endforeach()
  foreach(_dat ${_class_bbn_data})
    get_filename_component(_dat_name ${_dat} NAME)
    file(CREATE_LINK ${_dat} ${_class_pkg_dir}/data/bbn/${_dat_name} SYMBOLIC COPY_ON_ERROR)
  endforeach()
  foreach(_dat ${_class_hyrec_data})
    get_filename_component(_dat_name ${_dat} NAME)
    file(CREATE_LINK ${_dat} ${_class_pkg_dir}/data/hyrec/${_dat_name} SYMBOLIC COPY_ON_ERROR)
  endforeach()

  # 4) Wheel layout: the package <name>/ with __init__.py and the other modules, _core.<soabi>.so
  #    and data/{bbn,hyrec}. Each renamed build (classyref) carries its own data, so two of
  #    them never share a file.
  if(SKBUILD)
    install(TARGETS classy_core LIBRARY DESTINATION ${CLASS_PYTHON_MODULE_NAME})
    install(FILES ${_class_pkg_py} DESTINATION ${CLASS_PYTHON_MODULE_NAME})
    install(FILES ${_class_bbn_data} DESTINATION ${CLASS_PYTHON_MODULE_NAME}/data/bbn)
    install(FILES ${_class_hyrec_data} DESTINATION ${CLASS_PYTHON_MODULE_NAME}/data/hyrec)
  endif()
endif()
