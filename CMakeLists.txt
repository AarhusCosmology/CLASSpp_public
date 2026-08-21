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
# name `classyref`; the import name must therefore be configurable.
set(CLASS_PYTHON_MODULE_NAME classy CACHE STRING
  "Import name of the Python extension module")

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
  species/ncdm_base_species.cpp
  species/ncdm_family.cpp
  species/ncdm_interacting_species.cpp
  species/ncdm_species.cpp
  species/perturb_column_writer.cpp
  species/photons.cpp
  species/scalar_field.cpp
  species/scalar_field_potential.cpp
  species/species_collection.cpp
  species/species_input.cpp
  species/type3_species.cpp
  species/ultra_relativistic.cpp
  species/wdm_decay_product.cpp
)

set(CLASS_TOOLS_FILES
  tools/arrays.cpp
  tools/common.cpp
  tools/decay_transition_kernel.cpp
  tools/dei_rkck.cpp
  tools/evolver_etd.cpp
  tools/evolver_ndf15.cpp
  tools/evolver_rkck.cpp
  tools/evolver_rkdp45.cpp
  tools/nonfinite.cpp
  tools/exceptions.cpp
  tools/hyperspherical.cpp
  tools/parser.cpp
  tools/quadrature.cpp
  tools/reduced_collision_operator.cpp
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

  enable_testing()
  add_executable(test-parser tools/parser_test.cpp)
  add_executable(test-precision-evolver tools/precision_evolver_test.cpp)
  add_executable(test-errors tools/errors_test.cpp)
  add_executable(test-deposition source/energy_deposition_test.cpp)
  add_executable(test-bisection tools/bisection_test.cpp)
  add_executable(test-quadrature tools/quadrature_test.cpp)
  add_executable(test-decay-kernel tools/decay_kernel_test.cpp)
  add_executable(test-reduced-operator tools/reduced_operator_test.cpp)
  add_executable(test-evolver-etd tools/evolver_etd_test.cpp)
  add_executable(test-evolver-nonfinite tools/evolver_nonfinite_test.cpp)
  add_executable(test-photons species/photons_formula_test.cpp)
  add_executable(test-species-types species/species_type_name_test.cpp)
  add_executable(test-scf-potential species/scalar_field_potential_test.cpp)
  add_executable(test-cdm-coupled species/cdm_coupled_test.cpp)
  add_executable(test-scf-beta species/scalar_field_beta_test.cpp)
  add_executable(test-composite-classification species/composite_classification_test.cpp)
  add_executable(test-type3-gauge-guard species/type3_gauge_guard_test.cpp)
  add_executable(test-composite-layout species/composite_layout_test.cpp)
  add_executable(test-type3-coupling species/type3_coupling_test.cpp)
  add_executable(test-axion-scf-factory species/axion_scf_factory_test.cpp)
  add_executable(test-axion-ede-fluid species/axion_ede_fluid_test.cpp)
  add_executable(test-dcdm-wdm species/dcdm_wdm_test.cpp)
  add_executable(test-axion-ncdm species/axion_ncdm_test.cpp)
  add_executable(test-dncdm-switch-copy species/dncdm_switch_copy_test.cpp)
  add_executable(test-dncdm-inv species/dncdm_inv_test.cpp)
  add_executable(test-dncdm-proxy species/dncdm_proxy_test.cpp)
  add_executable(test-dr-psd species/dr_psd_test.cpp)
  add_executable(test-ncdm-family species/ncdm_family_test.cpp)
  add_executable(test-species-budget species/species_budget_test.cpp)
  foreach(_t IN ITEMS test-parser test-precision-evolver test-errors test-deposition test-bisection test-quadrature test-decay-kernel test-reduced-operator test-evolver-etd test-evolver-nonfinite test-photons test-species-types test-scf-potential test-cdm-coupled test-scf-beta test-composite-classification test-type3-gauge-guard test-composite-layout test-type3-coupling test-axion-scf-factory test-axion-ede-fluid test-dcdm-wdm test-axion-ncdm test-dncdm-switch-copy test-dncdm-inv test-dncdm-proxy test-dr-psd test-ncdm-family test-species-budget)
    target_link_libraries(${_t} PRIVATE classpp)
    # These tests assert on values; keep their assert()s live even under
    # RelWithDebInfo/Release (which add -DNDEBUG). -UNDEBUG re-enables asserts for
    # the test executables only; the classpp library stays optimized with NDEBUG.
    target_compile_options(${_t} PRIVATE -UNDEBUG)
    add_test(NAME ${_t} COMMAND ${_t} WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
  endforeach()
endif()

# ------------------------------------------------------------------- python --
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

  # 2) Cython transpile. Cython derives the module name (PyInit_*) from the
  #    source filename, so stage classy.pyx in the build dir under the
  #    requested module name. (cclassy.pxd is found via -I <build dir>.)
  set(_class_pyx ${CMAKE_CURRENT_BINARY_DIR}/${CLASS_PYTHON_MODULE_NAME}.pyx)
  add_custom_command(
    OUTPUT ${_class_pyx}
    COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_SOURCE_DIR}/classy.pyx ${_class_pyx}
    DEPENDS ${PROJECT_SOURCE_DIR}/classy.pyx
    COMMENT "Staging ${CLASS_PYTHON_MODULE_NAME}.pyx"
    VERBATIM
  )
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${CLASS_PYTHON_MODULE_NAME}.cpp
    COMMAND Python::Interpreter -m cython -3 --cplus
            -I ${CMAKE_CURRENT_BINARY_DIR}
            -o ${CMAKE_CURRENT_BINARY_DIR}/${CLASS_PYTHON_MODULE_NAME}.cpp
            ${_class_pyx}
    DEPENDS ${_class_pyx} ${CMAKE_CURRENT_BINARY_DIR}/cclassy.pxd
    COMMENT "Cythonizing ${CLASS_PYTHON_MODULE_NAME}.pyx"
    VERBATIM
  )

  # 3) The extension module.
  python_add_library(${CLASS_PYTHON_MODULE_NAME} MODULE
    ${CMAKE_CURRENT_BINARY_DIR}/${CLASS_PYTHON_MODULE_NAME}.cpp WITH_SOABI)
  target_link_libraries(${CLASS_PYTHON_MODULE_NAME} PRIVATE classpp Python::NumPy)
  target_compile_definitions(${CLASS_PYTHON_MODULE_NAME} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
  set_target_properties(${CLASS_PYTHON_MODULE_NAME} PROPERTIES
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
  )

  # 4) Wheel layout — identical to the old setup.py wheel:
  #    classy.<soabi>.so at the root, data under classy/{bbn,hyrec}.
  #    (The data stays under classy/ even for a renamed module — that matches
  #    the old classyref wheels and classy.pyx's runtime lookup.)
  if(SKBUILD)
    install(TARGETS ${CLASS_PYTHON_MODULE_NAME} LIBRARY DESTINATION .)
    file(GLOB _class_bbn_data ${PROJECT_SOURCE_DIR}/bbn/*.dat)
    file(GLOB _class_hyrec_data ${PROJECT_SOURCE_DIR}/hyrec/*.dat)
    install(FILES ${_class_bbn_data} DESTINATION classy/bbn)
    install(FILES ${_class_hyrec_data} DESTINATION classy/hyrec)
  endif()
endif()
