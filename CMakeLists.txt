cmake_minimum_required(VERSION 3.21)
project(CLASSpp LANGUAGES C CXX)

if(SKBUILD)
  set(_class_python_default ON)
else()
  set(_class_python_default OFF)
endif()
option(CLASS_BUILD_PYTHON "Build the classy Python extension" ${_class_python_default})

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
target_include_directories(classpp
  PUBLIC
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
  add_executable(test-bisection tools/bisection_test.cpp)
  add_executable(test-photons species/photons_formula_test.cpp)
  foreach(_t IN ITEMS test-parser test-bisection test-photons)
    target_link_libraries(${_t} PRIVATE classpp)
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
