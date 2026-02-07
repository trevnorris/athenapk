# AthenaPK - a performance portable block structured AMR MHD code
# Copyright (c) 2020-2023, Athena Parthenon Collaboration. All rights reserved.
# Licensed under the 3-Clause License (the "LICENSE");

cmake_minimum_required(VERSION 3.13)

# Imports machine-specific configuration
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/parthenon/CMakeLists.txt)
  include(external/parthenon/cmake/MachineCfg.cmake)
else()
  message(WARNING
    "Could not find MachineCfg.cmake to process machine config file."
    "If you're using Parthenon as submodule and see this message, please report this bug on GitLab."
    "Using a prebuilt Parthenon library with AthenaPK is currently not tested."
    "For further directions open an issue on GitLab.")
endif()

project(athenaPK VERSION 25.05 LANGUAGES CXX)

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)

# make paths in target_source() absolute
# used for sources in src/pgen folder
# new behavior introduces in CMake 3.13 that's why this is the min. version above
cmake_policy(SET CMP0076 NEW)

# Don't allow in-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR
      "You cannot build in a source directory (or any directory with a CMakeLists.txt file). "
      "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# If the user doesn't specify a build type, prefer Release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include(cmake/Format.cmake)

# Kokkos itself will actually be discovered by the Parthenon package.
# Just need to set a hint to use the AthenaPK submodule and not the Parthenon one.
set(Kokkos_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/external/Kokkos)

# The MPI package found (and it's variables) by the Parthenon package are not available
# in the parent scope here so we need to look for MPI here as it is used in the CMake
# logic for setting up tests.
set(NUM_MPI_PROC_TESTING "4" CACHE STRING "Number of mpi processors to use when running tests with MPI")
if (NOT PARTHENON_DISABLE_MPI)
  find_package(MPI COMPONENTS CXX)
  if (NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is required but couldn't be found. "
    "If you want to build Parthenon without MPI, please rerun CMake with -DPARTHENON_DISABLE_MPI=ON")
  endif()
endif()

option(AthenaPK_ENABLE_TESTING "Enable AthenaPK test" ON)
option(AthenaPK_ENABLE_4D_MODES "Enable experimental 4D x modes package scaffold" ON)
option(AthenaPK_ENABLE_4D_MODES_UNIT_TESTS "Build 4D modes unit-test executable" ON)
option(AthenaPK_ENABLE_4D_MODES_REGRESSION_TARGET "Add tuned 4D Harris regression target" ON)
set(PARTHENON_ENABLE_PYTHON_MODULE_CHECK ${AthenaPK_ENABLE_TESTING} CACHE BOOL "Check if local python version contains all modules required for running tests.")

set(PARTHENON_ENABLE_TESTING OFF CACHE BOOL "Disable Parthenon testing.")
set(PARTHENON_DISABLE_OPENMP ON CACHE BOOL "Disable OpenMP")
set(PARTHENON_DISABLE_EXAMPLES ON CACHE BOOL "Don't build Parthenon examples.")
set(PARTHENON_DISABLE_SPARSE ON CACHE BOOL "Disable sparse (not used in AthenaPK yet)")

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/external/parthenon/CMakeLists.txt)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/parthenon parthenon)
else()
  find_package(parthenon REQUIRED)
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_subdirectory(src)

if(AthenaPK_ENABLE_4D_MODES AND AthenaPK_ENABLE_4D_MODES_REGRESSION_TARGET)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_custom_target(
      modes4d_kk_coulomb_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/kk_coulomb_yukawa_regression.py
      USES_TERMINAL
      COMMENT "Running Gaussian KK Coulomb+Yukawa regression checks")
    add_custom_target(
      modes4d_controlled_parity_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_controlled_parity_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_controlled_parity_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_controlled_parity_regression
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running controlled-limit parity regression checks")
    add_custom_target(
      modes4d_controlled_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_controlled_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_controlled_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_controlled_regression
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running controlled-limit modes4d regression smoke checks")
    add_custom_target(
      modes4d_scalar_pulse_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_scalar_pulse_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/em4d_scalar_pulse_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_scalar_pulse_regression
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running scalar-channel 4D pulse regression checks")
    add_custom_target(
      modes4d_em_conservative_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_em_conservative_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_conservative_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_em_conservative_regression
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running conservative-EM transport regression checks (scalar + short Harris)")
    add_custom_target(
      modes4d_em_timelike_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_em_timelike_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_timelike_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_em_timelike_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_em_timelike_regression/outputs
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running targeted time-like EM mixed-coupling regression checks")
    add_custom_target(
      modes4d_em_conservative_timelike_diag_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_em_conservative_timelike_diag_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_conservative_timelike_diag_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_em_conservative_timelike_diag_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_em_conservative_timelike_diag_regression/outputs
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running conservative-path time-like EM diagnostics regression checks")
    add_custom_target(
      modes4d_em_conservative_transport_ledger_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_em_conservative_transport_ledger_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_conservative_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_em_conservative_transport_ledger_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_em_conservative_transport_ledger_regression/outputs
        --harris-tlim 0.05
        --harris-nlim 500
        --transport-norm-tol 1.0e-2
        --transport-abs-rate-tol 1.0e-10
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running strict conservative EM transport/source ledger regression checks")
    add_custom_target(
      modes4d_harris_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_regression_tuned.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_regression
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running tuned 4D Harris regression (closure + activity checks)")
    add_custom_target(
      modes4d_harris_pseudospectral_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_pseudospectral_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_regression_tuned.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_pseudospectral_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_harris_pseudospectral_regression/outputs
        --scan-controlled-arg modes4d/plasma_pseudospectral_transport=false
        --scan-full-arg modes4d/plasma_pseudospectral_transport=true
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running tuned 4D Harris regression with explicit pseudospectral transport gating")
    add_custom_target(
      modes4d_harris_nw_convergence_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_nw_convergence_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_nw_convergence_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_nw_convergence_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_harris_nw_convergence_regression/outputs
        --check-transport-closure
        --fail-on-check
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running tuned 4D Harris N_w-convergence regression")
    add_custom_target(
      modes4d_harris_pressure_experimental_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_pressure_experimental_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_regression_tuned.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_pressure_experimental_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_harris_pressure_experimental_regression/outputs
        --scan-full-arg modes4d/plasma_pressure_transport_gain=1.0e-3
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running experimental 4D Harris pressure-transport regression (full-case gain=1e-3)")
    add_custom_target(
      modes4d_harris_pressure_window_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_pressure_window_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_pressure_window_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_pressure_window_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_harris_pressure_window_regression/outputs
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running 4D Harris pressure-window regression (finite nonzero gains)")
    add_custom_target(
      modes4d_harris_wflux_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_wflux_regression
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_wflux_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_wflux_regression
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_harris_wflux_regression/outputs
        --enabled-gain 1.0e-1
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running tuned 4D Harris regression with projected two-fluid w-flux source enabled")
    add_custom_target(
      modes4d_harris_lundquist_scan
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_harris_lundquist_scan
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_lundquist_scan.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_harris_lundquist_scan
        --controlled-input ${CMAKE_CURRENT_SOURCE_DIR}/inputs/harris_4d_controlled.in
        --full-input ${CMAKE_CURRENT_SOURCE_DIR}/inputs/harris_4d_full.in
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_harris_lundquist_scan/outputs
        --fail-on-check
        --check-transport-closure
        --fail-on-scan-check
        --scan-min-abs-corr corr_logS_full_final_psi0_span=7.0e-1
        --scan-min-abs-corr corr_logS_full_final_s_leak_abs=5.0e-1
        --scan-min-abs-corr corr_logS_full_final_jw_ew_abs=5.0e-1
        --scan-min-abs-corr corr_logS_full_final_mixed_ew2=7.0e-1
        --scan-min-abs-corr corr_logS_full_final_em_leak_w_abs=2.0e-1
        --scan-min-abs-corr corr_logS_full_final_helicity_sub_abs=5.0e-1
        --scan-min-abs-corr corr_logS_full_final_edotb_sub_abs=5.0e-1
        --scan-max-corr corr_logS_full_final_helicity_sub_abs=-2.0e-1
        --scan-max-corr corr_logS_full_final_edotb_sub_abs=-2.0e-1
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running 4D Harris Lundquist scan (controlled vs full)")
    add_custom_target(
      modes4d_phase10_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/controlled_smoke
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_controlled_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/controlled_smoke
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/controlled_parity
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_controlled_parity_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/controlled_parity
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/kk_coulomb_yukawa_regression.py
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/scalar_pulse
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/em4d_scalar_pulse_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/scalar_pulse
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_timelike
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_timelike_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_timelike
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_timelike/outputs
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative_timelike_diag
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_conservative_timelike_diag_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative_timelike_diag
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative_timelike_diag/outputs
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_conservative_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative_transport_ledger
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/modes4d_em_conservative_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative_transport_ledger
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/em_conservative_transport_ledger/outputs
        --harris-tlim 0.05
        --harris-nlim 500
        --transport-norm-tol 1.0e-2
        --transport-abs-rate-tol 1.0e-10
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_regression_tuned.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris_pseudospectral
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_regression_tuned.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris_pseudospectral
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris_pseudospectral/outputs
        --scan-controlled-arg modes4d/plasma_pseudospectral_transport=false
        --scan-full-arg modes4d/plasma_pseudospectral_transport=true
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris_nw_convergence
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_nw_convergence_regression.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris_nw_convergence
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/harris_nw_convergence/outputs
        --check-transport-closure
        --fail-on-check
      COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/lundquist
      COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/scripts/harris_lundquist_scan.py
        --binary $<TARGET_FILE:athenaPK>
        --workdir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/lundquist
        --controlled-input ${CMAKE_CURRENT_SOURCE_DIR}/inputs/harris_4d_controlled.in
        --full-input ${CMAKE_CURRENT_SOURCE_DIR}/inputs/harris_4d_full.in
        --output-dir ${CMAKE_BINARY_DIR}/modes4d_phase10_regression/lundquist/outputs
        --fail-on-check
        --check-transport-closure
        --fail-on-scan-check
        --scan-min-abs-corr corr_logS_full_final_psi0_span=7.0e-1
        --scan-min-abs-corr corr_logS_full_final_s_leak_abs=5.0e-1
        --scan-min-abs-corr corr_logS_full_final_jw_ew_abs=5.0e-1
        --scan-min-abs-corr corr_logS_full_final_mixed_ew2=7.0e-1
        --scan-min-abs-corr corr_logS_full_final_em_leak_w_abs=2.0e-1
        --scan-min-abs-corr corr_logS_full_final_helicity_sub_abs=5.0e-1
        --scan-min-abs-corr corr_logS_full_final_edotb_sub_abs=5.0e-1
        --scan-max-corr corr_logS_full_final_helicity_sub_abs=-2.0e-1
        --scan-max-corr corr_logS_full_final_edotb_sub_abs=-2.0e-1
      DEPENDS athenaPK
      USES_TERMINAL
      COMMENT "Running Phase-10 4D validation bundle (controlled, scalar, KK, Harris, Lundquist)")
  else()
    message(STATUS "Python3 interpreter not found; skipping modes4d_harris_regression target.")
  endif()
endif()

if (AthenaPK_ENABLE_TESTING)
  include(CTest)
  add_subdirectory(tst/regression)
endif()
