# Copyright (c) 2026, AthenaPK 4d_modes contributors.
# Licensed under the 3-Clause License (the "LICENSE");

target_sources(athenaPK PRIVATE
    mode_tables.cpp
    mode_transform.cpp
    em4d_modes.cpp
    plasma4d_modes.cpp
    diagnostics4d.cpp
    pgen_harris4d.cpp
    package4d_modes.cpp
)

if(AthenaPK_ENABLE_4D_MODES_UNIT_TESTS)
  add_executable(modes4d_unit_tests
      mode_tables_tests.cpp
      mode_tables.cpp
      mode_transform.cpp
  )
  target_compile_features(modes4d_unit_tests PRIVATE cxx_std_17)
  if(BUILD_TESTING)
    add_test(NAME modes4d_unit_tests COMMAND modes4d_unit_tests)
  endif()
endif()
