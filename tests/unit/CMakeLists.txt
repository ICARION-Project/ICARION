# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Contributors

# Unit test for IonEnsemble (SoA data structure)

add_executable(test_ion_ensemble test_ion_ensemble.cpp)

target_link_libraries(test_ion_ensemble
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_ion_ensemble
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Register with CTest
add_test(NAME test_ion_ensemble COMMAND test_ion_ensemble)

# Add to test suite
set_tests_properties(test_ion_ensemble PROPERTIES
    LABELS "unit;soa"
    TIMEOUT 30
)
