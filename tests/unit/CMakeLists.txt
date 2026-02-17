# SPDX-License-Identifier: GPL-3.0-only
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

# Unit test for SimulationEngine SoA integration (Phase 2)

add_executable(test_simulation_engine_soa test_simulation_engine_soa.cpp)

target_link_libraries(test_simulation_engine_soa
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_simulation_engine_soa
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Register with CTest
add_test(NAME test_simulation_engine_soa COMMAND test_simulation_engine_soa)

# Add to test suite
set_tests_properties(test_simulation_engine_soa PROPERTIES
    LABELS "unit;soa;integration"
    TIMEOUT 60
)

# Benchmark for SoA performance (Phase 3)

add_executable(benchmark_soa_performance benchmark_soa_performance.cpp)

target_link_libraries(benchmark_soa_performance
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(benchmark_soa_performance
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Register with CTest (marked as benchmark, runs only when requested)
add_test(NAME benchmark_soa_performance COMMAND benchmark_soa_performance "[performance]")

set_tests_properties(benchmark_soa_performance PROPERTIES
    LABELS "benchmark;soa;performance"
    TIMEOUT 120
)

# Benchmark for OpenMP SoA scaling (Phase 6)

add_executable(benchmark_openmp_soa benchmark_openmp_soa.cpp)

target_link_libraries(benchmark_openmp_soa
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
        OpenMP::OpenMP_CXX
)

target_include_directories(benchmark_openmp_soa
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src
)

# Register with CTest (marked as benchmark, runs only when requested)
add_test(NAME benchmark_openmp_soa COMMAND benchmark_openmp_soa "[openmp]")

set_tests_properties(benchmark_openmp_soa PROPERTIES
    LABELS "benchmark;soa;openmp;performance"
    TIMEOUT 180
)

# Domain lookup test
add_executable(test_domain_manager_lookup test_domain_manager_lookup.cpp)
target_link_libraries(test_domain_manager_lookup PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_domain_manager_lookup PRIVATE ${PROJECT_SOURCE_DIR}/src)
add_test(NAME test_domain_manager_lookup COMMAND test_domain_manager_lookup)
set_tests_properties(test_domain_manager_lookup PROPERTIES LABELS "unit;domain" TIMEOUT 30)

# Domain lookup benchmark (not run by default)
add_executable(benchmark_domain_lookup benchmark_domain_lookup.cpp)
target_link_libraries(benchmark_domain_lookup PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(benchmark_domain_lookup PRIVATE ${PROJECT_SOURCE_DIR}/src)
add_test(NAME benchmark_domain_lookup COMMAND benchmark_domain_lookup "[benchmark][domain]")
set_tests_properties(benchmark_domain_lookup PROPERTIES LABELS "benchmark;domain" TIMEOUT 60)
