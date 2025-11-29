# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Integration strategy tests

add_executable(test_rk4_strategy test_rk4_strategy.cpp)
target_link_libraries(test_rk4_strategy PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_rk4_strategy PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_rk4_strategy COMMAND test_rk4_strategy)

add_executable(test_rk45_strategy test_rk45_strategy.cpp)
target_link_libraries(test_rk45_strategy PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_rk45_strategy PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_rk45_strategy COMMAND test_rk45_strategy)

add_executable(test_boris_strategy test_boris_strategy.cpp)
target_link_libraries(test_boris_strategy PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_boris_strategy PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_boris_strategy COMMAND test_boris_strategy)

# Phase 5A: SimulationEngine tests

add_executable(test_domain_manager test_domain_manager.cpp)
target_link_libraries(test_domain_manager PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_domain_manager PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_domain_manager COMMAND test_domain_manager)

add_executable(test_output_manager test_output_manager.cpp)
target_link_libraries(test_output_manager PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_output_manager PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_output_manager COMMAND test_output_manager)

add_executable(test_simulation_engine test_simulation_engine.cpp)
target_link_libraries(test_simulation_engine PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_simulation_engine PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_simulation_engine COMMAND test_simulation_engine)

# GPU Integration Tests (conditional on USE_GPU_ACCEL)

if (USE_GPU_ACCEL)
    add_executable(test_gpu_integration test_gpu_integration.cpp)
    target_link_libraries(test_gpu_integration PRIVATE icarion_core Catch2::Catch2WithMain)
    target_include_directories(test_gpu_integration PRIVATE ${CMAKE_SOURCE_DIR}/src)
    
    add_test(NAME test_gpu_integration COMMAND test_gpu_integration)
    set_tests_properties(test_gpu_integration PROPERTIES LABELS "gpu;integration")
    
    add_executable(test_cpu_gpu_parity test_cpu_gpu_parity.cpp)
    target_link_libraries(test_cpu_gpu_parity PRIVATE icarion_core Catch2::Catch2WithMain)
    target_include_directories(test_cpu_gpu_parity PRIVATE ${CMAKE_SOURCE_DIR}/src)
    
    add_test(NAME test_cpu_gpu_parity COMMAND test_cpu_gpu_parity)
    set_tests_properties(test_cpu_gpu_parity PROPERTIES LABELS "gpu;parity;integration")
    
    add_executable(test_gpu_field_interpolation test_gpu_field_interpolation.cpp)
    target_link_libraries(test_gpu_field_interpolation PRIVATE icarion_core Catch2::Catch2WithMain)
    target_include_directories(test_gpu_field_interpolation PRIVATE ${CMAKE_SOURCE_DIR}/src)
    
    add_test(NAME test_gpu_field_interpolation COMMAND test_gpu_field_interpolation)
    set_tests_properties(test_gpu_field_interpolation PROPERTIES 
        LABELS "gpu;fields;integration"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
