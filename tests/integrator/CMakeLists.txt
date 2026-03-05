# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: 2025 Christoph Schaefer

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

add_executable(test_simulation_engine_per_ion_dt test_simulation_engine_per_ion_dt.cpp)
target_link_libraries(test_simulation_engine_per_ion_dt PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_simulation_engine_per_ion_dt PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_simulation_engine_per_ion_dt COMMAND test_simulation_engine_per_ion_dt)

add_executable(test_simulation_engine_adaptive_sc test_simulation_engine_adaptive_sc.cpp)
target_link_libraries(test_simulation_engine_adaptive_sc PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_simulation_engine_adaptive_sc PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_simulation_engine_adaptive_sc COMMAND test_simulation_engine_adaptive_sc)

add_executable(test_domain_geometry test_domain_geometry.cpp)
target_link_libraries(test_domain_geometry PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_domain_geometry PRIVATE ${CMAKE_SOURCE_DIR}/src)

add_test(NAME test_domain_geometry COMMAND test_domain_geometry)

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
    target_compile_definitions(test_gpu_field_interpolation PRIVATE ICARION_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    
    add_test(NAME test_gpu_field_interpolation COMMAND test_gpu_field_interpolation)
    set_tests_properties(test_gpu_field_interpolation PROPERTIES 
        LABELS "gpu;fields;integration"
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    
    # GPU RK45 integration test
    add_executable(test_gpu_rk45 test_gpu_rk45.cpp)
    target_link_libraries(test_gpu_rk45 PRIVATE icarion_core)
    target_include_directories(test_gpu_rk45 PRIVATE ${CMAKE_SOURCE_DIR}/src)
    
    add_test(NAME test_gpu_rk45 COMMAND test_gpu_rk45)
    set_tests_properties(test_gpu_rk45 PROPERTIES LABELS "gpu;rk45;integration")
    
    # GPU Boris pusher integration test
    add_executable(test_gpu_boris test_gpu_boris.cpp)
    target_link_libraries(test_gpu_boris PRIVATE icarion_core)
    target_include_directories(test_gpu_boris PRIVATE ${CMAKE_SOURCE_DIR}/src)
    
    add_test(NAME test_gpu_boris COMMAND test_gpu_boris)
    set_tests_properties(test_gpu_boris PROPERTIES LABELS "gpu;boris;integration")
    
    # CPU/GPU parity test for RK45 and Boris
    add_executable(test_rk45_boris_parity test_rk45_boris_parity.cpp)
    target_link_libraries(test_rk45_boris_parity PRIVATE icarion_core Catch2::Catch2WithMain)
    target_include_directories(test_rk45_boris_parity PRIVATE ${CMAKE_SOURCE_DIR}/src)
    
    add_test(NAME test_rk45_boris_parity COMMAND test_rk45_boris_parity)
    set_tests_properties(test_rk45_boris_parity PROPERTIES LABELS "gpu;parity;rk45;boris")
endif()
