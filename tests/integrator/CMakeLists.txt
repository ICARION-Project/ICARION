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
