# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# GPU Boundary Checking Test
add_executable(test_gpu_boundaries test_gpu_boundaries.cpp)
target_link_libraries(test_gpu_boundaries PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_gpu_boundaries PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME test_gpu_boundaries COMMAND test_gpu_boundaries)

# GPU Space Charge P³M Test
add_executable(test_gpu_space_charge test_gpu_space_charge.cpp)
target_link_libraries(test_gpu_space_charge PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_gpu_space_charge PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME test_gpu_space_charge COMMAND test_gpu_space_charge)

# Adaptive Field Interpolation Test
add_executable(test_adaptive_interpolation test_adaptive_interpolation.cpp)
target_link_libraries(test_adaptive_interpolation PRIVATE icarion_core Catch2::Catch2WithMain)
target_include_directories(test_adaptive_interpolation PRIVATE ${CMAKE_SOURCE_DIR}/src)
add_test(NAME test_adaptive_interpolation COMMAND test_adaptive_interpolation)
