# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Collision handler tests

# Test: Collision Handler Factory
add_executable(test_collision_factory test_collision_factory.cpp)
target_link_libraries(test_collision_factory PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_collision_factory PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME CollisionFactory COMMAND test_collision_factory)
