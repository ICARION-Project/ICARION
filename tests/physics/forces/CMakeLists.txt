# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Physics Forces Tests
add_executable(test_force_registry
    test_force_registry.cpp
)

target_link_libraries(test_force_registry
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_force_registry
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME ForceRegistry COMMAND test_force_registry)
set_tests_properties(ForceRegistry PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
