# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Integration tests - test complete workflows

add_executable(test_main_integration test_main_integration.cpp)
target_link_libraries(test_main_integration
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)
target_include_directories(test_main_integration
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME test_main_integration COMMAND test_main_integration)
