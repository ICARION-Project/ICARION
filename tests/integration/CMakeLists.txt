# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: 2026 Christoph Schaefer

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

find_program(Python3_EXECUTABLE NAMES python3 python)
if(Python3_EXECUTABLE)
    add_test(NAME UserAnnotationSimulationCLI
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_annotation_cli.py
                --repo-root ${CMAKE_SOURCE_DIR}
                --icarion-bin $<TARGET_FILE:icarion_main>)
endif()
