# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Config system tests
add_executable(test_config_loader
    test_config_loader.cpp
)

target_link_libraries(test_config_loader
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_config_loader
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME ConfigLoader COMMAND test_config_loader)

# Legacy adapter tests
add_executable(test_legacy_adapter
    test_legacy_adapter.cpp
)

target_link_libraries(test_legacy_adapter
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_legacy_adapter
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME LegacyAdapter COMMAND test_legacy_adapter)
