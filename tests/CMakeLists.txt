# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Enable testing
enable_testing()

# Find Catch2
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found, config tests will not be built")
    return()
endif()

# Add test subdirectories
add_subdirectory(config)
add_subdirectory(io)
add_subdirectory(utils)
add_subdirectory(physics)
add_subdirectory(integrator)
add_subdirectory(integration)
add_subdirectory(instruments)
