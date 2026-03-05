# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: 2026 Christoph Schaefer

# Enable testing
enable_testing()

# Find Catch2
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    message(STATUS "Catch2 not found, config tests will not be built")
    return()
endif()

# Add test subdirectories
add_subdirectory(unit)
add_subdirectory(config)
add_subdirectory(io)
add_subdirectory(utils)
add_subdirectory(physics)
add_subdirectory(integrator)
add_subdirectory(integration)
add_subdirectory(instruments)

# GPU tests (if GPU acceleration enabled)
if(USE_GPU_ACCEL)
    add_subdirectory(gpu)
endif()
