# =====================================================================
#   ICARION - Ion Collision And Reaction IntegratiON
# =====================================================================
#   Modular C++ framework for simulating ion trajectories in
#   user-defined electric fields and background gas environments.
#
#   Build Configuration:
#     cmake -B build -DCMAKE_BUILD_TYPE=Release
#     cmake --build build -j$(nproc)
#
#   GPU Support:
#     cmake -B build -DUSE_GPU_ACCEL=ON
#
# =====================================================================

cmake_minimum_required(VERSION 3.16)

# Project name aligned with repository. Enable CUDA language so that
# check_language(CUDA) works reliably when -DUSE_GPU_ACCEL=ON is passed.
project(ICARION VERSION 1.0.0 LANGUAGES C CXX)

# ─────────────────────────────────────────────────────────────────────
# Build Configuration Options
# ─────────────────────────────────────────────────────────────────────

option(ICARION_BUILD_CORE_ONLY 
       "Build only ICARION Core (no FieldSolver or Optimizer)" 
       OFF)

option(USE_GPU_ACCEL 
       "Enable CUDA GPU acceleration" 
       OFF)

# Propagate core-only flag to compiler
if(ICARION_BUILD_CORE_ONLY)
    add_compile_definitions(ICARION_BUILD_CORE_ONLY)
    message(STATUS "ICARION Core-Only build enabled (no FieldSolver/Optimizer)")
endif()

# ─────────────────────────────────────────────────────────────────────
# Version and Git Information
# ─────────────────────────────────────────────────────────────────────

# Set version compile definition
add_compile_definitions(ICARION_VERSION="${PROJECT_VERSION}")

# Get git commit hash
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
    )
    if(NOT GIT_HASH)
        set(GIT_HASH "unknown")
    endif()
else()
    set(GIT_HASH "unknown")
endif()
add_compile_definitions(GIT_HASH="${GIT_HASH}")
message(STATUS "Git commit hash: ${GIT_HASH}")

# ─────────────────────────────────────────────────────────────────────
# CMake Modules (Modular Configuration)
# ─────────────────────────────────────────────────────────────────────

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompilerSettings.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Dependencies.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Targets.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CUDAConfig.cmake)

# ─────────────────────────────────────────────────────────────────────
# Global Include Paths
# ─────────────────────────────────────────────────────────────────────

include_directories(${PROJECT_SOURCE_DIR}/src)

# ─────────────────────────────────────────────────────────────────────
# Testing (must be at top level for CTest to work)
# ─────────────────────────────────────────────────────────────────────

enable_testing()

# ─────────────────────────────────────────────────────────────────────
# Build Subdirectories
# ─────────────────────────────────────────────────────────────────────

add_subdirectory(src)
add_subdirectory(tests)

# ─────────────────────────────────────────────────────────────────────
# Build Summary
# ─────────────────────────────────────────────────────────────────────

message(STATUS "")
message(STATUS "========================================")
message(STATUS "ICARION v${PROJECT_VERSION} Configuration")
message(STATUS "========================================")
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "Compiler:          ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ standard:      ${CMAKE_CXX_STANDARD}")
message(STATUS "Git commit:        ${GIT_HASH}")
message(STATUS "")
message(STATUS "Features:")
if(OpenMP_CXX_FOUND)
    message(STATUS "  OpenMP:          ${OpenMP_CXX_VERSION}")
else()
    message(STATUS "  OpenMP:          Not found")
endif()
if(USE_GPU_ACCEL)
    message(STATUS "  GPU (CUDA):      Enabled")
else()
    message(STATUS "  GPU (CUDA):      Disabled")
endif()
if(ICARION_BUILD_CORE_ONLY)
    message(STATUS "  Build mode:      Core-only")
else()
    message(STATUS "  Build mode:      Full")
endif()
message(STATUS "")
message(STATUS "Install prefix:    ${CMAKE_INSTALL_PREFIX}")
message(STATUS "========================================")
message(STATUS "")

