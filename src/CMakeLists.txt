# =====================================================================
#   ICARION v1.0 — Core Build
# =====================================================================

cmake_minimum_required(VERSION 3.16)

message(STATUS "Configuring ICARION Core...")

#
# Library: icarion_core
#
add_library(icarion_core STATIC)

target_compile_features(icarion_core PUBLIC cxx_std_20)

# ---------------------------------------------------------------------
# Source groups for the core library
# ---------------------------------------------------------------------

# Physics
file(GLOB_RECURSE ICARION_PHYSICS_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/physics/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/physics/*.h"
)

# Integrator
file(GLOB_RECURSE ICARION_INTEGRATOR_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/integrator/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/integrator/*.h"
)

# GPU (optional)
if(USE_GPU_ACCEL)
    file(GLOB_RECURSE ICARION_GPU_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/core/gpu/*.cu"
        "${CMAKE_CURRENT_SOURCE_DIR}/core/gpu/*.cuh"
        "${CMAKE_CURRENT_SOURCE_DIR}/core/gpu/*.cpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/core/gpu/*.h"
    )

    if(NOT ICARION_HAVE_NCCL)
        list(FILTER ICARION_GPU_SRC EXCLUDE REGEX "/multi_gpu_manager\\.cu$")
        message(STATUS "Excluding multi_gpu_manager.cu (NCCL not available)")
    endif()

    target_sources(icarion_core PRIVATE ${ICARION_GPU_SRC})

    if(ICARION_HAVE_NCCL)
        target_compile_definitions(icarion_core PUBLIC HAVE_NCCL)
        target_include_directories(icarion_core PRIVATE ${ICARION_NCCL_INCLUDE_DIR})
        target_link_libraries(icarion_core PUBLIC ${ICARION_NCCL_LIBRARY})
    endif()
endif()

# IO
file(GLOB_RECURSE ICARION_IO_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/io/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/io/*.h"
)

# Config
file(GLOB_RECURSE ICARION_CONFIG_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/config/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/config/*.h"
)

# Utils + Types
file(GLOB_RECURSE ICARION_UTILS_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/utils/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/utils/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/types/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/param/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/param/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.h"
)

# SimulationEngine (empty for now)
file(GLOB_RECURSE ICARION_SIMULATION_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/simulation/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/simulation/*.h"
)

# FieldSolver utilities (for field update API only)
set(ICARION_FIELDSOLVER_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/fieldsolver/utils/field_update_api.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/fieldsolver/utils/field_update_api.h"
)

# ---------------------------------------------------------------------
# Add sources to library
# ---------------------------------------------------------------------

target_sources(icarion_core PRIVATE
    ${ICARION_PHYSICS_SRC}
    ${ICARION_INTEGRATOR_SRC}
    ${ICARION_IO_SRC}
    ${ICARION_CONFIG_SRC}
    ${ICARION_UTILS_SRC}
    ${ICARION_SIMULATION_SRC}
    ${ICARION_FIELDSOLVER_SRC}
)

# Include root for headers
target_include_directories(icarion_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/core
)

# Dependencies (HDF5, jsoncpp, etc.)
target_link_libraries(icarion_core
    PUBLIC
        ${ICARION_CORE_DEPS}
)

# Enable warnings and sanitizers if configured
if(ICARION_ENABLE_SANITIZERS)
    add_sanitizers(icarion_core)
endif()

# ---------------------------------------------------------------------
# Build example main (optional)
# ---------------------------------------------------------------------

add_executable(icarion_main main/main.cpp)
target_link_libraries(icarion_main PRIVATE icarion_core)

message(STATUS "ICARION Core configuration completed.")
