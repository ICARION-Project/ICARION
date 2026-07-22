# =====================================================================
#   ICARION v1.1.1 - Core Build
# =====================================================================

cmake_minimum_required(VERSION 3.16)

message(STATUS "Configuring ICARION Core...")

#
# Library: icarion_core
#
add_library(icarion_core STATIC)

target_compile_features(icarion_core PUBLIC cxx_std_20)

# GCC 13 false-positive -Warray-bounds warning from spdlog/fmt headers.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
    target_compile_options(icarion_core PRIVATE -Wno-array-bounds -Wno-stringop-overflow)
endif()

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
    # Suppress pedantic line-marker noise from NVCC host compilation.
    target_compile_options(icarion_core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-Wno-pedantic>)

    # Link cuFFT for GPU Space Charge
    if(ICARION_CUFFT_LIBRARY)
        target_link_libraries(icarion_core PUBLIC ${ICARION_CUFFT_LIBRARY})
    endif()

    if(ICARION_HAVE_NCCL)
        target_compile_definitions(icarion_core PUBLIC HAVE_NCCL)
        target_include_directories(icarion_core PRIVATE ${ICARION_NCCL_INCLUDE_DIR})
        target_link_libraries(icarion_core PUBLIC ${ICARION_NCCL_LIBRARY})
    endif()
endif()

# IO + Logging
file(GLOB_RECURSE ICARION_IO_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/io/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/io/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/log/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/log/*.h"
)
# Exclude legacy configuration loader (replaced by ConfigLoader)
list(FILTER ICARION_IO_SRC EXCLUDE REGEX "configurationLoader\\.cpp$")

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
    "${CMAKE_CURRENT_SOURCE_DIR}/core/types/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/param/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/param/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.h"
)
# Exclude legacy paramUtils (removed)
list(FILTER ICARION_UTILS_SRC EXCLUDE REGEX "paramUtils\\.cpp$")
list(FILTER ICARION_UTILS_SRC EXCLUDE REGEX "paramUtils\\.h$")
if(NOT USE_GPU_ACCEL)
    # Do not compile GPU helper sources when CUDA is off (avoids cuda_runtime.h include errors)
    list(FILTER ICARION_UTILS_SRC EXCLUDE REGEX "/core/types/gpu/")
endif()

# SimulationEngine (empty for now)
file(GLOB_RECURSE ICARION_SIMULATION_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/core/simulation/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/simulation/*.h"
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
    # Explicit new field models to avoid glob reconfigure misses
    core/config/types/FieldProviderModel.cpp
    core/io/UserAnnotation.cpp
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

add_executable(icarion_main 
    main/main.cpp
    main/setup/PhysicsSetup.cpp
)
target_link_libraries(icarion_main PRIVATE icarion_core)

# CCS precompute tool
# CCS precompute tool
add_executable(ccs_precompute
    ../tools/ccs_precompute.cpp
)
target_link_libraries(ccs_precompute PRIVATE icarion_core)

# EHSS samples precompute tool
add_executable(ehss_samples_precompute
    ../tools/ehss_samples_precompute.cpp
)
target_link_libraries(ehss_samples_precompute PRIVATE icarion_core)

# EHSS offline sample precompute tool
add_executable(ehss_offline_precompute
    ../tools/ehss_offline_precompute.cpp
)
target_link_libraries(ehss_offline_precompute PRIVATE icarion_core)

# Interaction-potential offline sample precompute tool
add_executable(interaction_potential_precompute
    ../tools/interaction_potential_precompute.cpp
)
target_link_libraries(interaction_potential_precompute PRIVATE icarion_core)

# ---------------------------------------------------------------------
# Installable runtime package
# ---------------------------------------------------------------------

install(TARGETS
        icarion_main
        ccs_precompute
        ehss_samples_precompute
        ehss_offline_precompute
        interaction_potential_precompute
        RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
        COMPONENT runtime)

if(WIN32)
    set(ICARION_COMMAND_ALIAS "icarion.exe")
else()
    set(ICARION_COMMAND_ALIAS "icarion")
endif()

install(PROGRAMS "$<TARGET_FILE:icarion_main>"
        DESTINATION "${CMAKE_INSTALL_BINDIR}"
        RENAME "${ICARION_COMMAND_ALIAS}"
        COMPONENT runtime)

message(STATUS "ICARION Core configuration completed.")
