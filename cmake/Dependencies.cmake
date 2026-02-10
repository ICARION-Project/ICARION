# =====================================================================
# ICARION - Dependencies
# =====================================================================

message(STATUS "Loading Dependencies.cmake")

# JSON (jsoncpp) - minimum version 1.9.0
find_package(jsoncpp 1.9.0 CONFIG REQUIRED)
set(ICARION_CORE_DEPS jsoncpp_lib)

# HDF5 - minimum version 1.10.0
find_package(HDF5 1.10.0 COMPONENTS C CXX REQUIRED)
list(APPEND ICARION_CORE_DEPS HDF5::HDF5)

# Threads
find_package(Threads REQUIRED)
list(APPEND ICARION_CORE_DEPS Threads::Threads)

# Math / BLAS (optional, but harmless)
find_package(BLAS REQUIRED)
list(APPEND ICARION_CORE_DEPS BLAS::BLAS)

# cxxopts (CLI parsing)
find_package(cxxopts CONFIG QUIET)
if(NOT cxxopts_FOUND)
    message(STATUS "cxxopts not found, using FetchContent")
    include(FetchContent)
    FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG v3.1.1
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(cxxopts)
endif()
list(APPEND ICARION_CORE_DEPS cxxopts::cxxopts)

# spdlog (Logging library)
find_package(spdlog CONFIG QUIET)
if(NOT spdlog_FOUND)
    message(STATUS "spdlog not found, using FetchContent")
    include(FetchContent)
    FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.12.0
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(spdlog)
endif()

# GCC 13 can emit a known false-positive -Warray-bounds warning in bundled fmt
# used by spdlog (format.h bigint path). Keep warnings enabled globally, but
# suppress this one only for the third-party spdlog build target.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
    if(TARGET spdlog)
        target_compile_options(spdlog PRIVATE -Wno-array-bounds)
    endif()
endif()

list(APPEND ICARION_CORE_DEPS spdlog::spdlog)

# OpenSSL (for SHA256 file hashing)
find_package(OpenSSL REQUIRED)
list(APPEND ICARION_CORE_DEPS OpenSSL::Crypto)

# OpenMP (optional but highly recommended for performance) - minimum version 4.5
find_package(OpenMP 4.5)
if(OpenMP_CXX_FOUND)
    list(APPEND ICARION_CORE_DEPS OpenMP::OpenMP_CXX)
    add_compile_definitions(HAVE_OPENMP)
else()
    message(WARNING "OpenMP 4.5+ not found - simulations will be single-threaded")
endif()
