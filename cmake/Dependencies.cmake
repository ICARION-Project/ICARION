# =====================================================================
# ICARION - Dependencies
# =====================================================================

message(STATUS "Loading Dependencies.cmake")

# JSON (jsoncpp)
find_package(jsoncpp CONFIG REQUIRED)
set(ICARION_CORE_DEPS jsoncpp_lib)

# HDF5
find_package(HDF5 COMPONENTS C CXX REQUIRED)
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
list(APPEND ICARION_CORE_DEPS spdlog::spdlog)

# OpenSSL (for SHA256 file hashing)
find_package(OpenSSL REQUIRED)
list(APPEND ICARION_CORE_DEPS OpenSSL::Crypto)

# OpenMP (optional but highly recommended for performance)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    list(APPEND ICARION_CORE_DEPS OpenMP::OpenMP_CXX)
    add_compile_definitions(HAVE_OPENMP)
else()
    message(WARNING "OpenMP not found - simulations will be single-threaded")
endif()

