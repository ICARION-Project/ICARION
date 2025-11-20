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
