# =====================================================================
# ICARION - CUDA Configuration
# =====================================================================

message(STATUS "Loading CUDAConfig.cmake")

if (USE_GPU_ACCEL)
    find_package(CUDA REQUIRED)

    enable_language(CUDA)

    # Good defaults for ICARION kernels
    set(CMAKE_CUDA_STANDARD 17)

    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        # Default to common recent arch
        set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 90)
    endif()

    add_compile_definitions(ICARION_USE_GPU)

    # Detect NCCL for optional multi-GPU support
    set(_ICARION_HAVE_NCCL FALSE)
    find_path(_ICARION_NCCL_INCLUDE_DIR
              NAMES nccl.h
              HINTS $ENV{CUDA_PATH} $ENV{CUDA_HOME} ${CUDA_TOOLKIT_ROOT_DIR}
              PATH_SUFFIXES include include/nccl)
    find_library(_ICARION_NCCL_LIBRARY
                 NAMES nccl
                 HINTS $ENV{CUDA_PATH} $ENV{CUDA_HOME} ${CUDA_TOOLKIT_ROOT_DIR}
                 PATH_SUFFIXES lib lib64 lib/x64 lib/stubs lib64/stubs)

    if (_ICARION_NCCL_INCLUDE_DIR AND _ICARION_NCCL_LIBRARY)
        set(_ICARION_HAVE_NCCL TRUE)
        message(STATUS "NCCL detected: ${_ICARION_NCCL_LIBRARY}")
    else()
        message(STATUS "NCCL not found; multi-GPU support will be disabled.")
        if (NOT _ICARION_NCCL_INCLUDE_DIR)
            message(STATUS "  Searched include paths: $ENV{CUDA_PATH}, ${CUDA_TOOLKIT_ROOT_DIR}")
        endif()
    endif()

    set(ICARION_HAVE_NCCL ${_ICARION_HAVE_NCCL} CACHE INTERNAL "ICARION has NCCL available")
    if (ICARION_HAVE_NCCL)
        set(ICARION_NCCL_INCLUDE_DIR ${_ICARION_NCCL_INCLUDE_DIR} CACHE INTERNAL "ICARION NCCL include dir")
        set(ICARION_NCCL_LIBRARY ${_ICARION_NCCL_LIBRARY} CACHE INTERNAL "ICARION NCCL library")
    else()
        unset(ICARION_NCCL_INCLUDE_DIR CACHE)
        unset(ICARION_NCCL_LIBRARY CACHE)
    endif()

    # Find cuFFT library
    find_library(CUFFT_LIBRARY
                 NAMES cufft
                 HINTS $ENV{CUDA_PATH} $ENV{CUDA_HOME} ${CUDA_TOOLKIT_ROOT_DIR}
                 PATH_SUFFIXES lib lib64 lib/x64)
    
    if (CUFFT_LIBRARY)
        message(STATUS "cuFFT library found: ${CUFFT_LIBRARY}")
        set(ICARION_CUFFT_LIBRARY ${CUFFT_LIBRARY} CACHE INTERNAL "ICARION cuFFT library")
    else()
        message(FATAL_ERROR "cuFFT library not found (required for GPU Space Charge)")
    endif()
    
    message(STATUS "CUDA enabled. Architectures: ${CMAKE_CUDA_ARCHITECTURES}")
else()
    message(STATUS "CUDA disabled.")
    set(ICARION_HAVE_NCCL FALSE CACHE INTERNAL "ICARION has NCCL available")
    unset(ICARION_NCCL_INCLUDE_DIR CACHE)
    unset(ICARION_NCCL_LIBRARY CACHE)
endif()
