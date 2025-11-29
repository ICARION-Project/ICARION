# =====================================================================
# ICARION - Target Configuration Helpers
# =====================================================================

message(STATUS "Loading Targets.cmake")

# Function to add sanitizers
function(add_sanitizers target)
    if (ICARION_ENABLE_SANITIZERS AND CMAKE_BUILD_TYPE STREQUAL "Debug")
        if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
            target_compile_options(${target} PUBLIC -fsanitize=address,undefined)
            target_link_options(${target} PUBLIC -fsanitize=address,undefined)
        endif()
    endif()
endfunction()
