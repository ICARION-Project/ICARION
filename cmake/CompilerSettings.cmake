if (MSVC)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra)

    # GCC 13 false-positive warnings from spdlog/fmt (system + bundled headers).
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND
       CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 13)
        add_compile_options(-Wno-array-bounds -Wno-stringop-overflow)
    endif()
endif()
