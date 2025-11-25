# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Fast instrument ctests - basic instantiation and sanity checks
# Target: <100ms total runtime
add_executable(test_instrument_basic
    test_instrument_basic.cpp
)

target_link_libraries(test_instrument_basic
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_instrument_basic
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_instrument_basic PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME InstrumentBasic COMMAND test_instrument_basic)
set_tests_properties(InstrumentBasic PROPERTIES
    TIMEOUT 10
    LABELS "fast;instrument"
)
