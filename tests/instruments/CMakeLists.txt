# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

add_executable(test_instrument_physics
    ims_physics_test.cpp
    orbitrap_physics_test.cpp
)

target_link_libraries(test_instrument_physics
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_instrument_physics
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

add_test(NAME InstrumentPhysics COMMAND test_instrument_physics)
