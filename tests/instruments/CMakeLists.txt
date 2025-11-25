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

# IMS drift physics validation
add_executable(test_ims_drift
    test_ims_drift.cpp
)

target_link_libraries(test_ims_drift
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_ims_drift
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_ims_drift PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME IMSDrift COMMAND test_ims_drift)
set_tests_properties(IMSDrift PROPERTIES
    TIMEOUT 30
    LABELS "instrument;ims;physics"
)

# TOF flight time validation
add_executable(test_tof_flight_time
    test_tof_flight_time.cpp
)

target_link_libraries(test_tof_flight_time
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_tof_flight_time
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_tof_flight_time PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME TOFFlightTime COMMAND test_tof_flight_time)
set_tests_properties(TOFFlightTime PROPERTIES
    TIMEOUT 30
    LABELS "instrument;tof;physics"
)

# LQIT Mathieu stability validation
add_executable(test_lqit_stability
    test_lqit_stability.cpp
)

target_link_libraries(test_lqit_stability
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_lqit_stability
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_lqit_stability PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME LQITStability COMMAND test_lqit_stability)
set_tests_properties(LQITStability PROPERTIES
    TIMEOUT 60
    LABELS "instrument;lqit;physics"
)

# Orbitrap axial frequency validation
add_executable(test_orbitrap_frequency
    test_orbitrap_frequency.cpp
)

target_link_libraries(test_orbitrap_frequency
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_orbitrap_frequency
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_orbitrap_frequency PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME OrbitrapFrequency COMMAND test_orbitrap_frequency)
set_tests_properties(OrbitrapFrequency PROPERTIES
    TIMEOUT 30
    LABELS "instrument;orbitrap;physics"
)

# Domain transition validation
add_executable(test_domain_transition
    test_domain_transition.cpp
)

target_link_libraries(test_domain_transition
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_domain_transition
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_domain_transition PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME DomainTransition COMMAND test_domain_transition)
set_tests_properties(DomainTransition PROPERTIES
    TIMEOUT 30
    LABELS "instrument;transition;physics"
)

# Orbitrap frequency validation
add_executable(test_orbitrap_frequency
    test_orbitrap_frequency.cpp
)

target_link_libraries(test_orbitrap_frequency
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_orbitrap_frequency
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

set_target_properties(test_orbitrap_frequency PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME OrbitrapFrequency COMMAND test_orbitrap_frequency)
set_tests_properties(OrbitrapFrequency PROPERTIES
    TIMEOUT 30
    LABELS "instrument;orbitrap;physics"
)
