# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Physics Forces Tests

# Test: ForceRegistry
add_executable(test_force_registry
    test_force_registry.cpp
)

target_link_libraries(test_force_registry
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_force_registry
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME ForceRegistry COMMAND test_force_registry)
set_tests_properties(ForceRegistry PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: ElectricFieldForce
add_executable(test_electric_field_force
    test_electric_field_force.cpp
)

target_link_libraries(test_electric_field_force
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_electric_field_force
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME ElectricFieldForce COMMAND test_electric_field_force)
set_tests_properties(ElectricFieldForce PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: MagneticFieldForce + DampingForce
add_executable(test_magnetic_damping_forces
    test_magnetic_damping_forces.cpp
)

target_link_libraries(test_magnetic_damping_forces
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_magnetic_damping_forces
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME MagneticDampingForces COMMAND test_magnetic_damping_forces)
set_tests_properties(MagneticDampingForces PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: SpaceChargeForce
add_executable(test_space_charge_force
    test_space_charge_force.cpp
)

target_link_libraries(test_space_charge_force
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_space_charge_force
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME SpaceChargeForce COMMAND test_space_charge_force)
set_tests_properties(SpaceChargeForce PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: Force Integration (all forces combined)
add_executable(test_force_integration
    test_force_integration.cpp
)

target_link_libraries(test_force_integration
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_force_integration
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME ForceIntegration COMMAND test_force_integration)
set_tests_properties(ForceIntegration PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
