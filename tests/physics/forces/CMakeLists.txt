# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: 2026 Christoph Schaefer

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

# Test: AnalyticalFieldModel parity with ElectricFieldForce
add_executable(test_field_model_parity
    test_field_model_parity.cpp
)

target_link_libraries(test_field_model_parity
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_field_model_parity
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME FieldModelParity COMMAND test_field_model_parity)
set_tests_properties(FieldModelParity PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: FieldProviderModel integration path
add_executable(test_field_model_provider
    test_field_model_provider.cpp
)

target_link_libraries(test_field_model_provider
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_field_model_provider
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME FieldModelProvider COMMAND test_field_model_provider)
set_tests_properties(FieldModelProvider PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: MagneticFieldForce + DampingForce
add_executable(test_magnetic_and_damping_forces
    test_magnetic_damping_forces.cpp
)

target_link_libraries(test_magnetic_and_damping_forces
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_magnetic_and_damping_forces
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME MagneticAndDampingForces COMMAND test_magnetic_and_damping_forces)
set_tests_properties(MagneticAndDampingForces PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Test: SpaceChargeForce
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

# Test: Gas Flow Transport (physics validation)
add_executable(test_gas_flow_transport
    test_gas_flow_transport.cpp
)

target_link_libraries(test_gas_flow_transport
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_gas_flow_transport
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/tests
)

add_test(NAME GasFlowTransport COMMAND test_gas_flow_transport)
set_tests_properties(GasFlowTransport PROPERTIES
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
