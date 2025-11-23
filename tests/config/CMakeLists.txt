# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Config system tests
add_executable(test_config_loader
    test_config_loader.cpp
)

target_link_libraries(test_config_loader
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

target_include_directories(test_config_loader
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME ConfigLoader COMMAND test_config_loader)

# Species and Reaction loader smoke tests
add_executable(test_species_loader_smoke
    test_species_loader.cpp
)

target_link_libraries(test_species_loader_smoke
    PRIVATE
    icarion_core
)

target_include_directories(test_species_loader_smoke
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME SpeciesLoaderSmoke COMMAND test_species_loader_smoke)
set_tests_properties(SpeciesLoaderSmoke PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# Species loader unit tests
add_executable(test_species_loader_unit
    test_species_loader_unit.cpp
)

target_link_libraries(test_species_loader_unit
    PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)

target_include_directories(test_species_loader_unit
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME SpeciesLoaderUnit COMMAND test_species_loader_unit)

# Reaction loader unit tests
add_executable(test_reaction_loader_unit
    test_reaction_loader_unit.cpp
)

target_link_libraries(test_reaction_loader_unit
    PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)

target_include_directories(test_reaction_loader_unit
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME ReactionLoaderUnit COMMAND test_reaction_loader_unit)

# Reaction validation tests (order term validation rules #1-5)
add_executable(test_reaction_validation
    test_reaction_validation.cpp
)

target_link_libraries(test_reaction_validation
    PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)

target_include_directories(test_reaction_validation
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME ReactionValidation COMMAND test_reaction_validation)

# Database integration test (ConfigLoader with species/reactions)
add_executable(test_database_integration
    test_database_integration.cpp
)

target_link_libraries(test_database_integration
    PRIVATE
    icarion_core
)

target_include_directories(test_database_integration
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

add_test(NAME DatabaseIntegration COMMAND test_database_integration)

# Ion loader tests
add_executable(test_ion_loader
    test_ion_loader.cpp
)

target_link_libraries(test_ion_loader
    PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)

target_include_directories(test_ion_loader
    PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

add_test(NAME IonLoader COMMAND test_ion_loader)
