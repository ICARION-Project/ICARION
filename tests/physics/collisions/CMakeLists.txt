# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Collision handler tests

# Test: Collision Handler Factory
add_executable(test_collision_factory test_collision_factory.cpp)
target_link_libraries(test_collision_factory PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_collision_factory PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME CollisionFactory COMMAND test_collision_factory)

# Test: EHSS Collision Handler (thermalization)
add_executable(test_ehss_collision_handler test_ehss_collision_handler.cpp)
target_link_libraries(test_ehss_collision_handler PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_ehss_collision_handler PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME EHSSCollisionHandler COMMAND test_ehss_collision_handler)

# Test: HSS Collision Handler (thermalization)
add_executable(test_hss_collision_handler test_hss_collision_handler.cpp)
target_link_libraries(test_hss_collision_handler PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_hss_collision_handler PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME HSSCollisionHandler COMMAND test_hss_collision_handler)

# Test: OU Collision Handler (thermalization)
add_executable(test_ou_collision_handler test_ou_collision_handler.cpp)
target_link_libraries(test_ou_collision_handler PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_ou_collision_handler PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME OUCollisionHandler COMMAND test_ou_collision_handler)

# Test: Geometry Utils (conversion, loading)
add_executable(test_geometry_utils test_geometry_utils.cpp)
target_link_libraries(test_geometry_utils PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_geometry_utils PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME GeometryUtils COMMAND test_geometry_utils)

# Diagnostic: Energy distribution analysis
add_executable(test_collision_energy_diagnostic test_collision_energy_diagnostic.cpp)
target_link_libraries(test_collision_energy_diagnostic PRIVATE
    icarion_core
)
target_include_directories(test_collision_energy_diagnostic PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_collision_energy_diagnostic PRIVATE OpenMP::OpenMP_CXX)
endif()
