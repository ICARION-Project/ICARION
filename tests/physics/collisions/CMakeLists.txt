# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Collision handler tests

# Test: CollisionGeometry (core module - Phase 1)
add_executable(test_collision_geometry core/test_collision_geometry.cpp)
target_link_libraries(test_collision_geometry PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_collision_geometry PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME CollisionGeometry COMMAND test_collision_geometry)

# Test: VelocitySampling (core module - Phase 2)
add_executable(test_velocity_sampling core/test_velocity_sampling.cpp)
target_link_libraries(test_velocity_sampling PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_velocity_sampling PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME VelocitySampling COMMAND test_velocity_sampling)

# Test: CollisionKernels (core module - Phase 3)
add_executable(test_collision_kernels core/test_collision_kernels.cpp)
target_link_libraries(test_collision_kernels PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_collision_kernels PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME CollisionKernels COMMAND test_collision_kernels)

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

# Collision diagnostics (GPU only)
if(ICARION_USE_GPU)
    add_executable(test_collision_diagnostics test_collision_diagnostics.cpp)
    target_link_libraries(test_collision_diagnostics PRIVATE icarion_core)
    target_include_directories(test_collision_diagnostics PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
    )
endif()

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

# Test: Energy conservation in collision models
add_executable(test_collision_energy_conservation test_collision_energy_conservation.cpp)
target_link_libraries(test_collision_energy_conservation PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_collision_energy_conservation PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME CollisionEnergyConservation COMMAND test_collision_energy_conservation)

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

# Diagnostic: Thermalization convergence (time evolution)
add_executable(test_thermalization_convergence test_thermalization_convergence.cpp)
target_link_libraries(test_thermalization_convergence PRIVATE
    icarion_core
)
target_include_directories(test_thermalization_convergence PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_thermalization_convergence PRIVATE OpenMP::OpenMP_CXX)
endif()

# Test: Temperature scaling validation
add_executable(test_temperature_scaling test_temperature_scaling.cpp)
target_link_libraries(test_temperature_scaling PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_temperature_scaling PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_temperature_scaling PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME TemperatureScaling COMMAND test_temperature_scaling)

# Diagnostic: Realistic energy transfer with sampled neutrals
add_executable(test_realistic_transfer test_realistic_transfer.cpp)
target_link_libraries(test_realistic_transfer PRIVATE icarion_core)
target_include_directories(test_realistic_transfer PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

# Multi-gas collision handler
add_executable(test_multi_gas_collision test_multi_gas_collision.cpp)
target_link_libraries(test_multi_gas_collision PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_multi_gas_collision PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME MultiGasCollision COMMAND test_multi_gas_collision)

# GPU collision parity test (CPU vs GPU collision results)
if(ICARION_USE_GPU)
    add_executable(test_gpu_collision_parity test_gpu_collision_parity.cpp)
    target_link_libraries(test_gpu_collision_parity PRIVATE
        icarion_core
        Catch2::Catch2WithMain
    )
    target_include_directories(test_gpu_collision_parity PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
    )
    add_test(NAME GPUCollisionParity COMMAND test_gpu_collision_parity)
    
    # Quick standalone physics test
    add_executable(test_hss_physics_parity test_hss_physics_parity.cpp)
    target_link_libraries(test_hss_physics_parity PRIVATE icarion_core)
    target_include_directories(test_hss_physics_parity PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
    )
    
    # GPU thermalization convergence test
    add_executable(test_gpu_thermalization test_gpu_thermalization.cpp)
    target_link_libraries(test_gpu_thermalization PRIVATE icarion_core)
    target_include_directories(test_gpu_thermalization PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
    )
endif()

if(ICARION_USE_GPU)
    add_executable(test_gpu_only_therm test_gpu_only_therm.cpp)
    target_link_libraries(test_gpu_only_therm PRIVATE icarion_core)
    target_include_directories(test_gpu_only_therm PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
    )
endif()

