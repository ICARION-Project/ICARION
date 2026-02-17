# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Space charge validation tests

# Test executable for Poisson solver
add_executable(test_poisson_solver
    test_poisson_solver.cpp
)

target_include_directories(test_poisson_solver PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_poisson_solver PRIVATE
    Catch2::Catch2WithMain
    icarion_core
)

# Note: Space charge sources are compiled as part of icarion_core
# (via GLOB_RECURSE in src/CMakeLists.txt)

# Register with CTest
add_test(NAME PoissonSolver COMMAND test_poisson_solver)
set_tests_properties(PoissonSolver PROPERTIES
    LABELS "physics;spacecharge;poisson"
    TIMEOUT 60
)

# Test executable for charge deposition
add_executable(test_charge_deposition
    test_charge_deposition.cpp
)

target_include_directories(test_charge_deposition PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_charge_deposition PRIVATE
    Catch2::Catch2WithMain
    icarion_core
)

add_test(NAME ChargeDeposition COMMAND test_charge_deposition)
set_tests_properties(ChargeDeposition PROPERTIES
    LABELS "physics;spacecharge;deposition"
    TIMEOUT 30
)

# Test executable for space charge integration (CIC + Poisson + Forces)
add_executable(test_space_charge_integration
    test_space_charge_integration.cpp
)

target_include_directories(test_space_charge_integration PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_space_charge_integration PRIVATE
    Catch2::Catch2WithMain
    icarion_core
)

add_test(NAME SpaceChargeIntegration COMMAND test_space_charge_integration)
set_tests_properties(SpaceChargeIntegration PROPERTIES
    LABELS "physics;spacecharge;integration"
    TIMEOUT 60
)

# Direct model tests
add_executable(test_space_charge_model_direct
    test_space_charge_model_direct.cpp
)

target_include_directories(test_space_charge_model_direct PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_space_charge_model_direct PRIVATE
    Catch2::Catch2WithMain
    icarion_core
)

add_test(NAME SpaceChargeDirectModel COMMAND test_space_charge_model_direct)
set_tests_properties(SpaceChargeDirectModel PROPERTIES
    LABELS "physics;spacecharge;direct"
    TIMEOUT 10
)

# Grid model parity vs direct
add_executable(test_space_charge_model_parity
    test_space_charge_model_parity.cpp
)

target_include_directories(test_space_charge_model_parity PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_space_charge_model_parity PRIVATE
    Catch2::Catch2WithMain
    icarion_core
)

add_test(NAME SpaceChargeGridModelParity COMMAND test_space_charge_model_parity)
set_tests_properties(SpaceChargeGridModelParity PROPERTIES
    LABELS "physics;spacecharge;grid"
    TIMEOUT 15
)

# GPU model availability / stub test
add_executable(test_space_charge_gpu_model
    test_space_charge_gpu_model.cpp
)

target_include_directories(test_space_charge_gpu_model PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(test_space_charge_gpu_model PRIVATE
    Catch2::Catch2WithMain
    icarion_core
)

add_test(NAME SpaceChargeGPUModelTest COMMAND test_space_charge_gpu_model)
set_tests_properties(SpaceChargeGPUModelTest PROPERTIES
    LABELS "physics;spacecharge;gpu"
    TIMEOUT 10
)
