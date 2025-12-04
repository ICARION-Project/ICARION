# SPDX-License-Identifier: MIT
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# Reaction handler tests

# Test: Stochastic reaction handler
add_executable(test_stochastic_reaction_handler test_stochastic_reaction_handler.cpp)
target_link_libraries(test_stochastic_reaction_handler PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_stochastic_reaction_handler PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME StochasticReactionHandler COMMAND test_stochastic_reaction_handler)

# Test: Reaction handler factory
add_executable(test_reaction_factory test_reaction_factory.cpp)
target_link_libraries(test_reaction_factory PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_reaction_factory PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME ReactionFactory COMMAND test_reaction_factory)

# Test: Multi-gas mixture reaction weighting
add_executable(test_multi_gas_reaction test_multi_gas_reaction.cpp)
target_link_libraries(test_multi_gas_reaction PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_multi_gas_reaction PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME MultiGasReaction COMMAND test_multi_gas_reaction)

# AoS vs SoA parity for reactions
add_executable(test_reaction_soa_parity test_reaction_soa_parity.cpp)
target_link_libraries(test_reaction_soa_parity PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
target_include_directories(test_reaction_soa_parity PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)
add_test(NAME ReactionSoAParity COMMAND test_reaction_soa_parity)
