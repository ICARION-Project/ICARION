# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# HDF5Writer test
add_executable(test_hdf5_writer
    test_hdf5_writer.cpp
)

target_link_libraries(test_hdf5_writer
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
        ${HDF5_LIBRARIES}
)

target_include_directories(test_hdf5_writer
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src
        ${HDF5_INCLUDE_DIRS}
)

# Add to test suite
add_test(NAME HDF5Writer COMMAND test_hdf5_writer)
set_tests_properties(HDF5Writer PROPERTIES
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)
