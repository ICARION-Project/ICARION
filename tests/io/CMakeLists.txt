# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

# HDF5Writer v2 test
add_executable(test_hdf5_writer_v2
    test_hdf5_writer_v2.cpp
)

target_link_libraries(test_hdf5_writer_v2
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
        ${HDF5_LIBRARIES}
)

target_include_directories(test_hdf5_writer_v2
    PRIVATE
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_SOURCE_DIR}/src
        ${HDF5_INCLUDE_DIRS}
)

# Add to test suite
add_test(NAME HDF5WriterV2 COMMAND test_hdf5_writer_v2)
