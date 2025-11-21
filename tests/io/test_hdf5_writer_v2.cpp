// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @brief Smoke test for HDF5 C++ library
 * 
 * This is a minimal test to verify HDF5 C++ library is available and functional.
 * Full integration testing of HDF5Writer v2 requires field name alignment between
 * FullConfig types and IonState structures.
 * 
 * @todo Add full HDF5Writer v2 integration tests once resolved:
 * - IonState: domain_index → current_domain_index
 * - IonState: charge_C → ion_charge_C  
 * - SpeciesProperties: reduced_mobility_m2Vs → mobility conversion
 * - SpeciesProperties: ccs_m2 → CCS_m2
 * - Reaction: rate_constant → rate_constant_m3s
 * - DomainConfig: solver string → SolverType enum conversion
 */

#include <catch2/catch_test_macros.hpp>
#include <H5Cpp.h>
#include <filesystem>

TEST_CASE("HDF5 library is available and functional", "[hdf5][io]") {
    std::string test_file = "/tmp/test_hdf5_smoke.h5";
    
    if (std::filesystem::exists(test_file)) {
        std::filesystem::remove(test_file);
    }
    
    // Create file with basic structure
    {
        H5::H5File file(test_file, H5F_ACC_TRUNC);
        H5::Group group = file.createGroup("/test_group");
        
        // Write integer dataset
        hsize_t dims[1] = {5};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = group.createDataSet("test_data", 
                                                   H5::PredType::NATIVE_INT, 
                                                   dataspace);
        int data[5] = {1, 2, 3, 4, 5};
        dataset.write(data, H5::PredType::NATIVE_INT);
        
        dataset.close();
        group.close();
        file.close();
    }
    
    REQUIRE(std::filesystem::exists(test_file));
    
    // Verify read-back
    {
        H5::H5File file(test_file, H5F_ACC_RDONLY);
        REQUIRE(file.nameExists("/test_group"));
        REQUIRE(file.nameExists("/test_group/test_data"));
        
        H5::Group group = file.openGroup("/test_group");
        H5::DataSet dataset = group.openDataSet("test_data");
        
        int read_data[5];
        dataset.read(read_data, H5::PredType::NATIVE_INT);
        
        REQUIRE(read_data[0] == 1);
        REQUIRE(read_data[2] == 3);
        REQUIRE(read_data[4] == 5);
        
        dataset.close();
        group.close();
        file.close();
    }
    
    // Cleanup
    std::filesystem::remove(test_file);
}
