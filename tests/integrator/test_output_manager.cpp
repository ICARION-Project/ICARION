// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file test_output_manager.cpp
 * @brief Unit tests for OutputManager (HDF5 + text logging)
 * 
 * Tests:
 * 1. Basic initialization (metadata, species)
 * 2. Buffer management (accumulation, flush)
 * 3. Time-based write triggers
 * 4. Text logging (progress, finalization)
 * 5. Error handling (uninitialized, empty filename)
 */

#include "core/integrator/OutputManager.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "H5Cpp.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <filesystem>
#include <fstream>

using namespace ICARION::integrator;
using namespace ICARION::config;
using namespace ICARION;

namespace {
    // Test data directory
    const std::filesystem::path test_dir = "/tmp/icarion_test_output_manager";
    
    // Helper: Create minimal test config
    FullConfig create_test_config() {
        FullConfig config;
        config.simulation.total_time_s = 1e-6;
        config.simulation.dt_s = 1e-9;
        config.output.trajectory_file = "test_output.h5";
        
        // Empty domains list (valid for basic HDF5 creation test)
        return config;
    }
    
    // Helper: Create test ions (SoA)
    ICARION::core::IonEnsemble create_test_ions(size_t count) {
        std::vector<ICARION::core::IonState> ions;
        ions.reserve(count);
        
        for (size_t i = 0; i < count; ++i) {
            core::IonState ion;
            ion.mass_kg = 100.0 * 1.66053906660e-27;
            ion.ion_charge_C = 1.602176634e-19;
            ion.pos = {0.0, 0.0, static_cast<double>(i) * 0.001};
            ion.vel = {0.0, 0.0, 100.0};
            ion.active = true;
            ion.born = true;  // Must be born to be counted as active
            ion.species_id = "TestIon";
            ions.push_back(ion);
        }
        
        return core::IonEnsemble::from_legacy(ions);
    }
    
    // Helper: Read HDF5 attribute
    template<typename T>
    T read_hdf5_attribute(const std::string& filename, const std::string& attr_name, 
                          const H5::PredType& dtype) {
        H5::H5File file(filename, H5F_ACC_RDONLY);
        H5::Attribute attr = file.openAttribute(attr_name);
        T value;
        attr.read(dtype, &value);
        return value;
    }
}

TEST_CASE("OutputManager - Basic initialization", "[OutputManager]") {
    // Setup
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_init.h5";
    auto log_file = test_dir / "test_init.log";
    
    OutputManager manager(hdf5_file.string(), log_file.string());
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    // Initialize
    REQUIRE_NOTHROW(manager.initialize(config, ions));
    
    // Verify HDF5 file created
    REQUIRE(std::filesystem::exists(hdf5_file));
    
    // Verify text log created
    REQUIRE(std::filesystem::exists(log_file));
    
    // Cleanup
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - Buffer accumulation", "[OutputManager]") {
    // Setup
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_buffer.h5";
    
    OutputManager manager(hdf5_file.string(), "", 1e-3, 5);  // Buffer max = 5
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    manager.initialize(config, ions);
    
    // Log 3 steps (below buffer limit)
    manager.log_step(0.0, ions);
    manager.log_step(1e-6, ions);
    manager.log_step(2e-6, ions);
    
    REQUIRE(manager.buffer_size() == 3);
    
    // Log 2 more steps (reaches buffer limit → triggers flush)
    manager.log_step(3e-6, ions);
    manager.log_step(4e-6, ions);
    
    REQUIRE(manager.buffer_size() == 5);
    
    // Next log_step triggers auto-flush
    manager.log_step(5e-6, ions);
    
    REQUIRE(manager.buffer_size() == 1);  // Buffer cleared, new snapshot added
    
    // Cleanup
    manager.finalize(5e-6, ions);
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - Time-based write trigger", "[OutputManager]") {
    // Setup
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_time_trigger.h5";
    
    double write_interval = 1e-4;  // 100 μs
    OutputManager manager(hdf5_file.string(), "", write_interval, 100);  // Large buffer
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    manager.initialize(config, ions);
    
    // Log steps before time interval
    manager.log_step(0.0, ions);
    manager.log_step(5e-5, ions);  // 50 μs
    
    REQUIRE(manager.buffer_size() == 2);
    REQUIRE_FALSE(manager.should_write(5e-5));
    
    // Check that time interval would trigger write
    REQUIRE(manager.should_write(1.5e-4));  // 150 μs exceeds 100 μs interval
    
    // Log step after time interval (triggers auto-flush)
    manager.log_step(1.5e-4, ions);
    
    // Auto-flush triggered - buffer cleared, new snapshot added
    REQUIRE(manager.buffer_size() == 1);
    
    // Cleanup
    manager.finalize(1.5e-4, ions);
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - Text logging", "[OutputManager]") {
    // Setup
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_text_log.h5";
    auto log_file = test_dir / "test_text_log.log";
    
    OutputManager manager(hdf5_file.string(), log_file.string());
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    manager.initialize(config, ions);
    
    // Log progress
    manager.log_progress("Simulation started");
    manager.log_progress("50% completed (t = 0.5 ms)");
    
    // Finalize (writes completion summary)
    manager.finalize(1e-3, ions);
    
    // Verify log file contains progress messages
    std::ifstream log_stream(log_file);
    std::string log_content((std::istreambuf_iterator<char>(log_stream)),
                            std::istreambuf_iterator<char>());
    
    REQUIRE(log_content.find("Simulation started") != std::string::npos);
    REQUIRE(log_content.find("50% completed") != std::string::npos);
    
    // Cleanup
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - Finalization metadata", "[OutputManager]") {
    // Setup
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_finalize.h5";
    
    OutputManager manager(hdf5_file.string());
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    manager.initialize(config, ions);
    manager.log_step(0.0, ions);
    
    // Deactivate some ions
    auto* active = ions.active_data();
    active[3] = 0;
    active[7] = 0;
    
    manager.finalize(1e-3, ions);
    
    // Verify HDF5 completion metadata (now under /metadata/completion/ as datasets)
    H5::H5File file(hdf5_file.string(), H5F_ACC_RDONLY);
    H5::Group metadata_grp = file.openGroup("/metadata");
    H5::Group completion_grp = metadata_grp.openGroup("completion");
    
    H5::DataSet success_ds = completion_grp.openDataSet("success");
    int success_val;
    success_ds.read(&success_val, H5::PredType::NATIVE_INT);
    REQUIRE(success_val == 1);
    
    H5::DataSet active_ds = completion_grp.openDataSet("active_ions");
    int active_count;
    active_ds.read(&active_count, H5::PredType::NATIVE_INT);
    REQUIRE(active_count == 8);  // 10 - 2 deactivated
    
    completion_grp.close();
    metadata_grp.close();
    file.close();
    
    // Cleanup
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - Manual flush", "[OutputManager]") {
    // Setup
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_manual_flush.h5";
    
    OutputManager manager(hdf5_file.string(), "", 1.0, 100);  // Long interval, large buffer
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    manager.initialize(config, ions);
    
    // Log steps (won't auto-flush)
    manager.log_step(0.0, ions);
    manager.log_step(1e-4, ions);
    manager.log_step(2e-4, ions);
    
    REQUIRE(manager.buffer_size() == 3);
    
    // Manual flush
    manager.flush();
    
    REQUIRE(manager.buffer_size() == 0);
    
    // Cleanup
    manager.finalize(2e-4, ions);
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - Error handling: Empty filename", "[OutputManager]") {
    REQUIRE_THROWS_AS(
        OutputManager("", ""),
        std::invalid_argument
    );
}

TEST_CASE("OutputManager - Error handling: Uninitialized log_step", "[OutputManager]") {
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_uninit.h5";
    
    OutputManager manager(hdf5_file.string());
    auto ions = create_test_ions(10);
    
    // Try to log without initialization
    REQUIRE_THROWS_AS(
        manager.log_step(0.0, ions),
        std::runtime_error
    );
    
    std::filesystem::remove_all(test_dir);
}

TEST_CASE("OutputManager - No text log (disabled)", "[OutputManager]") {
    std::filesystem::create_directories(test_dir);
    auto hdf5_file = test_dir / "test_no_text_log.h5";
    
    OutputManager manager(hdf5_file.string());  // No log filename
    auto config = create_test_config();
    auto ions = create_test_ions(10);
    
    manager.initialize(config, ions);
    
    // Text logging disabled
    REQUIRE_FALSE(manager.has_text_log());
    REQUIRE(manager.get_log_filename().empty());
    
    // log_progress does nothing (no crash)
    REQUIRE_NOTHROW(manager.log_progress("This won't be written"));
    
    manager.finalize(1e-3, ions);
    std::filesystem::remove_all(test_dir);
}
