// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @brief Integration tests for HDF5Writer v2
 * 
 * Tests the complete HDF5 v2.0 output format with FullConfig.
 * Verifies file structure, metadata writing, trajectory appending, and finalization.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/io/hdf5Writer.h"
#include "core/config/types/FullConfig.h"
#include "core/config/loader/ConfigLoader.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include <H5Cpp.h>
#include <filesystem>
#include <fstream>

using namespace ICARION;
using Catch::Approx;

namespace {
    std::string get_test_file() {
        static std::string file = "/tmp/test_hdf5_v2.h5";
        if (std::filesystem::exists(file)) {
            std::filesystem::remove(file);
        }
        return file;
    }
    
    config::FullConfig create_minimal_config() {
        config::FullConfig cfg;
        
        // Simulation
        cfg.simulation.dt_s = 1e-9;
        cfg.simulation.total_time_s = 1e-6;
        cfg.simulation.total_steps = 1000;
        cfg.simulation.write_interval = 100;
        cfg.simulation.rng_seed = 42;
        cfg.simulation.integrator = "RK4";
        cfg.simulation.enable_gpu = false;
        
        // Physics
        cfg.physics.collision_model = config::CollisionModel::EHSS;
        cfg.physics.enable_reactions = false;
        cfg.physics.enable_space_charge = false;
        
        // Output
        cfg.output.trajectory_file = "test.h5";
        cfg.output.folder = "./output";
        
        // Species
        config::SpeciesProperties h3o;
        h3o.id = "H3O+";
        h3o.mass_amu = 19.02;
        h3o.charge = 1;
        h3o.mobility_cm2Vs = 2.8e-4;
        h3o.CCS_A2 = 110.0;
        h3o.convert_to_SI();  // Convert to SI units
        cfg.species_db.species["H3O+"] = h3o;
        
        config::SpeciesProperties no2;
        no2.id = "NO2+";
        no2.mass_amu = 46.0;
        no2.charge = 1;
        no2.mobility_cm2Vs = 2.5e-4;
        no2.CCS_A2 = 120.0;
        no2.convert_to_SI();
        cfg.species_db.species["NO2+"] = no2;
        
        // Domain
        config::DomainConfig domain;
        domain.name = "test_domain";
        domain.instrument = config::Instrument::IMS;
        domain.solver = config::SolverType::RK4;
        domain.domain_index = 0;
        domain.geometry.length_m = 0.1;
        domain.geometry.radius_m = 0.01;
        domain.environment.pressure_Pa = 200.0;
        domain.environment.temperature_K = 300.0;
        domain.environment.gas_species = "N2";
        cfg.domains.push_back(domain);
        
        return cfg;
    }
    
    std::vector<core::IonState> create_test_ions() {
        std::vector<core::IonState> ions;
        
        for (int i = 0; i < 10; ++i) {
            core::IonState ion;
            ion.species_id = (i % 2 == 0) ? "H3O+" : "NO2+";
            ion.pos = core::Vec3(0.001 * i, 0.0001 * i, 0);
            ion.vel = core::Vec3(100, 10, 0);
            ion.mass_kg = (i % 2 == 0) ? 19.02 * 1.66e-27 : 46.0 * 1.66e-27;
            ion.ion_charge_C = 1.602e-19;
            ion.birth_time_s = 0.0;
            ion.current_domain_index = 0;
            ion.active = true;
            ions.push_back(ion);
        }
        
        return ions;
    }

    core::IonEnsemble to_ensemble(const std::vector<core::IonState>& ions) {
        return core::IonEnsemble::from_legacy(ions);
    }
}

TEST_CASE("HDF5 library is available and functional", "[hdf5][io][smoke]") {
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

TEST_CASE("HDF5Writer v2 creates correct file structure", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    
    // Create file with metadata
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    io::HDF5Writer::create_file(test_file, config, ensemble, "test_git_hash", "gcc 11.4.0 -O3");
    
    REQUIRE(std::filesystem::exists(test_file));
    
    // Verify structure
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    
    // Check main groups
    REQUIRE(file.nameExists("/metadata"));
    REQUIRE(file.nameExists("/metadata/config"));
    REQUIRE(file.nameExists("/metadata/reproducibility"));
    REQUIRE(file.nameExists("/metadata/system"));
    REQUIRE(file.nameExists("/metadata/species"));
    REQUIRE(file.nameExists("/trajectory"));
    REQUIRE(file.nameExists("/ions"));
    REQUIRE(file.nameExists("/domains"));
    REQUIRE(file.nameExists("/domains/domain_0"));
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes simulation metadata correctly", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "abc123def", "gcc 11.4.0");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group cfg_group = file.openGroup("/metadata/config");
    
    // Check timestep
    H5::DataSet ds_dt = cfg_group.openDataSet("dt_s");
    double dt;
    ds_dt.read(&dt, H5::PredType::NATIVE_DOUBLE);
    REQUIRE(dt == Approx(1e-9));
    
    // Check total time
    H5::DataSet ds_time = cfg_group.openDataSet("total_time_s");
    double total_time;
    ds_time.read(&total_time, H5::PredType::NATIVE_DOUBLE);
    REQUIRE(total_time == Approx(1e-6));
    
    // Check integrator
    H5::DataSet ds_integ = cfg_group.openDataSet("integrator");
    H5::StrType str_type = ds_integ.getStrType();
    std::string integrator;
    ds_integ.read(integrator, str_type);
    REQUIRE(integrator == "RK4");
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes reproducibility metadata", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "abc123def", "gcc 11.4.0");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group repro = file.openGroup("/metadata/reproducibility");
    
    // Check RNG seed
    H5::DataSet ds_seed = repro.openDataSet("global_seed");
    unsigned int seed;
    ds_seed.read(&seed, H5::PredType::NATIVE_UINT);
    REQUIRE(seed == 42);
    
    // Check git hash
    H5::DataSet ds_git = repro.openDataSet("git_hash");
    H5::StrType str_type = ds_git.getStrType();
    std::string git_hash;
    ds_git.read(git_hash, str_type);
    REQUIRE(git_hash == "abc123def");
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes species metadata in tabular format", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "test", "test");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group species_group = file.openGroup("/metadata/species");
    
    // Check names dataset
    H5::DataSet ds_names = species_group.openDataSet("names");
    H5::DataSpace space = ds_names.getSpace();
    hsize_t dims[1];
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 2);  // H3O+ and NO2+
    
    // Check mass dataset
    H5::DataSet ds_mass = species_group.openDataSet("mass_kg");
    space = ds_mass.getSpace();
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 2);
    
    std::vector<double> masses(2);
    ds_mass.read(masses.data(), H5::PredType::NATIVE_DOUBLE);
    REQUIRE(masses[0] > 0.0);
    REQUIRE(masses[1] > 0.0);
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes domain configuration", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "test", "test");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group domain = file.openGroup("/domains/domain_0");
    
    // Check domain name
    H5::DataSet ds_name = domain.openDataSet("name");
    H5::StrType str_type = ds_name.getStrType();
    std::string name;
    ds_name.read(name, str_type);
    REQUIRE(name == "test_domain");
    
    // Check geometry group exists
    REQUIRE(file.nameExists("/domains/domain_0/geometry"));
    H5::Group geom = domain.openGroup("geometry");
    
    H5::DataSet ds_length = geom.openDataSet("length_m");
    double length;
    ds_length.read(&length, H5::PredType::NATIVE_DOUBLE);
    REQUIRE(length == Approx(0.1));
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes ion initial conditions", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "test", "test");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group ion_group = file.openGroup("/ions");
    
    // Check number of ions
    H5::DataSet ds_species = ion_group.openDataSet("initial_species_id");
    H5::DataSpace space = ds_species.getSpace();
    hsize_t dims[1];
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 10);
    
    // Check position data
    H5::DataSet ds_pos_x = ion_group.openDataSet("initial_pos_x");
    space = ds_pos_x.getSpace();
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 10);
    
    std::vector<double> pos_x(10);
    ds_pos_x.read(pos_x.data(), H5::PredType::NATIVE_DOUBLE);
    REQUIRE(pos_x[0] == Approx(0.0));
    REQUIRE(pos_x[9] == Approx(0.009));
    
    file.close();
}

TEST_CASE("HDF5Writer v2 appends trajectory data correctly", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "test", "test");
    
    // Append first timestep
    for (auto& ion : ions) {
        ion.pos.x += 0.0001;
        ion.vel.x += 1.0;
    }
    io::HDF5Writer::append_trajectory(test_file, 0.0, to_ensemble(ions));
    
    // Append second timestep
    for (auto& ion : ions) {
        ion.pos.x += 0.0001;
        ion.vel.x += 1.0;
    }
    io::HDF5Writer::append_trajectory(test_file, 1e-9, to_ensemble(ions));
    
    // Append third timestep
    for (auto& ion : ions) {
        ion.pos.x += 0.0001;
        ion.vel.x += 1.0;
    }
    io::HDF5Writer::append_trajectory(test_file, 2e-9, to_ensemble(ions));
    
    // Verify trajectory data
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group traj = file.openGroup("/trajectory");
    
    // Check time array
    H5::DataSet ds_time = traj.openDataSet("time");
    H5::DataSpace space = ds_time.getSpace();
    hsize_t dims[1];
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 3);  // 3 timesteps
    
    std::vector<double> times(3);
    ds_time.read(times.data(), H5::PredType::NATIVE_DOUBLE);
    REQUIRE(times[0] == Approx(0.0));
    REQUIRE(times[1] == Approx(1e-9));
    REQUIRE(times[2] == Approx(2e-9));
    
    // Check positions array [T × N × 3]
    H5::DataSet ds_pos = traj.openDataSet("positions");
    space = ds_pos.getSpace();
    hsize_t pos_dims[3];
    space.getSimpleExtentDims(pos_dims);
    REQUIRE(pos_dims[0] == 3);   // 3 timesteps
    REQUIRE(pos_dims[1] == 10);  // 10 ions
    REQUIRE(pos_dims[2] == 3);   // x,y,z
    
    // Check velocities array [T × N × 3]
    H5::DataSet ds_vel = traj.openDataSet("velocities");
    space = ds_vel.getSpace();
    hsize_t vel_dims[3];
    space.getSimpleExtentDims(vel_dims);
    REQUIRE(vel_dims[0] == 3);
    REQUIRE(vel_dims[1] == 10);
    REQUIRE(vel_dims[2] == 3);
    
    file.close();
}

TEST_CASE("HDF5Writer v2 finalization writes completion metadata", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    
    io::HDF5Writer::create_file(test_file, config, ensemble, "test", "test");
    
    // Simulate some trajectory steps
    io::HDF5Writer::append_trajectory(test_file, 0.0, to_ensemble(ions));
    io::HDF5Writer::append_trajectory(test_file, 1e-9, to_ensemble(ions));
    
    // Finalize
    io::HDF5Writer::finalize(test_file, true, 2e-9, 8);
    
    // Verify completion metadata
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    REQUIRE(file.nameExists("/metadata/completion"));
    
    H5::Group completion = file.openGroup("/metadata/completion");
    
    // Check success flag
    H5::DataSet ds_success = completion.openDataSet("success");
    hbool_t success;
    ds_success.read(&success, H5::PredType::NATIVE_HBOOL);
    REQUIRE(success);
    
    // Check final time
    H5::DataSet ds_time = completion.openDataSet("final_time_s");
    double final_time;
    ds_time.read(&final_time, H5::PredType::NATIVE_DOUBLE);
    REQUIRE(final_time == Approx(2e-9));
    
    // Check active ions
    H5::DataSet ds_active = completion.openDataSet("active_ions");
    int active;
    ds_active.read(&active, H5::PredType::NATIVE_INT);
    REQUIRE(active == 8);
    
    file.close();
}

TEST_CASE("HDF5Writer v2 handles empty reaction database", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    
    // reactions vector is already empty by default
    REQUIRE(config.reaction_db.reactions.empty());
    
    // Should not crash when writing empty reactions
    REQUIRE_NOTHROW(io::HDF5Writer::create_file(test_file, config, to_ensemble(ions), "test", "test"));
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    // Reactions group should still exist but be empty
    // (Implementation detail - may or may not create group for empty database)
    file.close();
}

TEST_CASE("HDF5Writer v2 writes reaction metadata correctly", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    
    // Enable reactions and add test reaction
    config.physics.enable_reactions = true;
    
    config::Reaction rxn;
    rxn.id = "charge_transfer_1";
    rxn.reactant = "H3O+";
    rxn.product = "NO2+";
    rxn.rate_constant = 1.5e-15;
    config.reaction_db.reactions.push_back(rxn);
    
    io::HDF5Writer::create_file(test_file, config, to_ensemble(ions), "test", "test");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    
    // Check reactions metadata exists
    REQUIRE(file.nameExists("/metadata/reactions"));
    H5::Group rxn_group = file.openGroup("/metadata/reactions");
    
    // Check reaction ID
    H5::DataSet ds_ids = rxn_group.openDataSet("id");
    H5::DataSpace space = ds_ids.getSpace();
    hsize_t dims[1];
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 1);  // 1 reaction
    
    // Check rate constant (stored with units in the name)
    H5::DataSet ds_rate = rxn_group.openDataSet("rate_constant_m3s");
    std::vector<double> rates(1);
    ds_rate.read(rates.data(), H5::PredType::NATIVE_DOUBLE);
    REQUIRE(rates[0] == Approx(1.5e-15));
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes multiple domains correctly", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    
    // Add second domain
    config::DomainConfig domain2;
    domain2.name = "drift_region_2";
    domain2.instrument = config::Instrument::TOF;
    domain2.solver = config::SolverType::RK45;
    domain2.domain_index = 1;
    domain2.geometry.length_m = 0.2;
    domain2.geometry.radius_m = 0.015;
    domain2.environment.pressure_Pa = 100.0;
    domain2.environment.temperature_K = 320.0;
    domain2.environment.gas_species = "He";
    config.domains.push_back(domain2);
    
    // Add third domain
    config::DomainConfig domain3;
    domain3.name = "detector";
    domain3.instrument = config::Instrument::NoFixedInstrument;
    domain3.solver = config::SolverType::RK4;
    domain3.domain_index = 2;
    domain3.geometry.length_m = 0.05;
    domain3.geometry.radius_m = 0.02;
    domain3.environment.pressure_Pa = 1e-3;
    domain3.environment.temperature_K = 300.0;
    domain3.environment.gas_species = "Vacuum";
    config.domains.push_back(domain3);
    
    io::HDF5Writer::create_file(test_file, config, to_ensemble(ions), "test", "test");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    
    // Check all domain groups exist
    REQUIRE(file.nameExists("/domains/domain_0"));
    REQUIRE(file.nameExists("/domains/domain_1"));
    REQUIRE(file.nameExists("/domains/domain_2"));
    
    // Verify domain 1 details
    H5::Group dom1 = file.openGroup("/domains/domain_1");
    H5::DataSet ds_name = dom1.openDataSet("name");
    H5::StrType str_type = ds_name.getStrType();
    std::string name;
    ds_name.read(name, str_type);
    REQUIRE(name == "drift_region_2");
    
    H5::DataSet ds_instr = dom1.openDataSet("instrument");
    std::string instrument;
    ds_instr.read(instrument, str_type);
    REQUIRE(instrument == "TOF");
    
    // Verify domain 2 geometry
    H5::Group dom2 = file.openGroup("/domains/domain_2");
    H5::Group geom2 = dom2.openGroup("geometry");
    H5::DataSet ds_length = geom2.openDataSet("length_m");
    double length;
    ds_length.read(&length, H5::PredType::NATIVE_DOUBLE);
    REQUIRE(length == Approx(0.05));
    
    file.close();
}

TEST_CASE("HDF5Writer v2 writes system metadata correctly", "[hdf5][io]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    auto ions = create_test_ions();
    
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    io::HDF5Writer::create_file(test_file, config, ensemble, "abc123", "gcc 11.4.0");
    
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group system = file.openGroup("/metadata/system");
    
    // Check hostname exists
    REQUIRE(file.nameExists("/metadata/system/hostname"));
    H5::DataSet ds_host = system.openDataSet("hostname");
    H5::StrType str_type = ds_host.getStrType();
    std::string hostname;
    ds_host.read(hostname, str_type);
    REQUIRE(hostname.length() > 0);
    
    // Check OS information exists
    REQUIRE(file.nameExists("/metadata/system/os"));
    REQUIRE(file.nameExists("/metadata/system/cpu_model"));
    
    file.close();
}

TEST_CASE("HDF5Writer v2 handles large ion ensembles", "[hdf5][io][performance]") {
    std::string test_file = get_test_file();
    auto config = create_minimal_config();
    
    // Create large ensemble (1000 ions)
    std::vector<core::IonState> ions;
    ions.reserve(1000);
    
    for (int i = 0; i < 1000; ++i) {
        core::IonState ion;
        ion.species_id = (i % 3 == 0) ? "H3O+" : "NO2+";
        ion.pos = core::Vec3(0.001 * i, 0.0001 * i, 0.00001 * i);
        ion.vel = core::Vec3(100 + i * 0.1, 10, 0);
        ion.mass_kg = 19.02 * 1.66e-27;
        ion.ion_charge_C = 1.602e-19;
        ion.birth_time_s = 0.0;
        ion.current_domain_index = 0;
        ion.active = true;
        ions.push_back(ion);
    }
    
    // Should handle large ensemble without issues
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    REQUIRE_NOTHROW(io::HDF5Writer::create_file(test_file, config, ensemble, "test", "test"));
    
    // Append trajectory snapshots
    for (int step = 0; step < 5; ++step) {
        for (auto& ion : ions) {
            ion.pos.x += 0.0001;
        }
        REQUIRE_NOTHROW(io::HDF5Writer::append_trajectory(test_file, step * 1e-9, to_ensemble(ions)));
    }
    
    // Verify data integrity
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group traj = file.openGroup("/trajectory");
    
    H5::DataSet ds_pos = traj.openDataSet("positions");
    H5::DataSpace space = ds_pos.getSpace();
    hsize_t dims[3];
    space.getSimpleExtentDims(dims);
    REQUIRE(dims[0] == 5);     // 5 timesteps
    REQUIRE(dims[1] == 1000);  // 1000 ions
    REQUIRE(dims[2] == 3);     // x,y,z
    
    file.close();
}

TEST_CASE("HDF5Writer v2 SHA256 hashing with ConfigLoader integration", "[hdf5][io][integration][hash]") {
    std::string test_file = get_test_file();
    
    // Create a real config file
    std::string config_path = "/tmp/test_sha256_config.json";
    // Create a simple ion cloud file first
    std::string ion_cloud_path = "/tmp/test_ion_cloud.json";
    std::ofstream ion_file(ion_cloud_path);
    ion_file << R"({"ions": [
        {"species": "test", "mass": 100, "charge": 1, "pos": [0,0,0], "vel": [0,0,0]}
    ]})";
    ion_file.close();
    
    std::ofstream ofs(config_path);
    ofs << R"({
        "simulation": {"dt_s": 1e-9, "total_time_s": 1e-6, "integrator": "RK4", "write_interval": 100},
        "physics": {"collision_model": "NoCollisions"},
        "output": {"folder": "./output", "trajectory_file": "test.h5"},
        "ion_cloud": "/tmp/test_ion_cloud.json",
        "domains": [{
            "name": "test_domain",
            "domain_index": 0,
            "instrument": "IMS",
            "solver": "rk45",
            "geometry": {"length_m": 0.1, "radius_m": 0.01},
            "environment": {"pressure_Pa": 101325.0, "temperature_K": 300.0, "gas_species": "He"},
            "fields": {"dc": {"axial_V": 0.0}}
        }]
    })";
    ofs.close();
    
    // Load config via ConfigLoader (should set config_file_path)
    auto config = config::ConfigLoader::load(config_path);
    
    // Verify config_file_path is set
    REQUIRE(!config.config_file_path.empty());
    REQUIRE(std::filesystem::exists(config.config_file_path));
    
    // Create test ions
    auto ions = create_test_ions();
    
    // Create HDF5 file with SHA256 hashing
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    io::HDF5Writer::create_file(test_file, config, ensemble, "abc123", "gcc 11.4.0");
    
    // Open and verify SHA256 hash was written
    H5::H5File file(test_file, H5F_ACC_RDONLY);
    H5::Group hash_group = file.openGroup("/metadata/reproducibility/input_hash");
    
    // Read SHA256 hash
    H5::DataSet ds_hash = hash_group.openDataSet("config_sha256");
    H5::StrType str_type = ds_hash.getStrType();
    std::string stored_hash;
    ds_hash.read(stored_hash, str_type);
    
    // Should NOT be "N/A" or "ERROR" (should be actual hash)
    REQUIRE(stored_hash != "N/A");
    REQUIRE(stored_hash != "ERROR");
    
    // Should be valid SHA256 (64 hex characters)
    REQUIRE(stored_hash.length() == 64);
    REQUIRE(stored_hash.find_first_not_of("0123456789abcdef") == std::string::npos);
    
    file.close();
    
    // Cleanup
    std::filesystem::remove(config_path);
    std::filesystem::remove(ion_cloud_path);
}
