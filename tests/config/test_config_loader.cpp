// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/config/loader/ConfigLoader.h"
#include <fstream>
#include <filesystem>
#include <cmath>

using namespace ICARION::config;
using Catch::Approx;

namespace {
    std::string create_temp_config(const std::string& content, const std::string& suffix = "") {
        std::string path = "/tmp/icarion_test_config" + suffix + ".json";
        std::ofstream file(path);
        file << content;
        file.close();
        return path;
    }
}

// ============================================================================
// Basic Config Loading
// ============================================================================

TEST_CASE("ConfigLoader loads minimal valid config", "[config][loader]") {
    std::string config = R"({
        "simulation": {
            "dt_s": 1e-9,
            "total_time_s": 1e-6,
            "integrator": "RK4",
            "write_interval": 100
        },
        "physics": {
            "collision_model": "NoCollisions"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 5,
                    "position": {
                        "type": "point",
                        "center": [0.0, 0.0, 0.001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": 300.0
                    }
                }
            ]
        },
        "domains": [
            {
                "name": "test_domain",
                "domain_index": 0,
                "instrument": "IMS",
                "solver": "rk45",
                "geometry": {
                    "length_m": 0.1,
                    "radius_m": 0.01
                },
                "environment": {
                    "pressure_Pa": 101325.0,
                    "temperature_K": 300.0,
                    "gas_species": "He"
                },
                "fields": {
                    "dc": {"axial_V": 0.0}
                }
            }
        ]
    })";
    
    std::string path = create_temp_config(config);
    FullConfig cfg = ConfigLoader::load(path);
    
    REQUIRE(cfg.simulation.dt_s == Approx(1e-9));
    REQUIRE(cfg.physics.collision_model == CollisionModel::NoCollisions);
    REQUIRE(cfg.output.folder == "./output");
    REQUIRE(cfg.domains.size() == 1);
    
    // Verify config_file_path is stored as absolute path (for SHA256 hashing)
    REQUIRE(!cfg.config_file_path.empty());
    REQUIRE(cfg.config_file_path[0] == '/');  // Absolute path starts with /
    REQUIRE(cfg.config_file_path.find(".json") != std::string::npos);
    REQUIRE(std::filesystem::exists(cfg.config_file_path));
}

TEST_CASE("ConfigLoader loads config with two domains", "[config][loader][domains]") {
    std::string config = R"({
        "simulation": {
            "dt_s": 5e-10,
            "total_time_s": 2e-6,
            "integrator": "RK4",
            "write_interval": 200,
            "rng_seed": 12345
        },
        "physics": {
            "collision_model": "EHSS",
            "enable_reactions": false,
            "enable_space_charge": false
        },
        "output": {
            "folder": "./test_output",
            "trajectory_file": "trajectories.h5",
            "print_progress": true
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 5,
                    "position": {
                        "type": "point",
                        "center": [0.0, 0.0, 0.001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": 300.0
                    }
                }
            ]
        },
        "domains": [
            {
                "name": "drift_region",
                "domain_index": 0,
                "instrument": "IMS",
                "solver": "rk45",
                "geometry": {
                    "length_m": 0.15,
                    "radius_m": 0.012,
                    "origin_m": [0.0, 0.0, 0.0]
                },
                "environment": {
                    "pressure_Pa": 101325.0,
                    "temperature_K": 300.0,
                    "gas_species": "He",
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "axial_V": 250.0,
                        "EN_Td": 10.0
                    },
                    "RF": {
                        "voltage_V": 0.0,
                        "frequency_MHz": 0.0
                    }
                }
            },
            {
                "name": "tof_region",
                "domain_index": 1,
                "instrument": "TOF",
                "solver": "rk4",
                "geometry": {
                    "length_m": 0.5,
                    "radius_m": 0.02,
                    "origin_m": [0.0, 0.0, 0.15]
                },
                "environment": {
                    "pressure_Pa": 1e-6,
                    "temperature_K": 300.0,
                    "gas_species": "He",
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "axial_V": 5000.0,
                        "EN_Td": 0.0
                    }
                }
            }
        ]
    })";
    
    std::string path = create_temp_config(config, "_two_domains");
    FullConfig cfg = ConfigLoader::load(path);
    
    // Simulation
    REQUIRE(cfg.simulation.dt_s == Approx(5e-10));
    REQUIRE(cfg.simulation.rng_seed == 12345);
    
    // Physics
    REQUIRE(cfg.physics.collision_model == CollisionModel::EHSS);
    REQUIRE_FALSE(cfg.physics.enable_reactions);
    REQUIRE_FALSE(cfg.physics.enable_space_charge);
    
    // Output
    REQUIRE(cfg.output.folder == "./test_output");
    REQUIRE(cfg.output.trajectory_file == "trajectories.h5");
    REQUIRE(cfg.output.print_progress);
    
    // Domains
    REQUIRE(cfg.domains.size() == 2);
    
    // Domain 0: IMS
    REQUIRE(cfg.domains[0].domain_index == 0);
    REQUIRE(cfg.domains[0].geometry.length_m == Approx(0.15));
    REQUIRE(cfg.domains[0].geometry.radius_m == Approx(0.012));
    REQUIRE(cfg.domains[0].environment.pressure_Pa == Approx(101325.0));
    REQUIRE(cfg.domains[0].environment.temperature_K == Approx(300.0));
    REQUIRE(cfg.domains[0].environment.gas_species == "He");
    REQUIRE(cfg.domains[0].fields.dc.axial_V.constant_value == Approx(250.0));
    REQUIRE(cfg.domains[0].fields.dc.EN_Td.constant_value == Approx(10.0));
    
    // Domain 1: TOF
    REQUIRE(cfg.domains[1].domain_index == 1);
    REQUIRE(cfg.domains[1].geometry.length_m == Approx(0.5));
    REQUIRE(cfg.domains[1].geometry.radius_m == Approx(0.02));
    REQUIRE(cfg.domains[1].geometry.origin_m.z == Approx(0.15));
    REQUIRE(cfg.domains[1].environment.pressure_Pa == Approx(1e-6));
    REQUIRE(cfg.domains[1].fields.dc.axial_V.constant_value == Approx(5000.0));
    
    // Verify derived properties were computed
    REQUIRE(cfg.domains[0].environment.particle_density_m_3 > 0.0);
    REQUIRE(cfg.domains[0].environment.mean_thermal_velocity_m_s > 0.0);
    REQUIRE(cfg.domains[1].environment.particle_density_m_3 > 0.0);
}

TEST_CASE("ConfigLoader supports pressure waveform in environment", "[config][loader][waveform][environment]") {
    std::string config = R"({
        "simulation": {
            "dt_s": 1e-8,
            "total_time_s": 1e-6,
            "integrator": "RK4",
            "write_interval": 100
        },
        "physics": {
            "collision_model": "NoCollisions"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "waveforms": {
            "pumpdown": {
                "type": "exponential",
                "offset": 1000.0,
                "amplitude": 1000.0,
                "rate_per_s": -1000.0,
                "start_time_s": 0.0,
                "end_time_s": 1.0,
                "clamp": true
            }
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 1,
                    "position": { "type": "point", "center": [0.0, 0.0, 0.001] },
                    "velocity": { "type": "thermal", "temperature_K": 300.0 }
                }
            ]
        },
        "domains": [
            {
                "name": "ims_domain",
                "instrument": "IMS",
                "geometry": {
                    "length_m": 0.05,
                    "radius_m": 0.01
                },
                "environment": {
                    "pressure_Pa": "@pumpdown",
                    "temperature_K": 300.0,
                    "gas_species": "He"
                },
                "fields": {
                    "DC": {"axial_V": 0.0}
                }
            }
        ]
    })";

    std::string path = create_temp_config(config, "_pressure_waveform");
    FullConfig cfg = ConfigLoader::load(path);

    REQUIRE(cfg.domains.size() == 1);
    REQUIRE(cfg.domains[0].environment.has_dynamic_pressure());
    REQUIRE(cfg.domains[0].environment.pressure_Pa == Approx(2000.0));

    cfg.domains[0].environment.update_time_dependent(0.001, cfg.waveforms);
    REQUIRE(cfg.domains[0].environment.pressure_Pa == Approx(1000.0 + 1000.0 * std::exp(-1.0)).epsilon(1e-9));
    REQUIRE(cfg.domains[0].environment.particle_density_m_3 > 0.0);
}

TEST_CASE("ConfigLoader handles file not found", "[config][loader]") {
    REQUIRE_THROWS_AS(
        ConfigLoader::load("/nonexistent/file.json"),
        std::runtime_error
    );
}

TEST_CASE("ConfigLoader requires dt_s and total_time_s", "[config][loader][validation]") {
    // Missing dt_s
    std::string config_no_dt = R"({
        "simulation": {
            "total_time_s": 1e-3
        },
        "physics": {
            "collision_model": "HSD"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "domains": [
            {
                "name": "test",
                "instrument": "IMS",
                "geometry": { "length_m": 0.1, "radius_m": 0.01, "origin_m": [0,0,0] },
                "environment": { "pressure_Pa": 101325, "temperature_K": 300, "gas_species": "He" }
            }
        ]
    })";
    
    std::string path_no_dt = create_temp_config(config_no_dt, "_no_dt");
    REQUIRE_THROWS_AS(ConfigLoader::load(path_no_dt), std::runtime_error);
    
    // Missing total_time_s
    std::string config_no_time = R"({
        "simulation": {
            "dt_s": 1e-9
        },
        "physics": {
            "collision_model": "HSD"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "domains": [
            {
                "name": "test",
                "instrument": "IMS",
                "geometry": { "length_m": 0.1, "radius_m": 0.01, "origin_m": [0,0,0] },
                "environment": { "pressure_Pa": 101325, "temperature_K": 300, "gas_species": "He" }
            }
        ]
    })";
    
    std::string path_no_time = create_temp_config(config_no_time, "_no_time");
    REQUIRE_THROWS_AS(ConfigLoader::load(path_no_time), std::runtime_error);
}
