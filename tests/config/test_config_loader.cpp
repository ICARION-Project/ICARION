// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "core/config/loader/ConfigLoader.h"
#include <fstream>
#include <chrono>
#include <filesystem>
#include <cmath>
#include <sstream>

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

    std::string make_lqit_config_with_fields(const std::string& fields_json) {
        return R"({
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
                        "count": 1,
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
                    "name": "lqit",
                    "instrument": "LQIT",
                    "geometry": {
                        "length_m": 0.1,
                        "radius_m": 0.005
                    },
                    "environment": {
                        "pressure_Pa": 2.0,
                        "temperature_K": 300.0,
                        "gas_species": "N2"
                    },
                    "fields": )" + fields_json + R"(
                }
            ]
        })";
    }
}

// ============================================================================
// Basic Config Loading
// ============================================================================

TEST_CASE("ConfigLoader resolves optional user annotations", "[config][annotation]") {
    Json::CharReaderBuilder builder;
    Json::Value root;
    std::string errors;
    std::istringstream input(make_lqit_config_with_fields(R"({"type":"uniform","electric_Vm":[0,0,0]})"));
    REQUIRE(Json::parseFromStream(builder, input, &root, &errors));

    SECTION("inline note") {
        root["metadata"]["note"] = "config note\nsecond line";
        auto config = ConfigLoader::load_from_json(root, "/tmp");
        REQUIRE(config.user_annotation.note == "config note\nsecond line");
        REQUIRE(config.user_annotation.source == "inline");
    }
    SECTION("file note relative to config directory") {
        const auto dir = std::filesystem::temp_directory_path() /
            ("icarion_config_note_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        std::filesystem::create_directories(dir);
        const auto path = dir / "note.md";
        std::ofstream(path, std::ios::binary) << "file annotation\n";
        root["metadata"]["note_file"] = path.filename().string();
        auto config = ConfigLoader::load_from_json(root, dir);
        REQUIRE(config.user_annotation.note == "file annotation\n");
        REQUIRE(config.user_annotation.source_filename == path.filename().string());
        std::filesystem::remove_all(dir);
    }
    SECTION("mutually exclusive and unknown fields fail") {
        root["metadata"]["note"] = "one";
        root["metadata"]["note_file"] = "two";
        REQUIRE_THROWS(ConfigLoader::load_from_json(root, "/tmp"));
        root.removeMember("metadata");
        root["metadata"]["unknown"] = true;
        REQUIRE_THROWS(ConfigLoader::load_from_json(root, "/tmp"));
    }
    SECTION("no metadata remains valid") {
        auto config = ConfigLoader::load_from_json(root, "/tmp");
        REQUIRE_FALSE(config.user_annotation.present);
    }
}

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

TEST_CASE("ConfigLoader parses space-charge model enum", "[config][loader][physics]") {
    std::string config = R"({
        "simulation": {
            "dt_s": 1e-9,
            "total_time_s": 1e-6,
            "integrator": "RK4",
            "write_interval": 100
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_space_charge": true,
            "space_charge_model": "direct"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 1,
                    "position": {"type": "point", "center": [0.0, 0.0, 0.001]},
                    "velocity": {"type": "thermal", "temperature_K": 300.0}
                }
            ]
        },
        "domains": [
            {
                "name": "test_domain",
                "instrument": "IMS",
                "geometry": {"length_m": 0.1, "radius_m": 0.01},
                "environment": {"pressure_Pa": 101325.0, "temperature_K": 300.0, "gas_species": "He"},
                "fields": {"dc": {"axial_V": 0.0}}
            }
        ]
    })";

    std::string path = create_temp_config(config, "_space_charge_model");
    FullConfig cfg = ConfigLoader::load(path);

    REQUIRE(cfg.physics.space_charge_model_type == SpaceChargeModel::Direct);
    REQUIRE(cfg.physics.space_charge_model == "direct");
}

TEST_CASE("ConfigLoader rejects unknown space-charge model", "[config][loader][physics]") {
    std::string config = R"({
        "simulation": {
            "dt_s": 1e-9,
            "total_time_s": 1e-6,
            "integrator": "RK4",
            "write_interval": 100
        },
        "physics": {
            "collision_model": "NoCollisions",
            "space_charge_model": "fmm"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 1,
                    "position": {"type": "point", "center": [0.0, 0.0, 0.001]},
                    "velocity": {"type": "thermal", "temperature_K": 300.0}
                }
            ]
        },
        "domains": [
            {
                "name": "test_domain",
                "instrument": "IMS",
                "geometry": {"length_m": 0.1, "radius_m": 0.01},
                "environment": {"pressure_Pa": 101325.0, "temperature_K": 300.0, "gas_species": "He"},
                "fields": {"dc": {"axial_V": 0.0}}
            }
        ]
    })";

    std::string path = create_temp_config(config, "_space_charge_model_bad");
    REQUIRE_THROWS_WITH(
        ConfigLoader::load(path),
        Catch::Matchers::ContainsSubstring("Unknown space_charge_model"));
}

TEST_CASE("ConfigLoader parses deep collision analysis output config", "[config][loader][output]") {
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
            "trajectory_file": "test.h5",
            "deep_analysis": {
                "mode": "sampled_events",
                "domain_filter_index": 0,
                "sample_every_n": 4,
                "max_events_per_ion": 7
            }
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 1,
                    "position": {"type": "point", "center": [0.0, 0.0, 0.001]},
                    "velocity": {"type": "thermal", "temperature_K": 300.0}
                }
            ]
        },
        "domains": [
            {
                "name": "test_domain",
                "domain_index": 0,
                "instrument": "IMS",
                "solver": "rk4",
                "geometry": {"length_m": 0.1, "radius_m": 0.01},
                "environment": {"pressure_Pa": 101325.0, "temperature_K": 300.0, "gas_species": "He"},
                "fields": {"dc": {"axial_V": 0.0}}
            }
        ]
    })";

    std::string path = create_temp_config(config, "_deep_analysis");
    FullConfig cfg = ConfigLoader::load(path);

    REQUIRE(cfg.output.deep_analysis.mode_type == DeepAnalysisMode::SampledEvents);
    REQUIRE(cfg.output.deep_analysis.mode == "sampled_events");
    REQUIRE(cfg.output.deep_analysis.domain_filter_index == 0);
    REQUIRE(cfg.output.deep_analysis.sample_every_n == 4);
    REQUIRE(cfg.output.deep_analysis.max_events_per_ion == 7);
}

TEST_CASE("ConfigLoader parses TIMS instrument, field program, and axial flow", "[config][loader][tims]") {
    std::string config = R"({
        "simulation": {
            "dt_s": 1e-8,
            "total_time_s": 1e-6,
            "integrator": "RK4",
            "write_interval": 10
        },
        "physics": {
            "collision_model": "HSS"
        },
        "output": {
            "folder": "./output",
            "trajectory_file": "tims_loader_test.h5"
        },
        "ions": {
            "species": [{
                "id": "H3O+",
                "count": 1,
                "position": {"type": "point", "center": [0.0, 0.0, 0.01]},
                "velocity": {"type": "thermal", "temperature_K": 300.0}
            }]
        },
        "domains": [{
            "name": "tims_domain",
            "instrument": "TIMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": 0.02,
                "radius_m": 0.005
            },
            "environment": {
                "pressure_Pa": 310.0,
                "temperature_K": 300.0,
                "gas_species": "N2",
                "flow_model": "axial_uniform",
                "flow_parameters": {
                    "axial_flow_velocity_m_s": 134.0
                }
            },
            "fields": {
                "RF": {
                    "voltage_V": 180.0,
                    "frequency_Hz": 2000000.0
                },
                "TIMS": {
                    "enabled": true,
                    "z_positions_m": [0.0, 0.02],
                    "axial_field_initial_profile_V_m": [0.0, -4495.0],
                    "axial_field_final_profile_V_m": [0.0, 0.0],
                    "ramp_start_s": 8e-4,
                    "ramp_end_s": 1.8e-3,
                    "ramp_mode": "linear"
                }
            }
        }]
    })";

    std::string path = create_temp_config(config, "_tims_loader");
    auto cfg = ConfigLoader::load(path);

    REQUIRE(cfg.domains.size() == 1);
    const auto& domain = cfg.domains[0];
    REQUIRE(domain.instrument == Instrument::TIMS);
    REQUIRE(domain.environment.flow_model == FlowModelKind::AxialUniform);
    REQUIRE(domain.environment.axial_flow_velocity_m_s == Approx(134.0));
    REQUIRE(domain.environment.gas_velocity_m_s.z == Approx(134.0));
    REQUIRE(domain.fields.tims.enabled);
    REQUIRE(domain.fields.tims.axial_field_initial_profile_V_m[1] == Approx(-4495.0));
    REQUIRE(domain.fields.tims.ramp_start_s == Approx(8e-4));
    REQUIRE(domain.fields.tims.ramp_end_s == Approx(1.8e-3));

    std::filesystem::remove(path);
}

TEST_CASE("ConfigLoader loads LQIT two-axis dipolar AC excitation", "[config][loader][lqit][ac]") {
    std::string config = make_lqit_config_with_fields(R"({
        "RF": {
            "voltage_V": 100.0,
            "frequency_Hz": 1e6
        },
        "AC": {
            "voltage_V": 10.0,
            "frequency_Hz": 2e5
        },
        "dipolar_excitation": {
            "x": {
                "enabled": true,
                "amplitude_V": 10.0,
                "frequency_Hz": 2e5,
                "phase_rad": 0.0
            },
            "y": {
                "enabled": true,
                "amplitude_V": 2.5,
                "frequency_Hz": 2e5,
                "phase_rad": 1.57079632679,
                "ramp": {
                    "type": "linear",
                    "start": 0.0,
                    "end": 1.0,
                    "start_time_s": 0.0,
                    "end_time_s": 1e-6
                }
            }
        }
    })");

    std::string path = create_temp_config(config, "_dipolar_ac");
    FullConfig cfg = ConfigLoader::load(path);

    const auto& ac = cfg.domains[0].fields.ac;
    REQUIRE(ac.voltage_V.constant_value.value() == Approx(10.0));
    REQUIRE(ac.frequency_Hz.constant_value.value() == Approx(2e5));
    REQUIRE(ac.dipolar_excitation_defined);
    REQUIRE(ac.dipolar_x.enabled);
    REQUIRE(ac.dipolar_y.enabled);
    REQUIRE(ac.dipolar_x.amplitude_V.constant_value.value() == Approx(10.0));
    REQUIRE(ac.dipolar_y.amplitude_V.constant_value.value() == Approx(2.5));
    REQUIRE(ac.dipolar_y.phase_rad == Approx(1.57079632679));
    REQUIRE(ac.dipolar_y.ramp.evaluate(0.5e-6, cfg.domains[0].fields.waveform_library) == Approx(0.5));
}

TEST_CASE("ConfigLoader rejects duplicate LQIT dipolar AC sources", "[config][loader][lqit][ac][ssot]") {
    std::string config = make_lqit_config_with_fields(R"({
        "RF": {
            "voltage_V": 100.0,
            "frequency_Hz": 1e6
        },
        "AC": {
            "voltage_V": 10.0,
            "frequency_Hz": 2e5,
            "x": {
                "enabled": true,
                "amplitude_V": 10.0,
                "frequency_Hz": 2e5
            }
        },
        "dipolar_excitation": {
            "x": {
                "enabled": true,
                "amplitude_V": 10.0,
                "frequency_Hz": 2e5
            }
        }
    })");

    std::string path = create_temp_config(config, "_dipolar_ac_duplicate_sources");
    REQUIRE_THROWS_AS(ConfigLoader::load(path), std::runtime_error);
}

TEST_CASE("ConfigLoader rejects duplicate LQIT dipolar axis amplitude aliases", "[config][loader][lqit][ac][ssot]") {
    std::string config = make_lqit_config_with_fields(R"({
        "RF": {
            "voltage_V": 100.0,
            "frequency_Hz": 1e6
        },
        "AC": {
            "voltage_V": 10.0,
            "frequency_Hz": 2e5
        },
        "dipolar_excitation": {
            "x": {
                "enabled": true,
                "amplitude_V": 10.0,
                "voltage_V": 10.0,
                "frequency_Hz": 2e5
            }
        }
    })");

    std::string path = create_temp_config(config, "_dipolar_ac_duplicate_axis_amplitude");
    REQUIRE_THROWS_AS(ConfigLoader::load(path), std::runtime_error);
}

TEST_CASE("ConfigLoader does not load global reactions when reactions are disabled", "[config][loader][reactions]") {
    const auto species_db_path = (std::filesystem::current_path() / "data" / "species_database_v1.json").string();
    std::string config = R"({
        "simulation": {
            "dt_s": 1e-9,
            "total_time_s": 1e-6,
            "integrator": "RK4",
            "write_interval": 100
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_reactions": false
        },
        "species_database": ")" + species_db_path + R"(",
        "output": {
            "folder": "./output",
            "trajectory_file": "test.h5"
        },
        "ions": {
            "species": [
                {
                    "id": "H3O+",
                    "count": 1,
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

    std::string path = create_temp_config(config, "_reactions_disabled");
    FullConfig cfg = ConfigLoader::load(path);

    REQUIRE_FALSE(cfg.physics.enable_reactions);
    REQUIRE(cfg.species_db.size() > 0);
    REQUIRE(cfg.reaction_db.size() == 0);
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
