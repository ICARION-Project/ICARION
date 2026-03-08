// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "ConfigLoader.h"
#include "SimulationConfigLoader.h"
#include "PhysicsConfigLoader.h"
#include "OutputConfigLoader.h"
#include "DomainConfigLoader.h"
#include "WaveformLoader.h"
#include <iostream>

namespace ICARION::config {

FullConfig ConfigLoader::load(const std::filesystem::path& config_path) {
    if (!std::filesystem::exists(config_path)) {
        throw std::runtime_error("Config file not found: " + config_path.string());
    }
    
    // Parse JSON
    Json::Value root = parse_json_file(config_path);
    
    // Base path for resolving relative paths
    auto base_path = config_path.parent_path();
    
    // Load from parsed JSON and store absolute config path
    auto config = load_from_json(root, base_path);
    config.config_file_path = std::filesystem::absolute(config_path).string();
    
    return config;
}

FullConfig ConfigLoader::load_from_json(const Json::Value& root, 
                                        const std::filesystem::path& base_path) {
    FullConfig config;
    
    // Load core sections
    if (root.isMember("simulation")) {
        config.simulation = SimulationConfigLoader::load(root["simulation"]);
    } else {
        throw std::runtime_error("Missing required 'simulation' section in config");
    }
    
    if (root.isMember("physics")) {
        config.physics = PhysicsConfigLoader::load(root["physics"]);
    } else {
        throw std::runtime_error("Missing required 'physics' section in config");
    }
    
    if (root.isMember("output")) {
        config.output = OutputConfigLoader::load(root["output"]);
    } else {
        throw std::runtime_error("Missing required 'output' section in config");
    }
    
    // Load global waveform library (SSOT)
    if (root.isMember("waveforms")) {
        config.waveforms = WaveformLoader::load_library(root["waveforms"]);
    }
    
    // Load domains with global integrator and global waveforms
    if (root.isMember("domains") && root["domains"].isArray()) {
        for (const auto& domain_json : root["domains"]) {
            config.domains.push_back(DomainConfigLoader::load(
                domain_json, 
                config.simulation.integrator,
                config.waveforms  // Pass global waveforms
            ));
        }
    } else {
        throw std::runtime_error("Missing or invalid 'domains' array in config");
    }
    
    // Load database paths
    if (root.isMember("species_database")) {
        config.species_database_path = resolve_path(
            root["species_database"].asString(), base_path);
    }
    
    if (root.isMember("reaction_database")) {
        config.reaction_database_path = resolve_path(
            root["reaction_database"].asString(), base_path);
    }
    
    // Ion configuration
    if (root.isMember("ions")) {
        config.ions = parse_ion_config(root["ions"], base_path);
    }
    
    // Optional metadata
    if (root.isMember("title")) {
        config.title = root["title"].asString();
    }
    
    // Note: config_file_path is set in load() function (not here, since we only have base_path)
    
    // Finalize all domains (compute derived quantities)
    config.finalize_all();
    
    // Load databases (species, reactions) if paths specified
    config.load_databases(base_path);
    
    // Validate complete configuration
    config.validate();
    
    std::cout << "[ConfigLoader] Loaded configuration: " 
              << config.domains.size() << " domains, "
              << config.species_db.size() << " species, "
              << config.reaction_db.size() << " reactions" << std::endl;
    
    return config;
}

Json::Value ConfigLoader::parse_json_file(const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open config file: " + path.string());
    }
    
    Json::Value root;
    Json::CharReaderBuilder builder;
    std::string errs;
    
    if (!Json::parseFromStream(builder, file, &root, &errs)) {
        throw std::runtime_error("JSON parse error in " + path.string() + ": " + errs);
    }
    
    return root;
}

std::string ConfigLoader::resolve_path(const std::string& path, 
                                       const std::filesystem::path& base_path) {
    if (path.empty()) {
        return "";
    }
    
    std::filesystem::path p(path);
    
    // If already absolute, return as-is
    if (p.is_absolute()) {
        return p.lexically_normal().string();
    }
    
    // Resolve relative to base path
    if (!base_path.empty()) {
        std::filesystem::path resolved = (base_path / p).lexically_normal();
        
        // Try canonical path if file exists (resolves symlinks and .. properly)
        std::error_code ec;
        auto canonical = std::filesystem::canonical(resolved, ec);
        if (!ec) {
            return canonical.string();
        }
        
        // Fallback to lexically normalized if file doesn't exist yet
        return resolved.string();
    }
    
    return path;
}

IonConfig ConfigLoader::parse_ion_config(const Json::Value& json,
                                         const std::filesystem::path& base_path) {
    IonConfig config;
    
    // Option 1: Load from file
    if (json.isMember("from_file")) {
        config.from_file = resolve_path(json["from_file"].asString(), base_path);
        return config;
    }
    
    // Option 2: Generate from species list
    if (json.isMember("species") && json["species"].isArray()) {
        for (const auto& spec_json : json["species"]) {
            IonSpeciesConfig spec;
            
            // Required fields (check both old 'id' and new 'species_id' formats)
            if ((!spec_json.isMember("id") && !spec_json.isMember("species_id")) || 
                !spec_json.isMember("count") ||
                !spec_json.isMember("position") || !spec_json.isMember("velocity")) {
                throw std::runtime_error("Ion species config missing required fields (species_id/id, count, position, velocity)");
            }
            
            // Support both old and new field names
            if (spec_json.isMember("species_id")) {
                spec.species_id = spec_json["species_id"].asString();
            } else {
                spec.species_id = spec_json["id"].asString();
            }
            spec.count = spec_json["count"].asInt();
            
            // Parse position distribution
            spec.position = parse_position_config(spec_json["position"]);
            
            // Parse velocity distribution
            spec.velocity = parse_velocity_config(spec_json["velocity"]);

            // Optional birth-time configuration
            // Fallback stays at t=0 s when nothing is provided.
            if (spec_json.isMember("birth_time_s")) {
                spec.birth_time_s = spec_json["birth_time_s"].asDouble();
            }

            if (spec_json.isMember("birth_time")) {
                const Json::Value& birth_json = spec_json["birth_time"];

                if (birth_json.isNumeric()) {
                    // Convenience: allow scalar shorthand
                    spec.birth_time_s = birth_json.asDouble();
                } else if (birth_json.isObject()) {
                    std::string type = birth_json.get("type", "gaussian").asString();

                    if (type == "fixed") {
                        if (birth_json.isMember("time_s")) {
                            spec.birth_time_s = birth_json["time_s"].asDouble();
                        } else if (birth_json.isMember("value_s")) {
                            spec.birth_time_s = birth_json["value_s"].asDouble();
                        } else {
                            throw std::runtime_error(
                                "birth_time.type='fixed' requires 'time_s' (or 'value_s')"
                            );
                        }
                    } else if (type == "gaussian" || type == "gaussian_truncated") {
                        if (!birth_json.isMember("t_min_s") || !birth_json.isMember("t_max_s")) {
                            throw std::runtime_error(
                                "birth_time gaussian requires 't_min_s' and 't_max_s'"
                            );
                        }

                        spec.use_birth_time_distribution = true;
                        spec.birth_time_distribution = IonSpeciesConfig::BirthTimeDistribution::GaussianTruncated;
                        spec.birth_time_min_s = birth_json["t_min_s"].asDouble();
                        spec.birth_time_max_s = birth_json["t_max_s"].asDouble();

                        if (spec.birth_time_max_s < spec.birth_time_min_s) {
                            throw std::runtime_error(
                                "birth_time requires t_max_s >= t_min_s"
                            );
                        }

                        const double window_s = spec.birth_time_max_s - spec.birth_time_min_s;
                        const double default_mean_s =
                            0.5 * (spec.birth_time_min_s + spec.birth_time_max_s);
                        const double default_std_s = (window_s > 0.0) ? (window_s / 6.0) : 0.0;

                        spec.birth_time_mean_s = birth_json.get("mean_s", default_mean_s).asDouble();
                        spec.birth_time_std_s = birth_json.get("std_s", default_std_s).asDouble();

                        if (spec.birth_time_std_s < 0.0) {
                            throw std::runtime_error("birth_time.std_s must be >= 0");
                        }
                    } else if (type == "uniform") {
                        if (!birth_json.isMember("t_min_s") || !birth_json.isMember("t_max_s")) {
                            throw std::runtime_error(
                                "birth_time uniform requires 't_min_s' and 't_max_s'"
                            );
                        }

                        spec.use_birth_time_distribution = true;
                        spec.birth_time_distribution = IonSpeciesConfig::BirthTimeDistribution::Uniform;
                        spec.birth_time_min_s = birth_json["t_min_s"].asDouble();
                        spec.birth_time_max_s = birth_json["t_max_s"].asDouble();

                        if (spec.birth_time_max_s < spec.birth_time_min_s) {
                            throw std::runtime_error(
                                "birth_time requires t_max_s >= t_min_s"
                            );
                        }

                        spec.birth_time_mean_s = 0.5 * (spec.birth_time_min_s + spec.birth_time_max_s);
                        spec.birth_time_std_s = 0.0;
                    } else {
                        throw std::runtime_error(
                            "Unknown birth_time type: " + type +
                            " (supported: fixed, gaussian, gaussian_truncated, uniform)"
                        );
                    }
                } else {
                    throw std::runtime_error(
                        "birth_time must be either a number or an object"
                    );
                }
            }
            
            config.species.push_back(spec);
        }
    }
    
    return config;
}

PositionConfig ConfigLoader::parse_position_config(const Json::Value& json) {
    PositionConfig config;
    
    if (!json.isMember("type")) {
        throw std::runtime_error("Position config missing 'type' field");
    }
    
    std::string type_str = json["type"].asString();
    
    if (type_str == "point") {
        config.type = PositionDistribution::Point;
        config.center = parse_vec3(json["center"]);
    }
    else if (type_str == "gaussian") {
        config.type = PositionDistribution::Gaussian;
        config.center = parse_vec3(json["center"]);
        config.std_dev = parse_vec3(json["std"]);
    }
    else if (type_str == "uniform_sphere") {
        config.type = PositionDistribution::UniformSphere;
        config.center = parse_vec3(json["center"]);
        config.radius = json["radius"].asDouble();
    }
    else if (type_str == "uniform_cylinder") {
        config.type = PositionDistribution::UniformCylinder;
        config.center = parse_vec3(json["center"]);
        config.cylinder_radius = json["radius"].asDouble();
        config.cylinder_length = json["length"].asDouble();
    }
    else if (type_str == "uniform_box") {
        config.type = PositionDistribution::UniformBox;
        config.box_min = parse_vec3(json["min"]);
        config.box_max = parse_vec3(json["max"]);
    }
    else {
        throw std::runtime_error("Unknown position distribution type: " + type_str);
    }
    
    return config;
}

VelocityConfig ConfigLoader::parse_velocity_config(const Json::Value& json) {
    VelocityConfig config;
    
    if (!json.isMember("type")) {
        throw std::runtime_error("Velocity config missing 'type' field");
    }
    
    std::string type_str = json["type"].asString();
    
    if (type_str == "fixed") {
        config.type = VelocityDistribution::Fixed;
        config.value = parse_vec3(json["value"]);
    }
    else if (type_str == "thermal") {
        config.type = VelocityDistribution::Thermal;
        config.temperature_K = json["temperature_K"].asDouble();
    }
    else if (type_str == "kinetic") {
        config.type = VelocityDistribution::Kinetic;
        config.energy_eV = json["energy_eV"].asDouble();
        config.direction = parse_vec3(json["direction"]);
        
        if (json.isMember("spread_angle_deg")) {
            config.spread_angle_deg = json["spread_angle_deg"].asDouble();
        }
    }
    else if (type_str == "gaussian") {
        config.type = VelocityDistribution::Gaussian;
        config.mean = parse_vec3(json["mean"]);
        config.std_dev = parse_vec3(json["std"]);
    }
    else {
        throw std::runtime_error("Unknown velocity distribution type: " + type_str);
    }
    
    return config;
}

Vec3 ConfigLoader::parse_vec3(const Json::Value& json) {
    if (!json.isArray() || json.size() != 3) {
        throw std::runtime_error("Vec3 must be array of 3 numbers");
    }
    
    return Vec3(
        json[0].asDouble(),
        json[1].asDouble(),
        json[2].asDouble()
    );
}

} // namespace ICARION::config
