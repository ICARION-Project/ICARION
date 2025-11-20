// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "ConfigLoader.h"
#include "SimulationConfigLoader.h"
#include "PhysicsConfigLoader.h"
#include "OutputConfigLoader.h"
#include "DomainConfigLoader.h"
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
    
    // Load from parsed JSON
    return load_from_json(root, base_path);
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
    
    // Load domains with global integrator as fallback
    if (root.isMember("domains") && root["domains"].isArray()) {
        for (const auto& domain_json : root["domains"]) {
            config.domains.push_back(DomainConfigLoader::load(domain_json, config.simulation.integrator));
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
    
    // Legacy ion_cloud path (deprecated, use "ions" instead)
    if (root.isMember("ion_cloud")) {
        config.ion_cloud_path = resolve_path(
            root["ion_cloud"].asString(), base_path);
    }
    
    // Modern ion configuration
    if (root.isMember("ions")) {
        config.ions = parse_ion_config(root["ions"], base_path);
    }
    
    // Optional metadata
    if (root.isMember("title")) {
        config.title = root["title"].asString();
    }
    
    // Store config file path
    config.config_file_path = base_path.string();
    
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
            
            // Required fields
            if (!spec_json.isMember("id") || !spec_json.isMember("count") ||
                !spec_json.isMember("position") || !spec_json.isMember("velocity")) {
                throw std::runtime_error("Ion species config missing required fields (id, count, position, velocity)");
            }
            
            spec.species_id = spec_json["id"].asString();
            spec.count = spec_json["count"].asInt();
            
            // Parse position distribution
            spec.position = parse_position_config(spec_json["position"]);
            
            // Parse velocity distribution
            spec.velocity = parse_velocity_config(spec_json["velocity"]);
            
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
