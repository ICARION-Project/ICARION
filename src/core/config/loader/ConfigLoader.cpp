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
    
    if (root.isMember("ion_cloud")) {
        config.ion_cloud_path = resolve_path(
            root["ion_cloud"].asString(), base_path);
    }
    
    // Optional metadata
    if (root.isMember("title")) {
        config.title = root["title"].asString();
    }
    
    // Store config file path
    config.config_file_path = base_path.string();
    
    // Finalize all domains (compute derived quantities)
    config.finalize_all();
    
    // Validate complete configuration
    config.validate();
    
    std::cout << "[ConfigLoader] Loaded configuration: " 
              << config.domains.size() << " domains" << std::endl;
    
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
        return path;
    }
    
    // Resolve relative to base path
    if (!base_path.empty()) {
        return (base_path / p).string();
    }
    
    return path;
}

} // namespace ICARION::config
