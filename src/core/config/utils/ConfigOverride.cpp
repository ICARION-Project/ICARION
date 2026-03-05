// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "ConfigOverride.h"
#include "../conversion/EnumMapper.h"
#include <iostream>
#include <algorithm>

namespace ICARION::config {

void ConfigOverride::apply(FullConfig& config, const std::map<std::string, std::string>& overrides) {
    if (overrides.empty()) {
        return;
    }
    
    std::cout << "[ConfigOverride] Applying " << overrides.size() << " override(s):\n";
    
    for (const auto& [key, value] : overrides) {
        std::cout << "  " << key << " = " << value << "\n";
        
        // === Simulation section ===
        if (key == "simulation.dt_s" || key == "simulation.timestep") {
            config.simulation.dt_s = parse_double(value, key);
        }
        else if (key == "simulation.total_time_s" || key == "simulation.total_time") {
            config.simulation.total_time_s = parse_double(value, key);
        }
        else if (key == "simulation.write_interval") {
            config.simulation.write_interval = parse_int(value, key);
        }
        else if (key == "simulation.rng_seed" || key == "simulation.seed") {
            config.simulation.rng_seed = parse_uint(value, key);
        }
        else if (key == "simulation.integrator") {
            config.simulation.integrator = value;
        }
        else if (key == "simulation.rk45_min_step_s") {
            config.simulation.rk45_min_step_s = parse_double(value, key);
        }
        else if (key == "simulation.enable_gpu") {
            config.simulation.enable_gpu = parse_bool(value, key);
        }
        else if (key == "simulation.enable_openmp") {
            config.simulation.enable_openmp = parse_bool(value, key);
        }
        
        // === Physics section ===
        else if (key == "physics.collision_model") {
            config.physics.collision_model = EnumMapper::parse_collision_model(value);
        }
        else if (key == "physics.enable_reactions") {
            config.physics.enable_reactions = parse_bool(value, key);
        }
        else if (key == "physics.enable_space_charge") {
            config.physics.enable_space_charge = parse_bool(value, key);
        }
        else if (key == "physics.enable_space_charge_gpu") {
            config.physics.enable_space_charge_gpu = parse_bool(value, key);
        }
        else if (key == "physics.enable_ou_thermalization") {
            config.physics.enable_ou_thermalization = parse_bool(value, key);
        }
        
        // === Output section ===
        else if (key == "output.folder") {
            config.output.folder = value;
        }
        else if (key == "output.trajectory_file" || key == "output.file") {
            config.output.trajectory_file = value;
        }
        else if (key == "output.print_progress") {
            config.output.print_progress = parse_bool(value, key);
        }
        else if (key == "output.buffer_byte_cap") {
            config.output.buffer_byte_cap = static_cast<size_t>(parse_uint64(value, key));
        }
        
        // === Database paths ===
        else if (key == "species_database" || key == "database.species") {
            config.species_database_path = value;
        }
        else if (key == "reaction_database" || key == "database.reactions") {
            config.reaction_database_path = value;
        }
        
        // === Unknown key ===
        else {
            std::cerr << "  Warning: Unknown config key '" << key << "' - ignored\n";
        }
    }
    
    std::cout << "[ConfigOverride] Overrides applied successfully\n";
}

// === Helper functions ===

double ConfigOverride::parse_double(const std::string& value, const std::string& key) {
    try {
        return std::stod(value);
    } catch (const std::exception& e) {
        throw std::runtime_error("Invalid double value for '" + key + "': " + value);
    }
}

int ConfigOverride::parse_int(const std::string& value, const std::string& key) {
    try {
        return std::stoi(value);
    } catch (const std::exception& e) {
        throw std::runtime_error("Invalid int value for '" + key + "': " + value);
    }
}

unsigned int ConfigOverride::parse_uint(const std::string& value, const std::string& key) {
    try {
        long val = std::stol(value);
        if (val < 0) {
            throw std::runtime_error("Negative value not allowed");
        }
        return static_cast<unsigned int>(val);
    } catch (const std::exception& e) {
        throw std::runtime_error("Invalid unsigned int value for '" + key + "': " + value);
    }
}

uint64_t ConfigOverride::parse_uint64(const std::string& value, const std::string& key) {
    try {
        long long val = std::stoll(value);
        if (val < 0) {
            throw std::runtime_error("Negative value not allowed");
        }
        return static_cast<uint64_t>(val);
    } catch (const std::exception& e) {
        throw std::runtime_error("Invalid uint64 value for '" + key + "': " + value);
    }
}

bool ConfigOverride::parse_bool(const std::string& value, const std::string& key) {
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    
    if (lower == "true" || lower == "1" || lower == "yes" || lower == "on") {
        return true;
    } else if (lower == "false" || lower == "0" || lower == "no" || lower == "off") {
        return false;
    } else {
        throw std::runtime_error("Invalid bool value for '" + key + "': " + value + 
                                " (expected: true/false, yes/no, on/off, 1/0)");
    }
}

} // namespace ICARION::config
