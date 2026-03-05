// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "IonLoader.h"
#include "utils/constants.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <json/json.h>

namespace ICARION::config {

IonGenerationResult IonLoader::generate_ions(
    const IonConfig& config,
    const SpeciesDatabase& species_db,
    std::mt19937& rng
) {
    IonGenerationResult result;
    
    if (!config.is_valid()) {
        result.validation.errors.push_back("Ion configuration invalid (no file or species specified)");
        result.validation.valid = false;
        return result;
    }
    
    // Option 1: Load from file
    if (config.from_file.has_value()) {
        return load_from_file(config.from_file.value(), species_db);
    }
    
    // Option 2: Generate from species list (each species has own boundaries)
    for (const auto& spec_config : config.species) {
        // Check species exists
        if (!species_db.has(spec_config.species_id)) {
            result.validation.errors.push_back("Species '" + spec_config.species_id + "' not found in database");
            result.validation.valid = false;
            continue;
        }
        
        try {
            auto species_ions = generate_species(spec_config, species_db, rng);
            result.ions.insert(result.ions.end(), species_ions.begin(), species_ions.end());
        } catch (const std::exception& e) {
            result.validation.errors.push_back("Failed to generate " + spec_config.species_id + ": " + e.what());
            result.validation.valid = false;
        }
    }
    
    if (result.validation.valid) {
        std::cout << "✓ Generated " << result.ions.size() << " ions from configuration\n";
    }
    
    return result;
}

IonGenerationResult IonLoader::load_from_file(
    const std::filesystem::path& filepath,
    const SpeciesDatabase& species_db
) {
    IonGenerationResult result;
    
    std::ifstream file(filepath);
    if (!file.is_open()) {
        result.validation.errors.push_back("Cannot open ion cloud file: " + filepath.string());
        result.validation.valid = false;
        return result;
    }
    
    Json::Value root;
    Json::CharReaderBuilder builder;
    std::string errs;
    if (!Json::parseFromStream(builder, file, &root, &errs)) {
        result.validation.errors.push_back("JSON parse error in " + filepath.string() + ": " + errs);
        result.validation.valid = false;
        return result;
    }
    
    if (!root.isMember("ions") || !root["ions"].isArray()) {
        result.validation.errors.push_back("Ion cloud file must have 'ions' array");
        result.validation.valid = false;
        return result;
    }
    
    const Json::Value& ions_array = root["ions"];
    result.ions.reserve(ions_array.size());
    
    for (const auto& ion_json : ions_array) {
        if (!ion_json.isMember("species")) {
            result.validation.warnings.push_back("Ion without 'species' field, skipping");
            continue;
        }
        
        std::string species_id = ion_json["species"].asString();
        
        if (!species_db.has(species_id)) {
            result.validation.warnings.push_back("Species '" + species_id + "' not in database, skipping ion");
            continue;
        }
        
        const auto& species = species_db.get(species_id);
        
        core::IonState ion;
        ion.species_id = species_id;
        ion.mass_kg = species.mass_kg;
        ion.ion_charge_C = species.charge_C;
        ion.reduced_mobility_cm2_Vs = species.mobility_m2Vs * 1e4;  // m²/Vs → cm²/Vs
        ion.CCS_m2 = species.CCS_m2;
        
        // Position
        if (ion_json.isMember("pos") && ion_json["pos"].isArray() && ion_json["pos"].size() >= 3) {
            const Json::Value& pos = ion_json["pos"];
            ion.pos.x = pos[0].asDouble();
            ion.pos.y = pos[1].asDouble();
            ion.pos.z = pos[2].asDouble();
        }
        
        // Velocity
        if (ion_json.isMember("vel") && ion_json["vel"].isArray() && ion_json["vel"].size() >= 3) {
            const Json::Value& vel = ion_json["vel"];
            ion.vel.x = vel[0].asDouble();
            ion.vel.y = vel[1].asDouble();
            ion.vel.z = vel[2].asDouble();
        }
        
        // Birth time
        ion.birth_time_s = ion_json.get("birth_time", 0.0).asDouble();
        ion.born = (ion.birth_time_s <= 0.0);
        
        ion.active = true;
        ion.history_index = -1;
        ion.current_domain_index = 0;
        
        result.ions.push_back(ion);
    }
    
    std::cout << "✓ Loaded " << result.ions.size() << " ions from " << filepath.filename() << "\n";
    
    return result;
}

std::vector<core::IonState> IonLoader::generate_species(
    const IonSpeciesConfig& spec_config,
    const SpeciesDatabase& species_db,
    std::mt19937& rng
) {
    // Caller already checked species exists
    const auto& species = species_db.get(spec_config.species_id);
    
    std::vector<core::IonState> ions;
    ions.reserve(spec_config.count);
    
    std::cout << "  Generating " << spec_config.count << " " << spec_config.species_id << " ions:\n";
    std::cout << "    Position: " << to_string(spec_config.position.type) << "\n";
    std::cout << "    Velocity: " << to_string(spec_config.velocity.type) << "\n";
    
    for (size_t i = 0; i < spec_config.count; ++i) {
        core::IonState ion;
        
        ion.species_id = spec_config.species_id;
        ion.mass_kg = species.mass_kg;
        ion.ion_charge_C = species.charge_C;
        ion.reduced_mobility_cm2_Vs = species.mobility_m2Vs * 1e4;
        ion.CCS_m2 = species.CCS_m2;
        
        // Sample position (within species-specific boundaries)
        ion.pos = sample_position(spec_config.position, rng);
        
        // Sample velocity
        ion.vel = sample_velocity(spec_config.velocity, species.mass_kg, rng);
        
        ion.birth_time_s = spec_config.birth_time_s;
        
        // SSOT: Birth logic in ONE place
        if (ion.birth_time_s <= 0.0) {
            // Immediate birth (t=0)
            ion.born = true;
            ion.active = true;
        } else {
            // Delayed birth (birth_time_s > 0)
            ion.born = false;
            ion.active = false;  // Will be activated by apply_ion_birth()
        }
        
        ion.history_index = -1;
        ion.current_domain_index = 0;
        
        ions.push_back(ion);
    }
    
    return ions;
}

Vec3 IonLoader::sample_position(const PositionConfig& config, std::mt19937& rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    switch (config.type) {
        case PositionDistribution::Point:
            return config.center;
        
        case PositionDistribution::Gaussian:
            return {
                config.center.x + normal(rng) * config.std_dev.x,
                config.center.y + normal(rng) * config.std_dev.y,
                config.center.z + normal(rng) * config.std_dev.z
            };
        
        case PositionDistribution::UniformSphere: {
            // Uniform in sphere: r³ distribution for radius, uniform angles
            double r = config.radius * std::cbrt(uniform(rng));
            double theta = 2.0 * M_PI * uniform(rng);
            double phi = std::acos(2.0 * uniform(rng) - 1.0);
            
            return {
                config.center.x + r * std::sin(phi) * std::cos(theta),
                config.center.y + r * std::sin(phi) * std::sin(theta),
                config.center.z + r * std::cos(phi)
            };
        }
        
        case PositionDistribution::UniformCylinder: {
            // Uniform in cylinder: r² distribution for radius, uniform z
            double r = config.cylinder_radius * std::sqrt(uniform(rng));
            double theta = 2.0 * M_PI * uniform(rng);
            double z = config.cylinder_length * (uniform(rng) - 0.5);
            
            return {
                config.center.x + r * std::cos(theta),
                config.center.y + r * std::sin(theta),
                config.center.z + z
            };
        }
        
        case PositionDistribution::UniformBox:
            return {
                config.box_min.x + uniform(rng) * (config.box_max.x - config.box_min.x),
                config.box_min.y + uniform(rng) * (config.box_max.y - config.box_min.y),
                config.box_min.z + uniform(rng) * (config.box_max.z - config.box_min.z)
            };
        
        default:
            throw std::runtime_error("Unknown position distribution type");
    }
}

Vec3 IonLoader::sample_velocity(const VelocityConfig& config, double ion_mass_kg, std::mt19937& rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    switch (config.type) {
        case VelocityDistribution::Fixed:
            return config.value;
        
        case VelocityDistribution::Thermal: {
            // Maxwell-Boltzmann: ALWAYS random directions (no directed drift)
            double sigma = std::sqrt(BOLTZMANN_CONSTANT * config.temperature_K / ion_mass_kg);
            return {
                normal(rng) * sigma,
                normal(rng) * sigma,
                normal(rng) * sigma
            };
        }
        
        case VelocityDistribution::Kinetic: {
            // Fixed energy, directed with optional spread
            double speed = std::sqrt(2.0 * config.energy_eV * ELEM_CHARGE_C / ion_mass_kg);
            
            // Normalize direction
            double dir_mag = std::sqrt(config.direction.x * config.direction.x +
                                     config.direction.y * config.direction.y +
                                     config.direction.z * config.direction.z);
            if (dir_mag < 1e-10) {
                // Use default z-direction if zero
                return {0, 0, speed};
            }
            
            Vec3 dir = {config.direction.x / dir_mag, 
                       config.direction.y / dir_mag, 
                       config.direction.z / dir_mag};
            
            // Apply spread angle (cone)
            if (config.spread_angle_deg > 0.0) {
                double max_theta = config.spread_angle_deg * M_PI / 180.0;
                double theta = max_theta * std::sqrt(uniform(rng));  // Uniform solid angle
                double phi = 2.0 * M_PI * uniform(rng);
                
                // Create perpendicular vectors for rotation
                Vec3 perp1, perp2;
                if (std::abs(dir.z) < 0.9) {
                    perp1 = {-dir.y, dir.x, 0.0};
                } else {
                    perp1 = {0.0, -dir.z, dir.y};
                }
                double perp1_mag = std::sqrt(perp1.x * perp1.x + perp1.y * perp1.y + perp1.z * perp1.z);
                perp1.x /= perp1_mag;
                perp1.y /= perp1_mag;
                perp1.z /= perp1_mag;
                
                // perp2 = dir × perp1
                perp2.x = dir.y * perp1.z - dir.z * perp1.y;
                perp2.y = dir.z * perp1.x - dir.x * perp1.z;
                perp2.z = dir.x * perp1.y - dir.y * perp1.x;
                
                // Rotate direction by theta around random axis in perp plane
                double cos_theta = std::cos(theta);
                double sin_theta = std::sin(theta);
                double cos_phi = std::cos(phi);
                double sin_phi = std::sin(phi);
                
                dir.x = dir.x * cos_theta + (perp1.x * cos_phi + perp2.x * sin_phi) * sin_theta;
                dir.y = dir.y * cos_theta + (perp1.y * cos_phi + perp2.y * sin_phi) * sin_theta;
                dir.z = dir.z * cos_theta + (perp1.z * cos_phi + perp2.z * sin_phi) * sin_theta;
            }
            
            return {dir.x * speed, dir.y * speed, dir.z * speed};
        }
        
        case VelocityDistribution::Gaussian:
            return {
                config.mean.x + normal(rng) * config.std_dev.x,
                config.mean.y + normal(rng) * config.std_dev.y,
                config.mean.z + normal(rng) * config.std_dev.z
            };
        
        default:
            throw std::runtime_error("Unknown velocity distribution type");
    }
}

} // namespace ICARION::config
