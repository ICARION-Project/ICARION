// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "SpeciesLoader.h"
#include <fstream>
#include <iostream>

namespace ICARION::config {

SpeciesDatabase SpeciesLoader::load(const std::filesystem::path& filepath) {
    // Open file
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open species database: " + filepath.string());
    }
    
    // Parse JSON
    Json::Value root;
    Json::CharReaderBuilder builder;
    std::string errs;
    if (!Json::parseFromStream(builder, file, &root, &errs)) {
        throw std::runtime_error("JSON parse error in " + filepath.string() + ": " + errs);
    }
    
    return load_from_json(root);
}

SpeciesDatabase SpeciesLoader::load_from_json(const Json::Value& json) {
    SpeciesDatabase db;
    
    // Expect top-level "species" object
    if (!json.isMember("species") || !json["species"].isObject()) {
        throw std::runtime_error("Species database must have top-level 'species' object");
    }
    
    const Json::Value& species_obj = json["species"];
    
    // Iterate over species dictionary
    for (auto it = species_obj.begin(); it != species_obj.end(); ++it) {
        std::string species_id = it.key().asString();
        const Json::Value& species_json = *it;
        
        SpeciesProperties species = parse_species(species_id, species_json);
        species.convert_to_SI();
        
        // Validate
        auto validation = species.validate();
        if (!validation.valid) {
            std::string error_msg = "Species '" + species_id + "' validation failed:\n";
            for (const auto& err : validation.errors) {
                error_msg += "  - " + err + "\n";
            }
            throw std::runtime_error(error_msg);
        }
        
        // Print warnings if any
        for (const auto& warning : validation.warnings) {
            std::cout << "⚠  Species '" + species_id + "': " << warning << "\n";
        }
        
        // Store
        db.species[species_id] = species;
    }
    
    // Validate database
    auto db_validation = db.validate();
    if (!db_validation.valid) {
        std::string error_msg = "Species database validation failed:\n";
        for (const auto& err : db_validation.errors) {
            error_msg += "  - " + err + "\n";
        }
        throw std::runtime_error(error_msg);
    }
    
    std::cout << "✓ Loaded " << db.size() << " species from database\n";
    
    return db;
}

SpeciesProperties SpeciesLoader::parse_species(const std::string& id, const Json::Value& json) {
    SpeciesProperties species;
    species.id = id;
    
    // Required fields
    if (!json.isMember("mass_amu") || !json["mass_amu"].isNumeric()) {
        throw std::runtime_error("Species '" + id + "': missing or invalid 'mass_amu'");
    }
    species.mass_amu = json["mass_amu"].asDouble();
    
    if (!json.isMember("charge") || !json["charge"].isInt()) {
        throw std::runtime_error("Species '" + id + "': missing or invalid 'charge'");
    }
    species.charge = json["charge"].asInt();
    
    // Optional fields (ions)
    species.mobility_cm2Vs = get_optional_double(json, "mobility_cm2Vs");
    species.CCS_A2 = get_optional_double(json, "CCS_A2");
    
    // Optional fields (neutrals)
    species.polarizability_A3 = get_optional_double(json, "polarizability_A3");
    
    // Optional metadata
    species.name = get_optional_string(json, "name");
    species.geometry_file = get_optional_string(json, "geometry_file");
    
    // Optional reference conditions
    species.reference_temperature_K = get_optional_double(json, "reference_temperature_K");
    species.reference_pressure_Pa = get_optional_double(json, "reference_pressure_Pa");
    species.ccs_method = get_optional_string(json, "ccs_method");

    return species;
}

std::optional<double> SpeciesLoader::get_optional_double(const Json::Value& json, const char* key) {
    if (json.isMember(key) && json[key].isNumeric()) {
        return json[key].asDouble();
    }
    return std::nullopt;
}

std::optional<std::string> SpeciesLoader::get_optional_string(const Json::Value& json, const char* key) {
    if (json.isMember(key) && json[key].isString()) {
        return json[key].asString();
    }
    return std::nullopt;
}

} // namespace ICARION::config
