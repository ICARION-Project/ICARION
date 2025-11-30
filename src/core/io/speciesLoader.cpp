// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "core/io/speciesLoader.h"
#include "core/log/Logger.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <json/json.h>

namespace ICARION {
namespace io {

// -----------------------------
// Species member functions
// -----------------------------

double Species::calculate_mobility(double temperature_K, double neutral_mass_kg) const {
    if (CCS_m2 <= 0.0) {
        throw std::runtime_error("Cannot calculate mobility: CCS not set for species '" + id + "'");
    }
    if (std::abs(charge_C) < 1e-30) {
        throw std::runtime_error("Cannot calculate mobility: species '" + id + "' is neutral");
    }
    
    // Mason-Schamp equation: K = (3*q)/(16*N) * sqrt(2*pi / (k*T*mu)) / CCS
    // where mu is reduced mass
    double reduced_mass = (mass_kg * neutral_mass_kg) / (mass_kg + neutral_mass_kg);
    double numerator = 3.0 * std::abs(charge_C) / (16.0 * LOSCHMIDT_CONSTANT);
    double sqrt_term = std::sqrt(2.0 * M_PI / (BOLTZMANN_CONSTANT * temperature_K * reduced_mass));
    
    return numerator * sqrt_term / CCS_m2;
}

double Species::calculate_CCS(double temperature_K, double neutral_mass_kg) const {
    if (mobility_m2Vs <= 0.0) {
        throw std::runtime_error("Cannot calculate CCS: mobility not set for species '" + id + "'");
    }
    if (std::abs(charge_C) < 1e-30) {
        throw std::runtime_error("Cannot calculate CCS: species '" + id + "' is neutral");
    }
    
    // Invert Mason-Schamp equation
    double reduced_mass = (mass_kg * neutral_mass_kg) / (mass_kg + neutral_mass_kg);
    double numerator = 3.0 * std::abs(charge_C) / (16.0 * LOSCHMIDT_CONSTANT);
    double sqrt_term = std::sqrt(2.0 * M_PI / (BOLTZMANN_CONSTANT * temperature_K * reduced_mass));
    
    return numerator * sqrt_term / mobility_m2Vs;
}

void Species::validate() const {
    if (id.empty()) {
        throw std::runtime_error("Species ID cannot be empty");
    }
    if (mass_u <= 0.0) {
        throw std::runtime_error("Species '" + id + "': mass must be positive, got " + std::to_string(mass_u));
    }
    if (std::abs(charge_e) > 10.0) {
        throw std::runtime_error("Species '" + id + "': charge seems unreasonable (|q| > 10), got " + std::to_string(charge_e));
    }
    if (CCS_m2 < 0.0) {
        throw std::runtime_error("Species '" + id + "': CCS cannot be negative, got " + std::to_string(CCS_m2));
    }
    if (mobility_m2Vs < 0.0) {
        throw std::runtime_error("Species '" + id + "': mobility cannot be negative, got " + std::to_string(mobility_m2Vs));
    }
}

// -----------------------------
// SpeciesDatabase member functions
// -----------------------------

void SpeciesDatabase::add(const Species& species) {
    if (has(species.id)) {
        throw std::runtime_error("Species '" + species.id + "' already exists in database");
    }
    
    species.validate();
    
    size_t index = species_list_.size();
    species_list_.push_back(species);
    id_to_index_[species.id] = index;
}

const Species& SpeciesDatabase::get(const std::string& id) const {
    auto it = id_to_index_.find(id);
    if (it == id_to_index_.end()) {
        throw std::runtime_error("Species '" + id + "' not found in database");
    }
    return species_list_[it->second];
}

bool SpeciesDatabase::has(const std::string& id) const {
    return id_to_index_.find(id) != id_to_index_.end();
}

void SpeciesDatabase::calculate_derived_quantities(double temperature_K, double neutral_mass_kg) {
    for (auto& species : species_list_) {
        // Calculate reduced mass
        species.reduced_mass_kg = (species.mass_kg * neutral_mass_kg) / 
                                 (species.mass_kg + neutral_mass_kg);
        
        // Skip neutrals for mobility/CCS calculation
        if (std::abs(species.charge_e) < 0.1) {
            continue;
        }
        
        // Calculate missing mobility or CCS
        if (species.CCS_m2 > 0.0 && species.mobility_m2Vs <= 0.0) {
            species.mobility_m2Vs = species.calculate_mobility(temperature_K, neutral_mass_kg);
        } else if (species.mobility_m2Vs > 0.0 && species.CCS_m2 <= 0.0) {
            species.CCS_m2 = species.calculate_CCS(temperature_K, neutral_mass_kg);
        }
    }
}

// -----------------------------
// Loader functions
// -----------------------------

SpeciesDatabase load_species(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open species file: " + filepath);
    }
    
    Json::Value root;
    try {
        file >> root;
    } catch (const Json::Exception& e) {
        throw std::runtime_error("JSON parse error in " + filepath + ": " + e.what());
    }
    
    if (!root.isMember("species")) {
        throw std::runtime_error("Species file missing 'species' array: " + filepath);
    }
    
    if (!root["species"].isArray()) {
        throw std::runtime_error("'species' must be an array in: " + filepath);
    }
    
    SpeciesDatabase db;
    const Json::Value& species_array = root["species"];
    
    for (Json::ArrayIndex i = 0; i < species_array.size(); ++i) {
        const Json::Value& j = species_array[i];
        
        // Required fields
        if (!j.isMember("id")) {
            throw std::runtime_error("Species entry " + std::to_string(i) + " missing required field 'id'");
        }
        if (!j.isMember("mass_u")) {
            throw std::runtime_error("Species '" + j["id"].asString() + "' missing required field 'mass_u'");
        }
        if (!j.isMember("charge_e")) {
            throw std::runtime_error("Species '" + j["id"].asString() + "' missing required field 'charge_e'");
        }
        
        Species species;
        species.id = j["id"].asString();
        species.name = j.get("name", species.id).asString();
        species.mass_u = j["mass_u"].asDouble();
        species.mass_kg = species.mass_u * AMU_TO_KG;
        species.charge_e = j["charge_e"].asDouble();
        species.charge_C = species.charge_e * ELEM_CHARGE_C;
        
        // Optional fields
        species.CCS_m2 = j.get("CCS_m2", 0.0).asDouble();
        species.mobility_m2Vs = j.get("mobility_m2Vs", 0.0).asDouble();
        
        if (j.isMember("geometry_file") && j["geometry_file"].isString()) {
            species.geometry_file = j["geometry_file"].asString();
        }
        
        // Validate that at least one of CCS or mobility is provided for ions
        if (std::abs(species.charge_e) > 0.1) {
            if (species.CCS_m2 <= 0.0 && species.mobility_m2Vs <= 0.0) {
                throw std::runtime_error("Species '" + species.id + "': must provide either CCS_m2 or mobility_m2Vs");
            }
        }
        
        db.add(species);
    }
    
    ICARION::log::Logger::config()->info("Loaded {} species from {}", db.size(), filepath);
    
    return db;
}

SpeciesDatabase load_species(const std::string& filepath, 
                             double temperature_K, 
                             double neutral_mass_kg) {
    auto db = load_species(filepath);
    db.calculate_derived_quantities(temperature_K, neutral_mass_kg);
    return db;
}

}  // namespace io
}  // namespace ICARION
