// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "core/io/reactionLoader.h"
#include "core/log/Logger.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <json/json.h>

namespace ICARION {
namespace io {

// -----------------------------
// Helper functions
// -----------------------------

ReactionType parse_reaction_type(const std::string& type_str) {
    std::string lower = type_str;
    std::transform(lower.begin(), lower.end(), lower.begin(), 
                  [](unsigned char c){ return std::tolower(c); });
    
    if (lower == "charge_transfer" || lower == "chargetransfer") 
        return ReactionType::ChargeTransfer;
    if (lower == "proton_transfer" || lower == "protontransfer") 
        return ReactionType::ProtonTransfer;
    if (lower == "association") 
        return ReactionType::Association;
    if (lower == "dissociation") 
        return ReactionType::Dissociation;
    if (lower == "switching") 
        return ReactionType::Switching;
    
    return ReactionType::Unknown;
}

std::string reaction_type_to_string(ReactionType type) {
    switch (type) {
        case ReactionType::ChargeTransfer: return "ChargeTransfer";
        case ReactionType::ProtonTransfer: return "ProtonTransfer";
        case ReactionType::Association: return "Association";
        case ReactionType::Dissociation: return "Dissociation";
        case ReactionType::Switching: return "Switching";
        default: return "Unknown";
    }
}

// -----------------------------
// Reaction member functions
// -----------------------------

double Reaction::effective_rate_constant(double temperature_K) const {
    if (activation_energy_eV <= 0.0) {
        return rate_constant_SI;
    }
    
    // Arrhenius equation: k(T) = k0 * exp(-Ea / (kB * T))
    double Ea_J = activation_energy_eV * ELEM_CHARGE_C;  // Convert eV to J
    double exponent = -Ea_J / (BOLTZMANN_CONSTANT * temperature_K);
    return rate_constant_SI * std::exp(exponent);
}

void Reaction::validate(const SpeciesDatabase& species_db) const {
    if (id.empty()) {
        throw std::runtime_error("Reaction ID cannot be empty");
    }
    
    if (reactants.empty()) {
        throw std::runtime_error("Reaction '" + id + "': must have at least one reactant");
    }
    
    if (products.empty()) {
        throw std::runtime_error("Reaction '" + id + "': must have at least one product");
    }
    
    if (rate_constant_SI <= 0.0) {
        throw std::runtime_error("Reaction '" + id + "': rate constant must be positive, got " + 
                                std::to_string(rate_constant_SI));
    }
    
    if (branching_ratio < 0.0 || branching_ratio > 1.0) {
        throw std::runtime_error("Reaction '" + id + "': branching ratio must be between 0 and 1, got " + 
                                std::to_string(branching_ratio));
    }
    
    // Validate all reactant species exist
    for (const auto& species_id : reactants) {
        if (!species_db.has(species_id)) {
            throw std::runtime_error("Reaction '" + id + "': reactant species '" + 
                                    species_id + "' not found in species database");
        }
    }
    
    // Validate all product species exist
    for (const auto& species_id : products) {
        if (!species_db.has(species_id)) {
            throw std::runtime_error("Reaction '" + id + "': product species '" + 
                                    species_id + "' not found in species database");
        }
    }
    
    // Validate order term species exist
    for (const auto& term : order_terms) {
        if (!species_db.has(term.species_id)) {
            throw std::runtime_error("Reaction '" + id + "': order term species '" + 
                                    term.species_id + "' not found in species database");
        }
    }
}

int Reaction::total_order() const {
    int order = 0;
    for (const auto& term : order_terms) {
        order += term.exponent;
    }
    return order;
}

// -----------------------------
// ReactionDatabase member functions
// -----------------------------

void ReactionDatabase::add(const Reaction& reaction) {
    if (has(reaction.id)) {
        throw std::runtime_error("Reaction '" + reaction.id + "' already exists in database");
    }
    
    size_t index = reactions_list_.size();
    reactions_list_.push_back(reaction);
    id_to_index_[reaction.id] = index;
    update_species_index(index);
}

const Reaction& ReactionDatabase::get(const std::string& id) const {
    auto it = id_to_index_.find(id);
    if (it == id_to_index_.end()) {
        throw std::runtime_error("Reaction '" + id + "' not found in database");
    }
    return reactions_list_[it->second];
}

bool ReactionDatabase::has(const std::string& id) const {
    return id_to_index_.find(id) != id_to_index_.end();
}

std::vector<const Reaction*> ReactionDatabase::get_reactions_for_species(const std::string& species_id) const {
    std::vector<const Reaction*> result;
    
    auto it = species_to_reactions_.find(species_id);
    if (it != species_to_reactions_.end()) {
        result.reserve(it->second.size());
        for (size_t idx : it->second) {
            result.push_back(&reactions_list_[idx]);
        }
    }
    
    return result;
}

void ReactionDatabase::validate_all(const SpeciesDatabase& species_db) const {
    for (const auto& reaction : reactions_list_) {
        reaction.validate(species_db);
    }
}

void ReactionDatabase::update_species_index(size_t reaction_idx) {
    const Reaction& rxn = reactions_list_[reaction_idx];
    
    // Add reaction to index for each reactant species
    for (const auto& species_id : rxn.reactants) {
        species_to_reactions_[species_id].push_back(reaction_idx);
    }
}

// -----------------------------
// Loader functions
// -----------------------------

ReactionDatabase load_reactions(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file) {
        throw std::runtime_error("Cannot open reactions file: " + filepath);
    }
    
    Json::Value root;
    try {
        file >> root;
    } catch (const Json::Exception& e) {
        throw std::runtime_error("JSON parse error in " + filepath + ": " + e.what());
    }
    
    if (!root.isMember("reactions")) {
        throw std::runtime_error("Reactions file missing 'reactions' array: " + filepath);
    }
    
    if (!root["reactions"].isArray()) {
        throw std::runtime_error("'reactions' must be an array in: " + filepath);
    }
    
    ReactionDatabase db;
    const Json::Value& reactions_array = root["reactions"];
    
    for (Json::ArrayIndex i = 0; i < reactions_array.size(); ++i) {
        const Json::Value& j = reactions_array[i];
        
        // Required fields
        if (!j.isMember("id")) {
            throw std::runtime_error("Reaction entry " + std::to_string(i) + " missing required field 'id'");
        }
        if (!j.isMember("reactants")) {
            throw std::runtime_error("Reaction '" + j["id"].asString() + "' missing required field 'reactants'");
        }
        if (!j.isMember("products")) {
            throw std::runtime_error("Reaction '" + j["id"].asString() + "' missing required field 'products'");
        }
        if (!j.isMember("rate_constant_SI")) {
            throw std::runtime_error("Reaction '" + j["id"].asString() + "' missing required field 'rate_constant_SI'");
        }
        
        Reaction reaction;
        reaction.id = j["id"].asString();
        reaction.rate_constant_SI = j["rate_constant_SI"].asDouble();
        reaction.activation_energy_eV = j.get("activation_energy_eV", 0.0).asDouble();
        reaction.branching_ratio = j.get("branching_ratio", 1.0).asDouble();
        
        // Parse reaction type
        std::string type_str = "unknown";
        if (j.isMember("type") && j["type"].isString()) {
            type_str = j["type"].asString();
        }
        reaction.type = parse_reaction_type(type_str);
        
        // Parse reactants array
        const Json::Value& reactants = j["reactants"];
        if (!reactants.isArray()) {
            throw std::runtime_error("Reaction '" + reaction.id + "': 'reactants' must be an array");
        }
        for (const auto& r : reactants) {
            reaction.reactants.push_back(r.asString());
        }
        
        // Parse products array
        const Json::Value& products = j["products"];
        if (!products.isArray()) {
            throw std::runtime_error("Reaction '" + reaction.id + "': 'products' must be an array");
        }
        for (const auto& p : products) {
            reaction.products.push_back(p.asString());
        }
        
        // Parse optional order terms
        if (j.isMember("order") && j["order"].isArray()) {
            const Json::Value& order = j["order"];
            for (const auto& term : order) {
                if (!term.isMember("species") || !term.isMember("exponent")) {
                    throw std::runtime_error("Reaction '" + reaction.id + 
                                            "': order term must have 'species' and 'exponent'");
                }
                ReactionOrderTerm rot;
                rot.species_id = term["species"].asString();
                rot.exponent = term["exponent"].asInt();
                reaction.order_terms.push_back(rot);
            }
        }
        
        db.add(reaction);
    }
    
    ICARION::log::Logger::config()->info("Loaded {} reactions from {}", db.size(), filepath);
    
    return db;
}

ReactionDatabase load_reactions(const std::string& filepath, 
                                const SpeciesDatabase& species_db) {
    auto db = load_reactions(filepath);
    db.validate_all(species_db);
    return db;
}

}  // namespace io
}  // namespace ICARION
