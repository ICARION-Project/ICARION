// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "ReactionLoader.h"
#include <fstream>
#include <iostream>

namespace ICARION::config {

ReactionDatabase ReactionLoader::load(const std::filesystem::path& filepath,
                                     const SpeciesDatabase* species_db) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open reaction database: " + filepath.string());
    }
    
    Json::Value root;
    Json::CharReaderBuilder builder;
    std::string errs;
    if (!Json::parseFromStream(builder, file, &root, &errs)) {
        throw std::runtime_error("JSON parse error in " + filepath.string() + ": " + errs);
    }
    
    return load_from_json(root, species_db);
}

ReactionDatabase ReactionLoader::load_from_json(const Json::Value& json,
                                                const SpeciesDatabase* species_db) {
    ReactionDatabase db;
    
    // Expect "reactions" array
    if (!json.isMember("reactions") || !json["reactions"].isArray()) {
        throw std::runtime_error("Reaction database must have 'reactions' array");
    }
    
    const Json::Value& reactions_array = json["reactions"];
    
    for (const auto& rxn_json : reactions_array) {
        Reaction rxn = parse_reaction(rxn_json);
        
        // Validate against species database if provided
        if (species_db) {
            auto validation = rxn.validate(*species_db);
            if (!validation.valid) {
                std::string error_msg = "Reaction '" + rxn.id + "' validation failed:\n";
                for (const auto& err : validation.errors) {
                    error_msg += "  - " + err + "\n";
                }
                throw std::runtime_error(error_msg);
            }
            
            // Print warnings
            for (const auto& warning : validation.warnings) {
                std::cout << "⚠  Reaction '" << rxn.id << "': " << warning << "\n";
            }
        }
        
        db.reactions.push_back(rxn);
    }
    
    // Validate database
    if (species_db) {
        auto db_validation = db.validate(*species_db);
        if (!db_validation.valid) {
            std::string error_msg = "Reaction database validation failed:\n";
            for (const auto& err : db_validation.errors) {
                error_msg += "  - " + err + "\n";
            }
            throw std::runtime_error(error_msg);
        }
    }
    
    std::cout << "✓ Loaded " << db.reactions.size() << " reactions from database\n";
    
    return db;
}

Reaction ReactionLoader::parse_reaction(const Json::Value& json) {
    Reaction rxn;
    
    // Required fields
    if (!json.isMember("id") || !json["id"].isString()) {
        throw std::runtime_error("Reaction missing 'id' field");
    }
    rxn.id = json["id"].asString();
    
    if (!json.isMember("reactant") || !json["reactant"].isString()) {
        throw std::runtime_error("Reaction '" + rxn.id + "': missing 'reactant' field");
    }
    rxn.reactant = json["reactant"].asString();
    
    if (!json.isMember("product") || !json["product"].isString()) {
        throw std::runtime_error("Reaction '" + rxn.id + "': missing 'product' field");
    }
    rxn.product = json["product"].asString();
    
    if (!json.isMember("rate_constant_m3s") || !json["rate_constant_m3s"].isNumeric()) {
        throw std::runtime_error("Reaction '" + rxn.id + "': missing or invalid 'rate_constant_m3s' field");
    }
    rxn.rate_constant_m3s = json["rate_constant_m3s"].asDouble();
    
    // === Optional: Temperature dependence ===
    
    // Rate model (default: Constant)
    if (json.isMember("rate_model") && json["rate_model"].isString()) {
        std::string model_str = json["rate_model"].asString();
        
        if (model_str == "Constant") {
            rxn.rate_model = RateModel::Constant;
        } else if (model_str == "Arrhenius") {
            rxn.rate_model = RateModel::Arrhenius;
        } else if (model_str == "ModifiedArrhenius") {
            rxn.rate_model = RateModel::ModifiedArrhenius;
        } else {
            throw std::runtime_error(
                "Reaction '" + rxn.id + "': unknown rate_model '" + model_str + 
                "' (allowed: Constant, Arrhenius, ModifiedArrhenius)"
            );
        }
    } else {
        rxn.rate_model = RateModel::Constant;  // Default
    }
    
    // Activation energy (optional)
    if (json.isMember("activation_energy_eV") && json["activation_energy_eV"].isNumeric()) {
        rxn.activation_energy_eV = json["activation_energy_eV"].asDouble();
        
        if (rxn.activation_energy_eV < 0.0) {
            throw std::runtime_error(
                "Reaction '" + rxn.id + "': activation_energy_eV must be non-negative"
            );
        }
    }
    
    // Temperature exponent (optional, modified Arrhenius only)
    if (json.isMember("temperature_exponent") && json["temperature_exponent"].isNumeric()) {
        rxn.temperature_exponent = json["temperature_exponent"].asDouble();
    }
    
    // Reference temperature (optional, modified Arrhenius only)
    if (json.isMember("reference_temperature_K") && json["reference_temperature_K"].isNumeric()) {
        rxn.reference_temperature_K = json["reference_temperature_K"].asDouble();
        
        if (rxn.reference_temperature_K <= 0.0) {
            throw std::runtime_error(
                "Reaction '" + rxn.id + "': reference_temperature_K must be positive"
            );
        }
    }
    
    // === Optional: Concentration dependence (order terms) ===
    if (json.isMember("order") && json["order"].isArray()) {
        for (const auto& term_json : json["order"]) {
            rxn.order_terms.push_back(parse_order_term(term_json));
        }
    }
    
    return rxn;
}

ReactionOrderTerm ReactionLoader::parse_order_term(const Json::Value& json) {
    ReactionOrderTerm term;
    
    if (!json.isMember("species") || !json["species"].isString()) {
        throw std::runtime_error("Order term missing 'species' field");
    }
    term.species = json["species"].asString();
    
    if (!json.isMember("exponent") || !json["exponent"].isInt()) {
        throw std::runtime_error("Order term for '" + term.species + "': missing or invalid 'exponent' field");
    }
    term.exponent = json["exponent"].asInt();
    
    // concentration_m3 is optional (can default to 0 and be overridden at runtime)
    if (json.isMember("concentration_m3") && json["concentration_m3"].isNumeric()) {
        term.concentration_m3 = json["concentration_m3"].asDouble();
    } else {
        term.concentration_m3 = 0.0;
    }
    
    return term;
}

} // namespace ICARION::config
