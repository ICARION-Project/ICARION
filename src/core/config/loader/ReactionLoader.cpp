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
    
    // Optional: order terms
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
