// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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
            
            // VALIDATION RULE #3: species exists in database (if not "neutral")
            for (const auto& term : rxn.order_terms) {
                if (term.species != "neutral" && !species_db->has(term.species)) {
                    throw std::runtime_error(
                        "Reaction '" + rxn.id + "': order term species '" + term.species + 
                        "' not found in species database"
                    );
                }
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
    
    if (!json.isMember("rate_constant") || !json["rate_constant"].isNumeric()) {
        throw std::runtime_error("Reaction '" + rxn.id + "': missing or invalid 'rate_constant' field");
    }
    rxn.rate_constant = json["rate_constant"].asDouble();
    
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
        std::unordered_map<std::string, int> species_count;  // Track duplicate species
        int neutral_fallback_count = 0;  // Track multiple neutral=-1 entries
        
        for (const auto& term_json : json["order"]) {
            ReactionOrderTerm term = parse_order_term(term_json);
            
            // VALIDATION RULE #4: No duplicate species in order terms
            if (species_count.count(term.species) > 0) {
                throw std::runtime_error(
                    "Reaction '" + rxn.id + "': duplicate order term for species '" + term.species + 
                    "'. Use exponent=2 instead of two separate terms."
                );
            }
            species_count[term.species] = 1;
            
            // VALIDATION RULE #5: Max one neutral fallback (concentration_m3 = -1)
            if (term.concentration_m3 == -1.0) {
                neutral_fallback_count++;
                if (neutral_fallback_count > 1) {
                    throw std::runtime_error(
                        "Reaction '" + rxn.id + "': multiple order terms with concentration_m3 = -1.0 " +
                        "(buffer gas fallback). Only one term can use buffer gas density."
                    );
                }
            }
            
            rxn.order_terms.push_back(term);
        }
        
        // ⚠️ DIMENSIONAL CONSISTENCY WARNING
        int total_order = 0;
        for (const auto& term : rxn.order_terms) {
            total_order += term.exponent;
        }
        
        // Warn if rate_constant dimensions likely mismatch order
        if (total_order == 0 && rxn.rate_constant > 1e-6) {
            std::cout << "⚠  Reaction '" << rxn.id << "': spontaneous decay (order=0) but k = " 
                      << rxn.rate_constant << " suggests [m³/s] units. Expected [s⁻¹] (~1e-3 to 1e6).\n";
        }
        
        if (total_order == 1 && (rxn.rate_constant < 1e-12 || rxn.rate_constant > 1e-6)) {
            std::cout << "⚠  Reaction '" << rxn.id << "': 2nd-order (total exponent=1) but k = " 
                      << rxn.rate_constant << " outside typical range [1e-12, 1e-6] m³/s.\n";
        }
        
        if (total_order == 2 && (rxn.rate_constant < 1e-30 || rxn.rate_constant > 1e-24)) {
            std::cout << "⚠  Reaction '" << rxn.id << "': 3rd-order (total exponent=2) but k = " 
                      << rxn.rate_constant << " outside typical range [1e-30, 1e-24] m⁶/s.\n";
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
    
    // VALIDATION RULE #1: exponent ∈ {0, 1, 2}
    if (term.exponent < 0 || term.exponent > 2) {
        throw std::runtime_error(
            "Order term for '" + term.species + "': exponent must be 0, 1, or 2 (got " + 
            std::to_string(term.exponent) + ")"
        );
    }
    
    // concentration_m3 is optional (default: -1.0 = use buffer gas)
    if (json.isMember("concentration_m3") && json["concentration_m3"].isNumeric()) {
        term.concentration_m3 = json["concentration_m3"].asDouble();
    } else {
        term.concentration_m3 = -1.0;  // Default: use buffer gas density
    }
    
    // VALIDATION RULE #2: concentration_m3 ≥ -1.0
    if (term.concentration_m3 < -1.0) {
        throw std::runtime_error(
            "Order term for '" + term.species + "': concentration_m3 must be ≥ -1.0 (got " + 
            std::to_string(term.concentration_m3) + "). Use -1.0 for buffer gas fallback."
        );
    }
    
    // Reject ambiguous negative values (e.g., -0.5)
    if (term.concentration_m3 > -1.0 && term.concentration_m3 < 0.0) {
        throw std::runtime_error(
            "Order term for '" + term.species + "': concentration_m3 cannot be between -1.0 and 0.0 (got " + 
            std::to_string(term.concentration_m3) + "). Use -1.0 for buffer gas or ≥ 0 for explicit value."
        );
    }
    
    return term;
}

} // namespace ICARION::config
