// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_FULL_CONFIG_H
#define ICARION_CONFIG_FULL_CONFIG_H

#include "SimulationConfig.h"
#include "PhysicsConfig.h"
#include "OutputConfig.h"
#include "DomainConfig.h"
#include <vector>
#include <string>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Complete ICARION configuration
 * 
 * Top-level aggregation of all config sections.
 * This is what ConfigLoader returns after parsing JSON.
 * 
 * Replaces the fragmented GlobalParams + vector<InstrumentDomain> pattern.
 */
struct FullConfig {
    // === Core sections ===
    SimulationConfig simulation;
    PhysicsConfig physics;
    OutputConfig output;
    std::vector<DomainConfig> domains;
    
    // === Database/file paths ===
    std::string species_database_path = "";     ///< Species properties database (includes geometry data)
    std::string reaction_database_path = "";    ///< Reaction rates database
    std::string ion_cloud_path = "";            ///< Initial ion cloud distribution
    
    // === Optional metadata ===
    std::string title = "";                     ///< Simulation title/description
    std::string config_file_path = "";          ///< Path to loaded config file (for reference)
    
    /**
     * @brief Finalize all domains
     * 
     * Calls finalize() on all domains to compute derived quantities.
     * Should be called after loading from JSON.
     */
    void finalize_all() {
        simulation.compute_derived();
        
        for (size_t i = 0; i < domains.size(); ++i) {
            domains[i].domain_index = static_cast<int>(i);
            domains[i].finalize();
        }
    }
    
    /**
     * @brief Validate complete configuration
     * 
     * @throws std::runtime_error if invalid
     */
    void validate() const {
        // Validate each section
        simulation.validate();
        physics.validate();
        output.validate();
        
        if (domains.empty()) {
            throw std::runtime_error("Configuration must have at least one domain");
        }
        
        for (const auto& domain : domains) {
            domain.validate();
        }
        
        // Cross-domain validation
        // Check for duplicate domain names
        for (size_t i = 0; i < domains.size(); ++i) {
            for (size_t j = i + 1; j < domains.size(); ++j) {
                if (domains[i].name == domains[j].name) {
                    throw std::runtime_error("Duplicate domain name: '" + domains[i].name + "'");
                }
            }
        }
        
        // Physics validation
        if (physics.enable_reactions && reaction_database_path.empty()) {
            throw std::runtime_error("Reactions enabled but no reaction database specified");
        }
        
        // Ion cloud validation
        if (ion_cloud_path.empty()) {
            throw std::runtime_error("No ion cloud file specified");
        }
    }
    
    /**
     * @brief Get domain by name
     * 
     * @param name Domain name to search for
     * @return Pointer to domain, or nullptr if not found
     */
    const DomainConfig* get_domain(const std::string& name) const {
        for (const auto& domain : domains) {
            if (domain.name == name) {
                return &domain;
            }
        }
        return nullptr;
    }
    
    /**
     * @brief Get domain by index
     * 
     * @param index Domain index (0-based)
     * @return Reference to domain
     * @throws std::out_of_range if index invalid
     */
    const DomainConfig& get_domain(size_t index) const {
        if (index >= domains.size()) {
            throw std::out_of_range("Domain index out of range: " + std::to_string(index));
        }
        return domains[index];
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_FULL_CONFIG_H
