// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "core/types/IonEnsemble.h"  // For IonReactionData view
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/CollisionTypes.h"  // PhysicsRng
#include <string>

namespace ICARION {
namespace physics {

/**
 * @brief Reaction statistics (for tracking/logging)
 */
struct ReactionStats {
    size_t total_reactions = 0;
    size_t failed_lookups = 0;  // Species not found in database
    double average_reaction_rate = 0.0;
};

/**
 * @brief Reaction handler interface
 * 
 * All reaction models implement this interface.
 * 
 * **SSOT Design:** Handler reads directly from databases and EnvironmentConfig.
 * No intermediate parameter structs (avoids temperature/density duplication).
 * 
 * **Design Pattern:** Strategy pattern for reaction handling
 * - StochasticReactionHandler: Monte Carlo discrete reactions
 * - NoReactionHandler: Null object (reactions disabled)
 * 
 * **Example Usage:**
 * @code
 * auto handler = ReactionHandlerFactory::create(config.physics);
 * 
 * auto view = ensemble.reaction_data(i);  // SoA view
 * bool reaction_occurred = handler->handle_reaction(
 *     view, dt, rng,
 *     config.reaction_db,    // Direct reference (SSOT!)
 *     config.species_db,     // Direct reference (SSOT!)
 *     domain.environment     // Contains temperature_K, particle_density_m_3
 * );
 * @endcode
 * 
 * @see ReactionHandlerFactory
 * @see StochasticReactionHandler
 */
class IReactionHandler {
public:
    virtual ~IReactionHandler() = default;
    
    /**
     * @brief Handle reactions using SoA view (cache-optimized hot path)
     * 
     * @param[in,out] view Ion reaction data view (species may be modified)
     * @param[in,out] cold_data Cold data arrays (CCS, mobility) for species updates
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] reaction_db Reaction database
     * @param[in] species_db Species database
     * @param[in] env Environment configuration
     * 
     * @return true if reaction occurred, false otherwise
     */
    virtual bool handle_reaction(
        core::IonReactionData& view,
        double dt,
        PhysicsRng& rng,
        const config::ReactionDatabase& reaction_db,
        const config::SpeciesDatabase& species_db,
        const config::EnvironmentConfig& env
    ) = 0;
    
    /**
     * @brief Get handler name
     * @return Human-readable handler name (e.g., \"Stochastic\", \"None\")
     */
    virtual std::string name() const = 0;
    
    /**
     * @brief Get reaction statistics
     * @return Statistics struct (total reactions, failed lookups, etc.)
     */
    virtual ReactionStats get_stats() const { return {}; }
    
    /**
     * @brief Reset statistics
     */
    virtual void reset_stats() {}
};

} // namespace physics
} // namespace ICARION
