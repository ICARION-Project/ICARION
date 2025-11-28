// IReactionHandler.h
// Interface for reaction handlers in ICARION particle simulation framework
//
// SSOT Design: Handler reads directly from ReactionDatabase, SpeciesDatabase, and EnvironmentConfig.
// No intermediate parameter structs (no ReactionContext!) to avoid parameter duplication.
//
// Created: 2025-11-22 (Phase 3 Refactor)

#pragma once

#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"  // For IonReactionData view
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/physics/collisions/collisionHelpers.h"  // EhssRng
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
 * bool reaction_occurred = handler->handle_reaction(
 *     ion, dt, rng,
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
     * @brief Handle reactions for single timestep
     * 
     * @param ion Ion state (species_id/mass/charge modified if reaction occurs)
     * @param dt Timestep [s]
     * @param rng Random number generator
     * @param reaction_db Reaction database (SSOT - direct reference!)
     * @param species_db Species database (SSOT - direct reference!)
     * @param env Environment config (contains temperature, density, etc.)
     * 
     * @return true if reaction occurred, false otherwise
     * 
     * **SSOT Compliance:** All parameters read directly from databases/config:
     * - Reaction rates: `reaction_db.reactions[i].rate_constant`
     * - Species properties: `species_db.species[id]`
     * - Temperature: `env.temperature_K`
     * - Density: `env.particle_density_m_3`
     * 
     * **No parameter copies, no intermediate structs!**
     * 
     * **Thread Safety:** Not thread-safe (designed for single-threaded integration loop).
     * For parallel execution, create one handler per thread.
     */
    virtual bool handle_reaction(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::ReactionDatabase& reaction_db,
        const config::SpeciesDatabase& species_db,
        const config::EnvironmentConfig& env
    ) = 0;
    
    /**
     * @brief Handle reactions using SoA view (Phase 3 - cache-optimized)
     * 
     * Zero-copy access to ion data via view struct.
     * Default implementation converts to IonState and calls handle_reaction().
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
     * 
     * @note Override for optimal SoA performance. Default wrapper provided for compatibility.
     */
    virtual bool handle_reaction_soa(
        core::IonReactionData& view,
        double* CCS_array,
        double* mobility_array,
        double dt,
        EhssRng& rng,
        const config::ReactionDatabase& reaction_db,
        const config::SpeciesDatabase& species_db,
        const config::EnvironmentConfig& env
    ) {
        // Default: convert to IonState and call legacy method
        IonState ion;
        ion.pos = view.kin.pos();
        ion.vel = view.kin.vel();
        ion.mass_kg = view.kin.get_mass();
        ion.ion_charge_C = view.kin.get_charge();
        ion.species_id = view.species_id();
        ion.CCS_m2 = CCS_array[view.kin.index];
        ion.reduced_mobility_cm2_Vs = mobility_array[view.kin.index];
        
        bool result = handle_reaction(ion, dt, rng, reaction_db, species_db, env);
        
        // Write back modified data
        view.kin.set_vel(ion.vel);
        view.kin.set_pos(ion.pos);
        // Note: species_id update requires ensemble method (handled by caller)
        CCS_array[view.kin.index] = ion.CCS_m2;
        mobility_array[view.kin.index] = ion.reduced_mobility_cm2_Vs;
        
        return result;
    }
    
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
