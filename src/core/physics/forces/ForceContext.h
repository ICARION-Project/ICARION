// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <vector>
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

// Forward declaration (IFieldProvider is in global namespace)
class IFieldProvider;

// Forward declare DomainConfig (in correct namespace)
namespace ICARION::config {
    struct DomainConfig;
}

namespace ICARION::physics {

/**
 * @brief Context for force computation
 * 
 * Contains shared data needed by multiple forces. This struct avoids parameter
 * explosion by grouping related data together.
 * 
 * Design Notes:
 * - Uses pointers for optional data (nullptr = not needed)
 * - No ownership (context is transient, data owned elsewhere)
 * - Lightweight (only pointers, cheap to copy)
 * 
 * Usage Pattern:
 * @code
 * ForceContext ctx;
 * ctx.field_provider = &my_field_provider;
 * ctx.domain = &domain_config;
 * ctx.all_ions = &ion_ensemble;
 * 
 * Vec3 force = force_registry.compute_total_force(ion, t, ctx);
 * @endcode
 */
struct ForceContext {
    // =========================================================================
    // Field Provider (for grid-based or analytical field sampling)
    // =========================================================================
    
    /**
     * @brief Optional field provider for E-field and B-field sampling
     * 
     * If non-null, forces can use this to sample fields at arbitrary positions.
     * If null, forces fall back to analytical field calculations.
     * 
     * @see IFieldProvider for field sampling interface
     */
    const ::IFieldProvider* field_provider = nullptr;
    
    // =========================================================================
    // Domain Configuration (geometry, instrument type, field parameters)
    // =========================================================================
    
    /**
     * @brief Optional domain configuration
     * 
     * Contains instrument type, geometry, field parameters, environment properties.
     * Forces can use this for instrument-specific calculations.
     * 
     * @see DomainConfig for complete structure
     */
    const ::ICARION::config::DomainConfig* domain = nullptr;
    
    // =========================================================================
    // Ion Ensemble (for N-body forces like space charge)
    // =========================================================================
    
    /**
     * @brief Optional ion ensemble for N-body force calculations
     * 
     * Required for space charge force (ion-ion interactions).
     * Can be null for single-particle forces (electric, magnetic, damping).
     * 
     * @note For space charge: exclude self-interaction (ion == current ion)
     */
    const std::vector<IonState>* all_ions = nullptr;
    
    // =========================================================================
    // SSOT Compliance Note
    // =========================================================================
    
    /**
     * **IMPORTANT: No duplicate environment parameters!**
     * 
     * All environment properties (temperature, pressure, gas density, etc.)
     * are accessed via `domain->environment`, NOT stored here.
     * 
     * **SSOT Pattern:**
     * ```cpp
     * // ✅ CORRECT: Read from domain
     * double T = ctx.domain->environment.temperature_K;
     * double P = ctx.domain->environment.pressure_Pa;
     * double n = ctx.domain->environment.particle_density_m_3;
     * 
     * // ❌ WRONG: Don't duplicate data in ForceContext!
     * // double temperature_K;  // SSOT violation!
     * ```
     * 
     * **Why?**
     * - Single Source of Truth: Config changes propagate automatically
     * - No synchronization bugs: Can't have out-of-sync values
     * - Smaller context: Less memory, faster copying
     * - Type safety: Compiler enforces correct config usage
     * 
     * **All environment data available via:**
     * - `ctx.domain->environment.temperature_K`
     * - `ctx.domain->environment.pressure_Pa`
     * - `ctx.domain->environment.particle_density_m_3`
     * - `ctx.domain->environment.gas_velocity_m_s` (Vec3)
     * - `ctx.domain->environment.neutral_mass_amu`
     * - `ctx.domain->environment.polarizability_m3`
     * - `ctx.domain->environment.collision_model`
     * 
     * @see EnvironmentConfig for complete structure
     */
};

} // namespace ICARION::physics
