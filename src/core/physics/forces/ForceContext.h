// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#pragma once

#include <vector>
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

namespace ICARION::physics {

// Forward declarations
class IFieldProvider;

namespace config {
    struct DomainConfig;
}

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
    const IFieldProvider* field_provider = nullptr;
    
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
    const config::DomainConfig* domain = nullptr;
    
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
    // Environment Properties (for damping forces)
    // =========================================================================
    
    /**
     * @brief Gas temperature [K]
     * 
     * Used by Langevin and friction forces for thermal velocity calculations.
     * Default: 300 K (room temperature)
     */
    double temperature_K = 300.0;
    
    /**
     * @brief Gas pressure [Pa]
     * 
     * Used by collision-based damping forces.
     * Default: 101325 Pa (1 atm)
     */
    double pressure_Pa = 101325.0;
    
    /**
     * @brief Gas flow velocity [m/s]
     * 
     * Used for drag force in IMS and drift tube instruments.
     * Default: {0, 0, 0} (no flow)
     */
    Vec3 gas_velocity_ms{0, 0, 0};
    
    /**
     * @brief Neutral particle density [m^-3]
     * 
     * Computed from pressure and temperature via ideal gas law.
     * Used by collision-based forces.
     * 
     * @note Can be computed as: n = P / (k_B * T)
     */
    double particle_density_m3 = 0.0;
    
    /**
     * @brief Neutral particle mass [kg]
     * 
     * Mass of background gas molecules (e.g., He, N2).
     * Used for reduced mass calculations in collision forces.
     */
    double neutral_mass_kg = 0.0;
};

} // namespace ICARION::physics
