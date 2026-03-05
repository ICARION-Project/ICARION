// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

#include <cstddef>
#include <vector>
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

// Forward declaration (IFieldProvider is in global namespace)
class IFieldProvider;

// Forward declare DomainConfig (in correct namespace)
namespace ICARION::config {
    struct DomainConfig;
    class IFieldModel;
}

// Forward declare IonEnsemble (SoA container)
namespace ICARION::core {
    class IonEnsemble;
}

namespace ICARION::physics {

/**
 * @brief Context for force computation (no ownership)
 * 
 * Holds optional pointers to shared data used by force implementations. Callers fill
 * what they need; nullptr means “not available”. Environment/geometry lives in the
 * DomainConfig referenced here.
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
    
    /**
     * @brief Optional analytical/map field model (new domain/field layer)
     *
     * If set, forces should prefer `field_model->E/B(...)` as SSOT. Typically
     * injected by PhysicsSetup (analytical vs. grid) and forwarded by
     * integrators; SimulationEngine only backfills from DomainManager as a
     * legacy fallback.
     */
    const ::ICARION::config::IFieldModel* field_model = nullptr;
    
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

    /**
     * @brief Optional SoA ensemble view for N-body forces in SoA path
     *
     * Used by space-charge forces to avoid AoS reconstruction in SoA integrators.
     * @code
     * ctx.ion_ensemble = &ensemble;
     * ctx.ion_index    = i;
     * @endcode
     */
    const ::ICARION::core::IonEnsemble* ion_ensemble = nullptr;

    /**
     * @brief Index of the current ion in ion_ensemble (SoA)
     *
     * Use together with `ion_ensemble` to identify the active ion in
     * space-charge computations. Ignored if `ion_ensemble` is null.
     */
    size_t ion_index = static_cast<size_t>(-1);
    
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
     * // CORRECT: Read from domain
     * double T = ctx.domain->environment.temperature_K;
     * double P = ctx.domain->environment.pressure_Pa;
     * double n = ctx.domain->environment.particle_density_m_3;
     * 
     * // WRONG: Don't duplicate data in ForceContext!
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
