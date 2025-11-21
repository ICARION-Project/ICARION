// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file EHSSCollisionHandler.h
 * @brief Explicit Hard-Sphere Scattering (EHSS) collision handler
 * 
 * Implements structure-resolved hard-sphere scattering using atom-centered spheres.
 * Uses molecular geometry (atom positions and radii) for accurate collision cross-sections.
 * 
 * **Physics:**
 * - Samples impact parameter and molecular orientation randomly
 * - Detects collision with first atom encountered along trajectory
 * - Applies elastic hard-sphere momentum transfer in center-of-mass frame
 * - Conserves momentum and energy
 * 
 * **SSOT Design:**
 * - Reads gas properties directly from `EnvironmentConfig` (temperature, pressure, etc.)
 * - Stores reference to geometry map (no copies!)
 * - No intermediate parameter structs
 * 
 * @date 2025-11-21
 * @version 1.0
 */

#pragma once

#include "ICollisionHandler.h"
#include "core/types/Vec3.h"
#include <vector>
#include <unordered_map>
#include <utility>

namespace ICARION::physics {

/**
 * @brief Geometry data: atom centers and radii for a molecular species
 * 
 * First element: atom positions [m] in molecule frame
 * Second element: atomic radii [m] for hard-sphere collision detection
 */
using GeometryData = std::pair<std::vector<Vec3>, std::vector<double>>;

/**
 * @brief Map of species ID to molecular geometry
 * 
 * Key: Species identifier (e.g., "H3O+", "N2")
 * Value: Geometry data (atom centers, radii)
 */
using GeometryMap = std::unordered_map<std::string, GeometryData>;

/**
 * @brief EHSS (Explicit Hard-Sphere Scattering) collision handler
 * 
 * Implements structure-resolved collision model using atom-centered hard spheres.
 * Provides accurate collision cross-sections for complex molecular geometries.
 * 
 * **Use cases:**
 * - Polyatomic ions with non-spherical geometries
 * - Accurate mobility calculations
 * - Cross-section validation studies
 * 
 * **Performance:**
 * - Slower than HSS (isotropic) due to geometry sampling
 * - O(N_atoms) collision detection per event
 * - Recommended for N_atoms < 50
 * 
 * **SSOT Pattern:**
 * ```cpp
 * // Correct: Store reference to geometry map
 * GeometryMap geometry;
 * geometry["H3O+"] = load_geometry("h3o+.xyz");
 * EHSSCollisionHandler handler(geometry);
 * 
 * // Handler reads from EnvironmentConfig directly
 * config::EnvironmentConfig env;
 * env.temperature_K = 300.0;
 * handler.handle_collision(ion, dt, rng, env);  // SSOT!
 * ```
 * 
 * @see HSSCollisionHandler for faster isotropic model
 * @see DampingForce for deterministic collision models
 */
class EHSSCollisionHandler : public ICollisionHandler {
public:
    /**
     * @brief Construct EHSS handler
     * 
     * @param geometry_map Map of species → (atom_centers, atom_radii) [REFERENCE stored, not copied!]
     * @param enable_logging Enable debug logging (writes CSV file with collision details)
     * 
     * ⚠️ **IMPORTANT:** `geometry_map` must outlive this handler instance (stores reference!)
     * 
     * @throws std::invalid_argument if geometry_map is empty
     */
    explicit EHSSCollisionHandler(
        const GeometryMap& geometry_map,
        bool enable_logging = false
    );
    
    /**
     * @brief Handle EHSS collision for single timestep
     * 
     * **Algorithm:**
     * 1. Look up geometry for ion species (fallback to ion.CCS_m2 if not found)
     * 2. Compute collision probability from mean free path
     * 3. If collision occurs:
     *    - Sample neutral velocity from Maxwell-Boltzmann distribution
     *    - Sample random molecular orientation
     *    - Detect collision with atom-centered spheres
     *    - Apply hard-sphere momentum transfer in COM frame
     * 
     * **SSOT:** Reads gas properties directly from `env`:
     * - env.temperature_K → thermal velocity
     * - env.pressure_Pa → particle density
     * - env.neutral_mass_kg → reduced mass
     * - env.gas_velocity_m_s → bulk flow
     * - env.neutral_radius_m → fallback CCS if no geometry
     * 
     * @param[in,out] ion Ion state (velocity modified if collision occurs)
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] env Environment configuration (SSOT!)
     * 
     * @return true if collision occurred, false otherwise
     */
    bool handle_collision(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::EnvironmentConfig& env
    ) override;
    
    std::string name() const override { return "EHSS"; }
    
    CollisionStats get_stats() const override { return stats_; }
    void reset_stats() override { stats_ = {}; }
    
private:
    // Reference to geometry map (no copy!)
    const GeometryMap& geometry_map_;
    
    // Debug logging
    bool enable_logging_;
    
    // Statistics
    mutable CollisionStats stats_;
    
    /**
     * @brief Compute effective collision cross-section from geometry
     * 
     * @param ion Ion state (contains species_id)
     * @param neutral_radius Neutral molecule radius [m] (from env.neutral_radius_m)
     * 
     * @return Effective CCS [m²]
     * 
     * Tries to find geometry for species. If found, computes CCS from atom-centered
     * spheres. Otherwise falls back to ion.CCS_m2 (stored effective cross-section).
     */
    double compute_effective_ccs(
        const IonState& ion,
        double neutral_radius
    ) const;
};

} // namespace ICARION::physics
