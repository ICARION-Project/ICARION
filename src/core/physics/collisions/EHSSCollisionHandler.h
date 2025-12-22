// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file EHSSCollisionHandler.h
 * @brief Explicit Hard-Sphere Scattering (EHSS) collision handler
 * 
 * Implements structure-resolved hard-sphere scattering using atom-centered spheres.
 * Uses molecular geometry (atom positions and radii) for collision detection; falls back
 * to CCS if geometry is missing.
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
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "core/config/types/SpeciesConfig.h"
#include <vector>
#include <unordered_map>
#include <utility>
#include <unordered_set>

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
 * Recommended when geometry is available and accurate scattering is needed; slower than\n+ * HSS due to per-atom checks.
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
     * @param geometry_map Map of species → (atom_centers, atom_radii) [COPIED to handler]
     * @param enable_logging Enable debug logging (writes CSV file with collision details)
     * 
     * @throws std::invalid_argument if geometry_map is empty
     */
    explicit EHSSCollisionHandler(
        GeometryMap geometry_map,
        bool enable_logging = false,
        const config::SpeciesDatabase* species_db = nullptr
    );
    
    /**
     * @brief Handle EHSS collision for single timestep (SoA)
     */
    bool handle_collision(
        core::IonCollisionData& view,
        double dt,
        PhysicsRng& rng,
        const config::EnvironmentConfig& env
    ) override;
    
    std::string name() const override { return "EHSS"; }
    
    CollisionStats get_stats() const override { return stats_; }
    void reset_stats() override { stats_ = {}; }
    
private:
    const GeometryMap geometry_map_;  ///< Copy of geometry map (owned by handler)
    bool enable_logging_;
    const config::SpeciesDatabase* species_db_ = nullptr;
    mutable CollisionStats stats_;
    mutable std::unordered_set<std::string> warned_missing_sigma_;
    
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
        double neutral_radius,
        const std::string& gas_id = ""
    ) const;
    
    /**
     * @brief Derive CCS for target gas from reference gas (HSS formula)
     * 
     * Uses hard-sphere approximation to extrapolate CCS from a known reference:
     * 1. Extract ion radius from reference CCS: r_ion = sqrt(σ_ref/π) - r_ref
     * 2. Compute CCS for target gas: σ_target = π(r_ion + r_target)²
     * 
     * This is the same formula used by tools/ccs_precompute.
     * 
     * @param sigma_ref_m2 Reference CCS [m²]
     * @param gas_ref Reference gas name (e.g., "He")
     * @param gas_target Target gas name (e.g., "N2")
     * @return Derived CCS [m²] or 0.0 if derivation not possible
     */
    double derive_ccs_for_target_gas(
        double sigma_ref_m2,
        const std::string& gas_ref,
        const std::string& gas_target
    ) const;
};

} // namespace ICARION::physics
