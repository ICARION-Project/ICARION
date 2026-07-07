// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file DomainManager.h
 * @brief Manages instrument domains: lookup, coordinate transforms, boundary conditions
 * 
 * Part of Phase 5A refactoring: SimulationEngine
 * 
 * Responsibilities:
 * - Find which domain contains an ion (spatial lookup)
 * - Transform coordinates between global and local frames
 * - Check aperture crossings (domain transitions)
 * - Update ion properties when transitioning domains
 * 
 * Extracts domain management logic from legacy integrate_trajectory().
 * Supports cylindrical and Orbitrap domains; other geometries are not handled.
 */

#pragma once

#include "core/config/types/DomainConfig.h"  // SSOT: config::DomainConfig
#include "core/config/types/IFieldModel.h"   // Field models (analytical/map)
#include "core/config/types/IDomainGeometry.h" // Geometry strategies
#include "core/types/Vec3.h"                // Vec3
#include "core/types/IonState.h"            // IonState
#include "boundaries/BoundaryAction.h"
#include "boundaries/BoundaryActionFactory.h"
#include <vector>
#include <memory>
#include <random>

namespace ICARION {
namespace integrator {

/**
 * @brief Floating-point tolerance for domain boundary checks [m]
 * 
 * Important for IMS/SIFDT: Ions starting at z=0 may have z≈-1e-16 due to FP roundoff.
 * Value: 1e-12 m = 1 picometer (negligible compared to typical domain sizes)
 */
constexpr double DOMAIN_BOUNDARY_EPSILON = 1e-12;

/**
 * @brief Numerical tolerance for ray-tracing calculations
 * 
 * Used to avoid division by zero in intersection computations.
 * Rays with direction magnitude < NUMERICAL_ZERO are treated as degenerate.
 * Value: 1e-15 (≈ sqrt(machine epsilon for double precision))
 */
constexpr double NUMERICAL_ZERO = 1e-15;

/**
 * @class DomainManager
 * @brief Manages spatial domains and coordinate transformations
 * 
 * Provides utilities for:
 * - Domain lookup (which domain contains a given position?)
 * - Coordinate transforms (global ↔ local)
 * - Aperture crossing detection (ion transitions between domains)
 * - Domain property updates (temperature, pressure, gas velocity)
 * 
 * **Floating-Point Tolerance:**
 * Domain boundary checks use epsilon tolerance (1e-12 m = 1 pm) to handle:
 * - Ions starting exactly at z=0 (IMS, SIFDT) may have z≈-1e-16 due to FP roundoff
 * - Integration errors accumulating near boundaries
 * - Coordinate transform precision loss
 * 
 * This prevents false-positive "ion left domain" for ions legitimately inside.
 * 
 * Thread-safe for read-only operations (OpenMP parallel integration).
 */
class DomainManager {
public:
    /**
     * @brief Construct domain manager
     * @param domains Vector of domain configs (reference, not owned)
     * @param rng_seed Random seed for boundary actions (diffuse/thermal reflection)
     * 
     * Stores a const reference to domains (owned by FullConfig).
     * Creates boundary actions for each domain based on config.
     * SSOT: Uses config::DomainConfig (modern format), not legacy InstrumentDomain.
     */
    explicit DomainManager(
        const std::vector<config::DomainConfig>& domains,
        unsigned int rng_seed = 42
    );
    
    /**
     * @brief Find which domain contains the ion position
     * @param pos Global position [m]
     * @return Domain index (0-based), or -1 if outside all domains
     * 
     * Uses internal geometry check for cylindrical/hyperbolic boundaries.
     * Returns index of first matching domain (order matters if domains overlap).
     */
    int find_domain_index(const Vec3& pos) const;
    
    /**
     * @brief Get domain by index
     * @param idx Domain index (must be valid: 0 <= idx < domains.size())
     * @return Const reference to domain (SSOT: config::DomainConfig)
     * @throws std::out_of_range if idx invalid
     */
    const config::DomainConfig& get_domain(int idx) const;
    
    /**
     * @brief Transform position from global to local coordinates
     * @param pos Global position [m]
     * @param domain_idx Domain index
     * @return Local position [m] (origin at domain origin, rotated to domain frame)
     * 
     * Formula: pos_local = R_g2l * (pos_global - origin)
     * where R_g2l is the rotation matrix from global to local frame.
     */
    Vec3 global_to_local_pos(const Vec3& pos, int domain_idx) const;
    
    /**
     * @brief Transform velocity from global to local coordinates
     * @param vel Global velocity [m/s]
     * @param domain_idx Domain index
     * @return Local velocity [m/s] (rotated to domain frame)
     * 
     * Formula: vel_local = R_g2l * vel_global
     */
    Vec3 global_to_local_vel(const Vec3& vel, int domain_idx) const;
    
    /**
     * @brief Transform position from local to global coordinates
     * @param pos_local Local position [m]
     * @param domain_idx Domain index
     * @return Global position [m]
     * 
     * Formula: pos_global = R_l2g * pos_local + origin
     */
    Vec3 local_to_global_pos(const Vec3& pos_local, int domain_idx) const;
    
    /**
     * @brief Transform velocity from local to global coordinates
     * @param vel_local Local velocity [m/s]
     * @param domain_idx Domain index
     * @return Global velocity [m/s]
     * 
     * Formula: vel_global = R_l2g * vel_local
     */
    Vec3 local_to_global_vel(const Vec3& vel_local, int domain_idx) const;
    
    /**
     * @brief Get number of domains
     */
    size_t num_domains() const { return domains_.size(); }

    /**
     * @brief Get field model for a domain (if constructed)
     * @param idx Domain index
     * @return Pointer to field model or nullptr
     */
    const config::IFieldModel* field_model(int idx) const;

    /**
     * @brief Access boundary action for a domain (may be null)
     */
    BoundaryAction* boundary_action(int idx) const;

    /**
     * @brief Compute inward surface normal for a global position on boundary.
     */
    Vec3 surface_normal_global(const Vec3& global_pos, int domain_idx) const;

    /**
     * @brief Compute boundary intersection for a segment in global coordinates.
     *
     * @param start_global Segment start
     * @param end_global Segment end
     * @param domain_idx Domain to test against
     * @param intersection_global Output intersection (global) if hit
     * @param normal_global Output inward normal (global) if hit
     * @return true if intersection detected
     */
    bool boundary_intersection_global(const Vec3& start_global,
                                      const Vec3& end_global,
                                      int domain_idx,
                                      Vec3& intersection_global,
                                      Vec3& normal_global) const;

    /**
     * @brief Preserve the existing axial bridge behavior when a step lands just beyond a domain end.
     */
    int forward_axial_bridge_domain(int current_domain_idx, const Vec3& pos_after) const;

    /**
     * @brief Preserve the existing shared-boundary handoff behavior near a downstream aperture.
     */
    int shared_boundary_handoff_domain(int current_domain_idx,
                                       int found_domain_idx,
                                       const Vec3& pos_after,
                                       double tolerance_m = 1e-9) const;

private:
    const std::vector<config::DomainConfig>& domains_;  ///< Reference to domain list (not owned, SSOT)
    std::vector<std::unique_ptr<BoundaryAction>> boundary_actions_;  ///< Boundary actions per domain
    std::vector<std::unique_ptr<config::IFieldModel>> field_models_; ///< Field models per domain (analytical fallback)
    std::vector<std::unique_ptr<config::IDomainGeometry>> geometries_; ///< Geometry strategy per domain
    std::mt19937 rng_;  ///< Random number generator for boundary actions

    struct DomainSpan {
        double z_min;
        double z_max;
        double radius;
        int idx;
    };
    std::vector<DomainSpan> spans_;  ///< Sorted by z_min for fast axial prefilter (cylindrical domains only)
    bool has_orbitrap_ = false;      ///< If any domain is orbitrap -> fallback to linear scan
    
};

}  // namespace integrator
}  // namespace ICARION
