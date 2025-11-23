// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

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
 * Extracts domain management logic from legacy integrate_trajectory()
 */

#pragma once

#include "core/config/types/DomainConfig.h"  // SSOT: config::DomainConfig
#include "core/types/Vec3.h"                // Vec3
#include "core/types/IonState.h"            // IonState
#include <vector>

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
     * 
     * Stores a const reference to domains (owned by FullConfig).
     * All methods are const (no state modification).
     * SSOT: Uses config::DomainConfig (modern format), not legacy InstrumentDomain.
     */
    explicit DomainManager(const std::vector<config::DomainConfig>& domains);
    
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
     * @brief Check if ion crossed aperture when exiting domain
     * @param ion Ion state (ion.active set to false if blocked)
     * @param domain_idx Current domain index
     * @param pos_before Local position before integration step [m]
     * @param pos_after Local position after integration step [m]
     * 
     * Detects if ion crossed the end plane (z = domain.length_m).
     * If crossing occurred, interpolates crossing point and checks radial distance.
     * If r_cross > end_aperture_m, ion is blocked (ion.active = false).
     * 
     * Called when ion exits a domain to enforce aperture constraints.
     */
    void check_aperture_crossing(IonState& ion, int domain_idx,
                                  const Vec3& pos_before, const Vec3& pos_after) const;
    
    /**
     * @brief Update ion domain-specific properties
     * @param ion Ion state (modified)
     * @param domain_idx New domain index
     * 
     * Updates:
     * - domain_neutral_mass_kg
     * - domain_temperature_K
     * - domain_particle_density_m3
     * - domain_gas_velocity_m_s
     * - current_domain_index
     * 
     * Called when ion transitions from one domain to another.
     */
    void update_domain_properties(IonState& ion, int domain_idx) const;
    
    /**
     * @brief Terminate ion at boundary (set position to intersection, velocity to 0)
     * @param ion Ion state (modified: pos corrected, vel=0, active=false)
     * @param domain_idx Current domain index
     * @param pos_before_local Local position before integration step [m]
     * @param pos_after_local Local position after integration step [m]
     * 
     * When an ion leaves the domain (detected via is_inside_domain check),
     * this method:
     * 1. Computes intersection point between trajectory line and domain boundary
     * 2. Sets ion position to intersection (not extrapolated pos_after)
     * 3. Sets velocity to zero (ion absorbed/stopped at boundary)
     * 4. Deactivates ion (ion.active = false)
     * 
     * **Why this matters:**
     * - Prevents unphysical positions in trajectory output (pos outside domain)
     * - Makes final position physically meaningful (where ion actually hit)
     * - Useful for analyzing loss patterns (e.g., radial losses in drift tubes)
     * 
     * **Boundary types handled:**
     * - Cylindrical instruments (IMS, TOF, LQIT, SLIM, etc.):
     *   * Radial boundary (r > radius_m): Ion hit cylindrical wall
     *   * Axial boundaries (z < 0 or z > length_m): Ion exited entrance/exit
     *   * Aperture (r > end_aperture_m at z=length_m): Ion blocked by aperture
     * - Orbitrap (hyperbolic electrodes):
     *   * Uses midpoint approximation (analytical ray-hyperbola intersection TODO)
     * 
     * Coordinates are in LOCAL frame (must transform to global after).
     */
    void terminate_ion_at_boundary(IonState& ion, int domain_idx,
                                    const Vec3& pos_before_local,
                                    const Vec3& pos_after_local) const;
    
    /**
     * @brief Get number of domains
     */
    size_t num_domains() const { return domains_.size(); }
    
private:
    const std::vector<config::DomainConfig>& domains_;  ///< Reference to domain list (not owned, SSOT)
    
    /**
     * @brief Check if position is inside domain (internal helper)
     * @param dom Domain to check (SSOT: config::DomainConfig)
     * @param pos Global position [m]
     * @return true if inside domain boundaries
     * 
     * Handles cylindrical geometry (most instruments) and hyperbolic (Orbitrap).
     * Replaces legacy isInsideDomain() from paramUtils.cpp.
     */
    bool is_inside_domain(const config::DomainConfig& dom, const Vec3& pos) const;
};

}  // namespace integrator
}  // namespace ICARION
