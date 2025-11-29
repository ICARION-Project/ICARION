// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#pragma once

#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "DomainManager.h"
#include <memory>

namespace ICARION {
namespace integrator {

/**
 * @brief RAII wrapper for ion state in local domain coordinates
 * 
 * Problem this solves:
 * - Eliminates manual copy-back of ion properties (SSOT violation)
 * - Automatic coordinate transformation on destruction
 * - No risk of forgetting to sync properties after collision/reaction
 * 
 * Usage:
 * ```cpp
 * {
 *     DomainContext ctx(ion, domain_idx, *domain_manager_);
 *     collision_handler->handle_collision(ctx.ion(), ...);
 *     reaction_handler->handle_reaction(ctx.ion(), ...);
 *     // Automatic sync on scope exit
 * }
 * ```
 * 
 * Design:
 * - ion is a reference (no copy, SSOT maintained)
 * - pos_local/vel_local are computed views
 * - Destructor transforms back to global coordinates
 */
class DomainContext {
public:
    /**
     * @brief Construct domain context and transform to local coordinates
     * 
     * @param ion Reference to ion state (will be modified in-place)
     * @param domain_idx Current domain index
     * @param domain_manager Manager for coordinate transformations
     */
    DomainContext(
        IonState& ion, 
        int domain_idx,
        const DomainManager& domain_manager
    ) : ion_(ion),
        domain_idx_(domain_idx),
        domain_manager_(domain_manager),
        pos_local_(domain_manager.global_to_local_pos(ion.pos, domain_idx)),
        vel_local_(domain_manager.global_to_local_vel(ion.vel, domain_idx))
    {}
    
    /**
     * @brief Destructor: Transform local coordinates back to global
     * 
     * This ensures ion.pos and ion.vel are always in sync with local state.
     * All other properties (species_id, mass, CCS, etc.) are already updated
     * in-place by handlers, so no manual copy-back needed!
     */
    ~DomainContext() {
        ion_.pos = domain_manager_.local_to_global_pos(pos_local_, domain_idx_);
        ion_.vel = domain_manager_.local_to_global_vel(vel_local_, domain_idx_);
    }
    
    // Prevent copying (would break RAII semantics)
    DomainContext(const DomainContext&) = delete;
    DomainContext& operator=(const DomainContext&) = delete;
    
    // Allow moving (for return optimization)
    DomainContext(DomainContext&&) = default;
    DomainContext& operator=(DomainContext&&) = default;
    
    // Accessors
    IonState& ion() { return ion_; }
    const IonState& ion() const { return ion_; }
    
    Vec3& pos_local() { return pos_local_; }
    const Vec3& pos_local() const { return pos_local_; }
    
    Vec3& vel_local() { return vel_local_; }
    const Vec3& vel_local() const { return vel_local_; }
    
    int domain_idx() const { return domain_idx_; }
    
    /**
     * @brief Update ion position/velocity from local coordinates
     * 
     * Use this if handlers modify pos_local/vel_local and you need
     * to sync ion state before destruction (e.g., for intermediate checks).
     */
    void sync_to_ion() {
        ion_.pos = pos_local_;
        ion_.vel = vel_local_;
    }
    
    /**
     * @brief Update local coordinates from ion state
     * 
     * Use this if ion properties were modified externally and you need
     * to recompute local coordinates.
     */
    void sync_from_ion() {
        pos_local_ = ion_.pos;
        vel_local_ = ion_.vel;
    }

private:
    IonState& ion_;                      ///< Reference to actual ion state (SSOT)
    int domain_idx_;                     ///< Current domain index
    const DomainManager& domain_manager_; ///< For coordinate transforms
    Vec3 pos_local_;                     ///< Position in local domain coordinates
    Vec3 vel_local_;                     ///< Velocity in local domain coordinates
};

}  // namespace integrator
}  // namespace ICARION
