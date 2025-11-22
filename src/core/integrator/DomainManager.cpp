// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "DomainManager.h"
#include "core/param/paramUtils.h"  // InstrumentDomain
#include "instrument/InstrumentTypes.h"  // Instrument enum
#include <stdexcept>
#include <cmath>

namespace ICARION {
namespace integrator {

DomainManager::DomainManager(const std::vector<InstrumentDomain>& domains)
    : domains_(domains)
{
    if (domains.empty()) {
        throw std::invalid_argument("DomainManager: domains vector is empty");
    }
}

int DomainManager::find_domain_index(const Vec3& pos) const {
    for (size_t i = 0; i < domains_.size(); ++i) {
        if (is_inside_domain(domains_[i], pos)) {
            return static_cast<int>(i);
        }
    }
    return -1;  // Outside all domains
}

const InstrumentDomain& DomainManager::get_domain(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(domains_.size())) {
        throw std::out_of_range("DomainManager::get_domain: invalid domain index " + 
                                std::to_string(idx));
    }
    return domains_[idx];
}

Vec3 DomainManager::global_to_local_pos(const Vec3& pos, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_global_to_local * (pos - dom.geom.origin_m);
}

Vec3 DomainManager::global_to_local_vel(const Vec3& vel, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_global_to_local * vel;
}

Vec3 DomainManager::local_to_global_pos(const Vec3& pos_local, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_local_to_global * pos_local + dom.geom.origin_m;
}

Vec3 DomainManager::local_to_global_vel(const Vec3& vel_local, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    return dom.rotation_local_to_global * vel_local;
}

void DomainManager::check_aperture_crossing(IonState& ion, int domain_idx,
                                             const Vec3& pos_before, const Vec3& pos_after) const {
    const auto& dom = get_domain(domain_idx);
    
    // No aperture constraint
    if (dom.geom.end_aperture_m <= 0.0) {
        return;  // Ion can pass freely
    }
    
    const double z_ap = dom.geom.length_m;
    
    // Check if ion crossed the aperture plane (forward or backward)
    const bool crossed_forward = (pos_before.z < z_ap && pos_after.z >= z_ap);
    const bool crossed_backward = (pos_before.z > z_ap && pos_after.z <= z_ap);
    
    if (crossed_forward || crossed_backward) {
        // Interpolate crossing point
        const double alpha = (z_ap - pos_before.z) / (pos_after.z - pos_before.z);
        const Vec3 cross_point = pos_before + (pos_after - pos_before) * alpha;
        
        // Check radial distance at crossing
        const double r_cross = std::sqrt(cross_point.x * cross_point.x + 
                                         cross_point.y * cross_point.y);
        
        if (r_cross > dom.geom.end_aperture_m) {
            // Ion blocked by aperture
            ion.active = false;
        }
    }
}

void DomainManager::update_domain_properties(IonState& ion, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    
    ion.domain_neutral_mass_kg = dom.env.neutral_mass_kg;
    ion.domain_temperature_K = dom.env.temperature_K;
    ion.domain_particle_density_m3 = dom.env.particle_density_m_3;
    ion.domain_gas_velocity_m_s = Vec3(dom.env.gas_velocity_m_s.x,
                                       dom.env.gas_velocity_m_s.y,
                                       dom.env.gas_velocity_m_s.z);
    ion.current_domain_index = dom.index;
}

bool DomainManager::is_inside_domain(const InstrumentDomain& dom, const Vec3& globalPos) const {
    // Transform to local coordinates
    Vec3 local = dom.rotation_global_to_local * (globalPos - dom.geom.origin_m);
    double r = std::sqrt(local.x * local.x + local.y * local.y);
    
    // Cylindrical geometry (most instruments)
    if (dom.instrument != Instrument::Orbitrap) {
        // Ion must be inside domain:
        // - z >= 0 (at or past entrance)
        // - z < length_m (before exit)
        // - r < radius_m (within radial boundary)
        return (local.z >= 0.0 && local.z < dom.geom.length_m) && (r < dom.geom.radius_m);
    }
    
    // Orbitrap: hyperbolic electrode geometry
    const double Rin = dom.geom.radius_in_m;
    const double Rout = dom.geom.radius_out_m;
    const double Rm = dom.geom.radius_char_m;
    const double z = std::fabs(local.z);
    
    // Compute allowed radial range for given z (hyperbolic surfaces)
    // r_inner(z) = sqrt(Rin^2 + (z/Rm)^2)
    // r_outer(z) = sqrt(Rout^2 + (z/Rm)^2)
    const double r_in_allowed  = std::sqrt(Rin * Rin + (z / Rm) * (z / Rm));
    const double r_out_allowed = std::sqrt(Rout * Rout + (z / Rm) * (z / Rm));
    
    // Inside domain if between hyperbolic surfaces
    return (r >= r_in_allowed) && (r <= r_out_allowed);
}

}  // namespace integrator
}  // namespace ICARION
