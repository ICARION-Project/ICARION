// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "DomainManager.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/config/types/AnalyticalFieldModel.h"
#include "core/config/types/CylindricalGeometry.h"
#include "core/config/types/OrbitrapGeometry.h"
#include <stdexcept>
#include <cmath>

namespace ICARION {
namespace integrator {

DomainManager::DomainManager(
    const std::vector<config::DomainConfig>& domains,
    unsigned int rng_seed
) : domains_(domains), rng_(rng_seed)
{
    if (domains.empty()) {
        throw std::invalid_argument("DomainManager: domains vector is empty");
    }
    
    // Create boundary actions, analytical field models, and geometry strategies for each domain.
    // PhysicsSetup will override field models with grid-backed versions when present.
    boundary_actions_.reserve(domains.size());
    field_models_.reserve(domains.size());
    geometries_.reserve(domains.size());
    for (const auto& domain : domains) {
        boundary_actions_.push_back(
            BoundaryActionFactory::create(domain.boundary, &rng_)
        );
        field_models_.push_back(
            std::make_unique<config::AnalyticalFieldModel>(domain)
        );
        if (domain.instrument == config::Instrument::Orbitrap) {
            geometries_.push_back(std::make_unique<config::OrbitrapGeometry>(domain));
        } else {
            geometries_.push_back(std::make_unique<config::CylindricalGeometry>(domain));
        }
    }
}

int DomainManager::find_domain_index(const Vec3& pos) const {
    for (size_t i = 0; i < domains_.size(); ++i) {
        if (is_inside_domain(domains_[i], pos)) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

const config::DomainConfig& DomainManager::get_domain(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(domains_.size())) {
        throw std::out_of_range("DomainManager::get_domain: invalid index " + 
                                std::to_string(idx));
    }
    return domains_[idx];
}

const config::IFieldModel* DomainManager::field_model(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(field_models_.size())) {
        return nullptr;
    }
    return field_models_[idx].get();
}

Vec3 DomainManager::global_to_local_pos(const Vec3& pos, int domain_idx) const {
    if (domain_idx < 0 || domain_idx >= static_cast<int>(geometries_.size())) {
        throw std::out_of_range("DomainManager::global_to_local_pos: invalid index");
    }
    return geometries_[domain_idx]->global_to_local_pos(pos);
}

Vec3 DomainManager::global_to_local_vel(const Vec3& vel, int domain_idx) const {
    if (domain_idx < 0 || domain_idx >= static_cast<int>(geometries_.size())) {
        throw std::out_of_range("DomainManager::global_to_local_vel: invalid index");
    }
    return geometries_[domain_idx]->global_to_local_vel(vel);
}

Vec3 DomainManager::local_to_global_pos(const Vec3& pos_local, int domain_idx) const {
    if (domain_idx < 0 || domain_idx >= static_cast<int>(geometries_.size())) {
        throw std::out_of_range("DomainManager::local_to_global_pos: invalid index");
    }
    return geometries_[domain_idx]->local_to_global_pos(pos_local);
}

Vec3 DomainManager::local_to_global_vel(const Vec3& vel_local, int domain_idx) const {
    if (domain_idx < 0 || domain_idx >= static_cast<int>(geometries_.size())) {
        throw std::out_of_range("DomainManager::local_to_global_vel: invalid index");
    }
    return geometries_[domain_idx]->local_to_global_vel(vel_local);
}

void DomainManager::check_aperture_crossing(IonState& ion, int domain_idx,
                                             const Vec3& pos_before, const Vec3& pos_after) const {
    const auto& dom = get_domain(domain_idx);
    
    double aperture = geometries_[domain_idx]->end_aperture();
    if (aperture <= 0.0) {
        return;  // No aperture constraint
    }
    
    Vec3 pb_local = global_to_local_pos(pos_before, domain_idx);
    Vec3 pa_local = global_to_local_pos(pos_after, domain_idx);
    
    double length = geometries_[domain_idx]->length();
    bool before_exit = (pb_local.z < length);
    bool after_exit = (pa_local.z >= length);
    
    if (before_exit && after_exit) {
        double t = (length - pb_local.z) / (pa_local.z - pb_local.z);
        Vec3 crossing = pb_local + (pa_local - pb_local) * t;
        double r_cross = std::sqrt(crossing.x*crossing.x + crossing.y*crossing.y);
        
        if (r_cross > aperture) {
            ion.active = false;
        }
    }
}

void DomainManager::update_domain_properties(IonState& ion, int domain_idx) const {
    const auto& dom = get_domain(domain_idx);
    
    // Domain cache fields removed from IonState - only update domain index
    ion.current_domain_index = dom.domain_index;
}

void DomainManager::terminate_ion_at_boundary(IonState& ion, int domain_idx,
                                                const Vec3& pos_before_local,
                                                const Vec3& pos_after_local) const {
    const auto& dom = get_domain(domain_idx);
    
    // ========== ORBITRAP ==========
    if (dom.instrument == config::Instrument::Orbitrap) {
        // Geometry strategy determined "outside"; terminate at midpoint for now.
        Vec3 intersection = (pos_before_local + pos_after_local) * 0.5;
        ion.pos = local_to_global_pos(intersection, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.death_time_s = ion.t;
        ion.active = false;
        return;
    }
    
    // ========== CYLINDRICAL GEOMETRY ==========
    Vec3 dir = pos_after_local - pos_before_local;
    double dir_len = std::sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);
    
    if (dir_len < NUMERICAL_ZERO) {
        ion.pos = local_to_global_pos(pos_before_local, domain_idx);
        ion.vel = {0.0, 0.0, 0.0};
        ion.death_time_s = ion.t;
        ion.active = false;
        return;
    }
    
    dir = dir * (1.0 / dir_len);
    
    double t_min = dir_len;
    Vec3 intersection = pos_after_local;
    
    // 1. Radial boundary
    const double R = geometries_[domain_idx]->radius();
    const double a = dir.x*dir.x + dir.y*dir.y;
    const double b = 2.0*(pos_before_local.x*dir.x + pos_before_local.y*dir.y);
    const double c = pos_before_local.x*pos_before_local.x + 
                     pos_before_local.y*pos_before_local.y - R*R;
    
    if (a > NUMERICAL_ZERO) {
        const double discriminant = b*b - 4.0*a*c;
        if (discriminant >= 0.0) {
            const double sqrt_disc = std::sqrt(discriminant);
            const double t1 = (-b + sqrt_disc) / (2.0*a);
            const double t2 = (-b - sqrt_disc) / (2.0*a);
            
            for (double t : {t1, t2}) {
                if (t > 0.0 && t < t_min) {
                    Vec3 candidate = pos_before_local + dir * t;
                    if (candidate.z >= -DOMAIN_BOUNDARY_EPSILON && 
                        candidate.z <= geometries_[domain_idx]->length()) {
                        t_min = t;
                        intersection = candidate;
                    }
                }
            }
        }
    }
    
    // 2. Entrance plane (z=0)
    if (std::abs(dir.z) > NUMERICAL_ZERO) {
        double t = (0.0 - pos_before_local.z) / dir.z;
        if (t > 0.0 && t < t_min) {
            Vec3 candidate = pos_before_local + dir * t;
            double r = std::sqrt(candidate.x*candidate.x + candidate.y*candidate.y);
            if (r <= R) {
                t_min = t;
                intersection = candidate;
            }
        }
    }
    
    // 3. Exit plane (z=length_m)
    if (std::abs(dir.z) > NUMERICAL_ZERO) {
        double t = (geometries_[domain_idx]->length() - pos_before_local.z) / dir.z;
        if (t > 0.0 && t < t_min) {
            Vec3 candidate = pos_before_local + dir * t;
            double r = std::sqrt(candidate.x*candidate.x + candidate.y*candidate.y);
            double aperture_limit = (geometries_[domain_idx]->end_aperture() > 0.0) ? 
                                    geometries_[domain_idx]->end_aperture() : R;
            if (r <= aperture_limit) {
                t_min = t;
                intersection = candidate;
            }
        }
    }
    
    // Compute surface normal at intersection (pointing inward)
    Vec3 normal = geometries_[domain_idx]->surface_normal(intersection);
    
    // Apply boundary action
    const auto& boundary_action = boundary_actions_[domain_idx];
    boundary_action->apply(
        ion,
        normal,
        local_to_global_pos(intersection, domain_idx),
        dom.environment.temperature_K,
        ion.t  // Current simulation time
    );
}

Vec3 DomainManager::compute_surface_normal(const Vec3& pos_local, int domain_idx) const {
    return geometries_[domain_idx]->surface_normal(pos_local);
}

bool DomainManager::is_inside_domain(const config::DomainConfig& dom, const Vec3& globalPos) const {
    if (dom.domain_index >= 0 && dom.domain_index < static_cast<int>(geometries_.size())) {
        return geometries_[dom.domain_index]->contains(globalPos);
    }
    return false;
}

}  // namespace integrator
}  // namespace ICARION
