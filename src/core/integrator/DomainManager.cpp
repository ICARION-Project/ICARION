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
        // Use geometry strategy directly; domain_index may be unset in tests
        if (geometries_[i]->contains(pos)) {
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

BoundaryAction* DomainManager::boundary_action(int idx) const {
    if (idx < 0 || idx >= static_cast<int>(boundary_actions_.size())) {
        return nullptr;
    }
    return boundary_actions_[idx].get();
}

Vec3 DomainManager::surface_normal_global(const Vec3& global_pos, int domain_idx) const {
    if (domain_idx < 0 || domain_idx >= static_cast<int>(geometries_.size())) {
        return Vec3{0.0, 0.0, 0.0};
    }
    Vec3 local = geometries_[static_cast<size_t>(domain_idx)]->global_to_local_pos(global_pos);
    Vec3 n_local = geometries_[static_cast<size_t>(domain_idx)]->surface_normal(local);
    return geometries_[static_cast<size_t>(domain_idx)]->local_to_global_vel(n_local); // rotation only
}

bool DomainManager::boundary_intersection_global(const Vec3& start_global,
                                                 const Vec3& end_global,
                                                 int domain_idx,
                                                 Vec3& intersection_global,
                                                 Vec3& normal_global) const {
    if (domain_idx < 0 || domain_idx >= static_cast<int>(geometries_.size())) {
        return false;
    }
    const auto& geom = geometries_[static_cast<size_t>(domain_idx)];
    Vec3 start_local = geom->global_to_local_pos(start_global);
    Vec3 end_local = geom->global_to_local_pos(end_global);
    Vec3 hit_local;
    if (!geom->first_boundary_intersection(start_local, end_local, hit_local)) {
        return false;
    }
    intersection_global = geom->local_to_global_pos(hit_local);
    Vec3 n_local = geom->surface_normal(hit_local);
    normal_global = geom->local_to_global_vel(n_local);
    return true;
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

}  // namespace integrator
}  // namespace ICARION
