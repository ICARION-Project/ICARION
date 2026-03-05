// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "DomainManager.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/config/types/AnalyticalFieldModel.h"
#include "core/config/types/CylindricalGeometry.h"
#include "core/config/types/OrbitrapGeometry.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

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
            has_orbitrap_ = true;
        } else {
            geometries_.push_back(std::make_unique<config::CylindricalGeometry>(domain));
            // Precompute axial span for cylindrical domains
            double z0 = domain.geometry.origin_m.z;
            double z1 = z0 + domain.geometry.length_m;
            spans_.push_back(DomainSpan{
                std::min(z0, z1) - DOMAIN_BOUNDARY_EPSILON,
                std::max(z0, z1) + DOMAIN_BOUNDARY_EPSILON,
                domain.geometry.radius_m + DOMAIN_BOUNDARY_EPSILON,
                static_cast<int>(spans_.size())  // will fix idx below
            });
        }
    }

    // Fix span indices to match domain ordering
    for (size_t i = 0, j = 0; i < domains_.size(); ++i) {
        if (domains_[i].instrument == config::Instrument::Orbitrap) continue;
        if (j < spans_.size()) {
            spans_[j].idx = static_cast<int>(i);
            ++j;
        }
    }
    std::sort(spans_.begin(), spans_.end(),
              [](const DomainSpan& a, const DomainSpan& b) { return a.z_min < b.z_min; });
}

int DomainManager::find_domain_index(const Vec3& pos) const {
    // Fallback to linear scan if index not applicable
    auto linear_lookup = [&]() -> int {
        for (size_t i = 0; i < domains_.size(); ++i) {
            if (geometries_[i]->contains(pos)) {
                return static_cast<int>(i);
            }
        }
        return -1;
    };

    if (has_orbitrap_ || spans_.empty()) {
        return linear_lookup();
    }

    // Axial prefilter: find spans whose z-range contains pos.z
    // spans_ sorted by z_min
    auto it = std::lower_bound(spans_.begin(), spans_.end(), pos.z,
                               [](const DomainSpan& span, double z) { return span.z_max < z; });
    for (; it != spans_.end(); ++it) {
        if (pos.z < it->z_min) {
            break;  // no further matches
        }
        // Optional radial check before full geometry
        if ((pos.x * pos.x + pos.y * pos.y) > it->radius * it->radius) {
            continue;
        }
        int idx = it->idx;
        if (idx >= 0 && idx < static_cast<int>(geometries_.size()) && geometries_[idx]->contains(pos)) {
            return idx;
        }
    }

    // Fallback if not found via spans
    return linear_lookup();

    // Original linear scan (kept for reference)
    // for (size_t i = 0; i < domains_.size(); ++i) {
    //     if (geometries_[i]->contains(pos)) {
    //         return static_cast<int>(i);
    //     }
    // }
    // return -1;
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
