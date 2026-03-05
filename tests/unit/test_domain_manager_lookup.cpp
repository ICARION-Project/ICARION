// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer

#include <catch2/catch_test_macros.hpp>
#include "core/integrator/DomainManager.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/Vec3.h"

using ICARION::integrator::DomainManager;
using ICARION::config::DomainConfig;
using ICARION::config::Instrument;
using ICARION::config::BoundaryActionType;
using ICARION::core::Vec3;

namespace {

DomainConfig make_cylinder(const std::string& name, double z0, double length, double radius) {
    DomainConfig d;
    d.name = name;
    d.instrument = Instrument::IMS;
    d.geometry.origin_m = {0.0, 0.0, z0};
    d.geometry.length_m = length;
    d.geometry.radius_m = radius;
    d.boundary.type = BoundaryActionType::Absorption;
    d.finalize();
    return d;
}

int expected_idx(const std::vector<DomainConfig>& domains, const Vec3& p) {
    for (size_t i = 0; i < domains.size(); ++i) {
        const auto& g = domains[i].geometry;
        double z0 = g.origin_m.z;
        double z1 = z0 + g.length_m;
        double z_min = std::min(z0, z1);
        double z_max = std::max(z0, z1);
        double r2 = p.x * p.x + p.y * p.y;
        if (p.z >= z_min && p.z <= z_max && r2 <= g.radius_m * g.radius_m) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

}  // namespace

TEST_CASE("DomainManager axial prefilter matches brute-force for cylinders", "[domain][lookup]") {
    std::vector<DomainConfig> domains;
    domains.push_back(make_cylinder("d0", 0.0, 0.1, 0.01));
    domains.push_back(make_cylinder("d1", 0.1, 0.1, 0.02));
    domains.push_back(make_cylinder("d2", 0.2, 0.05, 0.005));

    DomainManager mgr(domains, /*rng_seed=*/42);

    std::vector<Vec3> points = {
        {0.0, 0.0, 0.05},   // d0
        {0.0, 0.0, 0.15},   // d1
        {0.0, 0.0, 0.22},   // d2
        {0.03, 0.0, 0.15},  // outside d1 radius
        {0.0, 0.0, 0.5}     // outside all
    };

    for (const auto& p : points) {
        int idx = mgr.find_domain_index(p);
        int ref = expected_idx(domains, p);
        REQUIRE(idx == ref);
    }
}
