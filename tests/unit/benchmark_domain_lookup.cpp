// SPDX-License-Identifier: GPL-3.0-only
// Domain lookup micro-benchmark: prefilter vs linear scan for cylindrical domains.

#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <random>
#include "core/integrator/DomainManager.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/CylindricalGeometry.h"

using ICARION::integrator::DomainManager;
using ICARION::config::DomainConfig;
using ICARION::config::Instrument;
using ICARION::core::Vec3;

namespace {

DomainConfig make_cylinder(double z0, double length, double radius) {
    DomainConfig d;
    d.name = "d";
    d.instrument = Instrument::IMS;
    d.geometry.origin_m = {0.0, 0.0, z0};
    d.geometry.length_m = length;
    d.geometry.radius_m = radius;
    d.finalize();
    return d;
}

}  // namespace

TEST_CASE("Domain lookup benchmark (prefilter vs linear)", "[benchmark][domain][lookup]") {
    const size_t n_domains = 10;
    const size_t n_points = 100000;
    std::vector<DomainConfig> domains;
    domains.reserve(n_domains);
    for (size_t i = 0; i < n_domains; ++i) {
        double z0 = 0.1 * static_cast<double>(i);
        domains.push_back(make_cylinder(z0, 0.1, 0.02));
    }

    DomainManager mgr(domains, 123);

    // Build linear geometry list for reference timing
    std::vector<ICARION::config::CylindricalGeometry> geoms;
    for (const auto& d : domains) {
        geoms.emplace_back(d);
    }

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist_z(0.0, n_domains * 0.1);
    std::uniform_real_distribution<double> dist_r(-0.03, 0.03);

    std::vector<Vec3> points;
    points.reserve(n_points);
    for (size_t i = 0; i < n_points; ++i) {
        points.push_back({dist_r(rng), dist_r(rng), dist_z(rng)});
    }

    auto time_prefilter = [&]() {
        auto start = std::chrono::high_resolution_clock::now();
        size_t hits = 0;
        for (const auto& p : points) {
            int idx = mgr.find_domain_index(p);
            if (idx >= 0) ++hits;
        }
        auto end = std::chrono::high_resolution_clock::now();
        return std::make_pair(
            std::chrono::duration<double, std::milli>(end - start).count(),
            hits);
    };

    auto time_linear = [&]() {
        auto start = std::chrono::high_resolution_clock::now();
        size_t hits = 0;
        for (const auto& p : points) {
            int idx = -1;
            for (size_t i = 0; i < geoms.size(); ++i) {
                if (geoms[i].contains(p)) {
                    idx = static_cast<int>(i);
                    break;
                }
            }
            if (idx >= 0) ++hits;
        }
        auto end = std::chrono::high_resolution_clock::now();
        return std::make_pair(
            std::chrono::duration<double, std::milli>(end - start).count(),
            hits);
    };

    auto [ms_prefilter, hits_prefilter] = time_prefilter();
    auto [ms_linear, hits_linear] = time_linear();

    REQUIRE(hits_prefilter == hits_linear);
    INFO("prefilter_ms=" << ms_prefilter << " linear_ms=" << ms_linear);
}
