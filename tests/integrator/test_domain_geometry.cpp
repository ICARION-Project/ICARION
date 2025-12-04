// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/config/types/CylindricalGeometry.h"
#include "core/config/types/OrbitrapGeometry.h"
#include "core/config/types/DomainConfig.h"
#include "core/utils/mathUtils.h"

using ICARION::config::CylindricalGeometry;
using ICARION::config::OrbitrapGeometry;
using ICARION::config::DomainConfig;
using ICARION::core::Vec3;
using Catch::Approx;

namespace {

DomainConfig make_cyl_domain() {
    DomainConfig dom;
    dom.geometry.length_m = 0.10;
    dom.geometry.radius_m = 0.01;
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    return dom;
}

DomainConfig make_orbitrap_domain() {
    DomainConfig dom;
    dom.geometry.length_m = 0.20;
    dom.geometry.radius_in_m = 0.02;
    dom.geometry.radius_out_m = 0.04;
    dom.geometry.radius_char_m = 0.03;
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    return dom;
}

} // namespace

TEST_CASE("CylindricalGeometry contains/transform", "[geometry][cylindrical]") {
    auto dom = make_cyl_domain();
    CylindricalGeometry geom(dom);

    // Inside
    REQUIRE(geom.contains(Vec3{0.0, 0.0, 0.05}));
    REQUIRE(geom.contains(Vec3{0.005, -0.004, 0.01}));

    // Outside radial
    REQUIRE_FALSE(geom.contains(Vec3{0.02, 0.0, 0.05}));
    // Outside axial
    REQUIRE_FALSE(geom.contains(Vec3{0.0, 0.0, 0.11}));

    // Round-trip transform
    Vec3 local{0.003, -0.002, 0.04};
    Vec3 global = geom.local_to_global_pos(local);
    Vec3 back = geom.global_to_local_pos(global);
    REQUIRE(back.x == Approx(local.x));
    REQUIRE(back.y == Approx(local.y));
    REQUIRE(back.z == Approx(local.z));
}

TEST_CASE("OrbitrapGeometry corridor checks", "[geometry][orbitrap]") {
    auto dom = make_orbitrap_domain();
    OrbitrapGeometry geom(dom);

    // Center region inside corridor
    REQUIRE(geom.contains(Vec3{0.025, 0.0, 0.0}));
    // Near inner radius should be outside
    REQUIRE_FALSE(geom.contains(Vec3{0.0, 0.0, 0.0}));
    REQUIRE_FALSE(geom.contains(Vec3{0.015, 0.0, 0.0}));
    // Outside outer radius rejected
    REQUIRE_FALSE(geom.contains(Vec3{0.05, 0.0, 0.0}));
    // Axial beyond length rejected
    REQUIRE_FALSE(geom.contains(Vec3{0.025, 0.0, 0.15}));

    // Symmetry in z: point at mid-length still inside if radial fits corridor
    REQUIRE(geom.contains(Vec3{0.025, 0.0, 0.05}));
}
