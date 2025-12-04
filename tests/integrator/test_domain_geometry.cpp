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

// Helper replicated from OrbitrapGeometry for boundary expectation
double orbitrap_surface_residual(double r, double z, double R, double R_m) {
    const double term1 = z * z;
    const double term2 = 0.5 * (r * r - R * R);
    const double term3 = R_m * R_m * std::log(R / r);
    return term1 - term2 - term3;
}

double orbitrap_r_for_z(double z, double R, double R_m) {
    const double z_abs = std::fabs(z);
    const double eps = 1e-10;
    const int max_iter = 80;

    if (z_abs < eps) {
        return R;
    }

    double r_lo = 0.1 * R;
    double r_hi = R;

    double f_lo = orbitrap_surface_residual(r_lo, z_abs, R, R_m);
    double f_hi = orbitrap_surface_residual(r_hi, z_abs, R, R_m);

    int expand_iter = 0;
    while (f_lo * f_hi > 0.0 && expand_iter < 10) {
        r_lo *= 0.5;
        r_hi *= 1.5;
        f_lo = orbitrap_surface_residual(r_lo, z_abs, R, R_m);
        f_hi = orbitrap_surface_residual(r_hi, z_abs, R, R_m);
        expand_iter++;
    }

    double r_mid = R;
    for (int i = 0; i < max_iter; ++i) {
        r_mid = 0.5 * (r_lo + r_hi);
        double f_mid = orbitrap_surface_residual(r_mid, z_abs, R, R_m);
        if (f_lo * f_mid <= 0.0) {
            r_hi = r_mid;
            f_hi = f_mid;
        } else {
            r_lo = r_mid;
            f_lo = f_mid;
        }
        if (std::fabs(f_mid) < eps) break;
    }
    return r_mid;
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

TEST_CASE("OrbitrapGeometry near-boundary sampling", "[geometry][orbitrap][boundary]") {
    auto dom = make_orbitrap_domain();
    OrbitrapGeometry geom(dom);

    const double Rin = dom.geometry.radius_in_m;
    const double Rout = dom.geometry.radius_out_m;
    const double Rm = dom.geometry.radius_char_m;
    const double z_samples[] = {0.0, 0.02, 0.05, 0.09};

    for (double z : z_samples) {
        double r_in = orbitrap_r_for_z(z, Rin, Rm);
        double r_out = orbitrap_r_for_z(z, Rout, Rm);
        // If the corridor degenerates numerically, skip this sample
        if (r_out <= r_in) {
            continue;
        }
        const double delta_in = std::max(1e-4, 5e-3 * r_in);   // push clearly outside inner surface
        const double delta_out = std::max(1e-4, 5e-3 * r_out);  // push clearly outside outer surface
        const double delta_inside = std::max(1e-5, 1e-3 * r_out);
        // Mid-corridor should be inside
        double r_mid = 0.5 * (r_in + r_out);
        REQUIRE(geom.contains(Vec3{r_mid, 0.0, z}));
        // Just inside corridor
        REQUIRE(geom.contains(Vec3{r_in + delta_inside, 0.0, z}));
        REQUIRE(geom.contains(Vec3{r_out - delta_inside, 0.0, z}));
        // Just outside corridor
        REQUIRE_FALSE(geom.contains(Vec3{r_in - delta_in, 0.0, z}));
        REQUIRE_FALSE(geom.contains(Vec3{r_out + delta_out, 0.0, z}));
    }
}
