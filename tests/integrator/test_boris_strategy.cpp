// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/integrator/strategies/BorisStrategy.h"
#include "core/integrator/strategies/IntegrationStrategyFactory.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/IForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

#include <cmath>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;
using ICARION::core::IonEnsemble;
using Catch::Matchers::WithinAbs;

// =============================================================================
// TEST FORCE: Constant Electric Field
// =============================================================================
class ConstantElectricForce : public IForce {
public:
    explicit ConstantElectricForce(const Vec3& E_field) : E_(E_field) {}
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        double q = ensemble.charge_data()[ion_idx];
        return E_ * q;
    }

    Vec3 compute_soa(const ForceState& state, double t,
                     const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        return E_ * state.ion_charge_C;
    }
    
    std::string name() const override { return "ConstantElectric"; }
    
private:
    Vec3 E_;
};

// =============================================================================
// HELPER FUNCTIONS
// =============================================================================

DomainConfig create_test_domain() {
    DomainConfig domain;
    domain.geometry.length_m = 1.0;
    domain.geometry.radius_m = 0.1;
    domain.geometry.min_bound = Vec3{-0.1, -0.1, 0.0};
    domain.geometry.max_bound = Vec3{0.1, 0.1, 1.0};
    domain.fields.magnetic.enabled = false;
    return domain;
}

DomainConfig create_test_domain_with_B(double Bx, double By, double Bz) {
    DomainConfig domain = create_test_domain();
    domain.fields.magnetic.enabled = true;
    domain.fields.magnetic.field_strength_T = Vec3{Bx, By, Bz};
    return domain;
}

IonState create_test_ion(double mass_kg = 1.0, double charge_C = 1.0) {
    IonState ion;
    ion.mass_kg = mass_kg;
    ion.ion_charge_C = charge_C;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{0, 0, 0};
    return ion;
}

// =============================================================================
// TEST: Basic Functionality
// =============================================================================

TEST_CASE("BorisStrategy: Constructor and properties", "[integrator][boris]") {
    BorisStrategy strategy;
    
    SECTION("Name") {
        REQUIRE(strategy.name() == "Boris");
    }
    
    SECTION("Not adaptive") {
        REQUIRE(strategy.is_adaptive() == false);
    }
}

TEST_CASE("BorisStrategy: Factory creation", "[integrator][boris]") {
    auto strategy = IntegrationStrategyFactory::create("Boris");
    
    REQUIRE(strategy != nullptr);
    REQUIRE(strategy->name() == "Boris");
    REQUIRE(strategy->is_adaptive() == false);
}

// =============================================================================
// TEST: Pure Electric Field (No Magnetic Field)
// =============================================================================

TEST_CASE("BorisStrategy: Constant electric field acceleration", "[integrator][boris]") {
    // Setup: Constant E-field in x-direction
    auto domain = create_test_domain();
    auto ion = create_test_ion(1.0, 1.0);
    
    double E_x = 10.0;  // V/m
    Vec3 E_field{E_x, 0, 0};
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<ConstantElectricForce>(E_field));
    
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;
    
    BorisStrategy strategy;
    
    // Integrate for 1 second
    double t = 0.0;
    double dt = 0.001;  // 1 ms timestep
    int n_steps = 1000;
    
    for (int i = 0; i < n_steps; ++i) {
        strategy.step(ensemble, idx, t, dt, forces);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // Analytical solution: v = a*t, x = 0.5*a*t²
    IonState ion_initial = ensemble.to_legacy()[idx];  // mass/charge same as initial
    double a = E_x * ion_initial.ion_charge_C / ion_initial.mass_kg;
    double t_final = n_steps * dt;
    double vx_expected = a * t_final;
    double x_expected = 0.5 * a * t_final * t_final;
    
    REQUIRE_THAT(ion_out.vel.x, WithinAbs(vx_expected, 1e-4));
    REQUIRE_THAT(ion_out.pos.x, WithinAbs(x_expected, 1e-4));
}

// =============================================================================
// TEST: Pure Magnetic Field (Cyclotron Motion)
// =============================================================================

TEST_CASE("BorisStrategy: Cyclotron motion in uniform B-field", "[integrator][boris]") {
    // Setup: Uniform B-field in z-direction, particle moves in x-y plane
    double B_z = 1.0;  // 1 Tesla
    auto domain = create_test_domain_with_B(0, 0, B_z);
    
    auto ion = create_test_ion(1.0, 1.0);
    ion.vel = Vec3{1.0, 0, 0};  // Initial velocity in x-direction
    
    ForceRegistry forces(domain);  // No electric forces
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;
    
    BorisStrategy strategy;
    
    // Cyclotron frequency: ωc = q*B/m
    double omega_c = ion.ion_charge_C * B_z / ion.mass_kg;
    double T_cyclotron = 2.0 * M_PI / omega_c;  // Period
    
    // Integrate for one cyclotron period
    double t = 0.0;
    double dt = T_cyclotron / 100.0;  // 100 steps per period
    int n_steps = 100;
    
    // Save initial energy
    double v0_mag = std::sqrt(ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z);
    double E0 = 0.5 * ion.mass_kg * v0_mag * v0_mag;
    
    for (int i = 0; i < n_steps; ++i) {
        strategy.step(ensemble, idx, t, dt, forces);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // After one period, should return to approximately initial position
    // (Boris has small phase error, not exact)
    REQUIRE_THAT(ion_out.pos.x, WithinAbs(0.0, 0.05));
    REQUIRE_THAT(ion_out.pos.y, WithinAbs(0.0, 0.05));
    REQUIRE_THAT(ion_out.pos.z, WithinAbs(0.0, 1e-10));
    
    // Velocity magnitude should be conserved (energy conservation)
    double v_mag = std::sqrt(ion_out.vel.x * ion_out.vel.x + ion_out.vel.y * ion_out.vel.y + ion_out.vel.z * ion_out.vel.z);
    double E_final = 0.5 * ion_out.mass_kg * v_mag * v_mag;
    
    REQUIRE_THAT(v_mag, WithinAbs(v0_mag, 1e-6));
    REQUIRE_THAT(E_final, WithinAbs(E0, 1e-12));
}

// =============================================================================
// TEST: E×B Drift
// =============================================================================

TEST_CASE("BorisStrategy: E×B drift motion", "[integrator][boris]") {
    // Setup: E-field in x, B-field in z → drift in y-direction
    double E_x = 1.0;  // V/m
    double B_z = 0.01;  // T (weak field for fast drift)
    
    auto domain = create_test_domain_with_B(0, 0, B_z);
    
    auto ion = create_test_ion(1.0, 1.0);
    ion.vel = Vec3{0, 0, 0};  // Start at rest
    
    Vec3 E_field{E_x, 0, 0};
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<ConstantElectricForce>(E_field));
    
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;
    
    BorisStrategy strategy;
    
    // E×B drift velocity: v_drift = E×B / B²
    // E = (E_x, 0, 0), B = (0, 0, B_z)
    // E×B = (0, -E_x*B_z, 0)  [Correct cross product!]
    double v_drift_y = -E_x * B_z / (B_z * B_z);  // Negative!
    
    // Integrate for multiple cyclotron periods
    // ω_c = qB/m = 1.0 * 0.01 / 1.0 = 0.01 rad/s
    // T_c = 2π/ω_c ≈ 628 s → use 10 periods = 6280 s
    double t = 0.0;
    double dt = 0.01;  // Larger timestep for efficiency
    int n_steps = 100000;  // 1000 seconds
    
    for (int i = 0; i < n_steps; ++i) {
        strategy.step(ensemble, idx, t, dt, forces);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // Check that particle has drifted in y-direction
    double y_expected = v_drift_y * t;
    
    // After many cyclotron periods, drift should dominate
    // Allow 10% tolerance for numerical integration and cyclotron modulation
    REQUIRE(std::fabs(ion_out.pos.y - y_expected) / std::fabs(y_expected) < 0.1);
}

// =============================================================================
// TEST: Long-term Energy Stability (Symplectic Property)
// =============================================================================

TEST_CASE("BorisStrategy: Long-term energy conservation", "[integrator][boris][symplectic]") {
    // Boris pusher is symplectic → energy should not drift over long integration
    double B_z = 0.5;  // T
    auto domain = create_test_domain_with_B(0, 0, B_z);
    
    auto ion = create_test_ion(1.0, 1.0);
    ion.vel = Vec3{1.0, 0.5, 0};
    
    ForceRegistry forces(domain);
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;
    
    BorisStrategy strategy;
    
    // Initial energy
    double v0_sq = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
    double E0 = 0.5 * ion.mass_kg * v0_sq;
    
    // Integrate for 100 cyclotron periods
    double omega_c = ion.ion_charge_C * B_z / ion.mass_kg;
    double T_cyclotron = 2.0 * M_PI / omega_c;
    double t = 0.0;
    double dt = T_cyclotron / 50.0;  // 50 steps per period
    int n_periods = 100;
    int n_steps = n_periods * 50;
    
    for (int i = 0; i < n_steps; ++i) {
        strategy.step(ensemble, idx, t, dt, forces);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // Final energy
    double v_sq = ion_out.vel.x * ion_out.vel.x + ion_out.vel.y * ion_out.vel.y + ion_out.vel.z * ion_out.vel.z;
    double E_final = 0.5 * ion_out.mass_kg * v_sq;
    
    // Energy should be conserved to high precision (symplectic property)
    double rel_error = std::fabs(E_final - E0) / E0;
    REQUIRE(rel_error < 1e-10);  // <0.00000001% error after 100 orbits
}

// =============================================================================
// TEST: Comparison with RK4 for E-field only
// =============================================================================

TEST_CASE("BorisStrategy: Agreement with RK4 for E-field only", "[integrator][boris]") {
    // When B=0, Boris should give similar results to RK4
    auto domain = create_test_domain();  // No B-field
    
    auto ion_boris = create_test_ion(1.0, 1.0);
    auto ion_rk4 = ion_boris;
    
    ion_boris.vel = Vec3{0.1, 0.2, 0};
    ion_rk4.vel = ion_boris.vel;
    
    Vec3 E_field{5.0, 0, 0};
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<ConstantElectricForce>(E_field));
    
    IonEnsemble ensemble_boris = IonEnsemble::from_legacy({ion_boris});
    IonEnsemble ensemble_rk4 = IonEnsemble::from_legacy({ion_rk4});
    const size_t idx = 0;
    
    BorisStrategy boris;
    auto rk4 = IntegrationStrategyFactory::create("RK4");
    
    double t = 0.0;
    double dt = 0.001;
    int n_steps = 100;
    
    for (int i = 0; i < n_steps; ++i) {
        boris.step(ensemble_boris, idx, t, dt, forces);
        rk4->step(ensemble_rk4, idx, t, dt, forces);
        t += dt;
    }
    IonState ion_boris_out = ensemble_boris.to_legacy()[idx];
    IonState ion_rk4_out = ensemble_rk4.to_legacy()[idx];
    
    // Results should be similar (both 2nd order accurate, but different algorithms)
    REQUIRE_THAT(ion_boris_out.pos.x, WithinAbs(ion_rk4_out.pos.x, 1e-3));
    REQUIRE_THAT(ion_boris_out.pos.y, WithinAbs(ion_rk4_out.pos.y, 1e-3));
    REQUIRE_THAT(ion_boris_out.vel.x, WithinAbs(ion_rk4_out.vel.x, 1e-2));
    REQUIRE_THAT(ion_boris_out.vel.y, WithinAbs(ion_rk4_out.vel.y, 1e-2));
}
