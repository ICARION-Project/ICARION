// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   @file       test_rk4_strategy.cpp
 *   @brief      Unit tests for RK4Strategy
 *
 *   @details
 *   Tests RK4 integration against analytical solutions:
 *   - Free fall (constant acceleration)
 *   - Harmonic oscillator (periodic motion)
 *   - Exponential decay (damping)
 *
 *   Validates SSOT compliance (uses DomainConfig, ForceRegistry).
 *
 * =====================================================================
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/IntegrationStrategyFactory.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/IForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

#include <cmath>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;
using Catch::Matchers::WithinAbs;

// =============================================================================
// TEST FORCE: Constant Gravity (Free Fall)
// =============================================================================
class ConstantGravityForce : public IForce {
public:
    explicit ConstantGravityForce(double g = 9.81) : g_(g) {}
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override {
        (void)t; (void)ctx;  // Unused
        return Vec3{0, 0, -ion.mass_kg * g_};  // Downward force
    }
    
    std::string name() const override { return "ConstantGravity"; }
    
private:
    double g_;
};

// =============================================================================
// TEST FORCE: Harmonic Oscillator (Spring)
// =============================================================================
class HarmonicOscillatorForce : public IForce {
public:
    explicit HarmonicOscillatorForce(double k = 1.0) : k_(k) {}
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override {
        (void)t; (void)ctx;  // Unused
        // F = -k*x (Hooke's law)
        return ion.pos * (-k_);
    }
    
    std::string name() const override { return "HarmonicOscillator"; }
    
private:
    double k_;
};

// =============================================================================
// TEST FORCE: Linear Damping
// =============================================================================
class LinearDampingForce : public IForce {
public:
    explicit LinearDampingForce(double gamma = 0.1) : gamma_(gamma) {}
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override {
        (void)t; (void)ctx;  // Unused
        // F = -γ*m*v
        return ion.vel * (-gamma_ * ion.mass_kg);
    }
    
    std::string name() const override { return "LinearDamping"; }
    
private:
    double gamma_;
};

// =============================================================================
// HELPER: Create minimal DomainConfig
// =============================================================================
DomainConfig create_test_domain() {
    DomainConfig domain;
    domain.geometry.length_m = 1.0;
    domain.geometry.radius_m = 0.1;
    domain.geometry.min_bound = Vec3{-0.1, -0.1, 0.0};
    domain.geometry.max_bound = Vec3{0.1, 0.1, 1.0};
    return domain;
}

// =============================================================================
// TEST SUITE: RK4Strategy
// =============================================================================

TEST_CASE("RK4Strategy: Basic properties", "[integrator][rk4]") {
    RK4Strategy strategy;
    
    SECTION("Name") {
        REQUIRE(strategy.name() == "RK4");
    }
    
    SECTION("Not adaptive") {
        REQUIRE(strategy.is_adaptive() == false);
    }
}

TEST_CASE("RK4Strategy: Free fall (constant acceleration)", "[integrator][rk4]") {
    // Setup
    RK4Strategy strategy;
    ForceRegistry registry;
    registry.add_force(std::make_unique<ConstantGravityForce>(9.81));
    
    DomainConfig domain = create_test_domain();
    std::vector<IonState> all_ions;  // No space charge
    
    // Initial state: ion at rest at height 100m
    IonState ion;
    ion.pos = Vec3{0, 0, 100.0};
    ion.vel = Vec3{0, 0, 0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    
    double t = 0.0;
    double dt = 0.01;  // 10ms timestep
    
    // Integrate for 1 second
    for (int i = 0; i < 100; ++i) {
        strategy.step(ion, t, dt, registry, domain, all_ions);
        t += dt;
    }
    
    // Analytical solution: z(t) = z0 - 0.5*g*t^2
    double t_final = 1.0;
    double z_expected = 100.0 - 0.5 * 9.81 * t_final * t_final;
    double v_expected = -9.81 * t_final;
    
    SECTION("Position after 1s") {
        REQUIRE_THAT(ion.pos.z, WithinAbs(z_expected, 0.01));  // 1cm error
    }
    
    SECTION("Velocity after 1s") {
        REQUIRE_THAT(ion.vel.z, WithinAbs(v_expected, 0.01));  // 1cm/s error
    }
}

TEST_CASE("RK4Strategy: Harmonic oscillator (periodic motion)", "[integrator][rk4]") {
    // Setup: mass-spring system with ω = 1 rad/s
    RK4Strategy strategy;
    ForceRegistry registry;
    double k = 1.0;  // Spring constant
    registry.add_force(std::make_unique<HarmonicOscillatorForce>(k));
    
    DomainConfig domain = create_test_domain();
    std::vector<IonState> all_ions;
    
    // Initial state: displaced 1m, at rest
    IonState ion;
    ion.pos = Vec3{1.0, 0, 0};  // x0 = 1m
    ion.vel = Vec3{0, 0, 0};     // v0 = 0
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    
    double t = 0.0;
    double dt = 0.01;  // 10ms timestep
    double omega = std::sqrt(k / ion.mass_kg);  // ω = 1 rad/s
    
    // Integrate for one period (T = 2π/ω)
    double T = 2.0 * M_PI / omega;
    int steps = static_cast<int>(T / dt);
    
    for (int i = 0; i < steps; ++i) {
        strategy.step(ion, t, dt, registry, domain, all_ions);
        t += dt;
    }
    
    // After one period, should return to initial position
    SECTION("Position after one period") {
        REQUIRE_THAT(ion.pos.x, WithinAbs(1.0, 0.01));  // Return to x=1m
        REQUIRE_THAT(ion.vel.x, WithinAbs(0.0, 0.01));  // At rest again
    }
}

TEST_CASE("RK4Strategy: Exponential decay (damping)", "[integrator][rk4]") {
    // Setup: damped motion with γ = 0.5 s^-1
    RK4Strategy strategy;
    ForceRegistry registry;
    double gamma = 0.5;
    registry.add_force(std::make_unique<LinearDampingForce>(gamma));
    
    DomainConfig domain = create_test_domain();
    std::vector<IonState> all_ions;
    
    // Initial state: moving at 10 m/s
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{10.0, 0, 0};  // v0 = 10 m/s
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    
    double t = 0.0;
    double dt = 0.01;
    
    // Integrate for 2 seconds
    for (int i = 0; i < 200; ++i) {
        strategy.step(ion, t, dt, registry, domain, all_ions);
        t += dt;
    }
    
    // Analytical solution: v(t) = v0 * exp(-γ*t)
    double t_final = 2.0;
    double v_expected = 10.0 * std::exp(-gamma * t_final);
    
    SECTION("Velocity after 2s (exponential decay)") {
        REQUIRE_THAT(ion.vel.x, WithinAbs(v_expected, 0.01));
    }
}

TEST_CASE("RK4Strategy: Factory creation", "[integrator][factory]") {
    SECTION("Create RK4 via factory") {
        auto strategy = IntegrationStrategyFactory::create("RK4");
        REQUIRE(strategy != nullptr);
        REQUIRE(strategy->name() == "RK4");
        REQUIRE(strategy->is_adaptive() == false);
    }
    
    SECTION("Unknown strategy throws exception") {
        REQUIRE_THROWS_AS(
            IntegrationStrategyFactory::create("InvalidStrategy"),
            std::invalid_argument
        );
    }
    
    SECTION("RK45 not yet implemented") {
        REQUIRE_THROWS_AS(
            IntegrationStrategyFactory::create("RK45"),
            std::invalid_argument
        );
    }
}

TEST_CASE("RK4Strategy: SSOT compliance", "[integrator][ssot]") {
    // Verify RK4 uses DomainConfig (not GlobalParams)
    RK4Strategy strategy;
    ForceRegistry registry;
    registry.add_force(std::make_unique<ConstantGravityForce>());
    
    // Create DomainConfig (SSOT!)
    DomainConfig domain = create_test_domain();
    std::vector<IonState> all_ions;
    
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{0, 0, 0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    
    // This should compile and run (SSOT-compliant signature)
    REQUIRE_NOTHROW(strategy.step(ion, 0.0, 0.01, registry, domain, all_ions));
}
