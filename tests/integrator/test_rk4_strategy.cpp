// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/integrator/strategies/RK4Strategy.h"
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
// TEST FORCE: Constant Gravity (Free Fall)
// =============================================================================
class ConstantGravityForce : public IForce {
public:
    explicit ConstantGravityForce(double g = 9.81) : g_(g) {}
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        double m = ensemble.mass_data()[ion_idx];
        return Vec3{0, 0, -m * g_};
    }

    Vec3 compute_soa(const ForceState& state, double t,
                     const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        return Vec3{0, 0, -state.mass_kg * g_};
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
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        Vec3 pos = ensemble.get_pos(ion_idx);
        return pos * (-k_);
    }

    Vec3 compute_soa(const ForceState& state, double t,
                     const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        return state.pos * (-k_);
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
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        Vec3 vel = ensemble.get_vel(ion_idx);
        double mass = ensemble.mass_data()[ion_idx];
        return vel * (-gamma_ * mass);
    }

    Vec3 compute_soa(const ForceState& state, double t,
                     const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        return state.vel * (-gamma_ * state.mass_kg);
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
    DomainConfig domain = create_test_domain();
    ForceRegistry registry(domain);
    registry.add_force(std::make_unique<ConstantGravityForce>(9.81));
    
    // Initial state: ion at rest at height 100m
    IonState ion;
    ion.pos = Vec3{0, 0, 100.0};
    ion.vel = Vec3{0, 0, 0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;

    double t = 0.0;
    double dt = 0.001;  // smaller timestep for tighter parity
    
    // Integrate for 1 second
    for (int i = 0; i < 1000; ++i) {
        strategy.step(ensemble, idx, t, dt, registry);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // Analytical solution: z(t) = z0 - 0.5*g*t^2
    double t_final = 1.0;
    double z_expected = 100.0 - 0.5 * 9.81 * t_final * t_final;
    double v_expected = -9.81 * t_final;
    
    SECTION("Position after 1s") {
        REQUIRE_THAT(ion_out.pos.z, WithinAbs(z_expected, 0.01));  // 1cm error
    }
    
    SECTION("Velocity after 1s") {
        REQUIRE_THAT(ion_out.vel.z, WithinAbs(v_expected, 0.01));  // 1cm/s error
    }
}

TEST_CASE("RK4Strategy: Harmonic oscillator (periodic motion)", "[integrator][rk4]") {
    // Setup: mass-spring system with ω = 1 rad/s
    RK4Strategy strategy;
    DomainConfig domain = create_test_domain();
    ForceRegistry registry(domain);
    double k = 1.0;  // Spring constant
    registry.add_force(std::make_unique<HarmonicOscillatorForce>(k));
    
    // Initial state: displaced 1m, at rest
    IonState ion;
    ion.pos = Vec3{1.0, 0, 0};  // x0 = 1m
    ion.vel = Vec3{0, 0, 0};     // v0 = 0
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;
    
    double t = 0.0;
    double dt = 0.01;  // 10ms timestep
    double omega = std::sqrt(k / ion.mass_kg);  // ω = 1 rad/s
    
    // Integrate for one period (T = 2π/ω)
    double T = 2.0 * M_PI / omega;
    int steps = static_cast<int>(T / dt);
    
    for (int i = 0; i < steps; ++i) {
        strategy.step(ensemble, idx, t, dt, registry);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // After one period, should return to initial position
    SECTION("Position after one period") {
        REQUIRE_THAT(ion_out.pos.x, WithinAbs(1.0, 0.01));  // Return to x=1m
        REQUIRE_THAT(ion_out.vel.x, WithinAbs(0.0, 0.01));  // At rest again
    }
}

TEST_CASE("RK4Strategy: Exponential decay (damping)", "[integrator][rk4]") {
    // Setup: damped motion with γ = 0.5 s^-1
    RK4Strategy strategy;
    DomainConfig domain = create_test_domain();
    ForceRegistry registry(domain);
    double gamma = 0.5;
    registry.add_force(std::make_unique<LinearDampingForce>(gamma));
    
    // Initial state: moving at 10 m/s
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{10.0, 0, 0};  // v0 = 10 m/s
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    const size_t idx = 0;
    
    double t = 0.0;
    double dt = 0.01;
    
    // Integrate for 2 seconds
    for (int i = 0; i < 200; ++i) {
        strategy.step(ensemble, idx, t, dt, registry);
        t += dt;
    }
    IonState ion_out = ensemble.to_legacy()[idx];
    
    // Analytical solution: v(t) = v0 * exp(-γ*t)
    double t_final = 2.0;
    double v_expected = 10.0 * std::exp(-gamma * t_final);
    
    SECTION("Velocity after 2s (exponential decay)") {
        REQUIRE_THAT(ion_out.vel.x, WithinAbs(v_expected, 0.01));
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
    
    SECTION("RK45 implemented (Phase 4B)") {
        auto strategy = IntegrationStrategyFactory::create("RK45");
        REQUIRE(strategy != nullptr);
        REQUIRE(strategy->name() == "RK45");
        REQUIRE(strategy->is_adaptive() == true);
    }
    
    SECTION("Boris implemented (Phase 4B)") {
        auto strategy = IntegrationStrategyFactory::create("Boris");
        REQUIRE(strategy != nullptr);
        REQUIRE(strategy->name() == "Boris");
        REQUIRE(strategy->is_adaptive() == false);
    }
}

TEST_CASE("RK4Strategy: SSOT compliance", "[integrator][ssot]") {
    // Verify RK4 uses DomainConfig (not GlobalParams)
    RK4Strategy strategy;
    // Create DomainConfig (SSOT!)
    DomainConfig domain = create_test_domain();
    ForceRegistry registry(domain);
    registry.add_force(std::make_unique<ConstantGravityForce>());
    
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{0, 0, 0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1e-19;
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
    
    // This should compile and run (SSOT-compliant signature)
    REQUIRE_NOTHROW(strategy.step(ensemble, 0, 0.0, 0.01, registry));
}
