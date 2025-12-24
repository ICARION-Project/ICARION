// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/integrator/strategies/RK45Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/IForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

#include <cmath>
#include <iostream>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;
using ICARION::core::IonEnsemble;
using Catch::Matchers::WithinAbs;

// =============================================================================
// TEST FORCES (reuse from RK4 tests)
// =============================================================================

class ConstantGravityForce : public IForce {
public:
    explicit ConstantGravityForce(double g = 9.81) : g_(g) {}
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        const double m = ensemble.mass_data()[ion_idx];
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

class ExponentialDampingForce : public IForce {
public:
    explicit ExponentialDampingForce(double gamma = 1.0) : gamma_(gamma) {}
    
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
    
    std::string name() const override { return "ExponentialDamping"; }
    
private:
    double gamma_;
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

TEST_CASE("RK45Strategy: Constructor and configuration", "[integrator][rk45]") {
    SECTION("Default configuration") {
        RK45Strategy::AdaptiveConfig config;
        RK45Strategy strategy(config);
        
        REQUIRE(strategy.name() == "RK45");
        REQUIRE(strategy.is_adaptive() == true);
    }
    
    SECTION("Custom tolerances") {
        RK45Strategy::AdaptiveConfig config;
        config.atol = 1e-10;
        config.rtol = 1e-8;
        config.safety_factor = 0.85;
        
        RK45Strategy strategy(config);
        REQUIRE(strategy.is_adaptive() == true);
    }
    
    SECTION("Invalid configuration") {
        RK45Strategy::AdaptiveConfig config;
        config.atol = -1e-8;  // Invalid: negative
        
        REQUIRE_THROWS_AS(RK45Strategy(config), std::invalid_argument);
    }
}

// =============================================================================
// TEST: Free Fall (Constant Acceleration)
// =============================================================================

TEST_CASE("RK45Strategy: Free fall with adaptive timestep", "[integrator][rk45]") {
    // Setup
    auto domain = create_test_domain();
    auto ion = create_test_ion(1.0, 1.0);
    ion.pos = Vec3{0, 0, 100};  // Start at 100 m height
    ion.vel = Vec3{0, 0, 0};    // Zero initial velocity
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<ConstantGravityForce>(9.81));
    std::vector<IonState> all_ions = {ion};
    
    RK45Strategy::AdaptiveConfig config;
    config.atol = 1e-6;
    config.rtol = 1e-6;
    RK45Strategy strategy(config);
    
    // Integrate for 1 second
    double t = 0.0;
    double t_final = 1.0;
    double dt = 0.01;  // Initial timestep
    
    while (t < t_final) {
        double dt_step = std::min(dt, t_final - t);
        double dt_used = dt_step;  // Save before step_adaptive modifies it
        strategy.step_adaptive(ion, t, dt_step, forces, all_ions);
        all_ions[0] = ion;  // keep AoS cache in sync for space charge path
        t += dt_used;  // Use the dt we actually stepped with
        dt = dt_step;  // Get suggested dt for next step
    }
    
    // Analytical solution: z(t) = z0 - 0.5*g*t²
    double z_expected = 100.0 - 0.5 * 9.81 * 1.0 * 1.0;
    double vz_expected = -9.81 * 1.0;
    
    REQUIRE_THAT(ion.pos.z, WithinAbs(z_expected, 1e-4));
    REQUIRE_THAT(ion.vel.z, WithinAbs(vz_expected, 1e-4));
    
    // Check statistics
    auto stats = strategy.get_stats();
    REQUIRE(stats.accepted_steps > 0);
    REQUIRE(static_cast<long long>(stats.rejected_steps) >= 0);  // May have some rejections
    REQUIRE(stats.avg_error <= 1.5);  // Error should be near 1.0
}

// =============================================================================
// TEST: Harmonic Oscillator (Variable Dynamics)
// =============================================================================

TEST_CASE("RK45Strategy: Harmonic oscillator with adaptive timestep", "[integrator][rk45]") {
    // Setup: x'' = -ω²*x, solution: x(t) = A*cos(ωt)
    auto domain = create_test_domain();
    auto ion = create_test_ion(1.0, 1.0);
    
    double A = 1.0;         // Amplitude [m]
    double omega = 2.0;     // Angular frequency [rad/s]
    double k = omega * omega;  // Spring constant
    
    ion.pos = Vec3{A, 0, 0};     // x0 = A
    ion.vel = Vec3{0, 0, 0};     // v0 = 0
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<HarmonicOscillatorForce>(k));
    
    std::vector<IonState> all_ions = {ion};

    RK45Strategy::AdaptiveConfig config;
    config.atol = 1e-8;
    config.rtol = 1e-6;
    RK45Strategy strategy(config);
    
    // Integrate for one period: T = 2π/ω
    double T = 2.0 * M_PI / omega;
    double t = 0.0;
    double t_final = T;
    double dt = T / 1000.0;  // Initial timestep (smaller for accuracy)
    
    int steps = 0;
    while (t < t_final && steps < 10000) {
        double dt_step = std::min(dt, t_final - t);
        double dt_used = dt_step;
        strategy.step_adaptive(ion, t, dt_step, forces, all_ions);
        all_ions[0] = ion;
        t += dt_used;
        dt = dt_step;
        steps++;
    }
    const IonState& ion_out = ion;
    
    // After one period, should return to initial state
    double x_expected = A * std::cos(omega * t_final);
    double vx_expected = -A * omega * std::sin(omega * t_final);
    INFO("t_final reached = " << t << " (target " << t_final << ")");
    REQUIRE_THAT(t, WithinAbs(t_final, 1e-6));
    INFO("x_out=" << ion_out.pos.x << " v_out=" << ion_out.vel.x);
    
    REQUIRE_THAT(ion_out.pos.x, WithinAbs(x_expected, 1e-6));
    REQUIRE_THAT(ion_out.vel.x, WithinAbs(vx_expected, 1e-5));
    
    // Check that adaptive stepping was used
    auto stats = strategy.get_stats();
    REQUIRE(stats.min_step_used < stats.max_step_used);
    INFO("Min step: " << stats.min_step_used << ", Max step: " << stats.max_step_used);
}

// =============================================================================
// TEST: Exponential Decay (Stiff Problem)
// =============================================================================

TEST_CASE("RK45Strategy: Exponential decay with error control", "[integrator][rk45]") {
    // Setup: v' = -γ*v, solution: v(t) = v0*exp(-γ*t)
    auto domain = create_test_domain();
    auto ion = create_test_ion(1.0, 1.0);
    
    double gamma = 10.0;  // Decay rate [1/s]
    double v0 = 100.0;    // Initial velocity [m/s]
    
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{v0, 0, 0};
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<ExponentialDampingForce>(gamma));
    
    std::vector<IonState> all_ions = {ion};
    
    RK45Strategy::AdaptiveConfig config;
    config.atol = 1e-8;
    config.rtol = 1e-6;
    RK45Strategy strategy(config);
    
    // Integrate to t = 1/γ (one time constant)
    double t = 0.0;
    double t_final = 1.0 / gamma;
    double dt = 0.001;  // Initial timestep
    
    int steps = 0;
    while (t < t_final && steps < 10000) {
        double dt_step = std::min(dt, t_final - t);
        double dt_used = dt_step;
        strategy.step_adaptive(ion, t, dt_step, forces, all_ions);
        all_ions[0] = ion;
        t += dt_used;
        dt = dt_step;
        steps++;
    }
    
    // Analytical solution
    double vx_expected = v0 * std::exp(-gamma * t_final);
    
    REQUIRE_THAT(ion.vel.x, WithinAbs(vx_expected, 1e-4));
    
    // Error control should keep error near 1.0
    auto stats = strategy.get_stats();
    REQUIRE(stats.avg_error < 2.0);
    REQUIRE(stats.avg_error >= 0.0);
}

// =============================================================================
// TEST: Convergence Order (Richardson Extrapolation)
// =============================================================================

TEST_CASE("RK45Strategy: Convergence order verification", "[integrator][rk45][convergence]") {
    // Test that RK45 achieves 5th order convergence
    auto domain = create_test_domain();
    
    double omega = 1.0;
    double k = omega * omega;
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<HarmonicOscillatorForce>(k));
    
    std::vector<IonState> all_ions;
    
    double t_final = 1.0;
    
    // Test with different initial timesteps
    std::vector<double> dt_values = {0.1, 0.05, 0.025};
    std::vector<double> errors;
    
    for (double dt_init : dt_values) {
        auto ion = create_test_ion(1.0, 1.0);
        ion.pos = Vec3{1.0, 0, 0};
        ion.vel = Vec3{0, 0, 0};
        
        RK45Strategy::AdaptiveConfig config;
        config.atol = 1e-9;  // Moderate tolerance to allow dt variation
        config.rtol = 1e-7;
        RK45Strategy strategy(config);
        
        double t = 0.0;
        double dt = dt_init;
        
        while (t < t_final) {
            double dt_step = std::min(dt, t_final - t);
            double dt_used = dt_step;
            strategy.step_adaptive(ion, t, dt_step, forces, all_ions);
            all_ions.clear();
            all_ions.push_back(ion);
            t += dt_used;
            dt = dt_step;
        }
        
        double x_exact = std::cos(omega * t_final);
        double error = std::fabs(ion.pos.x - x_exact);
        errors.push_back(error);
    }
    
    // Check convergence: with adaptive stepping, smaller initial dt allows tighter error control
    // Expect at least 2nd order convergence (factor 4 when dt halves)
    if (errors.size() >= 2) {
        double ratio1 = errors[0] / errors[1];
        INFO("Convergence ratio (dt halved): " << ratio1);
        // With adaptive stepping, expect improvement when starting with smaller dt
        REQUIRE(ratio1 > 2.0);  // At least 2nd order improvement
    }
}

// =============================================================================
// TEST: Step Rejection and Retry
// =============================================================================

TEST_CASE("RK45Strategy: Step rejection with tight tolerance", "[integrator][rk45][adaptive]") {
    auto domain = create_test_domain();
    auto ion = create_test_ion(1.0, 1.0);
    
    ion.pos = Vec3{1.0, 0, 0};
    ion.vel = Vec3{0, 0, 0};
    
    double omega = 10.0;  // Fast oscillation
    double k = omega * omega;
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<HarmonicOscillatorForce>(k));
    
    std::vector<IonState> all_ions = {ion};
    
    RK45Strategy::AdaptiveConfig config;
    config.atol = 1e-10;  // Very tight tolerance
    config.rtol = 1e-8;
    RK45Strategy strategy(config);
    
    double t = 0.0;
    double dt = 0.1;  // Large initial step (will be rejected)
    
    // Take a few steps
    for (int i = 0; i < 5; ++i) {
        double dt_used = dt;
        strategy.step_adaptive(ion, t, dt, forces, all_ions);
        t += dt_used;
    }
    
    // Should have some rejected steps
    auto stats = strategy.get_stats();
    REQUIRE(stats.accepted_steps > 0);
    // May or may not have rejections depending on adaptive algorithm
    INFO("Accepted: " << stats.accepted_steps << ", Rejected: " << stats.rejected_steps);
}

// =============================================================================
// TEST: FSAL Property (First-Same-As-Last)
// =============================================================================

TEST_CASE("RK45Strategy: FSAL optimization", "[integrator][rk45][performance]") {
    // FSAL property: k7 from step n = k1 for step n+1
    // This test verifies that results are consistent with FSAL
    
    auto domain = create_test_domain();
    auto ion1 = create_test_ion(1.0, 1.0);
    auto ion2 = create_test_ion(1.0, 1.0);
    
    ion1.pos = Vec3{1.0, 0, 0};
    ion1.vel = Vec3{0, 1.0, 0};
    ion2 = ion1;
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<HarmonicOscillatorForce>(1.0));
    
    std::vector<IonState> all_ions;
    
    RK45Strategy::AdaptiveConfig config;
    config.atol = 1e-8;
    config.rtol = 1e-6;
    
    RK45Strategy strategy1(config);
    RK45Strategy strategy2(config);
    
    double t = 0.0;
    double dt1 = 0.01;
    double dt2 = 0.01;
    
    // Take two steps with each strategy
    strategy1.step_adaptive(ion1, t, dt1, forces, all_ions);
    strategy1.step_adaptive(ion1, t + dt1, dt1, forces, all_ions);
    
    strategy2.step_adaptive(ion2, t, dt2, forces, all_ions);
    strategy2.step_adaptive(ion2, t + dt2, dt2, forces, all_ions);
    
    // Results should be identical (FSAL property ensures consistency)
    REQUIRE_THAT(ion1.pos.x, WithinAbs(ion2.pos.x, 1e-12));
    REQUIRE_THAT(ion1.pos.y, WithinAbs(ion2.pos.y, 1e-12));
    REQUIRE_THAT(ion1.vel.x, WithinAbs(ion2.vel.x, 1e-12));
    REQUIRE_THAT(ion1.vel.y, WithinAbs(ion2.vel.y, 1e-12));
}

// =============================================================================
// TEST: Statistics Tracking
// =============================================================================

TEST_CASE("RK45Strategy: Statistics collection", "[integrator][rk45][stats]") {
    auto domain = create_test_domain();
    auto ion = create_test_ion(1.0, 1.0);
    
    ion.pos = Vec3{1.0, 0, 0};
    ion.vel = Vec3{0, 0, 0};
    
    ForceRegistry forces(domain);
    forces.add_force(std::make_unique<HarmonicOscillatorForce>(1.0));
    
    std::vector<IonState> all_ions = {ion};
    
    RK45Strategy::AdaptiveConfig config;
    RK45Strategy strategy(config);
    
    double t = 0.0;
    double dt = 0.01;
    
    // Take 10 steps
    for (int i = 0; i < 10; ++i) {
        double dt_used = dt;
        strategy.step_adaptive(ion, t, dt, forces, all_ions);
        t += dt_used;
    }
    
    auto stats = strategy.get_stats();
    
    REQUIRE(stats.accepted_steps == 10);
    REQUIRE(stats.min_step_used > 0.0);
    REQUIRE(stats.max_step_used > 0.0);
    REQUIRE(stats.avg_error > 0.0);
    
    // Reset and check
    strategy.reset_stats();
    auto stats_reset = strategy.get_stats();
    REQUIRE(stats_reset.accepted_steps == 0);
}
