// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/IForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

using namespace ICARION::physics;
using namespace ICARION::config;
using ICARION::core::IonEnsemble;
using Catch::Approx;

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Create minimal DomainConfig for testing
 */
static DomainConfig create_test_domain() {
    DomainConfig domain;
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.geometry.radius_m = 0.01;
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.environment.gas_species = "N2";
    domain.environment.compute_derived_properties();
    return domain;
}

// ============================================================================
// Mock Forces for Testing
// ============================================================================

/**
 * @brief Constant force (for testing)
 * 
 * Always returns the same force vector regardless of ion state.
 */
class ConstantForce : public IForce {
public:
    ConstantForce(const Vec3& force_vector) : force_(force_vector) {}
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)ensemble; (void)ion_idx; (void)t; (void)ctx;
        return force_;
    }
    
    std::string name() const override { return "ConstantForce"; }
    
private:
    Vec3 force_;
};

/**
 * @brief Conditional force (only applies to specific ions)
 * 
 * Returns force only if ion charge matches target charge.
 */
class ConditionalForce : public IForce {
public:
    ConditionalForce(double target_charge, const Vec3& force_vector)
        : target_charge_(target_charge), force_(force_vector) {}
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)ensemble; (void)ion_idx; (void)t; (void)ctx;
        return force_;
    }
    
    bool applies_to(const IonState& ion) const override {
        return std::abs(ion.ion_charge_C - target_charge_) < 1e-25;
    }
    
    std::string name() const override { return "ConditionalForce"; }
    
private:
    double target_charge_;
    Vec3 force_;
};

/**
 * @brief Gravity force (for testing realistic physics)
 * 
 * F = m * g (downward)
 */
class GravityForce : public IForce {
public:
    GravityForce(double g = 9.81) : g_(g) {}
    
    Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override {
        (void)t; (void)ctx;
        double m = ensemble.mass_data()[ion_idx];
        return Vec3{0, 0, -m * g_};
    }
    
    std::string name() const override { return "Gravity"; }
    
private:
    double g_;
};

// ============================================================================
// Unit Tests
// ============================================================================

TEST_CASE("ForceRegistry - Empty registry", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    SECTION("Empty registry returns zero force") {
        IonState ion;
        ion.mass_kg = 1e-26;
        ion.ion_charge_C = 1.602e-19;
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        REQUIRE(force.x == Approx(0.0));
        REQUIRE(force.y == Approx(0.0));
        REQUIRE(force.z == Approx(0.0));
    }
    
    SECTION("Empty registry properties") {
        REQUIRE(registry.empty());
        REQUIRE(registry.size() == 0);
        REQUIRE(registry.forces().size() == 0);
    }
}

TEST_CASE("ForceRegistry - Single force", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    // Add constant force: F = (1, 2, 3) N
    registry.add_force(std::make_unique<ConstantForce>(Vec3{1.0, 2.0, 3.0}));
    
    SECTION("Registry has one force") {
        REQUIRE_FALSE(registry.empty());
        REQUIRE(registry.size() == 1);
        REQUIRE(registry.forces().size() == 1);
    }
    
    SECTION("Computes correct force") {
        IonState ion;
        ion.mass_kg = 1e-26;
        ion.ion_charge_C = 1.602e-19;
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        REQUIRE(force.x == Approx(1.0));
        REQUIRE(force.y == Approx(2.0));
        REQUIRE(force.z == Approx(3.0));
    }
    
    SECTION("Force is time-independent") {
        IonState ion;
        ion.mass_kg = 1e-26;
        
        Vec3 f1 = registry.compute_total_force(ion, 0.0);
        Vec3 f2 = registry.compute_total_force(ion, 1.0);
        Vec3 f3 = registry.compute_total_force(ion, 100.0);
        
        REQUIRE(f1.x == f2.x);
        REQUIRE(f1.y == f2.y);
        REQUIRE(f1.z == f2.z);
        
        REQUIRE(f2.x == f3.x);
        REQUIRE(f2.y == f3.y);
        REQUIRE(f2.z == f3.z);
    }
}

TEST_CASE("ForceRegistry - Multiple forces", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    // Add three forces
    registry.add_force(std::make_unique<ConstantForce>(Vec3{1.0, 0.0, 0.0}));
    registry.add_force(std::make_unique<ConstantForce>(Vec3{0.0, 2.0, 0.0}));
    registry.add_force(std::make_unique<ConstantForce>(Vec3{0.0, 0.0, 3.0}));
    
    SECTION("Registry has three forces") {
        REQUIRE(registry.size() == 3);
    }
    
    SECTION("Forces are summed correctly") {
        IonState ion;
        ion.mass_kg = 1e-26;
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        // Should be (1, 2, 3) = sum of all forces
        REQUIRE(force.x == Approx(1.0));
        REQUIRE(force.y == Approx(2.0));
        REQUIRE(force.z == Approx(3.0));
    }
}

TEST_CASE("ForceRegistry - Conditional force", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    // Add conditional force (only applies to charge = 1.602e-19 C)
    double target_charge = 1.602e-19;
    registry.add_force(std::make_unique<ConditionalForce>(
        target_charge, 
        Vec3{10.0, 0.0, 0.0}
    ));
    
    SECTION("Force applies to matching ion") {
        IonState ion;
        ion.mass_kg = 1e-26;
        ion.ion_charge_C = target_charge;
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        REQUIRE(force.x == Approx(10.0));
        REQUIRE(force.y == Approx(0.0));
        REQUIRE(force.z == Approx(0.0));
    }
    
    SECTION("Force does not apply to non-matching ion") {
        IonState ion;
        ion.mass_kg = 1e-26;
        ion.ion_charge_C = 3.204e-19;  // Different charge
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        REQUIRE(force.x == Approx(0.0));
        REQUIRE(force.y == Approx(0.0));
        REQUIRE(force.z == Approx(0.0));
    }
}

TEST_CASE("ForceRegistry - Realistic physics (gravity)", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    // Earth gravity: g = 9.81 m/s^2
    registry.add_force(std::make_unique<GravityForce>(9.81));
    
    SECTION("Gravity force proportional to mass") {
        IonState ion1;
        ion1.mass_kg = 1.0;  // 1 kg
        
        IonState ion2;
        ion2.mass_kg = 2.0;  // 2 kg
        
        Vec3 f1 = registry.compute_total_force(ion1, 0.0);
        Vec3 f2 = registry.compute_total_force(ion2, 0.0);
        
        // F = m * g, so F2 should be twice F1
        REQUIRE(f2.z == Approx(2.0 * f1.z));
    }
    
    SECTION("Gravity force is downward (-z direction)") {
        IonState ion;
        ion.mass_kg = 1e-26;  // ~100 amu
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        REQUIRE(force.x == Approx(0.0));
        REQUIRE(force.y == Approx(0.0));
        REQUIRE(force.z < 0.0);  // Downward
        REQUIRE(force.z == Approx(-ion.mass_kg * 9.81));
    }
}

TEST_CASE("ForceRegistry - Clear functionality", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    // Add forces
    registry.add_force(std::make_unique<ConstantForce>(Vec3{1.0, 2.0, 3.0}));
    registry.add_force(std::make_unique<ConstantForce>(Vec3{4.0, 5.0, 6.0}));
    
    REQUIRE(registry.size() == 2);
    REQUIRE_FALSE(registry.empty());
    
    // Clear
    registry.clear();
    
    SECTION("Registry is empty after clear") {
        REQUIRE(registry.empty());
        REQUIRE(registry.size() == 0);
    }
    
    SECTION("Forces are zero after clear") {
        IonState ion;
        ion.mass_kg = 1e-26;
        
        Vec3 force = registry.compute_total_force(ion, 0.0);
        
        REQUIRE(force.x == Approx(0.0));
        REQUIRE(force.y == Approx(0.0));
        REQUIRE(force.z == Approx(0.0));
    }
}

TEST_CASE("ForceRegistry - Null force handling", "[forces][registry]") {
    ForceRegistry registry(create_test_domain());
    
    // Try to add null force (should be ignored)
    registry.add_force(nullptr);
    
    SECTION("Null force is not added") {
        REQUIRE(registry.empty());
        REQUIRE(registry.size() == 0);
    }
}

TEST_CASE("ForceRegistry - Force context (SSOT)", "[forces][registry]") {
    // Test that context is passed through correctly (SSOT-compliant)
    
    class ContextAwareForce : public IForce {
    public:
        Vec3 compute(const IonEnsemble& ensemble, size_t ion_idx, double t,
                     const ForceContext& ctx) const override {
            (void)ensemble; (void)ion_idx; (void)t;
            if (ctx.domain) {
                return Vec3{ctx.domain->environment.temperature_K, 0, 0};
            }
            return Vec3{0, 0, 0};
        }
        
        std::string name() const override { return "ContextAware"; }
    };
    
    ForceRegistry registry(create_test_domain());
    registry.add_force(std::make_unique<ContextAwareForce>());
    
    IonState ion;
    ion.mass_kg = 1e-26;
    
    SECTION("Context with domain (T = 300 K)") {
        DomainConfig domain;
        domain.environment.temperature_K = 300.0;
        
        ForceContext ctx;
        ctx.domain = &domain;
        
        Vec3 force = registry.compute_total_force(ion, 0.0, ctx);
        
        REQUIRE(force.x == Approx(300.0));
    }
    
    SECTION("Context with domain (T = 500 K)") {
        DomainConfig domain;
        domain.environment.temperature_K = 500.0;
        
        ForceContext ctx;
        ctx.domain = &domain;
        
        Vec3 force = registry.compute_total_force(ion, 0.0, ctx);
        
        REQUIRE(force.x == Approx(500.0));
    }
}
