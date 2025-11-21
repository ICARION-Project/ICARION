// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_collision_energy_conservation.cpp
 * @brief Unit tests for energy conservation in collision models
 * 
 * Tests that elastic collision models (HSS, EHSS) conserve total energy
 * (ion + neutral) in the collision process.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/physics/collisions/collisionHelpers.h"
#include "utils/constants.h"
#include <cmath>

using namespace ICARION;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// Helper: Calculate kinetic energy
double kinetic_energy_J(const Vec3& vel, double mass_kg) {
    double v2 = vel.x * vel.x + vel.y * vel.y + vel.z * vel.z;
    return 0.5 * mass_kg * v2;
}

TEST_CASE("Hard-sphere collision conserves total energy", "[collision][energy][conservation]") {
    // Setup collision parameters
    EHSSParams p;
    p.mi = 19.0 * AMU_TO_KG;  // H3O+
    p.mn = 4.0 * AMU_TO_KG;    // He
    p.kB = BOLTZMANN_CONSTANT;
    p.Tn = 300.0;
    p.ubx = 0.0;
    p.uby = 0.0;
    p.ubz = 0.0;
    p.sigma_eff = 25e-20;  // m^2
    p.Rn = std::sqrt(p.sigma_eff / M_PI);
    p.n = 1e20;  // m^-3
    p.dt = 1e-9;
    
    EhssRng rng(42);
    
    SECTION("Single collision with thermal neutral") {
        // Ion moving fast
        Vec3 v_ion_before{1000.0, 0.0, 0.0};
        
        // Sample neutral from thermal distribution
        Vec3 v_neutral_before = sample_neutral_velocity(p, rng);
        
        // Calculate energy before collision
        double E_before = kinetic_energy_J(v_ion_before, p.mi) + 
                         kinetic_energy_J(v_neutral_before, p.mn);
        
        // Apply collision
        Vec3 v_ion_after = collide_hs_cpu(v_ion_before, v_neutral_before, p, rng);
        
        // Calculate neutral velocity after (from momentum conservation)
        Vec3 momentum_before = v_ion_before * p.mi + v_neutral_before * p.mn;
        Vec3 momentum_after_ion = v_ion_after * p.mi;
        Vec3 v_neutral_after = (momentum_before - momentum_after_ion) * (1.0 / p.mn);
        
        // Calculate energy after collision
        double E_after = kinetic_energy_J(v_ion_after, p.mi) + 
                        kinetic_energy_J(v_neutral_after, p.mn);
        
        // Energy should be conserved to machine precision
        REQUIRE_THAT(E_after, WithinRel(E_before, 1e-10));
    }
    
    SECTION("Multiple collisions - energy conservation") {
        int N_trials = 1000;
        double max_rel_error = 0.0;
        
        for (int i = 0; i < N_trials; ++i) {
            // Random initial velocities
            Vec3 v_ion{rng.uniform01() * 2000 - 1000, 
                      rng.uniform01() * 2000 - 1000,
                      rng.uniform01() * 2000 - 1000};
            Vec3 v_neutral = sample_neutral_velocity(p, rng);
            
            double E_before = kinetic_energy_J(v_ion, p.mi) + 
                             kinetic_energy_J(v_neutral, p.mn);
            
            Vec3 v_ion_after = collide_hs_cpu(v_ion, v_neutral, p, rng);
            
            Vec3 momentum = v_ion * p.mi + v_neutral * p.mn;
            Vec3 v_neutral_after = (momentum - v_ion_after * p.mi) * (1.0 / p.mn);
            
            double E_after = kinetic_energy_J(v_ion_after, p.mi) + 
                            kinetic_energy_J(v_neutral_after, p.mn);
            
            double rel_error = std::abs((E_after - E_before) / E_before);
            max_rel_error = std::max(max_rel_error, rel_error);
        }
        
        INFO("Max relative energy error over " << N_trials << " trials: " << max_rel_error);
        REQUIRE(max_rel_error < 1e-9);
    }
    
    SECTION("Momentum conservation check") {
        Vec3 v_ion{500.0, -300.0, 200.0};
        Vec3 v_neutral{100.0, 50.0, -75.0};
        
        Vec3 p_before = v_ion * p.mi + v_neutral * p.mn;
        
        Vec3 v_ion_after = collide_hs_cpu(v_ion, v_neutral, p, rng);
        Vec3 v_neutral_after = (p_before - v_ion_after * p.mi) * (1.0 / p.mn);
        
        Vec3 p_after = v_ion_after * p.mi + v_neutral_after * p.mn;
        
        REQUIRE_THAT(p_after.x, WithinRel(p_before.x, 1e-10));
        REQUIRE_THAT(p_after.y, WithinRel(p_before.y, 1e-10));
        REQUIRE_THAT(p_after.z, WithinRel(p_before.z, 1e-10));
    }
}

TEST_CASE("EHSS collision with geometry conserves energy", "[collision][energy][ehss]") {
    // Setup for H3O+ with simple geometry
    EHSSParams p;
    p.mi = 19.0 * AMU_TO_KG;
    p.mn = 4.0 * AMU_TO_KG;
    p.kB = BOLTZMANN_CONSTANT;
    p.Tn = 300.0;
    p.ubx = 0.0;
    p.uby = 0.0;
    p.ubz = 0.0;
    p.sigma_eff = 25e-20;
    p.Rn = 1.3e-10;
    p.n = 1e20;
    p.dt = 1e-9;
    
    // Simple 4-atom geometry (H3O+)
    std::vector<Vec3> centers = {
        Vec3{0.0, 0.0, 0.0},        // O
        Vec3{0.96e-10, 0.0, 0.0},   // H1
        Vec3{-0.24e-10, 0.93e-10, 0.0},  // H2
        Vec3{-0.24e-10, -0.93e-10, 0.0}  // H3
    };
    std::vector<double> radii = {1.2e-10, 1.1e-10, 1.1e-10, 1.1e-10};
    
    EhssRng rng(123);
    
    SECTION("EHSS geometry collision conserves energy") {
        int N_trials = 100;
        double max_rel_error = 0.0;
        
        for (int i = 0; i < N_trials; ++i) {
            Vec3 v_ion{rng.uniform01() * 1000, 
                      rng.uniform01() * 1000,
                      rng.uniform01() * 1000};
            Vec3 v_neutral = sample_neutral_velocity(p, rng);
            
            double E_before = kinetic_energy_J(v_ion, p.mi) + 
                             kinetic_energy_J(v_neutral, p.mn);
            
            Vec3 v_ion_after = collide_ehss_cpu_geometry_given_neutral(
                v_ion, v_neutral, p, centers, radii, rng
            );
            
            Vec3 momentum = v_ion * p.mi + v_neutral * p.mn;
            Vec3 v_neutral_after = (momentum - v_ion_after * p.mi) * (1.0 / p.mn);
            
            double E_after = kinetic_energy_J(v_ion_after, p.mi) + 
                            kinetic_energy_J(v_neutral_after, p.mn);
            
            double rel_error = std::abs((E_after - E_before) / E_before);
            max_rel_error = std::max(max_rel_error, rel_error);
        }
        
        INFO("Max relative energy error (EHSS): " << max_rel_error);
        REQUIRE(max_rel_error < 1e-9);
    }
}

TEST_CASE("Center-of-mass frame energy conservation", "[collision][energy][com]") {
    double m_ion = 19.0 * AMU_TO_KG;
    double m_neutral = 4.0 * AMU_TO_KG;
    
    Vec3 v_ion{800.0, -200.0, 300.0};
    Vec3 v_neutral{50.0, 100.0, -50.0};
    
    // Calculate COM velocity
    double M = m_ion + m_neutral;
    Vec3 v_com = (v_ion * m_ion + v_neutral * m_neutral) * (1.0 / M);
    
    // Velocities in COM frame
    Vec3 v_ion_com = v_ion - v_com;
    Vec3 v_neutral_com = v_neutral - v_com;
    
    SECTION("COM frame has zero total momentum") {
        Vec3 p_com = v_ion_com * m_ion + v_neutral_com * m_neutral;
        
        REQUIRE_THAT(p_com.x, WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(p_com.y, WithinAbs(0.0, 1e-15));
        REQUIRE_THAT(p_com.z, WithinAbs(0.0, 1e-15));
    }
    
    SECTION("Relative kinetic energy in COM frame") {
        Vec3 v_rel = v_ion - v_neutral;
        double mu = (m_ion * m_neutral) / M;  // Reduced mass
        
        // E_rel = (1/2) * mu * v_rel^2
        double E_rel_direct = 0.5 * mu * (v_rel.x*v_rel.x + v_rel.y*v_rel.y + v_rel.z*v_rel.z);
        
        // E_rel = (1/2)*m_i*v_i_com^2 + (1/2)*m_n*v_n_com^2
        double E_rel_com = kinetic_energy_J(v_ion_com, m_ion) + 
                          kinetic_energy_J(v_neutral_com, m_neutral);
        
        REQUIRE_THAT(E_rel_com, WithinRel(E_rel_direct, 1e-10));
    }
}

TEST_CASE("Thermalization to correct equilibrium validates collision physics", "[collision][energy][thermalization]") {
    // CRITICAL DEBUG TEST: Validates that collision physics produces correct thermal equilibrium
    // 
    // This test is essential for catching bugs in:
    // - Collision probability calculation
    // - Neutral velocity sampling
    // - Energy transfer implementation
    // - Handler-level physics integration
    //
    // Theory: Ions colliding many times with Maxwell-Boltzmann neutrals at temperature T
    // should equilibrate to <E_ion> = (3/2)kT regardless of initial conditions.
    //
    // This test uses ensemble averaging (many ions) and validates that the equilibrium
    // energy matches thermal expectation to within ~10%.
    
    SECTION("H3O+ thermalizes to 300K in He buffer gas") {
        EHSSParams p;
        p.mi = 19.0 * AMU_TO_KG;
        p.mn = 4.0 * AMU_TO_KG;
        p.kB = BOLTZMANN_CONSTANT;
        p.Tn = 300.0;
        p.ubx = 0.0;
        p.uby = 0.0;
        p.ubz = 0.0;
        p.sigma_eff = 25e-20;
        p.Rn = std::sqrt(p.sigma_eff / M_PI);
        
        const double E_thermal = 1.5 * p.kB * p.Tn;  // (3/2)kT
        const int N_IONS = 500;
        const int N_COLLISIONS = 1000;
        
        double sum_E_final = 0.0;
        
        #pragma omp parallel for reduction(+:sum_E_final)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            Vec3 v_ion{3000.0, 0.0, 0.0};  // Start with high energy
            EhssRng rng(1000 + ion_idx);
            
            // Apply many collisions with properly sampled thermal neutrals
            for (int i = 0; i < N_COLLISIONS; ++i) {
                Vec3 v_neutral = sample_neutral_velocity(p, rng);
                v_ion = collide_hs_cpu(v_ion, v_neutral, p, rng);
            }
            
            double E_final = kinetic_energy_J(v_ion, p.mi);
            sum_E_final += E_final;
        }
        
        double avg_E_final = sum_E_final / N_IONS;
        double ratio = avg_E_final / E_thermal;
        
        INFO("Expected thermal energy: " << E_thermal / ELEM_CHARGE_C << " eV");
        INFO("Final average energy: " << avg_E_final / ELEM_CHARGE_C << " eV");
        INFO("Thermalization ratio: " << ratio);
        
        // CRITICAL: Should thermalize to within 10% of (3/2)kT
        // If this fails, there's likely a bug in collision physics or neutral sampling!
        REQUIRE_THAT(ratio, WithinRel(1.0, 0.10));
    }
    
    SECTION("He+ thermalizes to 300K in N2 buffer gas") {
        EHSSParams p;
        p.mi = 4.0 * AMU_TO_KG;
        p.mn = 28.0 * AMU_TO_KG;
        p.kB = BOLTZMANN_CONSTANT;
        p.Tn = 300.0;
        p.ubx = 0.0;
        p.uby = 0.0;
        p.ubz = 0.0;
        p.sigma_eff = 25e-20;
        p.Rn = std::sqrt(p.sigma_eff / M_PI);
        
        const double E_thermal = 1.5 * p.kB * p.Tn;
        const int N_IONS = 500;
        const int N_COLLISIONS = 1000;
        
        double sum_E_final = 0.0;
        
        #pragma omp parallel for reduction(+:sum_E_final)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            Vec3 v_ion{2500.0, 0.0, 0.0};
            EhssRng rng(2000 + ion_idx);
            
            for (int i = 0; i < N_COLLISIONS; ++i) {
                Vec3 v_neutral = sample_neutral_velocity(p, rng);
                v_ion = collide_hs_cpu(v_ion, v_neutral, p, rng);
            }
            
            double E_final = kinetic_energy_J(v_ion, p.mi);
            sum_E_final += E_final;
        }
        
        double avg_E_final = sum_E_final / N_IONS;
        double ratio = avg_E_final / E_thermal;
        
        INFO("He+ in N2 - Thermalization ratio: " << ratio);
        REQUIRE_THAT(ratio, WithinRel(1.0, 0.10));
    }
    
    SECTION("Equal mass thermalization (m_ion = m_neutral)") {
        EHSSParams p;
        p.mi = 20.0 * AMU_TO_KG;
        p.mn = 20.0 * AMU_TO_KG;
        p.kB = BOLTZMANN_CONSTANT;
        p.Tn = 300.0;
        p.ubx = 0.0;
        p.uby = 0.0;
        p.ubz = 0.0;
        p.sigma_eff = 25e-20;
        p.Rn = std::sqrt(p.sigma_eff / M_PI);
        
        const double E_thermal = 1.5 * p.kB * p.Tn;
        const int N_IONS = 500;
        const int N_COLLISIONS = 500;  // Faster thermalization with equal mass
        
        double sum_E_final = 0.0;
        
        #pragma omp parallel for reduction(+:sum_E_final)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            Vec3 v_ion{2000.0, 0.0, 0.0};
            EhssRng rng(3000 + ion_idx);
            
            for (int i = 0; i < N_COLLISIONS; ++i) {
                Vec3 v_neutral = sample_neutral_velocity(p, rng);
                v_ion = collide_hs_cpu(v_ion, v_neutral, p, rng);
            }
            
            double E_final = kinetic_energy_J(v_ion, p.mi);
            sum_E_final += E_final;
        }
        
        double avg_E_final = sum_E_final / N_IONS;
        double ratio = avg_E_final / E_thermal;
        
        INFO("Equal mass - Thermalization ratio: " << ratio);
        REQUIRE_THAT(ratio, WithinRel(1.0, 0.10));
    }
    
    SECTION("Temperature independence: 150K vs 600K") {
        // Verify thermalization works correctly across temperature range
        const int N_IONS = 300;
        const int N_COLLISIONS = 1000;
        
        for (double T_gas : {150.0, 300.0, 450.0, 600.0}) {
            EHSSParams p;
            p.mi = 19.0 * AMU_TO_KG;
            p.mn = 4.0 * AMU_TO_KG;
            p.kB = BOLTZMANN_CONSTANT;
            p.Tn = T_gas;
            p.ubx = 0.0;
            p.uby = 0.0;
            p.ubz = 0.0;
            p.sigma_eff = 25e-20;
            p.Rn = std::sqrt(p.sigma_eff / M_PI);
            
            const double E_thermal = 1.5 * p.kB * T_gas;
            double sum_E = 0.0;
            
            #pragma omp parallel for reduction(+:sum_E)
            for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
                Vec3 v_ion{3000.0, 0.0, 0.0};
                EhssRng rng(4000 + ion_idx * 10 + (int)T_gas);
                
                for (int i = 0; i < N_COLLISIONS; ++i) {
                    Vec3 v_neutral = sample_neutral_velocity(p, rng);
                    v_ion = collide_hs_cpu(v_ion, v_neutral, p, rng);
                }
                
                sum_E += kinetic_energy_J(v_ion, p.mi);
            }
            
            double ratio = (sum_E / N_IONS) / E_thermal;
            
            INFO("T = " << T_gas << " K: Ratio = " << ratio);
            REQUIRE_THAT(ratio, WithinRel(1.0, 0.10));
        }
    }
}


