// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// Debug test for SimulationEngine + Orbitrap bug

#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <iomanip>

#include "helpers/physics_sim_utils.h"
#include "core/integrator/SimulationEngine.h"
#include "core/config/types/FullConfig.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::tests;

TEST_CASE("SimulationEngine: Orbitrap bug reproducer", "[debug][orbitrap][engine]") {
    // Minimal Orbitrap config (same as working manual loop test)
    config::FullConfig cfg;
    cfg.simulation.dt_s = 1e-9;
    cfg.simulation.total_time_s = 1e-5;  // 10 µs
    cfg.simulation.write_interval = 100;
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    config::DomainConfig dom;
    dom.instrument = config::Instrument::Orbitrap;
    dom.name = "orbitrap";
    dom.domain_index = 0;
    
    double R_in = 6e-3;
    double R_out = 15e-3;
    double R_char = 22e-3;
    
    dom.geometry.radius_in_m = R_in;
    dom.geometry.radius_out_m = R_out;
    dom.geometry.radius_char_m = R_char;
    dom.geometry.length_m = 0.04;
    dom.geometry.origin_m = {0.0, 0.0, 0.0};
    
    dom.fields.dc.radial_V.constant_value = 3500.0;
    dom.fields.dc.axial_V.constant_value = 0.0;
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 0.0;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 0.0;
    dom.fields.ac.compute_derived();
    
    dom.environment.pressure_Pa = 1e-7;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    cfg.domains.push_back(dom);
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "debug_orbitrap.h5";
    cfg.output.print_progress = false;
    
    // Create ion
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {9e-3, 0.0, 6e-3};  // (r=9mm, z=6mm)
    
    double mass_amu = 200.0;
    double E_tang_eV = 1600.0;
    double m_kg = mass_amu * AMU_TO_KG;
    double v_tang = std::sqrt(2.0 * E_tang_eV * ELEM_CHARGE_C / m_kg);
    
    ion.vel = {0.0, v_tang, 0.0};
    ion.mass_kg = m_kg;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-19;
    ion.active = true;
    ion.born = true;
    ion.current_domain_index = 0;
    
    std::cout << "\n=== REPRODUCER: SimulationEngine + Orbitrap ===\n";
    std::cout << "Initial: r=" << 9.0 << "mm, z=" << 6.0 << "mm\n";
    std::cout << "Initial v_tang: " << v_tang << " m/s\n";
    std::cout << "Simulation: " << cfg.simulation.total_time_s*1e6 << " µs\n";
    
    SECTION("Test with run_simple_simulation(trace=false) - uses SimulationEngine") {
        std::cout << "\n--- Using SimulationEngine (trace=false) ---\n";
        auto result = run_simple_simulation(cfg, {ion}, false);
        
        REQUIRE(result.ions.size() == 1);
        const auto& final = result.ions[0];
        
        double r_final = std::sqrt(final.pos.x*final.pos.x + final.pos.y*final.pos.y);
        std::cout << "Final: r=" << r_final*1000 << "mm, z=" << final.pos.z*1000 << "mm\n";
        std::cout << "Active: " << final.active << "\n";
        std::cout << "Time: " << final.t*1e6 << " µs\n";
        
        // BUG: Ion should still be active but engine stops it immediately
        WARN("SimulationEngine BUG: Ion active = " << final.active << " (should be true)");
        CHECK(final.active);  // This will FAIL due to bug
    }
    
    SECTION("Test with run_simple_simulation(trace=true) - uses manual loop") {
        std::cout << "\n--- Using manual loop (trace=true) ---\n";
        auto result = run_simple_simulation(cfg, {ion}, true);
        
        REQUIRE(result.ions.size() == 1);
        const auto& final = result.ions[0];
        
        double r_final = std::sqrt(final.pos.x*final.pos.x + final.pos.y*final.pos.y);
        std::cout << "Final: r=" << r_final*1000 << "mm, z=" << final.pos.z*1000 << "mm\n";
        std::cout << "Active: " << final.active << "\n";
        std::cout << "Time: " << final.t*1e6 << " µs\n";
        
        // This works correctly
        REQUIRE(final.active);
    }
}
