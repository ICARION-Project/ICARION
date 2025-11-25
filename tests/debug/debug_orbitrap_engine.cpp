// Debug: Why does SimulationEngine stop immediately for Orbitrap?

#include <iostream>
#include <iomanip>
#include "core/integrator/SimulationEngine.h"
#include "core/config/types/FullConfig.h"
#include "core/IonState.h"
#include "utils/constants.h"

using namespace ICARION;

int main() {
    std::cout << "=== ORBITRAP SIMULATION ENGINE DEBUG ===\n\n";
    
    // Minimal Orbitrap config
    config::FullConfig cfg;
    cfg.simulation.dt_s = 1e-9;  // 1 ns
    cfg.simulation.total_time_s = 1e-5;  // 10 µs
    cfg.simulation.write_interval = 100;
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    // Orbitrap domain
    config::DomainConfig dom;
    dom.instrument = config::Instrument::Orbitrap;
    dom.name = "orbitrap";
    dom.domain_index = 0;
    
    // Geometry
    double R_in = 6e-3;
    double R_out = 15e-3;
    double R_char = 22e-3;
    
    dom.geometry.radius_in_m = R_in;
    dom.geometry.radius_out_m = R_out;
    dom.geometry.radius_char_m = R_char;
    dom.geometry.length_m = 0.04;  // ±20mm
    dom.geometry.origin_m = {0.0, 0.0, 0.0};
    
    // Fields
    dom.fields.dc.radial_V.constant_value = 3500.0;
    dom.fields.dc.axial_V.constant_value = 0.0;
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 0.0;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 0.0;
    dom.fields.ac.compute_derived();
    
    // Environment
    dom.environment.pressure_Pa = 1e-7;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    cfg.domains.push_back(dom);
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "debug_orbitrap.h5";
    cfg.output.print_progress = true;
    
    // Create ion at (r=9mm, z=6mm) with tangential velocity
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {9e-3, 0.0, 6e-3};
    
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
    
    std::cout << "Initial ion:\n";
    std::cout << "  Position: (" << ion.pos.x*1000 << ", " << ion.pos.y*1000 
              << ", " << ion.pos.z*1000 << ") mm\n";
    std::cout << "  Velocity: (" << ion.vel.x << ", " << ion.vel.y 
              << ", " << ion.vel.z << ") m/s\n";
    std::cout << "  v_tang = " << v_tang << " m/s\n";
    std::cout << "  m/z = " << mass_amu << "\n";
    std::cout << "  Active: " << ion.active << "\n\n";
    
    // Finalize config
    cfg.finalize_all();
    
    std::cout << "Running SimulationEngine...\n";
    std::cout << "Total steps: " << cfg.simulation.total_steps << "\n";
    std::cout << "dt = " << cfg.simulation.dt_s*1e9 << " ns\n\n";
    
    try {
        // Build force registries
        std::vector<std::shared_ptr<physics::ForceRegistry>> registries;
        for (const auto& d : cfg.domains) {
            auto reg = std::make_shared<physics::ForceRegistry>(d);
            registries.push_back(reg);
        }
        
        // Build integrator
        auto integrator = std::make_shared<integrator::RK4Strategy>();
        
        // Run engine
        integrator::SimulationEngine engine(cfg, registries, integrator, nullptr, nullptr);
        auto result = engine.run({ion});
        
        std::cout << "\n=== RESULT ===\n";
        std::cout << "Ions returned: " << result.size() << "\n";
        if (!result.empty()) {
            const auto& final = result[0];
            std::cout << "Final position: (" << final.pos.x*1000 << ", " 
                      << final.pos.y*1000 << ", " << final.pos.z*1000 << ") mm\n";
            std::cout << "Final velocity: (" << final.vel.x << ", " << final.vel.y 
                      << ", " << final.vel.z << ") m/s\n";
            std::cout << "Active: " << final.active << "\n";
            std::cout << "Time: " << final.t*1e6 << " µs\n";
            
            double r_final = std::sqrt(final.pos.x*final.pos.x + final.pos.y*final.pos.y);
            std::cout << "r_final = " << r_final*1000 << " mm\n";
            
            if (!final.active) {
                std::cout << "\n⚠️  ION DEACTIVATED!\n";
                std::cout << "Position barely changed → likely immediate deactivation\n";
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}
