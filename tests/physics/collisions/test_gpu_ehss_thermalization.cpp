// Test GPU EHSS thermalization with H3O+ geometry
#include "core/gpu/GPUContext.h"
#include "core/gpu/GPUCollisionHelper.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>

using namespace ICARION;
using namespace ICARION::config;

// Simple JSON parser for H3O+ molecule
struct Atom {
    Vec3 pos;
    double radius_angstrom;
};

std::vector<Atom> parse_h3o_geometry(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open " + filepath);
    }
    
    std::vector<Atom> atoms;
    std::string line;
    bool in_atoms = false;
    
    while (std::getline(file, line)) {
        // Find "pos": [x, y, z]
        if (line.find("\"pos\":") != std::string::npos) {
            size_t start = line.find('[');
            size_t end = line.find(']');
            if (start != std::string::npos && end != std::string::npos) {
                std::string coords = line.substr(start + 1, end - start - 1);
                std::istringstream iss(coords);
                double x, y, z;
                char comma;
                iss >> x >> comma >> y >> comma >> z;
                
                Atom atom;
                atom.pos = Vec3{x, y, z};  // Already in Angstrom from JSON
                atoms.push_back(atom);
            }
        }
        // Find "LJ_sigma_angstrom": value (used as atomic radius)
        if (line.find("\"LJ_sigma_angstrom\":") != std::string::npos && !atoms.empty()) {
            size_t start = line.find(':');
            if (start != std::string::npos) {
                double sigma;
                std::istringstream iss(line.substr(start + 1));
                iss >> sigma;
                atoms.back().radius_angstrom = sigma / 2.0;  // LJ sigma to radius
            }
        }
    }
    
    if (atoms.empty()) {
        throw std::runtime_error("No atoms found in " + filepath);
    }
    
    return atoms;
}

int main() {
    const int N_IONS = 5000;
    const double mass_kg = 19.0 * 1.66054e-27;  // H3O+ (1+1+1+16 = 19 AMU)
    const double T = 300.0;
    const double kB = 1.380649e-23;
    const double eV = 1.60218e-19;
    
    // Setup environment
    EnvironmentConfig env;
    env.temperature_K = T;
    env.pressure_Pa = 101325.0;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    
    GasMixtureComponent he;
    he.species = "He";
    he.mole_fraction = 1.0;
    he.mass_kg = 4.0026 * 1.66054e-27;
    he.radius_m = 1.4e-10;
    env.gas_mixture.push_back(he);
    
    // Initial velocity (10x thermal)
    double v_thermal = std::sqrt(3 * kB * T / mass_kg);
    double v_init = v_thermal * std::sqrt(10.0);
    
    std::vector<IonState> ions(N_IONS);
    for (auto& ion : ions) {
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.CCS_m2 = 24.9e-20;  // From H3O+.json
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v_init, 0.0, 0.0};
        ion.active = true;
    }
    
    // GPU context and helper (EHSS mode)
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    auto gpu_helper = icarion::gpu::GPUCollisionHelper::create(*gpu_ctx, 1000, "EHSS", 42);
    
    // Load H3O+ geometry from data/molecules/H3O+.json
    std::cout << "Loading H3O+ geometry...\n";
    auto atoms = parse_h3o_geometry("/home/chsch95/ICARION/data/molecules/H3O+.json");
    std::cout << "Loaded " << atoms.size() << " atoms\n";
    
    // Convert to geometry map format
    icarion::gpu::GPUCollisionHelper::GeometryMap geom_map;
    std::vector<Vec3> positions;
    std::vector<double> radii;
    for (const auto& atom : atoms) {
        positions.push_back(atom.pos * 1e-10);  // Angstrom -> meters
        radii.push_back(atom.radius_angstrom * 1e-10);  // Angstrom -> meters
    }
    geom_map["H3O+"] = {positions, radii};
    
    // Upload geometry to GPU
    gpu_helper->set_geometry(geom_map);
    std::cout << "Geometry uploaded to GPU\n";
    
    double dt = 1e-7;
    std::vector<int> steps = {100, 500, 1000, 2000, 5000};
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nGPU EHSS Thermalization (He gas, H3O+ geometry)\n";
    std::cout << "Expected thermal: " << 1.5 * kB * T / eV << " eV\n\n";
    std::cout << "Steps    KE(eV)   Ratio\n";
    std::cout << "-----  --------  ------\n";
    
    for (int N : steps) {
        // Reset ions
        for (auto& ion : ions) {
            ion.vel = Vec3{v_init, 0.0, 0.0};
        }
        
        // Process collisions
        for (int i = 0; i < N; ++i) {
            gpu_helper->process_collisions_batch(ions, dt, env);
        }
        
        // Compute mean KE
        double sum_v2 = 0.0;
        for (const auto& ion : ions) {
            sum_v2 += ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        }
        double KE = 0.5 * mass_kg * (sum_v2 / N_IONS) / eV;
        double ratio = KE / (1.5 * kB * T / eV);
        
        std::cout << std::setw(5) << N 
                  << std::setw(10) << KE 
                  << std::setw(8) << ratio << "\n";
    }
    
    std::cout << "\nExpected: Ratio ≈ 1.00 (100% thermalization)\n";
    
    return 0;
}
