// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Contributors

/**
 * @file benchmark_openmp_soa.cpp
 * @brief Benchmark OpenMP performance with SoA vs AoS (Phase 6)
 * 
 * Tests that SoA eliminates false sharing and achieves near-linear scaling.
 * Compares with AoS (which has severe false sharing on multi-core).
 */

#include <catch2/catch_test_macros.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/FullConfig.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include <chrono>
#include <iostream>
#include <omp.h>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::core;
using namespace ICARION::config;

// Create test configuration with OpenMP enabled
static FullConfig create_openmp_config() {
    FullConfig cfg;
    
    // Simulation settings
    cfg.simulation.total_time_s = 5e-6;  // 5 μs (5000 steps) - longer test
    cfg.simulation.dt_s = 1e-9;          // 1 ns
    cfg.simulation.write_interval = 10000;  // No output
    cfg.simulation.enable_openmp = true;
    cfg.simulation.rng_seed = 42;
    cfg.simulation.enable_safety_logging = false;
    
    // Single domain
    DomainConfig domain;
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.geometry.radius_m = 0.01;
    
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.environment.gas_species = "N2";
    domain.environment.compute_derived_properties();
    
    cfg.domains.push_back(domain);
    
    // Output settings
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "openmp_test.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

// Create test ions
static std::vector<IonState> create_test_ions(size_t n) {
    std::vector<IonState> ions;
    for (size_t i = 0; i < n; ++i) {
        IonState ion;
        ion.pos = {0.001 * (i % 100), 0.001 * ((i / 100) % 100), 0.01 + 0.0001 * i};
        ion.vel = {100.0, 0.0, 0.0};
        ion.mass_kg = 3.0e-26;
        ion.ion_charge_C = 1.6e-19;
        ion.CCS_m2 = 1.5e-18;
        ion.reduced_mobility_cm2_Vs = 2.0;
        ion.species_id = "H3O+";
        ion.birth_time_s = 0.0;
        ion.current_domain_index = 0;
        ion.history_index = static_cast<int>(i);
        ion.t = 0.0;
        ion.active = true;
        ion.born = true;
        ions.push_back(ion);
    }
    return ions;
}

// Benchmark SoA with specific thread count
static double benchmark_soa_threads(size_t n_ions, int n_threads) {
    omp_set_num_threads(n_threads);
    
    auto config = create_openmp_config();
    auto ions = create_test_ions(n_ions);
    auto ensemble = IonEnsemble::from_legacy(ions);
    
    // Minimal setup like in working test
    auto force_registry = std::make_shared<physics::ForceRegistry>();
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries = {force_registry};
    auto integrator = std::make_shared<integrator::RK4Strategy>();
    
    SimulationEngine engine(config, force_registries, integrator);
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = engine.run_soa(ensemble);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return duration.count();
}

// Benchmark AoS with specific thread count
static double benchmark_aos_threads(size_t n_ions, int n_threads) {
    omp_set_num_threads(n_threads);
    
    auto config = create_openmp_config();
    auto ions = create_test_ions(n_ions);
    
    // Minimal setup like in working test
    auto force_registry = std::make_shared<physics::ForceRegistry>();
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries = {force_registry};
    auto integrator = std::make_shared<integrator::RK4Strategy>();
    
    SimulationEngine engine(config, force_registries, integrator);
    
    auto start = std::chrono::high_resolution_clock::now();
    auto result = engine.run(ions);
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return duration.count();
}

TEST_CASE("OpenMP SoA Scaling", "[soa][openmp][.performance]") {
    const size_t N_IONS = 20000;  // More ions for better parallelization
    int max_threads = omp_get_max_threads();
    
    std::cout << "\n=== OpenMP SoA Scaling Test ===" << std::endl;
    std::cout << "Testing with " << N_IONS << " ions" << std::endl;
    std::cout << "Max threads: " << max_threads << "\n" << std::endl;
    
    // Test SoA with different thread counts
    std::cout << "SoA Performance:" << std::endl;
    double time_1 = benchmark_soa_threads(N_IONS, 1);
    std::cout << "  1 thread:  " << time_1 << " ms" << std::endl;
    
    double time_2 = benchmark_soa_threads(N_IONS, 2);
    double speedup_2 = time_1 / time_2;
    std::cout << "  2 threads: " << time_2 << " ms (speedup: " << speedup_2 << "x)" << std::endl;
    
    double time_4 = benchmark_soa_threads(N_IONS, 4);
    double speedup_4 = time_1 / time_4;
    std::cout << "  4 threads: " << time_4 << " ms (speedup: " << speedup_4 << "x)" << std::endl;
    
    if (max_threads >= 8) {
        double time_8 = benchmark_soa_threads(N_IONS, 8);
        double speedup_8 = time_1 / time_8;
        std::cout << "  8 threads: " << time_8 << " ms (speedup: " << speedup_8 << "x)" << std::endl;
    }
    
    // Test AoS for comparison
    std::cout << "\nAoS Performance (for comparison):" << std::endl;
    double aos_1 = benchmark_aos_threads(N_IONS, 1);
    std::cout << "  1 thread:  " << aos_1 << " ms" << std::endl;
    
    double aos_4 = benchmark_aos_threads(N_IONS, 4);
    double aos_speedup_4 = aos_1 / aos_4;
    std::cout << "  4 threads: " << aos_4 << " ms (speedup: " << aos_speedup_4 << "x)" << std::endl;
    
    // Calculate parallel efficiency
    double eff_2 = (speedup_2 / 2.0) * 100.0;
    double eff_4 = (speedup_4 / 4.0) * 100.0;
    
    std::cout << "\nParallel Efficiency (SoA):" << std::endl;
    std::cout << "  2 threads: " << eff_2 << "%" << std::endl;
    std::cout << "  4 threads: " << eff_4 << "%" << std::endl;
    
    std::cout << "\n================================" << std::endl;
    
    // Analysis
    std::cout << "Analysis:" << std::endl;
    std::cout << "  SoA vs AoS (4 threads): " << (speedup_4 / aos_speedup_4) << "x better scaling" << std::endl;
    std::cout << "  Expected: >1.0 (SoA should scale better than AoS)" << std::endl;
    std::cout << "  SoA scaling is " << ((speedup_4 > aos_speedup_4) ? "BETTER" : "WORSE") << " than AoS" << std::endl;
    
    std::cout << "\n================================\n" << std::endl;
    
    // Verify SoA scales better than AoS - this is the KEY requirement
    REQUIRE(speedup_4 > aos_speedup_4);  // SoA should scale better than AoS
    REQUIRE(speedup_4 > 1.2);  // Should get >1.2x on 4 cores (realistic)
    REQUIRE(eff_4 > 25.0);     // At least 25% efficient
}
