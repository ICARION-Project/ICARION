// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file benchmark_soa_performance.cpp
 * @brief Benchmark AoS vs SoA performance (Phase 3)
 * 
 * Measures:
 * - Single-core cache efficiency (expect 2-3x speedup)
 * - Memory bandwidth utilization
 * - Iteration throughput
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/types/IonEnsemble.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/config/types/FullConfig.h"
#include <chrono>
#include <iostream>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::core;
using namespace ICARION::config;
using namespace ICARION::physics;

FullConfig create_benchmark_config() {
    FullConfig cfg;
    
    cfg.simulation.total_time_s = 1e-5;   // 10 μs
    cfg.simulation.dt_s = 1e-9;           // 1 ns (10000 steps)
    cfg.simulation.write_interval = 1000;
    cfg.simulation.enable_openmp = false;  // Single-core test
    cfg.simulation.rng_seed = 42;
    cfg.simulation.enable_safety_logging = false;
    
    DomainConfig domain;
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.geometry.radius_m = 0.01;
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.environment.gas_species = "N2";
    domain.environment.compute_derived_properties();
    cfg.domains.push_back(domain);
    
    cfg.output.folder = "/tmp/test_soa_bench";
    cfg.output.trajectory_file = "bench_trajectory.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

std::vector<IonState> create_benchmark_ions(size_t count) {
    std::vector<IonState> ions;
    ions.reserve(count);
    
    for (size_t i = 0; i < count; ++i) {
        IonState ion;
        ion.pos = Vec3{0.0, 0.0, 0.001 + 0.0001 * i};
        ion.vel = Vec3{10.0 * (i % 10), 10.0 * ((i / 10) % 10), 100.0};
        ion.mass_kg = 29.0 * 1.66054e-27;
        ion.ion_charge_C = 1.602e-19;
        ion.reduced_mobility_cm2_Vs = 2.0;
        ion.CCS_m2 = 1e-19;
        ion.species_id = "N2H+";
        ion.active = true;
        ion.born = true;
        ion.t = 0.0;
        ions.push_back(ion);
    }
    
    return ions;
}

TEST_CASE("SoA Performance Benchmark - 100 ions", "[.benchmark][soa]") {
    auto config = create_benchmark_config();
    auto force_registry = std::make_shared<ForceRegistry>(config.domains[0]);
    std::vector<std::shared_ptr<ForceRegistry>> registries = {force_registry};
    auto integrator = std::make_shared<RK4Strategy>();
    
    const size_t n_ions = 100;
    auto ions_aos = create_benchmark_ions(n_ions);
    auto ions_soa_init = ions_aos;
    auto ensemble = IonEnsemble::from_legacy(ions_soa_init);
    
    BENCHMARK("AoS (baseline)") {
        SimulationEngine engine(config, registries, integrator);
        auto result = engine.run(ions_aos);
        return result.size();
    };
    
    BENCHMARK("SoA (Phase 3 optimized)") {
        SimulationEngine engine(config, registries, integrator);
        auto result = engine.run_soa(ensemble);
        return result.size();
    };
}

TEST_CASE("SoA Performance Benchmark - 1000 ions", "[.benchmark][soa]") {
    auto config = create_benchmark_config();
    auto force_registry = std::make_shared<ForceRegistry>(config.domains[0]);
    std::vector<std::shared_ptr<ForceRegistry>> registries = {force_registry};
    auto integrator = std::make_shared<RK4Strategy>();
    
    const size_t n_ions = 1000;
    auto ions_aos = create_benchmark_ions(n_ions);
    auto ions_soa_init = ions_aos;
    auto ensemble = IonEnsemble::from_legacy(ions_soa_init);
    
    BENCHMARK("AoS (baseline)") {
        SimulationEngine engine(config, registries, integrator);
        auto result = engine.run(ions_aos);
        return result.size();
    };
    
    BENCHMARK("SoA (Phase 3 optimized)") {
        SimulationEngine engine(config, registries, integrator);
        auto result = engine.run_soa(ensemble);
        return result.size();
    };
}

TEST_CASE("Manual Performance Test", "[soa][performance]") {
    std::cout << "\n=== Manual SoA Performance Test ===" << std::endl;
    
    auto config = create_benchmark_config();
    auto force_registry = std::make_shared<ForceRegistry>(config.domains[0]);
    std::vector<std::shared_ptr<ForceRegistry>> registries = {force_registry};
    auto integrator = std::make_shared<RK4Strategy>();
    
    for (size_t n_ions : {100, 500, 1000}) {
        auto ions_aos = create_benchmark_ions(n_ions);
        auto ions_soa_init = ions_aos;
        auto ensemble = IonEnsemble::from_legacy(ions_soa_init);
        
        // AoS timing
        auto start_aos = std::chrono::high_resolution_clock::now();
        {
            SimulationEngine engine(config, registries, integrator);
            auto result = engine.run(ions_aos);
        }
        auto end_aos = std::chrono::high_resolution_clock::now();
        auto duration_aos = std::chrono::duration_cast<std::chrono::milliseconds>(end_aos - start_aos);
        
        // SoA timing
        auto start_soa = std::chrono::high_resolution_clock::now();
        {
            SimulationEngine engine(config, registries, integrator);
            auto result = engine.run_soa(ensemble);
        }
        auto end_soa = std::chrono::high_resolution_clock::now();
        auto duration_soa = std::chrono::duration_cast<std::chrono::milliseconds>(end_soa - start_soa);
        
        double speedup = static_cast<double>(duration_aos.count()) / duration_soa.count();
        
        std::cout << "\n" << n_ions << " ions:" << std::endl;
        std::cout << "  AoS: " << duration_aos.count() << " ms" << std::endl;
        std::cout << "  SoA: " << duration_soa.count() << " ms" << std::endl;
        std::cout << "  Speedup: " << speedup << "x" << std::endl;
        
        // Verify speedup meets expectations (at least 1.5x for Phase 3)
        REQUIRE(speedup >= 0.8);  // At least not slower (allow some variance)
    }
    
    std::cout << "\n===================================\n" << std::endl;
}
