// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_SIMULATION_CONFIG_H
#define ICARION_CONFIG_SIMULATION_CONFIG_H

#include <string>
#include <vector>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Simulation control parameters
 * 
 * Time parameters and execution mode settings.
 * Replaces GlobalParams (simulation-related parts only).
 */
struct SimulationConfig {
    // === Time parameters ===
    double total_time_s = 1e-3;         ///< Total simulation time [s]
    double dt_s = 1e-9;                 ///< Timestep [s]
    int write_interval = 100;           ///< Steps between trajectory outputs
    
    // Derived (computed after load)
    int total_steps = 0;                ///< Computed: total_time_s / dt_s
    std::vector<double> t_eval;         ///< Time array for output snapshots
    
    // === Integrator ===
    std::string integrator = "RK4";     ///< Default integrator (can be overridden per-domain)
    
    // === Execution mode ===
    bool enable_gpu = false;            ///< GPU acceleration requested
    bool enable_openmp = false;         ///< OpenMP threading
    unsigned int rng_seed = 42;         ///< Random number generator seed
    
    // === Checkpointing ===
    std::string continue_from = "";     ///< HDF5 file to resume from
    double continue_time_s = 0.0;       ///< Additional time to simulate
    bool auto_continue_if_active = false; ///< Auto-detect incomplete runs
    
    /**
     * @brief Compute derived quantities (total_steps, t_eval)
     * 
     * Should be called after loading from JSON.
     */
    void compute_derived() {
        if (dt_s <= 0.0) {
            throw std::runtime_error("dt_s must be positive");
        }
        if (total_time_s <= 0.0) {
            throw std::runtime_error("total_time_s must be positive");
        }
        
        total_steps = static_cast<int>(total_time_s / dt_s);
        
        // Build output time array
        t_eval.clear();
        for (int i = 0; i <= total_steps; i += write_interval) {
            t_eval.push_back(i * dt_s);
        }
        // Ensure final time is included
        if (t_eval.empty() || t_eval.back() < total_time_s) {
            t_eval.push_back(total_time_s);
        }
    }
    
    /**
     * @brief Validate configuration parameters
     * 
     * @throws std::runtime_error if invalid
     */
    void validate() const {
        if (dt_s <= 0.0) {
            throw std::runtime_error("dt_s must be positive");
        }
        if (total_time_s <= 0.0) {
            throw std::runtime_error("total_time_s must be positive");
        }
        if (write_interval <= 0) {
            throw std::runtime_error("write_interval must be positive");
        }
        if (dt_s > total_time_s) {
            throw std::runtime_error("dt_s exceeds total_time_s");
        }
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_SIMULATION_CONFIG_H
