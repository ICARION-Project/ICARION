// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_CONFIG_SIMULATION_CONFIG_H
#define ICARION_CONFIG_SIMULATION_CONFIG_H

#include "../validation/ValidationResult.h"
#include <string>
#include <vector>
#include <stdexcept>
#include <optional>

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
    double rk45_min_step_s = 0.0;       ///< Absolute minimum dt for RK45 (0 = disabled, use relative factor)
    struct RK45RuntimeSettings {
        double atol = 0.0;
        double rtol = 0.0;
        double safety_factor = 0.0;
        double min_step_factor = 0.0;
        double max_step_factor = 0.0;
        double max_step_increase = 0.0;
        double max_step_decrease = 0.0;
        double absolute_min_step_s = 0.0;
    };
    std::optional<RK45RuntimeSettings> rk45_runtime_settings; ///< Populated at runtime from RK45Strategy (not user-configurable)
    
    // === Execution mode ===
    bool enable_gpu = false;            ///< GPU acceleration requested
    bool enable_openmp = false;         ///< OpenMP threading
    unsigned int rng_seed = 42;         ///< Random number generator seed
    
    // === Numerical Safety ===
    bool enable_safety_logging = false; ///< Enable detailed safety violation logging
    bool verbose_safety = false;        ///< Verbose safety output (performance impact)
    
    struct SafetyChecks {
        bool enable_nan_checks = true;       ///< Check for NaN/Inf values
        bool enable_bounds_checks = false;   ///< Check position/velocity bounds
        double max_position_m = 10.0;        ///< Maximum position magnitude [m]
        double max_velocity_ms = 1e6;        ///< Maximum velocity magnitude [m/s]
        double max_acceleration_ms2 = 1e12;  ///< Maximum acceleration magnitude [m/s²]
        bool throw_on_violation = false;     ///< Throw exception on violation (vs deactivate ion)
        bool attempt_recovery = false;       ///< Attempt to recover invalid values
    } safety_checks;
    
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
     * @brief Validate simulation parameters
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (dt_s <= 0.0) {
            result.add_error("dt_s must be positive");
        }
        if (total_time_s <= 0.0) {
            result.add_error("total_time_s must be positive");
        }
        if (write_interval <= 0) {
            result.add_error("write_interval must be positive");
        }
        if (dt_s > total_time_s) {
            result.add_error("dt_s exceeds total_time_s");
        }
        if (rk45_min_step_s < 0.0) {
            result.add_error("rk45_min_step_s must be non-negative");
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_SIMULATION_CONFIG_H
