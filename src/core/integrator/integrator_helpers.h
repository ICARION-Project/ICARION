// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file integrator_helpers.h
 * @brief Helper functions for trajectory integration
 * 
 * Internal utilities used by integrate_trajectory() for:
 * - Single timestep integration (RK4, RK45)
 * - Collision and reaction handling
 * - Trajectory buffering and HDF5 output
 * - Molecular geometry loading for EHSS collisions
 */
#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include "integrator.h"

namespace integrator_helpers {

/**
 * @brief Log trajectory snapshot and flush buffer if full
 * @param t Current simulation time [s]
 * @param t_end Final simulation time [s]
 * @param step_count Step counter (incremented)
 * @param next_write_time Next scheduled write time [s]
 * @param write_dt Time interval between trajectory snapshots [s]
 * @param ions Current ion states
 * @param times_buffer Buffered time points
 * @param trajectory_buffer Buffered ion states
 * @param buffer_max Maximum buffer size before flushing to HDF5
 * @param filename HDF5 output file path
 * 
 * Appends current state to buffer and writes to HDF5 when buffer is full.
 * Used to reduce I/O overhead during long simulations.
 */
void log_step(double t, double t_end, int& step_count, double& next_write_time, double write_dt,
             const std::vector<IonState>& ions, std::vector<double>& times_buffer,
             std::vector<std::vector<IonState>>& trajectory_buffer, size_t buffer_max,
             const std::string& filename);

/**
 * @brief Integrate single ion for one timestep
 * @param ion Ion state (modified in-place)
 * @param gParams Global simulation parameters
 * @param dom Current instrument domain
 * @param speciesDB Species database
 * @param reaction_list Available reactions
 * @param local_rng Random number generator for collisions
 * @param rk45 RK45 adaptive settings
 * @param current_ions All ions (for space-charge if enabled)
 * @param field_provider Field interpolation provider
 * 
 * Advances ion by one timestep using selected solver (RK4 or RK45).
 * Handles collisions (stochastic or deterministic) and reactions.
 * Updates ion position, velocity, and domain index.
 */
void integrate_one_step(
    IonState& ion,
    const GlobalParams& gParams,
    const InstrumentDomain& dom,
    const std::unordered_map<std::string, Species>& speciesDB,
    const std::vector<ReactionEntry>& reaction_list,
    EhssRng& local_rng,
    const RK45Settings& rk45,
    const std::vector<IonState>& current_ions
    , const IFieldProvider* field_provider = nullptr
);

// NOTE: load_geometry(), compute_CCS_from_geometry(), and handle_collision() removed
// in Phase 2D refactor. See src/core/physics/collisions/ for modern collision system:
// - load_geometry_from_file() in utils.h
// - CollisionHandlerFactory::create() for handler instantiation
// - ICollisionHandler::handle_collision() for collision processing

/**
 * @brief Handle ion-molecule reaction event
 * @param y Ion state (species and properties modified if reaction occurs)
 * @param rng Random number generator
 * @param dt Timestep [s]
 * @param gParams Global parameters
 * @param speciesDB Species database
 * @param reaction_list Available reactions
 * 
 * Checks for possible reactions based on current ion species.
 * Computes reaction probabilities from rate constants and neutral concentrations.
 * If reaction occurs, updates ion to product species with new mass/charge/mobility.
 */
void handle_reaction(IonState& y, EhssRng& rng, double dt, const GlobalParams& gParams,
                     const std::unordered_map<std::string, Species>& speciesDB,
                     const std::vector<ReactionEntry>& reaction_list);

/**
 * @brief Compute RK45 error estimate for adaptive timestep control
 * @param y5 5th-order RK solution
 * @param y4 4th-order RK solution (embedded)
 * @param y Current ion state
 * @param rk45 RK45 settings (tolerances)
 * @return Error norm (scaled by position and velocity tolerances)
 * 
 * Uses difference between 4th and 5th order solutions to estimate local error.
 * Error norm combines position and velocity errors with user-specified tolerances.
 */
double rk_error_norm(const IonState& y5, const IonState& y4, const IonState& y, const RK45Settings& rk45);

/**
 * @brief Single RK45 Dormand-Prince step
 * @param t Current time [s]
 * @param y Current ion state
 * @param h Timestep [s]
 * @param gParams Global parameters
 * @param dom Current domain
 * @param y5 Output: 5th-order solution
 * @param y4 Output: 4th-order solution (for error estimation)
 * @param current_ions All ions (for space-charge)
 * @param field_provider Field interpolation provider
 * 
 * Advances ion using 6-stage Dormand-Prince RK45 method.
 * Provides both 4th and 5th order solutions for adaptive timestep control.
 * Evaluates forces/accelerations at 6 intermediate points per step.
 */
void rk45_dp_step(double t, const IonState& y, double h, const GlobalParams& gParams,
                  const InstrumentDomain& dom, IonState& y5, IonState& y4, const std::vector<IonState>& current_ions,
                  const IFieldProvider* field_provider = nullptr);

} // namespace integrator_helpers
