// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *   under the influence of gas.
 *
 *   @file        integrator.h
 *   @brief       Integrates the ion trajectories with fixed RK4 or adaptive RK45
 *                implementation.
 *
 *   @details
 *   Advances the ion ensemble through time by numerically integrating Newton’s
 *   equations of motion using either a fixed-step RK4 or adaptive RK45 solver.
 *   The solver evaluates forces from electric fields, RF excitation, damping,
 *   and optional gas flow via the derivative function.

 *   Each ion’s trajectory is computed independently, with support for parallel
 *   execution if enabled. The solver handles boundary events, collision modeling,
 *   and optional reaction tracking at each step.

 *   Key features:
 *   Supports multiple instrument geometries (LQIT, SIFDT-MS, IMS, Orbitrap).
 *   Adapts to selected collision model (hard-sphere, Langevin, friction, EHSS, HSMC).
 *   Can record arrival times, reaction products, and full trajectory data.
 *   Compatible with time evaluation grid from `gParams.t_eval`.

 *   @note
 *   RK45 adapts step size based on local error estimates; RK4 uses fixed `dt_s`.
 *   Results can be exported to HDF5 for post-processing and visualization.
 *   Solver selection is based on instrument type and user configuration.
 *
 *
 *   @date        <2025-10-06>
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 *   @see         paramUtils/paramUtils.h
 *   @see         physics/collisionHelpers.h
 *   @see         types/SimulationResult.h
 *   @see         physics/ionMotion.h
 *   @see         physics/eventFunctions.h
 *   @see         physics/geometryReader.h
 *   @see         reactionUtils/reactionUtils.h
 *   @see         constants.h
 *   @see         io/hdf5Writer.h
 *
 * =====================================================================
 */

/**
 * =====================================================================
 *   under the influence of gas.
 *
 *   @file        integrator.h
 *   @brief       Integrates the ion trajectories with fixed RK4 or adaptive RK45
 *                implementation.
 *
 *   @details
 *   Advances the ion ensemble through time by numerically integrating Newton’s
 *   equations of motion using either a fixed-step RK4 or adaptive RK45 solver.
 *   The solver evaluates forces from electric fields, RF excitation, damping,
 *   and optional gas flow via the derivative function.
 *
 *   Each ion’s trajectory is computed independently, with support for parallel
 *   execution if enabled. The solver handles boundary events, collision modeling,
 *   and optional reaction tracking at each step.
 *
 *   Key features:
 *   Supports multiple instrument geometries (LQIT, SIFDT-MS, IMS, Orbitrap).
 *   Adapts to selected collision model (hard-sphere, Langevin, friction, EHSS, HSMC).
 *   Can record arrival times, reaction products, and full trajectory data.
 *   Compatible with time evaluation grid from `gParams.t_eval`.
 *
 *   @note
 *   RK45 adapts step size based on local error estimates; RK4 uses fixed `dt_s`.
 *   Results can be exported to HDF5 for post-processing and visualization.
 *   Solver selection is based on instrument type and user configuration.
 *
 *
 *   @date        <2025-10-06>
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 *   @see         paramUtils/paramUtils.h
 *   @see         physics/collisionHelpers.h
 *   @see         types/SimulationResult.h
 *   @see         physics/ionMotion.h
 *   @see         physics/eventFunctions.h
 *   @see         physics/geometryReader.h
 *   @see         reactionUtils/reactionUtils.h
 *   @see         constants.h
 *   @see         io/hdf5Writer.h
 *
 * =====================================================================
 */

#pragma once

#include <string>
#include <utility>
#include <vector>

#include "core/param/paramUtils.h"

#include "core/physics/collisions/collisionHelpers.h"
#include "core/physics/geometryReader.h"
#include "core/physics/computeAccelerations.h"

#include "core/io/hdf5Writer.h"

#include "H5Cpp.h"
#include "utils/constants.h"
#include "core/physics/reactions/reactionUtils.h"
#include "core/types/SimulationResult.h"
#include "core/io/logger.h"
// Context-based SimulationContext
#include "optimizer/SimulationContext.h"

#include "integrator/rk45_settings.h"

namespace ICARION {
namespace trajectory {

// -----------------------------
// Solver function
// -----------------------------

/**
 * @brief Integrate ion trajectories from t_start to t_end using RK4 or adaptive RK45
 * 
 * @param ions Vector of ion states (modified in-place during integration)
 * @param t_start Starting simulation time [s]
 * @param t_end Ending simulation time [s]
 * @param dt Fixed time step [s] (used for RK4; RK45 adapts internally)
 * @param gParams Global simulation parameters (fields, geometry, collision model, etc.)
 * @param speciesDB Map of species name to Species properties (mass, charge, collision cross-section)
 * @param reaction_list Vector of ion-molecule reactions with rate constants
 * @param instrumentDomains Spatial domains with different field/collision configurations
 * @param rk45 Adaptive RK45 settings (tolerances, min/max step sizes)
 * @param ctx Simulation context for optimizer integration (optional)
 * @param logger Run logger for output (optional)
 * 
 * @return SimulationResult containing arrival times, collision statistics, trajectory data
 * 
 * Integrates equations of motion for ion ensemble under electric fields, collisions,
 * and reactions. Supports multiple instrument types (LQIT, SIFDT-MS, IMS, Orbitrap)
 * and collision models (hard-sphere, Langevin, friction, EHSS, HSMC).
 * 
 * Writes trajectory snapshots to HDF5 at intervals defined by gParams.t_eval.
 * Handles boundary conditions (absorbing, reflecting, periodic) and detection events.
 */
SimulationResult integrate_trajectory(std::vector<IonState>& ions, double t_start, double t_end, double dt,
						   GlobalParams& gParams, 
						   const ICARION::io::SpeciesDatabase& speciesDB,
						   const std::vector<ReactionEntry>&               reaction_list,
						   const std::vector<InstrumentDomain>& instrumentDomains,
						   const RK45Settings& rk45 = RK45Settings(),
						   ICARION::io::RunLogger* logger = nullptr);

// Context-aware overload (modern API). Accepts a single SimulationContext
// containing all simulation state and settings and forwards to the
// canonical integrator implementation.
SimulationResult integrate_trajectory(optimization::SimulationContext& ctx,
                                     const RK45Settings& rk45 = RK45Settings(),
                                     ICARION::io::RunLogger* logger = nullptr);

// Legacy implementation is kept under a different name to avoid symbol
// collisions; the legacy definition lives in `src/solvers/integrator.cpp` and
// is available as `integrate_trajectory_legacy(...)` if needed for testing or
// reference.
SimulationResult integrate_trajectory_legacy(std::vector<IonState>& ions, double t_start, double t_end, double dt,
						   GlobalParams& gParams, 
						   const std::unordered_map<std::string, Species>& speciesDB,
						   const std::vector<ReactionEntry>&               reaction_list,
						   const std::vector<InstrumentDomain>& instrumentDomains,
						   const RK45Settings& rk45,
						   ICARION::io::RunLogger* logger = nullptr);

// -----------------------------
// RK helper functions
// -----------------------------

double rk_error_norm(const IonState& y5, const IonState& y4, const IonState& y, const RK45Settings& rk45);

void rk45_dp_step(double t, const IonState& y, double h, const GlobalParams& gParams, const InstrumentDomain& dom, IonState& y5,
				  IonState& y4, const std::vector<IonState>& current_ions);

void load_geometry(const std::string& filename, const std::string& targetName,
				   std::vector<Vec3>& h_centers, std::vector<double>& h_radii);

// Backwards-compatible overload that also returns an optional MobCal CCS value
// via the out parameter `ccs_mobcal` (0.0 when not present). Some code paths
// (notably the modern integrator in `src/integrator`) use this variant.
void load_geometry(const std::string& filename, const std::string& targetName,
				   std::vector<Vec3>& h_centers, std::vector<double>& h_radii,
				   double& ccs_mobcal);

// -----------------------------
// collision and reaction event helper functions
// -----------------------------                
void handle_collision(IonState& y, EhssRng& rng, double dt, const GlobalParams& gParams,
					  const std::unordered_map<std::string, std::pair<std::vector<Vec3>, std::vector<double>>>& geometry_map);

void handle_reaction(IonState& y, EhssRng& rng, double dt, const GlobalParams& gParams,
					 const ICARION::io::SpeciesDatabase& speciesDB,
					 const std::vector<ReactionEntry>& reaction_list);

// -----------------------------
// integration and loggin function
// -----------------------------        
void log_step(double t, double t_end, int& step_count, double& next_write_time, double write_dt,
					 const std::vector<IonState>& ions, std::vector<double>& times_buffer,
					 std::vector<std::vector<IonState>>& trajectory_buffer, size_t buffer_max,
					 const std::string& filename);

void integrate_one_step(
	IonState& ion,
	const GlobalParams& gParams,
	const InstrumentDomain& dom,
	const ICARION::io::SpeciesDatabase& speciesDB,
	const std::vector<ReactionEntry>& reaction_list,
	EhssRng& local_rng,
	const RK45Settings& rk45,
	const std::vector<IonState>& current_ions
);

// -----------------------------
// domain transition helper functions
// -----------------------------
bool check_aperture_crossing(IonState ion, const InstrumentDomain& dom,
								   const std::vector<InstrumentDomain>& domains,
								   const Vec3& local_old, const Vec3& local_new);

}  // namespace trajectory
}  // namespace ICARION

// Bring into global namespace for backward compatibility  
using ICARION::trajectory::integrate_trajectory;
using ICARION::trajectory::integrate_trajectory_legacy;
using ICARION::trajectory::rk_error_norm;
using ICARION::trajectory::rk45_dp_step;
using ICARION::trajectory::load_geometry;
using ICARION::trajectory::handle_collision;
using ICARION::trajectory::handle_reaction;
using ICARION::trajectory::log_step;
using ICARION::trajectory::integrate_one_step;
using ICARION::trajectory::check_aperture_crossing;