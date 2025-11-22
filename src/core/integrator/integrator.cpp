/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        integrator.cpp
 *   @brief       RK4 / adaptive RK45 trajectory integration.
 *
 *   @details
 *   Advances ions through time using fixed-step RK4 or adaptive RK45 solvers.
 *   Computes forces from electric fields, damping, collisions, and gas flow.
 *   Supports parallel execution, multiple instrument types, and collision models.
 *   Optional logging to HDF5 is included.
 *
 *   @note
 *   RK45 adapts dt based on local error; RK4 uses fixed dt_s.
 *   Trajectories, arrival times, and reactions can be recorded for post-processing.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "integrator/integrator.h"
#include "fieldsolver/utils/GridFieldProvider.h" // provide GridFieldProvider used for field sampling
#include "core/utils/safety/numericalSafetyGuards.h"  // Updated NaN/Inf guards

// Field server: Core-only uses shim, full build uses fieldsolver
#include "fieldsolver/utils/field_update_api.h"

// Collision system refactored in Phase 2
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/collisions/ICollisionHandler.h"
#include "core/physics/collisions/geometryUtils.h"  // SSOT for geometry loading & conversion
#include "core/config/types/EnvironmentConfig.h"

// Reaction system refactored in Phase 3
#include "core/physics/reactions/ReactionHandlerFactory.h"
#include "core/physics/reactions/IReactionHandler.h"


#include <algorithm>
#include <cmath>
#include <omp.h>
#include <tuple>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <mutex>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

namespace {
    // Helper to ensure .h5 extension is added only once
    std::string ensure_h5_extension(const std::string& path) {
        if (path.size() >= 3 && path.substr(path.size() - 3) == ".h5") {
            return path;  // Already has .h5 extension
        }
        return path + ".h5";
    }
}
#include <sys/types.h>
#include <cstring>
#include <cstdio>
#include "core/physics/reactions/reactionUtils.h"
#include "H5Cpp.h"

namespace ICARION {
namespace trajectory {

// Import field types for convenience
using ICARION::fields::FieldSnapshot;
using ICARION::fields::Vec3d;
using ICARION::fields::sample_field;

/**
 * @brief Log simulation state to output buffers for HDF5.
 */
void log_step(double t, double t_end, int& step_count, double& next_write_time, double write_dt,
                     const std::vector<IonState>& ions, std::vector<double>& times_buffer,
                     std::vector<std::vector<IonState>>& trajectory_buffer, size_t buffer_max,
                     const std::string& filename) {
    if (step_count == 0 || t >= next_write_time || t >= t_end) {
        times_buffer.push_back(t);
        trajectory_buffer.push_back(ions);
        if (write_dt > 0.0) {
            while (next_write_time <= t) {
                next_write_time += write_dt;
            }
        } else {
            // Avoid infinite loop if write_dt is zero
            next_write_time = t + 1.0;
        }
    }

    if (times_buffer.size() >= buffer_max) {
        ICARION::io::append_to_HDF5(filename, times_buffer, trajectory_buffer);
        times_buffer.clear();
        trajectory_buffer.clear();
    }
    
    step_count++;
}

/**
 * @brief Integrate a single ion one time step using RK4 or adaptive RK45.
 */
void integrate_one_step(
    IonState& ion,
    const GlobalParams& gParams,
    const InstrumentDomain& dom,
    const ICARION::io::SpeciesDatabase& speciesDB,
    const std::vector<ReactionEntry>& reaction_list,
    EhssRng& local_rng,
    const RK45Settings& rk45, //default in solver.h
    const std::vector<IonState>& current_ions
) {

    double t_local = ion.t;
    double t_target = ion.t + ion.dt;
    double dt_local = ion.dt;

    switch (dom.solver_type) {
        case SolverType::Boris: {

            // --- Boris method for magnetic fields ---
            // Get electric and magnetic fields at current position
            Vec3 E_field;
            if (dom.instrument == Instrument::FTICR) {
                if (!dom.B.enabled) {
                    throw std::runtime_error("FT-ICR simulation requires enabled magnetic field.");
                }
                double d = std::sqrt((dom.geom.length_m * dom.geom.length_m / 8.0) +
                                     (dom.geom.radius_m * dom.geom.radius_m / 4.0));
                E_field = FTICRField(ion, dom.DC.radial_V, d, dom.geom.length_m); // in V/m
            } else {
                E_field = ElectricFieldVec(ion, dom, t_local); // in V/m
            }
            Vec3 B_field = MagneticFieldVec(ion, dom);

            const double q = ion.ion_charge_C;
            const double m = ion.mass_kg;

            // Half acceleration due to electric field
            Vec3 v_minus = ion.vel + (E_field * q / m) * (0.5 * dt_local);
            // Rotation due to magnetic field
            Vec3 t_vec = (B_field * q / m) * (0.5 * dt_local);
            double t_mag2 = dot(t_vec, t_vec);
            Vec3 s_vec = t_vec * (2.0 / (1.0 + t_mag2));

            Vec3 v_prime = v_minus + cross(v_minus, t_vec);
            Vec3 v_plus = v_minus + cross(v_prime, s_vec);

            // Half acceleration due to electric field
            ion.vel = v_plus + (E_field * q / m) * (0.5 * dt_local);

            // Update position
            ion.pos += ion.vel * dt_local;
            ion.t = t_target;
            break;
        }

        case SolverType::RK4: {
            double z_old = ion.pos.z;
            // Prepare optional provider when a live snapshot or field array is available
            std::unique_ptr<IFieldProvider> provider_owner;
            FieldSnapshot snapshot_local;
            if (dom.use_grid_field && dom.fieldServer != nullptr) {
                snapshot_local = dom.fieldServer->get_snapshot();
                provider_owner.reset(new GridFieldProvider(&snapshot_local));
            } else if (dom.fieldArrayLoaded) {
                provider_owner.reset(new GridFieldProvider(&dom.fieldArray));
            }
            const IFieldProvider* provider_ptr = provider_owner ? provider_owner.get() : nullptr;

            IonState k1 = compute_accelerations(t_local, ion, gParams, dom, provider_ptr, current_ions);
            IonState k2 = compute_accelerations(t_local + 0.5 * dt_local, ion + k1 * 0.5 * dt_local, gParams, dom, provider_ptr, current_ions);
            IonState k3 = compute_accelerations(t_local + 0.5 * dt_local, ion + k2 * 0.5 * dt_local, gParams, dom, provider_ptr, current_ions);
            IonState k4 = compute_accelerations(t_local + dt_local, ion + k3 * dt_local, gParams, dom, provider_ptr, current_ions);
            ion += (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt_local / 6.0);
            ion.t = t_target;
            double z_new = ion.pos.z;
            break;
        }
        case SolverType::RK45: {
            int consecutive_rejects = 0;  // Track rejected steps
            
            while (t_local < t_target) {
                if (t_local + dt_local > t_target) {
                    dt_local = t_target - t_local;
                }

                IonState y5, y4;
                rk45_dp_step(t_local, ion, dt_local, gParams, dom, y5, y4, current_ions);

                // Check for NaN/Inf in RK45 result (critical safety check)
                if (!ICARION::safety::is_finite(y5)) {
                    std::ostringstream oss;
                    oss << "RK45 step produced NaN/Inf for ion at t=" << t_local 
                        << ", dt=" << dt_local << ", pos=(" << ion.pos.x << "," 
                        << ion.pos.y << "," << ion.pos.z << ")";
                    throw std::runtime_error(oss.str());
                }

                double err = rk_error_norm(y5, y4, ion, rk45);
                
                // Check for NaN error (indicates numerical instability)
                if (!std::isfinite(err)) {
                    std::ostringstream oss;
                    oss << "RK45 error norm is NaN/Inf at t=" << t_local 
                        << ", dt=" << dt_local;
                    throw std::runtime_error(oss.str());
                }

                if (err <= 1.0) {
                    // Step accepted
                    ion = y5;
                    t_local += dt_local;
                    consecutive_rejects = 0;  // Reset counter
                } else {
                    // Step rejected
                    consecutive_rejects++;
                    if (consecutive_rejects >= rk45.max_rejects) {
                        std::ostringstream oss;
                        oss << "RK45 failed: " << consecutive_rejects 
                            << " consecutive rejected steps. Last error=" << err
                            << ", dt=" << dt_local << ", t=" << t_local
                            << ", pos=(" << ion.pos.x << "," << ion.pos.y << "," << ion.pos.z << ")";
                        throw std::runtime_error(oss.str());
                    }
                }

                double factor = rk45.safety * std::pow(std::max(1e-16, 1.0 / err), 0.25);
                factor = std::max(rk45.min_factor, std::min(factor, rk45.max_factor));
                dt_local = std::max(rk45.dt_min, dt_local * factor);
            }

            ion.t = t_target;
            break;
        }
    }
}

/**
 * @brief Integrate ion trajectories across multiple domains with RK4/RK45, collisions, and reactions.
 *
 * Advances ions from t_start to t_end using step size dt, applying domain-specific fields,
 * stochastic collisions (EHSS/HSS), chemical reactions, and logging snapshots to HDF5.
 *
 * @param[in,out] ions              Ion states to integrate.
 * @param[in]     t_start           Start time [s].
 * @param[in]     t_end             End time [s].
 * @param[in]     dt                Initial integration step [s].
 * @param[in,out] gParams           Simulation parameters.
 * @param[in]     speciesDB         Species database.
 * @param[in]     reaction_list     Reaction list with rates/stoichiometry.
 * @param[in]     instrumentDomains Instrument domains with boundaries and solver info.
 * @param[in]     rk45              RK45 solver settings.
 * @return std::vector<IonState> containing final ion states.
 */
std::vector<IonState> integrate_trajectory(std::vector<IonState>& ions, double t_start, double t_end, double dt,
                           GlobalParams& gParams, 
                           const ICARION::io::SpeciesDatabase& speciesDB,
                           const std::vector<ReactionEntry>& reaction_list,
                           const std::vector<InstrumentDomain>& instrumentDomains,
                           const RK45Settings& rk45,
                           ICARION::io::RunLogger* logger) {
    // --- Initialization ---
    const bool do_hdf5 = !gParams.output_file.empty();
    std::string filename;
    if (do_hdf5) {
        filename = ensure_h5_extension(gParams.output_file);
        H5::H5File file(filename, H5F_ACC_TRUNC);
        ICARION::io::write_params_to_HDF5(file, gParams, instrumentDomains);
        
        // Write species metadata on file creation
        if (!ions.empty()) {
            ICARION::io::write_species_metadata(file, ions);
        }
        
        file.close();
    }

    // --- Logging buffers ---
    const int n_ions = static_cast<int>(ions.size());
    const size_t buffer_max = 50; // max timesteps in RAM
    std::vector<double>                times_buffer;
    std::vector<std::vector<IonState>> trajectory_buffer;
    if (do_hdf5) {
        times_buffer.reserve(buffer_max);
        trajectory_buffer.reserve(buffer_max);
    }

    // Space charge solver initialization
    if (gParams.enable_space_charge) {
        gParams.spaceChargeSolver = new SpaceChargeSolver(64,64,64,1e-3,1e-3,1e-3, Vec3{0,0,0});
    }
    
    // Field versioning for live field updates
    std::unordered_map<int, int> field_version_by_domain; // domain_index -> known_version
    for (const auto& dom : instrumentDomains) {
        if (dom.use_grid_field && dom.fieldServer != nullptr) {
            // Initialize with current version
            field_version_by_domain[dom.index] = dom.fieldServer->get_snapshot().version;
        }
    }
    
    // --- RNG per ion ---
    unsigned long long master_seed = 0xCAFEBABE12345678ULL;
    std::vector<EhssRng> rng_by_ion(n_ions);
    for (int i = 0; i < n_ions; ++i) rng_by_ion[i] = EhssRng(master_seed + i);
    if (logger) {
        std::ostringstream msg;
        msg << "RNG initialized for " << n_ions << " ions";
        logger->log(msg.str());
    }

    // Create PhysicsConfig from GlobalParams (used by collision and reaction handlers)
    ICARION::config::PhysicsConfig physics_config;
    physics_config.collision_model = gParams.collisionModel;
    physics_config.enable_reactions = gParams.enable_reactions;
    
    // Create collision handler using factory (Phase 2D refactor)
    std::unique_ptr<ICARION::physics::ICollisionHandler> collision_handler;
    if (gParams.collisionModel != CollisionModel::Friction) {
        // Load geometry map for EHSS if needed (Phase 2E: SSOT)
        ICARION::physics::GeometryMap geometry_map;
        if (gParams.collisionModel == CollisionModel::EHSS && !gParams.geometry_file.empty()) {
            // Extract unique species IDs from ions
            std::unordered_set<std::string> species_ids;
            for (const auto& ion : ions) {
                species_ids.insert(ion.species_id);
            }
            
            // SSOT: Use central geometry loading utility (Phase 2E)
            geometry_map = ICARION::physics::load_geometry_map(species_ids, gParams.geometry_file);
            
            if (logger && !geometry_map.empty()) {
                std::ostringstream msg;
                msg << "Loaded geometry for " << geometry_map.size() << " species";
                logger->log(msg.str());
            }
        }
        
        // Create handler from factory
        collision_handler = ICARION::physics::CollisionHandlerFactory::create(
            physics_config, &geometry_map);
        
        if (logger) {
            std::ostringstream msg;
            msg << "Collision handler created: model=" 
                << static_cast<int>(gParams.collisionModel);
            if (!geometry_map.empty()) {
                msg << ", loaded geometry for " << geometry_map.size() << " species";
            }
            logger->log(msg.str());
        }
    }
    
    // Create reaction handler using factory (Phase 3B refactor)
    // NOTE: Not yet wired into integrate_one_step due to database type mismatch
    // Will be connected in Phase 3C after SpeciesDatabase migration
    std::unique_ptr<ICARION::physics::IReactionHandler> reaction_handler;
    reaction_handler = ICARION::physics::ReactionHandlerFactory::create(physics_config, logger != nullptr);
    
    if (logger) {
        std::ostringstream msg;
        msg << "Reaction handler created: enable_reactions=" << gParams.enable_reactions;
        logger->log(msg.str());
    }

    // --- Solver setup ---
    double t_global = t_start;
    int step_count = 0;
    int log_step_count = 0;
    double write_dt = gParams.write_interval * dt;
    double next_write_time = t_start + write_dt;
    if (do_hdf5) {
        log_step(t_start, t_end, step_count, next_write_time, write_dt, ions,
                 times_buffer, trajectory_buffer, buffer_max, filename);
    }
    if (logger) {
        std::ostringstream msg;
        msg << "setup complete: t_start=" << t_start << " t_end=" << t_end << " dt=" << dt
            << " n_ions=" << n_ions << " collisionModel=" << static_cast<int>(gParams.collisionModel)
            << " enable_space_charge=" << gParams.enable_space_charge << " do_hdf5=" << do_hdf5;
        logger->log(msg.str());
    }
    while (t_global < t_end) {

        // quick diagnostic: print once at start of main loop
        static bool printed_loop_entry = false;
        if (!printed_loop_entry) {
            int active_count = 0;
            for (const auto& ion : ions) if (ion.born && ion.active) ++active_count;
            if (logger) {
                std::ostringstream msg;
                msg << "Entering main loop at t_global=" << t_global << " active_ions=" << active_count;
                logger->log(msg.str());
            }
            printed_loop_entry = true;
        }

        // Check for live field updates 
        for (const auto& dom : instrumentDomains) {
            if (dom.use_grid_field && dom.fieldServer != nullptr) {
                int known_version = field_version_by_domain[dom.index];
                if (dom.fieldServer->has_newer(known_version)) {
                    // Update known version for this domain
                    field_version_by_domain[dom.index] = dom.fieldServer->get_snapshot().version;
                    if (logger) {
                        logger->log("Field updated for domain " + std::to_string(dom.index) + 
                                    " at t=" + std::to_string(t_global) + " s");
                    }
                }
            }
        }

        //Activating newborn ions  
        for (auto& ion : ions) {
            if (!ion.born && t_global >= ion.birth_time_s) {
                ion.born = true;
                ion.active = true; 
                ion.t = ion.birth_time_s;
                ion.dt = dt; 
                ion.current_domain_index = find_domain_index(ion.pos, instrumentDomains);
            }
        }

        //parallel integrating over active ions
        #pragma omp parallel for if (gParams.parallelization)
        for (int idx = 0; idx < n_ions; ++idx) {
            IonState& y = ions[idx];
            EhssRng& local_rng = rng_by_ion[idx];
            if (!y.active || !y.born) continue;
            
            while(y.t < t_global + dt && y.active){
                //  Domain check 
                int idx = find_domain_index(y.pos, instrumentDomains);

                if (idx < 0) {
                    y.active = false;
                    continue; // Ion has left all domains
                }
                const InstrumentDomain& dom = instrumentDomains[idx];
                y.dt = dt;
                 // Update domain-specific properties
                if (y.current_domain_index != dom.index) {
                    y.domain_neutral_mass_kg = dom.env.neutral_mass_kg;
                    y.domain_temperature_K   = dom.env.temperature_K;
                    y.domain_particle_density_m3 = dom.env.particle_density_m_3;
                    y.domain_gas_velocity_m_s    = Vec3(dom.env.gas_velocity_m_s.x,
                                                    dom.env.gas_velocity_m_s.y,
                                                    dom.env.gas_velocity_m_s.z);
                    y.current_domain_index   = dom.index;
                }

                // Transform to local domain coordinates
                IonState y_loc = y;
                y_loc.pos = dom.rotation_global_to_local * (y.pos - dom.geom.origin_m);
                y_loc.vel = dom.rotation_global_to_local * y.vel;
                
                // Collision and reaction physics 
                if (collision_handler) {
                    // Build environment config from domain and ion state
                    ICARION::config::EnvironmentConfig env;
                    env.temperature_K = y_loc.domain_temperature_K;
                    env.pressure_Pa = y_loc.domain_particle_density_m3 * BOLTZMANN_CONSTANT * y_loc.domain_temperature_K;
                    env.gas_species = "He";  // TODO: Get from domain or ion state
                    env.gas_velocity_m_s = y_loc.domain_gas_velocity_m_s;
                    env.compute_derived_properties();
                    
                    collision_handler->handle_collision(y_loc, y.dt, local_rng, env);
                }
                handle_reaction(y_loc, local_rng, y.dt, gParams, speciesDB, reaction_list);
                Vec3 pos_before = y_loc.pos;
                integrate_one_step(y_loc, gParams, dom, speciesDB, reaction_list, local_rng, rk45, ions);
                Vec3 pos_after = y_loc.pos;

                // Check if ion transitions from domain to another and deactivate if it 
                // does not cross the end aperture of the old domain
                if (pos_before.z < dom.geom.length_m && pos_after.z >= dom.geom.length_m) {
                    // Ion is exiting the domain at the end
                    check_aperture_crossing(y_loc, dom, instrumentDomains, pos_before, pos_after);
                }

                // Back-transform  
                y.pos = dom.rotation_local_to_global * y_loc.pos + dom.geom.origin_m;
                y.vel = dom.rotation_local_to_global * y_loc.vel;
                y.t = y_loc.t;
                y.species_id              = y_loc.species_id;
                y.mass_kg                 = y_loc.mass_kg;
                y.reduced_mobility_cm2_Vs = y_loc.reduced_mobility_cm2_Vs;
                y.ion_charge_C            = y_loc.ion_charge_C;
                y.CCS_m2= y_loc.CCS_m2;
                }
        }

        // Exit check
        if (std::all_of(ions.begin(), ions.end(), [](const IonState& ion){ return ion.born && !ion.active; })) {
            if (logger) {
                std::ostringstream msg;
                msg << "All ions exited at t = " << t_global << " s";
                logger->log(msg.str());
            }
            break;
        }
        if (do_hdf5 && t_global >= next_write_time){
            log_step(t_global, t_end, step_count, next_write_time, write_dt, ions, times_buffer, trajectory_buffer, buffer_max, filename);
        }

        t_global += dt;

        // Update space charge solver with new ion positions
        if (gParams.enable_space_charge && gParams.spaceChargeSolver != nullptr) {
            gParams.spaceChargeSolver->update(ions);
        }

        // Progress logging
        if (logger && gParams.sim_time_steps > 0)
        {
            const int total_steps = gParams.sim_time_steps;
            const int log_interval = std::max(1, total_steps / 10);

            if (log_step_count % log_interval == 0 || log_step_count == total_steps - 1) {
                double percent = 100.0 * static_cast<double>(log_step_count) / static_cast<double>(total_steps);
                std::ostringstream msg;
                msg << std::fixed << std::setprecision(0)
                    << percent << "% completed (t = " << t_global * 1000.0 << " ms)";
                logger->log(msg.str());
            }
        }
    // Optional verbose position logging removed to keep test output clean
        log_step_count++;

    }

    // flush remaining logs
    if (do_hdf5 && !times_buffer.empty()) {
        if (times_buffer.back() < t_end) {
            times_buffer.push_back(t_end);
            trajectory_buffer.push_back(ions);
        }
        ICARION::io::append_to_HDF5(filename, times_buffer, trajectory_buffer);
    }
    if (logger)
        logger->log("Integration complete. Writing final data...");
    return ions;
}

/**
 * @brief Compute weighted RMS error norm between 5th- and 4th-order RK45 solutions.
 *
 * Used for adaptive step size control: error <= 1.0 means step is accepted.
 */
double rk_error_norm(const IonState& y5, const IonState& y4, const IonState& y, const RK45Settings& rk45) {
    double s = 0.0;
    s += std::pow((y5.pos.x - y4.pos.x) /
                          (rk45.absTol + rk45.relTol * std::max(std::fabs(y.pos.x), std::fabs(y5.pos.x))),
                  2);
    s += std::pow((y5.pos.y - y4.pos.y) /
                          (rk45.absTol + rk45.relTol * std::max(std::fabs(y.pos.y), std::fabs(y5.pos.y))),
                  2);
    s += std::pow((y5.pos.z - y4.pos.z) /
                          (rk45.absTol + rk45.relTol * std::max(std::fabs(y.pos.z), std::fabs(y5.pos.z))),
                  2);
    s += std::pow((y5.vel.x - y4.vel.x) /
                          (rk45.absTol + rk45.relTol * std::max(std::fabs(y.vel.x), std::fabs(y5.vel.x))),
                  2);
    s += std::pow((y5.vel.y - y4.vel.y) /
                          (rk45.absTol + rk45.relTol * std::max(std::fabs(y.vel.y), std::fabs(y5.vel.y))),
                  2);
    s += std::pow((y5.vel.z - y4.vel.z) /
                          (rk45.absTol + rk45.relTol * std::max(std::fabs(y.vel.z), std::fabs(y5.vel.z))),
                  2);
    return std::sqrt(s / 6.0);
}

/**
 * @brief Perform a single Dormand–Prince RK45 step.
 *
 * Computes both 5th-order (y5) and 4th-order embedded (y4) solutions.
 * Difference is used for local error estimation.
 */
void rk45_dp_step(double t, const IonState& y, double h, const GlobalParams& gParams, 
                  const InstrumentDomain& dom, IonState& y5,
                  IonState& y4, const std::vector<IonState>& current_ions) {
    
                    // Dormand–Prince coefficients               
    const double c2 = 1.0 / 5.0;
    const double c3 = 3.0 / 10.0;
    const double c4 = 4.0 / 5.0;
    const double c5 = 8.0 / 9.0;
    const double c6 = 1.0;
    const double c7 = 1.0;

    const double a21 = 1.0 / 5.0;
    const double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
    const double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
    const double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0, a53 = 64448.0 / 6561.0,
                 a54 = -212.0 / 729.0;
    const double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0,
                 a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
    const double a71 = 35.0 / 384.0, a72 = 0.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0,
                 a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;
    const double b1 = 35.0 / 384.0, b2 = 0.0, b3 = 500.0 / 1113.0, b4 = 125.0 / 192.0,
                 b5 = -2187.0 / 6784.0, b6 = 11.0 / 84.0, b7 = 0.0;
    const double bs1 = 5179.0 / 57600.0, bs2 = 0.0, bs3 = 7571.0 / 16695.0, bs4 = 393.0 / 640.0,
                 bs5 = -92097.0 / 339200.0, bs6 = 187.0 / 2100.0, bs7 = 1.0 / 40.0;

    // k1..k7
    // For RK45 steps prefer provider if possible (create local snapshot/provider)
    std::unique_ptr<IFieldProvider> provider_owner_all;
    FieldSnapshot prov_snapshot_all;
    if (dom.use_grid_field && dom.fieldServer != nullptr) {
        prov_snapshot_all = dom.fieldServer->get_snapshot();
        provider_owner_all.reset(new GridFieldProvider(&prov_snapshot_all));
    } else if (dom.fieldArrayLoaded) {
        provider_owner_all.reset(new GridFieldProvider(&dom.fieldArray));
    }
    const IFieldProvider* provider_all_ptr = provider_owner_all ? provider_owner_all.get() : nullptr;

    IonState k1 = compute_accelerations(t, y, gParams, dom, provider_all_ptr, current_ions) * h;
    IonState k2 = compute_accelerations(t + c2 * h, y + k1 * a21, gParams, dom, provider_all_ptr, current_ions) * h;
    IonState k3 = compute_accelerations(t + c3 * h, y + k1 * a31 + k2 * a32, gParams, dom, provider_all_ptr, current_ions) * h;
    IonState k4 = compute_accelerations(t + c4 * h, y + k1 * a41 + k2 * a42 + k3 * a43, gParams, dom, provider_all_ptr, current_ions) * h;
    IonState k5 = compute_accelerations(t + c5 * h, y + k1 * a51 + k2 * a52 + k3 * a53 + k4 * a54, gParams, dom, provider_all_ptr, current_ions) * h;
    IonState k6 = compute_accelerations(t + c6 * h, y + k1 * a61 + k2 * a62 + k3 * a63 + k4 * a64 + k5 * a65,
                             gParams, dom, provider_all_ptr, current_ions) *
                  h;
    IonState k7 = compute_accelerations(t + c7 * h,
                             y + k1 * a71 + k2 * a72 + k3 * a73 + k4 * a74 + k5 * a75 + k6 * a76,
                             gParams, dom, provider_all_ptr, current_ions) *
                  h;

    // solutions
    y5 = y + k1 * b1 + k2 * b2 + k3 * b3 + k4 * b4 + k5 * b5 + k6 * b6 + k7 * b7;
    y4 = y + k1 * bs1 + k2 * bs2 + k3 * bs3 + k4 * bs4 + k5 * bs5 + k6 * bs6 + k7 * bs7;
}

/**
 * @brief Load molecular geometry for a given molecule.
 *
 * Reads a geometry file, selects the first molecule matching @p targetName, 
 * and extracts atom positions and radii. Radii are computed as half of the 
 * Lennard–Jones sigma parameter.
 *
 * @param[in]  filename    Path to the geometry file.
 * @param[in]  targetName  Name of the molecule to load.
 * @param[out] h_centers   Vector to be filled with atom positions (meters).
 * @param[out] h_radii     Vector to be filled with atom radii (meters).
 *
 * @throws std::runtime_error if no molecules are found.
 *
 * @note This function only reads the first matching molecule from the file.
 */
// NOTE: load_geometry() removed in Phase 2D refactor
// Geometry loading now done via load_geometry_from_file() in CollisionHandlerFactory
// See src/core/physics/collisions/utils.h -> SSOT Violation?

// NOTE: handle_collision() removed in Phase 2D refactor
// Collision handling now done via CollisionHandlerFactory and ICollisionHandler interface
// See CollisionHandlerFactory::create() and handler->handle_collision() in integrate_trajectory()

/**
 * @brief Apply stochastic chemical reactions to an ion.
 *
 * Checks whether a reaction occurs for the given ion based on reaction rates
 * and timestep. Updates species, mass, mobility, charge, and CCS if a reaction occurs.
 *
 * @param[in,out] y             Ion state to update if a reaction occurs.
 * @param[in,out] rng           Random number generator for reproducibility.
 * @param[in]     dt            Current timestep [s].
 * @param[in]     gParams       Simulation parameters (enable_reactions, temperature, etc.).
 * @param[in]     speciesDB     Species database.
 * @param[in]     reaction_list List of reactions with rate constants and stoichiometry.
 *
 * @note
 * - Reactions are stochastic: P = 1 - exp(-k_eff * dt)
 * - Effective rate k_eff includes contributions from reaction orders and neutral concentration.
 * - Charge and CCS are updated if the product differs from the reactant.
 */
void handle_reaction(IonState& y, EhssRng& rng, double dt, const GlobalParams& gParams,
                     const ICARION::io::SpeciesDatabase& speciesDB,
                     const std::vector<ReactionEntry>& reaction_list) {
    if (!gParams.enable_reactions) return;

    const std::string current_species = y.species_id;

    for (const auto& rxn : reaction_list) {
        if (rxn.reactant != current_species) continue;
        
        // Calculate effective rate constant
        // k_eff starts as k in [m³/s] from JSON input (converted from cm³/s)
        double k_eff = rxn.rate_constant; // [m³/s]
        
        // Apply concentration dependence for reaction order
        // neutral_concentration is volume fraction (dimensionless)
        // For 2nd order: k_eff [m³/s] * n [m⁻³] = k_total [s⁻¹]
        double neutral_density_m3 = y.domain_particle_density_m3; // [m⁻³]
        
        for (const auto& term : rxn.order) {
            // rxn.neutral_concentration historically can be either
            // - a dimensionless volume fraction (0..1) which must be multiplied
            //   with the domain number density to get an absolute number density [m^-3], or
            // - already an absolute number density [m^-3] (tests/examples sometimes provide this).
            // Detect which by size: values >> 1 (e.g. >1e6) are considered number densities.
            double conc_density = 0.0; // [m^-3]
            if (rxn.neutral_concentration > 1e6) {
                // already a number density
                conc_density = rxn.neutral_concentration;
            } else {
                // treat as volume fraction (dimensionless)
                conc_density = rxn.neutral_concentration * neutral_density_m3; // [m^-3]
            }
            // guard against negative/zero
            if (conc_density < 0.0) conc_density = 0.0;
            k_eff *= std::pow(conc_density, term.exponent);
        }

        
        // Reaction probability: P = 1 - exp(-k_total * dt)
        // where k_total is in [s⁻¹] and dt in [s]
        double P = 1.0 - std::exp(-k_eff * dt);
        
        if (rng.uniform01() < P) {
            // Update ion properties
            if (!speciesDB.has(rxn.product)) {
                std::cerr << "Warning: Product species '" << rxn.product << "' not found in database\n";
                continue;
            }
            const auto& prod = speciesDB.get(rxn.product);
            y.species_id              = rxn.product;
            y.mass_kg                 = prod.mass_kg;
            y.reduced_mobility_cm2_Vs = prod.mobility_m2Vs * 1e4;  // m²/Vs to cm²/Vs
            y.ion_charge_C            = prod.charge_C;
            y.CCS_m2 = 3.0 / 16.0 / LOSCHMIDT_CONSTANT * y.ion_charge_C *
                       std::sqrt(2.0 * M_PI /
                                 (BOLTZMANN_CONSTANT * y.domain_temperature_K *
                                  (y.mass_kg * y.domain_neutral_mass_kg) /
                                  (y.mass_kg + y.domain_neutral_mass_kg))) *
                       1.0 / (y.reduced_mobility_cm2_Vs * 1e-4);
        }
    }
}

//--- Aperture-crossing check between domains ---
bool check_aperture_crossing(IonState ion, const InstrumentDomain& dom,
                                   const std::vector<InstrumentDomain>& domains,
                                   const Vec3& local_old, const Vec3& local_new)
{
    // no aperture defined
    if (dom.geom.end_aperture_m <= 0.0) return true;

    // Did we cross the aperture plane?
    const double z_ap = dom.geom.length_m;
    if ((local_old.z < z_ap && local_new.z >= z_ap) ||
        (local_old.z > z_ap && local_new.z <= z_ap))
    {
        double alpha = (z_ap - local_old.z) / (local_new.z - local_old.z);
        Vec3 cross_local = local_old + (local_new - local_old) * alpha;
        double r_cross = std::sqrt(cross_local.x * cross_local.x +
                                   cross_local.y * cross_local.y);

        if (r_cross <= dom.geom.end_aperture_m) {
            // allowed to pass
            return true;
        } else {
            // blocked
            ion.active = false;
            return false;
        }
    }

    // No crossing this step
    return true;
};

}  // namespace trajectory
}  // namespace ICARION
