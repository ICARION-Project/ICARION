// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "SimulationEngine.h"
#include "core/utils/safety/numericalSafetyGuards.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include "core/utils/mathUtils.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ICARION {
namespace integrator {

SimulationEngine::SimulationEngine(
    const config::FullConfig& config,
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries,
    std::shared_ptr<IIntegrationStrategy> integrator,
    std::shared_ptr<physics::ICollisionHandler> collision_handler,
    std::shared_ptr<physics::IReactionHandler> reaction_handler
) : config_(config),
    force_registries_(std::move(force_registries)),
    integrator_(integrator),
    collision_handler_(collision_handler),
    reaction_handler_(reaction_handler)
{
    // Validate force registries (must have one per domain)
    if (force_registries_.size() != config_.domains.size()) {
        throw std::invalid_argument(
            "SimulationEngine: force_registries size (" + std::to_string(force_registries_.size()) + 
            ") must match domains size (" + std::to_string(config_.domains.size()) + ")"
        );
    }
    
    for (size_t i = 0; i < force_registries_.size(); ++i) {
        if (!force_registries_[i]) {
            throw std::invalid_argument(
                "SimulationEngine: ForceRegistry[" + std::to_string(i) + "] cannot be null"
            );
        }
    }
    
    if (!integrator_) {
        throw std::invalid_argument("SimulationEngine: IntegrationStrategy cannot be null");
    }
    
    // Create domain manager
    domain_manager_ = std::make_unique<DomainManager>(config_.domains);
    
    // Create output manager
    std::string hdf5_path = config_.output.folder + "/" + config_.output.trajectory_file;
    std::string log_path = "";  // TODO: Add log file to OutputConfig
    if (config_.output.print_progress) {
        log_path = config_.output.folder + "/simulation.log";
    }
    
    output_manager_ = std::make_unique<OutputManager>(
        hdf5_path,
        log_path,
        config_.simulation.dt_s * config_.simulation.write_interval,  // Write interval in seconds
        50  // Buffer max
    );
    
    // Initialize numerical safety logger if enabled
    if (config_.simulation.enable_safety_logging) {
        std::string safety_log = config_.output.folder + "/numerical_safety.log";
        safety::NumericalSafetyLogger::getInstance(safety_log, config_.simulation.verbose_safety);
        
        // Configure logger (enable_logging, verbose_mode, buffer_size, max_history)
        safety::NumericalSafetyLogger::getInstance().configure(
            true,  // enable_logging
            config_.simulation.verbose_safety,  // verbose_mode
            1000,  // buffer_size (default)
            10000  // max_history (default)
        );
        
        output_manager_->log_progress("Numerical safety logging enabled: " + safety_log);
    }
}

void SimulationEngine::initialize(const std::vector<IonState>& ions) {
    // Initialize output system
    output_manager_->initialize(config_, ions);
    
    // Log initialization
    output_manager_->log_progress("Simulation engine initialized");
    
    std::ostringstream msg;
    msg << "Configuration: " << ions.size() << " ions, "
        << config_.domains.size() << " domains, "
        << "dt = " << config_.simulation.dt_s * 1e9 << " ns, "
        << "t_total = " << config_.simulation.total_time_s * 1e6 << " μs";
    output_manager_->log_progress(msg.str());
}

void SimulationEngine::apply_ion_birth(std::vector<IonState>& ions, double t) {
    for (auto& ion : ions) {
        if (!ion.born && ion.birth_time_s <= t) {
            ion.born = true;
            ion.active = true;
            ion.t = t;
        }
    }
}

void SimulationEngine::process_timestep(std::vector<IonState>& ions, double dt) {
    const int n_ions = static_cast<int>(ions.size());
    
    // Ion-based RNG for reproducibility (independent of OpenMP scheduling)
    // Each ion gets its own RNG seeded deterministically from: base_seed + ion_index
    // This ensures:
    // - Reproducible results regardless of thread count or scheduling
    // - Same ion always sees same random sequence
    // - Thread-safe (each thread accesses different ion RNG)
    std::vector<EhssRng> rng_by_ion;
    rng_by_ion.reserve(n_ions);
    for (int i = 0; i < n_ions; ++i) {
        uint64_t ion_seed = config_.simulation.rng_seed + static_cast<uint64_t>(i);
        rng_by_ion.emplace_back(ion_seed);
    }
    
    // Parallel ion processing (OpenMP if enabled)
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n_ions; ++i) {
            IonState& ion = ions[i];
            EhssRng& ion_rng = rng_by_ion[i];  // Ion-specific RNG
            
            // Skip inactive ions (still process ions waiting to be born)
            if (!ion.active) {
                continue;
            }
        
        // 1. Find current domain
        int domain_idx = domain_manager_->find_domain_index(ion.pos);
        if (domain_idx < 0) {
            ion.active = false;  // Ion left all domains
            continue;
        }
        
        const auto& domain_config = config_.domains[domain_idx];
        
        // 2. Update domain-specific properties (gas properties, environment)
        if (ion.current_domain_index != domain_idx) {
            domain_manager_->update_domain_properties(ion, domain_idx);
        }
        
        // 3. Transform to local domain coordinates
        Vec3 pos_local = domain_manager_->global_to_local_pos(ion.pos, domain_idx);
        Vec3 vel_local = domain_manager_->global_to_local_vel(ion.vel, domain_idx);
        
        // Store pre-integration position for aperture crossing check
        Vec3 pos_before = pos_local;
        
        // 4. Handle collisions (if enabled)
        if (collision_handler_) {
            IonState ion_local = ion;
            ion_local.pos = pos_local;
            ion_local.vel = vel_local;
            
            // Thread-safe RNG
            collision_handler_->handle_collision(ion_local, dt, ion_rng, domain_config.environment);
            
            pos_local = ion_local.pos;
            vel_local = ion_local.vel;
        }
        
        // 5. Handle reactions (if enabled)
        // SSOT: Direct access to config databases (no parameter copies!)
        if (reaction_handler_ && !config_.reaction_db.reactions.empty()) {
            IonState ion_local = ion;
            ion_local.pos = pos_local;
            ion_local.vel = vel_local;
            
            // Call reaction handler with SSOT databases
            bool reaction_occurred = reaction_handler_->handle_reaction(
                ion_local,
                dt,
                ion_rng,
                config_.reaction_db,  // ReactionDatabase (SSOT)
                config_.species_db,   // SpeciesDatabase (SSOT)
                domain_config.environment  // Temperature, density, etc.
            );
            
            // If reaction changed species, update ion properties
            if (reaction_occurred) {
                ion.species_id = ion_local.species_id;
                ion.mass_kg = ion_local.mass_kg;
                ion.ion_charge_C = ion_local.ion_charge_C;
                ion.CCS_m2 = ion_local.CCS_m2;
                ion.reduced_mobility_cm2_Vs = ion_local.reduced_mobility_cm2_Vs;
            }
            
            pos_local = ion_local.pos;
            vel_local = ion_local.vel;
        }
        
        // 6. Get domain-specific ForceRegistry
        const auto& force_registry = force_registries_[domain_idx];
        
        // 7. Compute forces (ForceRegistry with domain context)
        IonState ion_local = ion;
        ion_local.pos = pos_local;
        ion_local.vel = vel_local;
        
        physics::ForceContext ctx;  // TODO: Populate with field provider
        Vec3 total_force = force_registry->compute_total_force(ion_local, current_time_, ctx);
        
        // 8. Integrate trajectory (IIntegrationStrategy)
        // Note: IIntegrationStrategy::step() computes forces internally via ForceRegistry
        // ForceRegistry now knows its domain, so we don't need to pass domain_config
        integrator_->step(ion_local, current_time_, dt, *force_registry, ions);
        
        pos_local = ion_local.pos;
        vel_local = ion_local.vel;
        
        // 9. Check if ion left domain (boundary collision detection)
        Vec3 pos_after = pos_local;
        
        // First check aperture crossing (domain exit at z=length_m)
        if (pos_after.z >= domain_config.geometry.length_m && pos_before.z < domain_config.geometry.length_m) {
            domain_manager_->check_aperture_crossing(ion, domain_idx, pos_before, pos_after);
            
            // If blocked by aperture, terminate at boundary
            if (!ion.active) {
                domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
                continue;  // Ion absorbed, skip further processing
            }
            
            // Ion passed through aperture
            // Check if this is the last domain (no next domain available)
            bool is_last_domain = (domain_idx == static_cast<int>(config_.domains.size()) - 1);
            
            if (is_last_domain) {
                // No next domain - ion exits simulation
                ion.active = false;
                domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
                continue;
            }
            
            // Multi-domain: Transform to global and let next timestep find new domain
            ion.pos = domain_manager_->local_to_global_pos(pos_local, domain_idx);
            ion.vel = domain_manager_->local_to_global_vel(vel_local, domain_idx);
            ion.t += dt;
            continue;  // Skip remaining checks, next step will find new domain
        }
        
        // Then check all other boundaries (radial wall, entrance plane)
        // Note: Allow small tolerance at z_max for multi-domain transitions
        bool still_inside = (pos_after.z >= -DOMAIN_BOUNDARY_EPSILON);
        
        // Check z_max: strict for last domain, tolerant for others
        bool is_last_domain = (domain_idx == static_cast<int>(config_.domains.size()) - 1);
        if (is_last_domain) {
            // Last domain: strict boundary (no tolerance for exit)
            still_inside = still_inside && (pos_after.z < domain_config.geometry.length_m);
        } else {
            // Intermediate domain: allow small tolerance for aperture crossing
            still_inside = still_inside && (pos_after.z <= domain_config.geometry.length_m + DOMAIN_BOUNDARY_EPSILON);
        }
        
        // Check radial boundary
        if (still_inside) {
            double r = std::sqrt(pos_after.x*pos_after.x + pos_after.y*pos_after.y);
            still_inside = (r <= domain_config.geometry.radius_m + DOMAIN_BOUNDARY_EPSILON);
        }
        
        if (!still_inside) {
            // Ion left domain (hit wall or exited entrance/exit)
            domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
            continue;  // Ion absorbed, skip transform
        }
        
        // 9. Transform back to global coordinates
        ion.pos = domain_manager_->local_to_global_pos(pos_local, domain_idx);
        ion.vel = domain_manager_->local_to_global_vel(vel_local, domain_idx);
        
        // 10. Update ion time
        ion.t += dt;
        
        // 11. Numerical safety checks
        bool position_valid = ICARION::safety::is_finite(ion.pos);
        bool velocity_valid = ICARION::safety::is_finite(ion.vel);
        
        if (!position_valid || !velocity_valid) {
            // Log detailed violation if safety logging is enabled
            if (config_.simulation.enable_safety_logging) {
                safety::ViolationEvent event;
                event.type = !position_valid ? 
                    (std::isnan(ion.pos.x + ion.pos.y + ion.pos.z) ? 
                        safety::ViolationType::NAN_POSITION : safety::ViolationType::INF_POSITION) :
                    (std::isnan(ion.vel.x + ion.vel.y + ion.vel.z) ? 
                        safety::ViolationType::NAN_VELOCITY : safety::ViolationType::INF_VELOCITY);
                
                event.timestamp = std::chrono::steady_clock::now();
                event.ion_index = i;
                event.step_number = current_step_;
                event.simulation_time = ion.t;
                event.timestep = dt;
                event.position = ion.pos;
                event.velocity = ion.vel;
                event.violation_context = "Post-integration state check in domain " + std::to_string(domain_idx);
                event.violation_magnitude = !position_valid ? norm(ion.pos) : norm(ion.vel);
                event.recovery_attempted = false;
                event.recovery_successful = false;
                
                safety::NumericalSafetyLogger::getInstance().logViolation(event);
            }
            
            // Deactivate ion
            ion.active = false;
            
            // Log to standard output if not using safety logger
            if (!config_.simulation.enable_safety_logging) {
                std::cerr << "Warning: Ion " << i << " has invalid state (";
                if (!position_valid) std::cerr << "NaN/Inf in position";
                if (!position_valid && !velocity_valid) std::cerr << ", ";
                if (!velocity_valid) std::cerr << "NaN/Inf in velocity";
                std::cerr << ") at t = " << ion.t << " s, deactivating" << std::endl;
            }
        }
        
        // Optional bounds checking (if enabled)
        if (config_.simulation.safety_checks.enable_bounds_checks && 
            (position_valid && velocity_valid)) {
            
            double pos_mag = norm(ion.pos);
            double vel_mag = norm(ion.vel);
            
            bool bounds_violated = false;
            safety::ViolationType violation_type;
            
            if (pos_mag > config_.simulation.safety_checks.max_position_m) {
                bounds_violated = true;
                violation_type = safety::ViolationType::BOUNDS_POSITION;
            } else if (vel_mag > config_.simulation.safety_checks.max_velocity_ms) {
                bounds_violated = true;
                violation_type = safety::ViolationType::BOUNDS_VELOCITY;
            }
            
            if (bounds_violated) {
                if (config_.simulation.enable_safety_logging) {
                    safety::ViolationEvent event;
                    event.type = violation_type;
                    event.timestamp = std::chrono::steady_clock::now();
                    event.ion_index = i;
                    event.step_number = current_step_;
                    event.simulation_time = ion.t;
                    event.timestep = dt;
                    event.position = ion.pos;
                    event.velocity = ion.vel;
                    event.violation_context = "Bounds check exceeded in domain " + std::to_string(domain_idx);
                    event.violation_magnitude = (violation_type == safety::ViolationType::BOUNDS_POSITION) 
                        ? pos_mag : vel_mag;
                    event.recovery_attempted = false;
                    event.recovery_successful = false;
                    
                    safety::NumericalSafetyLogger::getInstance().logViolation(event);
                }
                
                if (config_.simulation.safety_checks.throw_on_violation) {
                    throw std::runtime_error("Bounds violation for ion " + std::to_string(i) + 
                                           " at t=" + std::to_string(ion.t));
                }
                
                ion.active = false;
            }
        }
        }  // End of parallel for loop
    }  // End of parallel region
}

bool SimulationEngine::should_continue(const std::vector<IonState>& ions, double t) const {
    // Stop if time exceeded
    if (t >= config_.simulation.total_time_s) {
        return false;
    }
    
    // Continue if any ion is active OR waiting to be born
    const double t_max = config_.simulation.total_time_s;
    bool any_active_or_waiting = std::any_of(ions.begin(), ions.end(), 
        [t_max](const IonState& ion) { 
            return (ion.active && ion.born) || (!ion.born && ion.birth_time_s <= t_max); 
        });
    
    return any_active_or_waiting;
}

void SimulationEngine::log_progress(double t) {
    const int total_steps = static_cast<int>(config_.simulation.total_time_s / config_.simulation.dt_s);
    const int log_interval = std::max(1, total_steps / 10);
    
    if (current_step_ % log_interval == 0 || current_step_ == total_steps - 1) {
        double percent = 100.0 * t / config_.simulation.total_time_s;
        
        std::ostringstream msg;
        msg << std::fixed << std::setprecision(0)
            << percent << "% completed (t = " << t * 1e3 << " ms, step " 
            << current_step_ << "/" << total_steps << ")";
        
        output_manager_->log_progress(msg.str());
    }
}

std::vector<IonState> SimulationEngine::run(std::vector<IonState>& ions) {
    // 1. Initialize subsystems
    initialize(ions);
    
    // 2. Main time loop
    const double dt = config_.simulation.dt_s;
    current_time_ = 0.0;
    current_step_ = 0;
    
    output_manager_->log_progress("Starting main simulation loop");
    
    while (should_continue(ions, current_time_)) {
        // Apply ion birth logic (delayed emission)
        apply_ion_birth(ions, current_time_);
        
        // Process one timestep
        process_timestep(ions, dt);
        
        // Log trajectory snapshot (auto-flush if needed)
        output_manager_->log_step(current_time_, ions);
        
        // Update time and step counter
        current_time_ += dt;
        current_step_++;
        
        // Progress logging (every 10%)
        log_progress(current_time_);
    }
    
    // 3. Finalization
    output_manager_->log_progress("Simulation completed");
    
    // Count active ions
    size_t active_count = std::count_if(ions.begin(), ions.end(),
        [](const IonState& ion) { return ion.active && ion.born; });
    
    std::ostringstream msg;
    msg << "Final state: " << active_count << "/" << ions.size() << " ions active";
    output_manager_->log_progress(msg.str());
    
    // Flush output and write completion metadata
    output_manager_->finalize(current_time_, ions);
    
    // Generate numerical safety report if logging was enabled
    if (config_.simulation.enable_safety_logging) {
        std::string report_file = config_.output.folder + "/numerical_safety_report.txt";
        safety::NumericalSafetyLogger::getInstance().generateSafetyReport(report_file);
        
        // Log summary statistics
        auto stats = safety::NumericalSafetyLogger::getInstance().getStatistics();
        std::ostringstream safety_msg;
        safety_msg << "Numerical safety: " << stats.total_violations << " violations detected";
        if (stats.recovery_attempts > 0) {
            safety_msg << ", " << stats.successful_recoveries << "/" 
                      << stats.recovery_attempts << " recoveries successful";
        }
        output_manager_->log_progress(safety_msg.str());
        output_manager_->log_progress("Safety report written: " + report_file);
    }
    
    return ions;
}

}  // namespace integrator
}  // namespace ICARION
