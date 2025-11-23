// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "SimulationEngine.h"
#include "core/utils/safety/numericalSafetyGuards.h"
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
    std::shared_ptr<physics::ForceRegistry> force_registry,
    std::shared_ptr<IIntegrationStrategy> integrator,
    std::shared_ptr<physics::ICollisionHandler> collision_handler,
    std::shared_ptr<physics::IReactionHandler> reaction_handler
) : config_(config),
    force_registry_(force_registry),
    integrator_(integrator),
    collision_handler_(collision_handler),
    reaction_handler_(reaction_handler)
{
    if (!force_registry_) {
        throw std::invalid_argument("SimulationEngine: ForceRegistry cannot be null");
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
    
    // Parallel ion processing (OpenMP if enabled)
    #pragma omp parallel for if(config_.simulation.enable_openmp) schedule(dynamic)
    for (int i = 0; i < n_ions; ++i) {
        IonState& ion = ions[i];
        
        // Skip inactive or unborn ions
        if (!ion.active || !ion.born) {
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
            
            // TODO: Thread-local RNG (currently using dummy - unsafe in parallel)
            static EhssRng dummy_rng(42);  // FIXME: Thread-unsafe!
            collision_handler_->handle_collision(ion_local, dt, dummy_rng, domain_config.environment);
            
            pos_local = ion_local.pos;
            vel_local = ion_local.vel;
        }
        
        // 5. Handle reactions (if enabled)
        if (reaction_handler_) {
            IonState ion_local = ion;
            ion_local.pos = pos_local;
            ion_local.vel = vel_local;
            
            // TODO: Thread-local RNG
            // reaction_handler_->handle_reaction(ion_local, dt, ...);
            
            pos_local = ion_local.pos;
            vel_local = ion_local.vel;
        }
        
        // 6. Compute forces (ForceRegistry)
        IonState ion_local = ion;
        ion_local.pos = pos_local;
        ion_local.vel = vel_local;
        
        physics::ForceContext ctx;  // TODO: Populate with field provider
        Vec3 total_force = force_registry_->compute_total_force(ion_local, current_time_, ctx);
        
        // 7. Integrate trajectory (IIntegrationStrategy)
        // Note: IIntegrationStrategy::step() computes forces internally via ForceRegistry
        // So we don't need to pass acceleration separately
        integrator_->step(ion_local, dt, current_time_, *force_registry_, domain_config, ions);
        
        pos_local = ion_local.pos;
        vel_local = ion_local.vel;
        
        // 8. Check if ion left domain (boundary collision detection)
        Vec3 pos_after = pos_local;
        
        // First check aperture crossing (domain exit at z=length_m)
        if (pos_after.z >= domain_config.geometry.length_m && pos_before.z < domain_config.geometry.length_m) {
            domain_manager_->check_aperture_crossing(ion, domain_idx, pos_before, pos_after);
            
            // If blocked by aperture, terminate at boundary
            if (!ion.active) {
                domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
                continue;  // Ion absorbed, skip further processing
            }
        }
        
        // Then check all other boundaries (radial wall, entrance plane)
        // Note: is_inside_domain checks domain AFTER transform back to global
        // So we check in local coordinates first
        bool still_inside = (pos_after.z >= -DOMAIN_BOUNDARY_EPSILON && 
                            pos_after.z < domain_config.geometry.length_m);
        if (still_inside) {
            double r = std::sqrt(pos_after.x*pos_after.x + pos_after.y*pos_after.y);
            still_inside = (r < domain_config.geometry.radius_m);
        }
        
        if (!still_inside) {
            // Ion left domain (hit wall or exited entrance)
            domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
            continue;  // Ion absorbed, skip transform
        }
        
        // 9. Transform back to global coordinates
        ion.pos = domain_manager_->local_to_global_pos(pos_local, domain_idx);
        ion.vel = domain_manager_->local_to_global_vel(vel_local, domain_idx);
        
        // 10. Update ion time
        ion.t += dt;
        
        // 11. Safety check (NaN/Inf detection)
        if (!ICARION::safety::is_finite(ion.pos) || !ICARION::safety::is_finite(ion.vel)) {
            std::cerr << "Warning: Ion " << i << " has invalid state (NaN/Inf) at t = " 
                      << ion.t << " s, deactivating" << std::endl;
            ion.active = false;
        }
    }
}

bool SimulationEngine::should_continue(const std::vector<IonState>& ions, double t) const {
    // Stop if time exceeded
    if (t >= config_.simulation.total_time_s) {
        return false;
    }
    
    // Stop if all ions inactive
    bool any_active = std::any_of(ions.begin(), ions.end(), 
        [](const IonState& ion) { return ion.active && ion.born; });
    
    return any_active;
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
    
    return ions;
}

}  // namespace integrator
}  // namespace ICARION
