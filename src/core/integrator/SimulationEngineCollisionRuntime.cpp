// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationEngine.h"

#include "core/utils/Profiler.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <sstream>
#include <stdexcept>

namespace ICARION {
namespace integrator {

void SimulationEngine::perform_collisions(core::IonEnsemble& ensemble,
                                          const std::vector<double>& dt_used_per_ion,
                                          const std::vector<int>& domain_indices) {
    if (!collision_handler_ || domain_indices.empty()) {
        return;
    }

    PROFILE_SCOPE_IF_ENABLED("Collision Handling");

    const bool has_batch = collision_handler_->supports_batch();
    const bool multi_event_mode = config_.physics.collision_multi_event_mode;
    const int configured_subcycles = std::max(1, config_.physics.collision_subcycles_per_step);
    const int min_multi_event_substeps = std::max(1, config_.physics.collision_max_events_per_step);
    const int subcycles = multi_event_mode
        ? std::max(configured_subcycles, min_multi_event_substeps)
        : configured_subcycles;
    const bool needs_cpu_split = subcycles > 1;
    const bool deep_collision_enabled = deep_collision_diagnostics_.enabled();
    const size_t domain_count = config_.domains.size();
    auto per_domain = group_active_indices_by_domain(ensemble, domain_indices, domain_count);

    for (size_t dom = 0; dom < domain_count; ++dom) {
        auto& indices = per_domain[dom];
        if (indices.empty()) {
            continue;
        }
        const auto& env = config_.domains[dom].environment;
        const bool batch_safe_flow = env.flow_model == config::FlowModelKind::Constant;
        bool handled = false;
        if (has_batch && !multi_event_mode && !needs_cpu_split && !deep_collision_enabled && batch_safe_flow) {
            // Use batch path only if all dt are equal to avoid per-ion timestep bias.
            const double dt_batch = dt_used_per_ion[indices.front()];
            bool uniform_dt = true;
            for (size_t idx : indices) {
                if (dt_used_per_ion[idx] != dt_batch) {
                    uniform_dt = false;
                    break;
                }
            }
            if (uniform_dt) {
                handled = collision_handler_->handle_batch(
                    ensemble, indices, dt_batch, env, rng_by_ion_);
                if (handled) {
                    collision_runtime_stats_.add_incomplete_batch_attempts(indices.size());
                }
            }
        }
        if (!handled) {
            handle_collisions_cpu(ensemble, dt_used_per_ion, indices, env, static_cast<int>(dom));
        }
    }
}

void SimulationEngine::handle_collisions_cpu(core::IonEnsemble& ensemble,
                                             const std::vector<double>& dt_used_per_ion,
                                             const std::vector<size_t>& indices,
                                             const config::EnvironmentConfig& env,
                                             int domain_index) {
    if (!collision_handler_ || indices.empty()) {
        return;
    }

    const auto* active = ensemble.active_data();
    const auto* born = ensemble.born_data();
    const bool use_omp = parallel_enabled_;
    const bool multi_event_mode = config_.physics.collision_multi_event_mode;
    const int configured_subcycles = std::max(1, config_.physics.collision_subcycles_per_step);
    const int min_multi_event_substeps = std::max(1, config_.physics.collision_max_events_per_step);
    const int collision_substeps = multi_event_mode
        ? std::max(configured_subcycles, min_multi_event_substeps)
        : configured_subcycles;
    uint64_t macro_attempts_local = 0;
    uint64_t substep_attempts_local = 0;
    uint64_t events_local = 0;
    std::atomic<bool> collision_exception{false};
    std::string collision_exception_msg;
    const config::GeometryConfig* geometry = nullptr;
    const Mat3* rotation_global_to_local = nullptr;
    const Mat3* rotation_local_to_global = nullptr;
    if (domain_index >= 0 && static_cast<size_t>(domain_index) < config_.domains.size()) {
        const auto& domain = config_.domains[static_cast<size_t>(domain_index)];
        geometry = &domain.geometry;
        rotation_global_to_local = &domain.rotation_global_to_local;
        rotation_local_to_global = &domain.rotation_local_to_global;
    }

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk) reduction(+:macro_attempts_local,substep_attempts_local,events_local)
        for (int k = 0; k < static_cast<int>(indices.size()); ++k) {
            if (collision_exception.load(std::memory_order_relaxed)) {
                continue;
            }
            const size_t ion_idx = indices[static_cast<size_t>(k)];
            if (!active[ion_idx] || !born[ion_idx]) {
                continue;
            }
            macro_attempts_local += 1;
            const double dt_total = dt_used_per_ion[ion_idx];
            double elapsed = 0.0;
            const double collision_dt = dt_total / static_cast<double>(collision_substeps);
            for (int substep = 0; substep < collision_substeps; ++substep) {
                if (collision_exception.load(std::memory_order_relaxed)) {
                    break;
                }
                substep_attempts_local += 1;
                auto view = ensemble.collision_data(ion_idx);
                    config::EnvironmentConfig local_env = env;
                    if (env.flow_model != config::FlowModelKind::Constant) {
                        local_env.gas_velocity_m_s = env.gas_velocity_at(
                            view.kin.pos(),
                            geometry,
                            rotation_global_to_local,
                            rotation_local_to_global);
                    }
                physics::CollisionEventDiagnostics event_diag{};
                const bool collect_deep = deep_collision_diagnostics_.enabled();
                const Vec3 pos_before = collect_deep ? view.kin.pos() : Vec3{};
                const Vec3 vel_before = collect_deep ? view.kin.vel() : Vec3{};
                bool collided = false;
                try {
                    collided = collision_handler_->handle_collision(
                        view,
                        collision_dt,
                        rng_by_ion_[ion_idx],
                        local_env,
                        collect_deep ? &event_diag : nullptr);
                } catch (const std::exception& e) {
                    if (!collision_exception.exchange(true)) {
                        std::ostringstream oss;
                        oss << "Collision handling failed in domain '"
                            << config_.domains[static_cast<size_t>(domain_index)].name
                            << "' for ion index " << ion_idx
                            << " (species='" << ensemble.species_id(ion_idx) << "'): "
                            << e.what();
                        #pragma omp critical(simengine_collision_exception)
                        {
                            collision_exception_msg = oss.str();
                        }
                    }
                    break;
                } catch (...) {
                    if (!collision_exception.exchange(true)) {
                        std::ostringstream oss;
                        oss << "Collision handling failed in domain '"
                            << config_.domains[static_cast<size_t>(domain_index)].name
                            << "' for ion index " << ion_idx
                            << " (species='" << ensemble.species_id(ion_idx) << "'): unknown exception";
                        #pragma omp critical(simengine_collision_exception)
                        {
                            collision_exception_msg = oss.str();
                        }
                    }
                    break;
                }
                if (collided) {
                    events_local += 1;
                }
                if (collect_deep) {
                    const double event_time = current_time_ + elapsed + 0.5 * collision_dt;
                    deep_collision_diagnostics_.note_collision(
                        ion_idx,
                        domain_index,
                        collided,
                        event_time,
                        view.kin.get_mass(),
                        pos_before,
                        event_diag.v_rel_before_m_s,
                        event_diag.sigma_mt_m2,
                        vel_before,
                        view.kin.vel(),
                        !active[ion_idx]);
                }
                elapsed += collision_dt;
            }
        }
    }

    collision_runtime_stats_.add_cpu_counts(macro_attempts_local, substep_attempts_local, events_local);
    if (collision_exception.load()) {
        if (collision_exception_msg.empty()) {
            throw std::runtime_error("Collision handling failed due to an exception in parallel section");
        }
        throw std::runtime_error(collision_exception_msg);
    }
}

void SimulationEngine::perform_reactions(core::IonEnsemble& ensemble,
                                         const std::vector<double>& dt_used_per_ion,
                                         const std::vector<int>& domain_indices) {
    if (!reaction_handler_ || config_.reaction_db.reactions.empty()) {
        return;
    }

    PROFILE_SCOPE_IF_ENABLED("Reaction Handling");

    const auto active_indices = collect_active_domain_indices(ensemble, domain_indices);

    if (reaction_handler_->supports_batch()) {
        bool uniform_dt = !active_indices.empty();
        const double dt_batch = uniform_dt ? dt_used_per_ion[active_indices.front()] : 0.0;
        for (size_t idx : active_indices) {
            if (dt_used_per_ion[idx] != dt_batch) {
                uniform_dt = false;
                break;
            }
        }
        if (uniform_dt) {
            const bool handled = reaction_handler_->handle_batch(
                ensemble,
                domain_indices,
                dt_batch,
                config_.reaction_db,
                config_.species_db,
                config_.domains,
                rng_by_ion_);
            if (handled) {
                return;
            }
        }
    }

    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int k = 0; k < static_cast<int>(active_indices.size()); ++k) {
            const size_t i = active_indices[static_cast<size_t>(k)];
            const int dom = domain_indices[i];
            auto reaction_view = ensemble.reaction_data(i);
            reaction_handler_->handle_reaction(
                reaction_view,
                dt_used_per_ion[i],
                rng_by_ion_[i],
                config_.reaction_db,
                config_.species_db,
                config_.domains[static_cast<size_t>(dom)].environment);
        }
    }
}

} // namespace integrator
} // namespace ICARION
