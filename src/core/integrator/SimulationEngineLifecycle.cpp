// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationEngine.h"

#include "core/utils/Profiler.h"
#include "core/utils/safety/numericalSafetyGuards.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>

namespace ICARION {
namespace integrator {

void SimulationEngine::update_space_charge_models(core::IonEnsemble& ensemble) {
    PROFILE_SCOPE_IF_ENABLED("Space Charge Update");
    std::vector<const physics::ISpaceChargeModel*> updated;
    updated.reserve(force_registries_.size());

    for (const auto& registry : force_registries_) {
        if (!registry) continue;
        auto* model = registry->space_charge_model();
        if (!model) continue;

        if (std::find(updated.begin(), updated.end(), model) != updated.end()) {
            continue;  // already updated (shared model)
        }
        model->update_fields(ensemble, current_time_);
        updated.push_back(model);
    }
}

void SimulationEngine::update_dynamic_environments(double t) {
    for (auto& domain : config_.domains) {
        if (!domain.environment.has_dynamic_pressure()) {
            continue;
        }

        try {
            domain.environment.update_time_dependent(t, config_.waveforms);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "Failed to evaluate pressure waveform for domain '" + domain.name + "': " + e.what());
        }
    }
}

void SimulationEngine::refresh_ensemble_environment_cache(core::IonEnsemble& ensemble) {
    const size_t n_ions = ensemble.size();
    for (size_t i = 0; i < n_ions; ++i) {
        const int dom = ensemble.domain_index(i);
        if (dom < 0 || static_cast<size_t>(dom) >= config_.domains.size()) {
            continue;
        }
        const auto& env = config_.domains[static_cast<size_t>(dom)].environment;
        ensemble.temperature_data()[i] = env.temperature_K;
        ensemble.gas_density_data()[i] = env.particle_density_m_3;
        ensemble.neutral_mass_data()[i] = env.gas_mass_kg;
    }
}

bool SimulationEngine::has_space_charge_model() const {
    for (const auto& registry : force_registries_) {
        if (registry && registry->space_charge_model()) {
            return true;
        }
    }
    return false;
}

std::vector<int> SimulationEngine::prepare_ions_for_integration(core::IonEnsemble& ensemble) {
    const size_t n_ions = ensemble.size();
    std::vector<int> integration_domains(n_ions, -1);

    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* active = ensemble.active_data();
    auto* born = ensemble.born_data();
    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n_ions); ++i) {
            if (!born[i] && current_time_ >= ensemble.birth_time(i)) {
                born[i] = 1;
                active[i] = 1;
            }

            if (!active[i] || !born[i]) {
                continue;
            }

            const Vec3 pos(pos_x[i], pos_y[i], pos_z[i]);
            int domain_idx = -1;
            {
                PROFILE_SCOPE_IF_ENABLED("Domain Finding");
                domain_idx = domain_manager_->find_domain_index(pos);
            }
            if (domain_idx < 0) {
                active[i] = 0;
                ensemble.set_death_time(i, current_time_);
                continue;
            }

            integration_domains[static_cast<size_t>(i)] = domain_idx;
            const auto& domain_config = config_.domains[static_cast<size_t>(domain_idx)];

            if (ensemble.domain_index(i) != domain_idx) {
                ensemble.update_domain_cache(i, domain_idx,
                    domain_config.environment.temperature_K,
                    domain_config.environment.particle_density_m_3,
                    domain_config.environment.gas_mass_kg);
            }
        }
    }

    return integration_domains;
}

void SimulationEngine::finalize_ions_after_integration(core::IonEnsemble& ensemble,
                                                       const std::vector<double>& dt_used_per_ion,
                                                       std::vector<double>& dt_next_per_ion) {
    const size_t n_ions = ensemble.size();
    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* active = ensemble.active_data();
    auto* born = ensemble.born_data();
    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n_ions); ++i) {
            if (!active[i] || !born[i]) {
                continue;
            }

            const Vec3 pos_after(pos_x[i], pos_y[i], pos_z[i]);
            int new_domain_idx;
            {
                PROFILE_SCOPE_IF_ENABLED("Boundary Checks");
                new_domain_idx = domain_manager_->find_domain_index(pos_after);
            }
            if (new_domain_idx < 0) {
                const int prev_dom = ensemble.domain_index(i);
                new_domain_idx = domain_manager_->forward_axial_bridge_domain(prev_dom, pos_after);
                if (new_domain_idx < 0) {
                    active[i] = 0;
                    ensemble.set_death_time(i, current_time_);
                    continue;
                }
            } else {
                const int cur_dom_idx = ensemble.domain_index(i);
                const int handoff_domain = domain_manager_->shared_boundary_handoff_domain(
                    cur_dom_idx, new_domain_idx, pos_after);
                if (handoff_domain >= 0) {
                    new_domain_idx = handoff_domain;
                }
            }
            if (new_domain_idx != ensemble.domain_index(i)) {
                const auto& new_dom = config_.domains[static_cast<size_t>(new_domain_idx)];
                ensemble.update_domain_cache(i, new_domain_idx,
                    new_dom.environment.temperature_K,
                    new_dom.environment.particle_density_m_3,
                    new_dom.environment.gas_mass_kg);
            }

            ensemble.set_time(i, ensemble.time(i) + dt_used_per_ion[static_cast<size_t>(i)]);
            if (dt_next_per_ion[static_cast<size_t>(i)] <= 0.0) {
                dt_next_per_ion[static_cast<size_t>(i)] =
                    dt_used_per_ion[static_cast<size_t>(i)] > 0.0
                        ? dt_used_per_ion[static_cast<size_t>(i)]
                        : config_.simulation.dt_s;
            }

            const Vec3 vel_check = ensemble.get_vel(i);
            const bool position_valid = ICARION::safety::is_finite(pos_after);
            const bool velocity_valid = ICARION::safety::is_finite(vel_check);

            if (!position_valid || !velocity_valid) {
                active[i] = 0;
                ensemble.set_death_time(i, current_time_);

                if (!config_.simulation.enable_safety_logging) {
                    std::cerr << "Warning: Ion " << i << " has invalid state at t = "
                              << ensemble.time(i) << " s, deactivating" << std::endl;
                }
            }
        }
    }
}

} // namespace integrator
} // namespace ICARION
