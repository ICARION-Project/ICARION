// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationEngine.h"

#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/utils/Profiler.h"

#include <algorithm>

namespace ICARION {
namespace integrator {

double SimulationEngine::perform_integration(core::IonEnsemble& ensemble,
                                             double t,
                                             const std::vector<double>& dt_per_ion,
                                             const std::vector<int>& domain_indices,
                                             std::vector<double>& dt_used_per_ion,
                                             std::vector<double>& dt_next_per_ion) {
    double max_dt_used = 0.0;

    // Only run batch if dt is uniform across active ions.
    bool uniform_dt = false;
    double dt_batch = 0.0;
    if (integrator_->supports_batch()) {
        uniform_dt = true;
        for (size_t i = 0; i < domain_indices.size(); ++i) {
            if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
            if (dt_batch == 0.0) {
                dt_batch = dt_per_ion[i];
            } else if (dt_per_ion[i] != dt_batch) {
                uniform_dt = false;
                break;
            }
        }
        if (uniform_dt && dt_batch > 0.0) {
            PROFILE_SCOPE_IF_ENABLED("Integrator Batch Step");
            if (integrator_->step_batch(ensemble, t, dt_batch, force_registries_, domain_indices)) {
                std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                return dt_batch;
            }
        }
    }

    if (has_space_charge_model()) {
        auto refresh_space_charge_stage = [&](double stage_time_s) {
            const double saved_time = current_time_;
            current_time_ = stage_time_s;
            update_space_charge_models(ensemble);
            current_time_ = saved_time;
        };

        if (auto* rk4 = dynamic_cast<RK4Strategy*>(integrator_.get())) {
            if (uniform_dt && dt_batch > 0.0) {
                PROFILE_SCOPE_IF_ENABLED("Integrator RK4 SC Batch");
                if (rk4->step_batch_with_stage_refresh(
                        ensemble, t, dt_batch, force_registries_, domain_indices,
                        refresh_space_charge_stage)) {
                    std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                    std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                    dt_per_ion_.assign(ensemble.size(), dt_batch);
                    return dt_batch;
                }
            }
            output_manager_->log_progress(
                "Warning: Space-charge present without uniform RK4 dt; falling back to per-ion integration with stale fields.");
        } else if (auto* rk45 = dynamic_cast<RK45Strategy*>(integrator_.get())) {
            double dt_seed = dt_batch;
            if (dt_seed <= 0.0) {
                for (size_t i = 0; i < ensemble.size(); ++i) {
                    if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                    dt_seed = dt_per_ion[i];
                    break;
                }
            }
            if (dt_seed <= 0.0) {
                dt_seed = config_.simulation.dt_s;
            }

            if (integrator_->is_adaptive()) {
                PROFILE_SCOPE_IF_ENABLED("Integrator RK45 SC Adaptive");
                const auto result = rk45->step_batch_adaptive(
                    ensemble, t, dt_seed, force_registries_, domain_indices,
                    refresh_space_charge_stage);
                if (result.accepted) {
                    for (size_t i = 0; i < ensemble.size(); ++i) {
                        if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                        dt_used_per_ion[i] = result.dt_used;
                        dt_next_per_ion[i] = result.dt_next;
                    }
                    dt_per_ion_.assign(ensemble.size(), result.dt_next);
                    return result.dt_used;
                }
            } else if (uniform_dt && dt_batch > 0.0) {
                PROFILE_SCOPE_IF_ENABLED("Integrator RK45 SC Batch");
                if (rk45->step_batch_fixed(
                        ensemble, t, dt_batch, force_registries_, domain_indices,
                        refresh_space_charge_stage)) {
                    std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                    std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                    dt_per_ion_.assign(ensemble.size(), dt_batch);
                    return dt_batch;
                }
            }
        }
    }

    auto* rk45 = integrator_->is_adaptive()
        ? dynamic_cast<RK45Strategy*>(integrator_.get())
        : nullptr;

    const size_t n = ensemble.size();
    for (size_t i = 0; i < n; ++i) {
        if (domain_indices[i] < 0 || !ensemble.is_active(i)) {
            continue;
        }
        const int dom = domain_indices[i];
        if (dom < 0 || dom >= static_cast<int>(force_registries_.size())) {
            continue;
        }
        const auto& registry = force_registries_[static_cast<size_t>(dom)];
        if (!registry) {
            continue;
        }
        const double dt = dt_per_ion[i];
        PROFILE_SCOPE_IF_ENABLED("Integrator Step");
        integrator_->step(ensemble, i, t, dt, *registry);

        if (rk45) {
            const double used = rk45->last_dt_used();
            const double next = rk45->last_dt_suggested();
            if (used > 0.0) {
                dt_used_per_ion[i] = used;
                max_dt_used = std::max(max_dt_used, used);
            }
            if (next > 0.0) {
                dt_next_per_ion[i] = next;
            }
        } else {
            dt_used_per_ion[i] = dt;
            dt_next_per_ion[i] = dt;
            max_dt_used = std::max(max_dt_used, dt);
        }
    }
    return max_dt_used;
}

} // namespace integrator
} // namespace ICARION
