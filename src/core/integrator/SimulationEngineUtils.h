// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/types/IonEnsemble.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace ICARION {
namespace integrator {

std::vector<size_t> collect_active_domain_indices(const core::IonEnsemble& ensemble,
                                                  const std::vector<int>& domain_indices);

std::vector<std::vector<size_t>> group_active_indices_by_domain(const core::IonEnsemble& ensemble,
                                                                const std::vector<int>& domain_indices,
                                                                size_t domain_count);

bool adaptive_space_charge_enabled();

double next_engine_time_after_step(const core::IonEnsemble& ensemble,
                                   double current_time,
                                   double total_time_s,
                                   double fallback_dt_s,
                                   double max_dt_used);

struct CollisionRuntimeStats {
    uint64_t macro_attempts_total = 0;
    uint64_t substep_attempts_total = 0;
    uint64_t events_total = 0;
    bool monitor_complete = true;

    void reset();
    void mark_batch_path_used();
    void add_batch_attempts(size_t ion_count);
    void add_incomplete_batch_attempts(size_t ion_count);
    void add_cpu_counts(uint64_t macro_attempts, uint64_t substep_attempts, uint64_t events);

    bool has_activity() const;
    bool should_warn_single_collision_timestep_load() const;
    std::string summary_message(int completed_steps) const;
    std::string single_collision_timestep_warning() const;

    double mean_events_per_step(int completed_steps) const;
    double event_fraction_per_ion_step() const;
    double event_fraction_per_substep() const;
};

}  // namespace integrator
}  // namespace ICARION
