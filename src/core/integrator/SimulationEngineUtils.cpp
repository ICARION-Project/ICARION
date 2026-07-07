// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "core/integrator/SimulationEngineUtils.h"
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace ICARION {
namespace integrator {

std::vector<size_t> collect_active_domain_indices(const core::IonEnsemble& ensemble,
                                                  const std::vector<int>& domain_indices) {
    std::vector<size_t> active_indices;
    active_indices.reserve(ensemble.size());
    const auto* active = ensemble.active_data();
    const auto* born = ensemble.born_data();
    for (size_t i = 0; i < ensemble.size(); ++i) {
        if (active[i] && born[i] && domain_indices[i] >= 0) {
            active_indices.push_back(i);
        }
    }
    return active_indices;
}

std::vector<std::vector<size_t>> group_active_indices_by_domain(const core::IonEnsemble& ensemble,
                                                                const std::vector<int>& domain_indices,
                                                                size_t domain_count) {
    std::vector<std::vector<size_t>> per_domain(domain_count);
    for (size_t idx : collect_active_domain_indices(ensemble, domain_indices)) {
        per_domain[static_cast<size_t>(domain_indices[idx])].push_back(idx);
    }
    return per_domain;
}

bool adaptive_space_charge_enabled() {
    const char* env = std::getenv("ICARION_ADAPTIVE_SC");
    if (!env) {
        return true;
    }
    std::string val(env);
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    return !(val == "0" || val == "false" || val == "off");
}

double next_engine_time_after_step(const core::IonEnsemble& ensemble,
                                   double current_time,
                                   double total_time_s,
                                   double fallback_dt_s,
                                   double max_dt_used) {
    double new_time = current_time;
    const double* time_ptr = ensemble.time_data();
    for (size_t i = 0; i < ensemble.size(); ++i) {
        new_time = std::max(new_time, time_ptr[i]);
    }

    if (new_time > current_time) {
        return new_time;
    }

    double next_birth_time = std::numeric_limits<double>::max();
    const auto* born_flags = ensemble.born_data();
    for (size_t i = 0; i < ensemble.size(); ++i) {
        if (born_flags[i]) {
            continue;
        }
        const double bt = ensemble.birth_time(i);
        if (bt > current_time && bt <= total_time_s) {
            next_birth_time = std::min(next_birth_time, bt);
        }
    }

    if (next_birth_time < std::numeric_limits<double>::max()) {
        return next_birth_time;
    }
    if (max_dt_used > 0.0) {
        return current_time + max_dt_used;
    }
    return current_time + fallback_dt_s;
}

void CollisionRuntimeStats::reset() {
    macro_attempts_total = 0;
    substep_attempts_total = 0;
    events_total = 0;
    monitor_complete = true;
}

void CollisionRuntimeStats::mark_batch_path_used() {
    monitor_complete = false;
}

void CollisionRuntimeStats::add_batch_attempts(size_t ion_count) {
    macro_attempts_total += static_cast<uint64_t>(ion_count);
    substep_attempts_total += static_cast<uint64_t>(ion_count);
}

void CollisionRuntimeStats::add_incomplete_batch_attempts(size_t ion_count) {
    mark_batch_path_used();
    add_batch_attempts(ion_count);
}

void CollisionRuntimeStats::add_cpu_counts(uint64_t macro_attempts,
                                           uint64_t substep_attempts,
                                           uint64_t events) {
    macro_attempts_total += macro_attempts;
    substep_attempts_total += substep_attempts;
    events_total += events;
}

bool CollisionRuntimeStats::has_activity() const {
    return macro_attempts_total > 0 || substep_attempts_total > 0;
}

bool CollisionRuntimeStats::should_warn_single_collision_timestep_load() const {
    return monitor_complete && event_fraction_per_ion_step() > 0.10;
}

std::string CollisionRuntimeStats::summary_message(int completed_steps) const {
    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "Collision load: avg collisions/step=" << mean_events_per_step(completed_steps)
        << ", event fraction per ion-step=" << (100.0 * event_fraction_per_ion_step()) << "%"
        << ", event fraction per substep=" << (100.0 * event_fraction_per_substep()) << "%"
        << " (events=" << events_total
        << ", ion_steps=" << macro_attempts_total
        << ", substeps=" << substep_attempts_total << ")";
    if (!monitor_complete) {
        msg << " [monitor incomplete: batch collision path used]";
    }
    return msg.str();
}

std::string CollisionRuntimeStats::single_collision_timestep_warning() const {
    std::ostringstream msg;
    msg << std::fixed << std::setprecision(2)
        << "WARNING: collision event fraction per ion-step is "
        << (100.0 * event_fraction_per_ion_step())
        << "% (>10%). Consider reducing dt or enabling collision micro-subcycling.";
    return msg.str();
}

double CollisionRuntimeStats::mean_events_per_step(int completed_steps) const {
    const int denom = (completed_steps > 0) ? completed_steps : 1;
    return static_cast<double>(events_total) / static_cast<double>(denom);
}

double CollisionRuntimeStats::event_fraction_per_ion_step() const {
    const uint64_t denom = (macro_attempts_total > 0) ? macro_attempts_total : 1;
    return static_cast<double>(events_total) / static_cast<double>(denom);
}

double CollisionRuntimeStats::event_fraction_per_substep() const {
    const uint64_t denom = (substep_attempts_total > 0) ? substep_attempts_total : 1;
    return static_cast<double>(events_total) / static_cast<double>(denom);
}

}  // namespace integrator
}  // namespace ICARION
