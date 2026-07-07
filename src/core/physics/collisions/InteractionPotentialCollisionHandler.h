// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file InteractionPotentialCollisionHandler.h
 * @brief Interaction-potential offline collision handler (formerly InteractionPotentialModel)
 */

#pragma once

#include "ICollisionHandler.h"
#include "InteractionPotentialOfflineSampleSet.h"
#include "core/config/types/SpeciesConfig.h"
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ICARION::physics {

/**
 * @brief InteractionPotential collision handler using offline samples
 *
 * Requires per-species offline sample tables (HDF5) referenced via
 * `ipm_samples_file` in the species database.
 */
class InteractionPotentialCollisionHandler : public ICollisionHandler {
public:
    explicit InteractionPotentialCollisionHandler(bool enable_logging = false,
                                                  const config::SpeciesDatabase* species_db = nullptr,
                                                  const std::string& orientation_mode = "random",
                                                  int fixed_orientation_index = 0,
                                                  const std::string& vrel_log_prefix = "",
                                                  const std::string& momentum_log_prefix = "");
    ~InteractionPotentialCollisionHandler() override;

    bool handle_collision(
        core::IonCollisionData& view,
        double dt,
        PhysicsRng& rng,
        const config::EnvironmentConfig& env,
        CollisionEventDiagnostics* diagnostics = nullptr
    ) override;

    std::string name() const override { return "InteractionPotentialModel"; }
    CollisionStats get_stats() const override;
    void reset_stats() override;

private:
    void flush_logs();

    struct VrelLog {
        std::vector<uint64_t> bin_counts;
        std::vector<double> logv_bins;
        double v_sum = 0.0;
        double v2_sum = 0.0;
        uint64_t samples = 0;
    };

    struct MomentumLog {
        double dp_par_sum = 0.0;
        double k_sum = 0.0;
        double v_sum = 0.0;
        double v2_sum = 0.0;
        uint64_t attempts = 0;
        uint64_t collisions = 0;
        double q_c = 0.0;
    };

    bool enable_logging_ = false;
    const config::SpeciesDatabase* species_db_ = nullptr;
    std::unordered_map<std::string, InteractionPotentialOfflineSampleSet> samples_;
    std::unordered_set<std::string> logged_use_;
    std::unordered_set<std::string> warned_low_coverage_;
    std::unordered_set<std::string> warned_high_coverage_;

    CollisionStats stats_;
    size_t rate_samples_ = 0;
    double rate_sum_ = 0.0;

    bool vrel_logging_enabled_ = false;
    std::string vrel_log_prefix_;
    std::unordered_map<std::string, VrelLog> vrel_logs_;

    bool momentum_logging_enabled_ = false;
    std::string momentum_log_prefix_;
    std::unordered_map<std::string, MomentumLog> momentum_logs_;
    mutable std::mutex state_mutex_;

    enum class OrientationMode { Random, Fixed };
    OrientationMode orientation_mode_ = OrientationMode::Random;
    size_t fixed_orientation_index_ = 0;

};
} // namespace ICARION::physics
