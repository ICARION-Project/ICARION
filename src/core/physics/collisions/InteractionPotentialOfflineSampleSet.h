// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace ICARION::physics {

inline constexpr const char* INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_FORMAT = "ipm_offline_samples";
inline constexpr const char* INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_UNITS = "logv, sigma_mt_m2, sigma_event_m2, b_max_m, dp_SI";
inline constexpr int INTERACTION_POTENTIAL_OFFLINE_SAMPLE_SET_VERSION = 1;

struct InteractionPotentialOfflineSampleSet {
    std::string gas;
    std::vector<double> logv_bins;    // size K
    std::vector<double> sigma_mt_m2;  // size N*K
    std::vector<double> sigma_event_m2; // size N*K
    std::vector<double> q1_m2;        // size N*K (optional; if absent, sigma_mt_m2 is used as Q1)
    std::vector<double> q2_m2;        // size N*K (optional)
    std::vector<double> q3_m2;        // size N*K (optional)
    std::vector<double> q12_m2;       // size N*K (optional; event-moment numerator <h1*h2>)
    std::vector<double> q13_m2;       // size N*K (optional; event-moment numerator <h1*h3>)
    std::vector<double> b_max_m;      // size N*K
    std::vector<long long> cdf_offsets; // size N*K
    std::vector<long long> cdf_counts;  // size N*K
    std::vector<double> cdf_values;     // size total_samples
    std::vector<double> dp_samples;     // size total_samples * 3
    std::vector<double> dp_stats;       // size N*K*4
    size_t n_orient = 0;
    size_t n_bins = 0;

    bool valid() const {
        if (n_orient == 0 || n_bins == 0) return false;
        if (logv_bins.size() != n_bins) return false;
        if (sigma_mt_m2.size() != n_orient * n_bins) return false;
        if (sigma_event_m2.size() != n_orient * n_bins) return false;
        if (!q1_m2.empty() && q1_m2.size() != n_orient * n_bins) return false;
        if (!q2_m2.empty() && q2_m2.size() != n_orient * n_bins) return false;
        if (!q3_m2.empty() && q3_m2.size() != n_orient * n_bins) return false;
        if (!q12_m2.empty() && q12_m2.size() != n_orient * n_bins) return false;
        if (!q13_m2.empty() && q13_m2.size() != n_orient * n_bins) return false;
        if (b_max_m.size() != n_orient * n_bins) return false;
        if (!cdf_offsets.empty() && cdf_offsets.size() != n_orient * n_bins) return false;
        if (!cdf_counts.empty() && cdf_counts.size() != n_orient * n_bins) return false;
        if (!dp_stats.empty() && dp_stats.size() != n_orient * n_bins * 4) return false;
        if (cdf_offsets.empty() != cdf_counts.empty()) return false;
        if (!dp_samples.empty() && (dp_samples.size() % 3 != 0)) return false;

        if (!cdf_offsets.empty()) {
            const size_t n = n_orient * n_bins;
            const size_t cdf_size = cdf_values.size();
            for (size_t i = 0; i < n; ++i) {
                const long long offset = cdf_offsets[i];
                const long long count = cdf_counts[i];
                if (count < 0) return false;

                // Backward compatibility: legacy precomputes use offset=-1 as
                // sentinel for empty bins (count=0, no CDF samples).
                if (count == 0) {
                    if (offset < -1) return false;
                    if (offset >= 0) {
                        const size_t off = static_cast<size_t>(offset);
                        if (off > cdf_size) return false;
                    }
                    continue;
                }

                if (offset < 0) return false;
                const size_t off = static_cast<size_t>(offset);
                const size_t cnt = static_cast<size_t>(count);
                if (off > cdf_size || cnt > (cdf_size - off)) return false;
                if (cnt > 0) {
                    double prev = cdf_values[off];
                    if (prev < 0.0 || prev > 1.0) return false;
                    for (size_t j = 1; j < cnt; ++j) {
                        const double v = cdf_values[off + j];
                        if (v < 0.0 || v > 1.0) return false;
                        if (v + 1e-12 < prev) return false;  // must be non-decreasing for lower_bound
                        prev = v;
                    }
                }
            }

            if (!cdf_values.empty()) {
                const size_t n_samples = dp_samples.size() / 3;
                if (n_samples < cdf_values.size()) return false;
            }
        }
        return true;
    }

    size_t index(size_t orient_idx, size_t bin_idx) const {
        return orient_idx * n_bins + bin_idx;
    }
};

bool load_interaction_potential_offline_sample_set_file(
    const std::filesystem::path& path,
    InteractionPotentialOfflineSampleSet& out,
    std::string* error_msg = nullptr
);

} // namespace ICARION::physics
