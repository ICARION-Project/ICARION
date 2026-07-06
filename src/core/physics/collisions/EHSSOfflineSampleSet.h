// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace ICARION::physics {

inline constexpr const char* EHSS_OFFLINE_SAMPLE_SET_FORMAT = "ehss_offline_samples";
inline constexpr const char* EHSS_OFFLINE_SAMPLE_SET_UNITS = "sigma_eff_m2,mu";

struct EHSSOfflineSampleSet {
    std::string gas;
    std::vector<double> sigma_eff_m2;    // size N
    std::vector<double> mu_samples_flat; // size N * mu_stride
    size_t mu_stride = 0;

    bool valid() const {
        return !sigma_eff_m2.empty()
            && mu_stride > 0
            && mu_samples_flat.size() == sigma_eff_m2.size() * mu_stride;
    }

    double mu_at(size_t orientation_idx, size_t mu_idx) const {
        return mu_samples_flat[orientation_idx * mu_stride + mu_idx];
    }
};

bool load_ehss_offline_sample_set_file(
    const std::filesystem::path& path,
    EHSSOfflineSampleSet& out,
    std::string* error_msg = nullptr
);

} // namespace ICARION::physics
