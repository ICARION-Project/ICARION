// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include <array>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <vector>

namespace ICARION::physics {

struct EHSSOrientationSamples {
    std::vector<std::array<double, 4>> orientations_quat;  // [w, x, y, z]
    std::unordered_map<std::string, std::vector<double>> areas_by_gas_m2;
};

bool load_ehss_samples_file(
    const std::filesystem::path& path,
    EHSSOrientationSamples& out,
    std::string* error_msg = nullptr
);

} // namespace ICARION::physics
