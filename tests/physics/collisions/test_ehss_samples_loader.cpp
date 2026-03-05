// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>

#include "core/physics/collisions/EHSSSamples.h"

using ICARION::physics::EHSSOrientationSamples;
using ICARION::physics::load_ehss_samples_file;

TEST_CASE("EHSS samples file loads", "[collision][ehss][samples]") {
    auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
    std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("ehss_samples_" + std::to_string(stamp) + ".json");

    std::ofstream ofs(path);
    REQUIRE(ofs.is_open());
    ofs << "{\n"
        << "  \"version\": 1,\n"
        << "  \"species_id\": \"X+\",\n"
        << "  \"n_orientations\": 2,\n"
        << "  \"orientations_quat\": [[1,0,0,0],[0,1,0,0]],\n"
        << "  \"areas_by_gas_m2\": {\n"
        << "    \"He\": [1.0e-18, 2.0e-18]\n"
        << "  }\n"
        << "}\n";
    ofs.close();

    EHSSOrientationSamples samples;
    std::string error;
    REQUIRE(load_ehss_samples_file(path, samples, &error));
    REQUIRE(samples.orientations_quat.size() == 2);
    REQUIRE(samples.areas_by_gas_m2.count("He") == 1);
    REQUIRE(samples.areas_by_gas_m2.at("He").size() == 2);

    std::filesystem::remove(path);
}
