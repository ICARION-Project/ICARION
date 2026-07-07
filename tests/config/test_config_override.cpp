// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer

#include "core/config/types/FullConfig.h"
#include "core/config/utils/ConfigOverride.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using namespace ICARION::config;

TEST_CASE("ConfigOverride applies collision subcycle and IPM controls", "[config][override][physics]") {
    FullConfig cfg;

    ConfigOverride::apply(cfg, {
        {"physics.collision_multi_event_mode", "true"},
        {"physics.collision_max_events_per_step", "8"},
        {"physics.collision_time_centered", "true"},
        {"physics.collision_time_randomized", "false"},
        {"physics.collision_subcycles_per_step", "4"},
        {"physics.ipm_orientation_mode", "fixed"},
        {"physics.ipm_fixed_orientation_index", "3"},
        {"physics.ipm_vrel_log_prefix", "vrel"},
        {"physics.ipm_momentum_log_prefix", "momentum"}
    });

    REQUIRE(cfg.physics.collision_multi_event_mode);
    REQUIRE(cfg.physics.collision_max_events_per_step == 8);
    REQUIRE(cfg.physics.collision_time_centered);
    REQUIRE_FALSE(cfg.physics.collision_time_randomized);
    REQUIRE(cfg.physics.collision_subcycles_per_step == 4);
    REQUIRE(cfg.physics.ipm_orientation_mode_type == IPMOrientationMode::Fixed);
    REQUIRE(cfg.physics.ipm_orientation_mode == "fixed");
    REQUIRE(cfg.physics.ipm_fixed_orientation_index == 3);
    REQUIRE(cfg.physics.ipm_vrel_log_prefix == "vrel");
    REQUIRE(cfg.physics.ipm_momentum_log_prefix == "momentum");
}

TEST_CASE("ConfigOverride rejects unknown IPM orientation mode", "[config][override][physics]") {
    FullConfig cfg;

    REQUIRE_THROWS_WITH(
        ConfigOverride::apply(cfg, {{"physics.ipm_orientation_mode", "uniform"}}),
        Catch::Matchers::ContainsSubstring("Unknown ipm_orientation_mode"));
}

TEST_CASE("ConfigOverride applies deep collision analysis mode", "[config][override][output]") {
    FullConfig cfg;

    ConfigOverride::apply(cfg, {{"output.deep_analysis.mode", "sampled_events"}});

    REQUIRE(cfg.output.deep_analysis.mode_type == DeepAnalysisMode::SampledEvents);
    REQUIRE(cfg.output.deep_analysis.mode == "sampled_events");
}

TEST_CASE("ConfigOverride rejects unknown deep collision analysis mode", "[config][override][output]") {
    FullConfig cfg;

    REQUIRE_THROWS_WITH(
        ConfigOverride::apply(cfg, {{"output.deep_analysis.mode", "verbose"}}),
        Catch::Matchers::ContainsSubstring("Unknown deep_analysis.mode"));
}
