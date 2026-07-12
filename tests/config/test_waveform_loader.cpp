// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/config/loader/WaveformLoader.h"
#include <json/json.h>

using namespace ICARION::config;
using Catch::Matchers::WithinAbs;

// Helper to parse JSON string
Json::Value parse_json(const std::string& json_str) {
    Json::CharReaderBuilder builder;
    Json::CharReader* reader = builder.newCharReader();
    Json::Value root;
    std::string errors;
    
    bool success = reader->parse(
        json_str.c_str(),
        json_str.c_str() + json_str.size(),
        &root,
        &errors
    );
    
    delete reader;
    
    if (!success) {
        throw std::runtime_error("JSON parse error: " + errors);
    }
    
    return root;
}

TEST_CASE("WaveformLoader: load constant waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "constant",
        "value": 123.45
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    REQUIRE(std::holds_alternative<ConstantWaveform>(w.data));
    const auto& cw = std::get<ConstantWaveform>(w.data);
    CHECK(cw.value == 123.45);
}

TEST_CASE("WaveformLoader: load linear waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "linear",
        "start": 0.0,
        "end": 500.0,
        "end_time_s": 0.501,
        "start_time_s": 0.5,
        "clamp": false
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    REQUIRE(std::holds_alternative<LinearWaveform>(w.data));
    const auto& lw = std::get<LinearWaveform>(w.data);
    CHECK(lw.start_value == 0.0);
    CHECK(lw.end_value == 500.0);
    CHECK(lw.end_time_s == 0.501);
    CHECK(lw.start_time_s == 0.5);
    CHECK(lw.clamp == false);
}

TEST_CASE("WaveformLoader: linear with defaults", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "linear",
        "start": 10,
        "end": 20,
        "end_time_s": 1.0
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    const auto& lw = std::get<LinearWaveform>(w.data);
    CHECK(lw.start_time_s == 0.0);  // Default
    CHECK(lw.clamp == true);        // Default
}

TEST_CASE("WaveformLoader: load quadratic waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "quadratic",
        "a": 10.0,
        "b": 5.0,
        "c": 2.0,
        "start_time_s": 1.0,
        "end_time_s": 2.0
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    REQUIRE(std::holds_alternative<QuadraticWaveform>(w.data));
    const auto& qw = std::get<QuadraticWaveform>(w.data);
    CHECK(qw.a == 10.0);
    CHECK(qw.b == 5.0);
    CHECK(qw.c == 2.0);
}

TEST_CASE("WaveformLoader: load exponential waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "exponential",
        "offset": 1.0,
        "amplitude": 4.0,
        "rate_per_s": -10.0,
        "start_time_s": 0.0,
        "end_time_s": 0.5,
        "clamp": true
    })";

    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);

    REQUIRE(std::holds_alternative<ExponentialWaveform>(w.data));
    const auto& ew = std::get<ExponentialWaveform>(w.data);
    CHECK(ew.offset == 1.0);
    CHECK(ew.amplitude == 4.0);
    CHECK(ew.rate_per_s == -10.0);
    CHECK(ew.start_time_s == 0.0);
    CHECK(ew.end_time_s == 0.5);
    CHECK(ew.clamp == true);
}

TEST_CASE("WaveformLoader: load sinusoidal waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "sinusoidal",
        "offset": 250.0,
        "amplitude": 50.0,
        "frequency_Hz": 100.0,
        "phase_rad": 1.57
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    REQUIRE(std::holds_alternative<SinusoidalWaveform>(w.data));
    const auto& sw = std::get<SinusoidalWaveform>(w.data);
    CHECK(sw.offset == 250.0);
    CHECK(sw.amplitude == 50.0);
    CHECK(sw.frequency_Hz == 100.0);
    CHECK_THAT(sw.phase_rad, WithinAbs(1.57, 1e-9));
}

TEST_CASE("WaveformLoader: load PWM waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "pwm",
        "low": 0.0,
        "high": 5.0,
        "frequency_Hz": 1000.0,
        "duty_cycle": 0.1,
        "phase_rad": 1.57,
        "start_time_s": 0.0,
        "end_time_s": 1.0,
        "clamp": true
    })";

    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);

    REQUIRE(std::holds_alternative<PWMWaveform>(w.data));
    const auto& pw = std::get<PWMWaveform>(w.data);
    CHECK(pw.low_value == 0.0);
    CHECK(pw.high_value == 5.0);
    CHECK(pw.frequency_Hz == 1000.0);
    CHECK(pw.duty_cycle == 0.1);
    CHECK_THAT(pw.phase_rad, WithinAbs(1.57, 1e-9));
    CHECK(pw.start_time_s == 0.0);
    CHECK(pw.end_time_s == 1.0);
    CHECK(pw.clamp == true);
}

TEST_CASE("WaveformLoader: load pulsed waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "pulsed",
        "low": 10.0,
        "high": 100.0,
        "pulse_start_s": 1.0,
        "pulse_width_s": 0.5
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    REQUIRE(std::holds_alternative<PulsedWaveform>(w.data));
    const auto& pw = std::get<PulsedWaveform>(w.data);
    CHECK(pw.low_value == 10.0);
    CHECK(pw.high_value == 100.0);
    CHECK(pw.pulse_start_s == 1.0);
    CHECK(pw.pulse_width_s == 0.5);
}

TEST_CASE("WaveformLoader: load arbitrary waveform", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "arbitrary",
        "times": [0.0, 1.0, 2.0, 3.0],
        "values": [0.0, 100.0, 50.0, 200.0],
        "interpolation": "linear"
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    REQUIRE(std::holds_alternative<ArbitraryWaveform>(w.data));
    const auto& aw = std::get<ArbitraryWaveform>(w.data);
    REQUIRE(aw.times_s.size() == 4);
    REQUIRE(aw.values.size() == 4);
    CHECK(aw.times_s[0] == 0.0);
    CHECK(aw.values[2] == 50.0);
    CHECK(aw.interp == ArbitraryWaveform::Interpolation::Linear);
}

TEST_CASE("WaveformLoader: arbitrary with step interpolation", "[waveform][loader]") {
    std::string json_str = R"({
        "type": "arbitrary",
        "times": [0, 1, 2],
        "values": [10, 20, 30],
        "interpolation": "step"
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    const auto& aw = std::get<ArbitraryWaveform>(w.data);
    CHECK(aw.interp == ArbitraryWaveform::Interpolation::Step);
}

TEST_CASE("WaveformLoader: load waveform library", "[waveform][loader]") {
    std::string json_str = R"({
        "ramp1": {
            "type": "linear",
            "start": 0,
            "end": 100,
            "end_time_s": 0.1
        },
        "sine1": {
            "type": "sinusoidal",
            "amplitude": 50,
            "frequency_Hz": 100
        }
    })";
    
    Json::Value json = parse_json(json_str);
    auto library = WaveformLoader::load_library(json);
    
    REQUIRE(library.size() == 2);
    REQUIRE(library.count("ramp1") == 1);
    REQUIRE(library.count("sine1") == 1);
    
    CHECK(library["ramp1"].id == "ramp1");
    CHECK(library["sine1"].id == "sine1");
    
    // Verify types
    CHECK(std::holds_alternative<LinearWaveform>(library["ramp1"].data));
    CHECK(std::holds_alternative<SinusoidalWaveform>(library["sine1"].data));
}

TEST_CASE("ValueOrWaveform: load static value", "[waveform][loader][value_or_waveform]") {
    Json::Value json = 123.45;
    std::map<std::string, Waveform> library;
    
    ValueOrWaveform val = WaveformLoader::load_value_or_waveform(json, library);
    
    REQUIRE(val.constant_value.has_value());
    CHECK(*val.constant_value == 123.45);
    CHECK_FALSE(val.is_time_varying());
}

TEST_CASE("ValueOrWaveform: load inline waveform", "[waveform][loader][value_or_waveform]") {
    std::string json_str = R"({
        "type": "linear",
        "start": 0,
        "end": 100,
        "end_time_s": 1.0
    })";
    
    Json::Value json = parse_json(json_str);
    std::map<std::string, Waveform> library;
    
    ValueOrWaveform val = WaveformLoader::load_value_or_waveform(json, library);
    
    REQUIRE(val.waveform.has_value());
    CHECK(val.is_time_varying());
    CHECK(std::holds_alternative<LinearWaveform>(val.waveform->data));
}

TEST_CASE("ValueOrWaveform: load waveform reference", "[waveform][loader][value_or_waveform]") {
    // Create library first
    std::map<std::string, Waveform> library;
    Waveform w;
    w.id = "my_ramp";
    w.data = LinearWaveform{0, 100, 0, 1.0, true};
    library["my_ramp"] = w;
    
    // Load reference
    Json::Value json = "@my_ramp";
    ValueOrWaveform val = WaveformLoader::load_value_or_waveform(json, library);
    
    REQUIRE(val.waveform.has_value());
    CHECK(val.waveform->id == "my_ramp");
    CHECK(std::holds_alternative<LinearWaveform>(val.waveform->data));
    CHECK(val.is_time_varying());
}

TEST_CASE("WaveformLoader: missing required field throws", "[waveform][loader][error]") {
    std::string json_str = R"({
        "type": "linear",
        "start": 0
    })";
    
    Json::Value json = parse_json(json_str);
    CHECK_THROWS_AS(WaveformLoader::load(json), std::runtime_error);
}

TEST_CASE("WaveformLoader: unknown type throws", "[waveform][loader][error]") {
    std::string json_str = R"({
        "type": "nonexistent_waveform",
        "value": 100
    })";
    
    Json::Value json = parse_json(json_str);
    CHECK_THROWS_AS(WaveformLoader::load(json), std::runtime_error);
}

TEST_CASE("WaveformLoader: invalid time range throws", "[waveform][loader][error]") {
    std::string json_str = R"({
        "type": "linear",
        "start": 0,
        "end": 100,
        "start_time_s": 1.0,
        "end_time_s": 0.5
    })";
    
    Json::Value json = parse_json(json_str);
    CHECK_THROWS_AS(WaveformLoader::load(json), std::runtime_error);
}

TEST_CASE("WaveformLoader: arbitrary with mismatched arrays throws", "[waveform][loader][error]") {
    std::string json_str = R"({
        "type": "arbitrary",
        "times": [0, 1, 2],
        "values": [10, 20]
    })";
    
    Json::Value json = parse_json(json_str);
    CHECK_THROWS_AS(WaveformLoader::load(json), std::runtime_error);
}

TEST_CASE("WaveformLoader: arbitrary with unsorted times throws", "[waveform][loader][error]") {
    std::string json_str = R"({
        "type": "arbitrary",
        "times": [0, 2, 1],
        "values": [10, 20, 30]
    })";
    
    Json::Value json = parse_json(json_str);
    CHECK_THROWS_AS(WaveformLoader::load(json), std::runtime_error);
}

TEST_CASE("ValueOrWaveform: invalid waveform reference throws", "[waveform][loader][error]") {
    std::map<std::string, Waveform> library;  // Empty library
    
    Json::Value json = "@nonexistent";
    CHECK_THROWS_AS(WaveformLoader::load_value_or_waveform(json, library), std::runtime_error);
}

TEST_CASE("ValueOrWaveform: reference without @ throws", "[waveform][loader][error]") {
    std::map<std::string, Waveform> library;
    
    Json::Value json = "my_waveform";  // Missing @ prefix
    CHECK_THROWS_AS(WaveformLoader::load_value_or_waveform(json, library), std::runtime_error);
}

TEST_CASE("WaveformLoader: loaded waveforms evaluate correctly", "[waveform][loader][integration]") {
    // Load linear waveform
    std::string json_str = R"({
        "type": "linear",
        "start": 0,
        "end": 100,
        "end_time_s": 1.0
    })";
    
    Json::Value json = parse_json(json_str);
    Waveform w = WaveformLoader::load(json);
    
    // Evaluate
    CHECK_THAT(w.evaluate(0.0), WithinAbs(0.0, 1e-9));
    CHECK_THAT(w.evaluate(0.5), WithinAbs(50.0, 1e-6));
    CHECK_THAT(w.evaluate(1.0), WithinAbs(100.0, 1e-9));
}
