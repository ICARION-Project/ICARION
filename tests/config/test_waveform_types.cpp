// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/config/types/WaveformConfig.h"

using namespace ICARION::config;
using Catch::Matchers::WithinAbs;

TEST_CASE("ConstantWaveform evaluates to constant", "[waveform][constant]") {
    ConstantWaveform w{100.0};
    
    CHECK(w.evaluate(0.0) == 100.0);
    CHECK(w.evaluate(1.0) == 100.0);
    CHECK(w.evaluate(-5.0) == 100.0);
    CHECK(w.evaluate(1e6) == 100.0);
}

TEST_CASE("LinearWaveform: before start time", "[waveform][linear]") {
    LinearWaveform w{0.0, 100.0, 1.0, 1.1, true};  // start=0, end=100, t0=1.0, t1=1.1
    
    CHECK(w.evaluate(0.0) == 0.0);
    CHECK(w.evaluate(0.5) == 0.0);
    CHECK(w.evaluate(0.999) == 0.0);
}

TEST_CASE("LinearWaveform: during ramp", "[waveform][linear]") {
    LinearWaveform w{0.0, 100.0, 1.0, 1.1, true};  // start=0, end=100, t0=1.0, t1=1.1
    
    CHECK_THAT(w.evaluate(1.0), WithinAbs(0.0, 1e-9));
    CHECK_THAT(w.evaluate(1.05), WithinAbs(50.0, 1e-6));
    CHECK_THAT(w.evaluate(1.1), WithinAbs(100.0, 1e-9));
}

TEST_CASE("LinearWaveform: after end time (clamped)", "[waveform][linear]") {
    LinearWaveform w{0.0, 100.0, 1.0, 1.1, true};  // start=0, end=100, t0=1.0, t1=1.1
    
    CHECK_THAT(w.evaluate(1.1), WithinAbs(100.0, 1e-9));
    CHECK_THAT(w.evaluate(2.0), WithinAbs(100.0, 1e-9));
    CHECK_THAT(w.evaluate(10.0), WithinAbs(100.0, 1e-9));
}

TEST_CASE("LinearWaveform: after end time (unclamped)", "[waveform][linear]") {
    LinearWaveform w{0.0, 100.0, 1.0, 1.1, false};  // start=0, end=100, t0=1.0, t1=1.1
    
    CHECK_THAT(w.evaluate(1.1), WithinAbs(0.0, 1e-9));
    CHECK_THAT(w.evaluate(2.0), WithinAbs(0.0, 1e-9));
}

TEST_CASE("LinearWaveform: negative ramp", "[waveform][linear]") {
    LinearWaveform w{100.0, 0.0, 0.0, 0.01, true};  // start=100, end=0, t0=0.0, t1=0.01
    
    CHECK_THAT(w.evaluate(0.0), WithinAbs(100.0, 1e-9));
    CHECK_THAT(w.evaluate(0.005), WithinAbs(50.0, 1e-6));
    CHECK_THAT(w.evaluate(0.01), WithinAbs(0.0, 1e-9));
}

TEST_CASE("QuadraticWaveform: before start time", "[waveform][quadratic]") {
    QuadraticWaveform w{10.0, 5.0, 2.0, 1.0, 2.0};
    
    CHECK(w.evaluate(0.0) == 10.0);
    CHECK(w.evaluate(0.5) == 10.0);
}

TEST_CASE("QuadraticWaveform: during interval", "[waveform][quadratic]") {
    QuadraticWaveform w{10.0, 5.0, 2.0, 1.0, 2.0};
    
    // y(t) = 10 + 5*t + 2*t²
    CHECK_THAT(w.evaluate(1.0), WithinAbs(17.0, 1e-9));  // 10 + 5*1 + 2*1² = 17
    CHECK_THAT(w.evaluate(1.5), WithinAbs(22.0, 1e-6));  // 10 + 5*1.5 + 2*1.5² = 10 + 7.5 + 4.5 = 22
}

TEST_CASE("QuadraticWaveform: after end time", "[waveform][quadratic]") {
    QuadraticWaveform w{10.0, 5.0, 2.0, 1.0, 2.0};
    
    // After end_time, quadratic waveform returns constant term a
    // But wait - at t=2.0 we're AT end_time, so still evaluate: 10 + 5*2 + 2*4 = 28
    // After end_time (t > 2.0), return a
    CHECK_THAT(w.evaluate(2.5), WithinAbs(10.0, 1e-9));
    CHECK(w.evaluate(5.0) == 10.0);
}

TEST_CASE("SinusoidalWaveform: at key phases", "[waveform][sinusoidal]") {
    // y = 250 + 50*sin(2π*100*t + 0)
    SinusoidalWaveform w{250.0, 50.0, 100.0, 0.0};
    
    CHECK_THAT(w.evaluate(0.0), WithinAbs(250.0, 1e-9));      // sin(0) = 0
    CHECK_THAT(w.evaluate(0.0025), WithinAbs(300.0, 1e-6));   // sin(π/2) = 1
    CHECK_THAT(w.evaluate(0.005), WithinAbs(250.0, 1e-6));    // sin(π) = 0
    CHECK_THAT(w.evaluate(0.0075), WithinAbs(200.0, 1e-6));   // sin(3π/2) = -1
}

TEST_CASE("SinusoidalWaveform: with phase offset", "[waveform][sinusoidal]") {
    // y = 0 + 1*sin(2π*1*t + π/2) = cos(2π*t)
    SinusoidalWaveform w{0.0, 1.0, 1.0, M_PI/2.0};
    
    CHECK_THAT(w.evaluate(0.0), WithinAbs(1.0, 1e-9));     // cos(0) = 1
    CHECK_THAT(w.evaluate(0.25), WithinAbs(0.0, 1e-6));    // cos(π/2) = 0
    CHECK_THAT(w.evaluate(0.5), WithinAbs(-1.0, 1e-6));    // cos(π) = -1
}

TEST_CASE("PulsedWaveform: before pulse", "[waveform][pulsed]") {
    PulsedWaveform w{10.0, 100.0, 1.0, 0.5};
    
    CHECK(w.evaluate(0.0) == 10.0);
    CHECK(w.evaluate(0.999) == 10.0);
}

TEST_CASE("PulsedWaveform: during pulse", "[waveform][pulsed]") {
    PulsedWaveform w{10.0, 100.0, 1.0, 0.5};
    
    CHECK(w.evaluate(1.0) == 100.0);
    CHECK(w.evaluate(1.25) == 100.0);
    CHECK(w.evaluate(1.499) == 100.0);
}

TEST_CASE("PulsedWaveform: after pulse", "[waveform][pulsed]") {
    PulsedWaveform w{10.0, 100.0, 1.0, 0.5};
    
    CHECK(w.evaluate(1.5) == 10.0);
    CHECK(w.evaluate(2.0) == 10.0);
    CHECK(w.evaluate(10.0) == 10.0);
}

TEST_CASE("ArbitraryWaveform: linear interpolation", "[waveform][arbitrary]") {
    ArbitraryWaveform w;
    w.times_s = {0.0, 1.0, 2.0, 3.0};
    w.values = {0.0, 100.0, 50.0, 200.0};
    w.interp = ArbitraryWaveform::Interpolation::Linear;
    
    // At exact points
    CHECK_THAT(w.evaluate(0.0), WithinAbs(0.0, 1e-9));
    CHECK_THAT(w.evaluate(1.0), WithinAbs(100.0, 1e-9));
    CHECK_THAT(w.evaluate(2.0), WithinAbs(50.0, 1e-9));
    CHECK_THAT(w.evaluate(3.0), WithinAbs(200.0, 1e-9));
    
    // Interpolated
    CHECK_THAT(w.evaluate(0.5), WithinAbs(50.0, 1e-6));
    CHECK_THAT(w.evaluate(1.5), WithinAbs(75.0, 1e-6));
    CHECK_THAT(w.evaluate(2.5), WithinAbs(125.0, 1e-6));
}

TEST_CASE("ArbitraryWaveform: step interpolation", "[waveform][arbitrary]") {
    ArbitraryWaveform w;
    w.times_s = {0.0, 1.0, 2.0};
    w.values = {10.0, 20.0, 30.0};
    w.interp = ArbitraryWaveform::Interpolation::Step;
    
    CHECK(w.evaluate(0.0) == 10.0);
    CHECK(w.evaluate(0.5) == 10.0);   // Hold previous (step returns v0)
    CHECK(w.evaluate(0.999) == 10.0);
    CHECK(w.evaluate(1.0) == 10.0);   // At boundary, step still returns previous v0
    CHECK(w.evaluate(1.5) == 20.0);   // Now in [1.0, 2.0) interval, returns v0=20
}

TEST_CASE("ArbitraryWaveform: out of bounds", "[waveform][arbitrary]") {
    ArbitraryWaveform w;
    w.times_s = {1.0, 2.0, 3.0};
    w.values = {100.0, 200.0, 300.0};
    w.interp = ArbitraryWaveform::Interpolation::Linear;
    
    // Before first point
    CHECK(w.evaluate(0.0) == 100.0);
    CHECK(w.evaluate(0.5) == 100.0);
    
    // After last point
    CHECK(w.evaluate(3.0) == 300.0);
    CHECK(w.evaluate(5.0) == 300.0);
}

TEST_CASE("Waveform variant: evaluate via std::visit", "[waveform][variant]") {
    Waveform w1;
    w1.data = ConstantWaveform{42.0};
    CHECK(w1.evaluate(0.0) == 42.0);
    
    Waveform w2;
    w2.data = LinearWaveform{0.0, 100.0, 0.0, 1.0, true};  // 0→100 over [0, 1]
    CHECK_THAT(w2.evaluate(0.5), WithinAbs(50.0, 1e-6));
    
    Waveform w3;
    w3.data = SinusoidalWaveform{250.0, 50.0, 100.0, 0.0};
    CHECK_THAT(w3.evaluate(0.0), WithinAbs(250.0, 1e-9));
}

TEST_CASE("ValueOrWaveform: constant value", "[waveform][value_or_waveform]") {
    ValueOrWaveform val;
    val.constant_value = 123.45;
    
    CHECK(val.is_valid());
    CHECK_FALSE(val.is_time_varying());
    
    std::map<std::string, Waveform> lib;
    CHECK(val.evaluate(0.0, lib) == 123.45);
    CHECK(val.evaluate(10.0, lib) == 123.45);
}

TEST_CASE("ValueOrWaveform: inline waveform", "[waveform][value_or_waveform]") {
    ValueOrWaveform val;
    Waveform w;
    w.data = LinearWaveform{0.0, 100.0, 0.0, 1.0, true};  // 0→100 over [0, 1]
    val.waveform = w;
    
    CHECK(val.is_valid());
    CHECK(val.is_time_varying());
    
    std::map<std::string, Waveform> lib;
    CHECK_THAT(val.evaluate(0.0, lib), WithinAbs(0.0, 1e-9));
    CHECK_THAT(val.evaluate(0.5, lib), WithinAbs(50.0, 1e-6));
    CHECK_THAT(val.evaluate(1.0, lib), WithinAbs(100.0, 1e-9));
}

TEST_CASE("ValueOrWaveform: waveform reference", "[waveform][value_or_waveform]") {
    ValueOrWaveform val;
    val.waveform_ref = "my_ramp";
    
    CHECK(val.is_valid());
    CHECK(val.is_time_varying());
    
    // Create library with referenced waveform
    std::map<std::string, Waveform> lib;
    Waveform w;
    w.id = "my_ramp";
    w.data = LinearWaveform{0.0, 200.0, 0.0, 2.0, true};  // 0→200 over [0, 2]
    lib["my_ramp"] = w;
    
    CHECK_THAT(val.evaluate(0.0, lib), WithinAbs(0.0, 1e-9));
    CHECK_THAT(val.evaluate(1.0, lib), WithinAbs(100.0, 1e-6));
    CHECK_THAT(val.evaluate(2.0, lib), WithinAbs(200.0, 1e-9));
}

TEST_CASE("ValueOrWaveform: invalid (no option set)", "[waveform][value_or_waveform]") {
    ValueOrWaveform val;
    
    CHECK_FALSE(val.is_valid());
}

TEST_CASE("ValueOrWaveform: invalid (multiple options set)", "[waveform][value_or_waveform]") {
    ValueOrWaveform val;
    val.constant_value = 100.0;
    val.waveform_ref = "some_waveform";
    
    CHECK_FALSE(val.is_valid());
}

TEST_CASE("ValueOrWaveform: missing waveform reference throws", "[waveform][value_or_waveform]") {
    ValueOrWaveform val;
    val.waveform_ref = "nonexistent";
    
    std::map<std::string, Waveform> lib;  // Empty library
    
    CHECK_THROWS_AS(val.evaluate(0.0, lib), std::runtime_error);
}

TEST_CASE("ArbitraryWaveform: empty arrays throw", "[waveform][arbitrary][error]") {
    ArbitraryWaveform w;
    w.times_s = {};
    w.values = {};
    
    CHECK_THROWS_AS(w.evaluate(0.0), std::runtime_error);
}

TEST_CASE("ArbitraryWaveform: mismatched array sizes throw", "[waveform][arbitrary][error]") {
    ArbitraryWaveform w;
    w.times_s = {0.0, 1.0, 2.0};
    w.values = {10.0, 20.0};  // One element missing
    
    CHECK_THROWS_AS(w.evaluate(0.0), std::runtime_error);
}
