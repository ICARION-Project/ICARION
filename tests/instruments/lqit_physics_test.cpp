// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <numeric>

#include "helpers/physics_test_utils.h"
#include "core/constants/constants.h"

using namespace ICARION;
using namespace ICARION::tests;

namespace {

Json::Value make_lqit_config(bool enable_gpu) {
    Json::Value cfg;
    cfg["title"] = "LQIT Instrument Test";

    Json::Value sim;
    sim["max_time_ns"] = 1.0e+05;
    sim["timestep_ns"] = 2.5;
    sim["output_interval"] = 100;
    sim["enable_gpu"] = enable_gpu;
    sim["random_seed"] = 777;
    cfg["simulation"] = sim;

    Json::Value instr;
    instr["type"] = "LQIT";
    instr["rf_voltage_V"] = 30.0;
    instr["rf_frequency_Hz"] = 1.5e6;
    instr["dc_axial_V"] = 1.0;
    instr["dc_quad_V"] = 0.0;
    instr["radius_m"] = 0.004;
    instr["length_m"] = 0.05;
    Json::Value ac;
    ac["voltage_V"] = 0.0;
    ac["frequency_Hz"] = 1.0e5;
    instr["ac_excitation"] = ac;
    cfg["instrument"] = instr;

    Json::Value ion;
    ion["mass_amu"] = 40.0;
    // Use SI charge key for GPU/CPU parity
    ion["charge_C"] = ELEM_CHARGE_C;
    ion["position_x"] = 5e-5;
    ion["position_y"] = 0.0;
    ion["position_z"] = 0.0;
    ion["velocity_x"] = 0.0;
    ion["velocity_y"] = 20.0;
    ion["velocity_z"] = 0.0;

    Json::Value ions(Json::arrayValue);
    ions.append(ion);
    cfg["initial_ions"] = ions;

    Json::Value output;
    Json::Value trajectory;

    trajectory["filename"] = enable_gpu ? "/tmp/lqit_physics_test_gpu.h5"
                   : "/tmp/lqit_physics_test_cpu.h5";

    output["trajectory"] = trajectory;
    cfg["output"] = output;

    cfg["collisions"]["model"] = "none";
    return cfg;
}

double compute_q_parameter(double charge_C,
                           double mass_kg,
                           double rf_voltage,
                           double rf_frequency_Hz,
                           double radius_m) {
    const double omega = 2.0 * M_PI * rf_frequency_Hz;
    return (4.0 * charge_C * rf_voltage) / (mass_kg * omega * omega * radius_m * radius_m);
}

double expected_secular_frequency(double charge_C,
                                  double mass_kg,
                                  double rf_voltage,
                                  double rf_frequency_Hz,
                                  double radius_m,
                                  double dc_quad_V) {
    const double omega = 2.0 * M_PI * rf_frequency_Hz;
    const double q_param = compute_q_parameter(charge_C, mass_kg, rf_voltage, rf_frequency_Hz, radius_m);
    const double a_param = (4.0 * charge_C * dc_quad_V) / (mass_kg * omega * omega * radius_m * radius_m);
    const double beta = std::sqrt(a_param + (q_param * q_param) / 2.0);
    return 0.5 * beta * omega;
}

double estimate_frequency_from_series(const TimeSeriesResult& capture) {
    std::vector<double> zero_crossings;
    if (capture.snapshots.size() < 2) return 0.0;
    for (size_t i = 1; i < capture.snapshots.size(); ++i) {
        const auto& prev = capture.snapshots[i - 1];
        const auto& curr = capture.snapshots[i];
        const double x_prev = prev.ions.front().pos.x;
        const double x_curr = curr.ions.front().pos.x;
        const double vx_curr = curr.ions.front().vel.x;
        if (x_prev <= 0.0 && x_curr > 0.0 && vx_curr > 0.0) {
            const double frac = x_prev / (x_prev - x_curr);
            zero_crossings.push_back(prev.time + frac * (curr.time - prev.time));
        }
    }
    if (zero_crossings.size() < 3) return 0.0;
    std::vector<double> periods;
    for (size_t i = 2; i < zero_crossings.size(); ++i) {
        periods.push_back(zero_crossings[i] - zero_crossings[i - 1]);
    }
    const double avg_period = std::accumulate(periods.begin(), periods.end(), 0.0) / periods.size();
    return 2.0 * M_PI / avg_period;
}

}  // namespace

TEST_CASE("LQIT secular frequency matches Mathieu theory", "[lqit][physics][cpu]") {
    auto config = make_lqit_config(false);
    config["initial_ions"][0]["mass_amu"] = 40.0;
    config["initial_ions"][0]["charge"] = 1.0;
    const auto capture = capture_cpu_time_series(config, 20);
    const auto summary = run_instrument_simulation(config);
    INFO("LQIT CPU ns/ion-step = " << summary.ns_per_ion_step);
    REQUIRE(capture.snapshots.size() > 3);

    const double freq_sim = estimate_frequency_from_series(capture);
    REQUIRE(freq_sim > 0.0);

    const double mass_amu = config["initial_ions"][0]["mass_amu"].asDouble();
    const double mass_kg = mass_amu * AMU_TO_KG;
    const double rf_voltage = config["instrument"]["rf_voltage_V"].asDouble();
    const double rf_freq = config["instrument"]["rf_frequency_Hz"].asDouble();
    const double radius = config["instrument"]["radius_m"].asDouble();
    const double dc_quad = config["instrument"]["dc_quad_V"].asDouble();
    const double freq_expected = expected_secular_frequency(ELEM_CHARGE_C, mass_kg, rf_voltage, rf_freq, radius, dc_quad);

    REQUIRE(std::fabs(freq_sim - freq_expected) / freq_expected < 0.05);
}

#ifdef USE_GPU_ACCEL
TEST_CASE("LQIT GPU vs CPU parity", "[lqit][gpu][parity]") {
    auto cpu_config = make_lqit_config(false);
    // Tight parity: reduce timestep
    cpu_config["simulation"]["timestep_ns"] = 0.05; // 0.05 ns
    auto gpu_config = make_lqit_config(true);
    gpu_config["simulation"]["timestep_ns"] = 0.05; // match CPU
    // Enable GPU validation to trigger adaptive download interval for parity testing
    gpu_config["debug"]["validate_gpu"] = true;

    auto cpu_run = run_instrument_simulation(cpu_config);
    auto gpu_run = run_instrument_simulation(gpu_config);

    INFO("LQIT CPU ns/ion-step = " << cpu_run.ns_per_ion_step);
    INFO("LQIT GPU ns/ion-step = " << gpu_run.ns_per_ion_step);

    REQUIRE(rms_error_pos(cpu_run.ions, gpu_run.ions) < 1e-8);
    REQUIRE(rms_error_vel(cpu_run.ions, gpu_run.ions) < 1e-5);
}
#else
TEST_CASE("LQIT GPU vs CPU parity (skipped)", "[lqit][gpu][parity]") {
    SUCCEED("GPU acceleration disabled at build time");
}
#endif
