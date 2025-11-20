// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include "helpers/physics_test_utils.h"
#include "core/constants/constants.h"

using namespace ICARION;
using namespace ICARION::tests;

namespace {

Json::Value make_quadrupole_config(bool enable_gpu, double rf_voltage, double total_time = 5e-5) {
    Json::Value cfg;
    cfg["title"] = "Quadrupole Instrument Test";

    Json::Value sim;
    // total_time is provided in seconds; convert to nanoseconds for schema
    sim["max_time_ns"] = total_time * 1e9;
    sim["timestep_ns"] = 0.5;
    sim["output_interval"] = 250;
    sim["enable_gpu"] = enable_gpu;
    sim["random_seed"] = 4242;
    cfg["simulation"] = sim;

    Json::Value env;
    env["pressure_Pa"] = 1e-10;
    env["temperature_K"] = 300.0;
    cfg["environment"] = env;

    Json::Value instr;
    instr["type"] = "Quadrupole";
    instr["rf_voltage_V"] = rf_voltage;
    instr["rf_frequency_Hz"] = 2e6;
    instr["dc_quad_V"] = 5.0;
    instr["dc_axial_V"] = 2.0;
    instr["radius_m"] = 0.004;
    instr["length_m"] = 0.1;
    cfg["instrument"] = instr;

    Json::Value ion;
    ion["mass_amu"] = 50.0;
    ion["charge_C"] = 1.0;
    ion["mobility_cm2_Vs"] = 1.5;
    ion["CCS_m2"] = 1e-18;
    ion["position_x"] = 5e-5;
    ion["position_y"] = 5e-5;
    ion["position_z"] = 0.0;
    ion["velocity_x"] = 0.0;
    ion["velocity_y"] = 0.0;
    ion["velocity_z"] = 100.0;

    Json::Value ions(Json::arrayValue);
    ions.append(ion);
    cfg["initial_ions"] = ions;

    Json::Value output;
    Json::Value trajectory;

    trajectory["filename"] = enable_gpu ? "/tmp/quadrupole_physics_test_gpu.h5"
                   : "/tmp/quadrupole_physics_test_cpu.h5";

    output["trajectory"] = trajectory;
    cfg["output"] = output;

    cfg["collisions"]["model"] = "none";
    return cfg;
}

double radial_distance(const core::IonState& ion) {
    return std::sqrt(ion.pos.x * ion.pos.x + ion.pos.y * ion.pos.y);
}

}  // namespace

TEST_CASE("Quadrupole stable trajectory", "[quad][physics][cpu]") {
    //for m/z 50 leads an amplitude of 250 V at 2 MHz to q parameter of 0.76 (stable)
    auto config = make_quadrupole_config(false, 250.0, 2.0e-5);
    config["initial_ions"][0]["mass_amu"] = 50.0;
    config["initial_ions"][0]["charge_C"] = ELEM_CHARGE_C;
    double radius = 0.004;
    auto capture = capture_cpu_time_series(config, 500);
    auto run = run_instrument_simulation(config);
    REQUIRE(!capture.summary.ions.empty());
    const double r_final = radial_distance(capture.summary.ions.front());
    INFO("Quadrupole CPU ns/ion-step = " << run.ns_per_ion_step);
    INFO("Quadrupole stable final radius = " << r_final);
    REQUIRE(r_final < radius);
}

TEST_CASE("Quadrupole unstable trajectory", "[quad][physics][cpu]") {
    //for m/z 50 leads an amplitude of 400 V at 2 MHz to q parameter of 1.22 (unstable)
    auto config = make_quadrupole_config(false, 400.0, 4.0e-5);
    config["initial_ions"][0]["mass_amu"] = 50.0;
    config["initial_ions"][0]["charge_C"] = ELEM_CHARGE_C;
    double radius = 0.004;
    auto capture = capture_cpu_time_series(config, 500);
    auto run = run_instrument_simulation(config);
    REQUIRE(!capture.summary.ions.empty());
    const double r_final = radial_distance(capture.summary.ions.front());
    INFO("Quadrupole (unstable) CPU ns/ion-step = " << run.ns_per_ion_step);
    INFO("Quadrupole unstable final radius = " << r_final);
    REQUIRE(r_final >= radius);
}

#ifdef USE_GPU_ACCEL
TEST_CASE("Quadrupole GPU vs CPU parity", "[quad][gpu][parity]") {
    // Use very short time for GPU parity test to avoid numerical drift
    auto cpu_config = make_quadrupole_config(false, 60.0);
    cpu_config["initial_ions"][0]["mass_amu"] = 50.0;
    cpu_config["initial_ions"][0]["charge_C"] = ELEM_CHARGE_C;
    auto gpu_config = make_quadrupole_config(true, 60.0);
    gpu_config["initial_ions"][0]["mass_amu"] = 50.0;
    gpu_config["initial_ions"][0]["charge_C"] = ELEM_CHARGE_C;
    // Enable GPU validation to trigger adaptive download interval for parity testing
    gpu_config["debug"]["validate_gpu"] = true;

    // Reduce timestep for higher numerical parity between CPU and GPU
    cpu_config["simulation"]["timestep_ns"] = 0.05; // 0.05 ns
    gpu_config["simulation"]["timestep_ns"] = 0.05;
    auto cpu_run = run_instrument_simulation(cpu_config);
    auto gpu_run = run_instrument_simulation(gpu_config);

    INFO("Quadrupole CPU ns/ion-step = " << cpu_run.ns_per_ion_step);
    INFO("Quadrupole GPU ns/ion-step = " << gpu_run.ns_per_ion_step);
    if (!cpu_run.ions.empty() && !gpu_run.ions.empty()) {
        const auto& c = cpu_run.ions.front();
        const auto& g = gpu_run.ions.front();
        INFO("CPU pos (x,y,z) = " << c.pos.x << ", " << c.pos.y << ", " << c.pos.z);
        INFO("GPU pos (x,y,z) = " << g.pos.x << ", " << g.pos.y << ", " << g.pos.z);
        INFO("CPU vel (x,y,z) = " << c.vel.x << ", " << c.vel.y << ", " << c.vel.z);
        INFO("GPU vel (x,y,z) = " << g.vel.x << ", " << g.vel.y << ", " << g.vel.z);
    }

    REQUIRE(rms_error_pos(cpu_run.ions, gpu_run.ions) < 1e-8);
    REQUIRE(rms_error_vel(cpu_run.ions, gpu_run.ions) < 1e-5);
}
#else
TEST_CASE("Quadrupole GPU vs CPU parity (skipped)", "[quad][gpu][parity]") {
    SUCCEED("GPU acceleration disabled at build time");
}
#endif
