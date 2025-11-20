// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "helpers/physics_test_utils.h"
#include "core/constants/constants.h"

using namespace ICARION;
using namespace ICARION::tests;

namespace {

constexpr double kIMSOutputDt = 1.0e-8;  // seconds
constexpr double kIMSRuntime = 2.0e-4;   // seconds
constexpr double kIMSElectricField = 200.0; // V/m
constexpr double kIMSPosTol = 1e-2;
constexpr double kIMSVelTol = 1e-5;

Json::Value make_ims_config(bool enable_gpu) {
    Json::Value cfg(Json::objectValue);

    cfg["metadata"]["description"] = "IMS Instrument Physics Test";
    cfg["metadata"]["purpose"] = "Regression test for IMS drift physics";

    auto& sim = cfg["simulation"];
    sim["timestep_ns"] = kIMSOutputDt * 1e9;
    sim["max_time_ns"] = kIMSRuntime * 1e9;
    sim["output_interval"] = 25;
    sim["enable_gpu"] = enable_gpu;
    sim["random_seed"] = 1337;

    cfg["integrator"]["type"] = "rk4";

    auto& instrument = cfg["instrument"];
    instrument["type"] = "IMS";
    instrument["length_m"] = 0.2;
    instrument["radius_m"] = 0.01;
    instrument["drift_field_V_m"] = kIMSElectricField;
    instrument["pressure_Pa"] = 200.0;
    instrument["temperature_K"] = 300.0;
    auto& geometry = instrument["geometry"];
    geometry["length_m"] = instrument["length_m"];
    geometry["radius_m"] = instrument["radius_m"];
    instrument["voltages"]["dc_axial_V"] = kIMSElectricField * instrument["length_m"].asDouble();

    auto& environment = cfg["environment"];
    environment["gas_species"] = "N2";
    environment["temperature_K"] = 300.0;
    environment["pressure_Pa"] = 200.0;
    auto& gas_velocity = environment["gas_velocity_m_s"];
    gas_velocity.append(0.0);
    gas_velocity.append(0.0);
    gas_velocity.append(0.0);

    auto& ions = cfg["ions"];
    Json::Value species(Json::arrayValue);
    Json::Value ion(Json::objectValue);
    ion["name"] = "IMS_Test_Ion+";
    ion["mass_amu"] = 40.0;
    ion["charge"] = 1;
    ion["initial_count"] = 1;
    ion["reduced_mobility_cm2_Vs"] = 1.5;
    ion["collision_cross_section_A2"] = 200.0;
    ion["ccs_m2"] = 2e-18;
    species.append(ion);
    ions["species"] = species;

    auto& distribution = ions["initial_distribution"];
    distribution["type"] = "point";
    auto& position = distribution["position"];
    position["type"] = "point";
    auto& center = position["center"];
    center.append(0.0);
    center.append(0.0);
    center.append(0.0);
    auto& velocity = distribution["velocity"];
    velocity["type"] = "fixed";
    auto& vel_value = velocity["value"];
    vel_value.append(0.0);
    vel_value.append(0.0);
    vel_value.append(0.0);

    auto& collisions = cfg["collisions"];
    collisions["enabled"] = true;
    collisions["model"] = "Friction";

    auto& output = cfg["output"];
    auto& trajectory_file = output["trajectory_file"];
    trajectory_file["filename"] = enable_gpu ? "/tmp/ims_physics_test_gpu.h5"
                                              : "/tmp/ims_physics_test_cpu.h5";

    return cfg;
}

double compute_expected_mobility(double total_time_s, double mass_amu) {
    const double charge_C = ELEM_CHARGE_C;
    const double mass_kg = mass_amu * AMU_TO_KG;
    return 0.5 * charge_C * total_time_s / mass_kg;
}

}  // namespace

TEST_CASE("IMS mobility CPU vs analytic", "[ims][physics][cpu]") {
    std::cout << "[IMS TEST] entering test body" << std::endl;
    auto config = make_ims_config(false);
    {
        Json::StreamWriterBuilder b; b["indentation"] = "";
        std::cout << "[IMS TEST] Config=" << Json::writeString(b, config) << std::endl;
    }
    auto run = run_instrument_simulation(config);

    REQUIRE(run.ions.size() == 1);

    // computes the simulated mobility via drifted distance over time and electric field
    const auto& ion = run.ions.front();
    const double delta_z = ion.pos.z;
    const double drift_velocity = delta_z / kIMSRuntime;
    const double particle_density_m_3 =
        (config["environment"]["pressure_Pa"].asDouble()) /
        (BOLTZMANN_CONSTANT * config["environment"]["temperature_K"].asDouble());
    std::cout << "Drifted distance: " << delta_z << " m over " << kIMSRuntime << " s" << std::endl;
    std::cout << "Drift velocity: " << drift_velocity << " m/s" << std::endl;
    std::cout << "Electric field: " << kIMSElectricField << " V/m" << std::endl;
    // Physics engine: ion_mobility = (reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / density
    // Config stores reduced mobility explicitly in ions.species[].reduced_mobility_cm2_Vs
    // Therefore: mobility_expected = (reduced_mobility_cm2_Vs * 1e-4) * LOSCHMIDT_CONSTANT / density
    const double mobility_sim = drift_velocity / kIMSElectricField;
    const double mobility_expected = (config["ions"]["species"][0]["reduced_mobility_cm2_Vs"].asDouble() * 1e-4) *
                                     LOSCHMIDT_CONSTANT / particle_density_m_3;
    std::cout << "Simulated mobility: " << mobility_sim << " m2/Vs" << std::endl;
    std::cout << "Expected mobility: " << mobility_expected << " m2/Vs" << std::endl;
    INFO("IMS CPU ns/ion-step = " << run.ns_per_ion_step);
    REQUIRE(std::fabs(mobility_sim - mobility_expected) / mobility_expected < 0.02);
}

#ifdef USE_GPU_ACCEL
TEST_CASE("IMS GPU vs CPU parity", "[ims][gpu][parity]") {
    auto cpu_config = make_ims_config(false);
    cpu_config["ions"]["species"][0]["reduced_mobility_cm2_Vs"] = 1.5;
    cpu_config["ions"]["species"][0]["ccs_m2"] = 2e-18;

    auto gpu_config = make_ims_config(true);
    gpu_config["ions"]["species"][0]["reduced_mobility_cm2_Vs"] = 1.5;
    gpu_config["ions"]["species"][0]["ccs_m2"] = 2e-18;
    // Enable GPU validation to trigger adaptive download interval for parity testing
    gpu_config["debug"]["validate_gpu"] = true;

    auto cpu_run = run_instrument_simulation(cpu_config);
    auto gpu_run = run_instrument_simulation(gpu_config);

    if (!cpu_run.ions.empty()) {
        const auto& cpu_ion = cpu_run.ions.front();
        INFO("CPU final pos=" << cpu_ion.pos.x << "," << cpu_ion.pos.y << "," << cpu_ion.pos.z
             << " vel=" << cpu_ion.vel.x << "," << cpu_ion.vel.y << "," << cpu_ion.vel.z);
    }
    if (!gpu_run.ions.empty()) {
        const auto& gpu_ion = gpu_run.ions.front();
        INFO("GPU final pos=" << gpu_ion.pos.x << "," << gpu_ion.pos.y << "," << gpu_ion.pos.z
             << " vel=" << gpu_ion.vel.x << "," << gpu_ion.vel.y << "," << gpu_ion.vel.z);
    }

    const double pos_rms = rms_error_pos(cpu_run.ions, gpu_run.ions);
    const double vel_rms = rms_error_vel(cpu_run.ions, gpu_run.ions);

    INFO("IMS CPU ns/ion-step = " << cpu_run.ns_per_ion_step);
    INFO("IMS GPU ns/ion-step = " << gpu_run.ns_per_ion_step);

    REQUIRE(pos_rms < kIMSPosTol);
    REQUIRE(vel_rms < kIMSVelTol);
}
#else
TEST_CASE("IMS GPU vs CPU parity (skipped)", "[ims][gpu][parity]") {
    SUCCEED("GPU acceleration disabled at build time");
}
#endif
