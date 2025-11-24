// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <vector>

#include "helpers/physics_sim_utils.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::tests;
using Catch::Approx;

namespace {

config::FullConfig make_orbitrap_base(double dt_s, double total_time_s) {
    config::FullConfig cfg;
    cfg.simulation.dt_s = dt_s;
    cfg.simulation.total_time_s = total_time_s;
    cfg.simulation.write_interval = 1000000;
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "orbitrap_physics.h5";
    cfg.output.print_progress = false;
    return cfg;
}

config::SpeciesProperties make_species(const std::string& id, double mass_amu, int charge_e) {
    config::SpeciesProperties sp;
    sp.id = id;
    sp.mass_amu = mass_amu;
    sp.charge = charge_e;
    sp.CCS_A2 = 150.0;
    sp.mobility_cm2Vs = 1.0;
    sp.convert_to_SI();
    return sp;
}

config::DomainConfig make_orbitrap_domain() {
    config::DomainConfig dom = make_single_domain(config::Instrument::Orbitrap);
    dom.geometry.length_m = 0.05;
    dom.geometry.radius_in_m = 0.010;
    dom.geometry.radius_out_m = 0.015;
    dom.geometry.radius_char_m = 0.014;  // ensure positive k denominator
    dom.fields.dc.radial_V = 1200.0;
    dom.environment.pressure_Pa = 1e-6;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "He";
    dom.environment.compute_derived_properties();
    return dom;
}

} // namespace

TEST_CASE("Orbitrap axial frequency matches analytic (20 oscillations)", "[instrument][orbitrap][field]") {
    // Base config (dt/total set after period calculation)
    config::FullConfig cfg = make_orbitrap_base(0.0, 0.0);
    cfg.domains.push_back(make_orbitrap_domain());

    auto sp = make_species("Orbitrap_Test_Ion", 200.0, 1);
    cfg.species_db.species[sp.id] = sp;

    const auto& dom = cfg.domains.front();
    const double r_char_sq = dom.geometry.radius_char_m * dom.geometry.radius_char_m;
    const double k = 2.0 * dom.fields.dc.radial_V /
        (r_char_sq * std::log(dom.geometry.radius_out_m / dom.geometry.radius_in_m)
         - 0.5 * (dom.geometry.radius_out_m * dom.geometry.radius_out_m
                  - dom.geometry.radius_in_m * dom.geometry.radius_in_m));
    const double omega_z = std::sqrt((sp.charge_C / sp.mass_kg) * k);
    const double period = 2.0 * M_PI / omega_z;

    cfg.simulation.total_time_s = 20.0 * period;
    cfg.simulation.dt_s = period / 200.0;
    cfg.simulation.write_interval = 1000000;

    // Ion state with axial displacement
    core::IonState ion;
    ion.species_id = sp.id;
    ion.mass_kg = sp.mass_kg;
    ion.ion_charge_C = sp.charge_C;
    ion.CCS_m2 = sp.CCS_m2;
    ion.pos = {cfg.domains.front().geometry.radius_char_m * 0.6, 0.0,
               0.5 * cfg.domains.front().geometry.length_m + 2e-3};
    // Launch with 1000 eV kinetic energy along x
    const double energy_J = 1000.0 * ELEM_CHARGE_C;
    const double v_mag = std::sqrt(2.0 * energy_J / ion.mass_kg);
    ion.vel = {v_mag, 0.0, -500.0};  // axial kick to enforce crossings

    auto sim = run_simple_simulation(cfg, {ion}, /*capture_trace=*/true);

    // Zero-crossing based period estimate
    const auto& t = sim.trace.times;
    const auto& z = sim.trace.z_positions;
    const double z_center_offset = 0.5 * cfg.domains.front().geometry.length_m;
    std::vector<double> crossings;
    for (size_t i = 1; i < z.size(); ++i) {
        const double z_prev = z[i-1] - z_center_offset;
        const double z_curr = z[i]   - z_center_offset;
        if (z_prev * z_curr < 0.0) {
            const double frac = z_prev / (z_prev - z_curr);
            const double tc = t[i-1] + frac * (t[i] - t[i-1]);
            crossings.push_back(tc);
        }
    }

    REQUIRE(crossings.size() >= 10);  // several oscillations observed
    std::vector<double> periods;
    for (size_t i = 2; i < crossings.size(); i += 2) {
        periods.push_back(crossings[i] - crossings[i-2]);  // two zero crossings per period
    }
    double period_sim = 0.0;
    for (double p : periods) period_sim += p;
    period_sim /= periods.size();

    REQUIRE(period_sim == Approx(period).margin(0.05 * period));
}

TEST_CASE("Orbitrap trapping keeps ion radially confined", "[instrument][orbitrap][confinement]") {
    config::FullConfig cfg = make_orbitrap_base(0.0, 0.0);
    cfg.domains.push_back(make_orbitrap_domain());

    auto sp = make_species("Orbitrap_Test_Ion", 200.0, 1);
    cfg.species_db.species[sp.id] = sp;
    config::IonSpeciesConfig ion_spec;
    ion_spec.species_id = sp.id;
    ion_spec.count = 1;
    cfg.ions.species.push_back(ion_spec);

    const auto& dom = cfg.domains.front();
    const double r_char_sq = dom.geometry.radius_char_m * dom.geometry.radius_char_m;
    const double k = 2.0 * dom.fields.dc.radial_V /
        (r_char_sq * std::log(dom.geometry.radius_out_m / dom.geometry.radius_in_m)
         - 0.5 * (dom.geometry.radius_out_m * dom.geometry.radius_out_m
                  - dom.geometry.radius_in_m * dom.geometry.radius_in_m));
    const double omega_z = std::sqrt((sp.charge_C / sp.mass_kg) * k);
    const double period = 2.0 * M_PI / omega_z;

    cfg.simulation.total_time_s = 10.0 * period;
    cfg.simulation.dt_s = period / 200.0;
    cfg.simulation.write_interval = 1000000;

    core::IonState ion;
    ion.species_id = sp.id;
    ion.mass_kg = sp.mass_kg;
    ion.ion_charge_C = sp.charge_C;
    ion.CCS_m2 = sp.CCS_m2;
    ion.pos = {dom.geometry.radius_char_m * 0.7, 0.0,
               0.5 * dom.geometry.length_m + 1e-3};
    const double energy_J = 500.0 * ELEM_CHARGE_C;
    const double v_mag = std::sqrt(2.0 * energy_J / ion.mass_kg);
    ion.vel = {v_mag, 0.0, -200.0};

    auto sim = run_simple_simulation(cfg, {ion}, /*capture_trace=*/false);

    REQUIRE(sim.ions.size() == 1);
    const auto& final = sim.ions.front();
    const double r_final = std::sqrt(final.pos.x * final.pos.x + final.pos.y * final.pos.y);
    REQUIRE(r_final < dom.geometry.radius_out_m);
    REQUIRE(final.active);
}
