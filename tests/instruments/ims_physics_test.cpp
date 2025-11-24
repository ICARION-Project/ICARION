// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <vector>
#include <chrono>

#include "helpers/physics_sim_utils.h"
#include "core/physics/collisions/geometryUtils.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::tests;

namespace {

config::FullConfig make_ims_base(config::CollisionModel model, double dt_s, double total_s) {
    config::FullConfig cfg;
    cfg.simulation.dt_s = dt_s;
    cfg.simulation.total_time_s = total_s;
    cfg.simulation.write_interval = 1000000;
    cfg.physics.collision_model = model;
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "ims_physics.h5";
    cfg.output.print_progress = false;
    return cfg;
}

config::SpeciesProperties make_species(const std::string& id, double mass_amu, int charge_e) {
    config::SpeciesProperties sp;
    sp.id = id;
    sp.mass_amu = mass_amu;
    sp.charge = charge_e;
    sp.CCS_A2 = 200.0;
    sp.mobility_cm2Vs = 1.5;
    sp.convert_to_SI();
    return sp;
}

config::DomainConfig make_ims_domain(double length_m, double radius_m, double E_field, double pressure, double temperature) {
    config::DomainConfig dom = make_single_domain(config::Instrument::IMS);
    dom.geometry.length_m = length_m;
    dom.geometry.radius_m = radius_m;
    // Initialize all field values
    dom.fields.dc.axial_V.constant_value = E_field * length_m;
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 0.0;
    dom.fields.dc.EN_Td.constant_value = 0.0;
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 0.0;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 0.0;
    dom.fields.ac.compute_derived();
    dom.environment.pressure_Pa = pressure;
    dom.environment.temperature_K = temperature;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    return dom;
}

core::IonState make_ims_ion(const config::SpeciesProperties& sp) {
    core::IonState ion;
    ion.species_id = sp.id;
    ion.mass_kg = sp.mass_kg;
    ion.ion_charge_C = sp.charge_C;
    ion.CCS_m2 = sp.CCS_m2;
    ion.reduced_mobility_cm2_Vs = sp.mobility_cm2Vs.value_or(1.5);
    ion.pos = {0.0, 0.0, 0.0};
    ion.vel = {0.0, 0.0, 0.0};
    return ion;
}

} // namespace

struct ModelCase {
    config::CollisionModel model;
    std::string name;
    bool enable_ou;
};

TEST_CASE("IMS drift-based mobility across collision models", "[instrument][ims][drift]") {
    constexpr double pressure_Pa = 200.0;
    constexpr double temperature_K = 300.0;
    constexpr double E_field = 200.0;
    constexpr double length_m = 0.05;
    constexpr double radius_m = 1.0;


    std::vector<ModelCase> cases = {
        {config::CollisionModel::Friction, "Friction+OU", false},
        {config::CollisionModel::Langevin, "Langevin+OU", false},
        {config::CollisionModel::HSD, "HSD+OU", false},
        {config::CollisionModel::HSS, "HSS", false},
        {config::CollisionModel::EHSS, "EHSS", false}
    };

    // Geometry for EHSS
    physics::GeometryMap geom_map;
    const std::filesystem::path data_dir =
        std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() / "data" / "molecules";
    geom_map = physics::load_geometry_map({"H3O+"}, data_dir.string());

    for (const auto& c : cases) {
        DYNAMIC_SECTION("Model=" << c.name) {
            config::FullConfig cfg = make_ims_base(c.model, 1e-9, 1e-6);
            cfg.output.folder = "/tmp";
            cfg.output.trajectory_file = "ims_" + c.name + ".h5";
            cfg.physics.enable_ou_thermalization = c.enable_ou;
            cfg.domains.push_back(make_ims_domain(length_m, radius_m, E_field, pressure_Pa, temperature_K));

            auto sp = make_species("H3O+", 19.0, 1);
            // Provide CCS for HSS/HSD/Langevin friction gamma calc
            sp.CCS_A2 = 24.9;
            sp.mobility_cm2Vs = 24.0;
            sp.convert_to_SI();
            cfg.species_db.species[sp.id] = sp;

            config::IonSpeciesConfig ion_spec;
            ion_spec.species_id = sp.id;
            ion_spec.count = 50;
            cfg.ions.species.push_back(ion_spec);

            // Ion state at entrance
            auto ion = make_ims_ion(sp);
            ion.reduced_mobility_cm2_Vs = sp.mobility_cm2Vs.value_or(2.0);
            ion.CCS_m2 = sp.CCS_m2;

            double gamma_for_ou = 0.0;
            if (c.enable_ou) {
                // Use friction formula for gamma: q/(K*m)
                const double n = cfg.domains.front().environment.particle_density_m_3;
                const double mobility = sp.mobility_cm2Vs.value_or(2.0) * 1e-4 * LOSCHMIDT_CONSTANT / n;
                gamma_for_ou = ion.ion_charge_C / (mobility * ion.mass_kg); // match damping to mobility
            }

            const physics::GeometryMap* gm_ptr = (c.model == config::CollisionModel::EHSS) ? &geom_map : nullptr;

            // Build ensemble
            std::vector<core::IonState> ions;
            ions.reserve(ion_spec.count);
            for (size_t i = 0; i < ion_spec.count; ++i) {
                ions.push_back(ion);
            }

            auto sim = run_simple_simulation(cfg, ions, false, gm_ptr, gamma_for_ou);

            REQUIRE(sim.ions.size() == ion_spec.count);
            double z_mean = 0.0;
            for (const auto& s : sim.ions) {
                z_mean += s.pos.z;
            }
            z_mean /= static_cast<double>(sim.ions.size());
            REQUIRE(z_mean > 0.0);
            const double drift_velocity = z_mean / cfg.simulation.total_time_s;
            const double mobility_sim = drift_velocity / E_field;

            const double n = cfg.domains.front().environment.particle_density_m_3;
            const double mobility_expected = sp.mobility_cm2Vs.value_or(2.0) * 1e-4 * LOSCHMIDT_CONSTANT / n;

            REQUIRE(mobility_sim == Catch::Approx(mobility_expected).margin(0.50 * mobility_expected));
        }
    }
}
