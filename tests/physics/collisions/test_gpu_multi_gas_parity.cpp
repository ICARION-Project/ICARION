// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/types/CollisionTypes.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"

#include <vector>
#include <cmath>

using namespace ICARION;
using namespace ICARION::config;
using namespace ICARION::physics;
using ICARION::core::IonEnsemble;
using ICARION::core::IonState;
using ICARION::core::Vec3;

#ifdef ICARION_USE_GPU

namespace {

constexpr size_t NUM_IONS = 2048;
constexpr double T_K = 300.0;
constexpr double MASS_KG = 29.0 * 1.66053906660e-27;  // Representative ion mass (H3O+)
constexpr double CCS_M2 = 45e-20;
constexpr double DT = 1e-7;
constexpr int SEED_BASE = 1337;

EnvironmentConfig make_mixture_env(double he_fraction) {
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 101325.0;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};

    GasMixtureComponent he;
    he.species = "He";
    he.mole_fraction = he_fraction;
    he.mass_kg = MOLAR_MASS_HE_KG;
    he.radius_m = RADIUS_HE_M;
    he.cross_section_m2 = 35e-20;
    env.gas_mixture.push_back(he);

    GasMixtureComponent n2;
    n2.species = "N2";
    n2.mole_fraction = 1.0 - he_fraction;
    n2.mass_kg = MOLAR_MASS_N2_KG;
    n2.radius_m = RADIUS_N2_M;
    n2.cross_section_m2 = 60e-20;
    env.gas_mixture.push_back(n2);

    env.compute_derived_properties();
    return env;
}

std::vector<IonState> make_ions(double velocity_scale) {
    std::vector<IonState> ions(NUM_IONS);
    const double v_thermal = std::sqrt(3.0 * BOLTZMANN_CONSTANT * T_K / MASS_KG);
    const double v_init = v_thermal * std::sqrt(velocity_scale);

    for (auto& ion : ions) {
        ion.species_id = "H3O+";
        ion.mass_kg = MASS_KG;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = CCS_M2;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v_init, 0.0, 0.0};
        ion.active = true;
    }
    return ions;
}

double compute_energy_ratio(const std::vector<IonState>& ions) {
    double sum_v2 = 0.0;
    for (const auto& ion : ions) {
        sum_v2 += ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
    }
    const double mean_v2 = sum_v2 / ions.size();
    const double KE_eV = 0.5 * MASS_KG * mean_v2 / ELEM_CHARGE_C;
    const double thermal_ke = 1.5 * BOLTZMANN_CONSTANT * T_K / ELEM_CHARGE_C;
    return KE_eV / thermal_ke;
}

double run_cpu_thermalization(const EnvironmentConfig& env,
                              double velocity_scale,
                              int n_steps) {
    HSSCollisionHandler handler(false, nullptr);
    std::vector<IonState> ions = make_ions(velocity_scale);

    for (size_t idx = 0; idx < ions.size(); ++idx) {
        IonState& ion = ions[idx];
        PhysicsRng rng(SEED_BASE + static_cast<int>(idx));

        for (int step = 0; step < n_steps; ++step) {
            std::vector<IonState> legacy{ion};
            auto ensemble = IonEnsemble::from_legacy(legacy);
            auto view = ensemble.collision_data(0);
            handler.handle_collision(view, DT, rng, env);
            ion = ensemble.ion_state(0);
        }
    }

    return compute_energy_ratio(ions);
}

double run_gpu_thermalization(icarion::gpu::GPUContext& gpu_ctx,
                              const EnvironmentConfig& env,
                              double velocity_scale,
                              int n_steps) {
    auto helper = icarion::gpu::GPUCollisionHelper::create(gpu_ctx, 256, "HSS", SEED_BASE);
    REQUIRE(helper != nullptr);

    std::vector<IonState> ions = make_ions(velocity_scale);
    for (int step = 0; step < n_steps; ++step) {
        bool used_gpu = helper->process_collisions_batch(ions, DT, env);
        REQUIRE(used_gpu);  // ensure GPU path taken
    }

    return compute_energy_ratio(ions);
}

}  // namespace

TEST_CASE("GPU multi-gas HSS thermalization parity", "[collision][gpu][multigas]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    if (!gpu_ctx) {
        SKIP("GPU not available");
    }

    struct Scenario {
        double he_fraction;
        double velocity_scale;
        int steps;
    };

    const std::vector<Scenario> scenarios = {
        {0.7, 10.0, 2000},
        {0.3, 0.1, 2000},
        {0.5, 1.0, 3000},
    };

    for (const auto& scenario : scenarios) {
        CAPTURE(scenario.he_fraction);
        CAPTURE(scenario.velocity_scale);
        CAPTURE(scenario.steps);

        auto env = make_mixture_env(scenario.he_fraction);
        double cpu_ratio = run_cpu_thermalization(env, scenario.velocity_scale, scenario.steps);
        double gpu_ratio = run_gpu_thermalization(*gpu_ctx, env, scenario.velocity_scale, scenario.steps);

        REQUIRE(cpu_ratio > 0.4);
        REQUIRE(cpu_ratio < 1.6);
        REQUIRE(gpu_ratio > 0.4);
        REQUIRE(gpu_ratio < 1.6);

        double diff_pct = std::abs(cpu_ratio - gpu_ratio) / cpu_ratio * 100.0;
        CHECK(diff_pct < 35.0);
    }
}

#else

TEST_CASE("GPU multi-gas parity (skipped)", "[collision][gpu][multigas]") {
    SKIP("ICARION_USE_GPU not enabled");
}

#endif  // ICARION_USE_GPU
