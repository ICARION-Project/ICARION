// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/config/types/DomainConfig.h"

using namespace ICARION;
using namespace ICARION::physics;
using Catch::Approx;

namespace {

std::vector<core::IonState> make_cloud() {
    std::vector<core::IonState> ions;
    ions.resize(3);

    ions[0].pos = {0.0, 0.0, 0.0};
    ions[1].pos = {1e-3, 0.0, 0.0};
    ions[2].pos = {0.0, 1e-3, 0.0};

    const double q = 1.602176634e-19;
    const double m = 1.0;
    for (auto& ion : ions) {
        ion.ion_charge_C = q;
        ion.mass_kg = m;
        ion.active = true;
        ion.born = true;
    }
    return ions;
}

}  // namespace

TEST_CASE("SpaceChargeDirect SoA matches AoS", "[forces][soa]") {
    config::DomainConfig domain{};
    ForceRegistry registry(domain);
    registry.add_force(std::make_unique<SpaceChargeDirect>(1e-10));

    auto ions = make_cloud();
    auto ensemble = core::IonEnsemble::from_legacy(ions);

    ForceContext ctx_aos{};
    ctx_aos.all_ions = &ions;

    ForceContext ctx_soa{};
    ctx_soa.ion_ensemble = &ensemble;

    for (size_t i = 0; i < ions.size(); ++i) {
        ctx_soa.ion_index = i;
        Vec3 F_aos = registry.compute_total_force(ions[i], /*t=*/0.0, ctx_aos);
        Vec3 F_soa = registry.compute_total_force(ensemble, i, /*t=*/0.0, ctx_soa);

        REQUIRE(F_soa.x == Approx(F_aos.x).margin(1e-12));
        REQUIRE(F_soa.y == Approx(F_aos.y).margin(1e-12));
        REQUIRE(F_soa.z == Approx(F_aos.z).margin(1e-12));
    }
}
