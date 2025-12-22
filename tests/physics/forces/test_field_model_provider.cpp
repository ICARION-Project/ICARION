// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include "core/config/types/FieldProviderModel.h"
#include <catch2/catch_approx.hpp>
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/utils/mathUtils.h"
#include "fieldsolver/utils/IFieldProvider.h"

using ICARION::config::FieldProviderModel;
using ICARION::config::DomainConfig;
using ICARION::config::Instrument;
using ICARION::physics::ElectricFieldForce;
using ICARION::physics::ForceContext;
using ICARION::core::Vec3;
using ICARION::core::IonState;
using Catch::Approx;

namespace {

class ConstantFieldProvider : public ::IFieldProvider {
public:
    explicit ConstantFieldProvider(const Vec3& E) : E_(E) {}
    Vec3 get_E(const Vec3& /*pos*/) const override { return E_; }
private:
    Vec3 E_;
};

DomainConfig make_dummy_domain() {
    DomainConfig dom;
    dom.instrument = Instrument::IMS;
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    dom.geometry.length_m = 1.0;
    dom.geometry.radius_m = 0.05;
    dom.fields.dc.axial_V = ICARION::config::ValueOrWaveform(100.0);
    dom.finalize();
    return dom;
}

} // namespace

TEST_CASE("FieldProviderModel feeds ElectricFieldForce", "[field][provider]") {
    // Constant provider field
    Vec3 E_const{1.0, -2.0, 3.5};
    auto provider = std::make_shared<ConstantFieldProvider>(E_const);
    FieldProviderModel model(provider);

    DomainConfig dom = make_dummy_domain();
    ElectricFieldForce force(dom);  // constructed with domain (analytic fallback)

    IonState ion;
    ion.pos = Vec3{0.0, 0.0, 0.0};
    ion.vel = Vec3{0.0, 0.0, 0.0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1.0;
    ion.reduced_mobility_cm2_Vs = 1.0;
    ion.CCS_m2 = 1.0;

    ForceContext ctx{};
    ctx.field_model = &model;   // prefer model path
    ctx.domain = &dom;

    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F = force.compute(ens, 0, 0.0, ctx);
    // Force = q * E, q=1
    REQUIRE(F.x == Approx(E_const.x));
    REQUIRE(F.y == Approx(E_const.y));
    REQUIRE(F.z == Approx(E_const.z));
}
