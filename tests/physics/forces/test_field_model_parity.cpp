// SPDX-License-Identifier: GPL-3.0-only
// SPDX-FileCopyrightText: 2026 Christoph Schaefer
//
// Parity test: AnalyticalFieldModel vs. ElectricFieldForce analytical path.
// Ensures the extracted field model returns identical E-fields to the legacy
// inline formulas for each supported instrument type.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/config/types/AnalyticalFieldModel.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/utils/mathUtils.h"

using ICARION::config::DomainConfig;
using ICARION::config::Instrument;
using ICARION::config::AnalyticalFieldModel;
using ICARION::physics::ElectricFieldForce;
using ICARION::physics::ForceContext;
using ICARION::core::IonState;
using ICARION::core::Vec3;
using ICARION::core::IonEnsemble;
using Catch::Approx;

namespace {

DomainConfig make_base_domain(Instrument inst) {
    DomainConfig dom;
    dom.instrument = inst;
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    dom.domain_index = 0;
    return dom;
}

void require_parity(const DomainConfig& dom, const Vec3& pos, double t) {
    IonState ion;
    ion.pos = pos;
    ion.vel = Vec3{0.0, 0.0, 0.0};
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1.0;
    ion.reduced_mobility_cm2_Vs = 1.0;
    ion.CCS_m2 = 1.0;
    ion.born = true;
    ion.active = true;

    AnalyticalFieldModel model(dom);
    ElectricFieldForce force(dom);

    ForceContext ctx_legacy;
    ctx_legacy.domain = &dom;

    ForceContext ctx_model;
    ctx_model.domain = &dom;
    ctx_model.field_model = &model;

    IonEnsemble ens = IonEnsemble::from_legacy({ion});
    ctx_legacy.ion_ensemble = &ens;
    ctx_model.ion_ensemble = &ens;
    ctx_legacy.ion_index = 0;
    ctx_model.ion_index = 0;

    Vec3 F_legacy = force.compute(ens, 0, t, ctx_legacy);
    Vec3 F_model = force.compute(ens, 0, t, ctx_model);

    REQUIRE(F_model.x == Approx(F_legacy.x).margin(1e-9));
    REQUIRE(F_model.y == Approx(F_legacy.y).margin(1e-9));
    REQUIRE(F_model.z == Approx(F_legacy.z).margin(1e-9));
}

} // namespace

TEST_CASE("AnalyticalFieldModel matches ElectricFieldForce", "[field][parity]") {
    // LQIT
    {
        DomainConfig dom = make_base_domain(Instrument::LQIT);
        dom.geometry.length_m = 0.10;
        dom.geometry.radius_m = 0.01;
        dom.fields.rf.voltage_V = ICARION::config::ValueOrWaveform(100.0);
        dom.fields.rf.frequency_Hz = ICARION::config::ValueOrWaveform(1.0e6);
        dom.fields.ac.voltage_V = ICARION::config::ValueOrWaveform(5.0);
        dom.fields.ac.frequency_Hz = ICARION::config::ValueOrWaveform(2.0e5);
        dom.fields.dc.quad_V = ICARION::config::ValueOrWaveform(10.0);
        dom.fields.dc.axial_V = ICARION::config::ValueOrWaveform(20.0);
        dom.finalize();
        require_parity(dom, Vec3{0.002, -0.0015, 0.03}, 1.0e-6);
    }

    // IMS
    {
        DomainConfig dom = make_base_domain(Instrument::IMS);
        dom.geometry.length_m = 1.0;
        dom.geometry.radius_m = 0.02;
        dom.fields.dc.axial_V = ICARION::config::ValueOrWaveform(100.0);
        dom.fields.dc.quad_V = ICARION::config::ValueOrWaveform(5.0);
        dom.fields.rf.voltage_V = ICARION::config::ValueOrWaveform(50.0);
        dom.fields.rf.frequency_Hz = ICARION::config::ValueOrWaveform(8.0e5);
        dom.finalize();
        require_parity(dom, Vec3{0.001, -0.001, 0.6}, 5.0e-7);
    }

    // TOF
    {
        DomainConfig dom = make_base_domain(Instrument::TOF);
        dom.geometry.length_m = 1.5;
        dom.geometry.acc_length_m = 0.2;
        dom.geometry.radius_m = 0.05;
        dom.fields.dc.axial_V = ICARION::config::ValueOrWaveform(1200.0);
        dom.finalize();
        require_parity(dom, Vec3{0.0, 0.0, 0.05}, 0.0);
        require_parity(dom, Vec3{0.0, 0.0, 0.3}, 0.0); // drift region: E=0
    }

    // FTICR
    {
        DomainConfig dom = make_base_domain(Instrument::FTICR);
        dom.geometry.length_m = 0.5;
        dom.geometry.radius_m = 0.05;
        dom.fields.dc.radial_V = ICARION::config::ValueOrWaveform(15.0);
        dom.finalize();
        require_parity(dom, Vec3{0.01, -0.005, 0.22}, 0.0);
    }

    // Orbitrap
    {
        DomainConfig dom = make_base_domain(Instrument::Orbitrap);
        dom.geometry.length_m = 0.2;
        dom.geometry.radius_in_m = 0.02;
        dom.geometry.radius_out_m = 0.04;
        dom.geometry.radius_char_m = 0.03;
        dom.fields.dc.radial_V = ICARION::config::ValueOrWaveform(1500.0);
        dom.finalize();
        require_parity(dom, Vec3{0.015, 0.0, 0.01}, 0.0);
    }

    // QuadrupoleRF
    {
        DomainConfig dom = make_base_domain(Instrument::QuadrupoleRF);
        dom.geometry.length_m = 0.6;
        dom.geometry.radius_m = 0.01;
        dom.fields.rf.voltage_V = ICARION::config::ValueOrWaveform(180.0);
        dom.fields.rf.frequency_Hz = ICARION::config::ValueOrWaveform(1.5e6);
        dom.fields.dc.quad_V = ICARION::config::ValueOrWaveform(12.0);
        dom.fields.dc.axial_V = ICARION::config::ValueOrWaveform(30.0);
        dom.finalize();
        require_parity(dom, Vec3{0.003, -0.002, 0.4}, 2.0e-6);
    }
}
