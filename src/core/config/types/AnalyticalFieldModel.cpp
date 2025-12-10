// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "AnalyticalFieldModel.h"
#include "core/utils/mathUtils.h"
#include <cmath>
#include <algorithm>

namespace {
    constexpr double MIN_VOLTAGE_THRESHOLD = 1e-12;

    inline double eval_value(const ICARION::config::ValueOrWaveform& val, double t,
                             const std::map<std::string, ICARION::config::Waveform>& lib) {
        if (val.constant_value.has_value()) {
            return val.constant_value.value();
        }
        return val.evaluate(t, lib);
    }
}

namespace ICARION::config {

Vec3 AnalyticalFieldModel::E(const Vec3& global_pos, double t) const {
    if (!domain_) return Vec3{0.0, 0.0, 0.0};
    using Inst = ICARION::instrument::InstrumentType;
    switch (domain_->instrument) {
        case Inst::LQIT:         return compute_lqit_field(global_pos, t);
        case Inst::IMS:          return compute_ims_field(global_pos, t);
        case Inst::TOF:          return compute_tof_field(global_pos);
        case Inst::FTICR:        return compute_fticr_field(global_pos);
        case Inst::Orbitrap:     return compute_orbitrap_field(global_pos, t);
        case Inst::QuadrupoleRF: return compute_quadrupole_rf_field(global_pos, t);
        default:                 return Vec3{0.0, 0.0, 0.0};
    }
}

Vec3 AnalyticalFieldModel::compute_lqit_field(const Vec3& pos, double t) const {
    Vec3 E_total{0.0, 0.0, 0.0};
    const auto& rf = domain_->fields.rf;
    const auto& dc = domain_->fields.dc;
    const auto& ac = domain_->fields.ac;
    const auto& geom = domain_->geometry;
    const auto& lib = domain_->fields.waveform_library;

    const double rf_voltage = (rf.voltage_V.constant_value.has_value() || rf.voltage_V.waveform_ref.has_value() || rf.voltage_V.waveform.has_value())
                             ? eval_value(rf.voltage_V, t, lib) : 0.0;
    const double rf_freq = (rf.frequency_Hz.constant_value.has_value() || rf.frequency_Hz.waveform_ref.has_value() || rf.frequency_Hz.waveform.has_value())
                          ? eval_value(rf.frequency_Hz, t, lib) : 0.0;
    const double ac_voltage = (ac.voltage_V.constant_value.has_value() || ac.voltage_V.waveform_ref.has_value() || ac.voltage_V.waveform.has_value())
                             ? eval_value(ac.voltage_V, t, lib) : 0.0;
    const double ac_freq = (ac.frequency_Hz.constant_value.has_value() || ac.frequency_Hz.waveform_ref.has_value() || ac.frequency_Hz.waveform.has_value())
                          ? eval_value(ac.frequency_Hz, t, lib) : 0.0;
    const double dc_quad = (dc.quad_V.constant_value.has_value() || dc.quad_V.waveform_ref.has_value() || dc.quad_V.waveform.has_value())
                          ? eval_value(dc.quad_V, t, lib) : 0.0;
    const double dc_axial = (dc.axial_V.constant_value.has_value() || dc.axial_V.waveform_ref.has_value() || dc.axial_V.waveform.has_value())
                           ? eval_value(dc.axial_V, t, lib) : 0.0;

    if (std::fabs(rf_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double r0_sq = geom.radius_m * geom.radius_m;
        const double omega = 2.0 * M_PI * rf_freq;
        const double U_eff = dc_quad + rf_voltage * std::cos(omega * t);
        E_total.x =  2.0 * pos.x * U_eff / r0_sq;
        E_total.y = -2.0 * pos.y * U_eff / r0_sq;
    }

    if (std::fabs(ac_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double omega_ac = 2.0 * M_PI * ac_freq;
        const double mag = -ac_voltage / geom.radius_m * std::cos(omega_ac * t);
        E_total.x += mag;
    }

    if (std::fabs(dc_axial) > MIN_VOLTAGE_THRESHOLD) {
        const double L = geom.length_m;
        const double alpha = dc_axial / (L * L);
        const double z_centered = pos.z - 0.5 * L;
        E_total.z = -2.0 * alpha * z_centered;
    }

    return E_total;
}

Vec3 AnalyticalFieldModel::compute_ims_field(const Vec3& pos, double t) const {
    Vec3 E_total{0.0, 0.0, 0.0};
    const auto& dc = domain_->fields.dc;
    const auto& rf = domain_->fields.rf;
    const auto& lib = domain_->fields.waveform_library;

    const double dc_axial = (dc.axial_V.constant_value.has_value() || dc.axial_V.waveform_ref.has_value())
                           ? eval_value(dc.axial_V, t, lib) : 0.0;
    const double dc_quad = (dc.quad_V.constant_value.has_value() || dc.quad_V.waveform_ref.has_value())
                          ? eval_value(dc.quad_V, t, lib) : 0.0;
    const double rf_voltage = (rf.voltage_V.constant_value.has_value() || rf.voltage_V.waveform_ref.has_value())
                             ? eval_value(rf.voltage_V, t, lib) : 0.0;
    const double rf_freq = (rf.frequency_Hz.constant_value.has_value() || rf.frequency_Hz.waveform_ref.has_value())
                          ? eval_value(rf.frequency_Hz, t, lib) : 0.0;

    E_total.z = dc_axial / domain_->geometry.length_m;

    if (std::fabs(rf_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double r0_sq = domain_->geometry.radius_m * domain_->geometry.radius_m;
        const double omega = 2.0 * M_PI * rf_freq;
        const double U_eff = dc_quad + rf_voltage * std::cos(omega * t);
        E_total.x =  2.0 * pos.x * U_eff / r0_sq;
        E_total.y = -2.0 * pos.y * U_eff / r0_sq;
    }
    return E_total;
}

Vec3 AnalyticalFieldModel::compute_tof_field(const Vec3& pos) const {
    const double z_local = pos.z - domain_->geometry.origin_m.z;
    const double acc_length = domain_->geometry.acc_length_m;

    if (z_local < acc_length) {
        const double dc_axial = eval_value(domain_->fields.dc.axial_V, 0.0, domain_->fields.waveform_library);
        return Vec3{0.0, 0.0, dc_axial / acc_length};
    }
    return Vec3{0.0, 0.0, 0.0};
}

Vec3 AnalyticalFieldModel::compute_fticr_field(const Vec3& pos) const {
    const double voltage = eval_value(domain_->fields.dc.radial_V, 0.0, domain_->fields.waveform_library);
    const double L = domain_->geometry.length_m;
    const double r = domain_->geometry.radius_m;
    const double d = std::sqrt((L * L / 8.0) + (r * r / 4.0));
    const double factor = voltage / (d * d);
    const double z_center = pos.z - 0.5 * L;
    return Vec3{
         factor * pos.x,
         factor * pos.y,
        -2.0 * factor * z_center
    };
}

Vec3 AnalyticalFieldModel::compute_orbitrap_field(const Vec3& global_pos, double t) const {
    (void)t; // Orbitrap field currently uses DC only
    // Transform position from global to local coordinates
    const Vec3 pos_local = domain_->rotation_global_to_local * (global_pos - domain_->geometry.origin_m);

    const double r_sq = std::max(1e-18, pos_local.x * pos_local.x + pos_local.y * pos_local.y);
    const double r_char = domain_->geometry.radius_char_m;
    const double r_char_sq = r_char * r_char;
    const double C = 1.0 - r_char_sq / r_sq;

    const double voltage = eval_value(domain_->fields.dc.radial_V, 0.0, domain_->fields.waveform_library);
    const double r_in = domain_->geometry.radius_in_m;
    const double r_out = domain_->geometry.radius_out_m;
    const double k = 2.0 * voltage / (r_char_sq * std::log(r_out / r_in)
                                     - 0.5 * (r_out * r_out - r_in * r_in));

    Vec3 E_local{
        0.5 * k * pos_local.x * C,
        0.5 * k * pos_local.y * C,
       -k * pos_local.z
    };
    return domain_->rotation_local_to_global * E_local;
}

Vec3 AnalyticalFieldModel::compute_quadrupole_rf_field(const Vec3& pos, double t) const {
    Vec3 E_total{0.0, 0.0, 0.0};
    const auto& rf = domain_->fields.rf;
    const auto& dc = domain_->fields.dc;
    const auto& geom = domain_->geometry;
    const auto& lib = domain_->fields.waveform_library;

    const double rf_voltage = (rf.voltage_V.constant_value.has_value() || rf.voltage_V.waveform_ref.has_value())
                             ? eval_value(rf.voltage_V, t, lib) : 0.0;
    const double rf_freq = (rf.frequency_Hz.constant_value.has_value() || rf.frequency_Hz.waveform_ref.has_value())
                          ? eval_value(rf.frequency_Hz, t, lib) : 0.0;
    const double dc_quad = (dc.quad_V.constant_value.has_value() || dc.quad_V.waveform_ref.has_value())
                          ? eval_value(dc.quad_V, t, lib) : 0.0;
    const double dc_axial = (dc.axial_V.constant_value.has_value() || dc.axial_V.waveform_ref.has_value())
                           ? eval_value(dc.axial_V, t, lib) : 0.0;

    if (std::fabs(rf_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double r0_sq = geom.radius_m * geom.radius_m;
        const double omega = 2.0 * M_PI * rf_freq;
        const double U_eff = dc_quad + rf_voltage * std::cos(omega * t);
        E_total.x = -2.0 * pos.x * U_eff / r0_sq;
        E_total.y = 2.0 * pos.y * U_eff / r0_sq;
    }

    if (std::fabs(dc_axial) > MIN_VOLTAGE_THRESHOLD) {
        E_total.z = dc_axial / geom.length_m;
    }

    return E_total;
}

} // namespace ICARION::config
