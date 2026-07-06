// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "IFieldModel.h"
#include "DomainConfig.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace ICARION::config {

/**
 * @brief TIMS axial field model.
 *
 * Evaluates E_z(z,t) = (1 - f(t)) * E_initial(z) + f(t) * E_final(z)
 * in local domain coordinates. Optional RF/DC radial confinement is added from
 * the standard IMS field settings.
 */
class TIMSAxialFieldModel : public IFieldModel {
public:
    explicit TIMSAxialFieldModel(const DomainConfig& domain)
        : domain_(&domain) {}

    Vec3 E(const Vec3& global_pos, double t) const override {
        if (!domain_) {
            return Vec3{0.0, 0.0, 0.0};
        }

        const Vec3 pos_local = domain_->rotation_global_to_local * (global_pos - domain_->geometry.origin_m);
        Vec3 e_local{0.0, 0.0, 0.0};

        const auto& tims = domain_->fields.tims;
        const double ramp_fraction = evaluate_ramp_fraction(t);
        const double e_initial = sample_profile(
            pos_local.z, tims.z_positions_m, tims.axial_field_initial_profile_V_m, tims.axial_field_initial_uniform_V_m);
        const double e_final = sample_profile(
            pos_local.z, tims.z_positions_m, tims.axial_field_final_profile_V_m, tims.axial_field_final_uniform_V_m);
        e_local.z = (1.0 - ramp_fraction) * e_initial + ramp_fraction * e_final;

        add_radial_confinement(pos_local, t, e_local);
        return domain_->rotation_local_to_global * e_local;
    }

private:
    const DomainConfig* domain_ = nullptr;

    static double eval_value(const ValueOrWaveform& value,
                             double t,
                             const std::map<std::string, Waveform>& waveform_library) {
        if (value.constant_value.has_value()) {
            return value.constant_value.value();
        }
        return value.evaluate(t, waveform_library);
    }

    static bool has_value(const ValueOrWaveform& value) {
        return value.constant_value.has_value() ||
               value.waveform_ref.has_value() ||
               value.waveform.has_value();
    }

    double evaluate_ramp_fraction(double t) const {
        const auto& tims = domain_->fields.tims;
        if (tims.ramp_fraction_defined) {
            const double raw = eval_value(tims.ramp_fraction, t, domain_->fields.waveform_library);
            return std::clamp(raw, 0.0, 1.0);
        }

        const double t0 = tims.ramp_start_s;
        const double t1 = tims.ramp_end_s;
        if (t <= t0) {
            return 0.0;
        }
        if (t >= t1) {
            return 1.0;
        }

        const double duration = t1 - t0;
        if (duration <= 0.0) {
            return 1.0;
        }

        const double rel = t - t0;
        if (tims.ramp_mode == "exponential") {
            const double tau = std::max(tims.ramp_tau_s, 1e-30);
            const double denom = 1.0 - std::exp(-duration / tau);
            if (std::abs(denom) < 1e-15) {
                return std::clamp(rel / duration, 0.0, 1.0);
            }
            return std::clamp((1.0 - std::exp(-rel / tau)) / denom, 0.0, 1.0);
        }

        return std::clamp(rel / duration, 0.0, 1.0);
    }

    static double sample_profile(double z_local_m,
                                 const std::vector<double>& z_nodes_m,
                                 const std::vector<double>& values_V_m,
                                 double uniform_fallback_V_m) {
        if (z_nodes_m.empty() || values_V_m.empty() || z_nodes_m.size() != values_V_m.size()) {
            return uniform_fallback_V_m;
        }
        if (z_local_m <= z_nodes_m.front()) {
            return values_V_m.front();
        }
        if (z_local_m >= z_nodes_m.back()) {
            return values_V_m.back();
        }

        auto it = std::lower_bound(z_nodes_m.begin(), z_nodes_m.end(), z_local_m);
        const size_t i1 = static_cast<size_t>(std::distance(z_nodes_m.begin(), it));
        const size_t i0 = i1 - 1;

        const double z0 = z_nodes_m[i0];
        const double z1 = z_nodes_m[i1];
        const double dz = z1 - z0;
        if (dz <= 0.0) {
            return values_V_m[i0];
        }

        const double alpha = (z_local_m - z0) / dz;
        return values_V_m[i0] + alpha * (values_V_m[i1] - values_V_m[i0]);
    }

    void add_radial_confinement(const Vec3& pos_local, double t, Vec3& e_local) const {
        const auto& rf = domain_->fields.rf;
        const auto& dc = domain_->fields.dc;
        const auto& lib = domain_->fields.waveform_library;

        const double rf_voltage = has_value(rf.voltage_V) ? eval_value(rf.voltage_V, t, lib) : 0.0;
        const double rf_freq = has_value(rf.frequency_Hz) ? eval_value(rf.frequency_Hz, t, lib) : 0.0;
        const double dc_quad = has_value(dc.quad_V) ? eval_value(dc.quad_V, t, lib) : 0.0;

        const double r0_sq = domain_->geometry.radius_m * domain_->geometry.radius_m;
        if (r0_sq <= 0.0) {
            return;
        }

        const double omega = 2.0 * M_PI * rf_freq;
        const double u_eff = dc_quad + rf_voltage * std::cos(omega * t);
        e_local.x += -2.0 * pos_local.x * u_eff / r0_sq;
        e_local.y += 2.0 * pos_local.y * u_eff / r0_sq;
    }
};

} // namespace ICARION::config
